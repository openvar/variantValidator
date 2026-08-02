import re
import copy
import logging
import vvhgvs
from vvhgvs.assemblymapper import AssemblyMapper
import vvhgvs.validator
from . import vvMixinInit
from . import seq_data
from . import hgvs_utils, hgvs_position_utils
from . import expanded_repeats
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from . import utils as fn
import sys
import json

from vvhgvs.exceptions import HGVSError, HGVSDataNotAvailableError, HGVSUnsupportedOperationError, \
     HGVSInvalidVariantError
from vvhgvs.enums import Datum # needed to handle r-> n mapping without re-parsing posedit
from VariantValidator.modules.hgvs_utils import hgvs_delins_parts_to_hgvs_obj, hgvs_dup_to_delins,\
        hgvs_obj_from_existing_edit
from VariantValidator.modules.transcript_map_data import TranscriptMapData

logger = logging.getLogger(__name__)


class AlleleSyntaxError(Exception):
    pass


class Mixin(vvMixinInit.Mixin):
    """
    Converters that use the Validator configuration.
    """

    def _expand_ref(self, ac, start, stop):
        """
        Fetch the first and last bases of a sequence interval.

        For intervals <= 1 kb, a single sequence fetch is faster than two
        separate 1 bp fetches. Coordinates are 0-based SeqRepo coordinates.
        """
        if stop - start > 1000:
            pre_base = self.sf.fetch_seq(ac, start, start + 1)
            post_base = self.sf.fetch_seq(ac, stop - 1, stop)
            return pre_base, post_base

        span = self.sf.fetch_seq(ac, start, stop)
        return span[0], span[-1]

    def coding(self, variant):
        """
        Return a transcript variant as a c. HGVS object where applicable.
        """
        if isinstance(variant, str):
            if ':c.' not in variant and ':n.' not in variant:
                return None
            variant = self.hp.parse_hgvs_variant(variant)

        if variant.type == 'n':
            return self.vm.n_to_c(variant)

        if variant.type == 'c':
            return variant

        return None

    def genomic(self, variant, evm, primary_assembly, vv_variant):
        """
        Return a variant as a genomic HGVS object where applicable.
        """
        hn = vv_variant.hn

        logger.info("Map %s to genomic position", variant)
        logger.info("Primary assembly: %s", primary_assembly)

        if not isinstance(variant, str):
            logger.info("Variant %s is not a string", variant)

            if variant.type in ('c', 'n'):
                try:
                    var_g = self.myevm_t_to_g(
                        variant,
                        evm,
                        primary_assembly,
                        hn,
                        vv_variant,
                    )
                except vvhgvs.exceptions.HGVSError as e:
                    logger.info("HGVS error: %s", e)
                    return 'error ' + str(e)

                logger.info("Variant %s mapped to %s", variant, var_g)
                return var_g

            if variant.type == 'g':
                return variant

            return None

        if ':c.' in variant or ':n.' in variant:
            hgvs_var = self.hp.parse_hgvs_variant(variant)

            try:
                return self.myevm_t_to_g(
                    hgvs_var,
                    evm,
                    primary_assembly,
                    hn,
                    vv_variant,
                )
            except vvhgvs.exceptions.HGVSError as e:
                return 'error ' + str(e)

        if ':g.' in variant:
            return self.hp.parse_hgvs_variant(variant)

        return None

    def myevm_t_to_g(
            self,
            hgvs_c,
            no_norm_evm,
            primary_assembly,
            hn,
            variant,
            reset_g_origin=False
    ):
        """
        Enhanced transcript-to-genome mapping using evm.

        Handles transcript positions affected by transcript/genome alignment gaps
        and attempts alternative genomic mappings when the preferred mapping is
        unavailable.

        Mapping preference:
            NC_ on requested assembly
            NC_ on alternative assembly
            NT_ on requested assembly
            NT_ on alternative assembly
            NW_ on requested assembly
            NW_ on alternative assembly
            NG_

        Requires a parsed c. or n. HGVS object and returns a parsed g. object.
        """
        alt_aln_method = self.alt_aln_method
        stored_hgvs_c = copy.deepcopy(hgvs_c)
        expand_out = False

        # ------------------------------------------------------------------
        # Local helpers
        # ------------------------------------------------------------------

        def gap_pre_tx_corrections(hgvs_cg):
            """
            Expand a non-intronic transcript variant before mapping across a
            gapped transcript/genome alignment.
            """
            if hgvs_cg.posedit.edit.type not in (
                    'identity', 'del', 'delins', 'dup', 'sub', 'ins', 'inv'
            ):
                return hgvs_cg

            # Gap handling is performed in n. coordinates where possible.
            if hgvs_cg.type == 'c':
                hgvs_cg = no_norm_evm.c_to_n(hgvs_cg)

            try:
                hn.normalize(hgvs_cg)

            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)

                if (
                        'intronic variant' not in error
                        and
                        'Length implied by coordinates must equal sequence deletion length'
                        in error
                        and hgvs_cg.ac.startswith('NR_')
                ):
                    hgvs_cg.posedit.pos.end.base = (
                            hgvs_c.posedit.pos.start.base
                            + len(hgvs_c.posedit.edit.ref)
                            - 1
                    )

            # Do not expand intronic variants.
            if not hgvs_position_utils.either_position_is_intronic(hgvs_cg):
                try:
                    hgvs_t = copy.deepcopy(hgvs_cg)
                    edit_type = hgvs_t.posedit.edit.type

                    if edit_type == 'inv':
                        inv_alt = self.revcomp(
                            hgvs_t.posedit.edit.ref
                        )

                        pre_base, post_base = self._expand_ref(
                            hgvs_t.ac,
                            hgvs_t.posedit.pos.start.base - 2,
                            hgvs_t.posedit.pos.end.base + 1
                        )

                        hgvs_t = hgvs_delins_parts_to_hgvs_obj(
                            hgvs_t.ac,
                            hgvs_t.type,
                            hgvs_t.posedit.pos.start.base - 1,
                            (
                                    pre_base
                                    + hgvs_t.posedit.edit.ref
                                    + post_base
                            ),
                            pre_base + inv_alt + post_base,
                            end=hgvs_t.posedit.pos.end.base + 1,
                            offset_pos=True
                        )

                    elif edit_type == 'dup':
                        pre_base, post_base = self._expand_ref(
                            hgvs_t.ac,
                            hgvs_t.posedit.pos.start.base - 2,
                            hgvs_t.posedit.pos.end.base + 1
                        )

                        hgvs_t = hgvs_delins_parts_to_hgvs_obj(
                            hgvs_t.ac,
                            hgvs_t.type,
                            hgvs_t.posedit.pos.start.base - 1,
                            (
                                    pre_base
                                    + hgvs_t.posedit.edit.ref
                                    + post_base
                            ),
                            (
                                    pre_base
                                    + hgvs_t.posedit.edit.ref
                                    + hgvs_t.posedit.edit.ref
                                    + post_base
                            ),
                            offset_pos=True
                        )

                    elif edit_type == 'ins':
                        # ins -> delins changes between-coordinate HGVS
                        # representation to an inclusive interval.
                        ins_ref = self.sf.fetch_seq(
                            hgvs_t.ac,
                            hgvs_t.posedit.pos.start.base - 2,
                            hgvs_t.posedit.pos.end.base + 1
                        )

                        ins_alt = (
                                ins_ref[:2]
                                + hgvs_t.posedit.edit.alt
                                + ins_ref[-2:]
                        )

                        hgvs_t = hgvs_delins_parts_to_hgvs_obj(
                            hgvs_t.ac,
                            hgvs_t.type,
                            hgvs_t.posedit.pos.start.base - 1,
                            ins_ref,
                            ins_alt,
                            end=hgvs_t.posedit.pos.end.base + 1,
                            offset_pos=True
                        )

                    else:
                        if hgvs_t.posedit.edit.alt is None:
                            hgvs_t.posedit.edit.alt = ''

                        pre_base, post_base = self._expand_ref(
                            hgvs_t.ac,
                            hgvs_t.posedit.pos.start.base - 2,
                            hgvs_t.posedit.pos.end.base + 1
                        )

                        hgvs_t = hgvs_delins_parts_to_hgvs_obj(
                            hgvs_t.ac,
                            hgvs_t.type,
                            hgvs_t.posedit.pos.start.base - 1,
                            (
                                    pre_base
                                    + hgvs_t.posedit.edit.ref
                                    + post_base
                            ),
                            (
                                    pre_base
                                    + hgvs_t.posedit.edit.alt
                                    + post_base
                            ),
                            end=hgvs_t.posedit.pos.end.base + 1,
                            offset_pos=True
                        )

                    hgvs_cg = copy.deepcopy(hgvs_t)

                # Historical defensive behaviour around malformed/gapped HGVS.
                except Exception:
                    pass

            # Convert back to c. where possible.
            try:
                hgvs_cg = no_norm_evm.n_to_c(hgvs_cg)
            except vvhgvs.exceptions.HGVSError:
                hgvs_cg = copy.deepcopy(stored_hgvs_c)

            # Ensure expansion has not crossed an exon/intron boundary.
            hgvs_check_boundaries = copy.deepcopy(hgvs_cg)

            try:
                hn.normalize(hgvs_check_boundaries)
            except vvhgvs.exceptions.HGVSError as e:
                if 'spanning the exon-intron boundary' in str(e):
                    hgvs_cg = copy.deepcopy(stored_hgvs_c)

            # Identity variants require an additional reference-only check.
            if hgvs_check_boundaries.posedit.edit.type == 'identity':
                hgvs_reform_ident = hgvs_delins_parts_to_hgvs_obj(
                    hgvs_cg.ac,
                    stored_hgvs_c.type,
                    hgvs_cg.posedit.pos,
                    hgvs_cg.posedit.edit.ref,
                    '',
                    offset_pos=True
                )

                try:
                    hn.normalize(hgvs_reform_ident)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)

                    if (
                            'spanning the exon-intron boundary' in error
                            or 'Normalization of intronic variants' in error
                    ):
                        hgvs_cg = copy.deepcopy(stored_hgvs_c)

            return hgvs_cg

        def map_to_genomic_ac(transcript_variant, genomic_ac):
            """
            Map to a specific genomic accession, applying transcript gap
            compensation where required.
            """
            if variant.map_dat.is_gapped_map(
                    hgvs_c.ac,
                    genomic_ac,
                    hdp=self.hdp
            ):
                logger.debug(
                    "gap_compensation_myevm enabled for %s against %s",
                    hgvs_c.ac,
                    genomic_ac
                )

                corrected = gap_pre_tx_corrections(transcript_variant)

                return super(
                    AssemblyMapper,
                    no_norm_evm
                ).t_to_g(
                    corrected,
                    genomic_ac,
                    alt_aln_method=alt_aln_method
                )

            return super(
                AssemblyMapper,
                no_norm_evm
            ).t_to_g(
                transcript_variant,
                genomic_ac,
                alt_aln_method=alt_aln_method
            )

        def rebuild_gap_variant(
                transcript_gap_n,
                transcript_gap_alt_n,
                genomic_ac
        ):
            """
            Reconstruct a transcript delins across a transcript/genome alignment
            gap and map it back to the genome.

            This replaces the two historically duplicated reconstruction blocks.
            """
            # Duplications do not expose alt in the same way as NARefAlt.
            try:
                alt = transcript_gap_alt_n.posedit.edit.alt
            except AttributeError:
                transcript_gap_n = hgvs_dup_to_delins(
                    transcript_gap_n
                )
                transcript_gap_alt_n = hgvs_dup_to_delins(
                    transcript_gap_alt_n
                )
                alt = transcript_gap_alt_n.posedit.edit.alt

            if alt is None:
                alternate_bases = (
                        ['X']
                        * len(transcript_gap_alt_n.posedit.edit.ref)
                )
            else:
                alternate_bases = list(alt)

            ref_start = transcript_gap_n.posedit.pos.start.base
            alt_start = transcript_gap_alt_n.posedit.pos.start.base

            ref_base_dict = {
                ref_start + index: base
                for index, base in enumerate(
                    transcript_gap_n.posedit.edit.ref
                )
            }

            alt_base_dict = {}

            for position in range(
                    transcript_gap_alt_n.posedit.pos.start.base,
                    transcript_gap_alt_n.posedit.pos.end.base + 1
            ):
                if position == alt_start:
                    alt_base_dict[position] = ''.join(
                        alternate_bases
                    )
                else:
                    alt_base_dict[position] = 'X'

            alternate_sequence_bases = []

            for position in range(
                    transcript_gap_n.posedit.pos.start.base,
                    transcript_gap_n.posedit.pos.end.base + 1
            ):
                if position in alt_base_dict:
                    alternate_sequence_bases.append(
                        alt_base_dict[position]
                    )
                elif position in ref_base_dict:
                    alternate_sequence_bases.append(
                        ref_base_dict[position]
                    )

            transcript_gap_n.posedit.edit.alt = ''.join(
                alternate_sequence_bases
            ).replace('X', '')

            try:
                transcript_gap_variant = self.vm.n_to_c(
                    transcript_gap_n
                )
            except vvhgvs.exceptions.HGVSError:
                transcript_gap_variant = transcript_gap_n

            try:
                mapped = self.vm.t_to_g(
                    transcript_gap_variant,
                    genomic_ac,
                    alt_aln_method
                )

                return hn.normalize(mapped)

            except vvhgvs.exceptions.HGVSError as e:
                if str(e) != 'base start position must be <= end position':
                    raise

            # Variant must be expanded one base at each side before it can
            # map back across the genomic gap.
            pre_base, post_base = self._expand_ref(
                transcript_gap_n.ac,
                transcript_gap_n.posedit.pos.start.base - 2,
                transcript_gap_n.posedit.pos.end.base + 1
            )

            transcript_gap_n.posedit.pos.start.base -= 1
            transcript_gap_n.posedit.pos.end.base += 1

            transcript_gap_n.posedit.edit.ref = (
                    pre_base
                    + transcript_gap_n.posedit.edit.ref
                    + post_base
            )

            transcript_gap_n.posedit.edit.alt = (
                    pre_base
                    + transcript_gap_n.posedit.edit.alt
                    + post_base
            )

            try:
                transcript_gap_variant = self.vm.n_to_c(
                    transcript_gap_n
                )
            except vvhgvs.exceptions.HGVSError:
                transcript_gap_variant = transcript_gap_n

            mapped = self.vm.t_to_g(
                transcript_gap_variant,
                genomic_ac,
                alt_aln_method
            )

            return hn.normalize(mapped)

        # ------------------------------------------------------------------
        # Determine available mappings
        # ------------------------------------------------------------------

        mapping_options = variant.map_dat.mapping_options(
            hgvs_c.ac,
            hdp=self.hdp
        )

        hgvs_genomic = None
        gap_corrected_hgvs_c = None
        attempted_mapping_errors = []

        if reset_g_origin or not hgvs_c.rel_ac:
            hgvs_c.rel_ac = ''

            for option in mapping_options:
                genomic_ac = option[1]

                if (
                        genomic_ac.startswith('NC_')
                        and seq_data.is_supported_for_mapping(
                    genomic_ac,
                    primary_assembly
                )
                ):
                    # Preserve historical behaviour: the final supported
                    # NC_ encountered becomes the relative accession.
                    hgvs_c.rel_ac = genomic_ac

        try:
            hn.normalize(hgvs_c)
        except vvhgvs.exceptions.HGVSError:
            pass

        # Preserve the preferred failed-normalisation mapping as a fallback.
        norm_f_hgvs_genomic = None

        # ------------------------------------------------------------------
        # First try the existing relative genomic origin
        # ------------------------------------------------------------------

        if hgvs_c.rel_ac:
            try:
                if variant.map_dat.is_gapped_map(
                        hgvs_c.ac,
                        hgvs_c.rel_ac,
                        hdp=self.hdp
                ):
                    logger.debug(
                        "gap_compensation_myevm enabled for %s against %s",
                        hgvs_c.ac,
                        hgvs_c.rel_ac
                    )

                    gap_corrected_hgvs_c = gap_pre_tx_corrections(
                        hgvs_c
                    )

                    hgvs_genomic = no_norm_evm.t_to_g(
                        gap_corrected_hgvs_c
                    )

                else:
                    hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)

                hn.normalize(hgvs_genomic)

            except vvhgvs.exceptions.HGVSError:
                norm_f_hgvs_genomic = hgvs_genomic
                hgvs_genomic = None

                if not mapping_options:
                    raise HGVSDataNotAvailableError(
                        "No alignment data between the specified transcript "
                        "reference sequence and any GRCh37 and GRCh38 genomic "
                        "reference sequences (including alternate chromosome "
                        "assemblies, patches and RefSeqGenes) are available."
                    )

        # ------------------------------------------------------------------
        # Search mapping options
        # ------------------------------------------------------------------

        for genomic_ac_type in ('NC_', 'NT_', 'NW_', 'NG_'):
            if hgvs_genomic is not None:
                break

            current_options = []

            for option in mapping_options:
                genomic_ac = option[1]
                alignment_method = option[2]

                if alignment_method.startswith('blat'):
                    continue

                if genomic_ac.startswith(genomic_ac_type):
                    current_options.append(genomic_ac)

            # Requested assembly first. NG_ is not chromosomal.
            if genomic_ac_type != 'NG_':
                for genomic_ac in current_options:
                    if not seq_data.is_supported_for_mapping(
                            genomic_ac,
                            primary_assembly
                    ):
                        continue

                    try:
                        hgvs_genomic = map_to_genomic_ac(
                            hgvs_c,
                            genomic_ac
                        )

                    except Exception as e:
                        attempted_mapping_errors.append(
                            f"{e}/{hgvs_c.ac}/{genomic_ac}~"
                        )
                        continue

                    try:
                        hn.normalize(hgvs_genomic)

                    except Exception:
                        if norm_f_hgvs_genomic is None:
                            norm_f_hgvs_genomic = hgvs_genomic

                        hgvs_genomic = None
                        continue

                    break

            if hgvs_genomic is not None:
                break

            # Then try mappings outside the requested assembly.
            for genomic_ac in current_options:
                if (
                        genomic_ac_type != 'NG_'
                        and seq_data.is_supported_for_mapping(
                    genomic_ac,
                    primary_assembly
                )
                ):
                    continue

                try:
                    hgvs_genomic = map_to_genomic_ac(
                        hgvs_c,
                        genomic_ac
                    )

                except Exception as e:
                    attempted_mapping_errors.append(
                        f"{e}/{hgvs_c.ac}/{genomic_ac}~"
                    )
                    continue

                try:
                    hn.normalize(hgvs_genomic)

                except Exception:
                    if norm_f_hgvs_genomic is None:
                        norm_f_hgvs_genomic = hgvs_genomic

                    hgvs_genomic = None
                    continue

                break

        if hgvs_genomic is None and norm_f_hgvs_genomic is not None:
            hgvs_genomic = norm_f_hgvs_genomic

        if hgvs_genomic is None:
            logger.debug("HGVS data not available error")

            raise HGVSDataNotAvailableError(
                ''.join(attempted_mapping_errors)
            )

        # ------------------------------------------------------------------
        # Determine whether gap expansion is required
        # ------------------------------------------------------------------

        gapped_mapping = variant.map_dat.is_gapped_map(
            hgvs_c.ac,
            hgvs_genomic.ac,
            hdp=self.hdp
        )

        expand_out = (
                gapped_mapping
                and not hgvs_position_utils.either_position_is_intronic(
            hgvs_c
        )
        )

        if (
                hgvs_c.posedit.edit.type == 'identity'
                and hgvs_genomic.posedit.edit.type == 'delins'
                and hgvs_genomic.posedit.edit.alt == ''
                and not expand_out
        ):
            hgvs_genomic.posedit.edit.alt = (
                hgvs_genomic.posedit.edit.ref
            )

        # ------------------------------------------------------------------
        # Correct malformed genomic insertion mappings across gaps
        # ------------------------------------------------------------------

        if hgvs_genomic.posedit.edit.type == 'ins' and gapped_mapping:
            try:
                hgvs_genomic = hn.normalize(hgvs_genomic)

            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)

                if error == 'insertion length must be 1':
                    ref = self.sf.fetch_seq(
                        hgvs_genomic.ac,
                        hgvs_genomic.posedit.pos.start.base - 1,
                        hgvs_genomic.posedit.pos.end.base
                    )

                    hgvs_genomic.posedit.edit.ref = ref
                    hgvs_genomic.posedit.edit.alt = (
                            ref[:1]
                            + hgvs_genomic.posedit.edit.alt
                            + ref[-1:]
                    )

                    hgvs_genomic = hn.normalize(hgvs_genomic)

                elif error == 'base start position must be <= end position':
                    start = hgvs_genomic.posedit.pos.start.base
                    end = hgvs_genomic.posedit.pos.end.base

                    hgvs_genomic.posedit.pos.start.base = end
                    hgvs_genomic.posedit.pos.end.base = start

                    hgvs_genomic = hn.normalize(hgvs_genomic)

        # ------------------------------------------------------------------
        # Restore references required by expanded gap descriptions
        # ------------------------------------------------------------------

        if (
                stored_hgvs_c.posedit.edit.ref in ('', None)
                and expand_out
        ):
            if stored_hgvs_c.type == 'c':
                stored_hgvs_n = self.vm.c_to_n(
                    stored_hgvs_c
                )
            else:
                stored_hgvs_n = stored_hgvs_c

            stored_hgvs_c.posedit.edit.ref = self.sf.fetch_seq(
                stored_hgvs_n.ac,
                stored_hgvs_n.posedit.pos.start.base - 1,
                stored_hgvs_n.posedit.pos.end.base
            )

        if (
                hgvs_genomic.posedit.edit.ref in ('', None)
                and expand_out
                and hgvs_genomic.posedit.edit.type == 'ins'
        ):
            stored_ref = self.sf.fetch_seq(
                hgvs_genomic.ac,
                hgvs_genomic.posedit.pos.start.base - 1,
                hgvs_genomic.posedit.pos.end.base
            )

            hgvs_genomic.posedit.edit.ref = stored_ref
            hgvs_genomic.posedit.edit.alt = (
                    stored_ref[:1]
                    + hgvs_genomic.posedit.edit.alt
                    + stored_ref[-1:]
            )

        # ------------------------------------------------------------------
        # Gap flank / within-gap reconstruction
        # ------------------------------------------------------------------

        if expand_out:
            nr_genomic = self.nr_vm.t_to_g(
                hgvs_c,
                hgvs_genomic.ac,
                alt_aln_method
            )

            try:
                hn.normalize(nr_genomic)

            except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                error_type_1 = str(e)

                if (
                        'Length implied by coordinates must equal sequence deletion length'
                        in error_type_1
                        or error_type_1
                        == 'base start position must be <= end position'
                ):
                    genomic_gap_variant = None

                    if (
                            'Length implied by coordinates must equal sequence deletion length'
                            in error_type_1
                    ):
                        logger.info(
                            "Variant is proximal to the flank of a genomic gap"
                        )

                        genomic_gap_variant = self.vm.t_to_g(
                            stored_hgvs_c,
                            hgvs_genomic.ac,
                            alt_aln_method
                        )

                        try:
                            hn.normalize(genomic_gap_variant)

                        except vvhgvs.exceptions.HGVSInvalidVariantError as e2:
                            if (
                                    'base start position must be <= end position'
                                    in str(e2)
                                    and
                                    'Length implied by coordinates must equal'
                                    in error_type_1
                            ):
                                make_gen_var = copy.copy(nr_genomic)

                                make_gen_var.posedit.edit.ref = (
                                    self.sf.fetch_seq(
                                        nr_genomic.ac,
                                        nr_genomic.posedit.pos.start.base - 1,
                                        nr_genomic.posedit.pos.end.base
                                    )
                                )

                                genomic_gap_variant = make_gen_var
                                error_type_1 = None

                        else:
                            genomic_gap_variant = self.nr_vm.t_to_g(
                                hgvs_c,
                                hgvs_genomic.ac,
                                alt_aln_method
                            )

                    if (
                            error_type_1
                            == 'base start position must be <= end position'
                    ):
                        logger.info(
                            "Variant is fully within a genomic gap"
                        )

                        genomic_gap_variant = self.vm.t_to_g(
                            stored_hgvs_c,
                            hgvs_genomic.ac,
                            alt_aln_method
                        )

                    try:
                        hn.normalize(genomic_gap_variant)

                    except Exception as gap_error:
                        gap_error_text = str(gap_error)

                        if (
                                gap_error_text
                                == 'base start position must be <= end position'
                        ):
                            gap_start = (
                                genomic_gap_variant.posedit.pos.end.base
                            )
                            gap_end = (
                                genomic_gap_variant.posedit.pos.start.base
                            )

                            genomic_gap_variant.posedit.pos.start.base = (
                                gap_start
                            )
                            genomic_gap_variant.posedit.pos.end.base = (
                                gap_end
                            )

                        if (
                                'Length implied by coordinates must equal sequence deletion length'
                                in gap_error_text
                        ):
                            logger.info(
                                "Variant is on the flank of a genomic gap "
                                "but not within the gap"
                            )

                            try:
                                try:
                                    norm_stored_c = hn.normalize(
                                        stored_hgvs_c
                                    )
                                except HGVSUnsupportedOperationError:
                                    norm_stored_c = stored_hgvs_c

                                if norm_stored_c.posedit.edit.type in (
                                        'sub',
                                        'identity'
                                ):
                                    flank_hgvs_genomic = self.vm.t_to_g(
                                        norm_stored_c,
                                        genomic_gap_variant.ac,
                                        alt_aln_method
                                    )

                                    self.vr.validate(
                                        flank_hgvs_genomic
                                    )

                                    # Preserve the historical special case:
                                    # a one-base transcript gap substitution
                                    # continues through gap reconstruction.
                                    transcript_gap_sub = (
                                            flank_hgvs_genomic.posedit.edit.type
                                            == 'sub'
                                            and
                                            norm_stored_c.posedit.edit.type
                                            == 'sub'
                                            and
                                            stored_hgvs_c.posedit.edit.type
                                            == 'sub'
                                            and
                                            len(
                                                genomic_gap_variant.posedit.edit.ref
                                            ) + 1
                                            == (
                                                    genomic_gap_variant.posedit.pos.end.base
                                                    - genomic_gap_variant.posedit.pos.start.base
                                                    + 1
                                            )
                                    )

                                    if not transcript_gap_sub:
                                        return flank_hgvs_genomic

                            except HGVSInvalidVariantError:
                                pass

                            genomic_gap_variant.posedit.pos.start.base -= 1
                            genomic_gap_variant.posedit.pos.end.base += 1
                            genomic_gap_variant.posedit.edit.ref = ''

                            stored_hgvs_c = copy.deepcopy(
                                hgvs_c
                            )

                        try:
                            genomic_gap_variant.posedit.edit.alt = ''
                        except AttributeError:
                            pass

                        genomic_gap_variant = hn.normalize(
                            genomic_gap_variant
                        )

                        transcript_gap_variant = self.vm.g_to_t(
                            genomic_gap_variant,
                            hgvs_c.ac,
                            alt_aln_method=alt_aln_method
                        )

                        if (
                                'Length implied by coordinates must equal sequence deletion length'
                                not in gap_error_text
                        ):
                            try:
                                transcript_gap_variant = hn.normalize(
                                    transcript_gap_variant
                                )
                            except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                                pass

                        try:
                            transcript_gap_n = no_norm_evm.c_to_n(
                                transcript_gap_variant
                            )
                            transcript_gap_alt_n = no_norm_evm.c_to_n(
                                stored_hgvs_c
                            )

                        except vvhgvs.exceptions.HGVSError:
                            transcript_gap_n = transcript_gap_variant
                            transcript_gap_alt_n = stored_hgvs_c

                        hgvs_genomic = rebuild_gap_variant(
                            transcript_gap_n,
                            transcript_gap_alt_n,
                            hgvs_genomic.ac
                        )

                        # Bypass the later expansion correction.
                        expand_out = False

        # ------------------------------------------------------------------
        # Correct expanded genomic descriptions
        # ------------------------------------------------------------------

        if (
                hgvs_c != stored_hgvs_c
                and expand_out
                and gapped_mapping
        ):
            genomic_ref = hgvs_genomic.posedit.edit.ref
            stored_ref = stored_hgvs_c.posedit.edit.ref

            if genomic_ref is None:
                if (
                        hgvs_genomic.posedit.edit.alt is not None
                        and len(hgvs_genomic.posedit.edit.alt) > 2
                ):
                    hgvs_genomic.posedit.edit.alt = (
                        hgvs_genomic.posedit.edit.alt[1:-1]
                    )

            elif len(genomic_ref) == len(stored_ref) + 2:
                hgvs_genomic.posedit.pos.start.base += 1
                hgvs_genomic.posedit.pos.end.base -= 1

                hgvs_genomic.posedit.edit.ref = genomic_ref[1:-1]

                if hgvs_genomic.posedit.edit.alt is not None:
                    hgvs_genomic.posedit.edit.alt = (
                        hgvs_genomic.posedit.edit.alt[1:-1]
                    )

            elif len(genomic_ref) == 2:
                hn.normalize(hgvs_genomic)

            elif len(genomic_ref) <= 1:
                genomic_gap_variant = self.vm.t_to_g(
                    stored_hgvs_c,
                    hgvs_genomic.ac,
                    alt_aln_method
                )

                try:
                    hn.normalize(genomic_gap_variant)

                except Exception as gap_error:
                    if (
                            str(gap_error)
                            == 'base start position must be <= end position'
                    ):
                        gap_start = (
                            genomic_gap_variant.posedit.pos.end.base
                        )
                        gap_end = (
                            genomic_gap_variant.posedit.pos.start.base
                        )

                        genomic_gap_variant.posedit.pos.start.base = (
                            gap_start
                        )
                        genomic_gap_variant.posedit.pos.end.base = (
                            gap_end
                        )

                    try:
                        genomic_gap_variant.posedit.edit.alt = ''
                    except AttributeError:
                        pass

                    genomic_gap_variant = hn.normalize(
                        genomic_gap_variant
                    )

                    transcript_gap_variant = self.vm.g_to_t(
                        genomic_gap_variant,
                        hgvs_c.ac,
                        alt_aln_method=alt_aln_method
                    )

                    transcript_gap_variant = hn.normalize(
                        transcript_gap_variant
                    )

                    try:
                        transcript_gap_n = no_norm_evm.c_to_n(
                            transcript_gap_variant
                        )
                        transcript_gap_alt_n = no_norm_evm.c_to_n(
                            stored_hgvs_c
                        )

                    except vvhgvs.exceptions.HGVSError:
                        transcript_gap_n = transcript_gap_variant
                        transcript_gap_alt_n = stored_hgvs_c

                    hgvs_genomic = rebuild_gap_variant(
                        transcript_gap_n,
                        transcript_gap_alt_n,
                        hgvs_genomic.ac
                    )

        # ------------------------------------------------------------------
        # Exon/exon insertion rescue
        # ------------------------------------------------------------------

        if (
                hgvs_c.posedit.edit.type == 'ins'
                and not hgvs_position_utils.either_position_is_intronic(
            hgvs_c
        )
        ):
            try:
                hn.normalize(hgvs_genomic)

            except vvhgvs.exceptions.HGVSError as e:
                if str(e) == 'insertion length must be 1':
                    if hgvs_c.type == 'c':
                        hgvs_t = self.vm.c_to_n(hgvs_c)
                    else:
                        hgvs_t = copy.copy(hgvs_c)

                    ins_ref = self.sf.fetch_seq(
                        hgvs_t.ac,
                        hgvs_t.posedit.pos.start.base - 1,
                        hgvs_t.posedit.pos.end.base
                    )

                    ins_alt = (
                            ins_ref[:1]
                            + hgvs_t.posedit.edit.alt
                            + ins_ref[-1:]
                    )

                    hgvs_t = hgvs_delins_parts_to_hgvs_obj(
                        hgvs_t.ac,
                        hgvs_t.type,
                        hgvs_t.posedit.pos.start.base,
                        ins_ref,
                        ins_alt,
                        end=hgvs_t.posedit.pos.end.base,
                        offset_pos=True
                    )

                    try:
                        hgvs_c = self.vm.n_to_c(hgvs_t)
                    except vvhgvs.exceptions.HGVSError:
                        hgvs_c = copy.copy(hgvs_t)

                    try:
                        hgvs_genomic = no_norm_evm.t_to_g(
                            hgvs_c
                        )
                    except vvhgvs.exceptions.HGVSError:
                        pass

        # ------------------------------------------------------------------
        # Non-gapped dup/inv reference correction
        # ------------------------------------------------------------------

        if (
                hgvs_c.posedit.edit.type in ('dup', 'inv')
                and not gapped_mapping
                and not hgvs_position_utils.either_position_is_intronic(
                    hgvs_c
                )
        ):
            old_gen_ref = hgvs_c.posedit.edit.ref

            # Normalisation will often populate a missing reference.
            if not old_gen_ref:
                try:
                    hn.normalize(hgvs_c)
                except vvhgvs.exceptions.HGVSError:
                    pass

                old_gen_ref = hgvs_c.posedit.edit.ref

            # If normalisation could not provide the reference, fetch it
            # directly from the transcript.
            if not old_gen_ref:
                if hgvs_c.type == 'c':
                    fix_n = self.vm.c_to_n(hgvs_c)
                else:
                    fix_n = hgvs_c

                old_gen_ref = self.sf.fetch_seq(
                    fix_n.ac,
                    fix_n.posedit.pos.start.base - 1,
                    fix_n.posedit.pos.end.base
                )

            strand = 1

            for option in mapping_options:
                if option[1] == hgvs_genomic.ac:
                    strand = int(option[4])
                    break

            if strand == -1:
                old_gen_ref = self.revcomp(old_gen_ref)

            if old_gen_ref != hgvs_genomic.posedit.edit.ref:
                if hgvs_c.posedit.edit.type == 'dup':
                    alt = old_gen_ref + old_gen_ref
                else:
                    alt = self.revcomp(old_gen_ref)

                hgvs_genomic = hgvs_delins_parts_to_hgvs_obj(
                    hgvs_genomic.ac,
                    hgvs_genomic.type,
                    hgvs_genomic.posedit.pos,
                    hgvs_genomic.posedit.edit.ref,
                    alt,
                    offset_pos=True
                )

                try:
                    hn.normalize(hgvs_genomic)
                except vvhgvs.exceptions.HGVSError:
                    pass

        return hgvs_genomic

    def noreplace_myevm_t_to_g(self, hgvs_c, variant):
        """
        USE WITH MAPPER THAT DOES NOT REPLACE THE REFERENCE GENOMIC BASES
        AND DOES NOT NORMALIZE.

        Enhanced transcript-to-genome mapping function using evm.

        Attempts to return a genomic mapping even when the transcript does not
        map directly to the specified primary assembly by trying the available
        UTA mapping options in order.

        Returns a parsed HGVS g. object.
        """
        alt_aln_method = self.alt_aln_method
        hgvs_genomic = None
        attempted_mapping_errors = []

        try:
            hgvs_genomic = variant.evm.t_to_g(hgvs_c)
            variant.hn.normalize(hgvs_genomic)

        # This can fail when multiple genomic references are available.
        except vvhgvs.exceptions.HGVSError:
            mapping_options = variant.map_dat.mapping_options(
                hgvs_c.ac,
                hdp=self.hdp
            )

            if not mapping_options:
                raise HGVSDataNotAvailableError(
                    "no g. mapping options available"
                )

            def search_in_options(
                    seqtype,
                    chr_num_val,
                    final=False
            ):
                for op in mapping_options:
                    if op[2].startswith('blat'):
                        continue

                    genomic_ac = op[1]

                    if not genomic_ac.startswith(seqtype):
                        continue

                    if not final:
                        chr_num = seq_data.is_supported_for_mapping(
                            genomic_ac,
                            variant.primary_assembly
                        )

                        if chr_num_val:
                            if chr_num == 'false':
                                continue
                        elif chr_num != 'false':
                            continue

                    try:
                        return self.vm.t_to_g(
                            hgvs_c,
                            genomic_ac,
                            alt_aln_method
                        )

                    except Exception as e:
                        attempted_mapping_errors.append(
                            f"{e}/{hgvs_c.ac}/{genomic_ac}~"
                        )

                return None

            # Preserve the historical mapping preference order.
            mapping_searches = (
                ('NC_', True, False),
                ('NC_', True, False),
                ('NC_', False, False),
                ('NT_', True, False),
                ('NT_', False, False),
                ('NW_', True, False),
                ('NW_', False, False),
                ('NG_', True, True),
            )

            for seqtype, chr_num_val, final in mapping_searches:
                candidate = search_in_options(
                    seqtype,
                    chr_num_val,
                    final=final
                )

                if candidate is None:
                    continue

                hgvs_genomic = candidate

                # NG_ is the final RefSeqGene fallback. Historically it was
                # returned without another normalization gate.
                if final:
                    break

                try:
                    variant.hn.normalize(hgvs_genomic)
                except vvhgvs.exceptions.HGVSError:
                    continue

                break

        if hgvs_genomic is None:
            raise HGVSDataNotAvailableError(
                'No available t_to_g liftover'
            )

        # Insertions at exon/exon boundaries can map badly. Convert these to
        # an equivalent delins when normalization identifies the problem.
        if (
                hgvs_c.posedit.edit.type == 'ins'
                and not hgvs_position_utils.either_position_is_intronic(hgvs_c)
        ):
            try:
                variant.hn.normalize(hgvs_genomic)

            except vvhgvs.exceptions.HGVSError as e:
                if str(e) == 'insertion length must be 1':
                    if hgvs_c.type == 'c':
                        hgvs_t = self.vm.c_to_n(hgvs_c)
                    else:
                        hgvs_t = copy.copy(hgvs_c)

                    ins_ref = self.sf.fetch_seq(
                        hgvs_t.ac,
                        hgvs_t.posedit.pos.start.base - 1,
                        hgvs_t.posedit.pos.end.base
                    )

                    ins_alt = (
                            ins_ref[:1]
                            + hgvs_t.posedit.edit.alt
                            + ins_ref[-1:]
                    )

                    hgvs_t = hgvs_delins_parts_to_hgvs_obj(
                        hgvs_t.ac,
                        hgvs_t.type,
                        hgvs_t.posedit.pos.start.base,
                        ins_ref,
                        ins_alt,
                        end=hgvs_t.posedit.pos.end.base,
                        offset_pos=True
                    )

                    try:
                        hgvs_c = self.vm.n_to_c(hgvs_t)
                    except vvhgvs.exceptions.HGVSError:
                        hgvs_c = copy.copy(hgvs_t)

                    try:
                        hgvs_genomic = variant.no_norm_evm.t_to_g(
                            hgvs_c
                        )
                    except vvhgvs.exceptions.HGVSError as e:
                        logger.info(
                            "Ins mapping error in noreplace_myevm_t_to_g %s",
                            e
                        )

        return hgvs_genomic

    def myevm_g_to_t(self, evm, hgvs_genomic, alt_ac):
        """
        Enhanced transcript to genome position on a specified genomic reference using vm
        Deals with mapping from transcript positions that do not exist in the genomic sequence
        i.e. the stated position aligns to a genomic gap!
        returns parsed hgvs g. object
        """
        hgvs_t = evm.g_to_t(hgvs_genomic, alt_ac)
        return hgvs_t

    def myvm_t_to_g(self, hgvs_c, alt_chr, no_norm_evm, hn, map_dat):
        # Store the input
        alt_aln_method = self.alt_aln_method
        stored_hgvs_c = copy.deepcopy(hgvs_c)
        expand_out = False

        utilise_gap_code = map_dat.is_gapped_map(
            hgvs_c.ac,
            alt_chr,
            hdp=self.hdp
        )

        # Warn gap code in use
        logger.debug("gap_compensation_mvm = %s", utilise_gap_code)

        if (
                utilise_gap_code
                and hgvs_c.posedit.edit.type
                in ('identity', 'del', 'delins', 'dup', 'sub', 'ins', 'inv')
        ):
            # If NM_ need the n. position
            if hgvs_c.type == 'c':
                hgvs_c = no_norm_evm.c_to_n(hgvs_c)

            # Check for intronic
            try:
                hn.normalize(hgvs_c)
            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)

                if 'intronic variant' in error:
                    logger.debug("Except passed, %s", e)

                elif (
                        'Length implied by coordinates must equal sequence deletion length'
                        in error
                        and hgvs_c.ac.startswith('NR_')
                ):
                    hgvs_c.posedit.pos.end.base = (
                            hgvs_c.posedit.pos.start.base
                            + len(hgvs_c.posedit.edit.ref)
                            - 1
                    )

            # Check again before continuing
            if not hgvs_position_utils.either_position_is_intronic(hgvs_c):
                try:
                    # For non-intronic sequence
                    hgvs_t = copy.deepcopy(hgvs_c)

                    # Handle inversions
                    if hgvs_t.posedit.edit.type == 'inv':
                        inv_alt = self.revcomp(hgvs_t.posedit.edit.ref)

                        pre_base, post_base = self._expand_ref(
                            hgvs_t.ac,
                            hgvs_t.posedit.pos.start.base - 2,
                            hgvs_t.posedit.pos.end.base + 1
                        )

                        hgvs_t = hgvs_delins_parts_to_hgvs_obj(
                            hgvs_t.ac,
                            hgvs_t.type,
                            hgvs_t.posedit.pos.start.base - 1,
                            pre_base + hgvs_t.posedit.edit.ref + post_base,
                            pre_base + inv_alt + post_base,
                            end=hgvs_t.posedit.pos.end.base + 1,
                            offset_pos=True
                        )

                    if hgvs_c.posedit.edit.type == 'dup':
                        pre_base, post_base = self._expand_ref(
                            hgvs_t.ac,
                            hgvs_t.posedit.pos.start.base - 2,
                            hgvs_t.posedit.pos.end.base + 1
                        )

                        hgvs_t = hgvs_delins_parts_to_hgvs_obj(
                            hgvs_t.ac,
                            hgvs_t.type,
                            hgvs_t.posedit.pos.start.base - 1,
                            pre_base + hgvs_t.posedit.edit.ref + post_base,
                            (
                                    pre_base
                                    + hgvs_t.posedit.edit.ref
                                    + hgvs_t.posedit.edit.ref
                                    + post_base
                            ),
                            offset_pos=True
                        )

                    elif hgvs_c.posedit.edit.type == 'ins':
                        ins_ref = self.sf.fetch_seq(
                            hgvs_t.ac,
                            hgvs_t.posedit.pos.start.base - 2,
                            hgvs_t.posedit.pos.end.base + 1
                        )

                        ins_alt = (
                                ins_ref[:2]
                                + hgvs_t.posedit.edit.alt
                                + ins_ref[-2:]
                        )

                        hgvs_t = hgvs_delins_parts_to_hgvs_obj(
                            hgvs_t.ac,
                            hgvs_t.type,
                            hgvs_t.posedit.pos.start.base - 1,
                            ins_ref,
                            ins_alt,
                            end=hgvs_t.posedit.pos.end.base + 1,
                            offset_pos=True
                        )

                    else:
                        if hgvs_t.posedit.edit.alt is None:
                            hgvs_t.posedit.edit.alt = ''

                        pre_base, post_base = self._expand_ref(
                            hgvs_t.ac,
                            hgvs_t.posedit.pos.start.base - 2,
                            hgvs_t.posedit.pos.end.base + 1
                        )

                        hgvs_t = hgvs_delins_parts_to_hgvs_obj(
                            hgvs_t.ac,
                            hgvs_t.type,
                            hgvs_t.posedit.pos.start.base - 1,
                            pre_base + hgvs_t.posedit.edit.ref + post_base,
                            pre_base + hgvs_t.posedit.edit.alt + post_base,
                            end=hgvs_t.posedit.pos.end.base + 1,
                            offset_pos=True
                        )

                    hgvs_c = copy.deepcopy(hgvs_t)

                    # Set expanded out test to true
                    expand_out = True

                except Exception:
                    pass

            # Convert back to c. position from n. position
            try:
                hgvs_c = no_norm_evm.n_to_c(hgvs_c)
            except vvhgvs.exceptions.HGVSError:
                hgvs_c = copy.deepcopy(stored_hgvs_c)

            # Ensure the altered c. variant has not crossed intron/exon boundaries
            hgvs_check_boundaries = copy.deepcopy(hgvs_c)

            try:
                hn.normalize(hgvs_check_boundaries)
            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)

                if 'spanning the exon-intron boundary' in error:
                    hgvs_c = copy.deepcopy(stored_hgvs_c)

            # Catch identity at the exon/intron boundary by trying to normalize ref only
            if hgvs_check_boundaries.posedit.edit.type == 'identity':
                hgvs_reform_ident = hgvs_delins_parts_to_hgvs_obj(
                    hgvs_c.ac,
                    stored_hgvs_c.type,
                    hgvs_c.posedit.pos,
                    hgvs_c.posedit.edit.ref,
                    '',
                    offset_pos=True
                )

                try:
                    hn.normalize(hgvs_reform_ident)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)

                    if (
                            'spanning the exon-intron boundary' in error
                            or 'Normalization of intronic variants' in error
                    ):
                        hgvs_c = copy.deepcopy(stored_hgvs_c)

        hgvs_genomic = self.vm.t_to_g(
            hgvs_c,
            alt_chr,
            alt_aln_method
        )

        if (
                hgvs_c.posedit.edit.type == 'identity'
                and hgvs_genomic.posedit.edit.type == 'delins'
                and hgvs_genomic.posedit.edit.alt == ''
                and not expand_out
        ):
            hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.ref

        if hgvs_genomic.posedit.edit.type == 'ins' and utilise_gap_code:

            if stored_hgvs_c.posedit.edit.type == "dup":
                stored_hgvs_c = hgvs_dup_to_delins(stored_hgvs_c)

            try:
                # Can move ins variants (and in doing so break
                # mid base == original bases assumption)
                pre_norm_genomic = copy.copy(hgvs_genomic)
                hgvs_genomic = hn.normalize(hgvs_genomic)

                if (
                        stored_hgvs_c.posedit.edit.alt
                        and len(stored_hgvs_c.posedit.edit.alt) + 2
                        == len(hgvs_c.posedit.edit.alt)
                        and hgvs_c.posedit.edit.alt
                        == pre_norm_genomic.posedit.edit.alt
                ):
                    pre_norm_genomic.posedit.edit.alt = (
                        pre_norm_genomic.posedit.edit.alt[1:-1]
                    )
                    hgvs_genomic = copy.copy(pre_norm_genomic)
                    hgvs_genomic = hn.normalize(hgvs_genomic)
                    hgvs_c.posedit.edit.alt = hgvs_c.posedit.edit.alt[1:-1]

            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)

                if error == 'insertion length must be 1':
                    ref = self.sf.fetch_seq(
                        hgvs_genomic.ac,
                        hgvs_genomic.posedit.pos.start.base - 1,
                        hgvs_genomic.posedit.pos.end.base
                    )
                    hgvs_genomic.posedit.edit.ref = ref
                    hgvs_genomic.posedit.edit.alt = (
                            ref[:1]
                            + hgvs_genomic.posedit.edit.alt
                            + ref[-1:]
                    )
                    hgvs_genomic = hn.normalize(hgvs_genomic)

                if error == 'base start position must be <= end position':
                    start = hgvs_genomic.posedit.pos.start.base
                    end = hgvs_genomic.posedit.pos.end.base
                    hgvs_genomic.posedit.pos.start.base = end
                    hgvs_genomic.posedit.pos.end.base = start
                    hgvs_genomic = hn.normalize(hgvs_genomic)

            except AttributeError as e:
                if "'Dup' object has no attribute 'alt'" in str(e):
                    logger.exception(
                        "Code triggered previously in very poor alignment so not "
                        "able to fully test, refer to test_inputs.py tests "
                        "test_alt_gapping_bug: hgvs_genomic: %s, "
                        "stored_hgvs_c: %s",
                        hgvs_genomic,
                        stored_hgvs_c
                    )
                    raise

        # Statements required to reformat stored_hgvs_c into a usable synonym
        if (
                (
                        stored_hgvs_c.posedit.edit.ref == ''
                        or stored_hgvs_c.posedit.edit.ref is None
                )
                and expand_out
        ):
            if stored_hgvs_c.type == 'c':
                stored_hgvs_n = self.vm.c_to_n(stored_hgvs_c)
            else:
                stored_hgvs_n = stored_hgvs_c

            stored_ref = self.sf.fetch_seq(
                stored_hgvs_n.ac,
                stored_hgvs_n.posedit.pos.start.base - 1,
                stored_hgvs_n.posedit.pos.end.base
            )
            stored_hgvs_c.posedit.edit.ref = stored_ref

        # First look for variants mapping to the flanks of gaps
        # either in the gap or on the flank but not fully within the gap
        if expand_out:
            nr_genomic = self.nr_vm.t_to_g(
                hgvs_c,
                hgvs_genomic.ac,
                alt_aln_method
            )

            try:
                hn.normalize(nr_genomic)

            except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                error_type_1 = str(e)

                if (
                        'Length implied by coordinates must equal sequence deletion length'
                        in error_type_1
                        or error_type_1
                        == 'base start position must be <= end position'
                ):
                    # This code is designed to use the fact that no-replace
                    # mappings don't adjust length to match the span to detect,
                    # and if needed handle, variants directly proximal to gap
                    # boundaries.
                    genomic_gap_variant = None

                    if (
                            'Length implied by coordinates must equal sequence deletion length'
                            in error_type_1
                    ):
                        logger.info(
                            'Variant is proximal to the flank of a genomic gap'
                        )

                        genomic_gap_variant = self.vm.t_to_g(
                            stored_hgvs_c,
                            hgvs_genomic.ac,
                            alt_aln_method
                        )

                        try:
                            hn.normalize(genomic_gap_variant)

                        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                            if (
                                    'base start position must be <= end position'
                                    in str(e)
                                    and 'Length implied by coordinates must equal'
                                    in error_type_1
                            ):
                                make_gen_var = copy.copy(nr_genomic)
                                make_gen_var.posedit.edit.ref = self.sf.fetch_seq(
                                    nr_genomic.ac,
                                    nr_genomic.posedit.pos.start.base - 1,
                                    nr_genomic.posedit.pos.end.base
                                )
                                genomic_gap_variant = make_gen_var
                                error_type_1 = None

                        else:
                            if (
                                    genomic_gap_variant.posedit.edit.ref is None
                                    and 'Length implied by coordinates must equal'
                                    in error_type_1
                            ):
                                # Handle gaps correctly for current expectations.
                                # With delGCTinsGGT when tx GCT maps to gen GT,
                                # output should be delG, not C>G, when normalized.
                                make_gen_var = copy.copy(nr_genomic)
                                make_gen_var.posedit.edit.ref = self.sf.fetch_seq(
                                    nr_genomic.ac,
                                    nr_genomic.posedit.pos.start.base - 1,
                                    nr_genomic.posedit.pos.end.base
                                )
                                genomic_gap_variant = make_gen_var
                                error_type_1 = None

                            else:
                                genomic_gap_variant = self.nr_vm.t_to_g(
                                    hgvs_c,
                                    hgvs_genomic.ac,
                                    alt_aln_method=alt_aln_method
                                )

                    if error_type_1 == 'base start position must be <= end position':
                        logger.info('Variant is fully within a genomic gap')
                        genomic_gap_variant = self.vm.t_to_g(
                            stored_hgvs_c,
                            hgvs_genomic.ac,
                            alt_aln_method
                        )

                    # Logic:
                    # We have checked that the variant does not cross boundaries,
                    # or is intronic, so it is likely mapping to a genomic gap.
                    try:
                        hn.normalize(genomic_gap_variant)

                    except Exception as ea1:
                        if str(ea1) == 'base start position must be <= end position':
                            # This will only happen when the variant is fully
                            # within the gap.
                            gap_start = genomic_gap_variant.posedit.pos.end.base
                            gap_end = genomic_gap_variant.posedit.pos.start.base
                            genomic_gap_variant.posedit.pos.start.base = gap_start
                            genomic_gap_variant.posedit.pos.end.base = gap_end

                        if (
                                'Length implied by coordinates must equal sequence deletion length'
                                in str(ea1)
                        ):
                            # Variant flanks the gap but is not inside it.
                            try:
                                try:
                                    norm_stored_c = hn.normalize(stored_hgvs_c)
                                except HGVSUnsupportedOperationError:
                                    norm_stored_c = stored_hgvs_c

                                if norm_stored_c.posedit.edit.type in (
                                        'sub',
                                        'identity'
                                ):
                                    flank_hgvs_genomic = self.vm.t_to_g(
                                        stored_hgvs_c,
                                        genomic_gap_variant.ac,
                                        alt_aln_method
                                    )
                                    init_flank_hgvs_genomic = copy.copy(
                                        flank_hgvs_genomic
                                    )

                                    # Handle genomic opening gap
                                    if (
                                            len(flank_hgvs_genomic.posedit.edit.ref)
                                            < len(stored_hgvs_c.posedit.edit.ref)
                                            and len(stored_hgvs_c.posedit.edit.ref)
                                            == len(
                                        genomic_gap_variant.posedit.edit.ref
                                    ) - 2
                                            and len(
                                        genomic_gap_variant.posedit.edit.ref
                                    )
                                            == len(
                                        genomic_gap_variant.posedit.edit.alt
                                    )
                                    ):
                                        n_flank_hgvs_genomic = hn.normalize(
                                            init_flank_hgvs_genomic
                                        )

                                        if (
                                                n_flank_hgvs_genomic
                                                == init_flank_hgvs_genomic
                                        ):
                                            return self.vm.t_to_g(
                                                norm_stored_c,
                                                genomic_gap_variant.ac,
                                                alt_aln_method
                                            )

                                        return hn.normalize(
                                            init_flank_hgvs_genomic
                                        )

                                    flank_hgvs_genomic = self.vm.t_to_g(
                                        norm_stored_c,
                                        genomic_gap_variant.ac,
                                        alt_aln_method
                                    )

                                    self.vr.validate(flank_hgvs_genomic)

                                    # Gap in the transcript e.g. NR2E3 tests
                                    if (
                                            len(
                                                init_flank_hgvs_genomic.posedit.edit.ref
                                            )
                                            > len(stored_hgvs_c.posedit.edit.ref)
                                            and len(
                                        stored_hgvs_c.posedit.edit.ref
                                    )
                                            == len(
                                        genomic_gap_variant.posedit.edit.ref
                                    ) - 2
                                    ):
                                        return hn.normalize(
                                            init_flank_hgvs_genomic
                                        )

                                    elif (
                                            flank_hgvs_genomic.posedit.edit.type
                                            == 'sub'
                                            and norm_stored_c.posedit.edit.type
                                            == 'sub'
                                            and stored_hgvs_c.posedit.edit.type
                                            == 'sub'
                                            and (
                                                    len(
                                                        genomic_gap_variant.posedit.edit.ref
                                                    ) + 1
                                                    == (
                                                            genomic_gap_variant.posedit.pos.end.base
                                                            - genomic_gap_variant.posedit.pos.start.base
                                                            + 1
                                                    )
                                            )
                                    ):
                                        pass

                                    else:
                                        return flank_hgvs_genomic

                            # Will occur if the variant still overlaps/is in gap
                            except HGVSInvalidVariantError:
                                pass

                            # If test fails, continue old processing
                            gap_start = (
                                    genomic_gap_variant.posedit.pos.start.base - 1
                            )
                            gap_end = (
                                    genomic_gap_variant.posedit.pos.end.base + 1
                            )
                            genomic_gap_variant.posedit.pos.start.base = gap_start
                            genomic_gap_variant.posedit.pos.end.base = gap_end
                            genomic_gap_variant.posedit.edit.ref = ''
                            stored_hgvs_c = copy.deepcopy(hgvs_c)

                        # Remove alt
                        try:
                            genomic_gap_variant.posedit.edit.alt = ''
                        except Exception as e:
                            logger.debug("Except passed, %s", e)

                        # Should be a delins so will normalize statically and
                        # replace the reference bases
                        genomic_gap_variant = hn.normalize(genomic_gap_variant)

                        # Static map to c. and static normalize
                        transcript_gap_variant = self.vm.g_to_t(
                            genomic_gap_variant,
                            hgvs_c.ac,
                            alt_aln_method=alt_aln_method
                        )

                        if (
                                'Length implied by coordinates must equal sequence deletion length'
                                not in str(ea1)
                        ):
                            try:
                                transcript_gap_variant = hn.normalize(
                                    transcript_gap_variant
                                )
                            except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                                logger.debug("Except passed, %s", e)

                        # If NM_ need the n. position
                        try:
                            transcript_gap_n = no_norm_evm.c_to_n(
                                transcript_gap_variant
                            )
                            transcript_gap_alt_n = no_norm_evm.c_to_n(
                                stored_hgvs_c
                            )
                        except vvhgvs.exceptions.HGVSError:
                            transcript_gap_n = transcript_gap_variant
                            transcript_gap_alt_n = stored_hgvs_c

                        # Ensure an ALT exists
                        try:
                            if transcript_gap_alt_n.posedit.edit.alt is None:
                                transcript_gap_alt_n.posedit.edit.alt = 'X'
                        except Exception as e:
                            if str(e) == "'Dup' object has no attribute 'alt'":
                                transcript_gap_n = hgvs_dup_to_delins(
                                    transcript_gap_n
                                )
                                transcript_gap_alt_n = hgvs_dup_to_delins(
                                    transcript_gap_alt_n
                                )

                        # Split reference/replacing ALT sequence into dictionaries
                        reference_bases = list(
                            transcript_gap_n.posedit.edit.ref
                        )

                        if transcript_gap_alt_n.posedit.edit.alt is not None:
                            alternate_bases = list(
                                transcript_gap_alt_n.posedit.edit.alt
                            )
                        else:
                            alternate_bases = ['X'] * len(
                                transcript_gap_alt_n.posedit.edit.ref
                            )

                        ref_start = transcript_gap_n.posedit.pos.start.base
                        alt_start = transcript_gap_alt_n.posedit.pos.start.base

                        ref_base_dict = {}

                        for base in reference_bases:
                            ref_base_dict[ref_start] = base
                            ref_start += 1

                        alt_base_dict = {}

                        # All variants forced into delete-insert format.
                        # Deleted ALT bases are represented by X.
                        for i in range(
                                transcript_gap_alt_n.posedit.pos.start.base,
                                transcript_gap_alt_n.posedit.pos.end.base + 1
                        ):
                            if i == alt_start:
                                alt_base_dict[i] = ''.join(alternate_bases)
                            else:
                                alt_base_dict[i] = 'X'

                        alternate_sequence_bases = []

                        for i in range(
                                transcript_gap_n.posedit.pos.start.base,
                                transcript_gap_n.posedit.pos.end.base + 1
                        ):
                            if i in alt_base_dict:
                                alternate_sequence_bases.append(
                                    alt_base_dict[i]
                                )
                            elif i in ref_base_dict:
                                alternate_sequence_bases.append(
                                    ref_base_dict[i]
                                )

                        alternate_sequence = ''.join(
                            alternate_sequence_bases
                        ).replace('X', '')

                        transcript_gap_n.posedit.edit.alt = alternate_sequence

                        try:
                            transcript_gap_variant = self.vm.n_to_c(
                                transcript_gap_n
                            )
                        except vvhgvs.exceptions.HGVSError:
                            transcript_gap_variant = transcript_gap_n

                        try:
                            hgvs_genomic = self.vm.t_to_g(
                                transcript_gap_variant,
                                hgvs_genomic.ac,
                                alt_aln_method
                            )
                            pre_norm_genomic = copy.copy(hgvs_genomic)
                            hgvs_genomic = hn.normalize(hgvs_genomic)

                        except Exception as e:
                            if str(e) == "base start position must be <= end position":
                                # Expansion out required to map back to genome
                                pre_base, post_base = self._expand_ref(
                                    transcript_gap_n.ac,
                                    transcript_gap_n.posedit.pos.start.base - 2,
                                    transcript_gap_n.posedit.pos.end.base + 1
                                )

                                transcript_gap_n.posedit.pos.start.base -= 1
                                transcript_gap_n.posedit.pos.end.base += 1

                                transcript_gap_n.posedit.edit.ref = (
                                        pre_base
                                        + transcript_gap_n.posedit.edit.ref
                                        + post_base
                                )
                                transcript_gap_n.posedit.edit.alt = (
                                        pre_base
                                        + transcript_gap_n.posedit.edit.alt
                                        + post_base
                                )

                                try:
                                    transcript_gap_variant = self.vm.n_to_c(
                                        transcript_gap_n
                                    )
                                except vvhgvs.exceptions.HGVSError:
                                    transcript_gap_variant = transcript_gap_n

                                hgvs_genomic = self.vm.t_to_g(
                                    transcript_gap_variant,
                                    hgvs_genomic.ac,
                                    alt_aln_method
                                )
                                pre_norm_genomic = copy.copy(hgvs_genomic)
                                hgvs_genomic = hn.normalize(hgvs_genomic)

                        # Bypass the next bit of gap code
                        expand_out = False

        # CASCADING STATEMENTS WHICH CAPTURE t to g MAPPING OPTIONS
        # Remove identity bases
        if hgvs_c == stored_hgvs_c:
            expand_out = False

        elif not expand_out or not utilise_gap_code:
            pass

        # Correct ref inside gap
        elif expand_out and hgvs_genomic.posedit.edit.ref is None:
            # Inserted entirely inside a gap in the genomic sequence.
            if (
                    hgvs_genomic.posedit.edit.alt is not None
                    and len(hgvs_genomic.posedit.edit.alt) > 2
            ):
                hgvs_genomic = pre_norm_genomic
                hgvs_genomic.posedit.edit.alt = (
                    hgvs_genomic.posedit.edit.alt[1:-1]
                )

                try:
                    hgvs_genomic = hn.normalize(hgvs_genomic)
                except vvhgvs.exceptions.HGVSError:
                    pass

        # Correct expansion ref + 2
        elif (
                expand_out
                and len(hgvs_genomic.posedit.edit.ref)
                == len(stored_hgvs_c.posedit.edit.ref) + 2
        ):
            hgvs_genomic.posedit.pos.start.base += 1
            hgvs_genomic.posedit.pos.end.base -= 1
            hgvs_genomic.posedit.edit.ref = hgvs_genomic.posedit.edit.ref[1:-1]

            try:
                if hgvs_genomic.posedit.edit.alt is not None:
                    hgvs_genomic.posedit.edit.alt = (
                        hgvs_genomic.posedit.edit.alt[1:-1]
                    )
            except AttributeError:
                pass

        elif (
                expand_out
                and len(hgvs_genomic.posedit.edit.ref)
                != len(stored_hgvs_c.posedit.edit.ref) + 2
        ):
            if len(hgvs_genomic.posedit.edit.ref) == 2:
                hn.normalize(hgvs_genomic)

            # Likely if start/end aligns to a gap in genomic sequence
            elif len(hgvs_genomic.posedit.edit.ref) <= 1:
                genomic_gap_variant = self.vm.t_to_g(
                    stored_hgvs_c,
                    hgvs_genomic.ac,
                    alt_aln_method
                )

                try:
                    hn.normalize(genomic_gap_variant)

                except Exception as e:
                    if str(e) == 'base start position must be <= end position':
                        gap_start = genomic_gap_variant.posedit.pos.end.base
                        gap_end = genomic_gap_variant.posedit.pos.start.base
                        genomic_gap_variant.posedit.pos.start.base = gap_start
                        genomic_gap_variant.posedit.pos.end.base = gap_end

                    # Remove alt
                    try:
                        genomic_gap_variant.posedit.edit.alt = ''
                    except Exception as e:
                        logger.debug("Except passed, %s", e)

                    genomic_gap_variant = hn.normalize(genomic_gap_variant)

                    transcript_gap_variant = self.vm.g_to_t(
                        genomic_gap_variant,
                        hgvs_c.ac,
                        alt_aln_method=alt_aln_method
                    )
                    transcript_gap_variant = hn.normalize(
                        transcript_gap_variant
                    )

                    # If NM_ need the n. position
                    try:
                        transcript_gap_n = no_norm_evm.c_to_n(
                            transcript_gap_variant
                        )
                        transcript_gap_alt_n = no_norm_evm.c_to_n(
                            stored_hgvs_c
                        )
                    except vvhgvs.exceptions.HGVSError:
                        transcript_gap_n = transcript_gap_variant
                        transcript_gap_alt_n = stored_hgvs_c

                    # Ensure an ALT exists
                    try:
                        if transcript_gap_alt_n.posedit.edit.alt is None:
                            transcript_gap_alt_n.posedit.edit.alt = 'X'
                    except Exception as e:
                        if str(e) == "'Dup' object has no attribute 'alt'":
                            transcript_gap_n = hgvs_dup_to_delins(
                                transcript_gap_n
                            )
                            transcript_gap_alt_n = hgvs_dup_to_delins(
                                transcript_gap_alt_n
                            )

                    reference_bases = list(
                        transcript_gap_n.posedit.edit.ref
                    )

                    if transcript_gap_alt_n.posedit.edit.alt is not None:
                        alternate_bases = list(
                            transcript_gap_alt_n.posedit.edit.alt
                        )
                    else:
                        alternate_bases = ['X'] * len(
                            transcript_gap_alt_n.posedit.edit.ref
                        )

                    ref_start = transcript_gap_n.posedit.pos.start.base
                    alt_start = transcript_gap_alt_n.posedit.pos.start.base

                    ref_base_dict = {}

                    for base in reference_bases:
                        ref_base_dict[ref_start] = base
                        ref_start += 1

                    alt_base_dict = {}

                    for i in range(
                            transcript_gap_alt_n.posedit.pos.start.base,
                            transcript_gap_alt_n.posedit.pos.end.base + 1
                    ):
                        if i == alt_start:
                            alt_base_dict[i] = ''.join(alternate_bases)
                        else:
                            alt_base_dict[i] = 'X'

                    alternate_sequence_bases = []

                    for i in range(
                            transcript_gap_n.posedit.pos.start.base,
                            transcript_gap_n.posedit.pos.end.base + 1
                    ):
                        if i in alt_base_dict:
                            alternate_sequence_bases.append(
                                alt_base_dict[i]
                            )
                        elif i in ref_base_dict:
                            alternate_sequence_bases.append(
                                ref_base_dict[i]
                            )

                    alternate_sequence = ''.join(
                        alternate_sequence_bases
                    ).replace('X', '')

                    transcript_gap_n.posedit.edit.alt = alternate_sequence

                    try:
                        transcript_gap_variant = self.vm.n_to_c(
                            transcript_gap_n
                        )
                    except vvhgvs.exceptions.HGVSError:
                        transcript_gap_variant = transcript_gap_n

                    try:
                        hgvs_genomic = self.vm.t_to_g(
                            transcript_gap_variant,
                            hgvs_genomic.ac,
                            alt_aln_method
                        )
                        hgvs_genomic = hn.normalize(hgvs_genomic)

                    except Exception as e:
                        if str(e) == "base start position must be <= end position":
                            pre_base, post_base = self._expand_ref(
                                transcript_gap_n.ac,
                                transcript_gap_n.posedit.pos.start.base - 2,
                                transcript_gap_n.posedit.pos.end.base + 1
                            )

                            transcript_gap_n.posedit.pos.start.base -= 1
                            transcript_gap_n.posedit.pos.end.base += 1

                            transcript_gap_n.posedit.edit.ref = (
                                    pre_base
                                    + transcript_gap_n.posedit.edit.ref
                                    + post_base
                            )
                            transcript_gap_n.posedit.edit.alt = (
                                    pre_base
                                    + transcript_gap_n.posedit.edit.alt
                                    + post_base
                            )

                            try:
                                transcript_gap_variant = self.vm.n_to_c(
                                    transcript_gap_n
                                )
                            except vvhgvs.exceptions.HGVSError:
                                transcript_gap_variant = transcript_gap_n

                            hgvs_genomic = self.vm.t_to_g(
                                transcript_gap_variant,
                                hgvs_genomic.ac,
                                alt_aln_method
                            )
                            hgvs_genomic = hn.normalize(hgvs_genomic)

        # Ins variants map badly - especially between c. exon/exon boundary
        if (
                hgvs_c.posedit.edit.type == 'ins'
                and not hgvs_position_utils.either_position_is_intronic(hgvs_c)
        ):
            try:
                hn.normalize(hgvs_genomic)

            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)

                if error == 'insertion length must be 1':
                    if hgvs_c.type == 'c':
                        hgvs_t = self.vm.c_to_n(hgvs_c)
                    else:
                        hgvs_t = copy.copy(hgvs_c)

                    ins_ref = self.sf.fetch_seq(
                        hgvs_t.ac,
                        hgvs_t.posedit.pos.start.base - 1,
                        hgvs_t.posedit.pos.end.base
                    )

                    ins_alt = (
                            ins_ref[:1]
                            + hgvs_t.posedit.edit.alt
                            + ins_ref[-1:]
                    )

                    hgvs_t = hgvs_delins_parts_to_hgvs_obj(
                        hgvs_t.ac,
                        hgvs_t.type,
                        hgvs_t.posedit.pos.start.base,
                        ins_ref,
                        ins_alt,
                        end=hgvs_t.posedit.pos.end.base,
                        offset_pos=True
                    )

                    try:
                        hgvs_c = self.vm.n_to_c(hgvs_t)
                    except vvhgvs.exceptions.HGVSError:
                        hgvs_c = copy.copy(hgvs_t)

                    try:
                        hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)
                    except vvhgvs.exceptions.HGVSError as e:
                        logger.info(
                            'Ins mapping error in myt_to_g %s',
                            e
                        )

        return hgvs_genomic

    def hgvs_r_to_c(self, hgvs_object):
        """
        Convert r. into c.
        """
        # Check for LRG_t with r.
        if 'LRG' in hgvs_object.ac:
            transcript_ac = self.db.get_refseq_transcript_id_from_lrg_transcript_id(
                hgvs_object.ac
            )
            if transcript_ac == 'none':
                raise HGVSDataNotAvailableError(
                    'Unable to identify a relevant transcript for '
                    + hgvs_object.ac
                )
            hgvs_object.ac = transcript_ac

        hgvs_object.type = 'c'
        edit = hgvs_object.posedit.edit

        # Uppercase and switch U to T.
        try:
            edit.ref = edit.ref.upper().replace('U', 'T')
        except AttributeError:
            pass

        try:
            edit.alt = edit.alt.upper().replace('U', 'T')
        except AttributeError:
            pass

        # Map positions onto coding-coordinate datums.
        if hgvs_position_utils.start_is_3_prime_utr(hgvs_object):
            hgvs_position_utils.set_start_as_3_prime_utr(hgvs_object)
        else:
            hgvs_object.posedit.pos.start.datum = Datum.CDS_START

        if hgvs_position_utils.end_is_3_prime_utr(hgvs_object):
            hgvs_position_utils.set_end_as_3_prime_utr(hgvs_object)
        else:
            hgvs_object.posedit.pos.end.datum = Datum.CDS_START

        return hgvs_object

    def relevant_transcripts(
            self,
            hgvs_genomic,
            evm,
            alt_aln_method,
            reverse_normalizer,
            select_transcripts
    ):
        """
        Automatically maps genomic positions onto all overlapping transcripts.
        """
        # The two region queries differ by one base at each boundary.
        # Combine both to avoid missing transcripts at either end.
        rts_list = self.hdp.get_tx_for_region(
            hgvs_genomic.ac,
            alt_aln_method,
            hgvs_genomic.posedit.pos.start.base - 1,
            hgvs_genomic.posedit.pos.end.base - 1
        )

        rts_dict = {
            tx_dat['tx_ac']: tx_dat['alt_strand']
            for tx_dat in rts_list
        }

        rts_list_2 = self.hdp.get_tx_for_region(
            hgvs_genomic.ac,
            self.alt_aln_method,
            hgvs_genomic.posedit.pos.start.base,
            hgvs_genomic.posedit.pos.end.base
        )

        for tx_dat in rts_list_2:
            rts_dict[tx_dat['tx_ac']] = tx_dat['alt_strand']

        rts = list(rts_dict)

        # Filter out transcripts that are not the latest versions.
        if select_transcripts in ("all", "None", None):
            rts = self.transcript_filter(rts)

        elif (
                select_transcripts not in ("all", "None", "raw")
                and "mane" not in select_transcripts
                and "select" not in select_transcripts
        ):
            rts = self.transcript_filter(rts, select_transcripts)

        # Prepare insertion as a forced delins for mappings where HGVS insertion
        # handling otherwise fails.
        hgvs_genomic_forced_delins = None

        if hgvs_genomic.posedit.edit.type == 'ins':
            start = hgvs_genomic.posedit.pos.start.base
            base = self.sf.fetch_seq(
                hgvs_genomic.ac,
                start_i=start - 1,
                end_i=start
            )
            alt = base + hgvs_genomic.posedit.edit.alt

            hgvs_genomic_forced_delins = (
                vvhgvs.sequencevariant.SequenceVariant(
                    ac=hgvs_genomic.ac,
                    type="g",
                    posedit=vvhgvs.posedit.PosEdit(
                        vvhgvs.location.Interval(
                            start=vvhgvs.location.SimplePosition(
                                base=start
                            ),
                            end=vvhgvs.location.SimplePosition(
                                base=start
                            ),
                            uncertain=hgvs_genomic.posedit.pos.uncertain
                        ),
                        vvhgvs.edit.NARefAlt(
                            ref=base,
                            alt=alt
                        )
                    )
                )
            )

        # Convert inversions to an equivalent delins without parsing.
        if hgvs_genomic.posedit.edit.type == 'inv':
            base = self.sf.fetch_seq(
                hgvs_genomic.ac,
                start_i=hgvs_genomic.posedit.pos.start.base - 1,
                end_i=hgvs_genomic.posedit.pos.end.base
            )

            alt = str(Seq(base).reverse_complement())

            hgvs_genomic = vvhgvs.sequencevariant.SequenceVariant(
                ac=hgvs_genomic.ac,
                type="g",
                posedit=vvhgvs.posedit.PosEdit(
                    vvhgvs.location.Interval(
                        start=vvhgvs.location.SimplePosition(
                            base=hgvs_genomic.posedit.pos.start.base
                        ),
                        end=vvhgvs.location.SimplePosition(
                            base=hgvs_genomic.posedit.pos.end.base
                        ),
                        uncertain=hgvs_genomic.posedit.pos.uncertain
                    ),
                    vvhgvs.edit.NARefAlt(
                        ref=base,
                        alt=alt
                    )
                )
            )

        # Project genomic variant onto overlapping transcripts.
        code_var = []

        for tx_ac in rts:
            try:
                variant = evm.g_to_t(
                    hgvs_genomic,
                    tx_ac
                )

            except vvhgvs.exceptions.HGVSError:
                curr_genomic = (
                    hgvs_genomic_forced_delins
                    if hgvs_genomic_forced_delins is not None
                    else hgvs_genomic
                )

                try:
                    variant = evm.g_to_t(
                        curr_genomic,
                        tx_ac
                    )
                except vvhgvs.exceptions.HGVSError:
                    continue

            except Exception as err:
                logger.info(
                    'non expected err type %s',
                    err
                )
                continue

            try:
                reverse_normalizer.normalize(variant)

            except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                if (
                        "Unsupported normalization of variants spanning the "
                        "exon-intron boundary" in str(e)
                        and variant.posedit.edit.type == "ins"
                ):
                    variant.posedit.pos.end.base = (
                        variant.posedit.pos.start.base
                    )
                    variant.posedit.pos.end.offset = 1

            # Corrective normalization of intronic and 3-prime UTR
            # descriptions in the antisense orientation.
            if (
                    hgvs_position_utils.either_position_is_intronic(variant)
                    or hgvs_position_utils.start_is_3_prime_utr(variant)
                    or hgvs_position_utils.end_is_3_prime_utr(variant)
            ):
                tx_ac = variant.ac

                try:
                    if rts_dict[tx_ac] < 0:
                        rev_hgvs_genomic = reverse_normalizer.normalize(
                            hgvs_genomic
                        )
                        variant = evm.g_to_t(
                            rev_hgvs_genomic,
                            tx_ac
                        )

                except vvhgvs.exceptions.HGVSInvalidIntervalError:
                    try:
                        variant = evm.g_to_t(
                            hgvs_genomic,
                            tx_ac
                        )
                    except vvhgvs.exceptions.HGVSInvalidIntervalError:
                        pass

            code_var.append(variant)

        return code_var

    def validateHGVS(self, query):
        """
        Take HGVS string, parse into hgvs object and validate
        """
        if type(query) is str:
            hgvs_input = self.hp.parse_hgvs_variant(query)
        else:
            hgvs_input = query

        try:
            self.vr.validate(hgvs_input)
        except vvhgvs.exceptions.HGVSError as e:
            return e
        else:
            return 'false'

    def entrez_efetch(self, db, id, rettype, retmode):
        """
        Search Entrez databases with efetch and SeqIO
        """
        Entrez.email = self.entrez_email
        Entrez.tool = 'VariantValidator'
        if self.entrez_api_key:
            Entrez.api_key = self.entrez_api_key
        # from Bio import SeqIO
        handle = Entrez.efetch(db=db, id=id, rettype=rettype, retmode=retmode)
        # Get record
        record = SeqIO.read(handle, "gb")
        # Place into text
        handle.close()
        return record

    def revcomp(self, bases):
        """
        Return the reverse complement of a nucleotide sequence.
        """
        return fn.simple_dna_revcomp(bases)

    def merge_hgvs_3pr(
            self,
            hgvs_variant_list,
            hn,
            genomic_reference=False,
            final_norm=True,
            hgvs_strict=False,
            map_dat=None
    ):
        """
        Merge multiple HGVS variants into a single delins using 3-prime
        normalization.

        Production paths are expected to supply parsed HGVS objects. Unit tests
        may supply HGVS strings, which are parsed at the testing boundary.
        """
        h_list = []
        store_ref_type = ""

        c_to_g_mapped = {
            "mapped": False,
            "ori": None,
            "transcript": None
        }

        tx_map_dat = map_dat
        if not tx_map_dat:
            tx_map_dat = TranscriptMapData(hdp=self.hdp)

        # Prepare and validate submitted HGVS variants.
        for hgvs_v in hgvs_variant_list:

            # Validate the HGVS object BEFORE converting c. coordinates to the
            # internal n. representation used by the merge.
            try:
                self.vr.validate(hgvs_v)

            except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                if 'Cannot validate sequence of an intronic variant' in str(e):
                    if genomic_reference is not False:
                        c_to_g_mapped["mapped"] = True
                        c_to_g_mapped["ori"] = tx_map_dat.mapped_exons(
                            hgvs_v.ac,
                            genomic_reference,
                            alt_aln_method=self.alt_aln_method
                        )[0][3]
                        c_to_g_mapped["transcript"] = hgvs_v.ac
                    else:
                        raise fn.mergeHGVSerror(
                            "AlleleSyntaxError: Intronic variants can only be "
                            "validated if a genomic/gene reference sequence is "
                            "also provided e.g. "
                            "NC_000017.11(NM_000088.3):c.589-1G>T"
                        )
                else:
                    raise

            except AssertionError:
                raise AlleleSyntaxError(
                    f"AlleleVariantError: {hgvs_v} is not a valid HGVS variant "
                    f"description. Please submit individually for additional "
                    f"guidance"
                )

            # Convert coding coordinates only AFTER validation.
            if hgvs_v.type == 'c':
                store_ref_type = 'c'

                try:
                    hgvs_v = self.vm.c_to_n(hgvs_v)
                except Exception:
                    raise fn.mergeHGVSerror(
                        "AlleleSyntaxError: Unable to map from c. position to "
                        "absolute position"
                    )

            elif hgvs_v.type not in ('g', 'n', 'm'):
                raise fn.mergeHGVSerror(
                    "AlleleSyntaxError: Unsupported HGVS reference type"
                )

            h_list.append(hgvs_v)

        # Map intronic transcript variants to the supplied genomic reference.
        if c_to_g_mapped["mapped"]:
            mapped_list = []

            for hgvs_v in h_list:
                try:
                    hgvs_v = self.vm.n_to_c(hgvs_v)
                except Exception:
                    pass

                logger.info("merge_hgvs_3pr before t_to_g: %s", hgvs_v)

                hgvs_v = self.vm.t_to_g(
                    hgvs_v,
                    genomic_reference,
                    alt_aln_method=self.alt_aln_method
                )

                logger.info("merge_hgvs_3pr after t_to_g: %s", hgvs_v)

                mapped_list.append(hgvs_v)

            h_list = mapped_list

            if c_to_g_mapped["ori"] == -1:
                h_list.reverse()

        # All downstream merge logic operates on the prepared HGVS objects,
        # rather than the original c. objects supplied by the caller.
        hgvs_variant_list = h_list

        # Identity preprocessing is required by gap handling.
        #
        # Do not simplify this independently of the gap-mapping tests. Identity
        # variants are re-created below, but their initial coordinate contributes
        # to construction of the merged interval.
        elec = 0
        cp_h_list = []

        for hgvs_v in h_list:
            if hgvs_v.posedit.edit.type == 'identity':
                if elec == 0:
                    hgvs_v.posedit.pos.end.base = (
                        hgvs_v.posedit.pos.start.base
                    )
                    hgvs_v.posedit.edit.ref = hgvs_v.posedit.edit.ref[0]
                    hgvs_v.posedit.edit.alt = hgvs_v.posedit.edit.alt[0]
                    cp_h_list.append(hgvs_v)
                continue

            cp_h_list.append(hgvs_v)

        h_list = cp_h_list

        # Construct the complete interval to merge.
        full_list = []
        accession = None
        seqtype = None
        merge_start_pos = None
        merge_end_pos = None

        for hgvs_v in h_list:

            # Intronic positions require a genomic reference.
            try:
                if (
                        hgvs_v.posedit.pos.start.offset != 0
                        and genomic_reference is False
                ):
                    raise fn.mergeHGVSerror(
                        "Base-offset position submitted"
                    )

                if (
                        hgvs_v.posedit.pos.end.offset != 0
                        and genomic_reference is False
                ):
                    raise fn.mergeHGVSerror(
                        "Base-offset position submitted"
                    )

            except AttributeError as e:
                logger.debug("Except passed, %s", e)

            try:
                hgvs_v = hn.normalize(hgvs_v)
            except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                pass

            if accession is None:
                accession = hgvs_v.ac
                seqtype = hgvs_v.type

            elif hgvs_v.ac != accession:
                raise fn.mergeHGVSerror(
                    "AlleleSyntaxError: More than one reference sequence submitted"
                )

            if merge_start_pos is None:
                merge_start_pos = hgvs_v.posedit.pos.start.base
                merge_end_pos = hgvs_v.posedit.pos.end.base
                full_list.append(hgvs_v)
                continue

            if hgvs_v.posedit.pos.start.base <= merge_end_pos:
                raise fn.mergeHGVSerror(
                    "AlleleSyntaxError: Submitted variants are out of order or "
                    "their ranges overlap"
                )

            # Re-create sequence between submitted variants as identity.
            if hgvs_v.posedit.pos.start.base > merge_end_pos + 1:
                identity_sequence = self.sf.fetch_seq(
                    hgvs_v.ac,
                    merge_end_pos,
                    hgvs_v.posedit.pos.start.base - 1
                )

                full_list.append(
                    hgvs_delins_parts_to_hgvs_obj(
                        hgvs_v.ac,
                        hgvs_v.type,
                        merge_end_pos + 1,
                        identity_sequence,
                        identity_sequence,
                        end=hgvs_v.posedit.pos.start.base - 1
                    )
                )

            merge_end_pos = hgvs_v.posedit.pos.end.base
            full_list.append(hgvs_v)

        # Strict HGVS merge-rule handling.
        check_frame_restore = False
        p_reference = None
        cp_hgvs_variant_list = hgvs_variant_list

        if hgvs_strict:
            merge_within_bases = 1

            # Preserve independent HGVS objects because this path may manipulate
            # them below.
            cp_hgvs_variant_list = copy.deepcopy(
                hgvs_variant_list
            )

            for vt in range(len(hgvs_variant_list) - 1):
                v1 = hgvs_variant_list[vt]
                v2 = hgvs_variant_list[vt + 1]

                vn1 = hn.normalize(v1)
                vn2 = self.reverse_hn.normalize(v2)

                distance = (
                        vn2.posedit.pos.start.base -
                        vn1.posedit.pos.end.base
                )

                if (
                        distance > merge_within_bases
                        and (
                        vn1.posedit.pos.start.base -
                        vn2.posedit.pos.end.base
                ) < 3
                        and store_ref_type == 'c'
                        and vn1.type != 'g'
                        and vn1.posedit.edit.type == 'sub'
                        and vn2.posedit.edit.type == 'sub'
                ):
                    tx_info = self.hdp.get_tx_identity_info(vn1.ac)
                    same_aa_span = False

                    for i in range(
                            int(tx_info[3] + 1),
                            int(tx_info[4]) + 1,
                            3
                    ):
                        if (
                                vn1.posedit.pos.end.base >= i
                                and vn2.posedit.pos.start.base <= i + 2
                        ):
                            same_aa_span = True
                            break

                    if not same_aa_span:
                        cp_hgvs_variant_list.remove(v1)

                elif distance > merge_within_bases:
                    cp_hgvs_variant_list.remove(v1)

            if len(cp_hgvs_variant_list) == 1:
                if (
                        store_ref_type == 'c'
                        and hgvs_variant_list[0].type != 'g'
                ):
                    cp_hgvs_variant_list = copy.deepcopy(
                        hgvs_variant_list
                    )

                    first_fs = None
                    last_fs = None

                    p_reference = self.hdp.get_pro_ac_for_tx_ac(
                        cp_hgvs_variant_list[0].ac
                    )

                    for index, n_variant in enumerate(
                            cp_hgvs_variant_list
                    ):
                        c_variant = self.vm.n_to_c(n_variant)
                        p_variant = self.vm.c_to_p(
                            c_variant,
                            p_reference
                        )

                        if p_variant.posedit.edit.type == 'fs':
                            if first_fs is None:
                                first_fs = index
                            else:
                                last_fs = index

                    if first_fs is not None and last_fs is not None:
                        hgvs_variant_list = cp_hgvs_variant_list[
                                            first_fs:last_fs + 1
                                            ]
                        check_frame_restore = True

                    else:
                        return False

                else:
                    hgvs_variant_list = cp_hgvs_variant_list

            else:
                hgvs_variant_list = cp_hgvs_variant_list

        # Build the merged alternate sequence directly from HGVS objects.
        alt_sequence = ''

        for hgvs_v in full_list:
            ref_alt = hgvs_utils.hgvs_ref_alt(
                hgvs_v,
                self.sf
            )
            alt_sequence += ref_alt['alt']

        reference_sequence = self.sf.fetch_seq(
            accession,
            merge_start_pos - 1,
            merge_end_pos
        )

        hgvs_delins = hgvs_delins_parts_to_hgvs_obj(
            accession,
            seqtype,
            merge_start_pos,
            reference_sequence,
            alt_sequence,
            end=merge_end_pos
        )

        try:
            hgvs_delins = self.vm.n_to_c(hgvs_delins)
        except Exception as e:
            logger.debug("Except passed, %s", e)

        if final_norm:
            try:
                hgvs_delins = hn.normalize(hgvs_delins)
            except HGVSUnsupportedOperationError as e:
                logger.debug("Except passed, %s", e)

        if hgvs_strict and len(cp_hgvs_variant_list) > 1:
            merge_these = []

            if c_to_g_mapped["mapped"]:
                transcript_variants = []

                if c_to_g_mapped["ori"] == -1:
                    hgvs_variant_list.reverse()

                for genomic_variant in hgvs_variant_list:
                    transcript_variant = self.vm.g_to_t(
                        genomic_variant,
                        c_to_g_mapped["transcript"],
                        alt_aln_method=self.alt_aln_method
                    )
                    transcript_variants.append(
                        transcript_variant
                    )

                hgvs_delins = self.vm.g_to_t(
                    hgvs_delins,
                    c_to_g_mapped["transcript"],
                    alt_aln_method=self.alt_aln_method
                )

                hgvs_variant_list = transcript_variants

            # Convert to strings only at the warning/output boundary.
            for variant in hgvs_variant_list:
                try:
                    variant = self.vm.n_to_c(variant)
                except Exception:
                    pass

                merge_these.append(
                    fn.valstr(variant).split('.')[2]
                )

            if not check_frame_restore:
                raise AlleleSyntaxError(
                    f"AlleleSyntaxError: Variants "
                    f"[{';'.join(merge_these)}] should be merged into "
                    f"{fn.valstr(hgvs_delins)}"
                )

            if hgvs_delins.type == 'c':
                hgvs_delins_p = self.vm.c_to_p(
                    hgvs_delins,
                    p_reference
                )

                if 'fs' not in str(hgvs_delins_p.posedit.edit):
                    raise AlleleSyntaxError(
                        f"AlleleSyntaxError: Merging variants "
                        f"[{';'.join(merge_these)}] restores the original "
                        f"reading frame, so should be described as "
                        f"{fn.valstr(hgvs_delins)}"
                    )

        return hgvs_delins

    def hgvs_alleles(self, my_variant, genomic_reference=False):
        """
        HGVS allele handling function which takes a single HGVS allele description
        and separates each allele into a list of HGVS variants.
        """
        logger.info(
            "HGVS allele handling function with variant %s and genomic reference "
            "set to %s",
            my_variant.quibble,
            genomic_reference
        )

        try:
            accession, remainder = my_variant.quibble.split(':')
            logger.info(
                "Accession: %s and remainder: %s",
                accession,
                remainder
            )

            if ("(" in accession or ")" in accession) and not genomic_reference:
                raise fn.alleleVariantError(
                    f"AlleleVariantError: Unexpected compound accession "
                    f"'{accession}' passed to hgvs_alleles()"
                )

            def _parse_allele_part(accession, var_type, pe):
                try:
                    if var_type == 'c':
                        posedit = self.hp.parse_c_posedit(pe)
                    elif var_type == 'g':
                        posedit = self.hp.parse_g_posedit(pe)
                    elif var_type == 'm':
                        posedit = self.hp.parse_m_posedit(pe)
                    elif var_type == 'n':
                        posedit = self.hp.parse_n_posedit(pe)
                    elif var_type == 'r':
                        logger.info(
                            "RNA variant %s identified with accession %s",
                            my_variant.quibble,
                            accession
                        )
                        raise fn.alleleVariantError(
                            "UnsupportedFormatError: RNA allele syntax variants "
                            "are not currently supported. Please submit "
                            "individually for additional guidance"
                        )
                    else:
                        raise fn.alleleVariantError(
                            f"UnsupportedFormatError: Unsupported HGVS reference "
                            f"type '{var_type}'"
                        )

                except vvhgvs.exceptions.HGVSError:
                    raise AlleleSyntaxError(
                        f"AlleleVariantError: {accession}:{var_type}.{pe} is not "
                        f"a valid HGVS variant description. Please submit "
                        f"individually for additional guidance"
                    )

                return vvhgvs.sequencevariant.SequenceVariant(
                    ac=accession,
                    type=var_type,
                    posedit=posedit
                )

            def _check_and_fix_for_ex_repeat(
                    accession,
                    var_type,
                    pe,
                    genomic_reference
            ):
                """
                Detect expanded repeat syntax within an allele and convert it to a
                normalised HGVS SequenceVariant for downstream processing.

                Returns
                -------
                tuple
                    (repeat_variant_or_None, genomic_reference)
                """
                if (
                        not pe.endswith("]")
                        or not re.search(r"[GATC]+\[\d+\]$", pe)
                ):
                    return None, genomic_reference

                logger.info(
                    "Checking allele %s for expanded repeats",
                    pe
                )

                expanded_variant = (
                    expanded_repeats.TandemRepeats.parse_repeat_variant(
                        f"{accession}:{var_type}.{pe}",
                        my_variant.primary_assembly,
                        "all",
                        self,
                    )
                )

                if expanded_variant is False:
                    return None, genomic_reference

                repeat_to_delins = expanded_variant.reformat(self)
                repeat_to_delins.posedit.expanded_rep = False

                try:
                    repeat_to_delins = self.hn.normalize(repeat_to_delins)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                    pass

                logger.info(
                    "Expanded repeat in allele normalised to %s",
                    repeat_to_delins
                )

                if (
                        genomic_reference is False
                        and repeat_to_delins.type in ("c", "n")
                        and (
                        repeat_to_delins.posedit.pos.start.offset != 0
                        or repeat_to_delins.posedit.pos.end.offset != 0
                )
                ):
                    logger.info(
                        "Looking up genomic reference for transcript %s",
                        repeat_to_delins.ac
                    )

                    for option in self.hdp.get_tx_mapping_options(
                            repeat_to_delins.ac
                    ):
                        genomic_ac = option[1]

                        if seq_data.get_chr_num_refseq(
                                genomic_ac,
                                my_variant.primary_assembly
                        ) is not None:
                            genomic_reference = genomic_ac

                            logger.info(
                                "Using genomic reference %s",
                                genomic_reference
                            )
                            break
                else:
                    logger.info(
                        "Created return using genomic reference %s",
                        genomic_reference
                    )

                return repeat_to_delins, genomic_reference

            def _parse_posedits(posedits, prefix=''):
                """
                Parse one semicolon-separated allele into HGVS objects.

                The returned list contains HGVS SequenceVariant objects only.
                """
                nonlocal genomic_reference

                current_allele = []

                for pe in posedits.split(';'):
                    if '?' in pe or pe == '0':
                        continue

                    pe = prefix + pe

                    tandem, genomic_reference = _check_and_fix_for_ex_repeat(
                        accession,
                        var_type,
                        pe,
                        genomic_reference
                    )

                    if tandem:
                        current_allele.append(tandem)
                    else:
                        current_allele.append(
                            _parse_allele_part(
                                accession,
                                var_type,
                                pe
                            )
                        )

                return current_allele

            def _validate_merges(alleles):
                """
                Apply strict HGVS merge validation to each non-empty allele.

                merge_hgvs_3pr() is called for its validation behaviour here; its
                returned merged variant is not required.
                """
                for each_allele in alleles:
                    if not each_allele:
                        continue

                    self.merge_hgvs_3pr(
                        each_allele,
                        my_variant.hn,
                        genomic_reference,
                        hgvs_strict=True
                    )

            # Shared-positions allele syntax:
            # NM_004006.2:c.2376[G>C];[G>C]
            if re.search(r'[gcn]\.\d+\[', remainder):
                var_type, remainder = remainder.split('.', 1)

                pos_match = re.match(r'\d+', remainder)
                pos = pos_match.group(0)

                # Remove the shared position and surrounding allele brackets.
                remainder = remainder[len(pos):]
                remainder = remainder[1:-1]

                alleles = remainder.split('];[')
                my_alleles = []

                for posedit in alleles:
                    # NM_004006.2:c.2376[G>C];[(G>C)]
                    if '(' in posedit:
                        continue

                    current_allele = _parse_posedits(
                        posedit,
                        prefix=pos
                    )
                    my_alleles.append(current_allele)

            else:
                var_type, remainder = remainder.split('.', 1)

                if '(;)' in remainder and '];' in remainder:
                    # NM_004006.2:c.[296T>G];[476T>C](;)1083A>C(;)1406del
                    pre_alleles = remainder.split('(;)')

                    pre_merges = []
                    alleles = []

                    for allele in pre_alleles:
                        if '[' in allele:
                            pre_merges.append(allele)
                        else:
                            alleles.append(allele)

                    my_alleles = []

                    # Unbracketed alleles.
                    for posedits in alleles:
                        my_alleles.append(
                            _parse_posedits(posedits)
                        )

                    # Bracketed alleles requiring merge validation.
                    merge_remainder = ';'.join(pre_merges)
                    merge_remainder = merge_remainder[1:-1]

                    for posedits in merge_remainder.split('];['):
                        my_alleles.append(
                            _parse_posedits(posedits)
                        )

                    _validate_merges(my_alleles)

                    # Preserve the existing behaviour: return the individual
                    # variants after strict merge validation.
                    my_alleles = [
                        [variant]
                        for each_allele in my_alleles
                        if each_allele
                        for variant in each_allele
                    ]

                elif '(;)' in remainder:
                    # Uncertain phase without bracketed allele syntax.
                    #
                    # NM_004006.2:c.2376G>C(;)3103del
                    # NM_000548.3:c.3623_3647del(;)3745_3756dup
                    if '[' in remainder:
                        raise fn.alleleVariantError(
                            'Unsupported format ' +
                            var_type +
                            '.' +
                            remainder
                        )

                    my_alleles = [
                        _parse_posedits(posedits)
                        for posedits in remainder.split('(;)')
                    ]

                else:
                    # Standard bracketed allele syntax.
                    #
                    # NM_004006.2:c.[2376G>C];[3103del]
                    # NM_004006.2:c.[296T>G;476C>T;1083A>C];
                    #                  [296T>G;1083A>C]
                    # NM_000548.3:c.[4358_4359del;4361_4372del]
                    if '(' in remainder:
                        raise fn.alleleVariantError(
                            'Unsupported format ' +
                            var_type +
                            '.' +
                            remainder
                        )

                    remainder = remainder[1:-1]

                    my_alleles = [
                        _parse_posedits(posedits)
                        for posedits in remainder.split('];[')
                    ]

                    _validate_merges(my_alleles)

                    # Preserve the existing behaviour: return individual variants,
                    # not the merged representation.
                    my_alleles = [
                        [variant]
                        for each_allele in my_alleles
                        if each_allele
                        for variant in each_allele
                    ]

            # String conversion belongs at the output boundary.
            return [
                str(allele)
                for alleles_l in my_alleles
                for allele in alleles_l
            ]

        except Exception as e:
            exc_type, exc_value, last_traceback = sys.exc_info()
            logger.error("%s %s", exc_type, exc_value)
            raise fn.alleleVariantError(str(e))

    def chr_to_rsg(self, hgvs_genomic, hn):
        """
        Convert a chromosomal HGVS description to RefSeqGene.
        """
        hgvs_genomic = hn.normalize(hgvs_genomic)

        chr_ac = hgvs_genomic.ac
        chr_start_pos = hgvs_genomic.posedit.pos.start.base
        chr_end_pos = hgvs_genomic.posedit.pos.end.base
        chr_edit = hgvs_genomic.posedit.edit

        all_info = self.db.get_g_to_g_info(
            gen_id=chr_ac,
            start=chr_start_pos,
            end=chr_end_pos
        )

        descriptions = []

        for line in all_info:
            if not (
                    chr_ac == line[1]
                    and chr_start_pos >= int(line[2])
                    and chr_end_pos <= int(line[3])
            ):
                continue

            rsg_ac = line[0]
            rsg_start = int(line[2])
            rsg_end = int(line[3])
            ori = line[4]
            gene = line[5]

            edit = copy.deepcopy(chr_edit)

            if ori == '+':
                start = chr_start_pos - rsg_start + 1
                end = chr_end_pos - rsg_start + 1

            elif ori == '-':
                if edit.type in (
                        'del', 'delins', 'dup', 'sub', 'inv', 'identity'
                ):
                    edit.ref = self.revcomp(edit.ref)

                if edit.type in (
                        'ins', 'delins', 'sub', 'identity'
                ):
                    edit.alt = self.revcomp(edit.alt)

                rsg_length = rsg_end - rsg_start

                start = (
                        rsg_length
                        - (chr_end_pos - rsg_start)
                        + 1
                )
                end = (
                        rsg_length
                        - (chr_start_pos - rsg_start)
                        + 1
                )

            else:
                continue

            hgvs_refseqgene = hgvs_obj_from_existing_edit(
                rsg_ac,
                'g',
                start,
                edit,
                end=end,
                offset_pos=True
            )

            try:
                hgvs_refseqgene = hn.normalize(hgvs_refseqgene)
            except vvhgvs.exceptions.HGVSError:
                descriptions.append({
                    'hgvs_refseqgene': hgvs_refseqgene,
                    'gene': gene,
                    'valid': 'Not in SeqRepo'
                })
                continue

            try:
                self.vr.validate(hgvs_refseqgene)

            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)

                if 'does not agree with reference sequence' in error:
                    match = re.findall(r'\(([GATC]+)\)', error)
                    hgvs_refseqgene.posedit.edit.ref = match[1]
                    error = 'true'

                descriptions.append({
                    'hgvs_refseqgene': hgvs_refseqgene,
                    'gene': gene,
                    'valid': error
                })

            else:
                descriptions.append({
                    'hgvs_refseqgene': hgvs_refseqgene,
                    'gene': gene,
                    'valid': 'true'
                })

        return descriptions

    def rsg_to_chr(self, hgvs_refseqgene, primary_assembly, hn):
        """
        Convert a RefSeqGene HGVS description to chromosomal HGVS.

        :param hgvs_refseqgene:
        :param primary_assembly:
        :param hn: HGVS Normalizer
        :return:
        """
        try:
            hgvs_refseqgene = hn.normalize(hgvs_refseqgene)
        except vvhgvs.exceptions.HGVSError as e:
            logger.debug("Except passed, %s", e)

        rsg_ac = hgvs_refseqgene.ac
        rsg_start_pos = hgvs_refseqgene.posedit.pos.start.base
        rsg_end_pos = hgvs_refseqgene.posedit.pos.end.base
        rsg_edit = hgvs_refseqgene.posedit.edit

        all_info = self.db.get_g_to_g_info(rsg_id=rsg_ac)

        descriptions = []

        for line in all_info:
            if not (
                    rsg_ac == line[0]
                    and primary_assembly == line[6]
            ):
                continue

            chr_ac = line[1]
            chr_start = int(line[2])
            chr_end = int(line[3])
            ori = line[4]
            gene = line[5]

            edit = copy.deepcopy(rsg_edit)

            if ori == '+':
                start = chr_start + rsg_start_pos - 1
                end = chr_start + rsg_end_pos - 1

            elif ori == '-':
                if edit.type in (
                        'del', 'delins', 'dup', 'sub', 'inv', 'identity'
                ):
                    edit.ref = self.revcomp(edit.ref)

                if edit.type in (
                        'ins', 'delins', 'dup', 'sub', 'inv', 'identity'
                ):
                    edit.alt = self.revcomp(edit.alt)

                chr_length = chr_end - chr_start

                start = (
                        chr_start
                        + chr_length
                        - rsg_end_pos
                        + 1
                )
                end = (
                        chr_start
                        + chr_length
                        - rsg_start_pos
                        + 1
                )

            else:
                continue

            hgvs_genomic = hgvs_obj_from_existing_edit(
                chr_ac,
                'g',
                start,
                edit,
                end=end,
                offset_pos=True
            )

            hgvs_genomic = hn.normalize(hgvs_genomic)

            try:
                self.vr.validate(hgvs_genomic)

            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)

                if 'does not agree with reference sequence' in error:
                    match = re.findall(r'\(([GATC]+)\)', error)
                    hgvs_genomic.posedit.edit.ref = match[1]
                    error = 'true'

                descriptions.append({
                    'hgvs_genomic': str(hgvs_genomic),
                    'gene': gene,
                    'valid': error
                })

            else:
                descriptions.append({
                    'hgvs_genomic': str(hgvs_genomic),
                    'gene': gene,
                    'valid': 'true'
                })

        return descriptions

    def transcript_filter(self, rts, select_transcripts=None):
        """
        Filter transcript lists to the latest accession versions, or return
        explicitly selected transcripts.
        """
        if self.testing or select_transcripts == 'raw':
            return rts

        if select_transcripts in (None, 'all'):
            latest_version = {}

            for tx_id in rts:
                if isinstance(tx_id, str):
                    # VV method: remove dud transcript IDs.
                    if "/" in tx_id or "_NG" in tx_id:
                        continue

                    accession, version = tx_id.split(".")

                    try:
                        version_int = int(version)
                    except ValueError:
                        logger.info(
                            "Transcript version error detected in %s.%s",
                            accession,
                            version
                        )
                        continue

                    if accession not in latest_version:
                        latest_version[accession] = {
                            "version": version,
                            "version_int": version_int,
                            "list": None
                        }
                    elif version_int > latest_version[accession]["version_int"]:
                        latest_version[accession] = {
                            "version": version,
                            "version_int": version_int,
                            "list": None
                        }

                else:
                    # VF method.
                    try:
                        accession, _sep, version = tx_id[0].partition(".")
                    except (IndexError, TypeError):
                        continue

                    try:
                        version_int = int(version)
                    except ValueError:
                        logger.info(
                            "Transcript version error detected in %s.%s",
                            accession,
                            version
                        )
                        continue

                    if accession not in latest_version:
                        latest_version[accession] = {
                            "version": version,
                            "version_int": version_int,
                            "list": tx_id[1:]
                        }
                    elif version_int > latest_version[accession]["version_int"]:
                        latest_version[accession] = {
                            "version": version,
                            "version_int": version_int,
                            "list": tx_id[1:]
                        }

            # Recreate list containing only the latest versions.
            filtered_rts = []

            for accession, data in latest_version.items():
                transcript = f"{accession}.{data['version']}"

                if data["list"] is None:
                    filtered_rts.append(transcript)
                else:
                    filtered_rts.append(
                        [transcript] + data["list"]
                    )

            return filtered_rts

        # If select_transcripts is JSON, decode it. Otherwise, treat the
        # supplied value as a single transcript selection.
        try:
            return list(json.loads(select_transcripts))
        except json.decoder.JSONDecodeError:
            return [select_transcripts]


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
