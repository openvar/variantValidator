import copy
import re
import vvhgvs.exceptions
from . import hgvs_utils, hgvs_position_utils
from VariantValidator.modules.hgvs_utils import hgvs_delins_parts_to_hgvs_obj, hgvs_dup_to_delins
from VariantValidator.modules.variant import TranscriptMapData
from VariantValidator.modules.utils import simple_dna_revcomp


# New functions for handling insertion type gapped mappings, relies on, and should only trigger with
# hgvs mapping improvements, needs to be used before the main GapMapper object
# no_norm_evm, and map_dat could be sourced from a variant object but will sill
def _expand_gap_insertion(
        insertion,
        gap_ref,
        start_offset,
        end_offset,
        reverse=False,
):
    """
    Add unchanged alignment-gap sequence to an insertion.

    start_offset and end_offset describe how much sequence from gap_ref must
    be restored on either side of the insertion. For reverse-strand mappings,
    gap_ref is reverse complemented and the offsets are exchanged so that
    they remain relative to the insertion sequence being constructed.
    """
    insertion = insertion or ""

    if reverse:
        gap_ref = simple_dna_revcomp(gap_ref)
        start_offset, end_offset = end_offset, start_offset

    if start_offset:
        insertion = gap_ref[:start_offset] + insertion

    if end_offset:
        insertion += gap_ref[-end_offset:]

    return insertion

def immediate_round_trip_gap_ins_handling(
        orig_hgvs_coding,
        orig_hgvs_genomic,
        variant_or_validator,
        map_dat=False,
):
    """
    Correct insertion sequence lost when mapping across a transcript/genome
    alignment gap.

    This function operates on fresh, unnormalised mappings. Normalisation can
    move transcript and genomic variants relative to one another and can also
    change the represented edit, so this correction must happen immediately
    after mapping.

    A typical problematic alignment is:

        TTTIIICIITTT
        GGG      GGG

    where T is normally aligned transcript sequence, I is transcript sequence
    inserted relative to the genome, and C is the actual sequence change.

    Mapping may represent the change as ``insC`` while omitting unchanged
    transcript sequence associated with the alignment gap. Where the mapping
    state can be established safely by a round trip, this function restores
    those omitted bases.

    Some branches deliberately retain support for genomic-to-transcript
    mapping states that current HGVS mapping does not reliably produce. These
    are retained as future-facing defensive behaviour rather than treated as
    dead code.

    Returns the original objects unchanged when no safe correction can be
    established.
    """
    no_norm_evm = variant_or_validator.no_norm_evm

    map_dat = getattr(
        variant_or_validator,
        "map_dat",
        map_dat,
    )
    if not map_dat:
        map_dat = TranscriptMapData(
            hdp=variant_or_validator.hdp
        )

    if hgvs_position_utils.either_position_is_intronic(
            orig_hgvs_coding
    ):
        return orig_hgvs_coding, orig_hgvs_genomic

    coding_edit = orig_hgvs_coding.posedit.edit
    genomic_edit = orig_hgvs_genomic.posedit.edit

    coding_edit_type = coding_edit.type
    genomic_edit_type = genomic_edit.type

    tx_ac = orig_hgvs_coding.ac
    genomic_ac = orig_hgvs_genomic.ac

    strand = map_dat.map_strand(
        tx_ac,
        genomic_ac,
    )
    reverse = strand < 0

    # Establish the effective alternate sequence on each reference.
    if genomic_edit_type == "dup":
        genomic_alt = genomic_edit.ref * 2
        if reverse:
            genomic_alt = simple_dna_revcomp(genomic_alt)

    elif coding_edit_type == "inv":
        if reverse:
            genomic_alt = genomic_edit.ref
        else:
            genomic_alt = simple_dna_revcomp(
                genomic_edit.ref
            )

    else:
        genomic_alt = genomic_edit.alt

    if coding_edit_type == "dup":
        coding_alt = coding_edit.ref * 2

    elif coding_edit_type == "inv":
        coding_alt = simple_dna_revcomp(
            coding_edit.ref
        )

    else:
        coding_alt = coding_edit.alt

    if coding_alt != genomic_alt:
        return orig_hgvs_coding, orig_hgvs_genomic

    # Genomic insertion, transcript edit is not an insertion.
    if (
            genomic_edit_type == "ins"
            and coding_edit_type != "ins"
    ):
        eq_hgvs_genomic = copy.copy(orig_hgvs_genomic)
        eq_hgvs_genomic.posedit.edit.ref = ""
        eq_hgvs_genomic.posedit.edit.alt = ""

        remap_hgvs_coding = no_norm_evm.g_to_n(
            eq_hgvs_genomic,
            tx_ac,
        )

        remap_ref = remap_hgvs_coding.posedit.edit.ref

        if len(coding_edit.ref) + 2 >= len(remap_ref):
            return orig_hgvs_coding, orig_hgvs_genomic

        if orig_hgvs_coding.type == "c":
            n_orig_hgvs_coding = no_norm_evm.c_to_n(
                orig_hgvs_coding
            )
        else:
            n_orig_hgvs_coding = orig_hgvs_coding

        start_offset = (
            n_orig_hgvs_coding.posedit.pos.start.base
            - remap_hgvs_coding.posedit.pos.start.base
            - 1
        )
        end_offset = (
            remap_hgvs_coding.posedit.pos.end.base
            - n_orig_hgvs_coding.posedit.pos.end.base
            - 1
        )

        new_ins = _expand_gap_insertion(
            coding_edit.alt,
            remap_ref[1:-1],
            start_offset,
            end_offset,
            reverse=reverse,
        )

        new_hgvs_genomic = copy.copy(orig_hgvs_genomic)
        new_hgvs_genomic.posedit.edit.alt = new_ins

        return orig_hgvs_coding, new_hgvs_genomic

    # Both mappings are insertions.
    if (
            genomic_edit_type == "ins"
            and coding_edit_type == "ins"
    ):
        eq_hgvs_genomic = copy.copy(orig_hgvs_genomic)
        eq_hgvs_genomic.posedit.edit.ref = ""
        eq_hgvs_genomic.posedit.edit.alt = ""

        remap_hgvs_coding = no_norm_evm.g_to_n(
            eq_hgvs_genomic,
            tx_ac,
        )

        if orig_hgvs_coding.type == "c":
            n_orig_hgvs_coding = no_norm_evm.c_to_n(
                orig_hgvs_coding
            )
        else:
            n_orig_hgvs_coding = orig_hgvs_coding

        stored_ac = False

        if n_orig_hgvs_coding.rel_ac != genomic_ac:
            stored_ac = n_orig_hgvs_coding.rel_ac
            n_orig_hgvs_coding.rel_ac = genomic_ac

        remap_hgvs_genomic = no_norm_evm.n_to_g(
            n_orig_hgvs_coding
        )

        if stored_ac:
            n_orig_hgvs_coding.rel_ac = stored_ac

        remap_genomic_edit = remap_hgvs_genomic.posedit.edit

        # Retained future-facing G->T handling.
        if remap_genomic_edit.type != "ins":
            start_offset = (
                orig_hgvs_genomic.posedit.pos.start.base
                - remap_hgvs_genomic.posedit.pos.start.base
            )
            end_offset = (
                remap_hgvs_genomic.posedit.pos.end.base
                - orig_hgvs_genomic.posedit.pos.end.base
            )

            new_ins = _expand_gap_insertion(
                genomic_edit.alt,
                remap_genomic_edit.ref,
                start_offset,
                end_offset,
                reverse=reverse,
            )

            new_hgvs_coding = copy.copy(
                orig_hgvs_coding
            )
            new_hgvs_coding.posedit.edit.alt = new_ins

            return new_hgvs_coding, orig_hgvs_genomic

        remap_ref = remap_hgvs_coding.posedit.edit.ref

        if len(remap_ref) > 2:
            start_offset = (
                n_orig_hgvs_coding.posedit.pos.start.base
                - remap_hgvs_coding.posedit.pos.start.base
                - 1
            )
            end_offset = (
                remap_hgvs_coding.posedit.pos.end.base
                - n_orig_hgvs_coding.posedit.pos.end.base
                - 1
            )

            new_ins = _expand_gap_insertion(
                coding_edit.alt,
                remap_ref[1:-1],
                start_offset,
                end_offset,
                reverse=reverse,
            )

            new_hgvs_genomic = copy.copy(
                orig_hgvs_genomic
            )
            new_hgvs_genomic.posedit.edit.alt = new_ins

            return orig_hgvs_coding, new_hgvs_genomic

        return orig_hgvs_coding, orig_hgvs_genomic

    # Transcript insertion, genomic mapping is not an insertion.
    # Retained as future-facing G->T handling.
    if coding_edit_type == "ins":
        if orig_hgvs_coding.type == "c":
            n_orig_hgvs_coding = no_norm_evm.c_to_n(
                orig_hgvs_coding
            )
        else:
            n_orig_hgvs_coding = orig_hgvs_coding

        stored_ac = False

        if n_orig_hgvs_coding.rel_ac != genomic_ac:
            stored_ac = n_orig_hgvs_coding.rel_ac
            n_orig_hgvs_coding.rel_ac = genomic_ac

        remap_hgvs_genomic = no_norm_evm.n_to_g(
            n_orig_hgvs_coding
        )

        if stored_ac:
            n_orig_hgvs_coding.rel_ac = stored_ac

        remap_ref = remap_hgvs_genomic.posedit.edit.ref

        if (
                not remap_ref
                or len(genomic_edit.ref) >= len(remap_ref)
        ):
            return orig_hgvs_coding, orig_hgvs_genomic

        start_offset = (
            orig_hgvs_genomic.posedit.pos.start.base
            - remap_hgvs_genomic.posedit.pos.start.base
        )
        end_offset = (
            remap_hgvs_genomic.posedit.pos.end.base
            - orig_hgvs_genomic.posedit.pos.end.base
        )

        new_ins = _expand_gap_insertion(
            genomic_edit.alt,
            remap_ref,
            start_offset,
            end_offset,
            reverse=reverse,
        )

        new_hgvs_coding = copy.copy(orig_hgvs_coding)
        new_hgvs_coding.posedit.edit.alt = new_ins

        return new_hgvs_coding, orig_hgvs_genomic

    return orig_hgvs_coding, orig_hgvs_genomic


def _position_offset(position):
    """Return an HGVS position offset, treating absent/None as zero."""
    return getattr(position, "offset", 0) or 0


def _start_offset(hgvs_variant):
    return _position_offset(hgvs_variant.posedit.pos.start)


def _end_offset(hgvs_variant):
    return _position_offset(hgvs_variant.posedit.pos.end)


def _has_any_offset(hgvs_variant):
    """True when either HGVS interval boundary carries an intronic offset."""
    return hgvs_position_utils.either_position_is_intronic(hgvs_variant)


def _has_both_offsets(hgvs_variant):
    """True when both HGVS interval boundaries carry intronic offsets."""
    return (
        hgvs_position_utils.start_position_is_intronic(hgvs_variant)
        and hgvs_position_utils.end_position_is_intronic(hgvs_variant)
    )


def _same_posedit(left, right):
    """Compare HGVS PosEdit objects without stringifying them."""
    return left.posedit == right.posedit


def _same_edit(left, right):
    """Compare HGVS edit objects without stringifying them."""
    return left.posedit.edit == right.posedit.edit


class GapMapper:

    def __init__(self, variant, validator):
        """
        Initialise the gap mapper.

        :param variant: variant.Variant()
        :param validator: Validator()
        """
        self.variant = variant
        self.validator = validator
        self.gapped_transcripts = ""
        self.auto_info = ""
        self.orientation = None
        self.hgvs_genomic_possibilities = []
        self.disparity_deletion_in = []
        self.hgvs_genomic_5pr = None
        self.tx_hgvs_not_delins = None

    def _ensure_map_data_provider(self):
        """Attach the Validator data provider to TranscriptMapData when needed."""
        if not self.variant.map_dat.hdp:
            self.variant.map_dat.hdp = self.validator.hdp

    def _transcript_is_selected(
            self,
            hgvs_transcript,
            select_transcripts_dict,
    ):
        """Return whether a mapped transcript passes the configured selector."""
        selection = self.validator.select_transcripts
        tx_ac = hgvs_transcript.ac

        if (
                selection not in ("all", "raw")
                and "select" not in selection
                and "mane" not in selection
                and "refseqgene" not in selection
        ):
            return tx_ac.split(".")[0] in select_transcripts_dict

        if selection == "select":
            annotation = self.validator.db.get_transcript_annotation(tx_ac)
            return any(
                value in annotation
                for value in (
                    '"select": "MANE"',
                    '"select": "RefSeq"',
                    '"select": "Ensembl"',
                )
            )

        annotation_key = {
            "mane": ('"mane_select": true', '"mane_plus_clinical": true'),
            "mane_select": ('"mane_select": true',),
            "refseq_select": ('"refseq_select": true',),
            "ensembl_select": ('"ensembl_select": true',),
        }.get(selection)

        if annotation_key is None:
            return True

        annotation = self.validator.db.get_transcript_annotation(tx_ac)
        return any(value in annotation for value in annotation_key)

    def _update_gap_warning(
            self,
            tx_ac,
            gen_ac,
            current_warning="",
            message=None,
    ):
        """Build a gap warning and merge useful automatic gap information."""
        gap_warnings = self.make_gap_warnings(
            tx_ac,
            gen_ac,
            self.variant.primary_assembly,
            message=message,
        )

        warning = gap_warnings["gapped_alignment_warning"]
        auto_info = gap_warnings["auto_info"]

        if auto_info and ("fewer" in auto_info or "extra" in auto_info):
            self.auto_info += auto_info

        return warning if warning is not None else current_warning

    def _adjust_transcript_gap_offsets(
            self,
            rn_tx_hgvs_not_delins,
            saved_hgvs_coding,
            hgvs_not_delins,
    ):
        """
        Remove transcript intronic offsets before gap-length comparison.

        Offset direction is read directly from the HGVS position objects;
        this replaces repeated string/regex inspection while preserving the
        existing movement rules.
        """
        start_offset = _start_offset(rn_tx_hgvs_not_delins)
        end_offset = _end_offset(rn_tx_hgvs_not_delins)

        if start_offset and end_offset:
            return (
                self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins),
                hgvs_not_delins,
            )

        if end_offset > 0:
            return self.move_tx_end_base_to_next_nonoffset(
                rn_tx_hgvs_not_delins,
                saved_hgvs_coding,
                back=False,
            )

        if start_offset > 0:
            return self.move_tx_start_base_to_previous_nonoffset(
                rn_tx_hgvs_not_delins,
                saved_hgvs_coding,
            )

        if end_offset < 0:
            return self.move_tx_end_base_to_next_nonoffset(
                rn_tx_hgvs_not_delins,
                saved_hgvs_coding,
            )

        if start_offset < 0:
            return self.move_tx_start_base_to_previous_nonoffset(
                rn_tx_hgvs_not_delins,
                saved_hgvs_coding,
                with_base_subtract=True,
            )

        return rn_tx_hgvs_not_delins, hgvs_not_delins

    def make_gap_warnings(
            self,
            tx_ac,
            gen_ac,
            primary_assembly,
            message=None,
    ):
        """
        Generate warnings describing gaps in a transcript/genome alignment.

        CIGAR strings for the transcript/genome exon alignments are inspected
        for insertions and deletions. Gaps are reported using coding
        coordinates for NM_ transcripts and non-coding coordinates for NR_
        transcripts.

        This function is only called with NM_ or NR_ transcript accessions.

        If message is supplied, it replaces the automatically generated gap
        description while retaining the standard warning structure.
        """
        map_dat = self.variant.map_dat

        # TranscriptMapData normally already has the Validator data provider,
        # but retain support for construction paths where it has not yet been
        # assigned.
        if not map_dat.hdp:
            map_dat.hdp = self.validator.hdp

        tx_exons = map_dat.mapped_exons(
            tx_ac,
            gen_ac,
            alt_aln_method=self.validator.alt_aln_method,
        )

        # Retain only exon alignments containing an insertion or deletion.
        # The reduced structure preserves the fields used by the historical
        # gap-position calculation below.
        gap_in_alignment = []

        for exon in tx_exons:
            cigar = exon[9]

            if "I" not in cigar and "D" not in cigar:
                continue

            gap_in_alignment.append([
                exon[0],
                exon[1],
                exon[3],
                int(exon[5]) + 1,
                int(exon[6]),
                int(exon[7]),
                int(exon[8]) + 1,
                cigar,
            ])

        gap_information = {
            "gapped_alignment_warning": "",
            "auto_info": "",
        }

        if not gap_in_alignment:
            return gap_information

        found_gaps = []

        # NM_ transcripts use coding coordinates; NR_ transcripts use
        # non-coding coordinates. These are the only transcript namespaces
        # accepted by this function.
        coordinate_type = (
            "c." if tx_ac.startswith("NM_") else "n."
        )

        for gap_loc in gap_in_alignment:
            cigar = gap_loc[-1]

            # Separate individual CIGAR operations while retaining each
            # operation character with its preceding length.
            cigar_parts = (
                cigar
                .replace("=", "=:")
                .replace("I", "I:")
                .replace("D", "D:")
                .replace("X", "X:")
                .split(":")
            )

            tx_exon_start = gap_loc[3]

            tx_annotation = (
                self.validator.hdp.get_tx_identity_info(
                    gap_loc[0]
                )
            )

            try:
                cds_start = int(tx_annotation[3]) + 1
                cds_end = int(tx_annotation[4])
                tx_position = tx_exon_start - cds_start
            except TypeError:
                # NR_ transcripts have no coding-region bounds, so their
                # positions remain transcript-relative.
                cds_start = None
                cds_end = None
                tx_position = tx_exon_start

            for cigar_part in cigar_parts:
                # X and = both consume transcript sequence without
                # representing an insertion/deletion in the alignment.
                cigar_part = cigar_part.replace("X", "=")

                if "=" in cigar_part:
                    match_length = int(
                        cigar_part.split("=", 1)[0]
                    )
                    tx_position += match_length
                    continue

                if "D" in cigar_part:
                    gap_length = int(
                        cigar_part.split("D", 1)[0]
                    )

                    gap_position = (
                        f"{coordinate_type}{tx_position}_"
                        f"{tx_position + gap_length + 1}"
                    )

                    found_gaps.append(
                        f"{gap_length} extra bases between "
                        f"{gap_position}"
                    )

                    tx_position += gap_length
                    continue

                if "I" in cigar_part:
                    gap_length = int(
                        cigar_part.split("I", 1)[0]
                    )

                    gap_position = (
                        f"{coordinate_type}{tx_position}_"
                        f"{tx_position + 1}"
                    )

                    found_gaps.append(
                        f"{gap_length} fewer bases between "
                        f"{gap_position}"
                    )

                    tx_position += gap_length

        # Convert NM_ gap positions beyond the CDS into the historical
        # 3-prime UTR representation used by VariantValidator warnings.
        if tx_ac.startswith("NM_"):
            converted_gaps = []

            for found_gap in found_gaps:
                coordinates = found_gap.split("c.")[-1]
                start = int(
                    coordinates.split("_", 1)[0]
                )

                if start + cds_start >= cds_end:
                    utr_3 = start - cds_end

                    utr_3_position = (
                        f"*{utr_3 + cds_start}_"
                        f"*{utr_3 + cds_start + 1}"
                    )

                    found_gap = found_gap.replace(
                        coordinates,
                        utr_3_position,
                    )

                converted_gaps.append(found_gap)

            found_gaps = converted_gaps

        gap_string = (
            message
            if message is not None
            else ", and ".join(found_gaps)
        )

        gap_information["gapped_alignment_warning"] = (
            "Submitted description does not represent a true variant because "
            f"it is an artefact of aligning {tx_ac} with {gen_ac} "
            f"(genome build {primary_assembly})"
        )

        gap_information["auto_info"] = (
            f"{tx_ac} contains {gap_string} than {gen_ac}"
        )

        return gap_information

    def gapped_g_to_c(self, rel_var, select_transcripts_dict):
        """
        Gap aware projection from g. to c.
        """

        # RefSeq or Ensembl?
        expanded_genomic_for_ensembl = False
        if self.validator.alt_aln_method == 'genebuild':
            # Expand the genomic variant to include the flanking bases as a delins
            reverse_normalized_hgvs_genomic = self.validator.reverse_hn.normalize(self.variant.hgvs_genomic)

            # VCF
            vcf_dict = hgvs_utils.hgvs2vcf(reverse_normalized_hgvs_genomic, self.variant.primary_assembly,
                                           self.variant.reverse_normalizer, self.validator.sf,
                                           extra_flank_bases=4)

            expanded_genomic_for_ensembl = hgvs_delins_parts_to_hgvs_obj(
                    reverse_normalized_hgvs_genomic.ac,
                    reverse_normalized_hgvs_genomic.type,
                    int(vcf_dict['pos']),
                    vcf_dict['ref'],
                    vcf_dict['alt'])

        # Set variables for problem specific warnings
        gapped_alignment_warning = ''
        corrective_action_taken = ''
        self.gapped_transcripts = ''
        self.auto_info = ''
        self.disparity_deletion_in = []

        # set map data provider
        self._ensure_map_data_provider()

        # Create a pseudo VCF so that normalization can be applied and a delins can be generated
        hgvs_genomic_variant = self.variant.hgvs_genomic

        # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic_variant)
        self.hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

        # VCF
        vcf_dict = hgvs_utils.hgvs2vcf(reverse_normalized_hgvs_genomic, self.variant.primary_assembly,
                                       self.variant.reverse_normalizer, self.validator.sf)
        stored_hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                self.hgvs_genomic_5pr.ac,
                self.hgvs_genomic_5pr.type,
                int(vcf_dict['pos']),vcf_dict['ref'],vcf_dict['alt'])

        # take a look at the input genomic variant for potential base salvage
        stash_ac = vcf_dict['chr']
        stash_input = self.variant.post_format_conversion
        if type(stash_input) is str:
            stash_input = self.validator.hp.parse_hgvs_variant(stash_input)
        # Re-Analyse genomic positions
        if 'NG_' in str(self.variant.hgvs_formatted):
            c = rel_var[0]
            if hasattr(c.posedit.edit, 'ref') and c.posedit.edit.ref is not None:
                c.posedit.edit.ref = c.posedit.edit.ref.upper()
            if hasattr(c.posedit.edit, 'alt') and c.posedit.edit.alt is not None:
                c.posedit.edit.alt = c.posedit.edit.alt.upper()
            stash_input = self.validator.myevm_t_to_g(c, self.variant.no_norm_evm,
                                                      self.variant.primary_assembly, self.variant.hn, self.variant)
        if stash_input.ac.startswith('NC_') or stash_input.ac.startswith('NT_') or stash_input.ac.startswith('NW_'):
            hgvs_stash = stash_input
            if hasattr(hgvs_stash.posedit.edit, 'ref') and hgvs_stash.posedit.edit.ref is not None:
                hgvs_stash.posedit.edit.ref = hgvs_stash.posedit.edit.ref.upper()
            if hasattr(hgvs_stash.posedit.edit, 'alt') and hgvs_stash.posedit.edit.alt is not None:
                hgvs_stash.posedit.edit.alt = hgvs_stash.posedit.edit.alt.upper()

            # MAKE A NO NORM HGVS2VCF
            vcf_dict = hgvs_utils.pos_lock_hgvs2vcf(hgvs_stash, self.variant.primary_assembly,
                                                      self.variant.reverse_normalizer, self.validator.sf)
            stash_ac = hgvs_stash.ac

        # Store a not real deletion insertion for stashed
        stash_hgvs_not_delins =  hgvs_delins_parts_to_hgvs_obj(
                stash_ac,
                self.hgvs_genomic_5pr.type,
                int(vcf_dict['pos']),vcf_dict['ref'],vcf_dict['alt'])
        pre_loop_stash_hgvs_not_delins = copy.copy(stash_hgvs_not_delins)

        # make an empty rel_var
        nw_rel_var = []

        # loop through rel_var and amend where required
        for var in rel_var:
            stash_hgvs_not_delins = pre_loop_stash_hgvs_not_delins

            # Blank the self.tx_hgvs_not_delins
            self.tx_hgvs_not_delins = None

            # Store the current hgvs:c. description
            try:
                saved_hgvs_coding = self.validator.hp.parse_hgvs_variant(var)
                original_var = copy.deepcopy(saved_hgvs_coding)
            except TypeError:
                saved_hgvs_coding = var
                original_var = var
            except vvhgvs.exceptions.HGVSInvalidVariantError:
                saved_hgvs_coding = var
                original_var = var

            if not self._transcript_is_selected(
                    saved_hgvs_coding,
                    select_transcripts_dict,
            ):
                continue

            ## Only apply to known gapped alignment mappings
            if not self.variant.map_dat.hdp:
                self.variant.map_dat.hdp = self.validator.hdp
            if self.variant.map_dat.is_gapped_map(saved_hgvs_coding.ac,hgvs_genomic_variant.ac):
                # Applies to ensembl only
                if expanded_genomic_for_ensembl is not False:
                    try:
                        hgvs_refreshed_variant = self.validator.vm.g_to_t(expanded_genomic_for_ensembl,
                                                                      saved_hgvs_coding.ac,
                                                                      alt_aln_method=self.validator.alt_aln_method)
                    except vvhgvs.exceptions.HGVSInvalidIntervalError as e:
                        if "start or end or both are beyond the bounds of transcript record" in str(e):
                            continue
                    try:
                        hgvs_refreshed_variant = self.validator.vm.n_to_c(hgvs_refreshed_variant)
                    except vvhgvs.exceptions.HGVSError:
                        pass

                    # Get the ref length difference
                    genomic_ref_len = len(expanded_genomic_for_ensembl.posedit.edit.ref)
                    transcript_ref_len = len(hgvs_refreshed_variant.posedit.edit.ref)
                    if genomic_ref_len != transcript_ref_len:
                        message = None
                        if genomic_ref_len > transcript_ref_len:
                            gap_length = genomic_ref_len - transcript_ref_len
                            message = f"{gap_length} fewer bases"
                        elif genomic_ref_len < transcript_ref_len:
                            gap_length = transcript_ref_len - genomic_ref_len
                            message = f"{gap_length} extra bases"
                        else:
                            message = None

                        if message is not None:
                            gap_warnings = self.make_gap_warnings(hgvs_refreshed_variant.ac,
                                                                  expanded_genomic_for_ensembl.ac,
                                                                  self.variant.primary_assembly,
                                                                  message=message)
                            gapped_alignment_warning = gap_warnings["gapped_alignment_warning"]
                            if not self.auto_info:
                                self.auto_info = gap_warnings["auto_info"]

                        # Will filter out intronic variants since intronic variants will not normalize
                        try:
                            hgvs_refreshed_variant = self.validator.genebuild_normalizer_cross.normalize(
                                hgvs_refreshed_variant)
                        except vvhgvs.exceptions.HGVSError:
                            nw_rel_var.append(saved_hgvs_coding)
                        else:
                            nw_rel_var.append(hgvs_refreshed_variant)
                        continue
                    else:
                        nw_rel_var.append(saved_hgvs_coding)

                # This next section is looking for exonic gaps so cannot be applied to intronic positions
                needs_a_push = False
                merged_variant = False

                if not hgvs_position_utils.either_position_is_intronic(
                        saved_hgvs_coding
                ):
                    """
                    Directly search for gaps using vcf hard_pushing left
                    we pre-normalise the input (at least as far as basic hgvs), to test for sharing
                    the n->g mapping, which can be the most time consuming step in the process
                    """
                    genomic_mapping = None
                    non_variant_genomic_ac = False # push left/right require False genomic ac when main variant
                    # is genomic
                    try:
                        hgvs_right = copy.copy(saved_hgvs_coding)
                        if hgvs_right.type == 'c':
                            hgvs_right = self.validator.vm.c_to_n(hgvs_right)
                            hgvs_left = copy.copy(hgvs_right)
                        else:
                            hgvs_left = copy.copy(saved_hgvs_coding)
                        hgvs_right = self.variant.hn.normalize(hgvs_right)
                        hgvs_left = self.variant.reverse_normalizer.normalize(hgvs_left)
                        if hgvs_right.type == 'n':
                            non_variant_genomic_ac = hgvs_genomic_variant.ac

                        if hgvs_right == hgvs_left and hgvs_right.type == 'n':
                            genomic_mapping = self.validator.vm.n_to_g(hgvs_right, hgvs_genomic_variant.ac,alt_aln_method=self.validator.alt_aln_method)
                        else:
                            genomic_mapping = False

                        vcf__dict = hgvs_utils.hard_left_hgvs2vcf(saved_hgvs_coding,
                                                                  self.variant.primary_assembly,
                                                                  self.variant.hn,
                                                                  self.variant.reverse_normalizer,
                                                                  self.validator.sf,
                                                                  saved_hgvs_coding.ac,
                                                                  self.variant.map_dat,
                                                                  self.validator.alt_aln_method,
                                                                  self.validator.hp,
                                                                  self.validator.vm,
                                                                  self.validator.merge_hgvs_3pr,
                                                                  genomic_ac=non_variant_genomic_ac,
                                                                  mapped_g=copy.copy(genomic_mapping),
                                                                  pre_norm=hgvs_left)

                        if vcf__dict['needs_a_push'] is True:
                            needs_a_push = True
                            merged_variant = vcf__dict['merged_variant']
                            if merged_variant is not False:
                                try:
                                    merged_variant = self.validator.vm.n_to_c(merged_variant)
                                except (TypeError,
                                        vvhgvs.exceptions.HGVSInvalidVariantError,
                                        vvhgvs.exceptions.HGVSUsageError):
                                    pass

                        vcf__dict = hgvs_utils.hard_right_hgvs2vcf(hgvs_right,
                                                                   self.variant.primary_assembly,
                                                                   self.variant.hn,
                                                                   self.variant.reverse_normalizer,
                                                                   self.validator.sf,
                                                                   saved_hgvs_coding.ac,
                                                                   self.variant.map_dat,
                                                                   self.validator.alt_aln_method,
                                                                   self.validator.hp,
                                                                   self.validator.vm,
                                                                   self.validator.merge_hgvs_3pr,
                                                                   genomic_ac=non_variant_genomic_ac,
                                                                   mapped_g=genomic_mapping,
                                                                   pre_norm=hgvs_right)

                        if vcf__dict['needs_a_push'] is True:
                            needs_a_push = True
                            merged_variant = vcf__dict['merged_variant']
                            if merged_variant is not False:
                                try:
                                    merged_variant = self.validator.vm.n_to_c(merged_variant)
                                except (TypeError, vvhgvs.exceptions.HGVSInvalidVariantError):
                                    pass

                    except (vvhgvs.exceptions.HGVSUnsupportedOperationError,
                            vvhgvs.exceptions.HGVSInvalidVariantError,
                            vvhgvs.exceptions.HGVSUsageError,
                            vvhgvs.exceptions.HGVSDataNotAvailableError,
                            ValueError):
                        pass

                    # Collect the hard_pushed variant information and adjust the variants accordingly
                    if needs_a_push is not False:
                        if merged_variant is not False:
                            saved_hgvs_coding = merged_variant
                            stash_hgvs_not_delins = self.validator.vm.t_to_g(saved_hgvs_coding,
                                                                             hgvs_genomic_variant.ac,
                                                                             alt_aln_method=self.validator.alt_aln_method)

                            # The merged variant may have created an ins or a del
                            if stash_hgvs_not_delins.posedit.edit.type == "del":
                                stash_hgvs_not_delins.posedit.edit.alt = ""
                            if stash_hgvs_not_delins.posedit.edit.type == "ins":
                                get_ref = copy.deepcopy(stash_hgvs_not_delins)
                                get_ref.posedit.edit.ref = ''
                                get_ref.posedit.edit.alt = ''
                                get_ref = self.variant.hn.normalize(get_ref)
                                ref_bases = get_ref.posedit.edit.ref
                                stash_hgvs_not_delins.posedit.edit.ref = ref_bases
                                stash_hgvs_not_delins.posedit.edit.alt = ref_bases[0] \
                                                                         + stash_hgvs_not_delins.posedit.edit.alt \
                                                                         + ref_bases[1]

                # Get orientation of the gene wrt genome and a list of exons mapped to the genome
                ori = self.variant.map_dat.tx_exons(
                        saved_hgvs_coding.ac, self.hgvs_genomic_5pr.ac,
                        self.validator.alt_aln_method,
                        hdp=self.validator.hdp)
                try:
                    self.orientation = int(ori[0]['alt_strand'])
                except TypeError:
                    continue
                # Set intronic params
                intronic_variant = 'false'
                hgvs_seek_var = self.get_hgvs_seek_var(self.variant.hgvs_genomic, saved_hgvs_coding)

                if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                        saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                        hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                        saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset):
                    pass
                else:
                    hgvs_seek_var = saved_hgvs_coding

                try:
                    self.variant.hn.normalize(hgvs_seek_var)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                    error = str(e)
                    if 'Normalization of intronic variants is not supported' in error or \
                            'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                            intronic_variant = 'hard_fail'
                        else:
                            # Double check to see whether the variant is actually intronic?
                            for exon in ori:
                                genomic_start = int(exon['alt_start_i'])
                                genomic_end = int(exon['alt_end_i'])
                                if genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end and \
                                        genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end:
                                    intronic_variant = 'false'
                                    break
                                else:
                                    intronic_variant = 'true'
                except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                    if "insertion length must be 1" in str(e):
                        pass

                if _has_any_offset(hgvs_seek_var):

                    # Double check to see whether the variant is actually intronic?
                    for exon in ori:
                        genomic_start = int(exon['alt_start_i'])
                        genomic_end = int(exon['alt_end_i'])
                        if (genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                                genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                            intronic_variant = 'false'
                            break
                        else:
                            intronic_variant = 'true'

                # If exonic, process
                if intronic_variant != 'true' and intronic_variant != 'hard_fail':

                    # Attempt to find gaps in reference sequence by catching disparity in genome length and
                    # overlapping transcript lengths
                    self.disparity_deletion_in = ['false', 'false']
                    hgvs_not_delins = ''

                    if stored_hgvs_not_delins != '':
                        # Refresh hgvs_not_delins from stored_hgvs_not_delins
                        try:
                            hgvs_not_delins = self.dup_ins_5prime_shift(stored_hgvs_not_delins, saved_hgvs_coding)
                        except vvhgvs.exceptions.HGVSInvalidIntervalError:
                            hgvs_not_delins = stored_hgvs_not_delins

                        try:
                            self.tx_hgvs_not_delins = self.variant.no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                      saved_hgvs_coding.ac)
                        except vvhgvs.exceptions.HGVSInvalidIntervalError as e:
                            if "start or end or both are beyond the bounds of transcript record" in str(e):
                                self.tx_hgvs_not_delins = saved_hgvs_coding
                            else:
                                self.tx_hgvs_not_delins = self.variant.no_norm_evm.g_to_n(self.hgvs_genomic_5pr,
                                                                                          saved_hgvs_coding.ac)
                        except vvhgvs.exceptions.HGVSError as e:
                            if str(e) == 'start or end or both are beyond the bounds of transcript record':
                                self.tx_hgvs_not_delins = saved_hgvs_coding

                        # Create normalized version of tx_hgvs_not_delins
                        rn_tx_hgvs_not_delins = copy.deepcopy(self.tx_hgvs_not_delins)

                        rn_tx_hgvs_not_delins, hgvs_not_delins = (
                            self._adjust_transcript_gap_offsets(
                                rn_tx_hgvs_not_delins,
                                saved_hgvs_coding,
                                hgvs_not_delins,
                            )
                        )

                        # Logic
                        if len(hgvs_not_delins.posedit.edit.ref) < len(rn_tx_hgvs_not_delins.posedit.edit.ref):
                            gap_length = len(rn_tx_hgvs_not_delins.posedit.edit.ref) - len(
                                hgvs_not_delins.posedit.edit.ref)
                            self.disparity_deletion_in = ['chromosome', gap_length]
                        elif len(hgvs_not_delins.posedit.edit.ref) > len(rn_tx_hgvs_not_delins.posedit.edit.ref):
                            gap_length = len(hgvs_not_delins.posedit.edit.ref) - len(
                                rn_tx_hgvs_not_delins.posedit.edit.ref)
                            self.disparity_deletion_in = ['transcript', gap_length]
                        else:
                            # store stash_hgvs_not_delins for restorstion after error below
                            restore_stash_hgvs_not_delins = copy.copy(stash_hgvs_not_delins)
                            try:
                                hgvs_stash_t = self.validator.vm.g_to_t(stash_hgvs_not_delins, saved_hgvs_coding.ac,
                                                                        alt_aln_method=self.validator.alt_aln_method)
                            except vvhgvs.exceptions.HGVSError as e:
                                if 'bounds' in str(e):
                                    stash_hgvs_not_delins = copy.copy(stored_hgvs_not_delins)
                                    try:
                                        hgvs_stash_t = self.validator.vm.g_to_t(stash_hgvs_not_delins,
                                                                                saved_hgvs_coding.ac,
                                                                                alt_aln_method=self.validator.alt_aln_method)
                                    except vvhgvs.exceptions.HGVSError:
                                        hgvs_stash_t = saved_hgvs_coding

                            if len(stash_hgvs_not_delins.posedit.edit.ref) > len(hgvs_stash_t.posedit.edit.ref):
                                try:
                                    self.variant.hn.normalize(hgvs_stash_t)
                                except Exception as e:
                                    pass
                                else:
                                    gap_length = len(stash_hgvs_not_delins.posedit.edit.ref) - len(
                                        hgvs_stash_t.posedit.edit.ref)
                                    self.disparity_deletion_in = ['transcript', gap_length]
                                    try:
                                        self.tx_hgvs_not_delins = self.validator.vm.c_to_n(hgvs_stash_t)
                                    except:
                                        self.tx_hgvs_not_delins = hgvs_stash_t
                                    hgvs_not_delins = stash_hgvs_not_delins
                            elif hgvs_stash_t.posedit.pos.start.offset != 0 or hgvs_stash_t.posedit.pos.end.offset != 0:
                                self.disparity_deletion_in = ['transcript', 'Requires Analysis']
                                try:
                                    self.tx_hgvs_not_delins = self.validator.vm.c_to_n(hgvs_stash_t)
                                except:
                                    self.tx_hgvs_not_delins = hgvs_stash_t
                                hgvs_not_delins = stash_hgvs_not_delins
                                self.hgvs_genomic_5pr = stash_hgvs_not_delins
                            else:
                                try:
                                    var_a = self.variant.hn.normalize(hgvs_stash_t)
                                    var_b = self.variant.hn.normalize(original_var)
                                except vvhgvs.exceptions.HGVSError:
                                    pass
                                else:
                                    if var_a.posedit.edit.type != var_b.posedit.edit.type:
                                        gap_warnings = self.make_gap_warnings(self.tx_hgvs_not_delins.ac,
                                                                              self.hgvs_genomic_5pr.ac,
                                                                              self.variant.primary_assembly)

                                        if gap_warnings["gapped_alignment_warning"] is not None \
                                                and gap_warnings["auto_info"] is not None:
                                            gapped_alignment_warning = gap_warnings["gapped_alignment_warning"]
                                            if ("fewer" in gap_warnings["auto_info"] or
                                                    "extra" in gap_warnings["auto_info"]):
                                                self.auto_info = self.auto_info + gap_warnings["auto_info"]

                            # Restore stash_hgvs_not_delins
                            stash_hgvs_not_delins = restore_stash_hgvs_not_delins

                    # Final sanity checks
                    try:
                        self.validator.vm.g_to_t(hgvs_not_delins, self.tx_hgvs_not_delins.ac,
                                                 alt_aln_method=self.validator.alt_aln_method)
                    except Exception as e:
                        if str(e) == 'start or end or both are beyond the bounds of transcript record':
                            hgvs_not_delins = saved_hgvs_coding
                            self.disparity_deletion_in = ['false', 'false']
                    try:
                        self.variant.hn.normalize(self.tx_hgvs_not_delins)
                    except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                        error = str(e)
                        if 'Normalization of intronic variants is not supported' in error or \
                                'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                            if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                                hgvs_not_delins = saved_hgvs_coding
                                self.disparity_deletion_in = ['false', 'false']
                            elif 'Normalization of intronic variants is not supported' in error:
                                # We know that this cannot be because of an intronic variant, so must be aligned to
                                # tx gap
                                self.disparity_deletion_in = ['transcript', 'Requires Analysis']

                    # Pre-processing of self.tx_hgvs_not_delins
                    try:
                        if self.tx_hgvs_not_delins.posedit.edit.alt is None:
                            self.tx_hgvs_not_delins.posedit.edit.alt = ''
                    except Exception as e:
                        if str(e) == "'Dup' object has no attribute 'alt'":
                            self.tx_hgvs_not_delins = hgvs_dup_to_delins(self.tx_hgvs_not_delins)

                    # GAP IN THE TRANSCRIPT DISPARITY DETECTED
                    if self.disparity_deletion_in[0] == 'transcript':
                        pass
                        # Check for issue https://github.com/openvar/variantValidator/issues/385 where the gap is
                        # being identified but oddly the vm is not compensating, likely due to odd sequence
                        try:
                            if len(self.tx_hgvs_not_delins.posedit.edit.ref) > \
                                    len(self.tx_hgvs_not_delins.posedit.edit.alt):
                                gen_len_difference = len(hgvs_not_delins.posedit.edit.ref) - \
                                                     len(hgvs_not_delins.posedit.edit.alt)
                                tx_len_difference = len(self.tx_hgvs_not_delins.posedit.edit.ref) - \
                                                    len(self.tx_hgvs_not_delins.posedit.edit.alt)
                            else:
                                gen_len_difference = len(hgvs_not_delins.posedit.edit.alt) - \
                                                     len(hgvs_not_delins.posedit.edit.ref)
                                tx_len_difference = len(self.tx_hgvs_not_delins.posedit.edit.alt) - \
                                                    len(self.tx_hgvs_not_delins.posedit.edit.ref)

                            # The logic here. Since there is a gap in the transcript,
                            # the actual length should be == gen_len_difference - 1 not == gen_len_difference
                            if tx_len_difference - self.disparity_deletion_in[1] == gen_len_difference:
                                # So here we know we need to knock off disparity_deletion_in[1] bases
                                if len(hgvs_not_delins.posedit.edit.alt) == len(
                                        self.tx_hgvs_not_delins.posedit.edit.alt):
                                    if self.orientation == 1:
                                        self.tx_hgvs_not_delins.posedit.edit.ref = hgvs_not_delins.posedit.ref
                                    else:
                                        replace_ref_bases = self.validator.revcomp(hgvs_not_delins.posedit.edit.ref)
                                        self.tx_hgvs_not_delins.posedit.edit.ref = replace_ref_bases
                                    self.tx_hgvs_not_delins.posedit.pos.end.offset = self.disparity_deletion_in[1]

                        except TypeError:
                            pass
                        except AttributeError:
                            pass

                        gap_warnings = self.make_gap_warnings(self.tx_hgvs_not_delins.ac,
                                                              self.hgvs_genomic_5pr.ac,
                                                              self.variant.primary_assembly)

                        if gap_warnings["gapped_alignment_warning"] is not None \
                                and gap_warnings["auto_info"] is not None:
                            gapped_alignment_warning = gap_warnings["gapped_alignment_warning"]
                            if ("fewer" in gap_warnings["auto_info"] or
                                    "extra" in gap_warnings["auto_info"]):
                                self.auto_info = self.auto_info + gap_warnings["auto_info"]

                        # ANY VARIANT WHOLLY WITHIN THE GAP
                        hgvs_refreshed_variant = self.transcript_disparity(
                            reverse_normalized_hgvs_genomic,
                            self.variant.hgvs_genomic,
                            1,
                        )

                        # Look for missed duplications into direct complte in-line repeats
                        try:
                            if (hasattr(hgvs_refreshed_variant.posedit.pos.start, "offset")
                             or hasattr(hgvs_refreshed_variant.posedit.pos.end, "offset")):
                                pass
                                if(hgvs_refreshed_variant.posedit.pos.start.base ==
                                        rn_tx_hgvs_not_delins.posedit.pos.start.base
                                        and
                                    hgvs_refreshed_variant.posedit.pos.end.base ==
                                        rn_tx_hgvs_not_delins.posedit.pos.end.base
                                        and
                                    hgvs_refreshed_variant.posedit.edit.alt ==
                                        rn_tx_hgvs_not_delins.posedit.edit.alt):
                                            pass
                                            hgvs_refreshed_variant = rn_tx_hgvs_not_delins
                                            try:
                                                hgvs_refreshed_variant = (self.validator.vm.
                                                                          n_to_c(hgvs_refreshed_variant))
                                            except vvhgvs.exceptions.HGVSError:
                                                pass
                                            pass
                                            hgvs_refreshed_variant = self.variant.hn.normalize(hgvs_refreshed_variant)
                                            pass
                        except Exception:
                            pass

                    # GAP IN THE CHROMOSOME
                    elif self.disparity_deletion_in[0] == 'chromosome':
                        # Set warning variables
                        gap_warnings = self.make_gap_warnings(self.tx_hgvs_not_delins.ac,
                                                              self.hgvs_genomic_5pr.ac,
                                                              self.variant.primary_assembly)

                        if gap_warnings["gapped_alignment_warning"] is not None \
                                and gap_warnings["auto_info"] is not None:
                            gapped_alignment_warning = gap_warnings["gapped_alignment_warning"]
                            if ("fewer" in gap_warnings["auto_info"] or
                                    "extra" in gap_warnings["auto_info"]):
                                self.auto_info = self.auto_info + gap_warnings["auto_info"]
                        hgvs_refreshed_variant = self.tx_hgvs_not_delins

                    else:
                        # Have we already had a hard push?
                        # if needs_a_push is not False:
                        if merged_variant is not False:
                            hgvs_refreshed_variant = saved_hgvs_coding
                        else:
                            try:
                                # Try the push to see if a gap is identified
                                hgvs_stash = copy.deepcopy(stash_hgvs_not_delins)
                                stash_ac = hgvs_stash.ac
                                genomic_mapping = None
                                genomic_ac_for_tx = False
                                hgvs_genomic_variant.ac
                                hgvs_right = copy.copy(hgvs_stash)
                                if hgvs_right.type == 'c':
                                    hgvs_right = self.validator.vm.c_to_n(hgvs_right)
                                    hgvs_left = copy.copy(hgvs_right)
                                else:
                                    hgvs_left = copy.copy(hgvs_stash)
                                hgvs_right = self.variant.hn.normalize(hgvs_right)
                                hgvs_left = self.variant.reverse_normalizer.normalize(hgvs_left)
                                if hgvs_left.type == 'n':
                                    genomic_ac_for_tx = hgvs_genomic_variant.ac
                                if hgvs_right == hgvs_left and hgvs_right.type == 'n':
                                    genomic_mapping = self.validator.vm.n_to_g(
                                            hgvs_right,
                                            hgvs_genomic_variant.ac,
                                            alt_aln_method=self.validator.alt_aln_method)

                                # Make a hard left and hard right not delins g.
                                stash_dict_right = hgvs_utils.hard_right_hgvs2vcf(hgvs_stash,
                                                                                  self.variant.primary_assembly,
                                                                                  self.variant.hn,
                                                                                  self.variant.reverse_normalizer,
                                                                                  self.validator.sf,
                                                                                  saved_hgvs_coding.ac,
                                                                                  self.variant.map_dat,
                                                                                  self.validator.alt_aln_method,
                                                                                  self.validator.hp,
                                                                                  self.validator.vm,
                                                                                  self.validator.merge_hgvs_3pr,
                                                                                  genomic_ac=genomic_ac_for_tx,
                                                                                  mapped_g=copy.copy(genomic_mapping),
                                                                                  pre_norm=hgvs_right)
                                stash_hgvs_not_delins_right = hgvs_delins_parts_to_hgvs_obj(
                                        stash_ac,
                                        hgvs_stash.type,
                                        int(stash_dict_right['pos']),
                                        stash_dict_right['ref'],
                                        stash_dict_right['alt'])

                                stash_dict_left = hgvs_utils.hard_left_hgvs2vcf(hgvs_stash,
                                                                                self.variant.primary_assembly,
                                                                                self.variant.hn,
                                                                                self.variant.reverse_normalizer,
                                                                                self.validator.sf,
                                                                                saved_hgvs_coding.ac,
                                                                                self.variant.map_dat,
                                                                                self.validator.alt_aln_method,
                                                                                self.validator.hp,
                                                                                self.validator.vm,
                                                                                self.validator.merge_hgvs_3pr,
                                                                                genomic_ac=genomic_ac_for_tx,
                                                                                mapped_g=copy.copy(genomic_mapping),
                                                                                pre_norm=hgvs_left)
                                stash_hgvs_not_delins_left =  hgvs_delins_parts_to_hgvs_obj(
                                        stash_ac,
                                        hgvs_stash.type,
                                        int(stash_dict_left['pos']),
                                        stash_dict_left['ref'],
                                        stash_dict_left['alt'])
                            except vvhgvs.exceptions.HGVSDataNotAvailableError:
                                continue

                            # Map in-situ to the transcript left and right
                            try:
                                tx_hard_right = self.validator.vm.g_to_t(stash_hgvs_not_delins_right,
                                                                         saved_hgvs_coding.ac,
                                                                         alt_aln_method=self.validator.alt_aln_method)
                            except Exception:
                                tx_hard_right = saved_hgvs_coding
                            else:
                                normalize_stash_right = self.variant.hn.normalize(stash_hgvs_not_delins_right)
                                if _same_posedit(normalize_stash_right, stash_hgvs_not_delins):
                                    tx_hard_right = saved_hgvs_coding
                            try:
                                tx_hard_left = self.validator.vm.g_to_t(stash_hgvs_not_delins_left,
                                                                        saved_hgvs_coding.ac,
                                                                        alt_aln_method=self.validator.alt_aln_method)
                            except Exception:
                                tx_hard_left = saved_hgvs_coding
                            else:
                                normalize_stash_left = self.variant.hn.normalize(stash_hgvs_not_delins_left)
                                if _same_posedit(normalize_stash_left, stash_hgvs_not_delins):
                                    tx_hard_left = saved_hgvs_coding

                            try:
                                # The Logic - Currently limited to genome gaps
                                if len(stash_hgvs_not_delins_right.posedit.edit.ref) < len(
                                        tx_hard_right.posedit.edit.ref) or \
                                        len(stash_hgvs_not_delins_right.posedit.edit.ref) > \
                                        len(tx_hard_right.posedit.edit.ref):
                                    tx_hard_right = self.variant.hn.normalize(tx_hard_right)

                                    gap_warnings = self.make_gap_warnings(self.tx_hgvs_not_delins.ac,
                                                                          self.hgvs_genomic_5pr.ac,
                                                                          self.variant.primary_assembly)

                                    if gap_warnings["gapped_alignment_warning"] is not None \
                                            and gap_warnings["auto_info"] is not None:
                                        gapped_alignment_warning = gap_warnings["gapped_alignment_warning"]
                                        if ("fewer" in gap_warnings["auto_info"] or
                                                "extra" in gap_warnings["auto_info"]):
                                            self.auto_info = self.auto_info + gap_warnings["auto_info"]
                                    hgvs_refreshed_variant = tx_hard_right

                                elif len(stash_hgvs_not_delins_left.posedit.edit.ref) < \
                                        len(tx_hard_left.posedit.edit.ref) or \
                                        len(stash_hgvs_not_delins_left.posedit.edit.ref) > \
                                        len(tx_hard_left.posedit.edit.ref):
                                    tx_hard_left = self.variant.hn.normalize(tx_hard_left)
                                    gap_warnings = self.make_gap_warnings(self.tx_hgvs_not_delins.ac,
                                                                          self.hgvs_genomic_5pr.ac,
                                                                          self.variant.primary_assembly)

                                    if gap_warnings["gapped_alignment_warning"] is not None \
                                            and gap_warnings["auto_info"] is not None:
                                        gapped_alignment_warning = gap_warnings["gapped_alignment_warning"]
                                        if ("fewer" in gap_warnings["auto_info"] or
                                                "extra" in gap_warnings["auto_info"]):
                                            self.auto_info = self.auto_info + gap_warnings["auto_info"]
                                    hgvs_refreshed_variant = tx_hard_left

                                else:
                                    # Keep the same by re-setting rel_var
                                    hgvs_refreshed_variant = saved_hgvs_coding

                            except TypeError:
                                # e.g. chr1:156561557G>GGGGTC (investigate at a later date)
                                hgvs_refreshed_variant = saved_hgvs_coding

                            except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                                # e.g. NG_005895.1:g.3684_44407del
                                hgvs_refreshed_variant = saved_hgvs_coding

                    # Edit the output
                    hgvs_refreshed_variant = self.edit_output(hgvs_refreshed_variant, saved_hgvs_coding)

                    nw_rel_var.append(hgvs_refreshed_variant)

                # Otherwise these variants need to be set
                else:
                    if not gapped_alignment_warning:
                        gapped_alignment_warning = ''
                    if not corrective_action_taken:
                        corrective_action_taken = ''
                    # Send to empty nw_rel_var
                    nw_rel_var.append(saved_hgvs_coding)

            # Otherwise these variants need to be set
            else:
                if not gapped_alignment_warning:
                    gapped_alignment_warning = ''
                if not corrective_action_taken:
                    corrective_action_taken = ''
                # Send to empty nw_rel_var
                nw_rel_var.append(saved_hgvs_coding)

        data = {
            'gapped_alignment_warning': gapped_alignment_warning,
            'corrective_action_taken': corrective_action_taken,
            'auto_info': self.auto_info,
            'disparity_deletion_in': self.disparity_deletion_in,
            'gapped_transcripts': self.gapped_transcripts
        }

        return data, nw_rel_var

    def g_to_t_compensation(self, ori, hgvs_coding, rec_var):
        self.orientation = int(ori[0]['alt_strand'])
        self.hgvs_genomic_possibilities = []
        hgvs_genomic = self.validator.myevm_t_to_g(hgvs_coding, self.variant.no_norm_evm, self.variant.primary_assembly,
                                                   self.variant.hn, self.variant)

        rn_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic)
        self.hgvs_genomic_possibilities.append([rn_hgvs_genomic, ['false', 'false']])

        try:
            if self.orientation != -1:
                chromosome_normalized_hgvs_coding = self.variant.reverse_normalizer.normalize(hgvs_coding)
            else:
                chromosome_normalized_hgvs_coding = self.variant.hn.normalize(hgvs_coding)
        except vvhgvs.exceptions.HGVSUnsupportedOperationError:
            chromosome_normalized_hgvs_coding = hgvs_coding

        most_3pr_hgvs_genomic = self.validator.myvm_t_to_g(chromosome_normalized_hgvs_coding, hgvs_genomic.ac,
                                                           self.variant.no_norm_evm, self.variant.hn,
                                                           self.variant.map_dat)
        self.hgvs_genomic_possibilities.append([most_3pr_hgvs_genomic, ['false', 'false']])

        # Push from side to side to try pick up odd placements
        # MAKE A NO NORM HGVS2VCF
        # First to the right
        hgvs_stash = copy.deepcopy(hgvs_coding)
        stash_tx_right = ''
        stash_tx_left = ''
        map_fail = False
        self._ensure_map_data_provider()

        try:
            if hgvs_stash.type == 'c':
                hgvs_stash = self.validator.vm.c_to_n(hgvs_stash)
            hgvs_right = copy.copy(hgvs_stash)
            hgvs_left = copy.copy(hgvs_stash)
            if not (
                    getattr(hgvs_right.posedit.pos.start,'offset',False) or
                    getattr(hgvs_right.posedit.pos.end,'offset',False)):
                hgvs_right = self.variant.hn.normalize(hgvs_right)
                hgvs_left = self.variant.reverse_normalizer.normalize(hgvs_left)
            if hgvs_right == hgvs_left and hgvs_right.type == 'n':
                genomic_mapping = self.validator.vm.n_to_g(
                        hgvs_right,
                        hgvs_genomic.ac,
                        alt_aln_method=self.validator.alt_aln_method)
            else:
                genomic_mapping = False
        except (vvhgvs.exceptions.HGVSUnsupportedOperationError,
                vvhgvs.exceptions.HGVSInvalidVariantError,
                vvhgvs.exceptions.HGVSUsageError,
                vvhgvs.exceptions.HGVSDataNotAvailableError):
            map_fail = True
        try:
            if map_fail:
                raise ValueError("Already failed in n->g mapping")
            stash_ac = hgvs_stash.ac
            stash_dict = hgvs_utils.hard_right_hgvs2vcf(hgvs_stash,
                                                        self.variant.primary_assembly,
                                                        self.variant.hn,
                                                        self.variant.reverse_normalizer,
                                                        self.validator.sf,
                                                        hgvs_coding.ac,
                                                        self.variant.map_dat,
                                                        self.validator.alt_aln_method,
                                                        self.validator.hp,
                                                        self.validator.vm,
                                                        self.validator.merge_hgvs_3pr,
                                                        genomic_ac=hgvs_genomic.ac,
                                                        mapped_g=copy.copy(genomic_mapping),
                                                        pre_norm=hgvs_right)
            # make a not real deletion insertion
            stash_hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                    stash_ac,
                    hgvs_stash.type,
                    int(stash_dict['pos']),stash_dict['ref'],stash_dict['alt'],
                    offset_pos=True)
            try:
                stash_hgvs_not_delins = self.variant.no_norm_evm.n_to_c(stash_hgvs_not_delins)
            except Exception as e:
                pass

            test_stash_tx_right = copy.deepcopy(stash_hgvs_not_delins)
            stash_genomic = self.validator.myvm_t_to_g(test_stash_tx_right, hgvs_genomic.ac, self.variant.no_norm_evm,
                                                       self.variant.hn,self.variant.map_dat)
            if len(test_stash_tx_right.posedit.edit.ref) == (
                    (stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                stash_tx_right = test_stash_tx_right
                if hasattr(test_stash_tx_right.posedit.edit,
                           'alt') and test_stash_tx_right.posedit.edit.alt is not None:
                    alt = test_stash_tx_right.posedit.edit.alt
                else:
                    alt = ''
                if hasattr(stash_genomic.posedit.edit,
                           'alt') and stash_genomic.posedit.edit.alt is not None:
                    g_alt = stash_genomic.posedit.edit.alt
                else:
                    g_alt = ''
                if (len(alt) - (
                        test_stash_tx_right.posedit.pos.end.base - test_stash_tx_right.posedit.pos.start.base) + 1) != (
                        len(g_alt) - (
                        stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])
                else:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
            elif test_stash_tx_right.posedit.edit.type == 'identity':
                hgvs_reform_ident = hgvs_delins_parts_to_hgvs_obj(
                        test_stash_tx_right.ac,
                        'c',
                        test_stash_tx_right.posedit.pos,
                        test_stash_tx_right.posedit.edit.ref,
                        '')
                try:
                    self.variant.hn.normalize(hgvs_reform_ident)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if "spanning the exon-intron boundary" in error:
                        stash_tx_right = test_stash_tx_right
                        self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_right = test_stash_tx_right
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])
            else:
                try:
                    self.variant.hn.normalize(test_stash_tx_right)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_right = test_stash_tx_right
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])

            # IDENTIFY GAP AND HARD SET
            if stash_dict['needs_a_push'] is True or stash_dict['identifying_g_variant'] is not False:
                # Look for merged variant from hard push
                if stash_dict['merged_variant'] is not False:
                    merged_variant = stash_dict['merged_variant']
                    try:
                        merged_variant = self.validator.vm.n_to_c(merged_variant)
                    except (TypeError, vvhgvs.exceptions.HGVSInvalidVariantError):
                        pass
                    identifying_variant = stash_dict['identifying_variant']
                    try:
                        identifying_variant = self.validator.vm.n_to_c(identifying_variant)
                    except (TypeError, vvhgvs.exceptions.HGVSInvalidVariantError):
                        pass
                    stash_genomic = self.validator.myvm_t_to_g(identifying_variant, stash_genomic.ac,
                                                               self.variant.no_norm_evm,self.variant.hn,
                                                               self.variant.map_dat)
                    stash_hgvs_not_delins = merged_variant
                    test_stash_tx_right = merged_variant

                elif stash_dict['identifying_g_variant'] is not False:
                    stash_hgvs_not_delins = stash_dict['identifying_g_variant']
                    stash_genomic = stash_dict['identifying_g_variant']

                # Look for gap info
                normalized_stash_genomic = self.variant.hn.normalize(stash_genomic)
                stash_tx_right = test_stash_tx_right
                if stash_hgvs_not_delins.posedit.edit.type == "ins":
                    len_tx = 2
                else:
                    len_tx = len(stash_hgvs_not_delins.posedit.edit.ref)
                if stash_genomic.posedit.edit.type == "ins":
                    len_gen = 2
                else:
                    len_gen = len(stash_genomic.posedit.edit.ref)
                if len_tx > len_gen:
                    gap_in = 'chromosome'
                    gap_len = len_tx - len_gen
                else:
                    gap_in = 'transcript'
                    gap_len = len_gen - len_tx

                # Set the options to a single option based on the results of pushing
                self.hgvs_genomic_possibilities = [[normalized_stash_genomic, [gap_in,
                                                                               gap_len,
                                                                               stash_hgvs_not_delins,
                                                                               stash_genomic]]]
        # Intronic positions not supported. Will cause a Value Error
        except (vvhgvs.exceptions.HGVSError, ValueError) as e:
            pass

        # Then to the left
        try:
            if map_fail:
                raise ValueError("Already failed in n->g mapping")
            stash_ac = hgvs_stash.ac
            stash_dict = hgvs_utils.hard_left_hgvs2vcf(hgvs_stash,
                                                       self.variant.primary_assembly,
                                                       self.variant.hn,
                                                       self.variant.reverse_normalizer,
                                                       self.validator.sf,
                                                       hgvs_coding.ac,
                                                       self.variant.map_dat,
                                                       self.validator.alt_aln_method,
                                                       self.validator.hp,
                                                       self.validator.vm,
                                                       self.validator.merge_hgvs_3pr,
                                                       genomic_ac=hgvs_genomic.ac,
                                                       mapped_g=copy.copy(genomic_mapping),
                                                       pre_norm=hgvs_left)
            # make a not real deletion insertion
            stash_hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                    stash_ac,
                    hgvs_stash.type,
                    int(stash_dict['pos']),stash_dict['ref'],stash_dict['alt'],
                    offset_pos=True)
            try:
                stash_hgvs_not_delins = self.variant.no_norm_evm.n_to_c(stash_hgvs_not_delins)
            except Exception as e:
                pass
                # Store a tx copy for later use
            test_stash_tx_left = copy.deepcopy(stash_hgvs_not_delins)
            stash_genomic = self.validator.myvm_t_to_g(test_stash_tx_left, hgvs_genomic.ac, self.variant.no_norm_evm,
                                                       self.variant.hn,self.variant.map_dat)

            if len(test_stash_tx_left.posedit.edit.ref) == ((stash_genomic.posedit.pos.end.base -
                                                             stash_genomic.posedit.pos.start.base) + 1):
                if hasattr(test_stash_tx_left.posedit.edit,
                           'alt') and test_stash_tx_left.posedit.edit.alt is not None:
                    alt = test_stash_tx_left.posedit.edit.alt
                else:
                    alt = ''
                if hasattr(stash_genomic.posedit.edit,
                           'alt') and stash_genomic.posedit.edit.alt is not None:
                    g_alt = stash_genomic.posedit.edit.alt
                else:
                    g_alt = ''

                if (len(alt) - (
                        test_stash_tx_left.posedit.pos.end.base - test_stash_tx_left.posedit.pos.start.base) + 1) != (
                        len(g_alt) - (
                        stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])
                else:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
            elif test_stash_tx_left.posedit.edit.type == 'identity':
                hgvs_reform_ident = hgvs_delins_parts_to_hgvs_obj(
                        test_stash_tx_left.ac,
                        'c',
                        test_stash_tx_left.posedit.pos,
                        test_stash_tx_left.posedit.edit.ref,
                        '')
                try:
                    self.variant.hn.normalize(hgvs_reform_ident)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if "spanning the exon-intron boundary" in error:
                        stash_tx_left = test_stash_tx_left
                        self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_left = test_stash_tx_left
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])
            else:
                try:
                    self.variant.hn.normalize(test_stash_tx_left)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_left = test_stash_tx_left
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])

            # IDENTIFY GAP AND HARD SET
            if stash_dict['needs_a_push'] is True or stash_dict['identifying_g_variant'] is not False:

                # Look for merged variant from hard push
                if stash_dict['merged_variant'] is not False:
                    merged_variant = stash_dict['merged_variant']
                    try:
                        merged_variant = self.validator.vm.n_to_c(merged_variant)
                    except (TypeError, vvhgvs.exceptions.HGVSInvalidVariantError):
                        pass
                    identifying_variant = stash_dict['identifying_variant']
                    try:
                        identifying_variant = self.validator.vm.n_to_c(identifying_variant)
                    except (TypeError, vvhgvs.exceptions.HGVSInvalidVariantError):
                        pass
                    stash_genomic = self.validator.myvm_t_to_g(identifying_variant, stash_genomic.ac,
                                                               self.variant.no_norm_evm,self.variant.hn,
                                                               self.variant.map_dat)
                    stash_hgvs_not_delins = merged_variant
                    test_stash_tx_left = merged_variant

                elif stash_dict['identifying_g_variant'] is not False:
                    stash_hgvs_not_delins = stash_dict['identifying_g_variant']
                    stash_genomic = stash_dict['identifying_g_variant']

                # Look for gap info
                normalized_stash_genomic = self.variant.hn.normalize(stash_genomic)
                stash_tx_left = test_stash_tx_left

                if stash_hgvs_not_delins.posedit.edit.type == "ins":
                    len_tx = 2
                else:
                    len_tx = len(stash_hgvs_not_delins.posedit.edit.ref)
                if stash_genomic.posedit.edit.type == "ins":
                    len_gen = 2
                else:
                    len_gen = len(stash_genomic.posedit.edit.ref)
                if len_tx > len_gen:
                    gap_in = 'chromosome'
                    gap_len = len_tx - len_gen
                else:
                    gap_in = 'transcript'
                    gap_len = len_gen - len_tx

                # Set the options to a single option based on the results of pushing
                # This code refers to https://github.com/openvar/variantValidator/issues/651
                if stash_dict["pre_merged_variant"] is not False and stash_dict["identifying_g_variant"] is False:
                    if ("transcript" in gap_in
                            and len_tx == len_gen
                            and stash_dict["pre_merged_variant"].type == "g")\
                            and normalized_stash_genomic.posedit.edit.type == "sub"\
                            and stash_dict["pre_merged_variant"].posedit.edit.type == "delins"\
                            and hgvs_coding.posedit.edit.type == "sub":
                        normalized_stash_genomic = stash_dict["pre_merged_variant"]

                self.hgvs_genomic_possibilities = [[normalized_stash_genomic, [gap_in,
                                                                               gap_len,
                                                                               stash_hgvs_not_delins,
                                                                               stash_genomic]]]
        # Intronic positions not supported. Will cause a Value Error
        except (vvhgvs.exceptions.HGVSError, ValueError)as e:
            pass

        # direct mapping from reverse_normalized transcript insertions in the delins format
        self.rev_norm_ins(hgvs_coding, hgvs_genomic)

        # Set variables for problem specific warnings
        self.gapped_transcripts = ''
        self.auto_info = ''

        # Mark as not disparity detected
        self.disparity_deletion_in = ['false', 'false']

        # Loop through to see if a gap can be located
        # Set the variables required for corrective normalization
        possibility_counter = 0
        suppress_c_normalization = 'false'  # Applies to boundary crossing normalization

        # If hard pushing identified a variant (in theory there can be only 1) then remove all other possibilities
        hard_possibility = []
        for hard_set_check in self.hgvs_genomic_possibilities:
            if hard_set_check[1] != ['false', 'false']:
                hard_possibility.append(hard_set_check)

        # Copy a version of hgvs_genomic_possibilities
        for a_possibility in self.hgvs_genomic_possibilities:
            possibility = a_possibility[0]
            disparity_info = a_possibility[1]
            possibility_counter = possibility_counter + 1

            # Loop out stash possibilities which will not spot gaps so are empty
            if possibility == '':
                continue

            # Use VCF generation code to push hgvs_genomic as for 5 prime as possible to uncover gaps
            hgvs_genomic_variant = copy.deepcopy(possibility)
            reverse_normalized_hgvs_genomic = ''
            # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
            try:
                reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic_variant)
            except vvhgvs.exceptions.HGVSError as e:
                # Strange error caused by gap in genomic
                error = str(e)
                if 'base start position must be <= end position' in error:
                    if hgvs_genomic.posedit.edit.type == 'delins':
                        start = hgvs_genomic.posedit.pos.start.base
                        end = hgvs_genomic.posedit.pos.end.base
                        lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                        rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                        hgvs_genomic.posedit.edit.ref = lhb + rhb
                        hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                        hgvs_genomic.posedit.pos.start.base = end
                        hgvs_genomic.posedit.pos.end.base = start
                        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic)
                    if hgvs_genomic.posedit.edit.type == 'del':
                        start = hgvs_genomic.posedit.pos.start.base
                        end = hgvs_genomic.posedit.pos.end.base
                        lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                        rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                        hgvs_genomic.posedit.edit.ref = lhb + rhb
                        hgvs_genomic.posedit.edit.alt = lhb + rhb
                        hgvs_genomic.posedit.pos.start.base = end
                        hgvs_genomic.posedit.pos.end.base = start
                        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic)
                if 'insertion length must be 1' in error:
                    if hgvs_genomic.posedit.edit.type == 'ins':
                        start = hgvs_genomic.posedit.pos.start.base
                        end = hgvs_genomic.posedit.pos.end.base
                        lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                        rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start, end)
                        hgvs_genomic.posedit.edit.ref = lhb + rhb
                        hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic)

            self.hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

            # Create VCF
            vcf_dict = hgvs_utils.hgvs2vcf(reverse_normalized_hgvs_genomic, self.variant.primary_assembly,
                                           self.variant.reverse_normalizer, self.validator.sf)
            # Store a not real deletion insertion to test for gapping
            stored_hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                    self.hgvs_genomic_5pr.ac,
                    self.hgvs_genomic_5pr.type,
                    int(vcf_dict['pos']),vcf_dict['ref'],vcf_dict['alt'])

            # Detect intronic variation using normalization
            intronic_variant = 'false'

            # Save a copy of current hgvs_coding
            try:
                saved_hgvs_coding = self.variant.no_norm_evm.g_to_t(stored_hgvs_not_delins,
                                                                    hgvs_coding.ac)
            except vvhgvs.exceptions.HGVSInvalidIntervalError as e:
                if str(e) == 'start or end or both are beyond the bounds of transcript record':
                    continue
                else:
                    saved_hgvs_coding = self.variant.no_norm_evm.g_to_n(reverse_normalized_hgvs_genomic,
                                                                        hgvs_coding.ac)

            # Look for normalized variant options that do not match hgvs_coding
            hgvs_seek_var = self.get_hgvs_seek_var(hgvs_genomic, hgvs_coding)

            if (
                    hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                    hgvs_coding.posedit.pos.start.base + hgvs_coding.posedit.pos.start.offset) and (
                    hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                    hgvs_coding.posedit.pos.end.base + hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
                pass
            else:
                hgvs_seek_var = saved_hgvs_coding

            try:
                self.variant.hn.normalize(hgvs_seek_var)
            except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                error = str(e)
                if 'Normalization of intronic variants is not supported' in error or \
                        'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                    if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        intronic_variant = 'hard_fail'
                    else:
                        # Double check to see whether the variant is actually intronic?
                        for exon in ori:
                            genomic_start = int(exon['alt_start_i'])
                            genomic_end = int(exon['alt_end_i'])
                            if (genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                                    genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                                intronic_variant = 'false'
                                break
                            else:
                                intronic_variant = 'true'

            if intronic_variant != 'hard_fail':
                if _has_any_offset(hgvs_seek_var):
                    # Double check to see whether the variant is actually intronic?
                    for exon in ori:
                        genomic_start = int(exon['alt_start_i'])
                        genomic_end = int(exon['alt_end_i'])
                        if (genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                                genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                            intronic_variant = 'false'
                            break
                        else:
                            intronic_variant = 'true'

            if re.search(r'\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search(r'\d+-', str(
                    hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+\+', str(
                    hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+-', str(hgvs_seek_var.posedit.pos)):
                # Double check to see whether the variant is actually intronic?
                for exon in ori:
                    genomic_start = int(exon['alt_start_i'])
                    genomic_end = int(exon['alt_end_i'])
                    if (genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                            genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                        intronic_variant = 'false'
                        break
                    else:
                        intronic_variant = 'true'

            if intronic_variant != 'true':
                # Flag RefSeqGene for amendment
                # amend_RefSeqGene = 'false'
                # Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping
                # transcript lengths
                hgvs_not_delins = ''
                if stored_hgvs_not_delins != '':
                    # Refresh hgvs_not_delins from stored_hgvs_not_delins
                    hgvs_not_delins = self.dup_ins_5prime_shift(stored_hgvs_not_delins, saved_hgvs_coding)

                    try:
                        self.tx_hgvs_not_delins = self.variant.no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                  saved_hgvs_coding.ac)
                    except vvhgvs.exceptions.HGVSInvalidIntervalError:
                        self.tx_hgvs_not_delins = self.variant.no_norm_evm.g_to_n(reverse_normalized_hgvs_genomic,
                                                                                  saved_hgvs_coding.ac)
                    # Create normalized version of tx_hgvs_not_delins
                    rn_tx_hgvs_not_delins = copy.deepcopy(self.tx_hgvs_not_delins)

                    # Check for +1 base and adjust
                    if '+' in str(rn_tx_hgvs_not_delins.posedit.pos.end) and '+' in str(
                            rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins = self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins)

                    elif '+' in str(rn_tx_hgvs_not_delins.posedit.pos.end):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_end_base_to_next_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding, back=False)

                    elif '+' in str(rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_start_base_to_previous_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding)

                    # Check for -ve base and adjust
                    elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.end) and '-' in str(
                            rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins = self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins)

                    elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.end):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_end_base_to_next_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding)

                    elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_start_base_to_previous_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding, with_base_subtract=True)

                    # Logic
                    hgvs_not_delins = self.logic_check(hgvs_not_delins, rn_tx_hgvs_not_delins, hgvs_coding)

                # 'At hgvs_genomic'
                # Final sanity checks
                try:
                    self.validator.vm.g_to_t(hgvs_not_delins, self.tx_hgvs_not_delins.ac,
                                             alt_aln_method=self.validator.alt_aln_method)
                except Exception as e:
                    if str(e) == 'start or end or both are beyond the bounds of transcript record':
                        continue
                try:
                    self.variant.hn.normalize(self.tx_hgvs_not_delins)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                    error = str(e)
                    if 'Normalization of intronic variants is not supported' in error or \
                            'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                            continue
                        elif 'Normalization of intronic variants is not supported' in error:
                            # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                            self.disparity_deletion_in = ['transcript', 'Requires Analysis']

                # Recreate hgvs_genomic
                if self.disparity_deletion_in[0] == 'transcript':
                    hgvs_genomic = hgvs_not_delins

                # Find oddly placed gaps where the tx variant is encompassed in the gap
                if self.disparity_deletion_in[0] == 'false' and (possibility_counter == 3 or possibility_counter == 4):
                    rg = self.variant.reverse_normalizer.normalize(hgvs_not_delins)
                    rtx = self.validator.vm.g_to_t(rg, self.tx_hgvs_not_delins.ac,
                                                   alt_aln_method=self.validator.alt_aln_method)
                    fg = self.variant.hn.normalize(hgvs_not_delins)
                    ftx = self.validator.vm.g_to_t(fg, self.tx_hgvs_not_delins.ac,
                                                   alt_aln_method=self.validator.alt_aln_method)
                    if (rtx.posedit.pos.start.offset == 0 and rtx.posedit.pos.end.offset == 0) and (
                            ftx.posedit.pos.start.offset != 0 and ftx.posedit.pos.end.offset != 0):
                        exons = self.variant.map_dat.mapped_exons(
                                ftx.ac, hgvs_not_delins.ac,
                                alt_aln_method=self.validator.alt_aln_method)
                        exonic = False
                        for ex_test in exons:
                            if ftx.posedit.pos.start.base in range(ex_test[6], ex_test[
                                    7]) and ftx.posedit.pos.end.base in range(ex_test[6], ex_test[7]):
                                exonic = True
                        if exonic is True:
                            hgvs_not_delins = fg
                            hgvs_genomic = fg
                            self.hgvs_genomic_5pr = fg
                            try:
                                self.tx_hgvs_not_delins = self.validator.vm.c_to_n(ftx)
                            except Exception:
                                self.tx_hgvs_not_delins = ftx
                            self.disparity_deletion_in = ['transcript', 'Requires Analysis']

                # Pre-processing of self.tx_hgvs_not_delins
                try:
                    if self.tx_hgvs_not_delins.posedit.edit.alt is None:
                        self.tx_hgvs_not_delins.posedit.edit.alt = ''
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        self.tx_hgvs_not_delins = hgvs_dup_to_delins(self.tx_hgvs_not_delins)


                # Has a hard set variant been identified from pushes?
                hard_set_outputs = False
                if disparity_info != ['false', 'false'] and len(disparity_info) == 4:
                    self.tx_hgvs_not_delins = disparity_info[2]
                    self.disparity_deletion_in = [disparity_info[0], disparity_info[1]]
                    hgvs_refreshed_variant = hgvs_coding
                    hgvs_genomic = possibility
                    suppress_c_normalization = 'true'
                    hard_set_outputs = True

                # GAP IN THE TRANSCRIPT DISPARITY DETECTED
                elif self.disparity_deletion_in[0] == 'transcript':
                    # Suppress intron boundary crossing due to non-intron intron based c. seq annotations
                    suppress_c_normalization = 'true'
                    # amend_RefSeqGene = 'true'
                    # ANY VARIANT WHOLLY WITHIN THE GAP
                    hgvs_refreshed_variant = self.transcript_disparity(
                        reverse_normalized_hgvs_genomic,
                        hgvs_genomic,
                        2,
                    )

                # GAP IN THE CHROMOSOME
                elif self.disparity_deletion_in[0] == 'chromosome':
                    suppress_c_normalization = 'true'
                    if possibility_counter == 3:
                        hgvs_refreshed_variant = stash_tx_right
                    elif possibility_counter == 4:
                        hgvs_refreshed_variant = stash_tx_left
                    else:
                        hgvs_refreshed_variant = chromosome_normalized_hgvs_coding

                else:
                    # Keep the same by re-setting rel_var
                    hgvs_refreshed_variant = hgvs_coding
                # amend_RefSeqGene = 'false'

                # Edit the output
                if 'NM_' in str(hgvs_refreshed_variant.ac) and not 'c' in str(hgvs_refreshed_variant.type):
                    hgvs_refreshed_variant = self.variant.no_norm_evm.n_to_c(hgvs_refreshed_variant)

                try:
                    self.variant.hn.normalize(hgvs_refreshed_variant)
                except Exception as e:
                    error = str(e)

                    # Ensure the final variant is not intronic nor does it cross exon boundaries
                    if 'Normalization of intronic variants is not supported' in error or \
                            'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        hgvs_refreshed_variant = saved_hgvs_coding
                    else:
                        pass
                        continue

                # Quick check to make sure the coding variant has not changed
                try:
                    to_test = self.variant.hn.normalize(hgvs_refreshed_variant)
                except:
                    to_test = hgvs_refreshed_variant
                if not _same_edit(to_test, hgvs_coding):
                    # Try the next available genomic option
                    if hgvs_coding.posedit.edit.type == 'identity' and to_test.posedit.edit.type == 'identity':
                        hgvs_coding = to_test
                    elif hard_set_outputs is True:
                        hgvs_coding = to_test
                    else:
                        continue

                if hard_set_outputs is not True:
                    # Update hgvs_genomic
                    hgvs_genomic = self.validator.myvm_t_to_g(hgvs_refreshed_variant, hgvs_genomic.ac,
                                                              self.variant.no_norm_evm, self.variant.hn,
                                                              self.variant.map_dat)
                    if hgvs_genomic.posedit.edit.type == 'identity':
                        re_c = self.validator.vm.g_to_t(hgvs_genomic, hgvs_refreshed_variant.ac,
                                                        alt_aln_method=self.validator.alt_aln_method)
                        if (self.variant.hn.normalize(re_c)) != (self.variant.hn.normalize(hgvs_refreshed_variant)):
                            shuffle_left_g = copy.copy(hgvs_genomic)
                            shuffle_left_g.posedit.edit.ref = ''
                            shuffle_left_g.posedit.edit.alt = ''
                            shuffle_left_g.posedit.pos.start.base = shuffle_left_g.posedit.pos.start.base - 1
                            shuffle_left_g.posedit.pos.end.base = shuffle_left_g.posedit.pos.end.base - 1
                            shuffle_left_g = self.variant.reverse_normalizer.normalize(shuffle_left_g)
                            re_c = self.validator.vm.g_to_t(shuffle_left_g, hgvs_refreshed_variant.ac,
                                                            alt_aln_method=self.validator.alt_aln_method)
                            if (self.variant.hn.normalize(re_c)) != (self.variant.hn.normalize(hgvs_refreshed_variant)):
                                hgvs_genomic = shuffle_left_g

            # Break if gap has been detected
            if self.disparity_deletion_in[0] != 'false':
                break

        # Normailse hgvs_genomic
        try:
            hgvs_genomic = self.variant.hn.normalize(hgvs_genomic)
        except vvhgvs.exceptions.HGVSError as e:
            # Strange error caused by gap in genomic

            if 'base start position must be <= end position' in error and self.disparity_deletion_in[0] == 'chromosome':
                if hgvs_genomic.posedit.edit.type == 'delins':
                    start = hgvs_genomic.posedit.pos.start.base
                    end = hgvs_genomic.posedit.pos.end.base
                    lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                    rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                    hgvs_genomic.posedit.edit.ref = lhb + rhb
                    hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                    hgvs_genomic.posedit.pos.start.base = end
                    hgvs_genomic.posedit.pos.end.base = start
                    hgvs_genomic = self.variant.hn.normalize(hgvs_genomic)
                if hgvs_genomic.posedit.edit.type == 'del':
                    start = hgvs_genomic.posedit.pos.start.base
                    end = hgvs_genomic.posedit.pos.end.base
                    lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                    rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                    hgvs_genomic.posedit.edit.ref = lhb + rhb
                    hgvs_genomic.posedit.edit.alt = lhb + rhb
                    hgvs_genomic.posedit.pos.start.base = end
                    hgvs_genomic.posedit.pos.end.base = start
                    hgvs_genomic = self.variant.hn.normalize(hgvs_genomic)

        pass

        return hgvs_genomic, suppress_c_normalization, hgvs_coding

    def g_to_t_gapped_mapping_stage2(self, ori, hgvs_coding, hgvs_genomic):
        pass

        pass

        hgvs_genomic_variant = hgvs_genomic
        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic_variant)
        self.hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)
        vcf_dict = hgvs_utils.hgvs2vcf(reverse_normalized_hgvs_genomic, self.variant.primary_assembly,
                                       self.variant.reverse_normalizer, self.validator.sf)

        stored_hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                self.hgvs_genomic_5pr.ac,
                self.hgvs_genomic_5pr.type,
                int(vcf_dict['pos']),vcf_dict['ref'],vcf_dict['alt'])
        self.orientation = int(ori[0]['alt_strand'])

        saved_hgvs_coding = copy.deepcopy(hgvs_coding)

        # is it in an exon?
        is_it_in_an_exon = 'no'
        for exon in ori:
            genomic_start = int(exon['alt_start_i'])
            genomic_end = int(exon['alt_end_i'])
            # Take from stored copy
            if (genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                    genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                is_it_in_an_exon = 'yes'
        if is_it_in_an_exon == 'yes':
            # map form reverse normalized g. to c.
            # Attempt to find gaps in reference sequence by catching disparity in genome length and
            # overlapping transcript lengths
            self.disparity_deletion_in = ['false', 'false']
            hgvs_not_delins = ''
            hard_fail = 'false'
            if stored_hgvs_not_delins != '':
                # Refresh hgvs_not_delins from stored_hgvs_not_delins
                hgvs_not_delins = self.dup_ins_5prime_shift(stored_hgvs_not_delins, saved_hgvs_coding)
                try:
                    self.tx_hgvs_not_delins = self.variant.no_norm_evm.g_to_n(hgvs_not_delins,
                                                                              saved_hgvs_coding.ac)
                except Exception as e:
                    if str(e) == 'start or end or both are beyond the bounds of transcript record':
                        self.tx_hgvs_not_delins = hgvs_coding
                        hard_fail = 'true'

                # Create normalized version of self.tx_hgvs_not_delins
                rn_tx_hgvs_not_delins = copy.deepcopy(self.tx_hgvs_not_delins)
                # Check for +ve base and adjust
                if '+' in str(rn_tx_hgvs_not_delins.posedit.pos.end) and '+' in \
                        str(rn_tx_hgvs_not_delins.posedit.pos.start):
                    rn_tx_hgvs_not_delins = self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins)

                elif '+' in str(rn_tx_hgvs_not_delins.posedit.pos.end):
                    rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_end_base_to_next_nonoffset(
                        rn_tx_hgvs_not_delins, saved_hgvs_coding, back=False)

                elif '+' in str(rn_tx_hgvs_not_delins.posedit.pos.start):
                    rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_start_base_to_previous_nonoffset(
                        rn_tx_hgvs_not_delins, saved_hgvs_coding)

                # Check for -ve base and adjust
                elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.end) and '-' in str(
                        rn_tx_hgvs_not_delins.posedit.pos.start):
                    rn_tx_hgvs_not_delins = self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins)

                elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.end):
                    rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_end_base_to_next_nonoffset(
                        rn_tx_hgvs_not_delins, saved_hgvs_coding)

                elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.start):
                    rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_start_base_to_previous_nonoffset(
                        rn_tx_hgvs_not_delins, saved_hgvs_coding, with_base_subtract=True)

                # Logic
                hgvs_not_delins = self.logic_check(hgvs_not_delins, rn_tx_hgvs_not_delins, hgvs_coding)

            # Final sanity checks
            try:
                self.validator.vm.g_to_t(hgvs_not_delins, self.tx_hgvs_not_delins.ac,
                                         alt_aln_method=self.validator.alt_aln_method)
            except Exception as e:
                if str(e) == 'start or end or both are beyond the bounds of transcript record':
                    pass
                    hard_fail == 'true'
            try:
                self.variant.hn.normalize(self.tx_hgvs_not_delins)
            except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                error = str(e)
                if 'Normalization of intronic variants is not supported' in error or \
                        'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                    if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        pass
                        hard_fail = 'true'
                    elif 'Normalization of intronic variants is not supported' in error:
                        # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                        self.disparity_deletion_in = ['transcript', 'Requires Analysis']

            if hard_fail == 'true':
                self.disparity_deletion_in = ['false', 'false']

            # Recreate hgvs_genomic
            if self.disparity_deletion_in[0] == 'transcript':
                hgvs_genomic = hgvs_not_delins

            # Pre-processing of tx_hgvs_not_delins
            try:
                if self.tx_hgvs_not_delins.posedit.edit.alt is None:
                    self.tx_hgvs_not_delins.posedit.edit.alt = ''
            except Exception as e:
                if str(e) == "'Dup' object has no attribute 'alt'":
                    self.tx_hgvs_not_delins = hgvs_dup_to_delins(self.tx_hgvs_not_delins)

            # GAP IN THE TRANSCRIPT DISPARITY DETECTED
            if self.disparity_deletion_in[0] == 'transcript':
                # ANY VARIANT WHOLLY WITHIN THE GAP
                hgvs_refreshed_variant = self.transcript_disparity(
                    reverse_normalized_hgvs_genomic,
                    hgvs_genomic,
                    3,
                )

                # Look for missed duplications into direct complete in-line repeats
                try:
                    check_refreshed_variant = (self.validator.vm.
                                              n_to_c(self.tx_hgvs_not_delins))
                except vvhgvs.exceptions.HGVSError:
                    check_refreshed_variant = self.tx_hgvs_not_delins
                pass

                try:
                    if (hasattr(check_refreshed_variant.posedit.pos.start, "offset")
                            or hasattr(check_refreshed_variant.posedit.pos.end, "offset")):
                        pass
                        if (hgvs_refreshed_variant.posedit.pos.end.base == check_refreshed_variant.posedit.pos.end.base
                            and (int(check_refreshed_variant.posedit.pos.start.base)) ==
                                (int(hgvs_refreshed_variant.posedit.pos.start.base)+1)
                                and hgvs_refreshed_variant.posedit.edit.alt[2:] ==
                                check_refreshed_variant.posedit.edit.alt):
                                    pass
                                    hgvs_refreshed_variant = check_refreshed_variant
                except Exception:
                    pass

            # GAP IN THE CHROMOSOME
            elif self.disparity_deletion_in[0] == 'chromosome':
                hgvs_refreshed_variant = self.tx_hgvs_not_delins

            else:
                # Keep the same by re-setting rel_var
                hgvs_refreshed_variant = saved_hgvs_coding

            # Edit the output
            hgvs_refreshed_variant = self.edit_output(hgvs_refreshed_variant, saved_hgvs_coding)

            # Sort out equality to equality c. events where the code will add 2 additional bases
            if hgvs_coding.posedit.edit.type == 'identity' and hgvs_refreshed_variant.posedit.edit.type == 'identity':
                pass
            else:
                hgvs_coding = copy.deepcopy(hgvs_refreshed_variant)

        pass
        return hgvs_coding

    def g_to_t_gap_compensation_version3(self, hgvs_alt_genomic, hgvs_coding, ori, alt_chr, rec_var):

        pass

        self.orientation = int(ori[0]['alt_strand'])
        hgvs_genomic = copy.deepcopy(hgvs_alt_genomic)

        pass

        rn_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_alt_genomic)
        self.hgvs_genomic_possibilities.append([rn_hgvs_genomic, ['false', 'false']])
        if self.orientation != -1:
            try:
                chromosome_normalized_hgvs_coding = self.variant.reverse_normalizer.normalize(
                    hgvs_coding)
            except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                chromosome_normalized_hgvs_coding = hgvs_coding
        else:
            try:
                chromosome_normalized_hgvs_coding = self.variant.hn.normalize(hgvs_coding)
            except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                chromosome_normalized_hgvs_coding = hgvs_coding

        most_3pr_hgvs_genomic = self.validator.myvm_t_to_g(chromosome_normalized_hgvs_coding, alt_chr,
                                                           self.variant.no_norm_evm, self.variant.hn,
                                                           self.variant.map_dat)
        self.hgvs_genomic_possibilities.append([most_3pr_hgvs_genomic, ['false', 'false']])

        # First to the right
        hgvs_stash = copy.deepcopy(hgvs_coding)
        stash_tx_right = ''
        stash_tx_left = ''
        self._ensure_map_data_provider()

        # Capture instances where variant merging hard-sets the outputs
        map_fail = False
        try:
            if hgvs_stash.type == 'c':
                hgvs_stash = self.validator.vm.c_to_n(hgvs_stash)
            hgvs_right = copy.copy(hgvs_stash)
            hgvs_left = copy.copy(hgvs_right)
            if not (
                    getattr(hgvs_right.posedit.pos.start,'offset',False) or
                    getattr(hgvs_right.posedit.pos.end,'offset',False)):
                hgvs_right = self.variant.hn.normalize(hgvs_right)
                hgvs_left = self.variant.reverse_normalizer.normalize(hgvs_left)
            #if hgvs_right.type == 'n':
            #    non_variant_genomic_ac = hgvs_genomic_variant.ac
            if hgvs_right == hgvs_left and hgvs_right.type == 'n':
                genomic_mapping = self.validator.vm.n_to_g(hgvs_right, hgvs_alt_genomic.ac,alt_aln_method=self.validator.alt_aln_method)
            else:
                genomic_mapping = False
        except (vvhgvs.exceptions.HGVSUnsupportedOperationError,
                vvhgvs.exceptions.HGVSInvalidVariantError,
                vvhgvs.exceptions.HGVSUsageError,
                vvhgvs.exceptions.HGVSDataNotAvailableError) as err1:
            map_fail = True
            pass

        try:
            if map_fail:
                raise ValueError("Already failed in n->g mapping")
            stash_ac = hgvs_stash.ac
            stash_dict = hgvs_utils.hard_right_hgvs2vcf(hgvs_stash,
                                                        self.variant.primary_assembly,
                                                        self.variant.hn,
                                                        self.variant.reverse_normalizer,
                                                        self.validator.sf,
                                                        hgvs_coding.ac,
                                                        self.variant.map_dat,
                                                        self.validator.alt_aln_method,
                                                        self.validator.hp,
                                                        self.validator.vm,
                                                        self.validator.merge_hgvs_3pr,
                                                        genomic_ac=hgvs_alt_genomic.ac,
                                                        mapped_g=copy.copy(genomic_mapping),
                                                        pre_norm=hgvs_right)

            # make a not real deletion insertion
            stash_hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                    stash_ac,
                    hgvs_stash.type,
                    int(stash_dict['pos']),stash_dict['ref'],stash_dict['alt'],
                    offset_pos=True)
            try:
                stash_hgvs_not_delins = self.variant.no_norm_evm.n_to_c(stash_hgvs_not_delins)
            except Exception as e:
                pass

            # Store a tx copy for later use
            test_stash_tx_right = copy.deepcopy(stash_hgvs_not_delins)
            stash_genomic = self.validator.myvm_t_to_g(test_stash_tx_right, hgvs_alt_genomic.ac,
                                                       self.variant.no_norm_evm, self.variant.hn,
                                                       self.variant.map_dat)

            if len(test_stash_tx_right.posedit.edit.ref) == ((stash_genomic.posedit.pos.end.base -
                                                              stash_genomic.posedit.pos.start.base) + 1):
                stash_tx_right = test_stash_tx_right
                alt = getattr(test_stash_tx_right.posedit.edit, 'alt', False)
                if not alt:
                    alt = ''
                g_alt = getattr(stash_genomic.posedit.edit, 'alt', False)
                if not g_alt:
                    g_alt = ''
                if (len(alt) - (
                        test_stash_tx_right.posedit.pos.end.base - test_stash_tx_right.posedit.pos.start.base) + 1) != (
                        len(g_alt) - (
                        stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])
                else:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
            elif test_stash_tx_right.posedit.edit.type == 'identity':
                hgvs_reform_ident = hgvs_delins_parts_to_hgvs_obj(
                        test_stash_tx_right.ac,
                        'c',
                        test_stash_tx_right.posedit.pos,
                        test_stash_tx_right.posedit.edit.ref,
                        '')
                try:
                    self.variant.hn.normalize(hgvs_reform_ident)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'spanning the exon-intron boundary' in error:
                        stash_tx_right = test_stash_tx_right
                        self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_right = test_stash_tx_right
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])
            else:
                try:
                    self.variant.hn.normalize(test_stash_tx_right)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_right = test_stash_tx_right
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])

            # IDENTIFY GAP AND HARD SET

            if stash_dict['needs_a_push'] is True or stash_dict['identifying_g_variant'] is not False:

                # Look for merged variant from hard push
                if stash_dict['merged_variant'] is not False:
                    merged_variant = stash_dict['merged_variant']
                    identifying_variant = stash_dict['identifying_variant']
                    try:
                        merged_variant = self.validator.vm.n_to_c(merged_variant)
                    except (TypeError, vvhgvs.exceptions.HGVSInvalidVariantError):
                        pass
                    try:
                        identifying_variant = self.validator.vm.n_to_c(identifying_variant)
                    except (TypeError, vvhgvs.exceptions.HGVSInvalidVariantError):
                        pass

                    stash_genomic = self.validator.myvm_t_to_g(identifying_variant, stash_genomic.ac,
                                                               self.variant.no_norm_evm,self.variant.hn,
                                                               self.variant.map_dat)
                    stash_hgvs_not_delins = merged_variant
                    test_stash_tx_right = merged_variant

                elif stash_dict['identifying_g_variant'] is not False:
                    stash_hgvs_not_delins = stash_dict['identifying_g_variant']
                    stash_genomic = stash_dict['identifying_g_variant']

                # Look for gap info
                normalized_stash_genomic = self.variant.hn.normalize(stash_genomic)
                stash_tx_right = test_stash_tx_right
                if stash_hgvs_not_delins.posedit.edit.type == "ins":
                    len_tx = 2
                else:
                    len_tx = len(stash_hgvs_not_delins.posedit.edit.ref)
                if stash_genomic.posedit.edit.type == "ins":
                    len_gen = 2
                else:
                    len_gen = len(stash_genomic.posedit.edit.ref)
                if len_tx > len_gen:
                    gap_in = 'chromosome'
                    gap_len = len_tx - len_gen
                else:
                    gap_in = 'transcript'
                    gap_len = len_gen - len_tx

                # Set the options to a single option based on the results of pushing
                self.hgvs_genomic_possibilities = [[normalized_stash_genomic, [gap_in,
                                                                               gap_len,
                                                                               stash_hgvs_not_delins,
                                                                               stash_genomic]]]

        except (vvhgvs.exceptions.HGVSError, ValueError) as e:
            pass

        # Then to the left
        try:
            if map_fail:
                raise ValueError("Already failed in n->g mapping")
            stash_ac = hgvs_stash.ac
            stash_dict = hgvs_utils.hard_left_hgvs2vcf(hgvs_stash,
                                                       self.variant.primary_assembly,
                                                       self.variant.hn,
                                                       self.variant.reverse_normalizer,
                                                       self.validator.sf,
                                                       hgvs_coding.ac,
                                                       self.variant.map_dat,
                                                       self.validator.alt_aln_method,
                                                       self.validator.hp,
                                                       self.validator.vm,
                                                       self.validator.merge_hgvs_3pr,
                                                       genomic_ac=hgvs_alt_genomic.ac,
                                                       mapped_g=copy.copy(genomic_mapping),
                                                       pre_norm=hgvs_left)

            # make a not real deletion insertion
            stash_hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                    stash_ac,
                    hgvs_stash.type,
                    int(stash_dict['pos']),stash_dict['ref'],stash_dict['alt'],
                    offset_pos=True)
            try:
                stash_hgvs_not_delins = self.variant.no_norm_evm.n_to_c(stash_hgvs_not_delins)
            except Exception as e:
                pass
                # Store a tx copy for later use
            test_stash_tx_left = copy.deepcopy(stash_hgvs_not_delins)
            stash_genomic = self.validator.myvm_t_to_g(test_stash_tx_left, hgvs_alt_genomic.ac,
                                                       self.variant.no_norm_evm, self.variant.hn,
                                                       self.variant.map_dat)

            if len(test_stash_tx_left.posedit.edit.ref) == ((stash_genomic.posedit.pos.end.base -
                                                             stash_genomic.posedit.pos.start.base) + 1):
                stash_tx_left = test_stash_tx_left
                alt = getattr(test_stash_tx_left.posedit.edit, 'alt', False)
                if not alt:
                    alt = ''
                g_alt = getattr(stash_genomic.posedit.edit, 'alt', False)
                if not g_alt:
                    g_alt = ''
                if (len(alt) - (
                        test_stash_tx_left.posedit.pos.end.base - test_stash_tx_left.posedit.pos.start.base) + 1) != (
                        len(g_alt) - (
                        stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])
                else:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
            elif test_stash_tx_left.posedit.edit.type == 'identity':
                hgvs_reform_ident = hgvs_delins_parts_to_hgvs_obj(
                        test_stash_tx_left.ac,
                        'c',
                        test_stash_tx_left.posedit.pos,
                        test_stash_tx_left.posedit.edit.ref,
                        '')
                try:
                    self.variant.hn.normalize(hgvs_reform_ident)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'spanning the exon-intron boundary' in error:
                        stash_tx_left = test_stash_tx_left
                        self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_left = test_stash_tx_left
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])

            else:
                try:
                    self.variant.hn.normalize(test_stash_tx_left)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_left = test_stash_tx_left
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])

            # IDENTIFY GAP AND HARD SET
            if stash_dict['needs_a_push'] is True or stash_dict['identifying_g_variant'] is not False:

                # Look for merged variant from hard push
                if stash_dict['merged_variant'] is not False:
                    merged_variant = stash_dict['merged_variant']
                    identifying_variant = stash_dict['identifying_variant']
                    try:
                        merged_variant = self.validator.vm.n_to_c(merged_variant)
                    except (TypeError, vvhgvs.exceptions.HGVSInvalidVariantError):
                        pass
                    try:
                        identifying_variant = self.validator.vm.n_to_c(identifying_variant)
                    except (TypeError, vvhgvs.exceptions.HGVSInvalidVariantError):
                        pass

                    stash_genomic = self.validator.myvm_t_to_g(identifying_variant, stash_genomic.ac,
                                                               self.variant.no_norm_evm,self.variant.hn,
                                                               self.variant.map_dat)
                    stash_hgvs_not_delins = merged_variant
                    test_stash_tx_left = merged_variant

                elif stash_dict['identifying_g_variant'] is not False:
                    stash_hgvs_not_delins = stash_dict['identifying_g_variant']
                    stash_genomic = stash_dict['identifying_g_variant']

                # Look for gap info
                normalized_stash_genomic = self.variant.hn.normalize(stash_genomic)
                stash_tx_left = test_stash_tx_left
                if stash_hgvs_not_delins.posedit.edit.type == "ins":
                    len_tx = 2
                else:
                    len_tx = len(stash_hgvs_not_delins.posedit.edit.ref)
                if stash_genomic.posedit.edit.type == "ins":
                    len_gen = 2
                else:
                    len_gen = len(stash_genomic.posedit.edit.ref)
                if len_tx > len_gen:
                    gap_in = 'chromosome'
                    gap_len = len_tx - len_gen
                else:
                    gap_in = 'transcript'
                    gap_len = len_gen - len_tx

                if stash_dict["pre_merged_variant"] is not False and stash_dict["identifying_g_variant"] is False:
                    if ("transcript" in gap_in
                            and len_tx == len_gen
                            and stash_dict["pre_merged_variant"].type == "g")\
                            and normalized_stash_genomic.posedit.edit.type == "sub"\
                            and stash_dict["pre_merged_variant"].posedit.edit.type == "delins"\
                            and hgvs_coding.posedit.edit.type == "sub":
                        normalized_stash_genomic = stash_dict["pre_merged_variant"]

                # Set the options to a single option based on the results of pushing
                self.hgvs_genomic_possibilities = [[normalized_stash_genomic, [gap_in,
                                                                               gap_len,
                                                                               stash_hgvs_not_delins,
                                                                               stash_genomic]]]

        except (vvhgvs.exceptions.HGVSError, ValueError) as e:
            pass

        # direct mapping from reverse_normalized transcript insertions in the delins format
        self.rev_norm_ins(hgvs_coding, hgvs_genomic)

        # Set variables for problem specific warnings
        self.gapped_transcripts = ''
        self.auto_info = ''

        # Mark as not disparity detected
        self.disparity_deletion_in = ['false', 'false']
        # Loop through to see if a gap can be located
        possibility_counter = 0

        # If hard pushing identified a variant (in theory there can be only 1) then remove all other possibilities
        hard_possibility = []
        for hard_set_check in self.hgvs_genomic_possibilities:
            if hard_set_check[1] != ['false', 'false']:
                hard_possibility.append(hard_set_check)
        if len(hard_possibility) >= 1:
            self.hgvs_genomic_possibilities = hard_possibility

        for a_possibility in self.hgvs_genomic_possibilities:
            possibility = a_possibility[0]
            disparity_info = a_possibility[1]
            possibility_counter = possibility_counter + 1
            # Loop out stash possibilities which will not spot gaps so are empty
            if possibility == '':
                continue

            # Use VCF generation code to push hgvs_genomic as for 5 prime as possible to uncover gaps
            hgvs_genomic_variant = possibility
            reverse_normalized_hgvs_genomic = ''

            # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
            try:
                reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(
                    hgvs_genomic_variant)
            except vvhgvs.exceptions.HGVSError as e:
                # Strange error caused by gap in genomic
                error = str(e)
                if 'base start position must be <= end position' in error:
                    if hgvs_genomic.posedit.edit.type == 'delins':
                        start = hgvs_genomic.posedit.pos.start.base
                        end = hgvs_genomic.posedit.pos.end.base
                        lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                        rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                        hgvs_genomic.posedit.edit.ref = lhb + rhb
                        hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                        hgvs_genomic.posedit.pos.start.base = end
                        hgvs_genomic.posedit.pos.end.base = start
                        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(
                            hgvs_genomic)
                    if hgvs_genomic.posedit.edit.type == 'del':
                        start = hgvs_genomic.posedit.pos.start.base
                        end = hgvs_genomic.posedit.pos.end.base
                        lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                        rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                        hgvs_genomic.posedit.edit.ref = lhb + rhb
                        hgvs_genomic.posedit.edit.alt = lhb + rhb
                        hgvs_genomic.posedit.pos.start.base = end
                        hgvs_genomic.posedit.pos.end.base = start
                        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(
                            hgvs_genomic)
                if 'insertion length must be 1' in error:
                    if hgvs_genomic.posedit.edit.type == 'ins':
                        start = hgvs_genomic.posedit.pos.start.base
                        end = hgvs_genomic.posedit.pos.end.base
                        ref_bases = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, end)
                        lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                        rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start, end)
                        hgvs_genomic.posedit.edit.ref = lhb + rhb
                        hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(
                            hgvs_genomic)

            self.hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)
            # Store a copy for later use

            # Make VCF
            vcf_dict = hgvs_utils.hgvs2vcf(reverse_normalized_hgvs_genomic, self.variant.primary_assembly,
                                           self.variant.reverse_normalizer, self.validator.sf)
            # Look for exonic gaps within transcript or chromosome
            # Store a not real deletion insertion to test for gapping
            stored_hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                    self.hgvs_genomic_5pr.ac,
                    self.hgvs_genomic_5pr.type,
                    int(vcf_dict['pos']),vcf_dict['ref'],vcf_dict['alt'])
            # Save a copy of current hgvs_coding
            saved_hgvs_coding = ''
            try:
                saved_hgvs_coding = self.variant.no_norm_evm.g_to_t(stored_hgvs_not_delins,
                                                                    hgvs_coding.ac)
            except Exception as e:
                if str(e) == 'start or end or both are beyond the bounds of transcript record':
                    continue

            # Detect intronic variation using normalization
            intronic_variant = 'false'
            # Look for normalized variant options that do not match hgvs_coding
            hgvs_seek_var = self.get_hgvs_seek_var(hgvs_genomic, hgvs_coding)
            if (
                    hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                    hgvs_coding.posedit.pos.start.base + hgvs_coding.posedit.pos.start.offset) and (
                    hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                    hgvs_coding.posedit.pos.end.base + hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
                pass
            else:
                hgvs_seek_var = saved_hgvs_coding

            try:
                self.variant.hn.normalize(hgvs_seek_var)
            except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                error = str(e)
                if 'Normalization of intronic variants is not supported' in error or \
                        'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                    if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        intronic_variant = 'hard_fail'
                    else:
                        # Double check to see whether the variant is actually intronic?
                        for exon in ori:
                            genomic_start = int(exon['alt_start_i'])
                            genomic_end = int(exon['alt_end_i'])
                            if (genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                                    genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                                intronic_variant = 'false'
                                break
                            else:
                                intronic_variant = 'true'

            if intronic_variant != 'hard_fail':
                if _has_any_offset(hgvs_seek_var):
                    # Double check to see whether the variant is actually intronic?
                    for exon in ori:
                        genomic_start = int(exon['alt_start_i'])
                        genomic_end = int(exon['alt_end_i'])
                        if (genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                                genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                            intronic_variant = 'false'
                            break
                        else:
                            intronic_variant = 'true'

            if intronic_variant != 'true':
                # Flag RefSeqGene for ammendment
                # amend_RefSeqGene = 'false'
                # Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping
                # transcript lengths
                hgvs_not_delins = ''
                if stored_hgvs_not_delins != '':
                    # Refresh hgvs_not_delins from stored_hgvs_not_delins
                    hgvs_not_delins = self.dup_ins_5prime_shift(stored_hgvs_not_delins, saved_hgvs_coding)

                    self.tx_hgvs_not_delins = self.variant.no_norm_evm.g_to_n(hgvs_not_delins,
                                                                              saved_hgvs_coding.ac)
                    # Create normalized version of tx_hgvs_not_delins
                    rn_tx_hgvs_not_delins = copy.deepcopy(self.tx_hgvs_not_delins)
                    # Check for +1 base and adjust
                    if '+' in str(rn_tx_hgvs_not_delins.posedit.pos.end) and '+' in str(
                            rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins = self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins)

                    elif '+' in str(rn_tx_hgvs_not_delins.posedit.pos.end):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_end_base_to_next_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding, back=False)

                    elif '+' in str(rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_start_base_to_previous_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding)

                    # Check for -ve base and adjust
                    elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.end) and '-' in str(
                            rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins = self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins)

                    elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.end):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_end_base_to_next_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding)

                    elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_start_base_to_previous_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding, with_base_subtract=True)

                    # Logic

                    hgvs_not_delins = self.logic_check(hgvs_not_delins, rn_tx_hgvs_not_delins, hgvs_coding,
                                                       do_continue=True, offset_check=True)

                # Final sanity checks
                try:
                    self.validator.vm.g_to_t(hgvs_not_delins,
                                             self.tx_hgvs_not_delins.ac,
                                             alt_aln_method=self.validator.alt_aln_method)
                except Exception as e:
                    if str(e) == 'start or end or both are beyond the bounds of transcript record':
                        continue
                try:
                    self.variant.hn.normalize(self.tx_hgvs_not_delins)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                    error = str(e)
                    if 'Normalization of intronic variants is not supported' in error or \
                            'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                            continue
                        elif 'Normalization of intronic variants is not supported' in error:
                            # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                            self.disparity_deletion_in = ['transcript', 'Requires Analysis']

                # Recreate hgvs_genomic
                if self.disparity_deletion_in[0] == 'transcript':
                    hgvs_genomic = hgvs_not_delins

                # Find oddly placed gaps where the tx variant is encompassed in the gap
                if self.disparity_deletion_in[0] == 'false' and (
                        possibility_counter == 3 or possibility_counter == 4):
                    rg = self.variant.reverse_normalizer.normalize(hgvs_not_delins)
                    rtx = self.validator.vm.g_to_t(rg, self.tx_hgvs_not_delins.ac,
                                                   alt_aln_method=self.validator.alt_aln_method)
                    fg = self.variant.hn.normalize(hgvs_not_delins)
                    ftx = self.validator.vm.g_to_t(fg, self.tx_hgvs_not_delins.ac,
                                                   alt_aln_method=self.validator.alt_aln_method)
                    if (rtx.posedit.pos.start.offset == 0 and rtx.posedit.pos.end.offset == 0) and (
                            ftx.posedit.pos.start.offset != 0 and ftx.posedit.pos.end.offset != 0):
                        exons = self.variant.map_dat.mapped_exons(
                                ftx.ac, hgvs_not_delins.ac,
                                alt_aln_method=self.validator.alt_aln_method)

                        exonic = False
                        for ex_test in exons:
                            if ftx.posedit.pos.start.base in range(ex_test[6], ex_test[7]) and \
                                    ftx.posedit.pos.end.base in range(ex_test[6], ex_test[7]):
                                exonic = True
                        if exonic is True:
                            hgvs_genomic = fg
                            self.hgvs_genomic_5pr = fg
                            try:
                                self.tx_hgvs_not_delins = self.validator.vm.c_to_n(ftx)
                            except Exception:
                                self.tx_hgvs_not_delins = ftx
                            self.disparity_deletion_in = ['transcript', 'Requires Analysis']

                # Pre-processing of tx_hgvs_not_delins
                try:
                    if self.tx_hgvs_not_delins.posedit.edit.alt is None:
                        self.tx_hgvs_not_delins.posedit.edit.alt = ''
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        self.tx_hgvs_not_delins = hgvs_dup_to_delins(self.tx_hgvs_not_delins)


                # Has a hard set variant been identified from pushes?
                hard_set_outputs = False
                if disparity_info != ['false', 'false'] and len(disparity_info) == 4:
                    self.tx_hgvs_not_delins = disparity_info[2]
                    self.disparity_deletion_in = [disparity_info[0], disparity_info[1]]
                    hgvs_refreshed_variant = hgvs_coding
                    hgvs_alt_genomic = possibility
                    hard_set_outputs = True

                elif self.disparity_deletion_in[0] == 'transcript':
                    # ANY VARIANT WHOLLY WITHIN THE GAP
                    pass
                    hgvs_refreshed_variant = self.transcript_disparity(
                        reverse_normalized_hgvs_genomic,
                        hgvs_genomic,
                        4,
                    )

                # GAP IN THE CHROMOSOME
                elif self.disparity_deletion_in[0] == 'chromosome':
                    # amend_RefSeqGene = 'true'
                    if possibility_counter == 3:
                        hgvs_refreshed_variant = stash_tx_right
                    elif possibility_counter == 4:
                        hgvs_refreshed_variant = stash_tx_left
                    else:
                        hgvs_refreshed_variant = chromosome_normalized_hgvs_coding
                else:
                    # Keep the same by re-setting rel_var
                    hgvs_refreshed_variant = hgvs_coding

                # Edit the output
                if 'NM_' in str(hgvs_refreshed_variant.ac) and 'c' not in str(hgvs_refreshed_variant.type):
                    hgvs_refreshed_variant = self.variant.no_norm_evm.n_to_c(hgvs_refreshed_variant)

                try:
                    self.variant.hn.normalize(hgvs_refreshed_variant)
                except Exception as e:
                    error = str(e)
                    # Ensure the final variant is not intronic nor does it cross exon boundaries
                    if 'Normalization of intronic variants is not supported' in error or \
                            'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        hgvs_refreshed_variant = saved_hgvs_coding
                    else:
                        continue

                # Quick check to make sure the coding variant has not changed UNLESS HARD SET
                try:
                    to_test = self.variant.hn.normalize(hgvs_refreshed_variant)
                except:
                    to_test = hgvs_refreshed_variant
                if not _same_edit(to_test, hgvs_coding):
                    # Try the next available genomic option
                    if hgvs_coding.posedit.edit.type == 'identity' and to_test.posedit.edit.type == 'identity':
                        hgvs_coding = to_test
                    elif hard_set_outputs is True:
                        hgvs_coding = to_test
                    else:
                        continue

                if hard_set_outputs is False:
                    # Update hgvs_genomic
                    hgvs_alt_genomic = self.validator.myvm_t_to_g(hgvs_refreshed_variant, alt_chr,
                                                                  self.variant.no_norm_evm, self.variant.hn,
                                                                  self.variant.map_dat)
                    if hgvs_alt_genomic.posedit.edit.type == 'identity':
                        re_c = self.validator.vm.g_to_t(hgvs_alt_genomic, hgvs_refreshed_variant.ac,
                                                        alt_aln_method=self.validator.alt_aln_method)
                        if (self.variant.hn.normalize(re_c)) != (self.variant.hn.normalize(hgvs_refreshed_variant)):
                            shuffle_left_g = copy.copy(hgvs_alt_genomic)
                            shuffle_left_g.posedit.edit.ref = ''
                            shuffle_left_g.posedit.edit.alt = ''
                            shuffle_left_g.posedit.pos.start.base = shuffle_left_g.posedit.pos.start.base - 1
                            shuffle_left_g.posedit.pos.end.base = shuffle_left_g.posedit.pos.end.base - 1
                            shuffle_left_g = self.variant.reverse_normalizer.normalize(shuffle_left_g)
                            re_c = self.validator.vm.g_to_t(shuffle_left_g, hgvs_refreshed_variant.ac,
                                                            alt_aln_method=self.validator.alt_aln_method)
                            if (self.variant.hn.normalize(re_c)) != (self.variant.hn.normalize(hgvs_refreshed_variant)):
                                hgvs_alt_genomic = shuffle_left_g

                                # If it is intronic, these vairables will not have been set

            # Break if gap has been detected
            if self.disparity_deletion_in[0] != 'false':
                break

        # Normailse hgvs_genomic
        try:
            hgvs_alt_genomic = self.variant.hn.normalize(hgvs_alt_genomic)
        except vvhgvs.exceptions.HGVSError as e:
            # Strange error caused by gap in genomic
            error = str(e)
            if 'base start position must be <= end position' in error and self.disparity_deletion_in[0] == 'chromosome':
                if hgvs_alt_genomic.posedit.edit.type == 'delins':
                    start = hgvs_alt_genomic.posedit.pos.start.base
                    end = hgvs_alt_genomic.posedit.pos.end.base
                    lhb = self.validator.sf.fetch_seq(str(hgvs_alt_genomic.ac), end - 1, end)
                    rhb = self.validator.sf.fetch_seq(str(hgvs_alt_genomic.ac), start - 1, start)
                    hgvs_alt_genomic.posedit.edit.ref = lhb + rhb
                    hgvs_alt_genomic.posedit.edit.alt = lhb + hgvs_alt_genomic.posedit.edit.alt + rhb
                    hgvs_alt_genomic.posedit.pos.start.base = end
                    hgvs_alt_genomic.posedit.pos.end.base = start
                    hgvs_alt_genomic = self.variant.hn.normalize(hgvs_alt_genomic)
                if hgvs_alt_genomic.posedit.edit.type == 'del':
                    start = hgvs_alt_genomic.posedit.pos.start.base
                    end = hgvs_alt_genomic.posedit.pos.end.base
                    lhb = self.validator.sf.fetch_seq(str(hgvs_alt_genomic.ac), end - 1, end)
                    rhb = self.validator.sf.fetch_seq(str(hgvs_alt_genomic.ac), start - 1, start)
                    hgvs_alt_genomic.posedit.edit.ref = lhb + rhb
                    hgvs_alt_genomic.posedit.edit.alt = lhb + rhb
                    hgvs_alt_genomic.posedit.pos.start.base = end
                    hgvs_alt_genomic.posedit.pos.end.base = start
                    hgvs_alt_genomic = self.variant.hn.normalize(hgvs_alt_genomic)

        # check for flanking substitutions which should be dels due to gap in transcript
        try:
            check_flank_genomic = self.validator.myvm_t_to_g(hgvs_refreshed_variant,
                                                             hgvs_alt_genomic.ac,
                                                             self.variant.no_norm_evm,
                                                             self.variant.hn,
                                                             self.variant.map_dat)

            if ((hgvs_alt_genomic.posedit.edit.type == hgvs_refreshed_variant.posedit.edit.type and
                    hgvs_alt_genomic.posedit.edit.type == "sub") and
                    check_flank_genomic.posedit.edit.type == "del"):

                hgvs_alt_genomic = check_flank_genomic
        except UnboundLocalError:
            pass

        pass
        return hgvs_alt_genomic, hgvs_coding


    def dup_ins_5prime_shift(
            self,
            stored_hgvs_not_delins,
            saved_hgvs_coding,
    ):
        """
        Adjust a single-position genomic variant after 5-prime shifting across
        an alignment gap.

        This path applies only to single-position variants represented as a
        duplication, insertion, or substitution. For duplications and insertions,
        expand the genomic interval by one base where required and refresh the
        reference sequence.

        The previous transcript-span comparison has been removed because both
        outcomes performed exactly the same adjustment.
        """
        hgvs_not_delins = copy.deepcopy(stored_hgvs_not_delins)

        position = hgvs_not_delins.posedit.pos
        edit = hgvs_not_delins.posedit.edit
        genomic_edit_type = self.hgvs_genomic_5pr.posedit.edit.type

        # This handling only applies to a single-position variant.
        if position.start.base != position.end.base:
            return hgvs_not_delins

        # Only duplications and insertion-containing edits require adjustment.
        if "dup" not in genomic_edit_type and "ins" not in genomic_edit_type:
            return hgvs_not_delins

        # Extend the genomic interval by one base.
        position.end.base = position.start.base + 1

        # A delins already has the required reference-spanning representation.
        if "ins" in genomic_edit_type and "del" in genomic_edit_type:
            return hgvs_not_delins

        # Duplications and pure insertions require the flanking reference bases
        # to reconstruct the expanded representation.
        start = position.start.base - 1
        end = position.end.base

        ref_bases = self.validator.sf.fetch_seq(
            hgvs_not_delins.ac,
            start,
            end,
        )

        edit.ref = ref_bases
        edit.alt = (
                ref_bases[:1]
                + edit.alt[1:]
                + ref_bases[1:]
        )

        return hgvs_not_delins

    def remove_offsetting_to_span_gap(
            self,
            rn_tx_hgvs_not_delins,
    ):
        """
        Remove transcript offsets so the variant spans the alignment gap.

        The end position is extended by one base to create a mappable interval,
        and the edit sequence is cleared before remapping.
        """
        position = rn_tx_hgvs_not_delins.posedit.pos
        edit = rn_tx_hgvs_not_delins.posedit.edit

        position.start.offset = 0
        position.end.offset = 0
        position.end.base += 1

        edit.ref = ""

        try:
            edit.alt = ""
        except AttributeError as error:
            # Some HGVS edit types may not expose a writable alt attribute.
            pass

        return rn_tx_hgvs_not_delins

    def move_tx_end_base_to_next_nonoffset(
            self,
            rn_tx_hgvs_not_delins,
            saved_hgvs_coding,
            back=True,
    ):
        """
        Move the transcript end position to the next non-offset base.

        The end offset is removed before rebuilding the variant across the
        transcript/genome alignment gap. When moving back across the gap, the
        reference base crossed by the move is appended to the alternate
        sequence.

        The adjusted transcript variant is then remapped transcript -> genome
        -> transcript using the non-normalising mapper. Normalisation must not
        be introduced here because it can reposition the variant relative to
        the alignment gap.
        """
        no_norm_evm = self.variant.no_norm_evm
        position = rn_tx_hgvs_not_delins.posedit.pos
        edit = rn_tx_hgvs_not_delins.posedit.edit

        position.end.offset = 0
        edit.ref = ""

        if back:
            # Preserve the transcript reference base crossed while removing
            # the offset.
            end = position.end.base
            edit.alt += self.validator.sf.fetch_seq(
                rn_tx_hgvs_not_delins.ac,
                end - 1,
                end,
            )
        else:
            # Move the end to the next available non-offset transcript base.
            position.end.base = (
                self.tx_hgvs_not_delins.posedit.pos.end.base + 1
            )

        # Coding RefSeq transcripts are represented as n. variants during
        # this gap-processing stage and must be returned to c. coordinates
        # before transcript-to-genome mapping.
        if rn_tx_hgvs_not_delins.ac.startswith("NM_"):
            test_tx_var = no_norm_evm.n_to_c(
                rn_tx_hgvs_not_delins
            )
        else:
            test_tx_var = rn_tx_hgvs_not_delins

        # Rebuild the genomic representation without normalising across the
        # alignment gap.
        hgvs_not_delins = self.validator.myevm_t_to_g(
            test_tx_var,
            no_norm_evm,
            self.variant.primary_assembly,
            self.variant.hn,
            self.variant,
        )

        # Map the rebuilt genomic representation back to transcript
        # coordinates.
        rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(
            hgvs_not_delins,
            saved_hgvs_coding.ac,
        )

        return rn_tx_hgvs_not_delins, hgvs_not_delins

    def move_tx_start_base_to_previous_nonoffset(
            self,
            rn_tx_hgvs_not_delins,
            saved_hgvs_coding,
            with_base_subtract=False,
    ):
        """
        Move the transcript start position to the previous non-offset base.

        The adjusted transcript variant is remapped transcript -> genome ->
        transcript using the non-normalising mapper. If removing the offset
        preserves the reference sequence across the round trip, sequence
        associated with the original offset is retained in the alternate.
        """
        no_norm_evm = self.variant.no_norm_evm

        # Preserve the original state because the round trip below may remove
        # information associated with the transcript offset.
        stored_tx_variant = copy.deepcopy(rn_tx_hgvs_not_delins)

        position = rn_tx_hgvs_not_delins.posedit.pos
        edit = rn_tx_hgvs_not_delins.posedit.edit

        position.start.offset = 0

        if with_base_subtract and position.start.base > 1:
            position.start.base -= 1

        edit.ref = ""

        # Coding RefSeq transcripts are represented as n. variants during this
        # gap-processing stage but require c. coordinates for T -> G mapping.
        if rn_tx_hgvs_not_delins.ac.startswith("NM_"):
            try:
                test_tx_var = no_norm_evm.n_to_c(
                    rn_tx_hgvs_not_delins
                )

            except vvhgvs.exceptions.HGVSInvalidVariantError as error:
                # Retain the historical recovery path only for the specific
                # state in which an n. variant was expected.
                if "Expected n. variant;" not in str(error):
                    raise

                rn_tx_hgvs_not_delins = self.validator.vm.c_to_n(
                    rn_tx_hgvs_not_delins
                )
                test_tx_var = no_norm_evm.n_to_c(
                    rn_tx_hgvs_not_delins
                )

        else:
            test_tx_var = rn_tx_hgvs_not_delins

        # Rebuild the genomic representation without normalising across the
        # alignment gap.
        hgvs_not_delins = self.validator.myevm_t_to_g(
            test_tx_var,
            no_norm_evm,
            self.variant.primary_assembly,
            self.variant.hn,
            self.variant,
        )

        try:
            rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(
                hgvs_not_delins,
                saved_hgvs_coding.ac,
            )

        except vvhgvs.exceptions.HGVSInvalidIntervalError:
            # Some gap-spanning intervals cannot be represented by the reverse
            # G -> N mapping. In that case retain the adjusted transcript
            # representation produced above.
            rn_tx_hgvs_not_delins = test_tx_var

        stored_start = stored_tx_variant.posedit.pos.start
        stored_edit = stored_tx_variant.posedit.edit

        remapped_start = rn_tx_hgvs_not_delins.posedit.pos.start
        remapped_edit = rn_tx_hgvs_not_delins.posedit.edit

        if (
                stored_start.offset != 0
                and remapped_start.offset == 0
                and stored_edit.ref == remapped_edit.ref
        ):
            # Preserve sequence associated with the offset when the round trip
            # removed the offset without changing the reference.
            remapped_edit.alt += stored_edit.ref[:stored_start.offset]

        else:
            remapped_start.offset = 0

        return rn_tx_hgvs_not_delins, hgvs_not_delins

    def c2_pos_edit(self, hgvs_genomic):
        """
        Refresh a transcript posedit where the gap affects the start position.

        Extend the transcript variant by one preceding non-offset base, map
        both transcript components to the genome, combine their reference and
        alternate sequences, and map the reconstructed genomic variant back
        to the transcript.
        """
        vm = self.validator.vm
        aln_method = self.validator.alt_aln_method
        genomic_ac = hgvs_genomic.ac

        try:
            c2 = vm.n_to_c(self.tx_hgvs_not_delins)
        except vvhgvs.exceptions.HGVSError:
            # NR_ transcripts cannot be converted to coding coordinates.
            c2 = self.tx_hgvs_not_delins

        c1 = copy.deepcopy(c2)

        c1_position = c1.posedit.pos
        c1_edit = c1.posedit.edit

        c1_position.start.base = c2.posedit.pos.start.base - 1
        c1_position.start.offset = 0
        c1_position.end = c2.posedit.pos.start
        c1_edit.ref = ""
        c1_edit.alt = ""

        if self.orientation != -1:
            g1 = vm.t_to_g(
                c1,
                genomic_ac,
                alt_aln_method=aln_method,
            )
            g2 = vm.t_to_g(
                c2,
                genomic_ac,
                alt_aln_method=aln_method,
            )

            # The preceding transcript component contributes reference
            # sequence only to the reconstructed genomic interval.
            g1.posedit.edit.alt = g1.posedit.edit.ref

        else:
            # Reverse-strand mappings require the genomic components in the
            # opposite order.
            g1 = vm.t_to_g(
                c2,
                genomic_ac,
                alt_aln_method=aln_method,
            )
            g2 = vm.t_to_g(
                c1,
                genomic_ac,
                alt_aln_method=aln_method,
            )

            g2.posedit.edit.alt = g2.posedit.edit.ref

        g1_edit = g1.posedit.edit
        g2_edit = g2.posedit.edit

        g3 = copy.deepcopy(g1)
        g3.posedit.pos.end.base = g2.posedit.pos.end.base
        g3.posedit.edit.ref = (
            g1_edit.ref
            + g2_edit.ref[1:]
        )
        g3.posedit.edit.alt = (
            g1_edit.alt
            + g2_edit.alt[1:]
        )

        return vm.g_to_t(
            g3,
            c1.ac,
            alt_aln_method=aln_method,
        )

    def c1_pos_edit(self, hgvs_genomic):
        """
        Refresh a transcript posedit where the gap affects the end position.

        Extend the transcript variant to the next non-offset base, map both
        transcript components to the genome, and reconstruct the complete
        reference sequence across the alignment gap.

        If the combined genomic variant cannot be mapped directly back to the
        transcript, reconstruct the transcript and intronic reference
        components separately.
        """
        vm = self.validator.vm
        aln_method = self.validator.alt_aln_method
        genomic_ac = hgvs_genomic.ac

        try:
            c1 = vm.n_to_c(self.tx_hgvs_not_delins)
        except vvhgvs.exceptions.HGVSError:
            # NR_ transcripts cannot be converted to coding coordinates.
            c1 = self.tx_hgvs_not_delins

        c2 = copy.deepcopy(c1)

        c1_position = c1.posedit.pos
        c2_position = c2.posedit.pos
        c2_edit = c2.posedit.edit

        c2_position.start = c1_position.end
        c2_position.end.base = c1_position.end.base + 1
        c2_position.end.offset = 0
        c2_edit.ref = ""
        c2_edit.alt = ""

        if self.orientation != -1:
            g1 = vm.t_to_g(
                c1,
                genomic_ac,
                alt_aln_method=aln_method,
            )
            g2 = vm.t_to_g(
                c2,
                genomic_ac,
                alt_aln_method=aln_method,
            )

            # The following transcript component contributes reference
            # sequence only to the reconstructed genomic interval.
            g2.posedit.edit.alt = g2.posedit.edit.ref

        else:
            # Reverse-strand mappings require the genomic components in the
            # opposite order.
            g1 = vm.t_to_g(
                c2,
                genomic_ac,
                alt_aln_method=aln_method,
            )
            g2 = vm.t_to_g(
                c1,
                genomic_ac,
                alt_aln_method=aln_method,
            )

            g1.posedit.edit.alt = g1.posedit.edit.ref

        g1_edit = g1.posedit.edit
        g2_edit = g2.posedit.edit

        g3 = copy.deepcopy(g1)
        g3.posedit.pos.end.base = g2.posedit.pos.end.base
        g3.posedit.edit.ref = (
            g1_edit.ref
            + g2_edit.ref[1:]
        )
        g3.posedit.edit.alt = (
            g1_edit.alt
            + g2_edit.alt[1:]
        )

        try:
            return vm.g_to_t(
                g3,
                c1.ac,
                alt_aln_method=aln_method,
            )

        except vvhgvs.exceptions.HGVSError:
            # A gap-spanning interval may not map directly back to the
            # transcript. Reconstruct its exonic and intronic reference
            # components separately instead.
            c_tx_part = copy.deepcopy(c1)
            tx_position = c_tx_part.posedit.pos
            tx_edit = c_tx_part.posedit.edit

            tx_position.end.offset = 0
            tx_edit.ref = ""
            tx_edit.alt = ""

            vm._replace_reference(c_tx_part)

            c_intronic_part = copy.deepcopy(c1)
            intronic_position = c_intronic_part.posedit.pos
            intronic_edit = c_intronic_part.posedit.edit

            intronic_position.start.base = intronic_position.end.base
            intronic_position.start.offset = 0
            intronic_edit.ref = ""
            intronic_edit.alt = ""

            g_intronic_ref_eq = vm.t_to_g(
                c_intronic_part,
                genomic_ac,
                alt_aln_method=aln_method,
            )

            intronic_ref = g_intronic_ref_eq.posedit.edit.ref

            if self.orientation == -1:
                intronic_ref = simple_dna_revcomp(intronic_ref)

            c3 = copy.deepcopy(c1)
            c3.posedit.edit.ref = (
                c_tx_part.posedit.edit.ref
                + intronic_ref[1:]
            )

            return c3

    def transcript_disparity(
            self,
            reverse_normalized_hgvs_genomic,
            hgvs_genomic,
            running_option,
    ):
        """
        Correct transcript/genome disparity caused by an alignment gap.

        This function handles variants whose transcript representation either
        contains offsets into an alignment gap or spans a known gap without
        retaining offsets.

        Mapping, normalisation, sequence fetching and gap reconstruction are kept
        in their established order because these operations can alter the HGVS
        representation and may fetch reference sequence.

        Further consolidation of the reconstruction paths should only be performed
        after the individual branches have been exercised with known variants.
        """
        vm = self.validator.vm
        aln_method = self.validator.alt_aln_method
        tx_variant = self.tx_hgvs_not_delins

        tx_start = tx_variant.posedit.pos.start
        tx_end = tx_variant.posedit.pos.end

        start_offset = tx_start.offset
        end_offset = tx_end.offset

        # Both ends of the transcript variant are offset into an alignment gap.
        if start_offset != 0 and end_offset != 0:
            self.gapped_transcripts += f" {tx_variant.ac}"

            # Work on an independent copy because the gap-filling representation
            # is modified extensively before mapping.
            tx_gap_fill_variant = copy.deepcopy(tx_variant)

            # Dup edits do not expose an alt attribute. Convert them to delins
            # when necessary so the gap-filling representation has an ALT.
            try:
                if tx_gap_fill_variant.posedit.edit.alt is None:
                    tx_gap_fill_variant.posedit.edit.alt = ""
            except Exception as error:
                if str(error) == "'Dup' object has no attribute 'alt'":
                    tx_gap_fill_variant = hgvs_dup_to_delins(
                        tx_gap_fill_variant
                    )

            gap_start = tx_gap_fill_variant.posedit.pos.start
            gap_end = tx_gap_fill_variant.posedit.pos.end
            gap_edit = tx_gap_fill_variant.posedit.edit

            # Move the transcript interval out of the offset region so that an
            # equivalent genomic interval spanning the alignment gap can be made.
            if gap_start.offset < 0:
                gap_start.base -= 1
                gap_start.offset = 0
                gap_end.offset = 0
                gap_edit.alt = ""
                gap_edit.ref = ""

            elif gap_start.offset > 0:
                gap_start.offset = 0
                gap_end.base += 1
                gap_end.offset = 0
                gap_edit.alt = ""
                gap_edit.ref = ""

            try:
                tx_gap_fill_variant = vm.n_to_c(
                    tx_gap_fill_variant
                )
            except Exception as error:
                # Retain the established fallback for transcript representations
                # that cannot be converted to coding coordinates.
                pass

            genomic_gap_fill_variant = vm.t_to_g(
                tx_gap_fill_variant,
                reverse_normalized_hgvs_genomic.ac,
                alt_aln_method=aln_method,
            )

            genomic_gap_fill_variant.posedit.edit.alt = (
                genomic_gap_fill_variant.posedit.edit.ref
            )

            try:
                c_tx_hgvs_not_delins = vm.n_to_c(
                    tx_variant
                )
            except Exception:
                c_tx_hgvs_not_delins = copy.copy(
                    tx_variant
                )

            genomic_gap_fill_variant_alt = self.validator.myvm_t_to_g(
                c_tx_hgvs_not_delins,
                self.hgvs_genomic_5pr.ac,
                self.variant.no_norm_evm,
                self.variant.hn,
                self.variant.map_dat,
            )

            # Ensure that the alternate genomic representation exposes an ALT.
            # Dup handling is retained here for later targeted investigation.
            try:
                if genomic_gap_fill_variant_alt.posedit.edit.alt is None:
                    genomic_gap_fill_variant_alt.posedit.edit.alt = "X"
            except Exception as error:
                if str(error) == "'Dup' object has no attribute 'alt'":
                    genomic_gap_fill_variant = hgvs_dup_to_delins(
                        genomic_gap_fill_variant
                    )
                    genomic_gap_fill_variant_alt = hgvs_dup_to_delins(
                        genomic_gap_fill_variant_alt
                    )

            alt_edit = genomic_gap_fill_variant_alt.posedit.edit
            alt_position = genomic_gap_fill_variant_alt.posedit.pos

            # Insertions need their flanking genomic reference bases added before
            # the replacement sequence can be projected across the gap interval.
            if alt_edit.type == "ins":
                append_ref = self.validator.sf.fetch_seq(
                    genomic_gap_fill_variant_alt.ac,
                    alt_position.start.base - 1,
                    alt_position.end.base,
                )
                alt_edit.alt = (
                        append_ref[0]
                        + alt_edit.alt
                        + append_ref[1]
                )

            gap_edit = genomic_gap_fill_variant.posedit.edit
            gap_position = genomic_gap_fill_variant.posedit.pos

            reference_bases = list(gap_edit.ref)

            if alt_edit.alt is not None:
                alternate_bases = list(alt_edit.alt)
            else:
                # A deletion has no inserted sequence. X is used internally to
                # mark deleted positions and is removed after reconstruction.
                alternate_bases = [
                    "X"
                    for _ in alt_edit.ref
                ]

            ref_start = gap_position.start.base
            alt_start = alt_position.start.base

            ref_base_dict = {
                ref_start + index: base
                for index, base in enumerate(reference_bases)
            }

            # Force the alternate representation into an interval-based
            # delete/insert form. Deleted positions are represented temporarily
            # by X so they can be removed after the complete sequence is built.
            alt_base_dict = {}

            for position in range(
                    alt_position.start.base,
                    alt_position.end.base + 1,
            ):
                if position == alt_start:
                    alt_base_dict[position] = "".join(
                        alternate_bases
                    )
                else:
                    alt_base_dict[position] = "X"

            alternate_sequence_bases = []

            for position in range(
                    gap_position.start.base,
                    gap_position.end.base + 1,
            ):
                if position in alt_base_dict:
                    alternate_sequence_bases.append(
                        alt_base_dict[position]
                    )
                else:
                    alternate_sequence_bases.append(
                        ref_base_dict[position]
                    )

            alternate_sequence = "".join(
                alternate_sequence_bases
            ).replace("X", "")

            gap_edit.alt = alternate_sequence

            hgvs_refreshed_variant = vm.g_to_t(
                genomic_gap_fill_variant,
                tx_gap_fill_variant.ac,
                alt_aln_method=aln_method,
            )

            # Record the number of reference bases represented by the alignment
            # disparity, excluding the two flanking bases.
            gap_size = str(len(gap_edit.ref) - 2)
            self.disparity_deletion_in[1] = [gap_size]

        else:
            # A transcript gap may be known even when neither endpoint retains an
            # offset. This can occur when the variant itself spans the gap.
            if start_offset == 0 and end_offset == 0:
                try:
                    c1 = vm.n_to_c(tx_variant)
                except Exception:
                    c1 = tx_variant

                g3 = self.validator.nr_vm.t_to_g(
                    c1,
                    hgvs_genomic.ac,
                    alt_aln_method=aln_method,
                )

                # Normalisation is intentionally retained here. The following
                # tests compare the remapped interval with the established
                # normalised genomic representation to determine whether the gap
                # has already been handled.
                hgvs_genomic_norm = self.variant.hn.normalize(
                    hgvs_genomic
                )

                g3_span = (
                        g3.posedit.pos.end.base
                        - g3.posedit.pos.start.base
                )
                genomic_span = (
                        hgvs_genomic_norm.posedit.pos.end.base
                        - hgvs_genomic_norm.posedit.pos.start.base
                )

                if (
                        g3_span > genomic_span
                        and hgvs_genomic_norm.posedit.edit.type == "del"
                        and (
                        g3.posedit.pos.end.base
                        == hgvs_genomic_norm.posedit.pos.end.base
                )
                ):
                    return tx_variant

                if (
                        g3_span > genomic_span
                        and hgvs_genomic_norm.posedit.edit.type == "del"
                        and (
                        g3.posedit.pos.start.base
                        < hgvs_genomic_norm.posedit.pos.start.base
                )
                        and (
                        g3.posedit.pos.end.base
                        + int(self.disparity_deletion_in[1])
                        < hgvs_genomic_norm.posedit.pos.end.base
                )
                ):
                    return tx_variant

                g3.posedit.pos.end.base = (
                        g3.posedit.pos.start.base
                        + len(g3.posedit.edit.ref)
                        - 1
                )

                try:
                    c2 = vm.g_to_t(
                        g3,
                        c1.ac,
                        alt_aln_method=aln_method,
                    )

                    if (
                            c2.posedit.pos.start.offset != 0
                            or c2.posedit.pos.end.offset != 0
                    ):
                        self.tx_hgvs_not_delins = c2

                        try:
                            self.tx_hgvs_not_delins = vm.c_to_n(
                                self.tx_hgvs_not_delins
                            )
                        except vvhgvs.exceptions.HGVSError as error:
                            pass

                except vvhgvs.exceptions.HGVSInvalidVariantError as error:
                    pass

                # The previous block may replace the stored transcript variant.
                tx_variant = self.tx_hgvs_not_delins
                tx_start = tx_variant.posedit.pos.start
                tx_end = tx_variant.posedit.pos.end

                start_offset = tx_start.offset
                end_offset = tx_end.offset

            if start_offset > 0 and end_offset <= 0:
                hgvs_refreshed_variant = self.c2_pos_edit(
                    hgvs_genomic
                )

            elif end_offset > 0 and start_offset <= 0:
                self.gapped_transcripts += f" {tx_variant.ac}"

                try:
                    hgvs_refreshed_variant = self.c1_pos_edit(
                        hgvs_genomic
                    )
                except vvhgvs.exceptions.HGVSDataNotAvailableError:
                    hgvs_refreshed_variant = tx_variant

            elif start_offset < 0 and end_offset >= 0:
                hgvs_refreshed_variant = self.c2_pos_edit(
                    hgvs_genomic
                )

            elif end_offset < 0 and start_offset >= 0:
                self.gapped_transcripts += f" {tx_variant.ac}"

                # This path handles variation associated with the first copy.
                if running_option == 1:
                    try:
                        c1 = vm.n_to_c(tx_variant)
                    except Exception:
                        c1 = tx_variant

                    c2 = copy.deepcopy(c1)
                    c2.posedit.pos.start = c1.posedit.pos.end
                    c2.posedit.pos.end.base = c1.posedit.pos.end.base
                    c2.posedit.pos.end.offset = 0
                    c2.posedit.edit.ref = ""
                    c2.posedit.edit.alt = ""

                    g2 = vm.t_to_g(
                        c2,
                        self.variant.hgvs_genomic.ac,
                        alt_aln_method=aln_method,
                    )

                    c2 = vm.g_to_t(
                        g2,
                        c2.ac,
                        alt_aln_method=aln_method,
                    )

                    alternate = (
                            c1.posedit.edit.alt
                            + c2.posedit.edit.ref[1:]
                    )

                    c3 = copy.deepcopy(c1)
                    c3.posedit.pos.end = c2.posedit.pos.end
                    c3.posedit.edit.ref = ""
                    c3.posedit.edit.alt = alternate

                    hgvs_refreshed_variant = c3

                else:
                    hgvs_refreshed_variant = self.c1_pos_edit(
                        hgvs_genomic
                    )

            else:
                # This path handles variation associated with the second copy.
                if running_option == 2:
                    tx_variant.posedit.pos.end.base = (
                            tx_variant.posedit.pos.start.base
                            + len(tx_variant.posedit.edit.ref)
                            - 1
                    )

                elif running_option != 4:
                    self.gapped_transcripts += f" {tx_variant.ac}"

                hgvs_refreshed_variant = tx_variant

        return hgvs_refreshed_variant

    def edit_output(self, hgvs_refreshed_variant, saved_hgvs_coding):
        """
        Normalise and tidy the refreshed transcript variant.

        If normalisation cannot be performed because the variant is intronic or
        spans an exon/intron or UTR/exon boundary, retain the original coding
        representation.
        """
        if (
                hgvs_refreshed_variant.ac.startswith("NM_")
                and hgvs_refreshed_variant.type != "c"
        ):
            hgvs_refreshed_variant = self.variant.evm.n_to_c(
                hgvs_refreshed_variant
            )

        try:
            hgvs_refreshed_variant = self.variant.hn.normalize(
                hgvs_refreshed_variant
            )

            pass

            edit = hgvs_refreshed_variant.posedit.edit
            pos = hgvs_refreshed_variant.posedit.pos

            if edit.type == "delins":
                if edit.ref[-1] == edit.alt[-1]:
                    edit.ref = edit.ref[:-1]
                    edit.alt = edit.alt[:-1]
                    pos.end.base -= 1

                    hgvs_refreshed_variant = self.variant.hn.normalize(
                        hgvs_refreshed_variant
                    )

                elif edit.ref[0] == edit.alt[0]:
                    edit.ref = edit.ref[1:]
                    edit.alt = edit.alt[1:]
                    pos.start.base += 1

                    hgvs_refreshed_variant = self.variant.hn.normalize(
                        hgvs_refreshed_variant
                    )

        except Exception as error:
            error_message = str(error)

            # Normalisation cannot safely process these transcript-coordinate
            # states. Retain the original coding representation instead.
            unsupported_normalisation = (
                    "Normalization of intronic variants is not supported"
                    in error_message
                    or
                    "Unsupported normalization of variants spanning the "
                    "exon-intron boundary"
                    in error_message
                    or
                    "Unsupported normalization of variants spanning the "
                    "UTR-exon boundary"
                    in error_message
            )

            if unsupported_normalisation:
                hgvs_refreshed_variant = saved_hgvs_coding

        return hgvs_refreshed_variant

    def logic_check(
            self,
            hgvs_not_delins,
            rn_tx_hgvs_not_delins,
            hgvs_coding,
            do_continue=False,
            offset_check=False,
    ):
        """
        Compare genomic and transcript reference lengths to identify alignment
        disparity.

        Where the initial comparison is inconclusive, inspect stored genomic
        possibilities and attempt to recover the corresponding transcript
        representation.
        """
        genomic_edit = hgvs_not_delins.posedit.edit
        tx_edit = rn_tx_hgvs_not_delins.posedit.edit

        if genomic_edit.ref is None:
            genomic_edit.ref = ""

        if tx_edit.ref is None:
            tx_edit.ref = ""

        genomic_ref_length = len(genomic_edit.ref)
        tx_ref_length = len(tx_edit.ref)

        if genomic_ref_length < tx_ref_length:
            gap_length = tx_ref_length - genomic_ref_length
            self.disparity_deletion_in = [
                "chromosome",
                gap_length,
            ]

        elif genomic_ref_length > tx_ref_length:
            gap_length = genomic_ref_length - tx_ref_length
            self.disparity_deletion_in = [
                "transcript",
                gap_length,
            ]

        else:
            re_capture_tx_variant = []

            for an_internal_possibility in self.hgvs_genomic_possibilities:
                try:
                    internal_possibility = an_internal_possibility[1][3]
                except IndexError:
                    internal_possibility = an_internal_possibility[0]

                if internal_possibility == "":
                    continue

                hgvs_t_possibility = self.validator.vm.g_to_t(
                    internal_possibility,
                    hgvs_coding.ac,
                    alt_aln_method=self.validator.alt_aln_method,
                )

                if hgvs_t_possibility.posedit.edit.type == "ins":
                    try:
                        hgvs_t_possibility = self.validator.vm.c_to_n(
                            hgvs_t_possibility
                        )
                    except Exception as error:
                        if do_continue:
                            continue
                        pass

                    if offset_check:
                        pos = hgvs_t_possibility.posedit.pos

                        if (
                                pos.start.offset != 0
                                or pos.end.offset != 0
                        ):
                            continue

                    pos = hgvs_t_possibility.posedit.pos

                    ins_ref = self.validator.sf.fetch_seq(
                        hgvs_t_possibility.ac,
                        pos.start.base - 1,
                        pos.start.base + 1,
                    )

                    try:
                        hgvs_t_possibility = self.validator.vm.n_to_c(
                            hgvs_t_possibility
                        )
                    except Exception as error:
                        if do_continue:
                            continue
                        pass

                    edit = hgvs_t_possibility.posedit.edit
                    edit.ref = ins_ref
                    edit.alt = ins_ref[0] + edit.alt + ins_ref[1]

                if internal_possibility.posedit.edit.type == "ins":
                    pos = internal_possibility.posedit.pos
                    edit = internal_possibility.posedit.edit

                    ins_ref = self.validator.sf.fetch_seq(
                        internal_possibility.ac,
                        pos.start.base - 1,
                        pos.end.base,
                    )

                    edit.ref = ins_ref
                    edit.alt = ins_ref[0] + edit.alt + ins_ref[1]

                if (
                        len(hgvs_t_possibility.posedit.edit.ref)
                        < len(internal_possibility.posedit.edit.ref)
                ):
                    gap_length = (
                            len(internal_possibility.posedit.edit.ref)
                            - len(hgvs_t_possibility.posedit.edit.ref)
                    )

                    re_capture_tx_variant = [
                        "transcript",
                        gap_length,
                        hgvs_t_possibility,
                    ]

                    hgvs_not_delins = internal_possibility
                    self.hgvs_genomic_5pr = internal_possibility
                    break

            if re_capture_tx_variant:
                try:
                    self.tx_hgvs_not_delins = self.validator.vm.c_to_n(
                        re_capture_tx_variant[2]
                    )
                except Exception:
                    self.tx_hgvs_not_delins = re_capture_tx_variant[2]

                self.disparity_deletion_in = re_capture_tx_variant[:-1]

        return hgvs_not_delins

    def get_hgvs_seek_var(
            self,
            hgvs_genomic,
            hgvs_coding,
            ori=None,
            with_query_genomic=False,
    ):
        """
        Position a genomic variant according to transcript orientation and map it
        back to the transcript to determine whether its representation moves.
        """
        if ori is None:
            ori = self.orientation

        if ori == -1:
            try:
                query_genomic = self.variant.reverse_normalizer.normalize(
                    hgvs_genomic
                )
            except Exception:
                query_genomic = hgvs_genomic
        else:
            # Position the genomic variant at its most 3-prime representation.
            try:
                query_genomic = self.variant.hn.normalize(
                    hgvs_genomic
                )
            except Exception:
                query_genomic = hgvs_genomic

        # Refresh an intronic transcript representation against the positioned
        # genomic variant when the caller requires query-genomic information.
        if (
                with_query_genomic
                and hgvs_coding.posedit.pos.start.offset != 0
        ):
            try:
                hgvs_coding = self.variant.evm.g_to_t(
                    query_genomic,
                    hgvs_coding.ac,
                )
            except vvhgvs.exceptions.HGVSInvalidIntervalError:
                pass

        try:
            hgvs_seek_var = self.variant.evm.g_to_t(
                query_genomic,
                hgvs_coding.ac,
            )
        except vvhgvs.exceptions.HGVSError:
            hgvs_seek_var = hgvs_coding

        if with_query_genomic:
            return (
                hgvs_seek_var,
                query_genomic,
                hgvs_coding,
            )

        return hgvs_seek_var

    def rev_norm_ins(self, hgvs_coding, hgvs_genomic):
        """
        Compare the genomic mappings of the most 3-prime and most 5-prime
        transcript representations of an insertion.

        Insertions are expanded to delins-like representations by adding their
        flanking reference bases before genomic mapping. The resulting genomic
        representations are normalised and retained when they indicate sequence
        disparity relative to the transcript.
        """
        try:
            if hgvs_coding.posedit.edit.type != "ins":
                return

            vm = self.validator.vm
            aln_method = self.validator.alt_aln_method

            most_5pr_tx = copy.deepcopy(hgvs_coding)
            most_3pr_tx = self.variant.reverse_normalizer.normalize(
                hgvs_coding
            )

            try:
                n_3pr = vm.c_to_n(most_3pr_tx)
                n_5pr = vm.c_to_n(most_5pr_tx)
            except Exception:
                n_3pr = most_3pr_tx
                n_5pr = most_5pr_tx

            # Expand each insertion using its flanking transcript reference bases.
            pr3_pos = n_3pr.posedit.pos
            pr5_pos = n_5pr.posedit.pos

            pr3_ref = self.validator.sf.fetch_seq(
                hgvs_coding.ac,
                pr3_pos.start.base - 1,
                pr3_pos.end.base,
            )
            pr5_ref = self.validator.sf.fetch_seq(
                hgvs_coding.ac,
                pr5_pos.start.base - 1,
                pr5_pos.end.base,
            )

            most_3pr_tx.posedit.edit.ref = pr3_ref
            most_5pr_tx.posedit.edit.ref = pr5_ref

            most_3pr_tx.posedit.edit.alt = (
                    pr3_ref[0]
                    + most_3pr_tx.posedit.edit.alt
                    + pr3_ref[1]
            )
            most_5pr_tx.posedit.edit.alt = (
                    pr5_ref[0]
                    + most_5pr_tx.posedit.edit.alt
                    + pr5_ref[1]
            )

            genomic_3pr = vm.t_to_g(
                most_3pr_tx,
                hgvs_genomic.ac,
                alt_aln_method=aln_method,
            )
            genomic_5pr = vm.t_to_g(
                most_5pr_tx,
                hgvs_genomic.ac,
                alt_aln_method=aln_method,
            )

            # If the variant spans a gap, normalisation should produce a stable
            # genomic representation. Some reverse-strand mappings can initially
            # return the interval in reverse order; retain the established repair.
            try:
                genomic_3pr = self.variant.hn.normalize(
                    genomic_3pr
                )
            except vvhgvs.exceptions.HGVSInvalidVariantError as error:
                if str(error) == "base start position must be <= end position":
                    start = genomic_3pr.posedit.pos.start.base
                    end = genomic_3pr.posedit.pos.end.base

                    genomic_3pr.posedit.pos.start.base = end
                    genomic_3pr.posedit.pos.end.base = start

                    genomic_3pr = self.variant.hn.normalize(
                        genomic_3pr
                    )

            try:
                genomic_5pr = self.variant.hn.normalize(
                    genomic_5pr
                )
            except vvhgvs.exceptions.HGVSInvalidVariantError as error:
                if str(error) == "base start position must be <= end position":
                    start = genomic_5pr.posedit.pos.start.base
                    end = genomic_5pr.posedit.pos.end.base

                    genomic_5pr.posedit.pos.start.base = end
                    genomic_5pr.posedit.pos.end.base = start

                    genomic_5pr = self.variant.hn.normalize(
                        genomic_5pr
                    )

            # Retain the historical Dup handling for now. This should eventually
            # be replaced by explicit edit-type handling rather than exception
            # message inspection.
            try:
                if genomic_3pr.posedit.edit.alt is None:
                    genomic_3pr.posedit.edit.alt = ""
            except Exception as error:
                if str(error) == "'Dup' object has no attribute 'alt'":
                    genomic_3pr = hgvs_dup_to_delins(genomic_3pr)

            try:
                if most_3pr_tx.posedit.edit.alt is None:
                    most_3pr_tx.posedit.edit.alt = ""
            except Exception as error:
                if str(error) == "'Dup' object has no attribute 'alt'":
                    most_3pr_tx = hgvs_dup_to_delins(most_3pr_tx)

            try:
                if genomic_5pr.posedit.edit.alt is None:
                    genomic_5pr.posedit.edit.alt = ""
            except Exception as error:
                if str(error) == "'Dup' object has no attribute 'alt'":
                    genomic_5pr = hgvs_dup_to_delins(genomic_5pr)

            try:
                if most_5pr_tx.posedit.edit.alt is None:
                    most_5pr_tx.posedit.edit.alt = ""
            except Exception as error:
                if str(error) == "'Dup' object has no attribute 'alt'":
                    most_5pr_tx = hgvs_dup_to_delins(most_5pr_tx)

            if (
                    len(genomic_3pr.posedit.edit.alt)
                    < len(most_3pr_tx.posedit.edit.alt)
            ):
                self.hgvs_genomic_possibilities.append(
                    [
                        genomic_3pr,
                        ["false", "false"],
                    ]
                )

            if (
                    len(genomic_5pr.posedit.edit.alt)
                    < len(most_5pr_tx.posedit.edit.alt)
            ):
                self.hgvs_genomic_possibilities.append(
                    [
                        genomic_5pr,
                        ["false", "false"],
                    ]
                )

        except vvhgvs.exceptions.HGVSUnsupportedOperationError:
            pass


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
