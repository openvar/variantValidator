# -*- coding: utf-8 -*-
"""
Liftover between genome builds is most accurate when mapping via a RefSeq transcript.
For intergenic regions, the process is more complex.

Step 1: attempt liftover through a common RefSeq transcript.
Step 2: fall back to PyLiftover.

PyLiftover:
Lift position > check bases > lift back and confirm the original position.
"""

import copy
import logging

import vvhgvs
import vvhgvs.exceptions
import vvhgvs.sequencevariant
from pyliftover import LiftOver
from vvhgvs.assemblymapper import AssemblyMapper

from . import hgvs_utils
from . import seq_data
from . import utils
from VariantValidator.modules.hgvs_utils import hgvs_delins_parts_to_hgvs_obj
from VariantValidator.modules.transcript_map_data import TranscriptMapData


vvhgvs.global_config.formatting.max_ref_length = 1000000

logger = logging.getLogger(__name__)

LO_CACHE = {}


def liftover(
        hgvs_genomic,
        build_from,
        build_to,
        hn,
        reverse_normalizer,
        evm,
        validator,
        specify_tx=False,
        liftover_level=None,
        g_to_g=False,
        gap_map=False,
        vfo=False,
        specified_tx_variant=False,
        genomic_data_w_vcf=False,
        force_pyliftover=False,
        map_dat=False,
):
    """
    Lift a genomic HGVS variant between genome builds.

    Step 1:
        Attempt transcript-mediated liftover using a common RefSeq transcript.

    Step 2:
        If transcript-mediated liftover cannot provide a result, fall back to
        PyLiftover and confirm the result by lifting the position back to the
        original assembly.

    :param hgvs_genomic:
    :param build_from:
    :param build_to:
    :param hn:
    :param reverse_normalizer:
    :param evm:
    :param validator: Validator object
    :param specify_tx: False or a specific transcript accession
    :param liftover_level: None, True or 'primary'
    :param g_to_g: True or False
    :param gap_map: False or VariantFormatter gap_map function
    :param vfo: False or VariantFormatter VFO object
    :param specified_tx_variant: False/None or specific HGVS transcript object
    :param genomic_data_w_vcf: False or existing validated genomic data
    :param force_pyliftover: Skip transcript-mediated liftover when True
    :param map_dat: False or TranscriptMapData object
    :return: dictionary containing genomic liftover data
    """

    if isinstance(hgvs_genomic, str):
        hgvs_genomic = validator.hp.parse(hgvs_genomic)

    # ------------------------------------------------------------------
    # Resolve source build
    # ------------------------------------------------------------------

    if build_from.startswith('GRC'):
        from_set = 'grc_chr'
        alt_from_set = 'ucsc_chr'

        if '37' in build_from:
            lo_from = 'hg19'
            alt_build_from = 'hg19'
        elif '38' in build_from:
            lo_from = 'hg38'
            alt_build_from = 'hg38'
        else:
            lo_from = ''
            alt_build_from = ''

    else:
        from_set = 'ucsc_chr'
        alt_from_set = 'grc_chr'

        if '19' in build_from:
            lo_from = 'hg19'
            alt_build_from = 'GRCh37'
        elif '38' in build_from:
            lo_from = 'hg38'
            alt_build_from = 'GRCh38'
        else:
            lo_from = ''
            alt_build_from = ''

    # ------------------------------------------------------------------
    # Resolve target build
    # ------------------------------------------------------------------

    if build_to.startswith('GRC'):
        to_set = 'grc_chr'
        alt_to_set = 'ucsc_chr'

        if '37' in build_to:
            lo_to = 'hg19'
            alt_build_to = 'hg19'
        elif '38' in build_to:
            lo_to = 'hg38'
            alt_build_to = 'hg38'
        else:
            lo_to = ''
            alt_build_to = ''

    else:
        to_set = 'ucsc_chr'
        alt_to_set = 'grc_chr'

        if '19' in build_to:
            lo_to = 'hg19'
            alt_build_to = 'GRCh37'
        elif '38' in build_to:
            lo_to = 'hg38'
            alt_build_to = 'GRCh38'
        else:
            lo_to = ''
            alt_build_to = ''

    # ------------------------------------------------------------------
    # Initialise response and source VCF
    # ------------------------------------------------------------------

    if genomic_data_w_vcf:
        lifted_response = {}

        for genome_build, genomic_data in genomic_data_w_vcf.items():
            lifted_response[genome_build] = {
                hgvs_genomic.ac: genomic_data
            }

        # Un-fix mitochondrial mapping and add the extra GRCh37 mapping
        # required for compatibility with raw genomic HGVS mapping.
        if (
                'hg19' in lifted_response
                and 'NC_001807.4' in lifted_response['hg19']
                and 'grch37' not in lifted_response
        ):
            lifted_response['grch37'] = copy.copy(
                lifted_response['hg19']
            )

        from_vcf = copy.copy(
            genomic_data_w_vcf[build_from.lower()]['vcf']
        )

        # Match the structure returned by report_hgvs2vcf().
        from_vcf[from_set] = from_vcf['chr']

    else:
        lifted_response = {}

        from_vcf = hgvs_utils.report_hgvs2vcf(
            hgvs_genomic,
            lo_from,
            reverse_normalizer,
            validator.sf,
        )

        lifted_response[build_from.lower()] = {
            hgvs_genomic.ac: {
                'hgvs_genomic_description': hgvs_genomic,
                'vcf': {
                    'chr': from_vcf[from_set],
                    'pos': str(from_vcf['pos']),
                    'ref': from_vcf['ref'],
                    'alt': from_vcf['alt'],
                },
            }
        }

        lifted_response[alt_build_from.lower()] = {
            hgvs_genomic.ac: {
                'hgvs_genomic_description': hgvs_genomic,
                'vcf': {
                    'chr': from_vcf[alt_from_set],
                    'pos': str(from_vcf['pos']),
                    'ref': from_vcf['ref'],
                    'alt': from_vcf['alt'],
                },
            }
        }

    # Ensure target-build dictionaries exist without overwriting anything
    # already present in genomic_data_w_vcf.
    lifted_response.setdefault(build_to.lower(), {})
    lifted_response.setdefault(alt_build_to.lower(), {})

    # ==================================================================
    # STEP 1: TRANSCRIPT-MEDIATED LIFTOVER
    # ==================================================================

    tx_list = []

    if not g_to_g:
        rts_dict = {
            tx_dat[0]: True
            for tx_dat in validator.hdp.get_tx_for_region(
                hgvs_genomic.ac,
                'splign',
                hgvs_genomic.posedit.pos.start.base - 1,
                hgvs_genomic.posedit.pos.end.base,
            )
        }

        if evm is not None:
            for tx in evm.relevant_transcripts(hgvs_genomic):
                rts_dict[tx] = True

        tx_list = list(rts_dict)

    # An explicitly specified transcript takes precedence.
    if specify_tx is not False:
        tx_list = [specify_tx]

    if (
            tx_list
            and not force_pyliftover
            and liftover_level is not None
    ):
        if not map_dat:
            map_dat = TranscriptMapData(hdp=validator.hdp)

        selected = []

        def accession_on_build(accession, build):
            if build.startswith('GRC'):
                return seq_data.to_chr_num_refseq(accession, build)

            if build.startswith('hg'):
                return seq_data.to_chr_num_ucsc(accession, build)

            return None

        for tx in tx_list:
            options = map_dat.mapping_options(
                tx,
                hdp=validator.hdp,
            )

            for option in options:
                accession = option[1]

                # Primary liftover is restricted to primary chromosomes.
                if liftover_level == 'primary':
                    if not accession.startswith('NC_'):
                        continue

                # Full liftover additionally permits alternate loci.
                elif not accession.startswith(('NC_', 'NT_', 'NW_')):
                    continue

                on_from_build = accession_on_build(
                    accession,
                    build_from,
                )
                on_to_build = accession_on_build(
                    accession,
                    build_to,
                )

                if on_from_build is not None:
                    selected.append(
                        [
                            option[0],
                            accession,
                            build_from,
                            alt_build_from,
                        ]
                    )

                elif on_to_build is not None:
                    selected.append(
                        [
                            option[0],
                            accession,
                            build_to,
                            alt_build_to,
                        ]
                    )

        # Remove duplicate genomic accessions while preserving insertion
        # order and the first usable transcript/accession combination.
        filtered = {}

        for tx, accession, mapped_build, mapped_alt_build in selected:
            if accession not in filtered:
                filtered[accession] = [
                    tx,
                    mapped_build,
                    mapped_alt_build,
                ]

        if filtered:
            added_data = False

            for accession, mapping_data in filtered.items():
                tx, mapped_build, mapped_alt_build = mapping_data
                am_i_gapped = None

                try:
                    # UTA may occasionally identify an overlapping transcript
                    # which cannot actually be mapped at this position.
                    hgvs_tx = validator.vm.g_to_t(
                        hgvs_genomic,
                        tx,
                    )

                    hgvs_alt_genomic = validator.vm.t_to_g(
                        hgvs_tx,
                        accession,
                    )

                    # ------------------------------------------------------
                    # VariantFormatter gap compensation
                    # ------------------------------------------------------

                    if (
                            gap_map is not False
                            and build_from not in mapping_data
                    ):
                        map_to_assembly = None

                        get_assembly = seq_data.supported_for_mapping(
                            accession,
                            "GRCh38",
                        )

                        if get_assembly is True:
                            map_to_assembly = "GRCh38"
                        else:
                            get_assembly = seq_data.supported_for_mapping(
                                accession,
                                "GRCh37",
                            )

                            if get_assembly is True:
                                map_to_assembly = "GRCh37"

                        if map_to_assembly is not None:
                            no_norm_evm = AssemblyMapper(
                                validator.hdp,
                                assembly_name=build_to,
                                alt_aln_method="splign",
                                normalize=False,
                                replace_reference=True,
                            )

                            check_current_alt_genomic = copy.copy(
                                hgvs_alt_genomic
                            )

                            try:
                                am_i_gapped = gap_map(
                                    specified_tx_variant,
                                    hgvs_alt_genomic,
                                    map_to_assembly,
                                    vfo,
                                )

                            except AttributeError:
                                if specified_tx_variant is None:
                                    try:
                                        am_i_gapped = gap_map(
                                            hgvs_tx,
                                            hgvs_alt_genomic,
                                            map_to_assembly,
                                            vfo,
                                        )
                                    except Exception:
                                        pass

                            except Exception:
                                pass

                            if am_i_gapped is not None:
                                hgvs_alt_genomic = (
                                    am_i_gapped["hgvs_genomic"]
                                )

                            else:
                                try:
                                    hgvs_alt_genomic = (
                                        validator.myvm_t_to_g(
                                            specified_tx_variant,
                                            hgvs_alt_genomic.ac,
                                            no_norm_evm,
                                            hn,
                                            map_dat,
                                        )
                                    )

                                except AttributeError:
                                    if specified_tx_variant is None:
                                        hgvs_alt_genomic = (
                                            validator.myvm_t_to_g(
                                                hgvs_tx,
                                                hgvs_alt_genomic.ac,
                                                no_norm_evm,
                                                hn,
                                                map_dat,
                                            )
                                        )

                                if (
                                        check_current_alt_genomic
                                        != hgvs_alt_genomic
                                ):
                                    am_i_gapped = {
                                        "gapped_alignment_warning":
                                            f"Variant "
                                            f"{utils.valstr(hgvs_alt_genomic)} "
                                            f"may be an artefact of alignment "
                                            f"with the selected transcript due "
                                            f"to an alignment gap in the "
                                            f"specified position",
                                        "gap_statement": "",
                                        "gap_position": "",
                                    }

                    # ------------------------------------------------------
                    # Target-build result
                    # ------------------------------------------------------

                    if (
                            mapped_build == build_to
                            or mapped_alt_build == alt_build_to
                    ):
                        alt_vcf = hgvs_utils.report_hgvs2vcf(
                            hgvs_alt_genomic,
                            build_to,
                            reverse_normalizer,
                            validator.sf,
                        )

                        if hgvs_alt_genomic.ac in (
                                'NC_012920.1',
                                'NC_001807.4',
                        ):
                            hgvs_alt_genomic.type = "m"

                        if mapped_build == build_to:
                            lifted_response[
                                build_to.lower()
                            ][hgvs_alt_genomic.ac] = {
                                'hgvs_genomic_description':
                                    hgvs_alt_genomic,
                                'vcf': {
                                    'chr': alt_vcf[to_set],
                                    'pos': str(alt_vcf['pos']),
                                    'ref': alt_vcf['ref'],
                                    'alt': alt_vcf['alt'],
                                },
                            }

                            added_data = True

                        if mapped_alt_build == alt_build_to:
                            lifted_response[
                                alt_build_to.lower()
                            ][hgvs_alt_genomic.ac] = {
                                'hgvs_genomic_description':
                                    hgvs_alt_genomic,
                                'vcf': {
                                    'chr': alt_vcf[alt_to_set],
                                    'pos': str(alt_vcf['pos']),
                                    'ref': alt_vcf['ref'],
                                    'alt': alt_vcf['alt'],
                                },
                            }

                            added_data = True

                    # ------------------------------------------------------
                    # Source-build/PAR information
                    # ------------------------------------------------------

                    if (
                            mapped_build == build_from
                            or mapped_alt_build == alt_build_from
                    ):
                        alt_vcf_from = hgvs_utils.report_hgvs2vcf(
                            hgvs_alt_genomic,
                            build_from,
                            reverse_normalizer,
                            validator.sf,
                        )

                        if mapped_build == build_from:
                            lifted_response[
                                build_from.lower()
                            ][hgvs_alt_genomic.ac] = {
                                'hgvs_genomic_description':
                                    hgvs_alt_genomic,
                                'vcf': {
                                    'chr': alt_vcf_from[to_set],
                                    'pos': str(alt_vcf_from['pos']),
                                    'ref': alt_vcf_from['ref'],
                                    'alt': alt_vcf_from['alt'],
                                },
                            }

                        if mapped_alt_build == alt_build_from:
                            lifted_response[
                                alt_build_from.lower()
                            ][hgvs_alt_genomic.ac] = {
                                'hgvs_genomic_description':
                                    hgvs_alt_genomic,
                                'vcf': {
                                    'chr':
                                        alt_vcf_from[alt_to_set],
                                    'pos':
                                        str(alt_vcf_from['pos']),
                                    'ref':
                                        alt_vcf_from['ref'],
                                    'alt':
                                        alt_vcf_from['alt'],
                                },
                            }

                    if am_i_gapped is not None:
                        lifted_response["am_i_gapped"] = am_i_gapped

                except vvhgvs.exceptions.HGVSError:
                    continue

            # A transcript-mediated target mapping was obtained.
            # Do not invoke PyLiftover.
            if added_data:
                return lifted_response

    # ==================================================================
    # STEP 2: PYLIFTOVER FALLBACK
    # ==================================================================

    genome_builds = [build_to]

    if lo_from not in LO_CACHE:
        LO_CACHE[lo_from] = {}

    if lo_to not in LO_CACHE[lo_from]:
        LO_CACHE[lo_from][lo_to] = LiftOver(
            lo_from,
            lo_to,
        )

    lo = LO_CACHE[lo_from][lo_to]

    # PyLiftover expects UCSC-style chromosome names.
    if from_vcf[from_set].startswith('chr'):
        source_chrom = from_vcf[from_set]
    else:
        source_chrom = 'chr' + from_vcf[from_set]

    liftover_list = lo.convert_coordinate(
        source_chrom,
        int(from_vcf['pos']),
    )

    for lifted in liftover_list:
        chrom = lifted[0]
        pos = lifted[1]
        orientation = lifted[2]

        lifted_ref_bases = from_vcf['ref']
        lifted_alt_bases = from_vcf['alt']

        if hgvs_genomic.posedit.edit.type == "dup":
            # Put the complete original sequence in REF and both copies
            # of the duplicated sequence in ALT.
            lifted_ref_bases += lifted_alt_bases[1:]
            lifted_alt_bases += lifted_alt_bases[1:]

        # Reverse-orientation mappings are not currently supported here.
        if orientation != '+':
            continue

        accession = seq_data.to_accession(
            chrom,
            lo_to,
        )

        if accession is None:
            logger.info(
                'Unable to identify an equivalent %s chromosome ID for %s',
                lo_to,
                chrom,
            )
            continue

        # --------------------------------------------------------------
        # Mitochondrial special cases
        # --------------------------------------------------------------

        # GRCh37 -> GRCh38
        if (
                "38" in build_to
                and "GRCh37" in build_from
                and accession == "NC_012920.1"
        ):
            mito_correction = True
            hgvs_lifted = hgvs_genomic

            if from_vcf[from_set].startswith('chr'):
                chrom = from_vcf[from_set]
            else:
                chrom = 'chr' + from_vcf[from_set]

            pos = int(from_vcf['pos'])

            liftback_list = [
                (chrom, pos, "+", "GRCh38"),
                (chrom, pos, "+", "GRCh37"),
            ]

        # GRCh38 -> GRCh37/hg19
        elif (
                build_to == "GRCh37"
                and "38" in build_from
                and accession == "NC_001807.4"
        ):
            mito_correction = True
            hgvs_lifted = hgvs_genomic

            lifted_data = list(lifted)
            lifted_data[-1] = "hg19"
            lifted = tuple(lifted_data)

            if from_vcf[from_set].startswith('chr'):
                chrom = from_vcf[from_set]
            else:
                chrom = 'chr' + from_vcf[from_set]

            pos = int(from_vcf['pos'])

            liftback_list = [
                (chrom, pos, "+", "GRCh37"),
                lifted,
            ]

            m19_hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                accession,
                'g',
                int(lifted[1]),
                '',
                lifted_alt_bases,
                end=(
                    int(lifted[1])
                    + len(lifted_ref_bases)
                    - 1
                ),
            )

            try:
                m19_hgvs_lifted = hn.normalize(
                    m19_hgvs_not_delins
                )
            except vvhgvs.exceptions.HGVSError:
                m19_hgvs_lifted = None

        # hg19/GRCh37 -> GRCh38 for old mitochondrial accession
        elif (
                build_to == "GRCh38"
                and ("37" in build_from or "19" in build_from)
                and hgvs_genomic.ac == "NC_001807.4"
        ):
            mito_correction = True
            m19_hgvs_lifted = hgvs_genomic

            lifted_data = list(lifted)
            lifted_data[-1] = "GRCh37"
            lifted = tuple(lifted_data)

            if from_vcf[from_set].startswith('chr'):
                chrom = from_vcf[from_set]
            else:
                chrom = 'chr' + from_vcf[from_set]

            pos = int(from_vcf['pos'])

            liftback_list = [
                lifted,
                (chrom, pos, "+", "hg19"),
            ]

            m38_hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                accession,
                'g',
                int(lifted[1]),
                '',
                lifted_alt_bases,
                end=(
                    int(lifted[1])
                    + len(lifted_ref_bases)
                    - 1
                ),
            )

            try:
                hgvs_lifted = hn.normalize(
                    m38_hgvs_not_delins
                )
            except vvhgvs.exceptions.HGVSError:
                continue

        # --------------------------------------------------------------
        # Standard genomic liftover
        # --------------------------------------------------------------

        else:
            mito_correction = False

            hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                accession,
                'g',
                pos,
                '',
                lifted_alt_bases,
                end=pos + len(lifted_ref_bases) - 1,
            )

            try:
                hgvs_lifted = hn.normalize(
                    hgvs_not_delins
                )

            except vvhgvs.exceptions.HGVSDataNotAvailableError:
                continue

            except vvhgvs.exceptions.HGVSInvalidVariantError:
                continue

            # Lift the result back to the original build. The mapping is
            # accepted only when it returns to the original position.
            if lo_to not in LO_CACHE:
                LO_CACHE[lo_to] = {}

            if lo_from not in LO_CACHE[lo_to]:
                LO_CACHE[lo_to][lo_from] = LiftOver(
                    lo_to,
                    lo_from,
                )

            reverse_lo = LO_CACHE[lo_to][lo_from]

            liftback_list = reverse_lo.convert_coordinate(
                chrom,
                pos,
            )

        # --------------------------------------------------------------
        # Validate lift-back
        # --------------------------------------------------------------

        for lifted_back in liftback_list:
            mito_build = False

            if mito_correction and len(lifted_back) > 3:
                mito_build = lifted_back[3]

            if mito_build == 'hg19':
                if m19_hgvs_lifted is None:
                    continue

                hgvs_lifted = m19_hgvs_lifted

            if lifted_back[0].startswith('chr'):
                lifted_back_chrom = lifted_back[0]
            else:
                lifted_back_chrom = 'chr' + lifted_back[0]

            if from_vcf[from_set].startswith('chr'):
                original_chrom = from_vcf[from_set]
            else:
                original_chrom = 'chr' + from_vcf[from_set]

            chromosome_matches = (
                lifted_back_chrom == original_chrom
            )

            position_matches = (
                lifted_back[1] == int(from_vcf['pos'])
                or mito_build in ("hg19", "GRCh37")
            )

            if not chromosome_matches or not position_matches:
                continue

            for build in genome_builds:
                vcf_dict = hgvs_utils.report_hgvs2vcf(
                    hgvs_lifted,
                    build,
                    reverse_normalizer,
                    validator.sf,
                )

                if hgvs_lifted.ac in (
                        'NC_012920.1',
                        'NC_001807.4',
                ):
                    hgvs_lifted.type = "m"

                # ------------------------------------------------------
                # Standard chromosome
                # ------------------------------------------------------

                if mito_build is False:
                    if build.startswith('GRC'):
                        target_chr = vcf_dict['grc_chr']
                        alternate_chr = vcf_dict['ucsc_chr']
                    else:
                        target_chr = vcf_dict['ucsc_chr']
                        alternate_chr = vcf_dict['grc_chr']

                    lifted_response[
                        build_to.lower()
                    ][hgvs_lifted.ac] = {
                        'hgvs_genomic_description': hgvs_lifted,
                        'vcf': {
                            'chr': target_chr,
                            'pos': str(vcf_dict['pos']),
                            'ref': vcf_dict['ref'],
                            'alt': vcf_dict['alt'],
                        },
                    }

                    lifted_response[
                        alt_build_to.lower()
                    ][hgvs_lifted.ac] = {
                        'hgvs_genomic_description': hgvs_lifted,
                        'vcf': {
                            'chr': alternate_chr,
                            'pos': str(vcf_dict['pos']),
                            'ref': vcf_dict['ref'],
                            'alt': vcf_dict['alt'],
                        },
                    }

                # ------------------------------------------------------
                # Mitochondrial corrected result
                # ------------------------------------------------------

                elif mito_build.startswith('GRC'):
                    lifted_response.setdefault(
                        mito_build.lower(),
                        {},
                    )

                    lifted_response[
                        mito_build.lower()
                    ][hgvs_lifted.ac] = {
                        'hgvs_genomic_description': hgvs_lifted,
                        'vcf': {
                            'chr': vcf_dict['grc_chr'],
                            'pos': str(vcf_dict['pos']),
                            'ref': vcf_dict['ref'],
                            'alt': vcf_dict['alt'],
                        },
                    }

                    # Preserve the historical GRCh38/hg38 mitochondrial
                    # representation expected by VariantValidator.
                    lifted_response.setdefault('grch38', {})
                    lifted_response.setdefault('hg38', {})

                    if not lifted_response['grch38']:
                        lifted_response[
                            'grch38'
                        ][hgvs_lifted.ac] = {
                            'hgvs_genomic_description': hgvs_lifted,
                            'vcf': {
                                'chr': vcf_dict['grc_chr'],
                                'pos': str(vcf_dict['pos']),
                                'ref': vcf_dict['ref'],
                                'alt': vcf_dict['alt'],
                            },
                        }

                        lifted_response[
                            'hg38'
                        ][hgvs_lifted.ac] = {
                            'hgvs_genomic_description': hgvs_lifted,
                            'vcf': {
                                'chr': 'chrM',
                                'pos': str(vcf_dict['pos']),
                                'ref': vcf_dict['ref'],
                                'alt': vcf_dict['alt'],
                            },
                        }

                else:
                    lifted_response.setdefault(
                        mito_build.lower(),
                        {},
                    )

                    lifted_response[
                        mito_build.lower()
                    ][hgvs_lifted.ac] = {
                        'hgvs_genomic_description': hgvs_lifted,
                        'vcf': {
                            'chr': vcf_dict['ucsc_chr'],
                            'pos': str(vcf_dict['pos']),
                            'ref': vcf_dict['ref'],
                            'alt': vcf_dict['alt'],
                        },
                    }

    # ------------------------------------------------------------------
    # Remove known invalid mitochondrial combinations
    # ------------------------------------------------------------------

    if (
            'hg19' in lifted_response
            and 'NC_012920.1' in lifted_response['hg19']
    ):
        lifted_response.pop('hg19')

    if (
            'grch37' in lifted_response
            and 'NC_001807.4' in lifted_response['grch37']
    ):
        lifted_response['grch37'].pop('NC_001807.4')

    return lifted_response


# <LICENSE>
# Copyright (C) 2016-2026 VariantValidator Contributors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>
