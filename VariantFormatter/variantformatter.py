# -*- coding: utf-8 -*-

"""
This module creates an initialization object.
This object connects to the hgvs Python library and associated databases.

The Initialization object is used by FormatVariant.
The FormatVariant object contains all HGVS descriptions available for a given
genomic variant, g_to_p.
"""

import collections
import copy
import json
import logging

import vvhgvs.exceptions

import VariantFormatter.formatter as formatter
import VariantValidator.modules.liftover as lo
from VariantValidator.modules import hgvs_utils, seq_data
import VariantValidator.modules.utils as fn


logger = logging.getLogger(
    f"VariantValidator.VariantFormatter."
    f"{__name__.removeprefix('VariantFormatter.')}"
)


# Custom Exceptions
class vcf2hgvsError(Exception):
    pass


class hgvs2VcfError(Exception):
    pass


class variableError(Exception):
    pass


class GenomicDescriptions:
    """
    Object contains genomic level sequence variant descriptions in the pseudo
    VCF (p_vcf), genomic HGVS (g_hgvs) and un-normalized g_hgvs.

    The reference bases are the HGVS description reference nucleotide sequence
    corresponding to the specified range.
    """

    def __init__(
            self,
            p_vcf,
            g_hgvs,
            un_norm_hgvs,
            hgvs_ref_bases,
            gen_error,
            genome_build,
            variant_description,
    ):
        if p_vcf == "None":
            p_vcf = None

        try:
            if g_hgvs == "None":
                g_hgvs = None
            elif (
                    g_hgvs.ac in ("NC_012920.1", "NC_001807.4")
                    and ":g." in variant_description
            ):
                gen_error = (
                    f"VariantTypeError: The given reference sequence "
                    f"({g_hgvs.ac}) does not match the DNA type (g). "
                    f"For {g_hgvs.ac}, please use (m). For g. variants, "
                    "please use a linear genomic reference sequence"
                )
        except AttributeError:
            if g_hgvs == "None":
                g_hgvs = None

        if un_norm_hgvs == "None":
            un_norm_hgvs = None

        if hgvs_ref_bases == "None":
            hgvs_ref_bases = None

        if gen_error == "None":
            gen_error = None

        try:
            accession = g_hgvs.ac
        except AttributeError:
            accession = None

        # Do not overwrite a more fundamental error already identified above.
        if gen_error is None:
            if (
                    accession == "NC_012920.1"
                    and genome_build == "hg19"
            ):
                gen_error = (
                    "GenomeBuildError: NC_012920.1 is not associated "
                    "with genome build hg19, instead use genome build "
                    "GRCh37"
                )
            elif (
                    accession == "NC_001807.4"
                    and genome_build == "GRCh37"
            ):
                gen_error = (
                    "GenomeBuildError: NC_001807.4 is not associated "
                    "with genome build GRCh37, instead use genome build "
                    "hg19"
                )

        self.p_vcf = p_vcf

        try:
            self.g_hgvs = formatter.remove_reference(g_hgvs)
        except AttributeError:
            self.g_hgvs = None

        self.un_norm_hgvs = un_norm_hgvs
        self.g_hgvs_ref = hgvs_ref_bases
        self.gen_error = gen_error
        self.gen_warnings = None
        self.selected_build = genome_build


class FormatVariant:
    """
    Return genomic, transcript and protein descriptions for a genomic variant.

    Parameters
    ----------
    variant_description : str
        Genomic HGVS or VCF-like variant description.
    genome_build : str
        Genome build.
    vfo
        VariantValidator object used by VariantFormatter.
    transcript_model : str, optional
        Transcript source selection.
    specify_transcripts : optional
        Transcript-selection mode or specified transcript accession(s).
    checkOnly : bool or str, optional
        Restrict processing where requested.
    liftover : bool, optional
        Enable liftover.
    legacy_genomic_structure : bool, optional
        Preserve the historical VariantFormatter genomic loci structure.
        Defaults to True for backwards compatibility. When False, genomic
        loci are returned using the VariantValidator structure.
    """

    def __init__(
        self,
        variant_description,
        genome_build,
        vfo,
        transcript_model=None,
        specify_transcripts=None,
        checkOnly=False,
        liftover=False,
        legacy_genomic_structure=True,
    ):
        self.variant_description = variant_description
        self.vfo = vfo
        self.warning_level = None
        self.liftover = liftover
        self.legacy_genomic_structure = legacy_genomic_structure
        self.direct_reformatting = {
            "instance": None,
            "reformat": None,
        }

        gen_error = None
        recovery_error = None

        if genome_build not in ("GRCh37", "GRCh38", "hg19", "hg38"):
            self.genomic_descriptions = GenomicDescriptions(
                None,
                None,
                None,
                None,
                (
                    "genome_build must be one of: "
                    "'GRCh37'; 'GRCh38'; 'hg19'; 'hg38'"
                ),
                genome_build,
                variant_description,
            )
            self.warning_level = "genomic_variant_warning"
            return

        self.genome_build = genome_build
        vfo.genome_build = genome_build

        if transcript_model is None:
            transcript_model = "all"

        if transcript_model not in ("ensembl", "refseq", "all"):
            self.genomic_descriptions = GenomicDescriptions(
                None,
                None,
                None,
                None,
                (
                    "transcript_model must be one of: "
                    "'ensembl'; 'refseq'; 'all'"
                ),
                genome_build,
                variant_description,
            )
            self.warning_level = "genomic_variant_warning"
            return

        self.transcript_model = transcript_model
        self.specify_transcripts = specify_transcripts

        # --------------------------------------------------------------
        # HGVS genomic input
        # --------------------------------------------------------------
        if self.variant_description.startswith(("NC_", "NT_", "NW_")):
            to_process, separator, reformat = variant_description.partition("|")

            if separator and reformat in ("l", "g", "m"):
                self.direct_reformatting["reformat"] = reformat
                self.direct_reformatting["instance"] = "methylation"
                variant_description = f"{to_process}="
                self.variant_description = variant_description

            try:
                hgvs_genomic = formatter.parse(
                    self.variant_description,
                    self.vfo,
                )
                vfo.vr.validate(hgvs_genomic)

            except Exception:
                validation = vfo.validate(
                    self.variant_description,
                    self.genome_build,
                    "all",
                    liftover_level=None,
                ).format_as_dict(test=True)

                reset_variant = None
                edit_warnings = None

                for val_val in validation.values():
                    try:
                        if "primary_assembly_loci" not in val_val:
                            continue

                        reset_variant = val_val[
                            "primary_assembly_loci"
                        ][genome_build.lower()][
                            "hgvs_genomic_description"
                        ]

                        validation_warned = val_val[
                            "validation_warnings"
                        ]

                        edit_warnings = [
                            warning
                            for warning in validation_warned
                            if "automapped to" in warning
                        ]
                        break

                    except AttributeError:
                        continue

                    except KeyError:
                        validation_warned = val_val[
                            "validation_warnings"
                        ]

                        edit_warnings = ", ".join(
                            validation_warned
                        )

                        self.genomic_descriptions = GenomicDescriptions(
                            None,
                            None,
                            None,
                            None,
                            edit_warnings,
                            genome_build,
                            variant_description,
                        )
                        self.warning_level = "genomic_variant_warning"
                        return

                if reset_variant is None:
                    self.genomic_descriptions = GenomicDescriptions(
                        None,
                        None,
                        None,
                        None,
                        (
                            f"InvalidSyntaxError: Unable to obtain a valid genomic HGVS description "
                            f"for {self.variant_description}"
                        ),
                        genome_build,
                        variant_description,
                    )
                    self.warning_level = "genomic_variant_warning"
                    return

                hgvs_genomic = formatter.parse(
                    reset_variant,
                    self.vfo,
                )

                if edit_warnings:
                    recovery_error = edit_warnings[0]

                self.warning_level = "genomic_variant_warning"

            # Check the parsed genomic reference against the selected build.
            #
            # This is deliberately performed after the normal and recovery
            # parsing paths converge so that both operate on an HGVS object.
            if self.genome_build.lower().startswith("grch"):
                seq_data_func = seq_data.to_chr_num_refseq
            else:
                seq_data_func = seq_data.to_chr_num_ucsc
            if seq_data_func(
                    hgvs_genomic.ac,
                    self.genome_build
            ) is None:
                gen_error = (
                    f"GenomeBuildError: chromosome ID {hgvs_genomic.ac} "
                    f"is not associated with genome build {self.genome_build}"
                )

                self.genomic_descriptions = GenomicDescriptions(
                    None,
                    hgvs_genomic,
                    None,
                    None,
                    gen_error,
                    genome_build,
                    variant_description,
                )
                self.warning_level = "genomic_variant_warning"
                return

            # Preserve the submitted/recovered HGVS object before normalization.
            un_norm_hgvs = copy.deepcopy(
                hgvs_genomic
            )

            # Normalize the HGVS object directly.
            #
            # Do not round-trip HGVS -> VCF -> HGVS. The VCF representation is
            # output only and must not become the source for reconstruction of an
            # already valid HGVS object.
            try:
                g_hgvs = vfo.splign_normalizer.normalize(
                    hgvs_genomic
                )

            except vvhgvs.exceptions.HGVSError as error:
                self.genomic_descriptions = GenomicDescriptions(
                    None,
                    None,
                    None,
                    None,
                    str(error),
                    genome_build,
                    variant_description,
                )
                self.warning_level = "genomic_variant_warning"
                return

            # Fetch reference sequence directly from the normalized HGVS object.
            hgvs_ref_bases = self.vfo.sf.fetch_seq(
                g_hgvs.ac,
                start_i=g_hgvs.posedit.pos.start.base - 1,
                end_i=g_hgvs.posedit.pos.end.base,
            )

            # Generate pseudo-VCF independently for output.
            try:
                vcf_dictionary = formatter.hgvs_genomic2vcf(
                    g_hgvs,
                    self.genome_build,
                    self.vfo,
                )

                if (
                        vcf_dictionary["grc_chr"] == "NC_001807.4"
                        and genome_build == "hg19"
                ):
                    chromosome = vcf_dictionary["ucsc_chr"]
                else:
                    chromosome = vcf_dictionary["grc_chr"]

                p_vcf = ":".join(
                    (
                        chromosome,
                        vcf_dictionary["pos"],
                        vcf_dictionary["ref"],
                        vcf_dictionary["alt"],
                    )
                )

            except Exception as error:
                error_message = str(error)

                if "Variant span is outside sequence bounds" in error_message:
                    accession = variant_description.split(":", 1)[0]
                    error_message = (
                        "The specified coordinate is outside the "
                        f"boundaries of reference sequence {accession}"
                    )

                try:
                    if "N" in g_hgvs.posedit.edit.ref:
                        error_message = (
                            "UncertainSequenceError: The submitted variant "
                            f"description "
                            f"{formatter.remove_reference(g_hgvs)} "
                            "refers to a genomic reference region with an "
                            "uncertain base composition (N)"
                        )
                except AttributeError:
                    pass

                self.genomic_descriptions = GenomicDescriptions(
                    None,
                    None,
                    None,
                    None,
                    error_message,
                    genome_build,
                    variant_description,
                )
                self.warning_level = "genomic_variant_warning"
                return

            # Report normalization only when the genomic coordinates move.
            if (
                    g_hgvs.posedit.pos.start.base
                    != un_norm_hgvs.posedit.pos.start.base
                    or g_hgvs.posedit.pos.end.base
                    != un_norm_hgvs.posedit.pos.end.base
            ):
                self.warning_level = "genomic_variant_warning"
                gen_error = (
                    f"AutoCorrectionWarning: {self.variant_description} updated to "
                    f"{formatter.remove_reference(g_hgvs)}"
                )

        # --------------------------------------------------------------
        # Unsupported HGVS reference types
        # --------------------------------------------------------------
        elif self.variant_description.startswith(("NG_", "LRG_")):
            if self.variant_description.startswith("NG_"):
                gen_error = (
                    f"UnsupportedFormatError: Variant description {self.variant_description} uses "
                    "the NG_ reference type, this is currently not "
                    "accepted through this tool"
                )
            else:
                gen_error = (
                    f"UnsupportedFormatError: Variant description {self.variant_description} uses "
                    "the LRG_ reference type, LRGs are no longer being "
                    "updated, and are not recommended, they are also not "
                    "currently accepted through this tool"
                )

            self.genomic_descriptions = GenomicDescriptions(
                None,
                None,
                None,
                None,
                gen_error,
                genome_build,
                variant_description,
            )
            self.warning_level = "submission_warning"
            return

        # --------------------------------------------------------------
        # VCF-like input
        # --------------------------------------------------------------
        elif self._is_vcf_like_description(
            self.variant_description
        ):
            try:
                genomic_level = formatter.vcf2hgvs_genomic(
                    self.variant_description,
                    self.genome_build,
                    self.vfo,
                )
            except Exception as error:
                raise vcf2hgvsError(str(error)) from error

            if genomic_level["error"]:
                self.genomic_descriptions = GenomicDescriptions(
                    None,
                    None,
                    None,
                    None,
                    genomic_level["error"],
                    genome_build,
                    variant_description,
                )
                self.warning_level = "genomic_variant_warning"
                return

            vcf_dictionary = formatter.hgvs_genomic2vcf(
                genomic_level["hgvs_genomic"],
                self.genome_build,
                self.vfo,
            )

            if (
                vcf_dictionary["grc_chr"] == "NC_001807.4"
                and genome_build == "hg19"
            ):
                chromosome = vcf_dictionary["ucsc_chr"]
            else:
                chromosome = vcf_dictionary["grc_chr"]

            p_vcf = "-".join(
                (
                    chromosome,
                    vcf_dictionary["pos"],
                    vcf_dictionary["ref"],
                    vcf_dictionary["alt"],
                )
            )

            g_hgvs = genomic_level["hgvs_genomic"]
            un_norm_hgvs = genomic_level[
                "un_normalized_hgvs_genomic"
            ]
            hgvs_ref_bases = genomic_level["ref_bases"]

        else:
            gen_error = (
                f"UnsupportedFormatError: Variant description {self.variant_description} is not "
                "in a supported format. This tool accepts vcf-like and "
                "HGVS genomic (g.) descriptions only"
            )

            self.genomic_descriptions = GenomicDescriptions(
                None,
                None,
                None,
                None,
                gen_error,
                genome_build,
                variant_description,
            )
            self.warning_level = "submission_warning"
            return

        if recovery_error is not None:
            gen_error = recovery_error

        self.genomic_descriptions = GenomicDescriptions(
            p_vcf,
            g_hgvs,
            un_norm_hgvs,
            hgvs_ref_bases,
            gen_error,
            genome_build,
            variant_description,
        )

        if checkOnly is True:
            return

        prelim_transcript_descriptions = {}

        if self.genome_build == "hg19":
            self.genome_build = "GRCh37"
        elif self.genome_build == "hg38":
            self.genome_build = "GRCh38"

        # --------------------------------------------------------------
        # Transcript selection
        # --------------------------------------------------------------
        if (
            self.specify_transcripts is not None
            and "select" not in self.specify_transcripts
            and "mane" not in self.specify_transcripts
            and "raw" not in self.specify_transcripts
        ):
            try:
                trans_list = json.loads(
                    self.specify_transcripts
                )
            except (json.decoder.JSONDecodeError, TypeError):
                trans_list = [
                    self.specify_transcripts
                ]

            transcript_list = [
                [transcript, ""]
                for transcript in trans_list
            ]

        else:
            transcript_list = formatter.fetch_aligned_transcripts(
                g_hgvs,
                self.transcript_model,
                self.vfo,
                genome_build,
            )

        transcript_dict = {}

        for tx in transcript_list:
            if "/" in tx[0]:
                continue

            transcript_dict.setdefault(tx[0], 1)

            if (
                len(tx) > 2
                and (
                    tx[1].startswith("NC_00")
                    or tx[1] in (
                        "NC_012920.1",
                        "NC_001807.4",
                    )
                )
            ):
                transcript_dict[tx[0]] += 1

        transcript_list = sorted(
            transcript_dict,
            key=transcript_dict.get,
            reverse=True,
        )

        g_to_g_lift = {}

        # --------------------------------------------------------------
        # Transcript processing
        # --------------------------------------------------------------
        for tx_id in transcript_list:
            try:
                annotation = vfo.db.get_transcript_annotation(
                    tx_id
                )
                annotation_dict = json.loads(annotation)

                gene_symbol = (
                    vfo.db.get_gene_symbol_from_transcript_id(
                        tx_id
                    )
                )

                gene_dict = {
                    "symbol": gene_symbol,
                    "hgnc_id": annotation_dict[
                        "db_xref"
                    ]["hgnc"],
                }

            except (
                json.decoder.JSONDecodeError,
                KeyError,
            ):
                continue

            select_dict = {
                key: True
                for key, value in annotation_dict.items()
                if value in ("true", True)
            }

            if self.specify_transcripts == "select":
                select_dict = {
                    key: value
                    for key, value in select_dict.items()
                    if "select" in key
                }

                if not select_dict:
                    continue

            elif self.specify_transcripts == "mane_select":
                select_dict = {
                    key: value
                    for key, value in select_dict.items()
                    if key == "mane_select"
                }

                if not select_dict:
                    continue

            elif self.specify_transcripts == "mane":
                select_dict = {
                    key: value
                    for key, value in select_dict.items()
                    if "mane" in key
                }

                if not select_dict:
                    continue

            if (
                self.specify_transcripts
                and any(
                    accession_type in self.specify_transcripts
                    for accession_type in (
                        "NM_",
                        "NR_",
                        "ENST",
                    )
                )
            ):
                overlapping_tx = (
                    formatter.fetch_aligned_transcripts(
                        g_hgvs,
                        self.transcript_model,
                        self.vfo,
                        genome_build,
                    )
                )

                overlapping_accessions = {
                    transcript[0]
                    for transcript in overlapping_tx
                }

                if tx_id not in overlapping_accessions:
                    continue

            hgvs_transcript_dict = (
                formatter.hgvs_genomic2hgvs_transcript(
                    g_hgvs,
                    tx_id,
                    self.vfo,
                )
            )

            hgvs_transcript = None
            hgvs_protein_tlc = None
            select_status = None

            logger.info(
                "VF transcript mapping: tx_id=%s, hgvs_transcript=%s, "
                "object_type=%s, error=%r",
                tx_id,
                hgvs_transcript_dict["hgvs_transcript"],
                type(hgvs_transcript_dict["hgvs_transcript"]).__name__,
                hgvs_transcript_dict["error"],
            )

            try:
                logger.info(
                    "VF before gap_checker: tx_id=%s, hgvs_transcript=%s, "
                    "edit_type=%s, ref=%r, alt=%r",
                    tx_id,
                    hgvs_transcript_dict["hgvs_transcript"],
                    getattr(
                        hgvs_transcript_dict["hgvs_transcript"].posedit.edit,
                        "type",
                        None,
                    ),
                    getattr(
                        hgvs_transcript_dict["hgvs_transcript"].posedit.edit,
                        "ref",
                        None,
                    ),
                    getattr(
                        hgvs_transcript_dict["hgvs_transcript"].posedit.edit,
                        "alt",
                        None,
                    ),
                )

                am_i_gapped = formatter.gap_checker(
                    hgvs_transcript_dict[
                        "hgvs_transcript"
                    ],
                    g_hgvs,
                    self.genome_build,
                    self.vfo,
                    transcript_model=self.transcript_model,
                )

            except Exception:
                logger.exception(
                    "VF gap_checker FAILED: tx_id=%s, input=%s",
                    tx_id,
                    hgvs_transcript_dict["hgvs_transcript"],
                )

                self.warning_level = "processing_error"

                if hgvs_transcript_dict["error"] == "":
                    hgvs_transcript_dict["error"] = None

                am_i_gapped = {
                    "hgvs_transcript": None,
                    "position_lock": False,
                    "gapped_alignment_warning": None,
                    "corrective_action": None,
                    "gap_position": None,
                    "transcript_accession": tx_id,
                    "error": hgvs_transcript_dict["error"],
                    "select_status": None,
                }

            else:
                if hgvs_transcript_dict["error"] == "":
                    hgvs_transcript_dict["error"] = None

                logger.info(
                    "VF after gap_checker: tx_id=%s, hgvs_transcript=%s, "
                    "object_type=%s",
                    tx_id,
                    am_i_gapped["hgvs_transcript"],
                    type(am_i_gapped["hgvs_transcript"]).__name__,
                )

                # Preserve the HGVS SequenceVariant returned by gap_checker.
                hgvs_transcript = am_i_gapped[
                    "hgvs_transcript"
                ]

                logger.info(
                    "VF selected transcript object: tx_id=%s, "
                    "hgvs_transcript=%s, hgvs_type=%s, edit_type=%s, "
                    "ref=%r, alt=%r",
                    tx_id,
                    hgvs_transcript,
                    getattr(hgvs_transcript, "type", None),
                    getattr(
                        hgvs_transcript.posedit.edit,
                        "type",
                        None,
                    ),
                    getattr(
                        hgvs_transcript.posedit.edit,
                        "ref",
                        None,
                    ),
                    getattr(
                        hgvs_transcript.posedit.edit,
                        "alt",
                        None,
                    ),
                )

                if (
                    isinstance(checkOnly, str)
                    and "tx" in checkOnly
                ):
                    hgvs_protein_tlc = None

                elif hgvs_transcript.type == "c":
                    try:
                        # Preserve the protein HGVS object.
                        hgvs_protein_tlc = (
                            formatter.hgvs_transcript2hgvs_protein(
                                hgvs_transcript,
                                self.genome_build,
                                self.vfo,
                            )
                        )

                        logger.info(
                            "VF protein mapping: tx_id=%s, protein=%s, "
                            "object_type=%s",
                            tx_id,
                            hgvs_protein_tlc,
                            type(hgvs_protein_tlc).__name__,
                        )

                    except (
                        NotImplementedError,
                        vvhgvs.exceptions.HGVSDataNotAvailableError,
                    ) as error:
                        logger.info(
                            "VF protein mapping failed: tx_id=%s, error=%s",
                            tx_id,
                            error,
                        )

                        hgvs_protein_tlc = None
                        hgvs_transcript_dict["error"] = str(
                            error
                        )

                am_i_gapped["error"] = (
                    hgvs_transcript_dict["error"]
                )
                select_status = select_dict

            if am_i_gapped["error"] == "":
                am_i_gapped["error"] = None

            # ----------------------------------------------------------
            # Output boundary
            #
            # Keep the canonical HGVS objects untouched. Make copies for
            # output-specific reference removal and serialize only here.
            # ----------------------------------------------------------
            order_my_tp = collections.OrderedDict()

            if hgvs_transcript is not None:
                logger.info(
                    "VF output before deepcopy/unset: tx_id=%s, "
                    "hgvs_transcript=%s, edit_type=%s, ref=%r, alt=%r",
                    tx_id,
                    hgvs_transcript,
                    getattr(
                        hgvs_transcript.posedit.edit,
                        "type",
                        None,
                    ),
                    getattr(
                        hgvs_transcript.posedit.edit,
                        "ref",
                        None,
                    ),
                    getattr(
                        hgvs_transcript.posedit.edit,
                        "alt",
                        None,
                    ),
                )

                output_hgvs_transcript = copy.deepcopy(
                    hgvs_transcript
                )
                output_hgvs_transcript = (
                    hgvs_utils.unset_hgvs_obj_ref(
                        output_hgvs_transcript, vf_mode=True
                    )
                )

                logger.info(
                    "VF output after unset: tx_id=%s, "
                    "output_hgvs_transcript=%s, edit_type=%s, "
                    "ref=%r, alt=%r",
                    tx_id,
                    output_hgvs_transcript,
                    getattr(
                        output_hgvs_transcript.posedit.edit,
                        "type",
                        None,
                    ),
                    getattr(
                        output_hgvs_transcript.posedit.edit,
                        "ref",
                        None,
                    ),
                    getattr(
                        output_hgvs_transcript.posedit.edit,
                        "alt",
                        None,
                    ),
                )

                order_my_tp["t_hgvs"] = fn.valstr(
                    output_hgvs_transcript
                )

                logger.info(
                    "VF final t_hgvs: tx_id=%s, t_hgvs=%r",
                    tx_id,
                    order_my_tp["t_hgvs"],
                )
            else:
                logger.info(
                    "VF output has no transcript: tx_id=%s",
                    tx_id,
                )
                order_my_tp["t_hgvs"] = None

            if hgvs_protein_tlc is not None:
                output_hgvs_protein = copy.deepcopy(
                    hgvs_protein_tlc
                )
                output_hgvs_protein = (
                    hgvs_utils.unset_hgvs_obj_ref(
                        output_hgvs_protein, vf_mode=True
                    )
                )

                order_my_tp["p_hgvs_tlc"] = fn.valstr(
                    output_hgvs_protein
                )

                order_my_tp["p_hgvs_slc"] = (
                    output_hgvs_protein.format(
                        {"p_3_letter": False}
                    )
                )
            else:
                order_my_tp["p_hgvs_tlc"] = None
                order_my_tp["p_hgvs_slc"] = None

            order_my_tp["select_status"] = (
                select_status
            )

            order_my_tp["gene_info"] = gene_dict

            order_my_tp["transcript_version_warning"] = (
                hgvs_transcript_dict["latest_version"]
            )

            order_my_tp["gapped_alignment_warning"] = (
                am_i_gapped["gapped_alignment_warning"]
            )

            order_my_tp["gap_statement"] = (
                am_i_gapped["gap_position"]
            )

            order_my_tp["transcript_variant_error"] = (
                am_i_gapped["error"]
            )

            # ----------------------------------------------------------
            # Liftover
            # ----------------------------------------------------------
            if self.liftover is not False:
                if (
                    self.genomic_descriptions.selected_build
                    in ("hg19", "GRCh37")
                ):
                    build_to = "GRCh38"
                else:
                    build_to = "GRCh37"

                if (
                    (
                        order_my_tp[
                            "transcript_variant_error"
                        ]
                        is not None
                        and not g_to_g_lift
                    )
                    or order_my_tp[
                        "transcript_variant_error"
                    ]
                    is None
                ):
                    # Reuse the HGVS object already produced above.
                    specified_tx_variant = hgvs_transcript

                    current_lift = lo.liftover(
                        self.genomic_descriptions.g_hgvs,
                        self.genomic_descriptions.selected_build,
                        build_to,
                        vfo.splign_normalizer,
                        vfo.reverse_splign_normalizer,
                        None,
                        vfo,
                        specify_tx=tx_id,
                        liftover_level=self.liftover,
                        gap_map=formatter.gap_checker,
                        vfo=self.vfo,
                        specified_tx_variant=specified_tx_variant,
                    )

                    if not current_lift[
                        build_to.lower()
                    ]:
                        direct_lift = lo.liftover(
                            self.genomic_descriptions.g_hgvs,
                            self.genomic_descriptions.selected_build,
                            build_to,
                            vfo.splign_normalizer,
                            vfo.reverse_splign_normalizer,
                            None,
                            vfo,
                            specify_tx=tx_id,
                            liftover_level=self.liftover,
                            gap_map=formatter.gap_checker,
                            vfo=self.vfo,
                            specified_tx_variant=specified_tx_variant,
                            force_pyliftover=True,
                        )

                        current_lift[
                            build_to.lower()
                        ] = direct_lift[
                            build_to.lower()
                        ]

                        if build_to == "GRCh37":
                            current_lift[
                                "hg19"
                            ] = direct_lift[
                                build_to.lower()
                            ]

                        elif build_to == "GRCh38":
                            current_lift[
                                "hg38"
                            ] = direct_lift[
                                build_to.lower()
                            ]

                    if "am_i_gapped" in current_lift:
                        lifted_gap = current_lift.pop(
                            "am_i_gapped"
                        )

                        if (
                            order_my_tp[
                                "gapped_alignment_warning"
                            ]
                            == ""
                        ):
                            order_my_tp[
                                "gapped_alignment_warning"
                            ] = lifted_gap[
                                "gapped_alignment_warning"
                            ]

                        if (
                            order_my_tp["gap_statement"]
                            == ""
                        ):
                            order_my_tp[
                                "gap_statement"
                            ] = lifted_gap[
                                "gap_position"
                            ]

                    if not g_to_g_lift:
                        g_to_g_lift = current_lift

                else:
                    current_lift = g_to_g_lift

                self._format_liftover_descriptions(
                    current_lift
                )

                primary_loci, alt_loci = (
                    self._split_liftover_loci(
                        current_lift
                    )
                )

                self._store_loci(
                    order_my_tp,
                    primary_loci,
                    alt_loci,
                )

            # ----------------------------------------------------------
            # No liftover
            # ----------------------------------------------------------
            else:
                try:
                    chromosome, pos, ref, alt = (
                        self.genomic_descriptions.p_vcf.split(
                            ":"
                        )
                    )
                except ValueError:
                    chromosome, pos, ref, alt = (
                        self.genomic_descriptions.p_vcf.split(
                            "-"
                        )
                    )

                accession = (
                    self.genomic_descriptions.g_hgvs.split(
                        ":",
                        1,
                    )[0]
                )

                primary = {
                    accession: {
                        "hgvs_genomic_description":
                            self.genomic_descriptions.g_hgvs,
                        "vcf": {
                            "chr": chromosome,
                            "pos": pos,
                            "ref": ref,
                            "alt": alt,
                        },
                    }
                }

                ucsc_chromosome = (
                    chromosome
                    if chromosome.startswith("chr")
                    else f"chr{chromosome}"
                )

                ucsc = {
                    accession: {
                        "hgvs_genomic_description":
                            self.genomic_descriptions.g_hgvs,
                        "vcf": {
                            "chr": ucsc_chromosome,
                            "pos": pos,
                            "ref": ref,
                            "alt": alt,
                        },
                    }
                }

                if (
                    self.genomic_descriptions.selected_build
                    == "GRCh38"
                ):
                    primary_loci = {
                        "grch38": primary,
                        "hg38": ucsc,
                    }
                else:
                    primary_loci = {
                        "grch37": primary,
                        "hg19": ucsc,
                    }

                self._store_loci(
                    order_my_tp,
                    primary_loci,
                    [],
                )

            prelim_transcript_descriptions[
                tx_id
            ] = order_my_tp

        self.t_and_p_descriptions = (
            prelim_transcript_descriptions
        )

    @staticmethod
    def _is_vcf_like_description(description):
        """
        Return True when input begins with a VCF-style reference token.
        """
        if not description:
            return False

        colon_index = description.find(":")
        hyphen_index = description.find("-")

        delimiter_indexes = [
            index
            for index in (
                colon_index,
                hyphen_index,
            )
            if index > 0
        ]

        if not delimiter_indexes:
            return False

        reference = description[
            :min(delimiter_indexes)
        ]

        if reference.lower().startswith("chr"):
            reference = reference[3:]

        return (
            bool(reference)
            and reference.replace("_", "").isalnum()
        )

    @staticmethod
    def _format_liftover_descriptions(
        current_lift,
    ):
        """
        Convert liftover HGVS descriptions to strings at the output boundary.
        """
        for mappings in current_lift.values():
            for locus in mappings.values():
                hgvs_description = locus[
                    "hgvs_genomic_description"
                ]

                if hasattr(hgvs_description, "format"):
                    locus[
                        "hgvs_genomic_description"
                    ] = hgvs_description.format(
                        {"max_ref_length": 0}
                    )

    def _split_liftover_loci(
        self,
        current_lift,
    ):
        """
        Split liftover results into primary and alternate loci.

        PAR chrY mappings are moved to alternate loci while chrX remains
        in primary_assembly_loci.
        """
        primary_loci = copy.deepcopy(
            current_lift
        )
        alt_loci = []

        for build, mappings in current_lift.items():
            is_par = (
                any(
                    "23." in accession
                    for accession in mappings
                )
                and any(
                    "24." in accession
                    for accession in mappings
                )
            )

            if is_par:
                self.genomic_descriptions.gen_warnings = (
                    "ParRegionWarning: Variant is located in a "
                    "pseudoautosomal region (PAR) of the X and Y "
                    "chromosomes, so the Y context description has "
                    "been moved to alt_genomic_loci"
                )

            for accession, locus in mappings.items():
                if accession.startswith("NC_"):
                    if (
                        is_par
                        and "24." in accession
                    ):
                        primary_loci[build].pop(
                            accession,
                            None,
                        )
                        alt_loci.append(
                            {build: locus}
                        )

                    continue

                primary_loci[build].pop(
                    accession,
                    None,
                )

                alt_loci.append(
                    {build: locus}
                )

        return primary_loci, alt_loci

    @staticmethod
    def _variantvalidator_primary_loci(
        primary_loci,
    ):
        """
        Flatten historical VariantFormatter primary loci into the
        VariantValidator primary_assembly_loci structure.
        """
        vv_primary = {}

        for build, mappings in primary_loci.items():
            if mappings:
                vv_primary[build] = next(
                    iter(mappings.values())
                )

        return vv_primary

    def _store_loci(
        self,
        output,
        primary_loci,
        alt_loci,
    ):
        """
        Store transcript-associated genomic loci.

        Legacy VariantFormatter structure is the default for backwards
        compatibility. VariantValidator structure is selected explicitly
        with legacy_genomic_structure=False.
        """
        if self.legacy_genomic_structure:
            output[
                "primary_assembly_loci"
            ] = primary_loci
        else:
            output[
                "primary_assembly_loci"
            ] = self._variantvalidator_primary_loci(
                primary_loci
            )

        output[
            "alt_genomic_loci"
        ] = alt_loci

    @staticmethod
    def _legacy_intergenic_alt_loci(
        alt_loci,
    ):
        """
        Restore the historical intergenic VariantFormatter alternate-loci
        dictionary from the internally normalised alternate-loci list.
        """
        legacy_alt_loci = {}

        for build_locus in alt_loci:
            for build, locus in build_locus.items():
                accession = locus[
                    "hgvs_genomic_description"
                ].split(":", 1)[0]

                legacy_alt_loci.setdefault(
                    build,
                    {},
                )[accession] = locus

        return legacy_alt_loci

    def _store_intergenic_loci(
        self,
        output,
        primary_loci,
        alt_loci,
    ):
        """
        Store intergenic genomic loci while preserving the historical
        VariantFormatter alternate-loci dictionary when legacy output is
        requested.
        """
        if self.legacy_genomic_structure:
            output[
                "primary_assembly_loci"
            ] = primary_loci

            output[
                "alt_genomic_loci"
            ] = self._legacy_intergenic_alt_loci(
                alt_loci
            )
            return

        output[
            "primary_assembly_loci"
        ] = self._variantvalidator_primary_loci(
            primary_loci
        )
        output[
            "alt_genomic_loci"
        ] = alt_loci

    def stucture_data(self):
        bring_order = collections.OrderedDict()

        bring_order[
            "p_vcf"
        ] = self.genomic_descriptions.p_vcf

        if self.genomic_descriptions.g_hgvs is not None:
            if self.genomic_descriptions.g_hgvs.startswith(
                ("NC_012920.1:", "NC_001807.4:")
            ):
                self.genomic_descriptions.g_hgvs = (
                    self.genomic_descriptions.g_hgvs.replace(
                        ":g.",
                        ":m.",
                    )
                )

                if self.genomic_descriptions.gen_error:
                    self.genomic_descriptions.gen_error = (
                        self.genomic_descriptions.gen_error.replace(
                            ":g.",
                            ":m.",
                        )
                    )

        bring_order[
            "g_hgvs"
        ] = self.genomic_descriptions.g_hgvs

        bring_order[
            "selected_build"
        ] = self.genomic_descriptions.selected_build

        bring_order[
            "genomic_variant_error"
        ] = self.genomic_descriptions.gen_error

        bring_order[
            "genomic_variant_warnings"
        ] = self.genomic_descriptions.gen_warnings

        try:
            if not self.t_and_p_descriptions:
                if self.legacy_genomic_structure:
                    intergenic = {
                        "alt_genomic_loci": None,
                    }
                else:
                    intergenic = {
                        "primary_assembly_loci": None,
                        "alt_genomic_loci": None,
                    }

                bring_order[
                    "hgvs_t_and_p"
                ] = {
                    "intergenic": intergenic
                }

                if self.liftover is not False:
                    if (
                        self.genomic_descriptions.selected_build
                        in ("hg19", "GRCh37")
                    ):
                        build_to = "GRCh38"
                    else:
                        build_to = "GRCh37"

                    current_lift = lo.liftover(
                        self.genomic_descriptions.g_hgvs,
                        self.genomic_descriptions.selected_build,
                        build_to,
                        self.vfo.splign_normalizer,
                        self.vfo.reverse_splign_normalizer,
                        None,
                        self.vfo,
                        specify_tx=False,
                        liftover_level=self.liftover,
                    )

                    self._format_liftover_descriptions(
                        current_lift
                    )

                    primary_loci, alt_loci = (
                        self._split_liftover_loci(
                            current_lift
                        )
                    )

                    self._store_intergenic_loci(
                        intergenic,
                        primary_loci,
                        alt_loci,
                    )

            else:
                bring_order[
                    "hgvs_t_and_p"
                ] = self.t_and_p_descriptions

        except AttributeError:
            bring_order[
                "hgvs_t_and_p"
            ] = None

        brought_order = {
            self.variant_description: bring_order
        }

        if (
            self.direct_reformatting["instance"]
            == "methylation"
        ):
            replace_json = json.dumps(
                brought_order
            )

            replace_json = replace_json.replace(
                '="',
                (
                    f'|'
                    f'{self.direct_reformatting["reformat"]}'
                    f'"'
                ),
            )

            brought_order = json.loads(
                replace_json
            )

        return brought_order

    def collect_metadata(self):
        meta = collections.OrderedDict()
        meta["api_version"] = self.vfo.version
        meta["hgvs_version"] = self.vfo.hgvsVersion
        meta["uta_schema"] = self.vfo.utaVersion
        meta["seqrepo_db"] = self.vfo.seqrepoVersion
        return meta


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
