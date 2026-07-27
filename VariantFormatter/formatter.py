# -*- coding: utf-8 -*-

"""
The core variantformatter functions:

vcf2hgvs_genomic
hgvs_genomic2vcf
hgvs_genomic2hgvs_transcript
hgvs_transcript2hgvs_protein
fetch_aligned_transcripts
remove_reference
gap_checker

A brief description is given above each function
"""

# Import Python modules
import copy

# Import hgvs modules
import vvhgvs.assemblymapper
import vvhgvs.exceptions

# Import VariantFormatter modules
import VariantValidator.modules.hgvs_utils as hgvs_utils
import VariantFormatter.gapGenes as gapGenes
from VariantValidator.modules import format_converters, utils
from VariantValidator.modules.variant import Variant


import logging

logger = logging.getLogger(
    f"VariantValidator.VariantFormatter.{__name__.removeprefix('VariantFormatter.')}"
)

"""
Simple function, parses an HGVS string into a hgvs (.py) object
NOTE: if the string is not valid, hgvs (py) parser will throw up an error
This function may be best used in a try/except statement
See hgvs (py) documentation for error types
"""


def parse(hgvs_string, vfo):
    hgvs_object = vfo.hp.parse_hgvs_variant(hgvs_string)
    return hgvs_object


"""
Function which takes a vcf string in the format Chr:Pos:Ref:Alt or Chr-Pos-Ref-Alt and the required genome build
and generates a hgvs python object. The object is validated to ensure the stated ref appears at the stated
positions in the RefSeq RefSeq reference sequence

Supported genome builds are GRCh38, GRCh37, hg19, hg38
"""


def vcf2hgvs_genomic(pseudo_vcf, genome_build, vfo):
    """
    Convert VCF-style genomic input to a validated and normalized genomic
    HGVS SequenceVariant.

    VariantValidator's VCF conversion stages are used in their normal order:

        vcf2hgvs_stage1
        vcf2hgvs_stage2
        vcf2hgvs_stage4

    HGVS variants are retained as objects once parsed or constructed.
    """

    result = {
        "error": "",
        "hgvs_genomic": "",
        "ref_bases": "",
        "un_normalized_hgvs_genomic": "",
    }

    variant = Variant(
        pseudo_vcf,
        quibble=pseudo_vcf,
        primary_assembly=genome_build,
    )

    # Required by the VV VCF conversion path.
    variant.hn = vfo.splign_normalizer

    batch_list = []

    # VCF type 1:
    # chr-pos-ref-alt -> chr:posRef>Alt
    toskip = format_converters.vcf2hgvs_stage1(
        variant,
        batch_list,
    )

    if toskip:
        if variant.warnings:
            result["error"] = variant.warnings[-1]
        return result

    # API type non-HGVS:
    # e.g. Chr16:2099572TC>T
    toskip = format_converters.vcf2hgvs_stage2(
        variant,
        vfo,
    )

    if toskip:
        if variant.warnings:
            result["error"] = variant.warnings[-1]
        return result

    # Convert non-substitution types such as GGGG>G to the appropriate
    # HGVS edit representation.
    toskip = format_converters.vcf2hgvs_stage4(
        variant,
        batch_list,
    )

    if toskip:
        if variant.warnings:
            result["error"] = variant.warnings[-1]
        return result

    # The conversion stages may already have produced an HGVS object.
    if isinstance(variant.quibble, str):
        try:
            hgvs_genomic = parse(
                variant.quibble,
                vfo,
            )
        except vvhgvs.exceptions.HGVSError as error:
            result["error"] = str(error)
            return result
    else:
        hgvs_genomic = variant.quibble

    # Validate the resulting genomic HGVS object.
    try:
        vfo.vr.validate(hgvs_genomic)

    except vvhgvs.exceptions.HGVSError as error:
        result["error"] = str(error)

    else:
        # Retain the variant before final normalization.
        result["un_normalized_hgvs_genomic"] = copy.deepcopy(
            hgvs_genomic
        )

        # Normalize once.
        hgvs_genomic = vfo.splign_normalizer.normalize(
            hgvs_genomic
        )

        result["hgvs_genomic"] = hgvs_genomic

        result["ref_bases"] = vfo.sf.fetch_seq(
            hgvs_genomic.ac,
            start_i=hgvs_genomic.posedit.pos.start.base - 1,
            end_i=hgvs_genomic.posedit.pos.end.base,
        )

    if result["error"]:
        errors = utils.normalise_warning_codes(
            [result["error"]]
        )
        result["error"] = errors[0]

    return result


"""
Function takes a genomic HGVS description and returns the component parts of a VCF
"""


def hgvs_genomic2vcf(hgvs_genomic, genome_build, vfo):
    vcf_dictionary = hgvs_utils.report_hgvs2vcf(hgvs_genomic, genome_build, vfo.reverse_splign_normalizer, vfo.sf)
    return vcf_dictionary


"""
Function which takes a hgvs Python genomic variant and maps to a specified transcript reference sequence. A transcript
level hgvs python object is returned. 
"""


def hgvs_genomic2hgvs_transcript(hgvs_genomic, tx_id, vfo):
    """
    Map a genomic HGVS SequenceVariant to a specified transcript.

    HGVS variants are retained as objects throughout mapping and
    normalization.
    """

    result = {
        "error": "",
        "hgvs_transcript": "",
        "ref_bases": "",
        "latest_version": tx_id,
    }

    # Configure mapping and normalization for the transcript source.
    if tx_id.startswith("ENST"):
        alt_aln_method = "genebuild"
        hn = vfo.genebuild_normalizer
        rhn = vfo.reverse_genebuild_normalizer

    elif tx_id.startswith(("NM_", "NR_")):
        alt_aln_method = "splign"
        hn = vfo.splign_normalizer
        rhn = vfo.reverse_splign_normalizer

    else:
        result["error"] = (
            f"TranscriptDataError: Unsupported transcript reference "
            f"sequence {tx_id}"
        )
        return result

    # Check for a more recent transcript version on this genomic alignment.
    tx_id_info = vfo.hdp.get_tx_identity_info(tx_id)
    uta_gene_symbol = tx_id_info[6]
    tx_for_gene = vfo.hdp.get_tx_for_gene(uta_gene_symbol)

    ac_root, ac_version = tx_id.rsplit(".", 1)
    current_version = int(ac_version)
    latest_version = current_version
    latest_accession = tx_id

    for accession in tx_for_gene:
        if hgvs_genomic.ac != accession[4]:
            continue

        candidate_ac = accession[3]

        if not candidate_ac.startswith(f"{ac_root}."):
            continue

        try:
            candidate_version = int(
                candidate_ac.rsplit(".", 1)[1]
            )
        except (ValueError, IndexError):
            continue

        if candidate_version > latest_version:
            latest_version = candidate_version
            latest_accession = candidate_ac

    if latest_accession != tx_id:
        result["latest_version"] = (
            f"TranscriptVersionWarning: A more recent version of the "
            f"selected reference sequence {tx_id} is available for "
            f"genome build {vfo.genome_build} ({latest_accession})"
        )
    else:
        result["latest_version"] = None

    # Obtain transcript orientation relative to the genomic accession.
    try:
        exon_alignments = vfo.hdp.get_tx_exons(
            tx_id,
            hgvs_genomic.ac,
            alt_aln_method,
        )
        orientation = int(
            exon_alignments[0]["alt_strand"]
        )

    except Exception:
        result["error"] = (
            f"TranscriptDataError: No alignment data available for "
            f"transcript {tx_id} and chromosome {hgvs_genomic.ac}"
        )
        return result

    # Normalize genomic variant in the appropriate direction.
    if orientation == -1:
        mapped_genomic = rhn.normalize(hgvs_genomic)
    else:
        mapped_genomic = hn.normalize(hgvs_genomic)

    logger.info(
        "VF g_to_t %s: mapped genomic = %s, edit type = %s",
        tx_id,
        mapped_genomic,
        mapped_genomic.posedit.edit.type,
    )

    # Map genomic variant directly to transcript coordinates.
    try:
        hgvs_tx = vfo.vm.g_to_t(
            mapped_genomic,
            tx_id,
            alt_aln_method=alt_aln_method,
        )

        logger.info(
            "VF g_to_t %s: raw transcript = %s, edit type = %s",
            tx_id,
            hgvs_tx,
            hgvs_tx.posedit.edit.type,
        )

    except vvhgvs.exceptions.HGVSError as error:
        logger.info(
            "VF g_to_t %s failed for %s: %s",
            tx_id,
            mapped_genomic,
            error,
        )
        result["error"] = str(error)

    else:
        # Ensure complete transcript normalization.
        try:
            hgvs_tx = hn.normalize(hgvs_tx)

            logger.info(
                "VF g_to_t %s: normalized transcript = %s, edit type = %s",
                tx_id,
                hgvs_tx,
                hgvs_tx.posedit.edit.type,
            )

        except vvhgvs.exceptions.HGVSError as error:
            logger.info(
                "VF g_to_t %s: transcript normalization failed for %s: %s",
                tx_id,
                hgvs_tx,
                error,
            )

            # Preserve the existing recovery for insertion normalization
            # on reverse-strand alignments.
            if (
                "insertion length must be 1" in str(error)
                and orientation == -1
            ):
                mapped_genomic = hn.normalize(
                    mapped_genomic
                )

                logger.info(
                    "VF g_to_t %s: retry mapped genomic = %s, "
                    "edit type = %s",
                    tx_id,
                    mapped_genomic,
                    mapped_genomic.posedit.edit.type,
                )

                hgvs_tx = vfo.vm.g_to_t(
                    mapped_genomic,
                    tx_id,
                    alt_aln_method=alt_aln_method,
                )

                logger.info(
                    "VF g_to_t %s: retry raw transcript = %s, "
                    "edit type = %s",
                    tx_id,
                    hgvs_tx,
                    hgvs_tx.posedit.edit.type,
                )

                try:
                    hgvs_tx = hn.normalize(hgvs_tx)

                    logger.info(
                        "VF g_to_t %s: retry normalized transcript = %s, "
                        "edit type = %s",
                        tx_id,
                        hgvs_tx,
                        hgvs_tx.posedit.edit.type,
                    )

                except vvhgvs.exceptions.HGVSInvalidVariantError as retry_error:
                    logger.info(
                        "VF g_to_t %s: retry transcript normalization "
                        "failed for %s: %s",
                        tx_id,
                        hgvs_tx,
                        retry_error,
                    )
                    pass

            # Other normalization failures are intentionally retained as
            # the mapped object. These may represent intronic variants.

        result["hgvs_transcript"] = hgvs_tx

        logger.info(
            "VF g_to_t %s: transcript stored in result = %s",
            tx_id,
            result["hgvs_transcript"],
        )

        try:
            result["ref_bases"] = hgvs_tx.posedit.edit.ref

        except vvhgvs.exceptions.HGVSError:
            sequence_hgvs = hgvs_tx

            if sequence_hgvs.type == "c":
                sequence_hgvs = vfo.vm.c_to_n(
                    sequence_hgvs
                )

            if not hgvs_utils.either_position_is_intronic(
                sequence_hgvs
            ):
                result["ref_bases"] = vfo.sf.fetch_seq(
                    sequence_hgvs.ac,
                    start_i=sequence_hgvs.posedit.pos.start.base - 1,
                    end_i=sequence_hgvs.posedit.pos.end.base,
                )
            else:
                result["ref_bases"] = ""

    # Normalize non-transcript-specific errors.
    if (
        result["error"]
        and "Transcript" not in result["error"]
    ):
        errors = utils.normalise_warning_codes(
            [result["error"]]
        )
        result["error"] = errors[0]

    logger.info(
        "VF g_to_t %s: final result hgvs_transcript = %s, error = %s",
        tx_id,
        result["hgvs_transcript"],
        result["error"],
    )

    return result


"""
Function which takes a hgvs Python transctript variant and maps to a specified protein reference sequence. A protein
level hgvs python object is returned.

Note the function currently assumes that the transcript description is correctly normalized having come from the 
previous g_to_t function
"""


def hgvs_transcript2hgvs_protein(hgvs_transcript, genome_build, vfo):
    """
    Map a transcript HGVS SequenceVariant to a protein HGVS
    SequenceVariant.

    Ensembl transcripts use the genebuild alignment method.
    RefSeq transcripts use the splign alignment method.

    Returns
    -------
    vvhgvs.sequencevariant.SequenceVariant
        Protein HGVS object.
    """

    if hgvs_transcript.ac.startswith("ENST"):
        alt_aln_method = "genebuild"
        rhn = vfo.reverse_genebuild_normalizer
    else:
        alt_aln_method = "splign"
        rhn = vfo.reverse_splign_normalizer

    # TODO: Reuse persistent AssemblyMapper instances from Validator
    # rather than creating a new mapper for each protein conversion.
    evm = vvhgvs.assemblymapper.AssemblyMapper(
        vfo.hdp,
        assembly_name=genome_build,
        alt_aln_method=alt_aln_method,
        normalize=True,
        replace_reference=True,
    )

    result = vfo.myc_to_p(
        hgvs_transcript,
        evm,
        False,
        rhn,
    )

    return result["hgvs_protein"]


"""
Return all aligned transcripts for a given genomic hgvs object
"""


def fetch_aligned_transcripts(
        hgvs_genomic,
        transcript_model,
        vfo,
        genome_build,
):
    """
    Return transcripts aligned to a genomic HGVS SequenceVariant.

    Ensembl transcripts use genebuild alignments and RefSeq transcripts
    use splign alignments. Both direct regional alignments and transcripts
    identified by AssemblyMapper are retained.

    The transcript-edge fallback is preserved for antisense mappings.
    """

    tx_list = []

    genome_build = (
        "GRCh38"
        if "38" in genome_build
        else "GRCh37"
    )

    start = hgvs_genomic.posedit.pos.start.base
    end = hgvs_genomic.posedit.pos.end.base

    def get_region_transcripts(alt_aln_method):
        transcripts = vfo.hdp.get_tx_for_region(
            hgvs_genomic.ac,
            alt_aln_method,
            start - 1,
            end,
        )

        # Preserve transcript-edge antisense fallback.
        if not transcripts:
            transcripts = vfo.hdp.get_tx_for_region(
                hgvs_genomic.ac,
                alt_aln_method,
                start,
                end - 1,
            )

        return transcripts

    def get_relevant_transcripts(alt_aln_method):
        # TODO: Reuse persistent AssemblyMapper instances from Validator
        # rather than creating a new mapper for each transcript search.
        evm = vvhgvs.assemblymapper.AssemblyMapper(
            vfo.hdp,
            assembly_name=genome_build,
            alt_aln_method=alt_aln_method,
            normalize=True,
            replace_reference=True,
        )

        return [
            [tx]
            for tx in evm.relevant_transcripts(hgvs_genomic)
        ]

    def has_version(accession):
        """
        Return True for a versioned transcript accession.
        """
        if "." not in accession:
            return False

        root, version = accession.rsplit(".", 1)

        return bool(root) and version.isdigit()

    if transcript_model in ("ensembl", "all"):
        enst_list = get_region_transcripts(
            "genebuild"
        )

        relevant_enst = get_relevant_transcripts(
            "genebuild"
        )

        tx_list.extend(
            transcript
            for transcript in enst_list
            if has_version(transcript[0])
        )

        tx_list.extend(
            transcript
            for transcript in relevant_enst
            if has_version(transcript[0])
        )

    if transcript_model in ("refseq", "all"):
        refseq_list = get_region_transcripts(
            "splign"
        )

        relevant_refseq = get_relevant_transcripts(
            "splign"
        )

        tx_list.extend(refseq_list)
        tx_list.extend(relevant_refseq)

    # Preserve explicitly requested/raw transcript sets. Otherwise filter
    # the collected transcripts to the current versions.
    select_transcripts = vfo.select_transcripts

    if (
        select_transcripts not in ("raw", "select")
        and not (
            select_transcripts
            and any(
                accession_type in select_transcripts
                for accession_type in (
                    "NR",
                    "NM_",
                    "ENST",
                )
            )
        )
    ):
        tx_list = vfo.transcript_filter(tx_list)

    return tx_list


"""
Format a nucleotide HGVS description without displaying reference bases.
"""


def remove_reference(hgvs_nucleotide):
    return hgvs_nucleotide.format(
        {"max_ref_length": 0}
    )


"""
Check transcript HGVS descriptions for alignment gaps.
"""


def gap_checker(
        hgvs_transcript,
        hgvs_genomic,
        genome_build,
        vfo,
        transcript_model="refseq",
):
    tx_id = hgvs_transcript.ac

    logger.info(
        "VF gap_checker ENTER %s: transcript=%s, edit_type=%s, genomic=%s",
        tx_id,
        hgvs_transcript,
        hgvs_transcript.posedit.edit.type,
        hgvs_genomic,
    )

    if tx_id.startswith("ENST"):
        hn = vfo.genebuild_normalizer
        rhn = vfo.reverse_genebuild_normalizer

    elif tx_id.startswith(("NM_", "NR_")):
        hn = vfo.splign_normalizer
        rhn = vfo.reverse_splign_normalizer

    else:
        raise ValueError(
            f"Unsupported transcript accession: {tx_id}"
        )

    checked = gapGenes.compensate_g_to_t(
        hgvs_transcript,
        hgvs_genomic,
        vfo.vm,
        hn,
        rhn,
        genome_build,
        vfo.hdp,
        vfo,
        transcript_model=transcript_model,
    )

    logger.info(
        "VF gap_checker EXIT %s: transcript=%s, edit_type=%s, "
        "warning=%s, corrective_action=%s, gap_position=%s",
        tx_id,
        checked.get("hgvs_transcript"),
        (
            checked["hgvs_transcript"].posedit.edit.type
            if checked.get("hgvs_transcript") is not None
            else None
        ),
        checked.get("gapped_alignment_warning"),
        checked.get("corrective_action"),
        checked.get("gap_position"),
    )

    return checked


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
