# -*- coding: utf-8 -*-

import logging

import vvhgvs.assemblymapper
import vvhgvs.exceptions

import VariantValidator.modules.gapped_mapping
import VariantValidator.modules.utils as fn
from VariantValidator.modules.transcript_map_data import TranscriptMapData
from VariantValidator.modules.variant import Variant


logger = logging.getLogger(
    f"VariantValidator.VariantFormatter.{__name__.removeprefix('VariantFormatter.')}"
)


def _get_alt_aln_method(
    transcript_accession,
    transcript_model="refseq",
    mixed_transcript_model=True,
):
    """
    Return the alignment method for a transcript.
    """
    if mixed_transcript_model:
        alt_aln_method = (
            "genebuild"
            if transcript_accession.startswith("ENST")
            else "splign"
        )
    else:
        alt_aln_method = (
            "genebuild"
            if transcript_model == "ensembl"
            else "splign"
        )

    logger.debug(
        "_get_alt_aln_method: transcript=%s, transcript_model=%s, "
        "mixed_transcript_model=%s, alt_aln_method=%s",
        transcript_accession,
        transcript_model,
        mixed_transcript_model,
        alt_aln_method,
    )

    return alt_aln_method


def compensate_g_to_t(
    hgvs_tx,
    hgvs_genomic,
    vm,
    hn,
    reverse_normalizer,
    primary_assembly,
    hdp,
    vfo,
    transcript_model="refseq",
    mixed_transcript_model=True,
):
    """
    Compensate genomic-to-transcript mapping across gapped alignments.
    """
    logger.info(
        "compensate_g_to_t ENTER: transcript=%s, genomic=%s, "
        "primary_assembly=%s, transcript_model=%s",
        hgvs_tx,
        hgvs_genomic,
        primary_assembly,
        transcript_model,
    )

    re_hash_hgvs_genomic = hgvs_genomic

    map_dat = TranscriptMapData(hdp=vfo.hdp)

    logger.info(
        "compensate_g_to_t %s: checking for gapped mapping against %s",
        hgvs_tx.ac,
        hgvs_genomic.ac,
    )

    gap_compensation = map_dat.is_gapped_map(
        hgvs_tx.ac,
        hgvs_genomic.ac,
    )

    logger.info(
        "compensate_g_to_t %s: gap_compensation=%s",
        hgvs_tx.ac,
        gap_compensation,
    )

    alt_aln_method = _get_alt_aln_method(
        hgvs_tx.ac,
        transcript_model=transcript_model,
        mixed_transcript_model=mixed_transcript_model,
    )

    logger.info(
        "compensate_g_to_t %s: alt_aln_method=%s",
        hgvs_tx.ac,
        alt_aln_method,
    )

    if not gap_compensation:
        logger.info(
            "compensate_g_to_t %s: no gapped mapping; "
            "calling fully_normalize",
            hgvs_tx.ac,
        )

        normalized_tx = fully_normalize(
            hgvs_tx,
            hgvs_genomic,
            hn,
            reverse_normalizer,
            vm,
            vfo,
            map_dat,
        )

        logger.info(
            "compensate_g_to_t %s: fully_normalize returned %s",
            hgvs_tx.ac,
            normalized_tx,
        )

        hgvs_tx_returns = [
            normalized_tx,
            False,
            None,
            None,
            None,
        ]

        gap_mapper = None

    else:
        logger.info(
            "compensate_g_to_t %s: gapped mapping detected; "
            "creating AssemblyMappers",
            hgvs_tx.ac,
        )

        no_norm_evm = vvhgvs.assemblymapper.AssemblyMapper(
            hdp,
            assembly_name=primary_assembly,
            alt_aln_method=alt_aln_method,
            normalize=False,
            replace_reference=True,
        )

        evm = vvhgvs.assemblymapper.AssemblyMapper(
            hdp,
            assembly_name=primary_assembly,
            alt_aln_method=alt_aln_method,
            normalize=True,
            replace_reference=True,
        )

        vfo.select_transcripts = "all"
        vfo.alt_aln_method = alt_aln_method
        vfo.reverse_hn = reverse_normalizer
        vfo.hn = hn

        # Variant currently accepts the submitted description as a string.
        # Keep this conversion at that API boundary only.
        logger.info(
            "compensate_g_to_t %s: creating Variant from genomic %s",
            hgvs_tx.ac,
            hgvs_genomic,
        )

        variant = Variant(fn.valstr(hgvs_genomic))

        variant.hgvs_genomic = hgvs_genomic
        variant.reverse_normalizer = reverse_normalizer
        variant.hn = hn
        variant.evm = evm
        variant.no_norm_evm = no_norm_evm
        variant.vm = vm
        variant.primary_assembly = primary_assembly
        variant.post_format_conversion = hgvs_genomic
        variant.map_dat = map_dat

        logger.info(
            "compensate_g_to_t %s: Variant configured: "
            "hgvs_genomic=%s, post_format_conversion=%s, "
            "primary_assembly=%s, alt_aln_method=%s",
            hgvs_tx.ac,
            variant.hgvs_genomic,
            variant.post_format_conversion,
            variant.primary_assembly,
            vfo.alt_aln_method,
        )

        gap_mapper = (
            VariantValidator.modules.gapped_mapping.GapMapper(
                variant,
                vfo,
            )
        )

        logger.info(
            "compensate_g_to_t %s: GapMapper created",
            hgvs_tx.ac,
        )

        # gapped_g_to_c currently expects transcript descriptions rather
        # than HGVS objects, so serialize only at this API boundary.
        logger.info(
            "compensate_g_to_t %s: calling gapped_g_to_c with %s",
            hgvs_tx.ac,
            hgvs_tx,
        )

        data, nw_rel_var = gap_mapper.gapped_g_to_c(
            [hgvs_tx],
            select_transcripts_dict={},
        )

        logger.info(
            "compensate_g_to_t %s: gapped_g_to_c returned "
            "nw_rel_var=%s, data=%s",
            hgvs_tx.ac,
            nw_rel_var,
            data,
        )

        logger.info(
            "compensate_g_to_t %s: obtaining exon orientation for %s/%s",
            hgvs_tx.ac,
            hgvs_tx.ac,
            hgvs_genomic.ac,
        )

        orientation = variant.map_dat.tx_exons(
            tx_ac=hgvs_tx.ac,
            alt_ac=hgvs_genomic.ac,
            alt_aln_method=vfo.alt_aln_method,
            hdp=vfo.hdp,
        )

        logger.info(
            "compensate_g_to_t %s: orientation=%s",
            hgvs_tx.ac,
            orientation,
        )

        try:
            logger.info(
                "compensate_g_to_t %s: calling g_to_t_compensation "
                "with transcript=%s",
                hgvs_tx.ac,
                nw_rel_var[0],
            )

            (
                re_hash_hgvs_genomic,
                _,
                hgvs_coding,
            ) = gap_mapper.g_to_t_compensation(
                orientation,
                nw_rel_var[0],
                "",
            )

            logger.info(
                "compensate_g_to_t %s: g_to_t_compensation returned "
                "genomic=%s, coding=%s",
                hgvs_tx.ac,
                re_hash_hgvs_genomic,
                hgvs_coding,
            )

        except vvhgvs.exceptions.HGVSInvalidVariantError as error:
            logger.info(
                "compensate_g_to_t %s: g_to_t_compensation raised "
                "HGVSInvalidVariantError: %s",
                hgvs_tx.ac,
                error,
            )

            if "insertion length must be 1" not in str(error):
                logger.info(
                    "compensate_g_to_t %s: error is not insertion-length "
                    "error; re-raising",
                    hgvs_tx.ac,
                )
                raise

            logger.info(
                "compensate_g_to_t %s: entering insertion-length "
                "validation fallback",
                hgvs_tx.ac,
            )

            validation = vfo.validate(
                fn.valstr(hgvs_genomic),
                primary_assembly,
                hgvs_tx.ac,
                liftover_level=None,
            ).format_as_dict()

            logger.info(
                "compensate_g_to_t %s: validation fallback returned %s",
                hgvs_tx.ac,
                validation,
            )

            hgvs_coding_description = next(
                (
                    key
                    for key in validation
                    if key.startswith(
                        ("NR_", "NM_", "ENST")
                    )
                ),
                None,
            )

            logger.info(
                "compensate_g_to_t %s: fallback coding description=%s",
                hgvs_tx.ac,
                hgvs_coding_description,
            )

            if hgvs_coding_description is not None:
                hgvs_coding = vfo.hp.parse_hgvs_variant(
                    hgvs_coding_description
                )
                nw_rel_var = [hgvs_coding]

                logger.info(
                    "compensate_g_to_t %s: fallback parsed coding=%s; "
                    "nw_rel_var=%s",
                    hgvs_tx.ac,
                    hgvs_coding,
                    nw_rel_var,
                )

        logger.info(
            "compensate_g_to_t %s: processing gap data",
            hgvs_tx.ac,
        )

        gapped_alignment_warning = data[
            "gapped_alignment_warning"
        ]

        is_alignment_artefact = (
            "does not represent a true variant"
            in gapped_alignment_warning
            or "may be an artefact of"
            in gapped_alignment_warning
        )

        logger.info(
            "compensate_g_to_t %s: "
            "gapped_alignment_warning=%s, is_alignment_artefact=%s",
            hgvs_tx.ac,
            gapped_alignment_warning,
            is_alignment_artefact,
        )

        logger.info(
            "compensate_g_to_t %s: constructing gap_compensated_tx "
            "from %s",
            hgvs_tx.ac,
            nw_rel_var[0],
        )

        gap_compensated_tx = [
            nw_rel_var[0],
            is_alignment_artefact,
            None,
            gapped_alignment_warning.replace(
                "the transcripts listed below",
                nw_rel_var[0].ac,
            ),
            data["auto_info"].replace("\n", ""),
        ]

        logger.info(
            "compensate_g_to_t %s: gap_compensated_tx=%s",
            hgvs_tx.ac,
            gap_compensated_tx,
        )

        if not gap_compensated_tx[1]:
            logger.info(
                "compensate_g_to_t %s: compensation is not an alignment "
                "artefact; calling fully_normalize",
                hgvs_tx.ac,
            )

            gap_compensated_tx[0] = fully_normalize(
                hgvs_tx,
                hgvs_genomic,
                hn,
                reverse_normalizer,
                vm,
                vfo,
                variant.map_dat,
            )

            logger.info(
                "compensate_g_to_t %s: fully_normalize returned %s",
                hgvs_tx.ac,
                gap_compensated_tx[0],
            )

        hgvs_tx_returns = gap_compensated_tx

    logger.info(
        "compensate_g_to_t %s: hgvs_tx_returns=%s",
        hgvs_tx.ac,
        hgvs_tx_returns,
    )

    hgvs_tx_dict = {
        "hgvs_transcript": hgvs_tx_returns[0],
        "gapped_alignment_warning": hgvs_tx_returns[3],
        "corrective_action": hgvs_tx_returns[2],
        "gap_position": hgvs_tx_returns[4],
        "transcript_accession": hgvs_tx_returns[0].ac,
        "hgvs_genomic": re_hash_hgvs_genomic,
    }

    logger.info(
        "compensate_g_to_t %s: initial output dictionary=%s",
        hgvs_tx.ac,
        hgvs_tx_dict,
    )

    if hgvs_tx_dict["gapped_alignment_warning"]:
        hgvs_tx_dict["gapped_alignment_warning"] = (
            "GappedAlignmentWarning: "
            f"{hgvs_tx_dict['gapped_alignment_warning']}"
        )

    if hgvs_tx_dict["corrective_action"]:
        hgvs_tx_dict["corrective_action"] = (
            "VariantMappingWarning"
            f"{hgvs_tx_dict['corrective_action']}"
        )

    if hgvs_tx_dict["gap_position"]:
        hgvs_tx_dict["gap_position"] = (
            "GappedAlignmentWarning: "
            f"{hgvs_tx_dict['gap_position']}"
        )

    # This branch is only relevant when the native gap mapper was created.
    if (
        gap_mapper is not None
        and not hgvs_tx_dict["gapped_alignment_warning"]
        and not hgvs_tx_dict["gap_position"]
    ):
        logger.info(
            "compensate_g_to_t %s: no gap warning returned; "
            "calling make_gap_warnings",
            hgvs_tx.ac,
        )

        gap_warnings = gap_mapper.make_gap_warnings(
            hgvs_tx.ac,
            hgvs_genomic.ac,
            primary_assembly,
        )

        logger.info(
            "compensate_g_to_t %s: make_gap_warnings returned %s",
            hgvs_tx.ac,
            gap_warnings,
        )

        hgvs_tx_dict["gapped_alignment_warning"] = (
            gap_warnings["gapped_alignment_warning"].replace(
                (
                    "Submitted description does not represent a true "
                    "variant because it is an artefact of aligning"
                ),
                (
                    "GappedAlignmentWarning: Variation described in "
                    "the context of an imperfect alignment of"
                ),
            )
        )

        hgvs_tx_dict["gap_position"] = (
            "GappedAlignmentWarning: "
            f"{gap_warnings['auto_info']}"
        )

    logger.info(
        "compensate_g_to_t EXIT %s: %s",
        hgvs_tx.ac,
        hgvs_tx_dict,
    )

    return hgvs_tx_dict


def fully_normalize(
    hgvs_tx,
    hgvs_genomic,
    hn,
    reverse_normalizer,
    vm,
    vfo,
    map_dat,
):
    """
    Fully normalize a transcript HGVS variant from its genomic variant.

    Used when gap-compensation output is not retained.
    """
    logger.info(
        "fully_normalize ENTER: transcript=%s, genomic=%s",
        hgvs_tx,
        hgvs_genomic,
    )

    tx_id = hgvs_tx.ac

    alt_aln_method = _get_alt_aln_method(
        tx_id,
        mixed_transcript_model=True,
    )

    logger.info(
        "fully_normalize %s: alt_aln_method=%s",
        tx_id,
        alt_aln_method,
    )

    exon_alignments = map_dat.tx_exons(
        tx_id,
        hgvs_genomic.ac,
        alt_aln_method,
    )

    logger.info(
        "fully_normalize %s: exon_alignments=%s",
        tx_id,
        exon_alignments,
    )

    orientation = int(
        exon_alignments[0]["alt_strand"]
    )

    logger.info(
        "fully_normalize %s: orientation=%s",
        tx_id,
        orientation,
    )

    if orientation == -1:
        logger.info(
            "fully_normalize %s: reverse-normalizing genomic %s",
            tx_id,
            hgvs_genomic,
        )

        hgvs_genomic = reverse_normalizer.normalize(
            hgvs_genomic
        )
    else:
        logger.info(
            "fully_normalize %s: forward-normalizing genomic %s",
            tx_id,
            hgvs_genomic,
        )

        hgvs_genomic = hn.normalize(
            hgvs_genomic
        )

    logger.info(
        "fully_normalize %s: normalized genomic=%s",
        tx_id,
        hgvs_genomic,
    )

    try:
        logger.info(
            "fully_normalize %s: calling vm.g_to_t",
            tx_id,
        )

        hgvs_tx = vm.g_to_t(
            hgvs_genomic,
            tx_id,
        )

        logger.info(
            "fully_normalize %s: vm.g_to_t returned %s",
            tx_id,
            hgvs_tx,
        )

    except vvhgvs.exceptions.HGVSError as error:
        logger.info(
            "fully_normalize %s: vm.g_to_t raised HGVSError: %s; "
            "retaining existing transcript=%s",
            tx_id,
            error,
            hgvs_tx,
        )

    try:
        logger.info(
            "fully_normalize %s: normalizing transcript %s",
            tx_id,
            hgvs_tx,
        )

        normalized_tx = hn.normalize(hgvs_tx)

        logger.info(
            "fully_normalize EXIT %s: normalized transcript=%s",
            tx_id,
            normalized_tx,
        )

        return normalized_tx

    except vvhgvs.exceptions.HGVSError as error:
        logger.info(
            "fully_normalize %s: transcript normalization failed "
            "for %s: %s",
            tx_id,
            hgvs_tx,
            error,
        )

        if "insertion length must be 1" not in str(error):
            logger.info(
                "fully_normalize EXIT %s: retaining unnormalized "
                "transcript=%s",
                tx_id,
                hgvs_tx,
            )
            return hgvs_tx

    logger.info(
        "fully_normalize %s: insertion-length recovery required",
        tx_id,
    )

    if orientation == -1:
        logger.info(
            "fully_normalize %s: recovery using forward normalizer",
            tx_id,
        )

        hgvs_genomic = hn.normalize(
            hgvs_genomic
        )
    else:
        logger.info(
            "fully_normalize %s: recovery using reverse normalizer",
            tx_id,
        )

        hgvs_genomic = reverse_normalizer.normalize(
            hgvs_genomic
        )

    logger.info(
        "fully_normalize %s: recovery genomic=%s",
        tx_id,
        hgvs_genomic,
    )

    hgvs_tx = vfo.vm.g_to_t(
        hgvs_genomic,
        tx_id,
        alt_aln_method=alt_aln_method,
    )

    logger.info(
        "fully_normalize %s: recovery g_to_t returned %s",
        tx_id,
        hgvs_tx,
    )

    try:
        normalized_tx = hn.normalize(hgvs_tx)

        logger.info(
            "fully_normalize EXIT %s: recovery normalized "
            "transcript=%s",
            tx_id,
            normalized_tx,
        )

        return normalized_tx

    except vvhgvs.exceptions.HGVSInvalidVariantError as error:
        logger.info(
            "fully_normalize EXIT %s: recovery normalization failed: "
            "%s; retaining transcript=%s",
            tx_id,
            error,
            hgvs_tx,
        )

        return hgvs_tx


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
