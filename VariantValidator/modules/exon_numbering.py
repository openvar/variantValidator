# Import vv_hgvs modules
import vvhgvs
import vvhgvs.exceptions

def _find_exon(position, exon_structure):
    """
    Find the exon or intron containing an HGVS transcript position.

    Exonic positions return the exon number.
    Intronic positions return the preceding exon number followed by ``i``.

    :param position: vvhgvs BaseOffsetPosition
    :param exon_structure: exon structure returned by gene2transcripts
    :return: exon/intron number as a string, or None if it cannot be found
    """

    # Exonic position
    if position.offset == 0:
        for exon in exon_structure:
            if (
                    exon["transcript_start"]
                    <= position.base
                    <= exon["transcript_end"]
            ):
                return str(exon["exon_number"])

        return None

    # Intronic position described relative to the end of the preceding exon,
    # e.g. n.100+1.
    if position.offset > 0:
        for exon in exon_structure:
            if position.base == exon["transcript_end"]:
                return f'{exon["exon_number"]}i'

        return None

    # Intronic position described relative to the start of the following exon,
    # e.g. n.101-1.
    if position.offset < 0:
        for exon in exon_structure:
            if position.base == exon["transcript_start"]:
                return f'{exon["exon_number"] - 1}i'

        return None
    return None

def finds_exon_number(variant, validator):
    """
    Find exon/intron numbering for the start and end positions of a variant.

    :param variant: VariantValidator variant object
    :param validator: VariantValidator validator object
    :return: dictionary containing start/end exon or intron numbers for each
             aligned chromosomal or gene reference sequence
    """
    response_dictionary = validator.gene2transcripts(
        variant,
        validator,
        bypass_web_searches=True
    )
    # Find the transcript record corresponding to the submitted transcript.
    transcript_info = None

    for transcript in response_dictionary["transcripts"]:
        if transcript["reference"] == variant.hgvs_coding.ac:
            transcript_info = transcript
            break

    if transcript_info is None:
        return {}

    exon_structure_dict = transcript_info["genomic_spans"]

    # Work in n. coordinates so that exon lookup uses transcript positions
    # directly and does not need to account manually for CDS start/end.
    try:
        hgvs_transcript = validator.vm.c_to_n(variant.hgvs_coding)
    except vvhgvs.exceptions.HGVSInvalidVariantError:
        hgvs_transcript = variant.hgvs_coding

    start_position = hgvs_transcript.posedit.pos.start
    end_position = hgvs_transcript.posedit.pos.end
    exon_start_and_end_positions = {}

    for accession, transcript_data in exon_structure_dict.items():
        exon_structure = transcript_data["exon_structure"]

        start_exon = _find_exon(
            start_position,
            exon_structure
        )

        end_exon = _find_exon(
            end_position,
            exon_structure
        )

        if start_exon is None:
            start_exon = "cannot be calculated"

        if end_exon is None:
            end_exon = "cannot be calculated"

        exon_start_and_end_positions[accession] = {
            "start_exon": start_exon,
            "end_exon": end_exon,
        }

    return exon_start_and_end_positions

# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
