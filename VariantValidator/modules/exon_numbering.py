"""
exon_numbering.py Module

Authors: Katie Williams (@kwi11iams) and Katherine Winfield (@kjwinfield)

This code will ultimately aim to provide exon numbering information for
VariantValidator

See exon_numbering.md markdown for a full description on how this
module operates.

Use exon_numbering_tests.py for automated testing of this module.
"""
import vvhgvs
import vvhgvs.exceptions
import re


def finds_exon_number(variant, validator):
    """
    :param variant: (obj): the variant object from VariantValidator
    :param validator (obj)
    :return: exon_start_end_positions (dict): a dictionary of the
                    exon/intron positions for the start and end of the given
                    variant for each aligned chromosomal or gene reference
                    sequence

    Function that finds and output exon numbering for a given variant
    """
    response_dictionary = validator.gene2transcripts(variant, validator, bypass_web_searches=True)

    # Filter out the response_dictionary for the variant transcript
    # This will find the exon structure dictionary for the given transcript
    # And select the coding start position number

    # Step 1 - Pull out all the attributes for the transcript in the context of all genome builds
    info_dict = {}
    for i in range(len(response_dictionary["transcripts"])):

        if response_dictionary["transcripts"][i]["reference"] == variant.hgvs_coding.ac:

            # Create record
            info_dict[(response_dictionary["transcripts"][i]["reference"])] = {}

            # Returns an exon structure dictionary
            info_dict[(response_dictionary["transcripts"][i]["reference"])]["exon_structure_dict"] = \
                response_dictionary["transcripts"][i]["genomic_spans"]

            # Returns the start of coding
            # (This is needed to correct the position)
            info_dict[(response_dictionary["transcripts"][i]["reference"])]['coding_start'] = \
                response_dictionary['transcripts'][i]["coding_start"]

    # Step 2 - Get the necessary variant information
    # Find the variant position from the variant nomenclature
    coordinates = str(variant.hgvs_coding.posedit.pos)

    # Identify start and end of variant from input coordinates
    if '_' in coordinates:
        start_position, end_position = coordinates.split('_')
    else:
        # If SNV, then start = end position
        start_position = coordinates
        end_position = coordinates

    # Create empty output dictionary
    exon_start_end_positions = {}

    # Create c_to_n varint
    try:
        to_n = validator.vm.c_to_n(variant.hgvs_coding)
    except vvhgvs.exceptions.HGVSInvalidVariantError:
        to_n = variant.hgvs_coding

    """
    This for loop identifies the exon/intron number for the transcript
    variant for each aligned chromosomal or gene reference sequence
    It populates output dictionary with the aligned chromosomal and gene
    reference sequences as keys
    Each of these keys has another dictionary as its value:
        keys: start_exon and end_exon
        values: start and position of variant in the reference sequence
    """
    exon_structure_dict = info_dict[variant.hgvs_coding.ac]["exon_structure_dict"]
    coding_start = info_dict[variant.hgvs_coding.ac]['coding_start']
    if coding_start is None:
        coding_start = 1

    for transcript in exon_structure_dict:
        for exon in exon_structure_dict[transcript]['exon_structure']:
            # For loop that runs to identify which exon/inton the variant is in
            # 'i' denotes introns, i.e. exon 2i is intron 2
            # Separated by start and end position of the variant as they may
            # be different if the variant is not a SNP.
            if ('+' not in str(start_position)
                    and not re.search('\d-\d', str(start_position))):
                # This works for positions in exons
                adj_start_position = to_n.posedit.pos.start.base
                if exon['transcript_start'] <= adj_start_position <= exon['transcript_end']:
                    start_exon = str(exon['exon_number'])

            elif re.match("-", str(start_position)):
                n_start_position = str(to_n.posedit.pos.start)
                if re.search("\d\+\d", str(n_start_position)):
                    nearest_exon_boundary = int(str(n_start_position).split('+')[0])
                    if nearest_exon_boundary == exon['transcript_end']:
                        start_exon = str(exon['exon_number']) + 'i'
                elif re.search("\d-\d", str(n_start_position)):
                    nearest_exon_boundary = int(str(n_start_position).split('-')[0])
                    if nearest_exon_boundary == exon['transcript_start']:
                        start_exon = str(exon['exon_number'] - 1) + 'i'

            elif re.match("\*", start_position) and "+" in start_position:
                n_start_position = str(to_n.posedit.pos.start)
                nearest_exon_boundary = int(str(n_start_position).split('+')[0])
                if nearest_exon_boundary == exon['transcript_end']:
                    start_exon = str(exon['exon_number']) + 'i'

            elif re.match("\*", start_position) and "-" in start_position:
                n_start_position = str(to_n.posedit.pos.start)
                nearest_exon_boundary = int(str(n_start_position).split('-')[0])
                if nearest_exon_boundary == exon['transcript_start']:
                    start_exon = str(exon['exon_number'] - 1) + 'i'

            elif '+' in start_position:
                # This works for positions that are + the exon boundary
                nearest_exon_boundary = start_position.split('+')[0]
                adj_nearest_exon_boundary = (int(nearest_exon_boundary)
                                             + coding_start - 1)
                if adj_nearest_exon_boundary == exon['transcript_end']:
                    start_exon = str(exon['exon_number']) + 'i'

            elif '-' in start_position:
                # This works for positions that are - the exon boundary
                nearest_exon_boundary = start_position.split('-')[0]
                adj_nearest_exon_boundary = (int(nearest_exon_boundary)
                                             + coding_start - 1)
                if adj_nearest_exon_boundary == exon['transcript_start']:
                    start_exon = str(exon['exon_number'] - 1) + 'i'

            # End position
            if '+' not in str(end_position) and not re.search('\d-\d', str(end_position)):
                # This works for positions in exons
                adj_end_position = to_n.posedit.pos.end.base
                if exon['transcript_start'] <= adj_end_position <= exon['transcript_end']:
                    end_exon = str(exon['exon_number'])

            elif re.match("-", end_position):
                n_end_position = str(to_n.posedit.pos.end)
                if re.search("\d\+\d", str(n_end_position)):
                    nearest_exon_boundary = int(str(n_end_position).split('+')[0])
                    if nearest_exon_boundary == exon['transcript_end']:
                        end_exon = str(exon['exon_number']) + 'i'
                elif re.search("\d-\d", str(n_end_position)):
                    nearest_exon_boundary = int(str(n_end_position).split('-')[0])
                    if nearest_exon_boundary == exon['transcript_start']:
                        end_exon = str(exon['exon_number'] -1) + 'i'

            elif re.match("\*", end_position) and "+" in end_position:
                n_end_position = str(to_n.posedit.pos.end)
                nearest_exon_boundary = int(str(n_end_position).split('+')[0])
                if nearest_exon_boundary == exon['transcript_end']:
                    end_exon = str(exon['exon_number']) + 'i'

            elif re.match("\*", end_position) and "-" in end_position:
                n_end_position = str(to_n.posedit.pos.end)
                nearest_exon_boundary = int(str(n_end_position).split('-')[0])
                if nearest_exon_boundary == exon['transcript_start']:
                    end_exon = str(exon['exon_number'] - 1) + 'i'

            elif '+' in end_position:
                # This works for positions that are + the exon boundary
                nearest_exon_boundary = end_position.split('+')[0]
                adj_nearest_exon_boundary = (int(nearest_exon_boundary)
                                             + coding_start - 1)
                if adj_nearest_exon_boundary == exon['transcript_end']:
                    end_exon = str(exon['exon_number']) + 'i'

            elif '-' in end_position:
                # This works for positions that are - the exon boundary
                if "+" in end_position:
                    adj_end_position = end_position.split("+")[0]
                else:
                    adj_end_position = end_position
                nearest_exon_boundary = adj_end_position.split('-')[0]
                adj_nearest_exon_boundary = (int(nearest_exon_boundary)
                                             + coding_start - 1)
                if adj_nearest_exon_boundary == exon['transcript_start']:
                    end_exon = str(exon['exon_number'] - 1) + 'i'
                elif adj_nearest_exon_boundary == exon['transcript_end']:
                    end_exon = str(exon['exon_number']) + 'i'

        try:
            exon_start_end_positions[transcript] = {"start_exon": start_exon,
                                                    "end_exon": end_exon}

        # This happens in genes like Shank 2 where there are problems with aligning exons to the genome build,
        # i.e. if there are exons in the transcript that are not aligned to the genome, it is impossible to calculate
        # the exon number
        except UnboundLocalError:
            exon_start_end_positions[transcript] = {"start_exon": "cannot be calculated",
                                                    "end_exon": "cannot be calculated"}

    return exon_start_end_positions

# <LICENSE>
# Copyright (C) 2016-2024 VariantValidator Contributors
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

