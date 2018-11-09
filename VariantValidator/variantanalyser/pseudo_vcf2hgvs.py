"""
psuedo_vcf2hgvs is a stripped down version of VariantValidator's vcf2hgvs functionality
The tool is used to convert only the pseudo VCF format e.g. chr-pos-ref-alt into hgvs
python objects. The two variations on the function provide the HGVS output as either a
3 prime normalized HGVS description or a 5 prime normalized HGVS description. 5 prime
normalization is primarily used in the process of merging several VCF calls into a single
HGVS description
"""
# Import  modules
import re
import copy
import hgvs
import hgvs.dataproviders
import hgvs.normalizer
import hgvs.parser
import supported_chromosome_builds as va_scb
from dbControls import data as va_dbCrl


# Error handling
class pseudoVCF2HGVSError(Exception):
    pass


# Set variables
hdp = hgvs.dataproviders.uta.connect(pooling=True)

# Reverse normalizer (5 prime)
reverse_normalize = hgvs.normalizer.Normalizer(hdp,
                                               cross_boundaries=False,
                                               shuffle_direction=5,
                                               alt_aln_method='splign'
                                               )

# normalizer (3 prime)
normalize = hgvs.normalizer.Normalizer(hdp,
                                       cross_boundaries=False,
                                       shuffle_direction=3,
                                       alt_aln_method='splign'
                                       )

# parser
hp = hgvs.parser.Parser()
# SeqFetcher
sf = hgvs.dataproviders.seqfetcher.SeqFetcher()


# pvcf is a pseudo_vcf string
# genome build is a build string e.g. GRCh37 hg19
# normalization direction an integer, 5 or 3.
def pvcf_to_hgvs(input, selected_assembly, normalization_direction):
    # Set normalizer
    if normalization_direction == 3:
        selected_normalizer = normalize
    if normalization_direction == 5:
        selected_normalizer = reverse_normalize

    # Gel stye pVCF
    input = input.replace(':', '-')
    
    # VCF type 1
    if re.search('-\d+-[GATC]+-[GATC]+', input):
        pre_input = copy.deepcopy(input)
        vcf_elements = pre_input.split('-')
        input = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], vcf_elements[3])
    elif re.search('-\d+-[GATC]+-', input):
        pre_input = copy.deepcopy(input)
        vcf_elements = pre_input.split('-')
        input = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], vcf_elements[2])
    else:
        raise pseudoVCF2HGVSError('Unsupported format: VCF specification 4.1 or later')

    # Chr16:2099572TC>T
    try:
        pre_input = copy.deepcopy(input)
        input_list = input.split(':')
        pos_ref_alt = str(input_list[1])
        positionAndEdit = input_list[1]
        if not re.match('N[CGWT]_', input) and not re.match('LRG_\d+$', input):
            chr_num = str(input_list[0])
            chr_num = chr_num.upper()
            chr_num = chr_num.strip()
            if re.match('CHR', chr_num):
                chr_num = chr_num.replace('CHR', '')
            # Use selected assembly
            accession = va_scb.to_accession(chr_num, selected_assembly)
            if accession is None:
                error = chr_num + ' is not part of genome build ' + selected_assembly + ' or is not supported'
                raise pseudoVCF2HGVSError(error)
        else:
            accession = input_list[0]

        # Assign reference sequence type
        ref_type = ':g.'
        if re.match('LRG_', accession):
            accession = va_dbCrl.get_RefSeqGeneID_from_lrgID(accession)

        # Reformat the variant
        input = str(accession) + ref_type + str(positionAndEdit)
    except Exception as e:
        error = str(e)
        raise pseudoVCF2HGVSError(error)

    # Find not_sub type in input e.g. GGGG>G
    not_sub = copy.deepcopy(input)
    not_sub_find = re.compile("([GATCgatc]+)>([GATCgatc]+)")
    if not_sub_find.search(not_sub):
        try:
            # If the length of either side of the substitution delimer (>) is >1
            matches = not_sub_find.search(not_sub)
            if len(matches.group(1)) > 1 or len(matches.group(2)) > 1 or re.search(
                    "([GATCgatc]+)>([GATCgatc]+),([GATCgatc]+)", input):
                # Search for and remove range
                range = re.compile("([0-9]+)_([0-9]+)")
                if range.search(not_sub):
                    m = not_sub_find.search(not_sub)
                    start = m.group(1)
                    delete = m.group(2)
                    beginning_string, middle_string = not_sub.split(':')
                    middle_string = middle_string.split('_')[0]
                    end_string = start + '>' + delete
                    not_sub = beginning_string + ':' + middle_string + end_string
                # Split description
                split_colon = not_sub.split(':')
                ref_ac = split_colon[0]
                remainder = split_colon[1]
                split_dot = remainder.split('.')
                ref_type = split_dot[0]
                remainder = split_dot[1]
                posedit = remainder
                split_greater = remainder.split('>')
                insert = split_greater[1]
                remainder = split_greater[0]
                # Split remainder using matches
                r = re.compile("([0-9]+)([GATCgatc]+)")
                try:
                    m = r.search(remainder)
                    start = m.group(1)
                    delete = m.group(2)
                    starts = posedit.split(delete)[0]
                    re_try = ref_ac + ':' + ref_type + '.' + starts + 'del' + delete[0] + 'ins' + insert
                    hgvs_re_try = hp.parse_hgvs_variant(re_try)
                    hgvs_re_try.posedit.edit.ref = delete
                    start_pos = str(hgvs_re_try.posedit.pos.start)
                    if re.search('\-', start_pos):
                        base, offset = start_pos.split('-')
                        new_offset = 0 - int(offset) + (len(delete))
                        end_pos = int(base)
                        hgvs_re_try.posedit.pos.end.base = int(end_pos)
                        hgvs_re_try.posedit.pos.end.offset = int(new_offset) - 1
                        not_delins = ref_ac + ':' + ref_type + '.' + start_pos + '_' + str(
                            hgvs_re_try.posedit.pos.end) + 'del' + delete + 'ins' + insert
                    elif re.search('\+', start_pos):
                        base, offset = start_pos.split('+')
                        end_pos = int(base) + (len(delete) - int(offset) - 1)
                        new_offset = 0 + int(offset) + (len(delete) - 1)
                        hgvs_re_try.posedit.pos.end.base = int(end_pos)
                        hgvs_re_try.posedit.pos.end.offset = int(new_offset)
                        not_delins = ref_ac + ':' + ref_type + '.' + start_pos + '_' + str(
                            hgvs_re_try.posedit.pos.end) + 'del' + delete + 'ins' + insert
                    else:
                        end_pos = int(start_pos) + (len(delete) - 1)
                        not_delins = ref_ac + ':' + ref_type + '.' + start_pos + '_' + str(
                            end_pos) + 'del' + delete + 'ins' + insert
                except:
                    not_delins = not_sub
                # Parse into hgvs object
                try:
                    hgvs_not_delins = hp.parse_hgvs_variant(not_delins)
                except hgvs.exceptions.HGVSError as e:
                    # Sort out multiple ALTS from VCF inputs
                    if re.search("([GATCgatc]+)>([GATCgatc]+),([GATCgatc]+)", not_delins):
                        # 						header,alts = not_delins.split('>')
                        # 						# Split up the alts into a list
                        # 						alt_list = alts.split(',')
                        # 						# Assemble and re-submit
                        # 						for alt in alt_list:
                        # 							validation['warnings'] = 'Multiple ALT sequences detected: auto-submitting all possible combinations'
                        # 							validation['write'] = 'false'
                        # 							refreshed_description = header + '>' + alt
                        # 							query = {'quibble' : refreshed_description, 'id' : validation['id'], 'warnings' : validation['warnings'], 'description' : '', 'coding' : '', 'coding_g' : '', 'genomic_r' : '', 'genomic_g' : '', 'protein' : '', 'write' : 'true', 'primary_assembly' : primary_assembly, 'order' : ordering}
                        # 							batch_list.append(query)
                        error = 'Multiple ALTs not supported by this function'
                        raise pseudoVCF2HGVSError(error)
                    else:
                        error = str(e)
                        raise pseudoVCF2HGVSError(error)

                # HGVS will deal with the errors
                hgvs_object = hgvs_not_delins
            else:
                hgvs_object = hp.parse_hgvs_variant(input)

        except Exception as e:
            error = str(e)
            raise pseudoVCF2HGVSError(error)
    else:
        hgvs_object = hp.parse_hgvs_variant(input)

    # Normalize
    hgvs_object = selected_normalizer.normalize(hgvs_object)
    # return
    return hgvs_object



# <LICENSE>
# Copyright (C) 2018  Peter Causey-Freeman, University of Leicester
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










