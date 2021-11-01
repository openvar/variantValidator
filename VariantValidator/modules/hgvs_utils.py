"""
A variety of functions that convert parser hgvs objects into VCF component parts
Each function has a slightly difference emphasis
"""

# Import modules
import re
import copy
from . import seq_data
from . import utils

# Import Biopython modules
from Bio.Seq import Seq
import vvhgvs
import vvhgvs.exceptions


# Database connections and hgvs objects are now passed from VariantValidator.py

# Error handling
class PseudoVCF2HGVSError(Exception):
    pass


def vcfcp_to_hgvsstr(vcf_dict, start_hgvs):
    """
    converts  vcf components to a string hgvs variant
    :param vcf_dict:
    :return: str(hgvs_variant) with no normalization
    """
    pos = int(vcf_dict['pos'])
    ref = vcf_dict['ref']
    alt = vcf_dict['alt']
    # Generate an end position
    end = str(pos + len(ref) - 1)
    str_hgvs = "%s:%s.%s_%sdel%sins%s" % (start_hgvs.ac, start_hgvs.type, pos, end, ref, alt)
    return str_hgvs


def hgvs_to_delins_hgvs(hgvs_object, hp, hn):
    """
    :param hgvs_object: parsed hgvs string
    :param hp: hgvs_parser
    :param hn: hgvs_normalizer (check function for hn vs reverse hn rules)
    :return: hgvs_object in delins format, see if statements for the details
    """

    # Duplications (alt = ref + ref)
    if hgvs_object.posedit.edit.type == "dup":
        v_pos = hgvs_object.posedit.pos.start.base
        v_ref = hgvs_object.posedit.edit.ref
        v_alt = v_ref + v_ref

    # Insertions (Generate the ref, then alt = ref[0] + insertion + ref[1]
    if hgvs_object.posedit.edit.type == "ins":
        alt_bs = hgvs_object.posedit.edit.alt
        hgvs_object.posedit.edit.alt = ""
        hgvs_object.posedit.edit.ref = ""
        hgvs_object = hn.normalize(hgvs_object)
        hgvs_object.posedit.edit.alt = \
            hgvs_object.posedit.edit.alt[0] + \
            alt_bs + \
            hgvs_object.posedit.edit.alt[1]
        # No stringing needed, return directly
        return hgvs_object

    # Deletions (Handles simple conversion by making alt = "")
    if hgvs_object.posedit.edit.type == "del":
        hgvs_object.posedit.edit.alt = ""
        return hgvs_object

    # Create the object via a string from the vcfcp_to_hgvsstr function
    hgvs_delins_string = vcfcp_to_hgvsstr({"pos": v_pos, "ref": v_ref, "alt": v_alt}, hgvs_object)
    hgvs_delins_object = hp.parse_hgvs_variant(hgvs_delins_string)
    return hgvs_delins_object


def pvcf_to_hgvs(query, selected_assembly, normalization_direction, reverse_normalizer, validator):
    """
    :param query: pseudo_vcf string
    :param selected_assembly:
    :param normalization_direction: normalization direction an integer, 5 or 3.
    :param reverse_normalizer:
    :param validator:
    :return:
    """
    # Set normalizer
    selected_normalizer = None
    if normalization_direction == 3:
        selected_normalizer = validator.hn
    if normalization_direction == 5:
        selected_normalizer = reverse_normalizer

    # Gel stye pVCF
    query = query.replace(':', '-')
    pre_input = copy.deepcopy(query)
    vcf_elements = pre_input.split('-')

    # VCF type 1
    if re.search(r'-\d+-[GATC]+-[GATC]+', query):
        query = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], vcf_elements[3])
    elif re.search(r'-\d+-[GATC]+-', query):
        query = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], vcf_elements[2])
    else:
        raise PseudoVCF2HGVSError('Unsupported format: VCF specification 4.1 or later')

    # Chr16:2099572TC>T
    try:
        input_list = query.split(':')
        position_and_edit = input_list[1]
        if not re.match(r'N[CGWT]_', query) and not re.match(r'LRG_\d+$', query):
            chr_num = str(input_list[0])
            chr_num = chr_num.upper()
            chr_num = chr_num.strip()
            if re.match(r'CHR', chr_num):
                chr_num = chr_num.replace('CHR', '')
            # Use selected assembly
            accession = seq_data.to_accession(chr_num, selected_assembly)
            if accession is None:
                error = chr_num + ' is not part of genome build ' + selected_assembly + ' or is not supported'
                raise PseudoVCF2HGVSError(error)
        else:
            accession = input_list[0]

        # Assign reference sequence type
        ref_type = ':g.'
        if 'LRG_' in accession:
            accession = validator.db.get_refseq_id_from_lrg_id(accession)

        # Reformat the variant
        query = str(accession) + ref_type + str(position_and_edit)
    except Exception as e:
        error = str(e)
        raise PseudoVCF2HGVSError(error)

    # Find not_sub type in input e.g. GGGG>G
    not_sub = copy.deepcopy(query)
    not_sub_find = re.compile(r"([GATCgatc]+)>([GATCgatc]+)")
    if not_sub_find.search(not_sub):
        try:
            # If the length of either side of the substitution delimer (>) is >1
            matches = not_sub_find.search(not_sub)
            if len(matches.group(1)) > 1 or len(matches.group(2)) > 1 or re.search(
                    r"([GATCgatc]+)>([GATCgatc]+),([GATCgatc]+)", query):
                # Search for and remove range
                range = re.compile(r"([0-9]+)_([0-9]+)")
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
                r = re.compile(r"([0-9]+)([GATCgatc]+)")
                try:
                    m = r.search(remainder)
                    delete = m.group(2)
                    starts = posedit.split(delete)[0]
                    re_try = ref_ac + ':' + ref_type + '.' + starts + 'del' + delete[0] + 'ins' + insert
                    hgvs_re_try = validator.hp.parse_hgvs_variant(re_try)
                    hgvs_re_try.posedit.edit.ref = delete
                    start_pos = str(hgvs_re_try.posedit.pos.start)
                    if '-' in start_pos:
                        base, offset = start_pos.split('-')
                        new_offset = 0 - int(offset) + (len(delete))
                        end_pos = int(base)
                        hgvs_re_try.posedit.pos.end.base = int(end_pos)
                        hgvs_re_try.posedit.pos.end.offset = int(new_offset) - 1
                        not_delins = ref_ac + ':' + ref_type + '.' + start_pos + '_' + str(
                            hgvs_re_try.posedit.pos.end) + 'del' + delete + 'ins' + insert
                    elif '+' in start_pos:
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
                    hgvs_not_delins = validator.hp.parse_hgvs_variant(not_delins)
                except vvhgvs.exceptions.HGVSError as e:
                    # Sort out multiple ALTS from VCF inputs
                    if re.search("([GATCgatc]+)>([GATCgatc]+),([GATCgatc]+)", not_delins):
                        # header,alts = not_delins.split('>')
                        # # Split up the alts into a list
                        # alt_list = alts.split(',')
                        # # Assemble and re-submit
                        # for alt in alt_list:
                        # 	validation['warnings'] = 'Multiple ALT sequences detected: auto-submitting all possible combinations'
                        # 	validation['write'] = 'false'
                        # 	refreshed_description = header + '>' + alt
                        # 	query = {'quibble' : refreshed_description, 'id' : validation['id'], 'warnings' : validation['warnings'], 'description' : '', 'coding' : '', 'coding_g' : '', 'genomic_r' : '', 'genomic_g' : '', 'protein' : '', 'write' : 'true', 'primary_assembly' : primary_assembly, 'order' : ordering}
                        # 	batch_list.append(query)
                        error = 'Multiple ALTs not supported by this function'
                        raise PseudoVCF2HGVSError(error)
                    else:
                        error = str(e)
                        raise PseudoVCF2HGVSError(error)

                # HGVS will deal with the errors
                hgvs_object = hgvs_not_delins
            else:
                hgvs_object = validator.hp.parse_hgvs_variant(query)

        except Exception as e:
            error = str(e)
            raise PseudoVCF2HGVSError(error)
    else:
        hgvs_object = validator.hp.parse_hgvs_variant(query)

    # Normalize
    hgvs_object = selected_normalizer.normalize(hgvs_object)
    # return
    return hgvs_object


def hgvs2vcf(hgvs_genomic, primary_assembly, reverse_normalizer, sf):
    """
    Simple conversionwhich ensures identity is as 5 prime as possible by adding an extra 5
    prime base. Necessary for most gap handling situations

    :param hgvs_genomic:
    :param primary_assembly:
    :param reverse_normalizer:
    :param sf:
    :return:
    """
    hgvs_genomic_variant = hgvs_genomic
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic_variant)
    # hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

    # Chr
    chr = seq_data.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac, primary_assembly)
    if chr is not None:
        pass
    else:
        chr = reverse_normalized_hgvs_genomic.ac

    # Identity
    if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity':
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'ins':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        alt_start = start - 1  #
        # Recover sequences
        ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), alt_start, end - 1)
        ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        # Assemble
        pos = start
        ref = ref_seq
        alt = ref_seq + ins_seq

    # Substitutions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'sub':
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'del':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, start)
        # Assemble
        pos = str(start)
        ref = pre_base + hgvs_del_seq
        alt = pre_base

    # inv
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'inv':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
                # Recover sequences
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = ins_seq
        if 'inv' in str(reverse_normalized_hgvs_genomic.posedit):
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'delins':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
                # Recover sequences
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = vcf_del_seq[:1] + ins_seq

    # Duplications
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'dup':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)  #
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2  #
        start = start - 1  #
        # Recover sequences
        dup_seq = reverse_normalized_hgvs_genomic.posedit.edit.ref
        vcf_ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_ref_seq
        alt = vcf_ref_seq + dup_seq
    else:
        chr = ''
        ref = ''
        alt = ''
        pos = ''

    # ensure as 5' as possible
    if chr != '' and pos != '' and ref != '' and alt != '':
        if len(ref) > 1:
            if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity':
                pos = int(pos) - 1
                prev = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), pos - 1, pos)
                pos = str(pos)
                ref = prev + ref
                alt = prev + alt

    # Dictionary the VCF
    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': reverse_normalized_hgvs_genomic}
    return vcf_dict


def report_hgvs2vcf(hgvs_genomic, primary_assembly, reverse_normalizer, sf):
    """
    Used to report the Most true representation of the VCF i.e. 5 prime normalized but no
    additional bases added. NOTE: no gap handling capabilities

    :param hgvs_genomic:
    :param primary_assembly:
    :param reverse_normalizer:
    :param sf:
    :return:
    """

    import time
    hgvs_genomic_variant = hgvs_genomic

    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic_variant)
    # hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

    ucsc_pa = ''
    grc_pa = ''
    # Sort the primary assemblies
    if 'GRC' in primary_assembly:
        if '37' in primary_assembly:
            ucsc_pa = 'hg19'
            grc_pa = primary_assembly
        if '38' in primary_assembly:
            ucsc_pa = 'hg38'
            grc_pa = primary_assembly
    else:
        if '19' in primary_assembly:
            ucsc_pa = primary_assembly
            grc_pa = 'GRCh37'
        if '38' in primary_assembly:
            ucsc_pa = primary_assembly
            grc_pa = 'GRCh38'

    # UCSC Chr
    ucsc_chr = seq_data.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac, ucsc_pa)
    if ucsc_chr is not None:
        pass
    else:
        ucsc_chr = reverse_normalized_hgvs_genomic.ac

    # GRC Chr
    grc_chr = seq_data.to_chr_num_refseq(reverse_normalized_hgvs_genomic.ac, grc_pa)
    if grc_chr is not None:
        pass
    else:
        grc_chr = reverse_normalized_hgvs_genomic.ac

    # Identity
    if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity':
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'ins':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        alt_start = start - 1  #
        # Recover sequences
        ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), alt_start, end - 1)
        ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        # Assemble
        pos = start
        ref = ref_seq
        alt = ref_seq + ins_seq

    # Substitutions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'sub':
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'del':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, start)
        # Assemble
        pos = str(start)
        ref = pre_base + hgvs_del_seq
        alt = pre_base

    # inv
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'inv':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        # pos = str(start-1)
        # ref = bs + vcf_del_seq
        ref = vcf_del_seq
        alt = ins_seq
        if reverse_normalized_hgvs_genomic.posedit.edit.type == 'inv':
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'delins':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        # pos = str(start)
        # ref = vcf_del_seq
        # alt = vcf_del_seq[:1] + ins_seq
        pos = str(start + 1)
        ref = vcf_del_seq[1:]
        alt = ins_seq

    # Duplications
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'dup':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)  #
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2  #
        start = start - 1  #
        # Recover sequences
        vcf_ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_ref_seq[0]
        alt = vcf_ref_seq
    else:
        ref = ''
        alt = ''
        pos = ''

    # Dictionary the VCF
    vcf_dict = {'pos': str(pos), 'ref': ref, 'alt': alt, 'ucsc_chr': ucsc_chr, 'grc_chr': grc_chr,
                'normalized_hgvs': reverse_normalized_hgvs_genomic}
    return vcf_dict


def pos_lock_hgvs2vcf(hgvs_genomic, primary_assembly, reverse_normalizer, sf):
    """
    No normalization at all. No additional bases added. Simply returns an in-situ VCF
    :param hgvs_genomic:
    :param primary_assembly:
    :param reverse_normalizer:
    :param sf:
    :return:
    """
    # Replace reference manually
    if hgvs_genomic.posedit.edit.ref == '':
        hgvs_genomic.posedit.edit.ref = sf.fetch_seq(str(hgvs_genomic.ac), hgvs_genomic.posedit.pos.start.base - 1,
                                                     hgvs_genomic.posedit.pos.end.base)

    reverse_normalized_hgvs_genomic = hgvs_genomic
    if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity' and len(
            reverse_normalized_hgvs_genomic.posedit.edit.ref) == 0:
        reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(reverse_normalized_hgvs_genomic)

    # hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

    # Chr
    chr = seq_data.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac, primary_assembly)
    if chr is not None:
        pass
    else:
        chr = reverse_normalized_hgvs_genomic.ac

    # Identity
    if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity':
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'ins':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        alt_start = start - 1  #
        # Recover sequences
        ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), alt_start, end - 1)
        ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        # Assemble
        pos = start
        ref = ref_seq
        alt = ref_seq + ins_seq

    # Substitutions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'sub':
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'del':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, start)
        # Assemble
        pos = str(start)
        ref = pre_base + hgvs_del_seq
        alt = pre_base

    # inv
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'inv':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        # pos = str(start-1)
        # ref = bs + vcf_del_seq
        ref = vcf_del_seq
        alt = ins_seq
        if 'inv' in str(reverse_normalized_hgvs_genomic.posedit):
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'delins':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = vcf_del_seq[:1] + ins_seq

    # Duplications
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'dup':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)  #
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2  #
        start = start - 1  #
        # Recover sequences
        dup_seq = reverse_normalized_hgvs_genomic.posedit.edit.ref
        vcf_ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_ref_seq
        alt = vcf_ref_seq + dup_seq
    else:
        chr = ''
        ref = ''
        alt = ''
        pos = ''

    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': reverse_normalized_hgvs_genomic}
    return vcf_dict


def hard_right_hgvs2vcf(hgvs_genomic, primary_assembly, hn, reverse_normalizer, sf, tx_ac, hdp, alt_aln_method, hp, vm,
                        mrg, genomic_ac=False):
    """
    Designed specifically for gap handling.
    hard right pushes as 3 prime as possible and adds additional bases
    :param hgvs_genomic:
    :param primary_assembly:
    :param hn:
    :param reverse_normalizer:
    :param sf:
    :param tx_ac:
    :param hdp:
    :param alt_aln_method:
    :param hp:
    :param vm:
    :param genomic_ac:
    :return:
    """

    # c. must be in n. format
    try:
        hgvs_genomic = vm.c_to_n(hgvs_genomic)  # Need in n. context
    except TypeError:
        pass
    except vvhgvs.exceptions.HGVSInvalidVariantError:
        pass

    hgvs_genomic_variant = hgvs_genomic
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    normalized_hgvs_genomic = hn.normalize(hgvs_genomic_variant)

    # Variants in/on a genomic gap that cause issues need sorting by making them span so coordinates do not reverse
    if hgvs_genomic.type is not "g":
        hgvs_genomic_g = vm.n_to_g(normalized_hgvs_genomic, genomic_ac)
        try:
            hn.normalize(hgvs_genomic_g)
        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
            if "base start position must be <= end position" in str(e):
                hgvs_genomic_g_identity = copy.deepcopy(hgvs_genomic_g)
                # First, get ins and del into delins
                if hgvs_genomic_g_identity.posedit.edit.type == "dup":
                    # Duplications have no "alt" in the object structure so convert to delins
                    stb = hgvs_genomic_g_identity.posedit.pos.end.base
                    edb = hgvs_genomic_g_identity.posedit.pos.start.base
                    hgvs_genomic_g_identity.posedit.pos.start.base = stb
                    hgvs_genomic_g_identity.posedit.pos.end.base = edb
                    hgvs_genomic_g_identity.posedit.edit.ref = sf.fetch_seq(hgvs_genomic_g_identity.ac, stb - 1, edb)
                    hgvs_genomic_g_identity = hgvs_to_delins_hgvs(hgvs_genomic_g_identity, hp, hn)
                else:
                    adjust_s = hgvs_genomic_g_identity.posedit.pos.end.base
                    adjust_e = hgvs_genomic_g_identity.posedit.pos.start.base
                    hgvs_genomic_g_identity.posedit.pos.start.base = adjust_s
                    hgvs_genomic_g_identity.posedit.pos.end.base = adjust_e
                hgvs_genomic_g_identity.posedit.edit.ref = ''
                hgvs_genomic_g_identity.posedit.edit.alt = ''
                hgvs_genomic_g_identity = hn.normalize(hgvs_genomic_g_identity)
                hgvs_genomic_n_gap = vm.g_to_n(hgvs_genomic_g_identity, hgvs_genomic.ac)
                hgvs_genomic_n_identity = copy.deepcopy(hgvs_genomic_n_gap)
                hgvs_genomic_n_identity.posedit.edit.ref = ""
                hgvs_genomic_n_identity.posedit.edit.alt = ""
                hgvs_genomic_n_identity = hn.normalize(hgvs_genomic_n_identity)

                """
                At this stage we have the gap position at g. in hgvs_genomic_g_identity
                # The impact in hgvs_genomic_n_gap
                # And the range the variant needs to span to cover the gap in hgvs_genomic_n_identity
                """
                # First, get ins and del into delins
                if hgvs_genomic.posedit.edit.type == "del":
                    normalized_hgvs_genomic = hgvs_to_delins_hgvs(normalized_hgvs_genomic, hp, hn)
                elif normalized_hgvs_genomic.posedit.edit.type == "ins":
                    normalized_hgvs_genomic = hgvs_to_delins_hgvs(normalized_hgvs_genomic, hp, hn)
                elif normalized_hgvs_genomic.posedit.edit.type == "dup":
                    normalized_hgvs_genomic = hgvs_to_delins_hgvs(normalized_hgvs_genomic, hp, hn)

                # Now create the variant
                if hgvs_genomic_n_identity.posedit.pos.start.base < \
                        normalized_hgvs_genomic.posedit.pos.start.base:
                    v1 = copy.deepcopy(hgvs_genomic_n_identity)
                    v1.posedit.pos.end.base = normalized_hgvs_genomic.posedit.pos.start.base - 1
                    v1.posedit.edit.ref = ""
                    v1.posedit.edit.alt = ""
                    v1 = hn.normalize(v1)
                    if (hgvs_genomic_n_identity.posedit.pos.end.base
                            > normalized_hgvs_genomic.posedit.pos.end.base):
                        v3 = copy.deepcopy(hgvs_genomic_n_identity)
                        v3.posedit.pos.start.base = normalized_hgvs_genomic.posedit.pos.end.base + 1
                        v3.posedit.edit.ref = ""
                        v3.posedit.edit.alt = ""
                        v3 = hn.normalize(v3)
                    else:
                        v3 = False

                    # Assemble
                    hgvs_genomic_n_assembled = copy.deepcopy(hgvs_genomic_n_identity)
                    if v3 is not False:
                        hgvs_genomic_n_assembled.posedit.pos.end.base = v3.posedit.pos.end.base
                        ass_ref = (v1.posedit.edit.ref +
                                   normalized_hgvs_genomic.posedit.edit.ref +
                                   v3.posedit.edit.ref)
                        ass_alt = (v1.posedit.edit.alt +
                                   normalized_hgvs_genomic.posedit.edit.alt +
                                   v3.posedit.edit.alt)
                    else:
                        hgvs_genomic_n_assembled.posedit.pos.end.base = \
                            normalized_hgvs_genomic.posedit.pos.end.base
                        ass_ref = (v1.posedit.edit.ref +
                                   normalized_hgvs_genomic.posedit.edit.ref)
                        ass_alt = (v1.posedit.edit.alt +
                                   normalized_hgvs_genomic.posedit.edit.alt)
                    hgvs_genomic_n_assembled.posedit.edit.ref = ass_ref
                    hgvs_genomic_n_assembled.posedit.edit.alt = ass_alt
                    normalized_hgvs_genomic = copy.deepcopy(hgvs_genomic_n_assembled)
                else:
                    v3 = copy.deepcopy(hgvs_genomic_n_identity)
                    v3.posedit.pos.end.base = normalized_hgvs_genomic.posedit.pos.start.base - 1
                    v3.posedit.edit.ref = ""
                    v3.posedit.edit.alt = ""
                    if (hgvs_genomic_n_identity.posedit.pos.end.base
                            > normalized_hgvs_genomic.posedit.pos.end.base):
                        v3 = copy.deepcopy(hgvs_genomic_n_identity)
                        v3.posedit.pos.start.base = normalized_hgvs_genomic.posedit.pos.end.base + 1
                        v3.posedit.edit.ref = ""
                        v3.posedit.edit.alt = ""
                        v3 = hn.normalize(v3)
                    else:
                        v3 = False

                    # Assemble
                    hgvs_genomic_n_assembled = copy.deepcopy(hgvs_genomic_n_identity)
                    if v3 is not False:
                        hgvs_genomic_n_assembled.posedit.pos.end.base = v3.posedit.pos.end.base
                        ass_ref = (normalized_hgvs_genomic.posedit.edit.ref +
                                   v3.posedit.edit.ref)
                        ass_alt = (normalized_hgvs_genomic.posedit.edit.alt +
                                   v3.posedit.edit.alt)
                    else:
                        hgvs_genomic_n_assembled.posedit.pos.end.base = \
                            normalized_hgvs_genomic.posedit.pos.end.base
                        ass_ref = normalized_hgvs_genomic.posedit.edit.ref
                        ass_alt = normalized_hgvs_genomic.posedit.edit.alt
                    hgvs_genomic_n_assembled.posedit.edit.ref = ass_ref
                    hgvs_genomic_n_assembled.posedit.edit.alt = ass_alt
                    normalized_hgvs_genomic = copy.deepcopy(hgvs_genomic_n_assembled)

    # Chr
    chr = seq_data.to_chr_num_ucsc(normalized_hgvs_genomic.ac, primary_assembly)
    if chr is not None:
        pass
    else:
        chr = normalized_hgvs_genomic.ac

    # identity
    if normalized_hgvs_genomic.posedit.edit.type == 'identity':
        pos = str(normalized_hgvs_genomic.posedit.pos.start)
        ref = normalized_hgvs_genomic.posedit.edit.ref
        alt = normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif normalized_hgvs_genomic.posedit.edit.type == 'ins':
        end = int(normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(normalized_hgvs_genomic.posedit.pos.start.base)
        alt_start = start - 1  #
        # Recover sequences
        ref_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), alt_start, end - 1)
        ins_seq = normalized_hgvs_genomic.posedit.edit.alt
        # Assemble
        pos = start
        ref = ref_seq
        alt = ref_seq + ins_seq

    # Substitutions
    elif normalized_hgvs_genomic.posedit.edit.type == 'sub':
        ref = normalized_hgvs_genomic.posedit.edit.ref
        alt = normalized_hgvs_genomic.posedit.edit.alt
        pos = str(normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif normalized_hgvs_genomic.posedit.edit.type == 'del':
        end = int(normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), start, end)
        pre_base = sf.fetch_seq(str(normalized_hgvs_genomic.ac), adj_start, start)
        # Assemble
        pos = str(start)
        ref = pre_base + hgvs_del_seq
        alt = pre_base

    # inv
    elif normalized_hgvs_genomic.posedit.edit.type == 'inv':
        end = int(normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 1
        start = start
        try:
            ins_seq = normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        vcf_del_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        # pos = str(start-1)
        # ref = bs + vcf_del_seq
        ref = vcf_del_seq
        alt = ins_seq
        if normalized_hgvs_genomic.posedit.edit.type == 'inv':
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif normalized_hgvs_genomic.posedit.edit.type == 'delins':
        end = int(normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(normalized_hgvs_genomic.posedit.pos.start.base - 1)
        adj_start = start - 1
        start = start
        try:
            ins_seq = normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        vcf_del_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = vcf_del_seq[:1] + ins_seq

    # Duplications
    elif normalized_hgvs_genomic.posedit.edit.type == 'dup':
        end = int(normalized_hgvs_genomic.posedit.pos.end.base)  #
        start = int(normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2  #
        start = start - 1  #
        # Recover sequences
        dup_seq = normalized_hgvs_genomic.posedit.edit.ref
        vcf_ref_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_ref_seq
        alt = vcf_ref_seq + dup_seq
    else:
        chr = ''
        ref = ''
        alt = ''
        pos = ''

    # ADD SURROUNDING BASES
    # If possible, capture and alt variant that spans the gap
    merged_variant = False
    pre_merged_variant = False
    identifying_variant = False

    if chr != '' and pos != '' and ref != '' and alt != '':

        # Set exon boundary
        if genomic_ac is False:
            # Find the boundaries at the genomic level for the current exon
            exon_set = hdp.get_tx_exons(tx_ac, hgvs_genomic.ac, alt_aln_method)
            exon_end_genomic = None
            for exon in exon_set:
                if int(exon[7]) + 1 <= int(pos) <= int(exon[8]):
                    exon_end_genomic = int(exon[8])
                    break
        else:
            # Trick the system using transcript positions
            exon_set = hdp.get_tx_exons(hgvs_genomic.ac, genomic_ac, alt_aln_method)
            exon_end_genomic = None
            for exon in exon_set:
                if int(exon[5]) + 1 <= int(pos) <= int(exon[6]):
                    exon_end_genomic = int(exon[6])
                    break

        # Set loop variables for extending the push
        push_ref = ref
        push_alt = alt
        working_pos = int(pos) + len(ref)
        needs_a_push = False
        if genomic_ac is False:
            genomic_ac = hgvs_genomic.ac
        # Clear staging_loop
        staging_loop = 0

        # Loop and add bases - up to the range defined below - unless we go into an intron/past the transcript
        for push in range(50):
            post = sf.fetch_seq(str(normalized_hgvs_genomic.ac), working_pos-1, working_pos)
            push_ref = push_ref + post
            push_alt = push_alt + post

            # Create a not_delins for normalisation checking
            normlize_check_variant = hp.parse_hgvs_variant(hgvs_genomic.ac + ':' + hgvs_genomic.type + '.' +
                                                           str(pos) + '_' + str(working_pos) + 'del' + push_ref
                                                           + 'ins' + push_alt)

            # Check to see of we end up spanning a gap
            try:
                if hgvs_genomic.type != "g":
                    normlize_check_mapped = vm.n_to_g(normlize_check_variant, genomic_ac)
                else:
                    normlize_check_mapped = vm.g_to_n(normlize_check_variant, tx_ac)

            # Catch out-of-bounds errors
            except vvhgvs.exceptions.HGVSInvalidIntervalError:
                needs_a_push = False
                break

            """
            Break out from loop parameters
            """
            if normlize_check_mapped.posedit.pos.start.base > normlize_check_mapped.posedit.pos.end.base:
                needs_a_push = False
                break
            if len(normlize_check_mapped.posedit.edit.ref) <= 1:
                staging_loop = staging_loop + 1

            # Check here for the gap (Has it been crossed?) Note: if gap in tx, we have the whole gap spanned
            if (((len(normlize_check_mapped.posedit.edit.ref) != len(normlize_check_variant.posedit.edit.ref) and
                 len(normlize_check_mapped.posedit.edit.ref) > 1))
                    or
                    (normlize_check_variant.posedit.edit.type == 'identity')
                    and len(normlize_check_mapped.posedit.edit.alt) != len(normlize_check_variant.posedit.edit.ref)):

                # Add the identifying variant
                identifying_variant = normlize_check_variant

                if push == 0:  # Already crossing the gap so return original vcf
                    end_seq_check_variant = copy.copy(normlize_check_variant)
                    # end_seq_check_variant.posedit.edit.alt = end_seq_check_variant.posedit.edit.ref
                else:
                    # Look to see if the gap has been identified by addition of bases in sequence
                    end_seq_check_variant = hp.parse_hgvs_variant(hgvs_genomic.ac
                                                                  + ':'
                                                                  + hgvs_genomic.type
                                                                  + '.'
                                                                  + str(normlize_check_variant.posedit.pos.end.base - 1
                                                                        - staging_loop)
                                                                  + "_"
                                                                  + str(normlize_check_variant.posedit.pos.end.base)
                                                                  + "del"
                                                                  + push_ref[-2 - staging_loop:]
                                                                  + "ins"
                                                                  + push_ref[-2 - staging_loop])

                # Check to see of we end up spanning a gap at the last 2 bases
                if hgvs_genomic.type != "g":
                    end_seq_check_mapped = vm.n_to_g(end_seq_check_variant, genomic_ac)
                else:
                    end_seq_check_mapped = vm.g_to_n(end_seq_check_variant, tx_ac)

                # For genomic_variant mapped onto gaps, we end up with an offset
                start_offset = False
                end_offset = False
                try:
                    end_seq_check_mapped.posedit.pos.start.offset
                except AttributeError:
                    start_offset = False
                else:
                    if end_seq_check_mapped.posedit.pos.start.offset != 0:
                        start_offset = True
                try:
                    end_seq_check_mapped.posedit.pos.end.offset
                except AttributeError:
                    end_offset = False
                else:
                    if end_seq_check_mapped.posedit.pos.end.offset != 0:
                        end_offset = True
                if start_offset is True or end_offset is True:

                    # To identify the gap, we need to span it before mapping back
                    if end_offset is True:
                        end_seq_check_mapped.posedit.pos.end.base = end_seq_check_mapped.posedit.pos.start.base + 1
                        end_seq_check_mapped.posedit.pos.end.offset = 0
                        end_seq_check_mapped.posedit.edit.ref = ''
                        norml_end_seq_check_mapped = end_seq_check_mapped
                    elif start_offset is True:
                        end_seq_check_mapped.posedit.pos.start.base = end_seq_check_mapped.posedit.pos.end.base - 1
                        end_seq_check_mapped.posedit.pos.start.offset = 0
                        end_seq_check_mapped.posedit.edit.ref = ''
                        norml_end_seq_check_mapped = end_seq_check_mapped

                    # now map back onto original reference sequence
                    try:
                        norml_end_seq_check_mapped = vm.c_to_n(norml_end_seq_check_mapped)  # Need in n. context
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass
                    if hgvs_genomic.type == "g":
                        map_back = vm.n_to_g(norml_end_seq_check_mapped, genomic_ac)
                    else:
                        map_back = vm.g_to_n(norml_end_seq_check_mapped, tx_ac)

                    # Normalize variants, original and the gap induced variant (note, variant pre-normalized)
                    map_back = hn.normalize(map_back)  # gap is left so normalize right
                    map_back_rn = reverse_normalizer.normalize(map_back)
                    try:
                        map_back = vm.c_to_n(map_back)  # Need in n. context
                        map_back_rn = vm.c_to_n(map_back_rn)
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass

                    # Can the variants be normalized together
                    if ((
                            (map_back.posedit.pos.end.base >=
                             normalized_hgvs_genomic.posedit.pos.start.base - 1)
                            and
                            (map_back.posedit.pos.end.base <=
                             normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back_rn.posedit.pos.end.base >=
                                     normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.end.base <=
                                     normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back.posedit.pos.start.base >=
                                     normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back.posedit.pos.start.base <=
                                     normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back_rn.posedit.pos.start.base >=
                                     normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.start.base <=
                                     normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )):

                        # Create a variant that reflects the impact of the gap.
                        # This uses variant merging
                        # We merge the "gap" variant and the variant itself
                        v1 = hgvs_genomic
                        v2 = map_back
                        if v2.posedit.edit.type == "identity":
                            needs_a_push = True  # Return new vcf only
                            break
                        if "g" not in hgvs_genomic.type:
                            v1 = vm.n_to_g(hgvs_genomic, genomic_ac)
                            v2 = vm.n_to_g(map_back, genomic_ac)
                        try:
                            v1 = hn.normalize(v1)
                            v2 = hn.normalize(v2)
                        except vvhgvs.exceptions.HGVSInvalidVariantError:
                            needs_a_push = True  # Return new vcf only
                            break
                        else:
                            try:
                                if v1.posedit.pos.start.base < v2.posedit.pos.start.base:
                                    pre_merged_variant = mrg([v1, v2], reverse_normalizer, final_norm=False)
                                else:
                                    pre_merged_variant = mrg([v2, v1], reverse_normalizer, final_norm=False)
                                if "g" in pre_merged_variant.type:
                                    merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                else:
                                    merged_variant = pre_merged_variant
                            except utils.mergeHGVSerror as e:
                                needs_a_push = True  # Return new vcf only
                                break
                            except vvhgvs.exceptions.HGVSParseError:
                                needs_a_push = True  # Return new vcf only
                                break

                            # Ensure merged variant is not in a "non-intron" if mapped back to n.
                            if merged_variant is not False:
                                try:
                                    if (merged_variant.posedit.pos.start.offset != 0
                                            or merged_variant.posedit.pos.start.offset != 0):
                                        # Try from normalized genomic
                                        pre_merged_variant = hn.normalize(pre_merged_variant)
                                        test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                        if (test_merged_variant.posedit.pos.start.offset == 0
                                                and test_merged_variant.posedit.pos.start.offset == 0):
                                            merged_variant = pre_merged_variant
                                        else:
                                            pre_merged_variant = reverse_normalizer.normalize(pre_merged_variant)
                                            test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                            if (test_merged_variant.posedit.pos.start.offset == 0
                                                    and test_merged_variant.posedit.pos.start.offset == 0):
                                                merged_variant = pre_merged_variant
                                    # Map back to n.
                                    if "g" in merged_variant.type:
                                        merged_variant = vm.g_to_n(merged_variant, tx_ac)
                                except AttributeError:
                                    pass

                            needs_a_push = True  # Keep the new vcf
                            break
                    else:
                        needs_a_push = False  # Restore old vcf
                        break

                # Or we have identified the gap again at the expected position
                if len(end_seq_check_mapped.posedit.edit.ref) != len(end_seq_check_variant.posedit.edit.ref):

                    """
                    At this stage, we have done the following, illustrated by a  gap in transcript

                    g. NNNNNNNNNN
                    n. NNNNNNN--N

                    We forced the gap to be projected by making the end_seq_check_variant n.=

                             NN  Deletion in g.
                             |
                    g. NNNNNNNN
                    n. NNNNNNNN                   

                    So we need to make the g. == again before mapping back, which will make an ins in the n.
                    """

                    # Now normalize the variants to see if they meet
                    norml_end_seq_check_mapped = copy.deepcopy(end_seq_check_mapped)
                    norml_end_seq_check_mapped.posedit.edit.alt = norml_end_seq_check_mapped.posedit.edit.ref

                    # now map back onto original reference sequence
                    try:
                        norml_end_seq_check_mapped = vm.c_to_n(norml_end_seq_check_mapped)  # Need in n. context
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass
                    if hgvs_genomic.type == "g":
                        map_back = vm.n_to_g(norml_end_seq_check_mapped, genomic_ac)
                    else:
                        map_back = vm.g_to_n(norml_end_seq_check_mapped, tx_ac)

                    # In transcript gaps, this can push us fully into the gap
                    try:
                        if map_back.posedit.pos.start.offset != 0 and map_back.posedit.pos.start.offset != 0:
                            needs_a_push = False
                            break
                    except AttributeError:
                        pass

                    # Normalize variants, original and the gap induced variant (note, variant pre-normalized)
                    try:
                        map_back = hn.normalize(map_back)  # gap is left so normalize right
                        map_back_rn = reverse_normalizer.normalize(map_back)
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        needs_a_push = False  # Restore old vcf
                        break
                    try:
                        map_back = vm.c_to_n(map_back)  # Need in n. context
                        map_back_rn = vm.c_to_n(map_back_rn)
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass

                    # Is the gap variant the same as the incoming variant?
                    if normalized_hgvs_genomic == map_back:
                        needs_a_push = True
                        push_ref = end_seq_check_variant.posedit.edit.ref
                        push_alt = end_seq_check_variant.posedit.edit.alt
                        pos = end_seq_check_variant.posedit.pos.start.base
                        break

                    # Can the variants be normalized together
                    if ((
                            (map_back.posedit.pos.end.base >=
                             normalized_hgvs_genomic.posedit.pos.start.base - 1)
                            and
                            (map_back.posedit.pos.end.base <=
                             normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back_rn.posedit.pos.end.base >=
                                     normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.end.base <=
                                     normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back.posedit.pos.start.base >=
                                     normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back.posedit.pos.start.base <=
                                     normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back_rn.posedit.pos.start.base >=
                                     normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.start.base <=
                                     normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )):

                        # Create a variant that reflects the impact of the gap.
                        # This uses variant merging
                        # We merge the "gap" variant and the variant itself
                        v1 = hgvs_genomic
                        v2 = map_back
                        if v2.posedit.edit.type == "identity":
                            needs_a_push = True  # Return new vcf only
                            break
                        if "g" not in hgvs_genomic.type:
                            v1 = vm.n_to_g(hgvs_genomic, genomic_ac)
                            v2 = vm.n_to_g(map_back, genomic_ac)
                        try:
                            v1 = hn.normalize(v1)
                            v2 = hn.normalize(v2)
                        except vvhgvs.exceptions.HGVSInvalidVariantError:
                            needs_a_push = True  # Return new vcf only
                            break
                        else:
                            try:
                                if v1.posedit.pos.start.base < v2.posedit.pos.start.base:
                                    pre_merged_variant = mrg([v1, v2], reverse_normalizer, final_norm=False)
                                else:
                                    pre_merged_variant = mrg([v2, v1], reverse_normalizer, final_norm=False)
                                if "g" in pre_merged_variant.type:
                                    merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                else:
                                    merged_variant = pre_merged_variant
                            except utils.mergeHGVSerror as e:
                                needs_a_push = True  # Return new vcf only
                                break
                            except vvhgvs.exceptions.HGVSParseError:
                                needs_a_push = True  # Return new vcf only
                                break

                            # Ensure merged variant is not in a "non-intron" if mapped back to n.
                            if merged_variant is not False:
                                try:
                                    if (merged_variant.posedit.pos.start.offset != 0
                                            or merged_variant.posedit.pos.start.offset != 0):
                                        # Try from normalized genomic
                                        pre_merged_variant = hn.normalize(pre_merged_variant)
                                        test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                        if (test_merged_variant.posedit.pos.start.offset == 0
                                                and test_merged_variant.posedit.pos.start.offset == 0):
                                            merged_variant = pre_merged_variant
                                        else:
                                            pre_merged_variant = reverse_normalizer.normalize(pre_merged_variant)
                                            test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                            if (test_merged_variant.posedit.pos.start.offset == 0
                                                    and test_merged_variant.posedit.pos.start.offset == 0):
                                                merged_variant = pre_merged_variant
                                    # Map back to n.
                                    if "g" in merged_variant.type:
                                        merged_variant = vm.g_to_n(merged_variant, tx_ac)
                                except AttributeError:
                                    pass

                            needs_a_push = True  # Keep the new vcf
                            break
                    else:
                        needs_a_push = False  # Restore old vcf
                        break

                else:
                    # Everything missed, assume no push required
                    needs_a_push = False
                    break

            # exon boundary hit. Break before intron
            elif working_pos == exon_end_genomic:
                break

            # Continue looping
            else:
                working_pos = working_pos + 1
                continue

        # Clear staging_loop
        staging_loop = 0

        # Create vcf dict
        if needs_a_push is True:
            # Re-sep pos-ref-alt (pos remains equal)
            ref = push_ref
            alt = push_alt
        # else:
        #     # Old behavior
        #     # Hard addition of 2 post bases for sequence mismatches in gaps
        #     pos = int(pos)
        #     pre_end_pos = pos + len(ref)
        #     end_pos = pre_end_pos + 1
        #     if pre_end_pos <= exon_end_genomic:
        #         post = sf.fetch_seq(str(normalized_hgvs_genomic.ac), pre_end_pos - 1, end_pos)
        #         ref = ref + post
        #         alt = alt + post

    # Dictionary the VCF
    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': normalized_hgvs_genomic,
                'merged_variant': merged_variant, 'identifying_variant': identifying_variant,
                'pre_merged_variant': pre_merged_variant}
    str_hgvs = vcfcp_to_hgvsstr(vcf_dict, hgvs_genomic)
    vcf_dict['str_hgvs'] = str_hgvs
    vcf_dict['needs_a_push'] = needs_a_push
    return vcf_dict


def hard_left_hgvs2vcf(hgvs_genomic, primary_assembly, hn, reverse_normalizer, sf, tx_ac, hdp, alt_aln_method,
                       hp, vm, mrg, genomic_ac=False):
    """
    Designed specifically for gap handling.
    hard left pushes as 5 prime as possible and adds additional bases
    :param hgvs_genomic:
    :param primary_assembly:
    :param hn:
    :param reverse_normalizer:
    :param sf:
    :param tx_ac:
    :param hdp:
    :param alt_aln_method:
    :param hp:
    :param vm:
    :param genomic_ac:
    :return:
    """
    # c. must be in n. format
    try:
        hgvs_genomic = vm.c_to_n(hgvs_genomic)  # Need in n. context
    except TypeError:
        pass
    except vvhgvs.exceptions.HGVSInvalidVariantError:
        pass

    hgvs_genomic_variant = hgvs_genomic

    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic_variant)

    # Variants in/on a genomic gap that cause issues need sorting by making them span so coordinates do not reverse
    if hgvs_genomic.type is not "g":
        hgvs_genomic_g = vm.n_to_g(reverse_normalized_hgvs_genomic, genomic_ac)
        try:
            hn.normalize(hgvs_genomic_g)
        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
            if "base start position must be <= end position" in str(e):
                hgvs_genomic_g_identity = copy.deepcopy(hgvs_genomic_g)
                # First, get ins and del into delins
                if hgvs_genomic_g_identity.posedit.edit.type == "dup":
                    # Duplications have no "alt" in the object structure so convert to delins
                    stb = hgvs_genomic_g_identity.posedit.pos.end.base
                    edb = hgvs_genomic_g_identity.posedit.pos.start.base
                    hgvs_genomic_g_identity.posedit.pos.start.base = stb
                    hgvs_genomic_g_identity.posedit.pos.end.base = edb
                    hgvs_genomic_g_identity.posedit.edit.ref = sf.fetch_seq(hgvs_genomic_g_identity.ac, stb - 1, edb)
                    hgvs_genomic_g_identity = hgvs_to_delins_hgvs(hgvs_genomic_g_identity, hp, hn)
                else:
                    adjust_s = hgvs_genomic_g_identity.posedit.pos.end.base
                    adjust_e = hgvs_genomic_g_identity.posedit.pos.start.base
                    hgvs_genomic_g_identity.posedit.pos.start.base = adjust_s
                    hgvs_genomic_g_identity.posedit.pos.end.base = adjust_e
                hgvs_genomic_g_identity.posedit.edit.ref = ''
                hgvs_genomic_g_identity.posedit.edit.alt = ''
                hgvs_genomic_g_identity = hn.normalize(hgvs_genomic_g_identity)
                hgvs_genomic_n_gap = vm.g_to_n(hgvs_genomic_g_identity, hgvs_genomic.ac)
                hgvs_genomic_n_identity = copy.deepcopy(hgvs_genomic_n_gap)
                hgvs_genomic_n_identity.posedit.edit.ref = ""
                hgvs_genomic_n_identity.posedit.edit.alt = ""
                hgvs_genomic_n_identity = hn.normalize(hgvs_genomic_n_identity)

                """
                At this stage we have the gap position at g. in hgvs_genomic_g_identity
                # The impact in hgvs_genomic_n_gap
                # And the range the variant needs to span to cover the gap in hgvs_genomic_n_identity
                """
                # First, get ins and del into delins
                if reverse_normalized_hgvs_genomic.posedit.edit.type == "del":
                    reverse_normalized_hgvs_genomic = hgvs_to_delins_hgvs(reverse_normalized_hgvs_genomic, hp, hn)
                elif reverse_normalized_hgvs_genomic.posedit.edit.type == "ins":
                    reverse_normalized_hgvs_genomic = hgvs_to_delins_hgvs(reverse_normalized_hgvs_genomic, hp, hn)
                elif reverse_normalized_hgvs_genomic.posedit.edit.type == "dup":
                    reverse_normalized_hgvs_genomic = hgvs_to_delins_hgvs(reverse_normalized_hgvs_genomic, hp, hn)

                # Now create the variant
                if hgvs_genomic_n_identity.posedit.pos.start.base < \
                        reverse_normalized_hgvs_genomic.posedit.pos.start.base:
                    v1 = copy.deepcopy(hgvs_genomic_n_identity)
                    v1.posedit.pos.end.base = reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1
                    v1.posedit.edit.ref = ""
                    v1.posedit.edit.alt = ""
                    v1 = hn.normalize(v1)
                    if (hgvs_genomic_n_identity.posedit.pos.end.base
                            > reverse_normalized_hgvs_genomic.posedit.pos.end.base):
                        v3 = copy.deepcopy(hgvs_genomic_n_identity)
                        v3.posedit.pos.start.base = reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1
                        v3.posedit.edit.ref = ""
                        v3.posedit.edit.alt = ""
                        v3 = hn.normalize(v3)
                    else:
                        v3 = False

                    # Assemble
                    hgvs_genomic_n_assembled = copy.deepcopy(hgvs_genomic_n_identity)
                    if v3 is not False:
                        hgvs_genomic_n_assembled.posedit.pos.end.base = v3.posedit.pos.end.base
                        ass_ref = (v1.posedit.edit.ref +
                                   reverse_normalized_hgvs_genomic.posedit.edit.ref +
                                   v3.posedit.edit.ref)
                        ass_alt = (v1.posedit.edit.alt +
                                   reverse_normalized_hgvs_genomic.posedit.edit.alt +
                                   v3.posedit.edit.alt)
                    else:
                        hgvs_genomic_n_assembled.posedit.pos.end.base = \
                            reverse_normalized_hgvs_genomic.posedit.pos.end.base
                        ass_ref = (v1.posedit.edit.ref +
                                   reverse_normalized_hgvs_genomic.posedit.edit.ref)
                        ass_alt = (v1.posedit.edit.alt +
                                   reverse_normalized_hgvs_genomic.posedit.edit.alt)
                    hgvs_genomic_n_assembled.posedit.edit.ref = ass_ref
                    hgvs_genomic_n_assembled.posedit.edit.alt = ass_alt
                    reverse_normalized_hgvs_genomic = copy.deepcopy(hgvs_genomic_n_assembled)
                else:
                    v3 = copy.deepcopy(hgvs_genomic_n_identity)
                    v3.posedit.pos.end.base = reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1
                    v3.posedit.edit.ref = ""
                    v3.posedit.edit.alt = ""
                    if (hgvs_genomic_n_identity.posedit.pos.end.base
                            > reverse_normalized_hgvs_genomic.posedit.pos.end.base):
                        v3 = copy.deepcopy(hgvs_genomic_n_identity)
                        v3.posedit.pos.start.base = reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1
                        v3.posedit.edit.ref = ""
                        v3.posedit.edit.alt = ""
                        v3 = hn.normalize(v3)
                    else:
                        v3 = False

                    # Assemble
                    hgvs_genomic_n_assembled = copy.deepcopy(hgvs_genomic_n_identity)
                    if v3 is not False:
                        hgvs_genomic_n_assembled.posedit.pos.end.base = v3.posedit.pos.end.base
                        ass_ref = (reverse_normalized_hgvs_genomic.posedit.edit.ref +
                                   v3.posedit.edit.ref)
                        ass_alt = (reverse_normalized_hgvs_genomic.posedit.edit.alt +
                                   v3.posedit.edit.alt)
                    else:
                        hgvs_genomic_n_assembled.posedit.pos.end.base = \
                            reverse_normalized_hgvs_genomic.posedit.pos.end.base
                        ass_ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
                        ass_alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
                    hgvs_genomic_n_assembled.posedit.edit.ref = ass_ref
                    hgvs_genomic_n_assembled.posedit.edit.alt = ass_alt
                    reverse_normalized_hgvs_genomic = copy.deepcopy(hgvs_genomic_n_assembled)

    # Chr
    chr = seq_data.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac, primary_assembly)
    if chr is not None:
        pass
    else:
        chr = reverse_normalized_hgvs_genomic.ac

    # Identity
    if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity':
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'ins':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        alt_start = start - 1  #
        # Recover sequences
        ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), alt_start, end - 1)
        ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        # Assemble
        pos = start
        ref = ref_seq
        alt = ref_seq + ins_seq

    # Substitutions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'sub':
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'del':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, start)
        # Assemble
        pos = str(start)
        ref = pre_base + hgvs_del_seq
        alt = pre_base

    # inv
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'inv':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        # pos = str(start-1)
        # ref = bs + vcf_del_seq
        ref = vcf_del_seq
        alt = ins_seq
        if reverse_normalized_hgvs_genomic.posedit.edit.type == 'inv':
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'delins':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = vcf_del_seq[:1] + ins_seq

    # Duplications
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'dup':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)  #
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2  #
        start = start - 1  #
        # Recover sequences
        dup_seq = reverse_normalized_hgvs_genomic.posedit.edit.ref
        vcf_ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_ref_seq
        alt = vcf_ref_seq + dup_seq
    else:
        chr = ''
        ref = ''
        alt = ''
        pos = ''

    # ADD SURROUNDING BASES
    # If possible, capture and alt variant that spans the gap
    merged_variant = False
    pre_merged_variant = False
    identifying_variant = False

    if chr != '' and pos != '' and ref != '' and alt != '':

        # Set exon boundary
        if genomic_ac is False:
            # Find the boundaries at the genomic level for the current exon
            exon_set = hdp.get_tx_exons(tx_ac, hgvs_genomic.ac, alt_aln_method)
            exon_start_genomic = None
            for exon in exon_set:
                if int(exon[7]) + 1 <= int(pos) <= int(exon[8]):
                    exon_start_genomic = int(exon[7] + 1)
                    break
        else:
            # Trick the system using transcript positions
            exon_set = hdp.get_tx_exons(hgvs_genomic.ac, genomic_ac, alt_aln_method)
            exon_start_genomic = None
            for exon in exon_set:
                if int(exon[5]) + 1 <= int(pos) <= int(exon[6]):
                    exon_start_genomic = int(exon[5] + 1)
                    break

        # Set loop variables for extending the push
        push_ref = ref
        push_alt = alt
        push_pos_by = 1
        needs_a_push = False
        staging_loop = 0
        if genomic_ac is False:
            genomic_ac = hgvs_genomic.ac

        # Loop and add bases - up to the range defined below - unless we go into an intron/past the transcript
        for push in range(50):
            pre_pos = int(pos) - push_pos_by
            prev = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), pre_pos - 1, pre_pos)
            push_ref = prev + push_ref
            push_alt = prev + push_alt

            # Create a not_delins for normalisation checking
            var_end = str(pre_pos + len(push_ref) - 1)
            normlize_check_variant = hp.parse_hgvs_variant(hgvs_genomic.ac + ':' + hgvs_genomic.type + '.' +
                                                           str(pre_pos) + '_' + var_end + 'del' + push_ref
                                                           + 'ins' + push_alt)
            # Check to see of we end up spanning a gap
            try:
                if hgvs_genomic.type != "g":
                    normlize_check_mapped = vm.n_to_g(normlize_check_variant, genomic_ac)
                else:
                    normlize_check_mapped = vm.g_to_n(normlize_check_variant, tx_ac)
            # Catch out-of-bounds errors
            except vvhgvs.exceptions.HGVSInvalidIntervalError:
                needs_a_push = False
                break

            """
            Break out from loop parameters
            """
            if normlize_check_mapped.posedit.pos.start.base > normlize_check_mapped.posedit.pos.end.base:
                needs_a_push = False
                break
            if len(normlize_check_mapped.posedit.edit.ref) <= 1:
                staging_loop = staging_loop + 1

            # Check here for the gap (Has it been crossed?) Note: if gap in tx, we have the whole gap spanned
            if (((len(normlize_check_mapped.posedit.edit.ref) != len(normlize_check_variant.posedit.edit.ref) and
                 len(normlize_check_mapped.posedit.edit.ref) > 1))
                    or
                    (normlize_check_variant.posedit.edit.type == 'identity')
                    and len(normlize_check_mapped.posedit.edit.alt) != len(normlize_check_variant.posedit.edit.ref)):

                # Add the identifying variant
                identifying_variant = normlize_check_variant

                if push == 0:  # Already crossing the gap so return original vcf
                    end_seq_check_variant = copy.deepcopy(normlize_check_variant)
                    # end_seq_check_variant.posedit.edit.alt = end_seq_check_variant.posedit.edit.ref

                else:
                    # Look to see if the gap has been identified by addition of bases in sequence
                    end_seq_check_variant = hp.parse_hgvs_variant(hgvs_genomic.ac
                                                                  + ':'
                                                                  + hgvs_genomic.type
                                                                  + '.'
                                                                  + str(normlize_check_variant.posedit.pos.start.base)
                                                                  + "_"
                                                                  + str(normlize_check_variant.posedit.pos.start.base
                                                                        + 1 + staging_loop)
                                                                  + "del"
                                                                  + push_ref[0:2 + staging_loop]
                                                                  + "ins"
                                                                  + push_ref[0:2 + staging_loop])

                # Check to see of we end up spanning a gap at the last 2 bases
                if hgvs_genomic.type != "g":
                    end_seq_check_mapped = vm.n_to_g(end_seq_check_variant, genomic_ac)
                else:
                    end_seq_check_mapped = vm.g_to_n(end_seq_check_variant, tx_ac)

                # For genomic_variant mapped onto gapps, we end up with an offset
                start_offset = False
                end_offset = False
                try:
                    end_seq_check_mapped.posedit.pos.start.offset
                except AttributeError:
                    start_offset = False
                else:
                    if end_seq_check_mapped.posedit.pos.start.offset != 0:
                        start_offset = True
                try:
                    end_seq_check_mapped.posedit.pos.end.offset
                except AttributeError:
                    end_offset = False
                else:
                    if end_seq_check_mapped.posedit.pos.end.offset != 0:
                        end_offset = True
                if start_offset is True or end_offset is True:

                    # To identify the gap, we need to span it before mapping back
                    if end_offset is True:
                        end_seq_check_mapped.posedit.pos.end.base = end_seq_check_mapped.posedit.pos.start.base + 1
                        end_seq_check_mapped.posedit.pos.end.offset = 0
                        end_seq_check_mapped.posedit.edit.ref = ''
                        norml_end_seq_check_mapped = end_seq_check_mapped
                    elif start_offset is True:
                        end_seq_check_mapped.posedit.pos.start.base = end_seq_check_mapped.posedit.pos.end.base - 1
                        end_seq_check_mapped.posedit.pos.start.offset = 0
                        end_seq_check_mapped.posedit.edit.ref = ''
                        norml_end_seq_check_mapped = end_seq_check_mapped

                    # now map back onto original reference sequence
                    try:
                        norml_end_seq_check_mapped = vm.c_to_n(norml_end_seq_check_mapped)  # Need in n. context
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass
                    if hgvs_genomic.type == "g":
                        map_back = vm.n_to_g(norml_end_seq_check_mapped, genomic_ac)
                    else:
                        map_back = vm.g_to_n(norml_end_seq_check_mapped, tx_ac)

                    # Normalize variants, original and the gap induced variant (note, variant pre-normalized)
                    if map_back.posedit.pos.start.base > map_back.posedit.pos.end.base:
                        needs_a_push = False
                        break
                    map_back = hn.normalize(map_back)  # gap is left so normalize right
                    map_back_rn = reverse_normalizer.normalize(map_back)
                    try:
                        map_back = vm.c_to_n(map_back)  # Need in n. context
                        map_back_rn = vm.c_to_n(map_back_rn)
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass

                    # Can the variants be normalized together
                    if ((
                            (map_back.posedit.pos.end.base >=
                             reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                            and
                            (map_back.posedit.pos.end.base <=
                             reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back_rn.posedit.pos.end.base >=
                                     reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.end.base <=
                                     reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back.posedit.pos.start.base >=
                                     reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back.posedit.pos.start.base <=
                                     reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back_rn.posedit.pos.start.base >=
                                     reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.start.base <=
                                     reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )):

                        # Create a variant that reflects the impact of the gap.
                        # This uses variant merging
                        # We merge the "gap" variant and the variant itself
                        v1 = hgvs_genomic
                        v2 = map_back
                        # if v2.posedit.edit.type == "identity":
                        #     needs_a_push = True  # Return new vcf only
                        #     push_pos_by = push_pos_by + 1
                        #     break
                        if "g" not in hgvs_genomic.type:
                            v1 = vm.n_to_g(hgvs_genomic, genomic_ac)
                            v2 = vm.n_to_g(map_back, genomic_ac)
                        try:
                            v1 = reverse_normalizer.normalize(v1)
                            v2 = reverse_normalizer.normalize(v2)
                        except vvhgvs.exceptions.HGVSInvalidVariantError:
                            needs_a_push = True  # Restore old vcf
                            push_pos_by = push_pos_by + 1
                            break
                        else:
                            try:
                                if v1.posedit.pos.start.base < v2.posedit.pos.start.base:
                                    pre_merged_variant = mrg([v1, v2], reverse_normalizer, final_norm=False)
                                else:
                                    pre_merged_variant = mrg([v2, v1], reverse_normalizer, final_norm=False)
                                if "g" in pre_merged_variant.type:
                                    merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                else:
                                    merged_variant = pre_merged_variant
                            except utils.mergeHGVSerror as e:
                                needs_a_push = True  # Return new vcf only
                                push_pos_by = push_pos_by + 1
                                break
                            except vvhgvs.exceptions.HGVSParseError as e:
                                needs_a_push = True  # Return new vcf only
                                push_pos_by = push_pos_by + 1
                                break

                            # Ensure merged variant is not in a "non-intron" if mapped back to n.
                            if merged_variant is not False:
                                try:
                                    if (merged_variant.posedit.pos.start.offset != 0
                                            or merged_variant.posedit.pos.start.offset != 0):
                                        # Try from normalized genomic
                                        pre_merged_variant = hn.normalize(pre_merged_variant)
                                        test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                        if (test_merged_variant.posedit.pos.start.offset == 0
                                                and test_merged_variant.posedit.pos.start.offset == 0):
                                            merged_variant = pre_merged_variant
                                        else:
                                            pre_merged_variant = reverse_normalizer.normalize(pre_merged_variant)
                                            test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                            if (test_merged_variant.posedit.pos.start.offset == 0
                                                    and test_merged_variant.posedit.pos.start.offset == 0):
                                                merged_variant = pre_merged_variant
                                    # Map back to n.
                                    if "g" in merged_variant.type:
                                        merged_variant = vm.g_to_n(merged_variant, tx_ac)
                                except AttributeError:
                                    pass

                            needs_a_push = True  # Keep the new vcf
                            push_pos_by = push_pos_by + 1
                            break
                    else:
                        needs_a_push = False  # Restore old vcf
                        push_pos_by = push_pos_by + 1
                        break

                # Or we have identified the gap again at the expected position
                if len(end_seq_check_mapped.posedit.edit.ref) != len(end_seq_check_variant.posedit.edit.ref):

                    """
                    At this stage, we have done the following, illustrated by a  gap in transcript
                    
                    g. NNNNNNNNNN
                    n. NNNNNNN--N
                    
                    We forced the gap to be projected by making the end_seq_check_variant n.=
                    
                             NN  Deletion in g.
                             |
                    g. NNNNNNNN
                    n. NNNNNNNN                   
                    
                    So we need to make the g. == again before mapping back, which will make an ins in the n.
                    """

                    # Now normalize the variants to see if they meet
                    norml_end_seq_check_mapped = copy.deepcopy(end_seq_check_mapped)
                    norml_end_seq_check_mapped.posedit.edit.alt = norml_end_seq_check_mapped.posedit.edit.ref

                    # now map back onto original reference sequence
                    try:
                        norml_end_seq_check_mapped = vm.c_to_n(norml_end_seq_check_mapped)  # Need in n. context
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass

                    if hgvs_genomic.type == "g":
                        map_back = vm.n_to_g(norml_end_seq_check_mapped, genomic_ac)
                    else:
                        map_back = vm.g_to_n(norml_end_seq_check_mapped, tx_ac)

                    # In transcript gaps, this can push us fully into the gap
                    try:
                        if map_back.posedit.pos.start.offset != 0 and map_back.posedit.pos.start.offset != 0:
                            needs_a_push = False
                            push_pos_by = push_pos_by + 1
                            break
                    except AttributeError:
                        pass

                    # Normalize variants, original and the gap induced variant (note, variant pre-normalized)
                    if map_back.posedit.pos.start.base > map_back.posedit.pos.end.base:
                        needs_a_push = False
                        break
                    try:
                        map_back = hn.normalize(map_back)  # gap is left so normalize right
                        map_back_rn = reverse_normalizer.normalize(map_back)
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        needs_a_push = False  # Restore old vcf
                        push_pos_by = push_pos_by + 1
                        break
                    try:
                        map_back = vm.c_to_n(map_back)  # Need in n. context
                        map_back_rn = vm.c_to_n(map_back_rn)
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass

                    # Is the gap variant the same as the incoming variant?
                    if reverse_normalized_hgvs_genomic == map_back_rn:
                        needs_a_push = True
                        push_pos_by = 1
                        push_ref = end_seq_check_variant.posedit.edit.ref
                        push_alt = end_seq_check_variant.posedit.edit.alt
                        pos = end_seq_check_variant.posedit.pos.start.base
                        break

                    # Can the variants be normalized together
                    if ((
                            (map_back.posedit.pos.end.base >=
                             reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                            and
                            (map_back.posedit.pos.end.base <=
                             reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back_rn.posedit.pos.end.base >=
                                     reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.end.base <=
                                     reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back.posedit.pos.start.base >=
                                     reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back.posedit.pos.start.base <=
                                     reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back_rn.posedit.pos.start.base >=
                                     reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.start.base <=
                                     reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )):

                        # Create a variant that reflects the impact of the gap.
                        # This uses variant merging
                        # We merge the "gap" variant and the variant itself
                        v1 = hgvs_genomic
                        v2 = map_back
                        # if v2.posedit.edit.type == "identity":
                        #     needs_a_push = True  # Return new vcf only
                        #     push_pos_by = push_pos_by + 1
                        #     break

                        if "g" not in hgvs_genomic.type:
                            v1 = vm.n_to_g(hgvs_genomic, genomic_ac)
                            v2 = vm.n_to_g(map_back, genomic_ac)

                        try:
                            v1 = reverse_normalizer.normalize(v1)
                            v2 = reverse_normalizer.normalize(v2)
                        except vvhgvs.exceptions.HGVSInvalidVariantError:
                            needs_a_push = True  # Restore old vcf
                            push_pos_by = push_pos_by + 1
                            break
                        else:
                            try:
                                if v1.posedit.pos.start.base < v2.posedit.pos.start.base:
                                    pre_merged_variant = mrg([v1, v2], reverse_normalizer, final_norm=False)
                                else:
                                    pre_merged_variant = mrg([v2, v1], reverse_normalizer, final_norm=False)
                                if "g" in pre_merged_variant.type:
                                    merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                else:
                                    merged_variant = pre_merged_variant
                            except utils.mergeHGVSerror as e:
                                needs_a_push = True  # Return new vcf only
                                push_pos_by = push_pos_by + 1
                                break
                            except vvhgvs.exceptions.HGVSParseError as e:
                                needs_a_push = True  # Return new vcf only
                                break

                            # Ensure merged variant is not in a "non-intron" if mapped back to n.
                            if merged_variant is not False:
                                try:
                                    if (merged_variant.posedit.pos.start.offset != 0
                                            or merged_variant.posedit.pos.start.offset != 0):
                                        # Try from normalized genomic
                                        pre_merged_variant = hn.normalize(pre_merged_variant)
                                        test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                        if (test_merged_variant.posedit.pos.start.offset == 0
                                                and test_merged_variant.posedit.pos.start.offset == 0):
                                            merged_variant = pre_merged_variant
                                        else:
                                            pre_merged_variant = reverse_normalizer.normalize(pre_merged_variant)
                                            test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                            if (test_merged_variant.posedit.pos.start.offset == 0
                                                    and test_merged_variant.posedit.pos.start.offset == 0):
                                                merged_variant = pre_merged_variant
                                    # Map back to n.
                                    if "g" in merged_variant.type:
                                        merged_variant = vm.g_to_n(merged_variant, tx_ac)
                                except AttributeError:
                                    pass

                            needs_a_push = True  # Keep the new vcf
                            push_pos_by = push_pos_by + 1
                            break
                    else:
                        needs_a_push = False  # Restore old vcf
                        push_pos_by = push_pos_by + 1
                        break

                else:
                    # Everything missed, assume no push required
                    needs_a_push = False
                    push_pos_by = push_pos_by + 1
                    break

            # exon boundary hit. Break before intron
            elif pre_pos == exon_start_genomic:
                push_pos_by = push_pos_by + 1
                break

            # Continue looping
            else:
                push_pos_by = push_pos_by + 1
                continue

        # Clear staging_loop
        staging_loop = 0

        # Populate vcf dict
        if needs_a_push is True:
            # Re-sep pos-ref-alt
            pos = (int(pos) - (push_pos_by-1))
            ref = push_ref
            alt = push_alt
        # else:
        #     # Old behavior
        #     # Hard addition of 2 post bases for sequence mismatches in gaps
        #     pre_pos = int(pos) - 1
        #     if pre_pos >= exon_start_genomic:
        #         prev = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), pre_pos - 1, pre_pos)
        #         pos = str(pre_pos)
        #         ref = prev + ref
        #         alt = prev + alt

    # Dictionary the VCF
    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': reverse_normalized_hgvs_genomic,
                'merged_variant': merged_variant, 'identifying_variant': identifying_variant,
                'pre_merged_variant': pre_merged_variant}
    str_hgvs = vcfcp_to_hgvsstr(vcf_dict, hgvs_genomic)
    vcf_dict['str_hgvs'] = str_hgvs
    vcf_dict['needs_a_push'] = needs_a_push
    return vcf_dict


def hgvs_ref_alt(hgvs_variant, sf):
    if re.search(r'[GATC]+=', str(hgvs_variant.posedit)):
        ref = hgvs_variant.posedit.edit.ref
        alt = hgvs_variant.posedit.edit.ref

    # Insertions
    elif 'ins' in str(hgvs_variant.posedit) and 'del' not in str(hgvs_variant.posedit):
        end = int(hgvs_variant.posedit.pos.end.base)
        start = int(hgvs_variant.posedit.pos.start.base)
        alt_start = start - 1  #
        # Recover sequences
        ref_seq = sf.fetch_seq(str(hgvs_variant.ac), alt_start, end)
        ins_seq = hgvs_variant.posedit.edit.alt
        # Assemble
        ref = ref_seq
        alt = ref_seq[:1] + ins_seq + ref_seq[-1:]

    # Substitutions
    elif '>' in str(hgvs_variant.posedit):
        ref = hgvs_variant.posedit.edit.ref
        alt = hgvs_variant.posedit.edit.alt

    # Deletions
    elif 'del' in str(hgvs_variant.posedit) and 'ins' not in str(hgvs_variant.posedit):
        ref = hgvs_variant.posedit.edit.ref
        alt = ''

    # inv
    elif 'inv' in str(hgvs_variant.posedit):
        ref = hgvs_variant.posedit
        my_seq = Seq(ref)
        alt = str(my_seq.reverse_complement())

    # Delins
    elif 'del' in str(hgvs_variant.posedit) and 'ins' in str(hgvs_variant.posedit):
        ref = hgvs_variant.posedit.edit.ref
        alt = hgvs_variant.posedit.edit.alt

    # Duplications
    elif 'dup' in str(hgvs_variant.posedit):
        ref = hgvs_variant.posedit.edit.ref
        alt = hgvs_variant.posedit.edit.ref + hgvs_variant.posedit.edit.ref
    else:
        ref = ''
        alt = ''

    ref_alt_dict = {'ref': ref, 'alt': alt}
    return ref_alt_dict

# <LICENSE>
# Copyright (C) 2016-2021 VariantValidator Contributors
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
