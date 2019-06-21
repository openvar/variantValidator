"""
A variety of functions that convert parser hgvs objects into VCF component parts
Each function has a slightly difference emphasis
"""

# Import modules
import re
import copy
from . import seq_data

# Import Biopython modules
from Bio.Seq import Seq
import hgvs
import hgvs.exceptions


# Database connections and hgvs objects are now passed from VariantValidator.py

# Error handling
class PseudoVCF2HGVSError(Exception):
    pass


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
                except hgvs.exceptions.HGVSError as e:
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

    if re.search(r'[GATC]+=', str(reverse_normalized_hgvs_genomic.posedit)):
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif 'ins' in str(reverse_normalized_hgvs_genomic.posedit) and 'del' not in str(
            reverse_normalized_hgvs_genomic.posedit):
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
    elif '>' in str(reverse_normalized_hgvs_genomic.posedit):
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif 'del' in str(reverse_normalized_hgvs_genomic.posedit) and 'ins' not in str(
            reverse_normalized_hgvs_genomic.posedit):
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
    elif 'inv' in str(reverse_normalized_hgvs_genomic.posedit):
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
    elif 'del' in str(reverse_normalized_hgvs_genomic.posedit) and 'ins' in str(
            reverse_normalized_hgvs_genomic.posedit):
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
    elif 'dup' in str(reverse_normalized_hgvs_genomic.posedit):
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

    if re.search(r'[GATC]+=', str(reverse_normalized_hgvs_genomic.posedit)):
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif 'ins' in str(reverse_normalized_hgvs_genomic.posedit) and 'del' not in str(
            reverse_normalized_hgvs_genomic.posedit):
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
    elif '>' in str(reverse_normalized_hgvs_genomic.posedit):
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif 'del' in str(reverse_normalized_hgvs_genomic.posedit) and 'ins' not in str(
            reverse_normalized_hgvs_genomic.posedit):
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
    elif 'inv' in str(reverse_normalized_hgvs_genomic.posedit):
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
        if re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif 'del' in str(reverse_normalized_hgvs_genomic.posedit) and 'ins' in str(
            reverse_normalized_hgvs_genomic.posedit):
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
    elif 'dup' in str(reverse_normalized_hgvs_genomic.posedit):
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

    if re.search(r'[GATC]+=', str(reverse_normalized_hgvs_genomic.posedit)):
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif 'ins' in str(reverse_normalized_hgvs_genomic.posedit) and 'del' not in str(
            reverse_normalized_hgvs_genomic.posedit):
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
    elif '>' in str(reverse_normalized_hgvs_genomic.posedit):
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif 'del' in str(reverse_normalized_hgvs_genomic.posedit) and 'ins' not in str(
            reverse_normalized_hgvs_genomic.posedit):
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
    elif 'inv' in str(reverse_normalized_hgvs_genomic.posedit):
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
    elif 'del' in str(reverse_normalized_hgvs_genomic.posedit) and 'ins' in str(
            reverse_normalized_hgvs_genomic.posedit):
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
    elif 'dup' in str(reverse_normalized_hgvs_genomic.posedit):
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


def hard_right_hgvs2vcf(hgvs_genomic, primary_assembly, hn, sf):
    """
    Designed specifically for gap handling.
    hard right pushes as 3 prime as possible and adds additional bases
    """
    hgvs_genomic_variant = hgvs_genomic
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    normalized_hgvs_genomic = hn.normalize(hgvs_genomic_variant)

    # Chr
    chr = seq_data.to_chr_num_ucsc(normalized_hgvs_genomic.ac, primary_assembly)
    if chr is not None:
        pass
    else:
        chr = normalized_hgvs_genomic.ac

    if re.search(r'[GATC]+=', str(normalized_hgvs_genomic.posedit)):
        pos = str(normalized_hgvs_genomic.posedit.pos.start)
        ref = normalized_hgvs_genomic.posedit.edit.ref
        alt = normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif 'ins' in str(normalized_hgvs_genomic.posedit) and 'del' not in str(
            normalized_hgvs_genomic.posedit):
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
    elif '>' in str(normalized_hgvs_genomic.posedit):
        ref = normalized_hgvs_genomic.posedit.edit.ref
        alt = normalized_hgvs_genomic.posedit.edit.alt
        pos = str(normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif 'del' in str(normalized_hgvs_genomic.posedit) and 'ins' not in str(normalized_hgvs_genomic.posedit):
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
    elif 'inv' in str(normalized_hgvs_genomic.posedit):
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
        if re.search('inv', str(normalized_hgvs_genomic.posedit)):
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif 'del' in str(normalized_hgvs_genomic.posedit) and 'ins' in str(normalized_hgvs_genomic.posedit):
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
    elif 'dup' in str(normalized_hgvs_genomic.posedit):
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
    if chr != '' and pos != '' and ref != '' and alt != '':
        # Add 2 post bases
        pos = int(pos)
        pre_end_pos = pos + len(ref)
        end_pos = pre_end_pos + 1
        post = sf.fetch_seq(str(normalized_hgvs_genomic.ac), pre_end_pos - 1, end_pos)
        ref = ref + post
        alt = alt + post

    # Dictionary the VCF
    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': normalized_hgvs_genomic}
    return vcf_dict


def hard_left_hgvs2vcf(hgvs_genomic, primary_assembly, reverse_normalizer, sf):
    """
    Designed specifically for gap handling.
    hard left pushes as 5 prime as possible and adds additional bases

    :param hgvs_genomic:
    :param primary_assembly:
    :param reverse_normalizer:
    :param sf:
    :return:
    """
    hgvs_genomic_variant = hgvs_genomic
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic_variant)

    # Chr
    chr = seq_data.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac, primary_assembly)
    if chr is not None:
        pass
    else:
        chr = reverse_normalized_hgvs_genomic.ac

    if re.search(r'[GATC]+=', str(reverse_normalized_hgvs_genomic.posedit)):
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif 'ins' in str(reverse_normalized_hgvs_genomic.posedit) and 'del' not in str(
            reverse_normalized_hgvs_genomic.posedit):
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
    elif '>' in str(reverse_normalized_hgvs_genomic.posedit):
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif 'del' in str(reverse_normalized_hgvs_genomic.posedit) and 'ins' not in str(
            reverse_normalized_hgvs_genomic.posedit):
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
    elif 'inv' in str(reverse_normalized_hgvs_genomic.posedit):
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
        if re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif 'del' in str(reverse_normalized_hgvs_genomic.posedit) and 'ins' in str(
            reverse_normalized_hgvs_genomic.posedit):
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
    elif 'dup' in str(reverse_normalized_hgvs_genomic.posedit):
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
    if chr != '' and pos != '' and ref != '' and alt != '':
        pre_pos = int(pos) - 1
        prev = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), pre_pos - 1, pre_pos)
        pos = str(pre_pos)
        ref = prev + ref
        alt = prev + alt

    # Dictionary the VCF
    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': reverse_normalized_hgvs_genomic}
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