import requests
import functools
import logging
import re
import copy
from VariantValidator.modules import seq_data

logger = logging.getLogger(__name__)

# Translation table derived from Extended translation table in
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec27
PROT_TRANSLATION_DICT = {
        'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
        'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
        'TGT':'C', 'TGC':'C',
        'TGG': 'W',
        'GAA':'E', 'GAG':'E',
        'GAT':'D', 'GAC':'D',
        'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
        'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
        'AAT':'N', 'AAC':'N',
        'ATG':'M',
        'AAA':'K', 'AAG':'K',
        'TAT':'Y', 'TAC':'Y',
        'ATT':'I', 'ATC':'I', 'ATA':'I',
        'CAA':'Q', 'CAG':'Q',
        'TTT':'F', 'TTC':'F',
        'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
        'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
        'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
        'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
        'CAT':'H', 'CAC':'H',
        'TAA':'*', 'TAG':'*', 'TGA':'*'
        }

PROT_TRANSLATION_DICT_SEL = copy.copy(PROT_TRANSLATION_DICT)
PROT_TRANSLATION_DICT_SEL['TGA'] = 'U'



def handleCursor(func):
    """
    Decorator function for handling opening and closing cursors.
    """
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        out = func(self, *args, **kwargs)
        return out
    return wrapper


def hgnc_rest(path):
    data = {
        'record': '',
        'error': 'false'
    }
    # HGNC server
    headers = {
        'Accept': 'application/json',
    }
    domain = 'http://rest.genenames.org'
    url = domain + path
    r = requests.get(url, headers=headers)

    if r.status_code == 200:
        data['record'] = r.json()
    else:
        my_error = "Problem encountered while connecting genenames.org: URL=%s: Status=%s" % (url, str(r.status_code))
        data['error'] = my_error
    return data


def ensembl_rest(id, endpoint, genome, options=False):
    """
    fires requests to the Ensembl APIs
    :param id: Usually a transcript ID (base accession minus the version)
    :param endpoint: See https://rest.ensembl.org/
    :param genome: Genome build, grch37 or grch38
    :param options: set of options for additional data, see https://rest.ensembl.org/
    :return: json of the requested data
    """
    data = {
        'record': '',
        'error': 'false'
    }

    id = id.split('.')[0]

    # Set base URL
    if genome == 'GRCh37':
        base_url = 'https://grch37.rest.ensembl.org'
    if genome == 'GRCh38':
        base_url = 'https://rest.ensembl.org'

    headers = {
        'Accept': 'application/json',
    }

    if options is False:
        url = '%s%s%s?content-type=application/json' % (base_url, endpoint, id)
    else:
        url = '%s%s%s?%s;content-type=application/json' % (base_url, endpoint, id, options)

    # Request info
    r = requests.get(url, headers=headers)

    if r.status_code == 200:
        data['record'] = r.json()
    else:
        my_error = "Problem encountered while connecting Ensembl REST: URL=%s: Status=%s" % (url, str(r.status_code))
        data['error'] = my_error
    return data


def ensembl_tark(id, endpoint, options=False):
    """
    fires requests to the Ensembl APIs
    :param id: Usually a transcript ID
    :param endpoint: See http://dev-tark.ensembl.org/
    :param options: set of options for additional data, see http://dev-tark.ensembl.org/
    :return: json of the requested data
    """
    data = {
        'record': '',
        'error': 'false'
    }

    # Set base URL
    base_url = 'https://tark.ensembl.org'


    headers = {
        'Accept': 'application/json',
    }

    if options is False:
        url = '%s%s?stable_id_with_version=%s&content-type=application/json' % (base_url, endpoint, id)
    else:
        url = '%s%s?stable_id_with_version=%s&content-type=application/json' % (base_url, endpoint, id)

    # Request info
    try:
        r = requests.get(url, headers=headers)
    except requests.exceptions.InvalidSchema as e:
        my_error = str(e)
        data['error'] = my_error
        return data

    if r.status_code == 200:
        data['record'] = r.json()
    else:
        my_error = "Problem encountered while connecting Ensembl REST: URL=%s: Status=%s" % (url, str(r.status_code))
        data['error'] = my_error
    return data


def valstr(hgvs_variant):
    """
    Required for final validation and stringifying parsed hgvs variants prior to printing/passing to html.
    Function to ensure the required number of reference bases are displayed in descriptions
    """
    cp_hgvs_variant = copy.deepcopy(hgvs_variant)
    if cp_hgvs_variant.posedit.edit.type == 'identity':
        if len(cp_hgvs_variant.posedit.edit.ref) > 0:
            cp_hgvs_variant = remove_reference(cp_hgvs_variant)
        cp_hgvs_variant = str(cp_hgvs_variant)
    else:
        cp_hgvs_variant = remove_reference(cp_hgvs_variant)
        cp_hgvs_variant = str(cp_hgvs_variant)
    return cp_hgvs_variant


def single_letter_protein(hgvs_protein):
    """
    format protein description into single letter aa code
    """
    return hgvs_protein.format({'p_3_letter': False})


def remove_reference(hgvs_nucleotide):
    """
    format nucleotide descriptions to not display reference base
    """
    hgvs_nucleotide_refless = hgvs_nucleotide.format({'max_ref_length': 0})
    return hgvs_nucleotide_refless


def remove_reference_string(variant_string):

    """
    format stringified nucleotide descriptions to not display reference base
    """
    # deletions
    del_match = re.search('del[GATC]+$', variant_string)
    if del_match:
        variant_string = variant_string.replace(del_match[0], 'del')
    # inversions
    inv_match = re.search('inv[GATC]+$', variant_string)
    if inv_match:
        variant_string = variant_string.replace(inv_match[0], 'inv')
    # delins
    delins_match = re.search('del[GATC]+ins[GATC+]$', variant_string)
    if delins_match:
        delins_match_b = delins_match[0].split('ins')[1]
        delins_match_b = 'delins' + delins_match_b
        variant_string = variant_string.replace(delins_match[0], delins_match_b)
    return variant_string


def user_input(query):
    """
    user_input
    collect the input from the form and convert to a hgvs readable string
        Removes brackets and contained information -if given
        Identifies variant type (p. c. etc)
        Returns a dictionary containing a formated input string which is optimal for hgvs
        parsing and the variant type
        Accepts c, g, n, r currently. And now P also 15.07.15
    """
    raw_variant = query.strip()

    # Set regular expressions for if statements
    pat_g = re.compile(r":g\.")  # Pattern looks for :g.
    pat_gene = re.compile(r'\(.+?\)')  # Pattern looks for (....)
    pat_c = re.compile(r":c\.")  # Pattern looks for :c.
    pat_r = re.compile(r":r\.")  # Pattern looks for :r.
    pat_n = re.compile(r":n\.")  # Pattern looks for :n.
    pat_p = re.compile(r":p\.")  # Pattern looks for :p.
    pat_m = re.compile(r":m\.")  # Pattern looks for :m.
    pat_est = re.compile(r"\d:\d")  # Pattern looks for number:number

    # If statements
    if pat_g.search(raw_variant):  # If the :g. pattern is present in the raw_variant, g_in is linked to the raw_variant
        if pat_gene.search(raw_variant):  # If pat gene is present in the raw_variant
            variant = pat_gene.sub('', raw_variant)  # variant is set to the raw_variant string with the pattern (...) substituted out
            formated = {'variant': variant, 'type': ':g.'}
            return formated
        else:
            variant = raw_variant  # Otherwise it is set to raw_variant
            formated = {'variant': variant, 'type': ':g.'}
            return formated

    elif pat_r.search(raw_variant):
        if pat_gene.search(raw_variant):
            variant = pat_gene.sub('', raw_variant)
            formated = {'variant': variant, 'type': ':r.'}
            return formated
        else:
            variant = raw_variant
            formated = {'variant': variant, 'type': ':r.'}
            return formated

    elif pat_n.search(raw_variant):
        if pat_gene.search(raw_variant):
            variant = pat_gene.sub('', raw_variant)
            formated = {'variant': variant, 'type': ':n.'}
            return formated
        else:
            variant = raw_variant
            formated = {'variant': variant, 'type': ':n.'}
            return formated

    elif pat_c.search(raw_variant):
        if pat_gene.search(raw_variant):
            variant = pat_gene.sub('', raw_variant)
            formated = {'variant': variant, 'type': ':c.'}
            return formated
        else:
            variant = raw_variant
            formated = {'variant': variant, 'type': ':c.'}
            return formated

    elif pat_p.search(raw_variant):
        variant = raw_variant
        formated = {'variant': variant, 'type': ':p.'}
        return formated

    elif pat_m.search(raw_variant):
        variant = raw_variant
        formated = {'variant': variant, 'type': ':m.'}
        return formated

    elif pat_est.search(raw_variant):
        variant = raw_variant
        formated = {'variant': variant, 'type': 'est'}
        return formated

    else:
        formatted = 'invalid'
        return formatted


def pro_inv_info(prot_ref_seq, prot_var_seq):
    """
    Function which predicts the protein effect of c. inversions
    """
    info = {
        'variant': 'true',
        'prot_del_seq': '',
        'prot_ins_seq': '',
        'edit_start': 0,
        'edit_end': 0,
        'terminate': 'false',
        'ter_pos': 0,
        'error': 'false'
    }

    # Is there actually any variation?
    if prot_ref_seq == prot_var_seq:
        info['variant'] = 'false'
        info['variant'] = 'identity'
        return info
    else:
        # Deal with terminations
        if '*' in prot_var_seq:
            # Set the termination reporter to true
            info['terminate'] = 'true'
            # The termination position will be equal to the length of the variant sequence because it's a TERMINATOR!!!
            info['ter_pos'] = len(prot_var_seq)
            # cut the ref sequence to == size
            prot_ref_seq = prot_ref_seq[0:info['ter_pos']]
            prot_var_seq = prot_var_seq[0:info['ter_pos']]

            # Whether terminated or not, the sequences should now be the same length
            # Unless the termination codon has been disrupted
            if len(prot_var_seq) < len(prot_ref_seq):
                info['error'] = 'true'
                return info
            else:
                # Set the counter
                aa_counter = 0

                # Make list copies of the sequences to gather the required info
                ref = list(prot_ref_seq)
                var = list(prot_var_seq)

                # Loop through ref list to find the first missmatch position
                for aa in ref:
                    if ref[aa_counter] == var[aa_counter]:
                        aa_counter = aa_counter + 1
                    else:
                        break

                # Enter the start position
                info['edit_start'] = aa_counter + 1
                # Remove those elements form the list
                del ref[0:aa_counter]
                del var[0:aa_counter]

                # the sequences should now be the same length
                # Except if the termination codon was removed
                if len(ref) > len(var):
                    info['error'] = 'true'
                    return info
                else:
                    # Reset the aa_counter but to go backwards
                    aa_counter = 0
                    # reverse the lists
                    ref = ref[::-1]
                    var = var[::-1]
                    # Reverse loop through ref list to find the first missmatch position
                    for aa in ref:
                        if var[aa_counter] == r'\*':
                            break
                        if aa == var[aa_counter]:
                            aa_counter = aa_counter + 1
                        else:
                            break
                    # Remove those elements form the list
                    del ref[0:aa_counter]
                    del var[0:aa_counter]
                    # re-reverse the lists
                    ref = ref[::-1]
                    var = var[::-1]

                    # If the var is > ref, the ter has been removed, need to re-add ter to each
                    if len(ref) < len(var):
                        ref.append('*')
                        if prot_var_seq[-1] == '*':
                            var.append('*')
                    # the sequences should now be the same length
                    # Except if the ter was removed
                    if len(ref) > len(var):
                        info['error'] = 'true'
                        return info
                    else:
                        # Enter the sequences
                        info['prot_del_seq'] = ''.join(ref)
                        info['prot_ins_seq'] = ''.join(var)
                        info['edit_end'] = info['edit_start'] + len(ref) - 1
                        return info


def pro_delins_info(prot_ref_seq, prot_var_seq, in_frame=False):
    info = {
            'variant': 'true',
            'prot_del_seq': '',
            'prot_ins_seq': '',
            'edit_start': 0,
            'edit_end': 0,
            'terminate': 'false',
            'ter_pos': 0,
            'error': 'false'
            }

    # Is there actually any variation?
    if prot_ref_seq == prot_var_seq:
        info['variant'] = 'false'
        info['variant'] = 'identity'
        return info
    else:
        # Deal with terminations (Cannot be used as a marker for the delins pathway because in frame deletions have Ter
        if '*' in prot_var_seq:
            # Set the termination reporter to true
            info['terminate'] = 'true'

            # Set the terminal pos dependant on the shortest sequence
            # This is where we look for in-frame deletions / delins that can be shortened to a simple del/delins
            if len(prot_var_seq) <= len(prot_ref_seq):

                # Look for early termination rather than just deletions. These params may need to be altered.
                if in_frame is not False and in_frame == (len(prot_var_seq) - len(prot_ref_seq)):
                    info['ter_pos'] = len(prot_ref_seq)

                else:
                    # This code deals with the early termination out of frame variants
                    if prot_var_seq[-1] == "*":
                        info['ter_pos'] = len(prot_var_seq)
                    # Otherwise, if no termination, we carry on as normal
                    else:
                        info['ter_pos'] = len(prot_ref_seq)
            else:
                info['ter_pos'] = len(prot_var_seq)

            # cut the ref sequence to == size
            prot_ref_seq = prot_ref_seq[0:info['ter_pos']]
            prot_var_seq = prot_var_seq[0:info['ter_pos']]

            # Set the counter
            aa_counter = 0
            # Make list copies of the sequences to gather the required info
            ref = list(prot_ref_seq)
            var = list(prot_var_seq)
            # Loop through ref list to find the first missmatch position
            for aa in ref:
                if ref[aa_counter] == var[aa_counter]:
                    aa_counter = aa_counter + 1
                else:
                    break

            # Enter the start position
            info['edit_start'] = aa_counter + 1

            # Remove those elements form the list
            del ref[0:aa_counter]
            del var[0:aa_counter]

            # Reset the aa_counter but to go backwards
            aa_counter = 0
            # reverse the lists
            ref = ref[::-1]
            var = var[::-1]
            # Reverse loop through ref list to find the first missmatch position
            for aa in ref:
                try:
                    if var[aa_counter] == r'\*':
                        break
                except IndexError:
                    break
                if aa == var[aa_counter]:
                    aa_counter = aa_counter + 1
                else:
                    break
            # Remove those elements form the list
            del ref[0:aa_counter]
            del var[0:aa_counter]
            # re-reverse the lists
            ref = ref[::-1]
            var = var[::-1]

            # Enter the sequences
            info['prot_del_seq'] = ''.join(ref)
            info['prot_ins_seq'] = ''.join(var)
            info['edit_end'] = info['edit_start'] + len(ref) - 1
            return info


def translate(ed_seq, cds_start, modified_aa=None, tolerate_no_stop_cds = False):
    """
    Translate c. reference sequences, including those that have been modified
    must have the CDS in the specified position
    """
    ed_seq = ed_seq.strip()
    # Ensure the starting codon is in the correct position
    met = ed_seq[cds_start:cds_start + 3]
    met = met.upper() #this should be redundant with all inputs upper case
    """
    >>> mito_table.start_codons
    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']
    """
    if met not in ['ATG', 'TTG', 'CTG', 'GTG', 'ATT', 'ATC', 'ATA', 'ACG']:
        translation = 'error'
        return translation

    # Remove the 5 prime UTR
    coding_sequence = ed_seq[cds_start:].upper()
    if modified_aa == "Sec":
        use_dict = PROT_TRANSLATION_DICT_SEL
        stops = ['TAA', 'TAG']
    else:
        use_dict = PROT_TRANSLATION_DICT
        stops = ['TAA', 'TAG', 'TGA']

    # Translate
    if len(coding_sequence) % 3:
        last_codon_end = int(len(coding_sequence)/3) * 3
    else:
        last_codon_end = len(coding_sequence)
    codon_list = [coding_sequence[i:i+3] for i in range(0, last_codon_end, 3)]
    translation = []
    for codon in codon_list:
        translation.append(use_dict[codon])
        if codon in stops:
            break
    if translation[-1] != '*':
        if not tolerate_no_stop_cds:
            raise IndexError('No stop CDS')
        translation.append('X')

    return "".join(translation)


def one_to_three(seq):
    """
    Convert single letter amino acid code to 3 letter code
    """
    aacode = {
        'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu',
        'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
        'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
        'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
        'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
        '*': 'Ter', 'U': 'Sec'}

    oned = list(seq)
    out = []
    for aa in oned:
        get_value = aacode.get(aa)
        out.append(get_value)

    threed_up = ''.join(out)
    return threed_up


def three_to_one(seq):

    aacode = {
        'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E',
        'Phe': 'F', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Lys': 'K', 'Leu': 'L', 'Met': 'M', 'Asn': 'N',
        'Pro': 'P', 'Gln': 'Q', 'Arg': 'R', 'Ser': 'S',
        'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y',
        'Ter': '*'}

    threed = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    out = []

    for aa in threed:
        get_value = aacode.get(aa)
        out.append(get_value)

    oned_up = ''.join(out)
    return oned_up


# n. Inversions - This comes from VariantValidator, not validation!!!!
def n_inversion(ref_seq, del_seq, inv_seq, interval_start, interval_end):
    """
    Takes a reference sequence and inverts the specified position
    """
    # Use string indexing to check whether the sequences are the same
    test = ref_seq[interval_start - 1:interval_end]
    if test == del_seq:
        sequence = ref_seq[0:interval_start - 1] + inv_seq + ref_seq[interval_end:]
        return sequence
    else:
        sequence = 'error'
        return sequence


def hgvs_dup2indel(hgvs_seq):
    """Will convert hgvs variant object dup into a string with del and ins"""
    string = "%s:%s.%s_%sdel%sins%s%s" % (
        hgvs_seq.ac,
        hgvs_seq.type,
        hgvs_seq.posedit.pos.start.base,
        hgvs_seq.posedit.pos.end.base,
        hgvs_seq.posedit.edit.ref,
        hgvs_seq.posedit.edit.ref,
        hgvs_seq.posedit.edit.ref
        )
    return string


def get_exon_boundary_list(variant, validator):
    """
    Function to get the exon boundaries of a transcript
    """
    # Get the transcript
    transcript = variant.quibble.split(':')[0]
    if transcript.startswith('NM_' or 'NR_' or 'ENST'):
        # Get alignment options and identify the relevant primary assembly chrom
        mapping_options = validator.hdp.get_tx_mapping_options(transcript)
        chromosome_reference = None
        for option in mapping_options:
            is_in_assembly = seq_data.to_chr_num_refseq(option[1], variant.primary_assembly)
            if is_in_assembly is not None:
                chromosome_reference = option[1]
                break

        # Set the offsets for CDS
        transcript_info = validator.hdp.get_tx_identity_info(transcript)
        try:
            cds_offset = int(transcript_info[3]) + 1
        except ValueError:
            cds_offset = 0
        except TypeError:
            cds_offset = 0
        try:
            cds_end = int(transcript_info[4])
        except ValueError:
            cds_end = 0
        except TypeError:
            cds_end = 0

        # Get the exon boundaries of the transcript
        exons = validator.hdp.get_tx_exons(transcript, chromosome_reference, validator.alt_aln_method)

        # Extract the exon boundaries
        exon_boundaries = []
        exon_boundaries.append(transcript)
        exon_boundaries.append(chromosome_reference)

        for exon in exons:
            # Set the exon boundaries
            if exon[5] > cds_end:
                exon_boundaries.append(f"*{str(exon[5]+1 - cds_end)}")
            else:
                exon_boundaries.append(str(exon[5]+1 - (cds_offset-1)))
            if exon[6] > cds_end:
                exon_boundaries.append(f"*{str(exon[6] - cds_end)}")
            else:
                exon_boundaries.append(str(exon[6] - (cds_offset-1)))

        return exon_boundaries
    else:
        raise ExonMappingError(f"{transcript} is not a valid transcript reference sequence ID")


# Custom Exceptions
class VariantValidatorError(Exception):
    pass


class mergeHGVSerror(Exception):
    pass


class alleleVariantError(Exception):
    pass


class DatabaseConnectionError(Exception):
    pass


class ObsoleteSeqError(Exception):
    pass


class ExonMappingError(Exception):
    pass

# <LICENSE>
# Copyright (C) 2016-2025 VariantValidator Contributors
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
