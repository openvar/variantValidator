from Bio import Entrez,SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import httplib2 as http
import json
from urlparse import urlparse #Python 2
import functools
import traceback
import sys
from vvLogging import logger
import re
import copy
import mysql

#from urllib.parse import urlparse #Python 3

def handleCursor(func):
    #Decorator function for handling opening and closing cursors.
    @functools.wraps(func)
    def wrapper(self,*args,**kwargs):
        self.db.pool=mysql.connector.pooling.MySQLConnectionPool(pool_size=10, **self.db.dbConfig)
        self.db.conn=self.db.pool.get_connection()
        self.db.cursor = self.db.conn.cursor(buffered=True)
        out=func(self,*args,**kwargs)
        if self.db.cursor:
            self.db.cursor.close()
        if self.db.conn:
            self.db.conn.close()
        #self.cursor=None
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
    uri = 'http://rest.genenames.org'
    target = urlparse(uri + path)
    method = 'GET'
    body = ''
    h = http.Http()
    # collect the response
    response, content = h.request(
        target.geturl(),
        method,
        body,
        headers)
    if response['status'] == '200':
        # assume that content is a json reply
        # parse content with the json module
        data['record'] = json.loads(content)
    else:
        data['error'] = "Unable to contact the HGNC database: Please try again later"
    return data

# method for final validation and stringifying parsed hgvs variants prior to printing/passing to html
def valstr(hgvs_variant):
    """
    Function to ensure the required number of reference bases are displayed in descriptions
    """
    cp_hgvs_variant = copy.deepcopy(hgvs_variant)
    if cp_hgvs_variant.posedit.edit.type == 'identity':
        if len(cp_hgvs_variant.posedit.edit.ref) > 1:
            cp_hgvs_variant = remove_reference(cp_hgvs_variant)
        cp_hgvs_variant = str(cp_hgvs_variant)
    else:
        cp_hgvs_variant = remove_reference(cp_hgvs_variant)
        cp_hgvs_variant = str(cp_hgvs_variant)
    return cp_hgvs_variant

# From output_formatter
"""
format protein description into single letter aa code
"""
def single_letter_protein(hgvs_protein):
    return hgvs_protein.format({'p_3_letter': False})
"""
format nucleotide descriptions to not display reference base
"""
def remove_reference(hgvs_nucleotide):
    hgvs_nucleotide_refless = hgvs_nucleotide.format({'max_ref_length': 0})
    return hgvs_nucleotide_refless

def exceptPass(validation=None):
    exc_type, exc_value, last_traceback = sys.exc_info()
    te = traceback.format_exc()
    tbk = [str(exc_type), str(exc_value), str(te)]
    er = str('\n'.join(tbk))
    if last_traceback:
        logger.warning(
            "Except pass for " + str(exc_type) + " " + str(exc_value) + " at line " + str(last_traceback.tb_lineno))
    else:
        logger.warning("Except pass for " + str(exc_type) + " " + str(exc_value))
    logger.debug(er)

# From functions.py
"""
user_input
collect the input from the form and convert to a hgvs readable string
    Removes brackets and contained information -if given
    Identifies variant type (p. c. etc)
    Returns a dictionary containing a formated input string which is optimal for hgvs 
    parsing and the variant type
    Accepts c, g, n, r currently. And now P also 15.07.15
"""
def user_input(input):
    raw_variant = input.strip()

    # Set regular expressions for if statements
    pat_g = re.compile("\:g\.")  # Pattern looks for :g.
    pat_gene = re.compile('\(.+?\)')  # Pattern looks for (....)
    pat_c = re.compile("\:c\.")  # Pattern looks for :c.
    pat_r = re.compile("\:r\.")  # Pattern looks for :r.
    pat_n = re.compile("\:n\.")  # Pattern looks for :n.
    pat_p = re.compile("\:p\.")  # Pattern looks for :p.
    pat_m = re.compile("\:m\.")  # Pattern looks for :m.
    pat_est = re.compile("\d\:\d")  # Pattern looks for number:number

    # If statements
    if pat_g.search(raw_variant):  # If the :g. pattern is present in the raw_variant, g_in is linked to the raw_variant
        if pat_gene.search(raw_variant):  # If pat gene is present in the raw_variant
            variant = pat_gene.sub('',
                                   raw_variant)  # variant is set to the raw_variant string with the pattern (...) substituted out
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

# From links.py
"""
Function which predicts the protein effect of c. inversions
"""

def pro_inv_info(prot_ref_seq, prot_var_seq):
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
    else:
        # Deal with terminations
        term = re.compile("\*")
        if term.search(prot_var_seq):
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
                        if var[aa_counter] == '\*':
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

def pro_delins_info(prot_ref_seq, prot_var_seq):
    info = {
            'variant' : 'true',
            'prot_del_seq' : '',
            'prot_ins_seq' : '',
            'edit_start' : 0,
            'edit_end' : 0,
            'terminate' : 'false',
            'ter_pos' : 0,
            'error' : 'false'
            }

    # Is there actually any variation?
    if prot_ref_seq == prot_var_seq:
        info['variant'] = 'false'
    else:
        # Deal with terminations
        term = re.compile("\*")
        if term.search(prot_var_seq):
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
                        if var[aa_counter] == '\*':
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
#                   if len(ref) < len(var):
#                       ref.append('*')
#                       if prot_var_seq[-1] == '*':
#                           var.append('*')

                    # the sequences should now be the same length
                    # Except if the ter was removed
                    if len(ref) > len(var):
                        info['error'] = 'true'
                        return info
                    else:
                        # Enter the sequences
                        info['prot_del_seq'] = ''.join(ref)
                        info['prot_ins_seq'] = ''.join(var)
                        info['edit_end'] = info['edit_start'] + len(ref) -1
                        return info

"""
Translate c. reference sequences, including those that have been modified 
must have the CDS in the specified position
"""
def translate(ed_seq, cds_start):
    # ed_seq = ed_seq.replace('\n', '')
    ed_seq = ed_seq.strip()
    # Ensure the starting codon is in the correct position
    met = ed_seq[cds_start:cds_start + 3]
    if (met == 'ATG') or (met == 'atg'):
        # Remove the 5 prime UTR
        sequence = ed_seq[cds_start:]
        coding_dna = Seq(str(sequence), IUPAC.unambiguous_dna)
        # Translate
        trans = coding_dna.translate()
        aain = list(trans)
        aaout = []
        count = 0
        while aain:
            if aain[count] != '*':
                aaout.append(aain[count])
                count = count + 1
            else:
                aaout.append(aain[count])
                break
        translation = ''.join(aaout)
        # Apply a width of 60 characters to the string output
        # translation = textwrap.fill(translation, width=60)
        return translation
    else:
        translation = 'error'
        return translation

"""
Convert single letter amino acid code to 3 letter code
"""
def one_to_three(seq):
    aacode = {
        'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu',
        'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
        'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
        'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
        'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
        '*': 'Ter'}

    oned = list(seq)
    out = []
    for aa in oned:
        get_value = aacode.get(aa)
        out.append(get_value)

    threed_up = ''.join(out)

    return threed_up


""" 
Takes a reference sequence and inverts the specified position
"""
# n. Inversions - This comes from VariantValidator, not validation!!!!
def n_inversion(ref_seq, del_seq, inv_seq, interval_start, interval_end):
    sequence = ''
    # Use string indexing to check whether the sequences are the same
    test = ref_seq[interval_start - 1:interval_end]
    if test == del_seq:
        sequence = ref_seq[0:interval_start - 1] + inv_seq + ref_seq[interval_end:]
        return sequence
    else:
        sequence = 'error'
        return sequence


# Custom Exceptions
class VariantValidatorError(Exception):
    pass
