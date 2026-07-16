import requests
import functools
import logging
import re
import copy
from VariantValidator.modules import seq_data
import time

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

DNA_TRANS_TBL = str.maketrans("ACTG", "TGAC")

def simple_dna_revcomp(dna):
    """
    Simplest possible reverse compliment, for use on validated input and
    DNA (or cDNA) of known origin.

    Multiple online published performance benchmarks put this method at
    the top of the performance comparison for simple DNA reverse
    compliment, VS loop based and dict lookup among others.
    Biopython does(did?) the same internally, but adds extra checks.
    """
    return dna.upper().translate(DNA_TRANS_TBL)[::-1]

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
    """
    Fires requests to the HGNC REST API.
    """

    data = {
        'record': '',
        'error': 'false'
    }

    headers = {
        'Accept': 'application/json',
    }

    domain = 'http://rest.genenames.org'
    url = domain + path

    max_retries = 3
    retry_status = {429, 500, 502, 503, 504}

    last_status = None

    for attempt in range(max_retries):
        try:
            r = requests.get(url, headers=headers, timeout=60)
            last_status = r.status_code

            if r.status_code == 200:
                data['record'] = r.json()
                return data

            if r.status_code not in retry_status:
                data['error'] = (
                    "Problem encountered while connecting genenames.org: "
                    "URL=%s: Status=%s" % (url, r.status_code)
                )
                logger.warning(data['error'])
                return data

        except (
            requests.exceptions.Timeout,
            requests.exceptions.ConnectionError,
            requests.exceptions.SSLError,
        ) as e:
            if attempt == max_retries - 1:
                data['error'] = (
                    "Problem encountered while connecting genenames.org: "
                    "URL=%s: %s" % (url, str(e))
                )
                logger.warning(data['error'])
                return data

        if attempt < max_retries - 1:
            time.sleep(0.5 * (2 ** attempt))

    data['error'] = (
        "Problem encountered while connecting genenames.org: "
        "URL=%s: Status=%s" % (
            url,
            last_status if last_status is not None else "Unknown"
        )
    )
    logger.warning(data['error'])
    return data


def ensembl_rest(id, endpoint, genome, options=False):
    """
    fires requests to the Ensembl APIs
    :param id: Usually a transcript ID (base accession minus the version)
    :param endpoint: See https://rest.ensembl.org/
    :param genome: Genome build, grch37 or grch38
    :param options: set of options for additional data
    :return: json of the requested data
    """

    data = {
        'record': '',
        'error': 'false'
    }

    id = id.split('.')[0]

    if genome == 'GRCh37':
        base_url = 'https://grch37.rest.ensembl.org'
    elif genome == 'GRCh38':
        base_url = 'https://rest.ensembl.org'
    else:
        data['error'] = "Unknown genome build '%s'" % genome
        logger.warning(data['error'])
        return data

    headers = {
        'Accept': 'application/json',
    }

    if options:
        url = '%s%s%s?%s;content-type=application/json' % (
            base_url,
            endpoint,
            id,
            options
        )
    else:
        url = '%s%s%s?content-type=application/json' % (
            base_url,
            endpoint,
            id
        )

    max_retries = 3
    retry_status = {429, 500, 502, 503, 504}

    last_status = None

    for attempt in range(max_retries):
        try:
            r = requests.get(url, headers=headers, timeout=60)
            last_status = r.status_code

            if r.status_code == 200:
                data['record'] = r.json()
                return data

            if r.status_code not in retry_status:
                data['error'] = (
                    "Problem encountered while connecting Ensembl REST: "
                    "URL=%s: Status=%s" % (url, r.status_code)
                )
                logger.warning(data['error'])
                return data

        except (
            requests.exceptions.Timeout,
            requests.exceptions.ConnectionError,
            requests.exceptions.SSLError,
        ) as e:
            if attempt == max_retries - 1:
                data['error'] = (
                    "Problem encountered while connecting Ensembl REST: "
                    "URL=%s: %s" % (url, str(e))
                )
                logger.warning(data['error'])
                return data

        if attempt < max_retries - 1:
            time.sleep(0.5 * (2 ** attempt))

    data['error'] = (
        "Problem encountered while connecting Ensembl REST: "
        "URL=%s: Status=%s" % (
            url,
            last_status if last_status is not None else "Unknown"
        )
    )
    logger.warning(data['error'])
    return data


def ensembl_tark(id, endpoint, options=False):
    """
    fires requests to the Ensembl TARK API
    """

    data = {
        'record': '',
        'error': 'false'
    }

    base_url = 'https://tark.ensembl.org'

    headers = {
        'Accept': 'application/json',
    }

    if options:
        url = (
            "%s%s?stable_id_with_version=%s&%s&content-type=application/json"
            % (base_url, endpoint, id, options)
        )
    else:
        url = (
            "%s%s?stable_id_with_version=%s&content-type=application/json"
            % (base_url, endpoint, id)
        )

    max_retries = 3
    retry_status = {429, 500, 502, 503, 504}

    last_status = None

    for attempt in range(max_retries):
        try:
            r = requests.get(url, headers=headers, timeout=60)
            last_status = r.status_code

            if r.status_code == 200:
                data['record'] = r.json()
                return data

            if r.status_code not in retry_status:
                data['error'] = (
                    "Problem encountered while connecting Ensembl REST: "
                    "URL=%s: Status=%s" % (url, r.status_code)
                )
                logger.warning(data['error'])
                return data

        except (
            requests.exceptions.InvalidSchema,
            requests.exceptions.Timeout,
            requests.exceptions.ConnectionError,
            requests.exceptions.SSLError,
        ) as e:
            if isinstance(e, requests.exceptions.InvalidSchema):
                data['error'] = str(e)
                return data

            if attempt == max_retries - 1:
                data['error'] = (
                    "Problem encountered while connecting Ensembl REST: "
                    "URL=%s: %s" % (url, str(e))
                )
                logger.warning(data['error'])
                return data

        if attempt < max_retries - 1:
            time.sleep(0.5 * (2 ** attempt))

    data['error'] = (
        "Problem encountered while connecting Ensembl REST: "
        "URL=%s: Status=%s" % (
            url,
            last_status if last_status is not None else "Unknown"
        )
    )
    logger.warning(data['error'])
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
    logger.info("pro_inv_info function called")
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
    logger.info(f"pro_delins_info function called")
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


def translate(ed_seq, cds_start, modified_aa=None, tolerate_no_stop_cds=False):
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
            # Add Polyadenylation stop codon completing bases to relevant
            # transcripts
            spare_end = len(coding_sequence) % 3
            if spare_end and coding_sequence[-spare_end:] in ['T','TA']:
                translation.append('*')
            else:
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
        mapping_options = variant.map_dat.mapping_options(transcript,hdp=validator.hdp)
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
        exons = variant.map_dat.mapped_exons(
                transcript, chromosome_reference,
                alt_aln_method=validator.alt_aln_method,
                hdp=validator.hdp)

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

WARNING_CODE_MAP = {

    # ========================================================================
    # Syntax
    # ========================================================================

    "Unable to identify a colon":
        ("VariantSyntaxError", None),

    "Unable to identify a dot":
        ("VariantSyntaxError", None),

    "Removing redundant gene symbol":
        ("VariantSyntaxError", None),

    "Removing redundant reference bases":
        ("VariantSyntaxError", None),

    "This not a valid HGVS description, due to characters being in the wrong case":
        ("InvalidCaseError",
         "This is not a valid HGVS description because characters are in the wrong case. "
         "Please check the use of upper- and lowercase characters."),

    "Trailing digits are not permitted in":
        ("VariantSyntaxError", None),

    "Refer to http://varnomen.hgvs.org":
        ("VariantSyntaxError", None),

    # ========================================================================
    # Reference sequences / genome build
    # ========================================================================

    "no reference sequence ID has been provided":
        ("ReferenceSequenceError", None),

    "Reference sequence ":
        ("ReferenceSequenceError", None),

    "Reference type incorrectly stated":
        ("ReferenceTypeError", None),

    "invalid reference sequence identifier":
        ("ReferenceSequenceError", None),

    "is not a valid, and is also not a valid gene symbol":
        ("InvalidReferenceError", None),

    "HGVS variant nomenclature does not allow the use of a gene symbol":
        ("InvalidReferenceError", None),

    "A transcript reference sequence has not been provided":
        ("TranscriptReferenceError", None),

    "Multiple genomic reference sequences have been provided":
        ("ReferenceSequenceError", None),

    "apparent Transcript Reference ID was not recognised":
        ("TranscriptReferenceError", None),

    "is not part of genome build":
        ("GenomeReferenceWarning", None),

    "cannot be mapped directly to genome build":
        ("GenomeMismatchWarning", None),

    "Validation will fail if the selected chromosome reference sequence":
        ("GenomeMismatchWarning",
         "The selected chromosome reference sequence may not correspond to the selected genome build."),

    "did you mean GRCh":
        ("GenomeBuildWarning", None),

    "See alternative genomic loci":
        ("GenomeReferenceWarning", None),

    # ========================================================================
    # RNA
    # ========================================================================

    "RNA sequence must be lower-case":
        ("RnaAlphabetError", None),

    "RNA sequence contains Uracil":
        ("RnaAlphabetError",
         "RNA sequence contains thymine (T). RNA descriptions must use uracil (U)."),

    "The IUPAC RNA alphabet dictates":
        ("RnaAlphabetError", None),

    "The variant type for an RNA description must be r.":
        ("VariantSyntaxError", None),

    "Invalid variant type for non-coding transcript. Instead use n.":
        ("NonCodingTranscriptError", None),

    "Intronic descriptions are only valid in the context of a c. description":
        ("IntronicVariantError", None),

    # ========================================================================
    # Transcript
    # ========================================================================

    "not in our database":
        ("TranscriptMissingError", None),

    "versions of the requested transcript are available":
        ("TranscriptVersionWarning", None),

    "A more recent version of the selected reference sequence":
        ("TranscriptVersionWarning", None),

    "No individual transcripts have been identified":
        ("TranscriptIdentificationWarning", None),

    "Mapping unavailable for RefSeqGene":
        ("TranscriptMappingError", None),

    "Transcript ":
        ("TranscriptDataError", None),

    "Required information for ":
        ("TranscriptDataError", None),

    "Universal Transcript Archive":
        ("TranscriptDataError", None),

    "Query gene2transcripts":
        ("TranscriptSuggestionWarning", None),

    "None of the specified transcripts":
        ("TranscriptSelectionError", None),

    "Transcripts were found but the current transcript type limitation":
        ("TranscriptSelectionWarning", None),

    # ========================================================================
    # Mapping / normalisation
    # ========================================================================

    "automapped to genome position":
        ("GenomeMappingWarning", None),

    "automapped to":
        ("VariantMappingWarning", None),

    "normalized to":
        ("VariantNormalizationWarning", None),

    "mapped to":
        ("ReferenceMismatchWarning", None),

    "No relevant genomic mapping options":
        ("TranscriptMappingError", None),

    "Full alignment data between the specified transcript reference sequence":
        ("AlignmentDataWarning", None),

    "Alignment is incomplete":
        ("AlignmentDataWarning", None),

    "Suspected incomplete alignment between transcript":
        ("AlignmentDataWarning", None),

    "This coding sequence variant description spans at least one intron":
        ("IntronSpanningWarning", None),

    "Automap is unable to correct the input exon/intron boundary coordinates":
        ("ExonBoundaryError", None),

    # ========================================================================
    # Alleles
    # ========================================================================

    "The alleleic description is in the correct syntax":
        ("AlleleExtractionWarning",
         "The allelic description is syntactically correct and all possible variant descriptions have been extracted."),

    "Each variant is validated independently":
        ("AlleleValidationWarning", None),

    # ========================================================================
    # Protein
    # ========================================================================

    "Protein level variant descriptions are not fully supported":
        ("ProteinSupportWarning", None),

    "Cannot identify an in-frame Termination codon in the reference mRNA sequence":
        ("TranscriptTypeError", None),

    "Cannot identify an in-frame Termination codon in the variant mRNA sequence":
        ("ProteinTranslationWarning", None),

    "is HGVS compliant and contains a valid reference amino acid description":
        ("ProteinSupportWarning", None),

    "contains a valid reference amino acid description":
        ("ProteinTranslationInfo", None),

    "The amino acid at position":
        ("AminoMismatchError", None),

    "affects the initiation amino acid":
        ("InitiationCodonWarning", None),

    # ========================================================================
    # Insertions
    # ========================================================================

    "The inserted sequence must be provided":
        ("InsertionSequenceError", None),

    "Insertion length must be 1":
        ("InsertionLengthError", None),

    "An insertion must be provided with the two positions between which the insertion has taken place":
        ("InsertionLengthError", None),

    # ========================================================================
    # Expanded repeats
    # ========================================================================

    "The coordinates for the repeat region are stated incorrectly":
        ("ExpandedRepeatError", None),

    "should only be used as an annotation for the core HGVS descriptions provided":
        ("ExpandedRepeatWarning", None),

    # ========================================================================
    # VCF
    # ========================================================================

    "Conversions are no longer valid HGVS Sequence Variant Descriptions":
        ("LegacySyntaxError", None),

    "Insufficient or incorrect VCF elements provided":
        ("VcfFormatError", None),

    "Not stating ALT bases is ambiguous":
        ("AmbiguousVcfWarning", None),

    "VariantValidator has output both alternatives":
        ("AmbiguousVcfWarning", None),

    "CNV identified, and mapped to":
        ("VcfConversionWarning", None),

    "Multiple ALT sequences detected":
        ("MultipleAlleleWarning", None),

    # ========================================================================
    # LRG
    # ========================================================================

    "updated to equivalent RefSeq record":
        ("LrgMappingWarning", None),

    "updated to RefSeq record":
        ("LrgMappingWarning", None),

    "is pending therefore changes may be made":
        ("LrgStatusWarning", None),

    # ========================================================================
    # Mitochondrial
    # ========================================================================

    "is not associated with genome build hg19":
        ("MitochondrialBuildError", None),

    "is not associated with genome build GRCh37":
        ("MitochondrialBuildError", None),

    "does not match the DNA type (g).":
        ("MitochondrialReferenceError", None),

    # ========================================================================
    # HGVS normalizer
    # ========================================================================

    "Unsupported normalization of protein level variants":
        ("ProteinNormalizationError", None),

    "Unsupported normalization of conversion variants":
        ("ConversionNormalizationError", None),

    "Normalization of intronic variants is not supported":
        ("IntronicVariantError", None),

    "Unsupported normalization of variants spanning the exon-intron boundary":
        ("ExonBoundaryError", None),

    "No mapping info available for":
        ("TranscriptMappingError", None),

    "No identity info available for":
        ("TranscriptDataError", None),

    "Variant span is outside sequence bounds":
        ("OutOfBoundsError", None),

    # ========================================================================
    # HGVS validator
    # ========================================================================

    "Cannot validate sequence of an intronic variant":
        ("IntronicValidationWarning", None),

    "Variant reference (":
        ("ReferenceMismatchError", None),

    "does not agree with reference sequence":
        ("ReferenceMismatchError", None),

    "Variant coordinate is out of the bound of CDS region":
        ("CDSBoundaryError", None),

    "No transcript data for accession":
        ("TranscriptDataError", None),

    # ========================================================================
    # HGVS VariantMapper
    # ========================================================================

    "Expected a g. variant; got":
        ("VariantTypeError", None),

    "Expected a c. or n. variant; got":
        ("VariantTypeError", None),

    "Expected a cDNA (c.); got":
        ("VariantTypeError", None),

    "Expected a cDNA (c.) variant; got":
        ("VariantTypeError", None),

    "Expected n. variant; got":
        ("VariantTypeError", None),

    "Only NARefAlt/Dup/Inv types are currently implemented":
        ("UnsupportedEditError", None),

    "Can only update references for type c, g, m, n, r":
        ("ReferenceUpdateError", None),

    "Getting altered sequence for":
        ("SequenceGenerationError", None),

    # ========================================================================
    # HGVS AssemblyMapper
    # ========================================================================

    "Expected a coding (c.) or non-coding (n.) variant; got":
        ("VariantTypeError", None),

    "No alignments for":
        ("TranscriptMappingError", None),

    "non-pseudoautosomal region":
        ("MultipleAlignmentError", None),

    "likely pseudoautosomal region":
        ("PseudoautosomalRegionWarning", None),

    "in_par_assume=":
        ("PseudoautosomalRegionError", None),

    # ========================================================================
    # HGVS TranscriptMapper
    # ========================================================================

    "No transcript info":
        ("TranscriptDataError", None),

    "No transcript exons":
        ("TranscriptMappingError", None),

    "No transcript identity info":
        ("TranscriptDataError", None),

    "CDS start_i and end_i must be both defined or both undefined":
        ("TranscriptDataError", None),

    "CDS is undefined for":
        ("NonCodingTranscriptError", None),

    "The given coordinate is outside the bounds of the reference sequence":
        ("OutOfBoundsError", None),

    # ========================================================================
    # Internal
    # ========================================================================

    "If the following error message does not address the issue":
        ("InternalValidationError", None),

    # ========================================================================
    # Uncertain / fuzzy positions
    # ========================================================================

    "Uncertain positions are not fully supported, however the syntax is valid":
        ("UncertainPositionWarning", None),

    "Uncertain positions are not fully supported, however the start position is > the end position":
        ("InvalidRangeError", None),

    "Uncertain positions are not fully supported, however the provided positions are out of order":
        ("InvalidRangeError", None),

    "Selected transcript does not span the entire range of the genomic variation":
        ("TranscriptRangeWarning", None),

    "Only a single transcript can be processed, updating to select":
        ("TranscriptSelectionWarning", None),

    "Only a single transcript can be processed, updating to Select":
        ("TranscriptSelectionWarning", None),

    # ========================================================================
    # Fuzzy position exceptions
    # ========================================================================

    "Fuzzy/unknown variant start and end positions":
        ("FuzzyRangeError", None),

    "Fuzzy/unknown variant start position":
        ("FuzzyPositionError", None),

    "Fuzzy/unknown variant end position":
        ("FuzzyPositionError", None),

    "Invalid range submitted, missing underscore":
        ("InvalidRangeError", None),

    "is an invalid range for accession":
        ("InvalidRangeError", None),

    "Length implied by coordinates must equal":
        ("InvalidRangeError", None),

    "exon boundary ":
        ("ExonBoundaryError", None),

    "is not known to be compatible with variant type":
        ("IncompatibleTypeError", None),

    "base start position must be <= end position":
        ("IntervalOrderError", None),



    # ========================================================================
    # Use checking
    # ========================================================================

    "VariantValidator operates on variant descriptions, but":
        ("InvalidVariantError", None),

    "is a concatenation of":
        ("VariantSyntaxError", None),

    "HGVS descriptions contain a single colon":
        ("VariantSyntaxError", None),

    "Illegal addition of the invalid characters":
        ("VariantSyntaxError", None),

    "The format(s)":
        ("VariantSyntaxError", None),

    "lacks the . character between":
        ("VariantSyntaxError", None),

    "Stripping unnecessary characters":
        ("VariantSyntaxError", None),

    "is not in an accepted format":
        ("InvalidVariantError", None),

    "An insertion must be provided with the two positions":
        ("InsertionPositionError", None),

    "The length of the variant is not formatted following the HGVS guidelines":
        ("InsertionLengthError", None),

    "may also be written as":
        ("AlternativeRepresentationWarning", None),

    "Base substitution (>) submitted with a reference sequence range":
        ("VariantSyntaxError", None),

    "Transcript reference sequence input as genomic":
        ("ReferenceTypeError", None),

    "Non-coding transcript reference sequence input as coding":
        ("ReferenceTypeError", None),

    "Protein reference sequence input as Nucleotide":
        ("ReferenceTypeError", None),

    "Coding transcript reference sequence input as non-coding":
        ("ReferenceTypeError", None),

    "Using a nucleotide reference sequence (NM_ NR_ NG_ NC_)":
        ("ReferenceTypeError", None),

    "NG_:c.PositionVariation descriptions should not be used":
        ("ReferenceTypeError", None),

    "The variant positions are valid but we cannot normalize variants spanning the origin of circular reference sequences":
        ("CircularReferenceWarning", None),

    "Using a transcript reference sequence to specify a variant position that lies outside of the reference sequence":
        ("TranscriptCoordinateError", None),

    "Cannot map":
        ("TranscriptMappingError", None),

    "Interval start position":
        ("IntervalOrderError", None),

    "Interval end position":
        ("IntervalOrderError", None),

    "The given coordinate is outside the boundaries of reference sequence":
        ("OutOfBoundsError", None),

    "UncertainSequenceError:":
        ("UncertainSequenceError", None),

    # ========================================================================
    # Additional use checking
    # ========================================================================

    "contains only numeric characters":
        ("NumericInputError", None),

    "contains only alphanumeric characters":
        ("IncompleteVariantError", None),

    "auto-corrected":
        ("AutoCorrectionWarning", None),

    "updated to":
        ("AutoCorrectionWarning", None),

    "contains uncertainty":
        ("UncertainVariantWarning", None),

    "outside the reference sequence is not HGVS-compliant":
        ("TranscriptCoordinateError", None),

    "spanning the origin of circular reference sequences":
        ("CircularReferenceWarning", None),

    "auto-mapped to":
        ("VariantMappingWarning", None),

    # ========================================================================
    # Mixin converters
    # ========================================================================

    "No available t_to_g liftover":
        ("TranscriptMappingError", None),

    "no g. mapping options available":
        ("TranscriptMappingError", None),

    "Unable to identify a relevant transcript for":
        ("TranscriptIdentificationError", None),

    "AlleleVariantError:":
        ("AlleleVariantError", None),

    "AlleleSyntaxError:":
        ("AlleleSyntaxError", None),

    # ========================================================================
    # Expanded repeats
    # ========================================================================

    "RepeatSyntaxError:":
        ("RepeatSyntaxError", None),

    "ExonBoundaryError:":
        ("ExonBoundaryError", None),

    # ========================================================================
    # Overlapping position
    # ========================================================================
    "is > or overlaps":
        ("OverlappingPositionError", None),

    }


# Compile once at import
_ALREADY_CODED = re.compile(r"^[A-Z][A-Za-z0-9]+(?:Error|Warning|Info): ")

# Longest strings first
_sorted = sorted(
    WARNING_CODE_MAP.items(),
    key=lambda x: len(x[0]),
    reverse=True,
)

# Heuristic:
# Warnings beginning with an uppercase letter are assumed to match the
# start of the warning string. Warnings beginning with a lowercase letter
# are assumed to occur within the warning text.
_PREFIX_LOOKUPS = tuple(
    item for item in _sorted
    if item[0][0].isupper()
)

_SUBSTRING_LOOKUPS = tuple(
    item for item in _sorted
    if item[0][0].islower()
)


def normalise_warning_codes(warnings):
    """Apply standard VV warning/error codes.

    As a compromise for performance, warning lookups are pre-split at import
    into prefix and substring searches. This reduces the number of substring
    (`in`) comparisons during normal operation while avoiding the maintenance
    overhead of manually maintaining two lookup tables.
    """

    output = []

    for warning in warnings:
        warning = str(warning)

        # Already coded
        if _ALREADY_CODED.match(warning):
            output.append(warning)
            continue

        # Fast prefix lookups
        for search, (code, replacement) in _PREFIX_LOOKUPS:
            if warning.startswith(search):
                if replacement is None:
                    output.append(f"{code}: {warning}")
                else:
                    output.append(f"{code}: {replacement}")
                break
        else:
            # Slower substring lookups
            for search, (code, replacement) in _SUBSTRING_LOOKUPS:
                if search in warning:
                    if replacement is None:
                        output.append(f"{code}: {warning}")
                    else:
                        output.append(f"{code}: {replacement}")
                    break
            else:
                output.append(warning)

    return output


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
