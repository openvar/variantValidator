# -*- coding: utf-8 -*-
"""
functions.py

Module containing VariantValidator sub-functions. The majoirty of these functions require
hgvs Python package top-level functions or sub-functions contained in uta.py and
seqfetcher.py
"""

# IMPORT REQUIRED PYTHON MODULES
import re
import os
import sys
import copy
from vvLogging import logger

# Setup functions

# Config Section Mapping function
def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                logger.warning("skip: %s" % option)
        except:
            logger.warning("exception on %s!" % option)
            dict1[option] = None
    return dict1


# Set up paths
# FUNCTIONS_ROOT = os.path.dirname(os.path.abspath(__file__))
ENTREZ_ID = os.environ.get('ENTREZ_ID')
if ENTREZ_ID is None:
    from configparser import ConfigParser

    CONF_ROOT = os.environ.get('CONF_ROOT')
    Config = ConfigParser()
    Config.read(os.path.join(CONF_ROOT, 'config.ini'))
    ENTREZ_ID = ConfigSectionMap("EntrezID")['entrezid']

# IMPORT HGVS MODULES and create instances
import hgvs
import hgvs.exceptions
import hgvs.sequencevariant

# Error types
from hgvs.exceptions import HGVSError, HGVSDataNotAvailableError, HGVSUnsupportedOperationError


class mergeHGVSerror(Exception):
    pass


class alleleVariantError(Exception):
    pass

# # Connect to UTA
# hdp = hgvs.dataproviders.uta.connect(pooling=True)
# # Create normalizer
# hn = hgvs.normalizer.Normalizer(hdp,
#                                 cross_boundaries=False,
#                                 shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
#                                 alt_aln_method='splign'
#                                 )
# reverse_hn = hgvs.normalizer.Normalizer(hdp,
#                                         cross_boundaries=False,
#                                         shuffle_direction=5,
#                                         alt_aln_method='splign'
#                                         )
#
# # Create normalizer
# merge_normalizer = hgvs.normalizer.Normalizer(hdp,
#                                               cross_boundaries=False,
#                                               shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
#                                               alt_aln_method='splign',
#                                               validate=False
#                                               )
# reverse_merge_normalizer = hgvs.normalizer.Normalizer(hdp,
#                                                       cross_boundaries=False,
#                                                       shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
#                                                       alt_aln_method='splign',
#                                                       validate=False
#                                                       )
#
# # Validator
# vr = hgvs.validator.Validator(hdp)
# # parser
# hp = hgvs.parser.Parser()
# # Variantmapper
# vm = hgvs.variantmapper.VariantMapper(hdp, replace_reference=True)  # , normalize=False)
# nr_vm = hgvs.variantmapper.VariantMapper(hdp, replace_reference=False)
# # SeqFetcher
# sf = hgvs.dataproviders.seqfetcher.SeqFetcher()
#
# #create no_norm_evm
# no_norm_evm_38 = hgvs.assemblymapper.AssemblyMapper(hdp,
#                                                     assembly_name='GRCh38',
#                                                     alt_aln_method='splign',
#                                                     normalize=False,
#                                                     replace_reference=True
#                                                     )
#
# no_norm_evm_37 = hgvs.assemblymapper.AssemblyMapper(hdp,
#                                                     assembly_name='GRCh37',
#                                                     alt_aln_method='splign',
#                                                     normalize=False,
#                                                     replace_reference=True
#                                                     )

# variantanalyser modules
import dbControls
import supported_chromosome_builds
import hgvs2vcf
import pseudo_vcf2hgvs
import gap_genes
import links

# BioPython
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq

# HGNC rest variables
import httplib2 as http
import json

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

"""
usr_input
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
    pat_g = re.compile(":g\.")  # Pattern looks for :g.
    pat_gene = re.compile('\(.+?\)')  # Pattern looks for (....)
    pat_c = re.compile(":c\.")  # Pattern looks for :c.
    pat_r = re.compile(":r\.")  # Pattern looks for :r.
    pat_n = re.compile(":n\.")  # Pattern looks for :n.
    pat_p = re.compile(":p\.")  # Pattern looks for :p.
    pat_m = re.compile(":m\.")  # Pattern looks for :m.
    pat_est = re.compile("\d:\d")  # Pattern looks for number:number

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


"""
r_to_c
parses r. variant strings into hgvs object and maps to the c. equivalent. 

Marked for removal
"""

# def r_to_c(variant, evm, hp):
#     # convert the input string into a hgvs object by parsing
#     var_r = hp.parse_hgvs_variant(variant)
#     # map to the coding sequence
#     var_c = evm.r_to_c(var_r)  # coding level variant
#     variant = str(var_c)
#     c_from_r = {'variant': variant, 'type': ':c.'}
#     return c_from_r


""" 
Maps transcript variant descriptions onto specified RefSeqGene reference sequences
Return an hgvs object containing the genomic sequence variant relative to the RefSeqGene 
acession
refseq_ac = RefSeqGene ac

Marked for removal
"""


# def refseq(variant, vm, refseq_ac, hp, hdp, no_norm_evm, primary_assembly, vr, sf, nr_vm, hn):
#     # parse the variant into hgvs object
#     var_c = hp.parse_hgvs_variant(variant)
#     # map to the genomic co-ordinates using the easy variant mapper set to alt_aln_method = alt_aln_method
#     var_g = myevm_t_to_g(var_c, hdp, no_norm_evm, primary_assembly, vm, hp, hn, sf, nr_vm)
#     # Get overlapping transcripts - forcing a splign alignment
#     start_i = var_g.posedit.pos.start.base
#     end_i = var_g.posedit.pos.end.base
#     alt_ac = var_g.ac
#     alt_aln_method = 'splign'
#     transcripts = hdp.get_tx_for_region(alt_ac, alt_aln_method, start_i - 1, end_i)
#     # Take the first transcript
#     for trans in transcripts:
#         tx_ac = trans[0]
#         try:
#             ref_c = vm.g_to_t(var_g, tx_ac, alt_aln_method='splign')
#         except:
#             continue
#         else:
#             # map the variant co-ordinates to the refseq Gene accession using vm
#             ref_g_dict = {
#                 'ref_g': '',
#                 'error': 'false'
#             }
#             try:
#                 ref_g_dict['ref_g'] = vm.t_to_g(ref_c, alt_ac=refseq_ac, alt_aln_method='splign')
#             except:
#                 e = sys.exc_info()[0]
#                 ref_g_dict['error'] = e
#             try:
#                 vr.validate(ref_g_dict['ref_g'])
#             except:
#                 e = sys.exc_info()[0]
#                 ref_g_dict['error'] = e
#             if ref_g_dict['error'] == 'false':
#                 return ref_g_dict
#             else:
#                 continue
#     # Return as an error if all fail
#     return ref_g_dict


"""
Parses genomic variant strings into hgvs objects
Maps genomic hgvs object into a coding hgvs object if the c accession string is provided
returns a c. variant description string

Marked for removal
"""

# def g_to_c(var_g, tx_ac, hp, evm):
#     pat_g = re.compile(":g\.")  # Pattern looks for :g.
#     # If the :g. pattern is present in the input variant
#     if pat_g.search(var_g):
#         # convert the input string into a hgvs object by parsing
#         var_g = hp.parse_hgvs_variant(var_g)
#         # Map to coding variant
#         var_c = str(evm.g_to_c(var_g, tx_ac))
#         return var_c


"""
Parses genomic variant strings into hgvs objects
Maps genomic hgvs object into a non-coding hgvs object if the n accession string is provided
returns a n. variant description string

Marked for removal
"""

# def g_to_n(var_g, tx_ac, hp, evm):
#     pat_g = re.compile(":g\.")  # Pattern looks for :g.
#     # If the :g. pattern is present in the input variant
#     if pat_g.search(var_g):
#         # convert the input string into a hgvs object by parsing
#         var_g = hp.parse_hgvs_variant(var_g)
#         # Map to coding variant
#         var_n = str(evm.g_to_n(var_g, tx_ac))
#         return var_n


"""
Ensures variant strings are transcript c. or n.
returns parsed hgvs c. or n. object
"""


def coding(variant, hp):
    # If the :c. pattern is present in the input variant
    if re.search(':c.', variant) or re.search(':n.', variant):
        # convert the input string into a hgvs object
        var_c = hp.parse_hgvs_variant(variant)
        return var_c


"""
Mapping transcript to genomic position from a HGVS string rather than an hgvs (py) parsed object 
Interfaces with myevm t_to_g
Ensures variant strings are transcript c. or n.
returns parsed hgvs g. object
"""


def genomic(variant, no_norm_evm, hp, hdp, primary_assembly, vm, hn, sf, nr_vm):
    # Set regular expressions for if statements
    pat_g = re.compile(":g\.")  # Pattern looks for :g.
    pat_n = re.compile(":n\.")
    pat_c = re.compile(":c\.")  # Pattern looks for :c.

    # If the :c. pattern is present in the input variant
    if pat_c.search(variant) or pat_n.search(variant):
        error = 'false'
        hgvs_var = hp.parse_hgvs_variant(variant)
        try:
            var_g = myevm_t_to_g(hgvs_var, hdp, no_norm_evm, primary_assembly, vm, hp, hn, sf, nr_vm)
        except hgvs.exceptions.HGVSError as e:
            error = e
        if error != 'false':
            var_g = 'error ' + str(e)
        return var_g

    # If the :g. pattern is present in the input variant
    elif (pat_g.search(variant)):  # or (pat_n.search(variant)):
        # convert the input string into a hgvs object
        var_g = hp.parse_hgvs_variant(variant)
        return var_g


"""


Mapping transcript to protein prediction
Accepts a variant string rather than a parsed hgvs_object
Ensures variant strings are transcript c.
returns parsed hgvs p. object

Replaced by myc_to_p and marked for removal
"""

# def protein(variant, evm, hp, hdp):
#     # Set regular expressions for if statements
#     pat_c = re.compile(":c\.")  # Pattern looks for :c. Note (gene) has been removed
#
#     # If the :c. pattern is present in the input variant
#     if pat_c.search(variant):
#         # convert the input string into a hgvs object
#         var_c = hp.parse_hgvs_variant(variant)
#         # Does the edit affect the start codon?
#         if ((var_c.posedit.pos.start.base >= 1 and var_c.posedit.pos.start.base <= 3 and var_c.posedit.pos.start.offset == 0) or (
#                 var_c.posedit.pos.end.base >= 1 and var_c.posedit.pos.end.base <= 3 and var_c.posedit.pos.end.offset == 0)) and not re.search('\*', str(
#                 var_c.posedit.pos)):
#             ass_prot = hdp.get_pro_ac_for_tx_ac(var_c.ac)
#             if str(ass_prot) == 'None':
#                 cod = str(var_c)
#                 cod = cod.replace('inv', 'del')
#                 cod = hp.parse_hgvs_variant(cod)
#                 p = evm.c_to_p(cod)
#                 ass_prot = p.ac
#             var_p = hgvs.sequencevariant.SequenceVariant(ac=ass_prot, type='p', posedit='(Met1?)')
#         else:
#             var_p = evm.c_to_p(var_c)
#         return var_p
#     if re.search(':n.', variant):
#         var_p = hp.parse_hgvs_variant(variant)
#         var_p.ac = 'Non-coding transcript'
#         var_p.posedit = ''
#         return var_p

"""
Function which takes a NORMALIZED hgvs Python transcript variant and maps to a specified protein reference sequence. A protein
level hgvs python object is returned.

Note the function currently assumes that the transcript description is correctly normalized having come from the 
previous g_to_t function
"""


def myc_to_p(hgvs_transcript, evm, hdp, hp, hn, vm, sf, re_to_p):
    # Create dictionary to store the information
    hgvs_transcript_to_hgvs_protein = {'error': '', 'hgvs_protein': '', 'ref_residues': ''}

    # Collect the associated protein
    if hgvs_transcript.type == 'c':
        associated_protein_accession = hdp.get_pro_ac_for_tx_ac(hgvs_transcript.ac)
        # This method sometimes fails
        if str(associated_protein_accession) == 'None':
            cod = str(hgvs_transcript)
            cod = cod.replace('inv', 'del')
            cod = hp.parse_hgvs_variant(cod)
            p = evm.c_to_p(cod)
            associated_protein_accession = p.ac
    else:
        pass

        # Check for non-coding transcripts
    if hgvs_transcript.type == 'c':
        # Handle non inversions with simple c_to_p mapping

        if (hgvs_transcript.posedit.edit.type != 'inv') and (hgvs_transcript.posedit.edit.type != 'delins') and (
                re_to_p is False):
            # Does the edit affect the start codon?
            if ((
                        hgvs_transcript.posedit.pos.start.base >= 1 and hgvs_transcript.posedit.pos.start.base <= 3 and hgvs_transcript.posedit.pos.start.offset == 0) or (
                        hgvs_transcript.posedit.pos.end.base >= 1 and hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset == 0)) \
                    and not re.search('\*', str(
                hgvs_transcript.posedit.pos)):
                residue_one = sf.fetch_seq(associated_protein_accession, start_i=1-1,end_i=1)
                threed_residue_one = links.one_to_three(residue_one)
                r_one_report = '(%s1?)' % threed_residue_one
                hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                    type='p', posedit=r_one_report)
            else:
                try:
                    hgvs_protein = evm.c_to_p(hgvs_transcript)
                except IndexError as e:
                    error = str(e)
                    if re.search('string index out of range', error) and re.search('dup', str(hgvs_transcript)):
                        hgvs_ins = hp.parse_hgvs_variant(str(hgvs_transcript))
                        hgvs_ins = hn.normalize(hgvs_ins)
                        inst = hgvs_ins.ac + ':c.' + str(hgvs_ins.posedit.pos.start.base - 1) + '_' + str(
                            hgvs_ins.posedit.pos.start.base) + 'ins' + hgvs_ins.posedit.edit.ref
                        hgvs_transcript = hp.parse_hgvs_variant(inst)
                        hgvs_protein = evm.c_to_p(hgvs_transcript)

            try:
                hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                return hgvs_transcript_to_hgvs_protein
            except UnboundLocalError:
                hgvs_transcript_to_hgvs_protein = myc_to_p(hgvs_transcript, evm, hdp, hp, hn, vm, sf, re_to_p=True)
                return hgvs_transcript_to_hgvs_protein

        else:
            # Additional code required to process inversions
            # Note, this code was developed for VariantValidator and is not native to the biocommons hgvs Python package
            # Convert positions to n. position
            hgvs_naughty = vm.c_to_n(hgvs_transcript)

            # Collect the deleted sequence using fetch_seq
            del_seq = sf.fetch_seq(str(hgvs_naughty.ac), start_i=hgvs_naughty.posedit.pos.start.base - 1,
                                   end_i=hgvs_naughty.posedit.pos.end.base)

            # Make the inverted sequence
            my_seq = Seq(del_seq)

            if hgvs_transcript.posedit.edit.type == 'inv':
                inv_seq = my_seq.reverse_complement()
            else:
                inv_seq = hgvs_transcript.posedit.edit.alt
                if inv_seq is None:
                    inv_seq = ''

            # Look for p. delins or del
            not_delins = True
            if hgvs_transcript.posedit.edit.type != 'inv':
                try:
                    shifts = evm.c_to_p(hgvs_transcript)
                    if re.search('del', shifts.posedit.edit.type):
                        not_delins = False
                except Exception:
                    not_delins = False
            else:
                not_delins = False

            # Use inv delins code?
            if not_delins == False:
                # Collect the associated protein
                associated_protein_accession = hdp.get_pro_ac_for_tx_ac(hgvs_transcript.ac)

                # Intronic inversions are marked as uncertain i.e. p.?
                if re.search('\d+\-', str(hgvs_transcript.posedit.pos)) or re.search('\d+\+', str(
                        hgvs_transcript.posedit.pos)) or re.search('\*', str(hgvs_transcript.posedit.pos)) or re.search(
                        '[cn].\-', str(hgvs_transcript)):
                    if ((
                                hgvs_transcript.posedit.pos.start.base >= 1 and hgvs_transcript.posedit.pos.start.base <= 3 and hgvs_transcript.posedit.pos.start.offset == 0)
                        or
                        (
                                hgvs_transcript.posedit.pos.end.base >= 1 and hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset == 0)) \
                            and not re.search('\*', str(hgvs_transcript.posedit.pos)):
                        residue_one = sf.fetch_seq(associated_protein_accession, start_i=1 - 1, end_i=1)
                        threed_residue_one = links.one_to_three(residue_one)
                        r_one_report = '(%s1?)' % threed_residue_one # was (MET1?)
                        hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                            type='p', posedit=r_one_report)
                    else:
                        # Make the variant
                        hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p',
                                                                            posedit='?')
                    hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                    return hgvs_transcript_to_hgvs_protein
                else:
                    # Need to obtain the cds_start
                    inf = hdp.get_tx_identity_info(hgvs_transcript.ac)
                    cds_start = inf[3]

                    # Extract the reference coding sequence from SeqRepo
                    try:
                        ref_seq = sf.fetch_seq(str(hgvs_naughty.ac))
                    except Exception as e:
                        error = str(e)
                        hgvs_transcript_to_hgvs_protein['error'] = error
                        return hgvs_transcript_to_hgvs_protein

                    # Create the variant coding sequence
                    var_seq = links.n_inversion(ref_seq, del_seq, inv_seq,
                                                hgvs_naughty.posedit.pos.start.base,
                                                hgvs_naughty.posedit.pos.end.base)
                    # Translate the reference and variant proteins
                    prot_ref_seq = links.translate(ref_seq, cds_start)

                    try:
                        prot_var_seq = links.translate(var_seq, cds_start)
                    except IndexError:
                        hgvs_transcript_to_hgvs_protein[
                            'error'] = 'Cannot identify an in-frame Termination codon in the variant mRNA sequence'
                        hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p',
                                                                            posedit='?')
                        hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                        return hgvs_transcript_to_hgvs_protein

                    if prot_ref_seq == 'error':
                        error = 'Unable to generate protein variant description'
                        hgvs_transcript_to_hgvs_protein['error'] = error
                        return hgvs_transcript_to_hgvs_protein
                    elif prot_var_seq == 'error':
                        # Does the edit affect the start codon?
                        if ((
                                    hgvs_transcript.posedit.pos.start.base >= 1 and hgvs_transcript.posedit.pos.start.base <= 3 and hgvs_transcript.posedit.pos.start.offset == 0)
                            or
                            (
                                    hgvs_transcript.posedit.pos.end.base >= 1 and hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset == 0)) \
                                and not re.search('\*', str(hgvs_transcript.posedit.pos)):
                            residue_one = sf.fetch_seq(associated_protein_accession, start_i=1 - 1, end_i=1)
                            threed_residue_one = links.one_to_three(residue_one)
                            r_one_report = '(%s1?)' % threed_residue_one # was (MET1?)
                            hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                                type='p', posedit=r_one_report)

                            hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                            return hgvs_transcript_to_hgvs_protein
                        else:
                            error = 'Unable to generate protein variant description'
                            hgvs_transcript_to_hgvs_protein['error'] = error
                            return hgvs_transcript_to_hgvs_protein
                    else:
                        # Gather the required information regarding variant interval and sequences
                        if hgvs_transcript.posedit.edit.type != 'delins':
                            pro_inv_info = links.pro_inv_info(prot_ref_seq, prot_var_seq)
                        else:
                            pro_inv_info = links.pro_delins_info(prot_ref_seq, prot_var_seq)

                        # Error has occurred
                        if pro_inv_info['error'] == 'true':
                            error = 'Translation error occurred, please contact admin'
                            hgvs_transcript_to_hgvs_protein['error'] = error
                            return hgvs_transcript_to_hgvs_protein

                        # The Nucleotide variant has not affected the protein sequence i.e. synonymous
                        elif pro_inv_info['variant'] != 'true':
                            # Make the variant
                            hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                                type='p',
                                                                                posedit='=')
                            hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                            return hgvs_transcript_to_hgvs_protein

                        else:
                            # Early termination i.e. stop gained
                            # if pro_inv_info['terminate'] == 'true':
                            #     end = 'Ter' + str(pro_inv_info['ter_pos'])
                            #     pro_inv_info['prot_ins_seq'].replace('*', end)

                            # Complete variant description
                            # Recode the single letter del and ins sequences into three letter amino acid codes
                            del_thr = links.one_to_three(pro_inv_info['prot_del_seq'])
                            ins_thr = links.one_to_three(pro_inv_info['prot_ins_seq'])

                            # Write the HGVS position and edit
                            del_len = len(del_thr)
                            from_aa = del_thr[0:3]
                            to_aa = del_thr[del_len - 3:]

                            # Handle a range of amino acids
                            if pro_inv_info['edit_start'] != pro_inv_info['edit_end']:
                                if len(ins_thr) > 0:
                                    if re.search('Ter', del_thr) and ins_thr[-3:] != 'Ter':
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + str(
                                            pro_inv_info['edit_end']) + 'delins' + ins_thr + '?)'
                                    else:
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + str(
                                            pro_inv_info['edit_end']) + 'delins' + ins_thr + ')'
                                else:
                                    if re.search('Ter', del_thr) and ins_thr[-3:] != 'Ter':
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + str(
                                            pro_inv_info['edit_end']) + 'del?)'
                                    else:
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + str(
                                            pro_inv_info['edit_end']) + 'del)'
                            else:
                                # Handle extended proteins i.e. stop_lost
                                if del_thr == 'Ter' and (len(ins_thr) > len(del_thr)):
                                    # Nucleotide variant range aligns to the Termination codon
                                    if ins_thr[-3:] == 'Ter':
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + str(
                                            ins_thr[:3]) + 'ext' + str(ins_thr[-3:]) + str((len(ins_thr) / 3) - 1) + ')'
                                    # Nucleotide variant range spans the Termination codon
                                    else:
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + str(
                                            ins_thr[:3]) + 'ext?)'

                                # Nucleotide variation has not affected the length of the protein thus substitution or del
                                else:
                                    if len(ins_thr) == 3:
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + ins_thr + ')'
                                    elif len(ins_thr) == 0:
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + 'del)'
                                    else:
                                        posedit = '(' + from_aa + str(
                                            pro_inv_info['edit_start']) + 'delins' + ins_thr + ')'

                            # Complete the variant
                            hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                                type='p',
                                                                                posedit=posedit)

                            hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein

            else:
                hgvs_transcript_to_hgvs_protein['hgvs_protein'] = shifts

            # Return
            return hgvs_transcript_to_hgvs_protein


    # Handle non-coding transcript and non transcript descriptions
    elif hgvs_transcript.type == 'n':
        # non-coding transcripts
        hgvs_protein = copy.deepcopy(hgvs_transcript)
        hgvs_protein.ac = 'Non-coding '
        hgvs_protein.posedit = ''
        hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
        return hgvs_transcript_to_hgvs_protein
    else:
        hgvs_transcript_to_hgvs_protein['error'] = 'Unable to map %s to %s' % (
            hgvs_transcript.ac, associated_protein_accession)
        return hgvs_transcript_to_hgvs_protein


"""
Marked for removal
"""
# Return an hgvs object containing the rna sequence variant
# def rna(variant, evm, hp):
#   Set regular expressions for if statements
#   pat_c = re.compile(":c\.")         # Pattern looks for :c. Note (gene) has been removed
#   If the :c. pattern is present in the input variant
#   if  pat_c.search(variant):
#       convert the input string into a hgvs object
#       var_c = hp.parse_hgvs_variant(variant)
#       map to the genomic sequence
#       var_r = evm.c_to_n(var_c)   # rna level variant
#       return var_r

"""
Marked for removal
"""
# def hgvs_rna(variant, hp):
#   # Set regular expressions for if statements
#   pat_r = re.compile(":n\.")         # Pattern looks for :n. Note (gene) has been removed
#   # If the :r. pattern is present in the input variant
#   if  pat_r.search(variant):
#       # convert the input string into a hgvs object
#       var_r = hp.parse_hgvs_variant(variant)
#       return var_r


"""
Ensures variant strings are g.
returns parsed hgvs g. object

Marked for removal
"""

# def hgvs_genomic(variant, hp):
#     # Set regular expressions for if statements
#     pat_g = re.compile(":g\.")  # Pattern looks for :g. Note (gene) has been removed
#     # If the :g. pattern is present in the input variant
#     if pat_g.search(variant):
#         # convert the input string into a hgvs object
#         var_g = hp.parse_hgvs_variant(variant)
#         return var_g


"""
Enhanced transcript to genome position mapping function using evm
Deals with mapping from transcript positions that do not exist in the genomic sequence
i.e. the stated position aligns to a genomic gap!
Trys to ensure that a genomic position is always returned even if the c. or n. transcript
will not map to the specified genome build primary assembly.
Deals with transcript mapping to several genomic assemblies
Order 
Map to a single NC_ for the specified genome build primary assembly
Map to a single NC_ for an alternate genome build primary assembly
Map to an NT_ from the specified genome build
Map to an NT_ from an alternative genome build
Map to an NW_ from the specified genome build
Map to an NW_ from an alternative genome buildRequires parsed c. or n. object
returns parsed hgvs g. object
"""


def myevm_t_to_g(hgvs_c, hdp, no_norm_evm, primary_assembly, vm, hp, hn, sf, nr_vm):

    # store the input
    stored_hgvs_c = copy.deepcopy(hgvs_c)
    expand_out = 'false'
    utilise_gap_code = True

    # Gap gene black list
    try:
        gene_symbol = dbControls.data.get_gene_symbol_from_transcriptID(hgvs_c.ac)
    except Exception:
        utilise_gap_code = False
    else:
        # If the gene symbol is not in the list, the value False will be returned
        utilise_gap_code = gap_genes.gap_black_list(gene_symbol)
    # Warn gap code in use
    logger.warning("gap_compensation_myevm = " + str(utilise_gap_code))

    if utilise_gap_code is True and (
            hgvs_c.posedit.edit.type == 'identity' or hgvs_c.posedit.edit.type == 'del' or hgvs_c.posedit.edit.type == 'delins' or hgvs_c.posedit.edit.type == 'dup' or hgvs_c.posedit.edit.type == 'sub' or hgvs_c.posedit.edit.type == 'ins' or hgvs_c.posedit.edit.type == 'inv'):

        # if NM_ need the n. position
        if re.match('NM_', str(hgvs_c.ac)):
            hgvs_c = no_norm_evm.c_to_n(hgvs_c)

        # Check for intronic
        try:
            hn.normalize(hgvs_c)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if re.search('intronic variant', error):
                pass
            elif re.search('Length implied by coordinates must equal sequence deletion length', error) and re.match(
                    'NR_', hgvs_c.ac):
                hgvs_c.posedit.pos.end.base = hgvs_c.posedit.pos.start.base + len(hgvs_c.posedit.edit.ref) - 1

        # Check again before continuing
        if re.search('\d+\+', str(hgvs_c.posedit.pos)) or re.search('\d+\-', str(hgvs_c.posedit.pos)) or re.search(
                '\*\d+\+', str(hgvs_c.posedit.pos)) or re.search('\*\d+\-', str(hgvs_c.posedit.pos)):
            pass

        else:
            try:
                # For non-intronic sequence
                hgvs_t = copy.deepcopy(hgvs_c)
                if hgvs_t.posedit.edit.type == 'inv':
                    inv_alt = revcomp(hgvs_t.posedit.edit.ref)
                    t_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(
                        hgvs_t.posedit.pos.end.base) + 'del' + hgvs_t.posedit.edit.ref + 'ins' + inv_alt
                    hgvs_t_delins = hp.parse_hgvs_variant(t_delins)
                    pre_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                            hgvs_t.posedit.pos.start.base - 1)
                    post_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
                                             hgvs_t.posedit.pos.end.base + 1)
                    hgvs_t.posedit.edit.ref = pre_base + hgvs_t.posedit.edit.ref + post_base
                    inv_alt = pre_base + inv_alt + post_base
                    hgvs_t.posedit.pos.start.base = hgvs_t.posedit.pos.start.base - 1
                    start = hgvs_t.posedit.pos.start.base
                    hgvs_t.posedit.pos.start.base = start + 1
                    hgvs_t.posedit.pos.end.base = hgvs_t.posedit.pos.end.base + 1
                    end = hgvs_t.posedit.pos.end.base
                    hgvs_t.posedit.pos.start.base = start
                    hgvs_t.posedit.pos.end.base = end
                    hgvs_str = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(start) + '_' + str(
                        end) + 'del' + hgvs_t.posedit.edit.ref + 'ins' + inv_alt
                    hgvs_t = hp.parse_hgvs_variant(hgvs_str)
                elif hgvs_c.posedit.edit.type == 'dup':
                    pre_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                            hgvs_t.posedit.pos.start.base - 1)
                    post_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
                                             hgvs_t.posedit.pos.end.base + 1)
                    alt = pre_base + hgvs_t.posedit.edit.ref + hgvs_t.posedit.edit.ref + post_base
                    ref = pre_base + hgvs_t.posedit.edit.ref + post_base
                    dup_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(
                        hgvs_t.posedit.pos.start.base - 1) + '_' + str(
                        (hgvs_t.posedit.pos.start.base + len(ref)) - 2) + 'del' + ref + 'ins' + alt
                    hgvs_t = hp.parse_hgvs_variant(dup_to_delins)
                elif hgvs_c.posedit.edit.type == 'ins':
                    ins_ref = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                           hgvs_t.posedit.pos.end.base + 1)
                    ins_alt = ins_ref[:2] + hgvs_t.posedit.edit.alt + ins_ref[-2:]
                    ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(
                        hgvs_t.posedit.pos.start.base - 1) + '_' + str(
                        hgvs_t.posedit.pos.end.base + 1) + 'del' + ins_ref + 'ins' + ins_alt
                    hgvs_t = hp.parse_hgvs_variant(ins_to_delins)
                else:
                    if str(hgvs_t.posedit.edit.alt) == 'None':
                        hgvs_t.posedit.edit.alt = ''
                    pre_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                            hgvs_t.posedit.pos.start.base - 1)
                    post_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
                                             hgvs_t.posedit.pos.end.base + 1)
                    hgvs_t.posedit.edit.ref = pre_base + hgvs_t.posedit.edit.ref + post_base
                    hgvs_t.posedit.edit.alt = pre_base + hgvs_t.posedit.edit.alt + post_base
                    hgvs_t.posedit.pos.start.base = hgvs_t.posedit.pos.start.base - 1
                    start = hgvs_t.posedit.pos.start.base
                    hgvs_t.posedit.pos.start.base = start + 1
                    hgvs_t.posedit.pos.end.base = hgvs_t.posedit.pos.end.base + 1
                    end = hgvs_t.posedit.pos.end.base
                    hgvs_t.posedit.pos.start.base = start
                    hgvs_t.posedit.pos.end.base = end
                    hgvs_str = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(start) + '_' + str(end) + str(
                        hgvs_t.posedit.edit)
                    hgvs_t = hp.parse_hgvs_variant(hgvs_str)
                hgvs_c = copy.deepcopy(hgvs_t)

                # Set expanded out test to true
                expand_out = 'true'

            except Exception:
                hgvs_c = hgvs_c

        if re.match('NM_', str(hgvs_c.ac)):
            try:
                hgvs_c = no_norm_evm.n_to_c(hgvs_c)
            except hgvs.exceptions.HGVSError as e:
                hgvs_c = copy.deepcopy(stored_hgvs_c)

        # Ensure the altered c. variant has not crossed intro exon boundaries
        hgvs_check_boundaries = copy.deepcopy(hgvs_c)
        try:
            h_variant = hn.normalize(hgvs_check_boundaries)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if re.search('spanning the exon-intron boundary', error):
                hgvs_c = copy.deepcopy(stored_hgvs_c)
        # Catch identity at the exon/intron boundary by trying to normalize ref only
        if hgvs_check_boundaries.posedit.edit.type == 'identity':
            reform_ident = str(hgvs_c).split(':')[0]
            reform_ident = reform_ident + ':' + stored_hgvs_c.type + '.' + str(hgvs_c.posedit.pos) + 'del' + str(
                hgvs_c.posedit.edit.ref)  # + 'ins' + str(hgvs_c.posedit.edit.alt)
            hgvs_reform_ident = hp.parse_hgvs_variant(reform_ident)
            try:
                hn.normalize(hgvs_reform_ident)
            except hgvs.exceptions.HGVSError as e:
                error = str(e)
                if re.search('spanning the exon-intron boundary', error) or re.search(
                        'Normalization of intronic variants', error):
                    hgvs_c = copy.deepcopy(stored_hgvs_c)
    try:
        hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)
        hn.normalize(hgvs_genomic)  # Check the validity of the mapping
        # This will fail on multiple refs for NC_
    except hgvs.exceptions.HGVSError as e:
        # Recover all available mapping options from UTA
        mapping_options = hdp.get_tx_mapping_options(hgvs_c.ac)

        if mapping_options == []:
            raise HGVSDataNotAvailableError(
                "No alignment data between the specified transcript reference sequence and any GRCh37 and GRCh38 genomic reference sequences (including alternate chromosome assemblies, patches and RefSeqGenes) are available.")

        # Capture errors from attempted mappings
        attempted_mapping_error = ''

        for option in mapping_options:
            if re.match('blat', option[2]):
                continue
            if re.match('NC_', option[1]):
                chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
                if chr_num != 'false':
                    try:
                        hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                        break
                    except Exception as e:
                        attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + option[
                            1] + '~'
                        print e
                        continue

        # If not mapped, raise error
        try:
            hn.normalize(hgvs_genomic)
        except:
            for option in mapping_options:
                if re.match('blat', option[2]):
                    continue
                if re.match('NC_', option[1]):
                    chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
                    if chr_num == 'false':
                        try:
                            hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                            break
                        except Exception as e:
                            if re.search(option[1], attempted_mapping_error):
                                pass
                            else:
                                attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + option[
                                    1] + '~'
                            print e
                            continue
            try:
                hn.normalize(hgvs_genomic)
            except:
                for option in mapping_options:
                    if re.match('blat', option[2]):
                        continue
                    if re.match('NT_', option[1]):
                        chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
                        if chr_num != 'false':
                            try:
                                hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                break
                            except Exception as e:
                                attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + option[
                                    1] + '~'
                                print e
                                continue
                try:
                    hn.normalize(hgvs_genomic)
                except:
                    for option in mapping_options:
                        if re.match('blat', option[2]):
                            continue
                        if re.match('NT_', option[1]):
                            chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]),
                                                                                        primary_assembly)
                            if chr_num == 'false':
                                try:
                                    hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                    break
                                except Exception as e:
                                    attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + \
                                                              option[
                                                                  1] + '~'
                                    print e
                                    continue
                    try:
                        hn.normalize(hgvs_genomic)
                    except:
                        for option in mapping_options:
                            if re.match('blat', option[2]):
                                continue
                            if re.match('NW_', option[1]):
                                chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]),
                                                                                            primary_assembly)
                                if chr_num != 'false':
                                    try:
                                        hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                        break
                                    except Exception as e:
                                        attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + \
                                                                  option[1] + '~'
                                        print e
                                        continue
                        try:
                            hn.normalize(hgvs_genomic)
                        except:
                            for option in mapping_options:
                                if re.match('blat', option[2]):
                                    continue
                                if re.match('NW_', option[1]):
                                    chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]),
                                                                                                primary_assembly)
                                    if chr_num == 'false':
                                        try:
                                            hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                            break
                                        except Exception as e:
                                            attempted_mapping_error = attempted_mapping_error + str(
                                                e) + "/" + hgvs_c.ac + "/" + \
                                                                      option[1] + '~'
                                            print e
                                            continue

                            # Only a RefSeqGene available
                            try:
                                hn.normalize(hgvs_genomic)
                            except:
                                for option in mapping_options:
                                    if re.match('blat', option[2]):
                                        continue
                                    if re.match('NG_', option[1]):
                                        try:
                                            hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                            break
                                        except Exception as e:
                                            attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + \
                                                                      option[1] + '~'
                                            print e
                                            continue

    # If not mapped, raise error
    try:
        hgvs_genomic
    except Exception:
        raise HGVSDataNotAvailableError(attempted_mapping_error)

    if hgvs_c.posedit.edit.type == 'identity' and hgvs_genomic.posedit.edit.type == 'delins' and hgvs_genomic.posedit.edit.alt == '' and expand_out != 'true':
        hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.ref
    if hgvs_genomic.posedit.edit.type == 'ins' and utilise_gap_code is True:
        try:
            hgvs_genomic = hn.normalize(hgvs_genomic)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if error == 'insertion length must be 1':
                ref = sf.fetch_seq(str(hgvs_genomic.ac), hgvs_genomic.posedit.pos.start.base - 1,
                                   hgvs_genomic.posedit.pos.end.base)
                hgvs_genomic.posedit.edit.ref = ref
                hgvs_genomic.posedit.edit.alt = ref[0:1] + hgvs_genomic.posedit.edit.alt + ref[-1:]
                hgvs_genomic = hn.normalize(hgvs_genomic)
            if error == 'base start position must be <= end position':
                start = hgvs_genomic.posedit.pos.start.base
                end = hgvs_genomic.posedit.pos.end.base
                hgvs_genomic.posedit.pos.start.base = end
                hgvs_genomic.posedit.pos.end.base = start
                hgvs_genomic = hn.normalize(hgvs_genomic)

    # Statements required to reformat the stored_hgvs_c into a useable synonym
    if (stored_hgvs_c.posedit.edit.ref == '' or stored_hgvs_c.posedit.edit.ref is None) and expand_out != 'false':
        if stored_hgvs_c.type == 'c':
            stored_hgvs_n = vm.c_to_n(stored_hgvs_c)
        else:
            stored_hgvs_n = stored_hgvs_c
        stored_ref = sf.fetch_seq(str(stored_hgvs_n.ac), stored_hgvs_n.posedit.pos.start.base - 1,
                                  stored_hgvs_n.posedit.pos.end.base)
        stored_hgvs_c.posedit.edit.ref = stored_ref

    if (hgvs_genomic.posedit.edit.ref == '' or hgvs_genomic.posedit.edit.ref is None) and expand_out != 'false':
        if hgvs_genomic.posedit.edit.type == 'ins':
            stored_ref = sf.fetch_seq(str(hgvs_genomic.ac), hgvs_genomic.posedit.pos.start.base - 1,
                                      hgvs_genomic.posedit.pos.end.base)
            stored_alt = stored_ref[:1] + hgvs_genomic.posedit.edit.alt + stored_ref[-1:]
            hgvs_genomic.posedit.edit.ref = stored_ref
            hgvs_genomic.posedit.edit.alt = stored_alt

    # First look for variants mapping to the flanks of gaps
    # either in the gap or on the flank but not fully within the gap
    if expand_out == 'true':

        nr_genomic = nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)

        try:
            hn.normalize(nr_genomic)
        except hgvs.exceptions.HGVSInvalidVariantError as e:
            error_type_1 = str(e)
            if re.match('Length implied by coordinates must equal sequence deletion length', str(e)) or str(
                    e) == 'base start position must be <= end position':
                # Effectively, this code is designed to handle variants that are directly proximal to
                # gap BOUNDARIES, but in some cases the replace reference function of hgvs mapping has removed bases due to
                # the deletion length being > the specified range.

                # Warn of variant location wrt the gap
                if re.match('Length implied by coordinates must equal sequence deletion length', str(e)):
                    logger.warning('Variant is proximal to the flank of a genomic gap')
                    genomic_gap_variant = vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)
                    try:
                        hn.normalize(genomic_gap_variant)
                    # Still a problem
                    except hgvs.exceptions.HGVSInvalidVariantError as e:
                        if 'base start position must be <= end position' in str(e) and \
                                'Length implied by coordinates must equal' in error_type_1:
                            make_gen_var = copy.copy(nr_genomic)
                            make_gen_var.posedit.edit.ref = sf.fetch_seq(nr_genomic.ac,
                                                                         nr_genomic.posedit.pos.start.base - 1,
                                                                         nr_genomic.posedit.pos.end.base)
                            genomic_gap_variant = make_gen_var

                            error_type_1 = None
                    else:
                        genomic_gap_variant = nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)

                if error_type_1 == 'base start position must be <= end position':
                    logger.warning('Variant is fully within a genomic gap')
                    genomic_gap_variant = vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)

                # Logic
                # We have checked that the variant does not cross boundaries, or is intronic
                # So is likely mapping to a genomic gap
                try:
                    hn.normalize(genomic_gap_variant)
                except Exception as e:
                    if str(e) == 'base start position must be <= end position':
                        # This will only happen when the variant is fully within the gap
                        gap_start = genomic_gap_variant.posedit.pos.end.base
                        gap_end = genomic_gap_variant.posedit.pos.start.base
                        genomic_gap_variant.posedit.pos.start.base = gap_start
                        genomic_gap_variant.posedit.pos.end.base = gap_end
                    if re.match('Length implied by coordinates must equal sequence deletion length', str(e)):
                        # This will only happen if the variant is flanking the gap but is
                        # not inside the gap
                        logger.warning('Variant is on the flank of a genomic gap but not within the gap')
                        gap_start = genomic_gap_variant.posedit.pos.start.base - 1
                        gap_end = genomic_gap_variant.posedit.pos.end.base + 1
                        genomic_gap_variant.posedit.pos.start.base = gap_start
                        genomic_gap_variant.posedit.pos.end.base = gap_end
                        genomic_gap_variant.posedit.edit.ref = ''
                        stored_hgvs_c = copy.deepcopy(hgvs_c)

                    # Remove alt
                    try:
                        genomic_gap_variant.posedit.edit.alt = ''
                    except Exception as e:
                        if str(e) == "'Dup' object has no attribute 'alt'":
                            pass

                    # Should be a delins so will normalize statically and replace the reference bases
                    genomic_gap_variant = hn.normalize(genomic_gap_variant)
                    # Static map to c. and static normalize
                    transcript_gap_variant = vm.g_to_t(genomic_gap_variant, hgvs_c.ac)
                    stored_transcript_gap_variant = transcript_gap_variant

                    if not re.match('Length implied by coordinates must equal sequence deletion length', str(e)):
                        try:
                            transcript_gap_variant = hn.normalize(transcript_gap_variant)
                        except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                            if ' Unsupported normalization of variants spanning the UTR-exon boundary' in str(e):
                                pass

                    # if NM_ need the n. position
                    if re.match('NM_', str(hgvs_c.ac)):
                        transcript_gap_n = no_norm_evm.c_to_n(transcript_gap_variant)
                        transcript_gap_alt_n = no_norm_evm.c_to_n(stored_hgvs_c)
                    else:
                        transcript_gap_n = transcript_gap_variant
                        transcript_gap_alt_n = stored_hgvs_c

                    # Ensure an ALT exists
                    try:
                        if transcript_gap_alt_n.posedit.edit.alt is None:
                            transcript_gap_alt_n.posedit.edit.alt = 'X'
                    except Exception as e:
                        if str(e) == "'Dup' object has no attribute 'alt'":
                            transcript_gap_n_delins_from_dup = transcript_gap_n.ac + ':' + transcript_gap_n.type + '.' + str(
                                transcript_gap_n.posedit.pos.start.base) + '_' + str(
                                transcript_gap_n.posedit.pos.end.base) + 'del' + transcript_gap_n.posedit.edit.ref + 'ins' + transcript_gap_n.posedit.edit.ref + transcript_gap_n.posedit.edit.ref
                            transcript_gap_n = hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)
                            transcript_gap_alt_n_delins_from_dup = transcript_gap_alt_n.ac + ':' + transcript_gap_alt_n.type + '.' + str(
                                transcript_gap_alt_n.posedit.pos.start.base) + '_' + str(
                                transcript_gap_alt_n.posedit.pos.end.base) + 'del' + transcript_gap_alt_n.posedit.edit.ref + 'ins' + transcript_gap_alt_n.posedit.edit.ref + transcript_gap_alt_n.posedit.edit.ref
                            transcript_gap_alt_n = hp.parse_hgvs_variant(transcript_gap_alt_n_delins_from_dup)

                    # Split the reference and replacing alt sequence into a dictionary
                    reference_bases = list(transcript_gap_n.posedit.edit.ref)
                    if transcript_gap_alt_n.posedit.edit.alt is not None:
                        alternate_bases = list(transcript_gap_alt_n.posedit.edit.alt)
                    else:
                        # Deletions with no ins
                        pre_alternate_bases = list(transcript_gap_alt_n.posedit.edit.ref)
                        alternate_bases = []
                        for base in pre_alternate_bases:
                            alternate_bases.append('X')

                    # Create the dictionaries
                    ref_start = transcript_gap_n.posedit.pos.start.base
                    alt_start = transcript_gap_alt_n.posedit.pos.start.base
                    ref_base_dict = {}
                    for base in reference_bases:
                        ref_base_dict[ref_start] = str(base)
                        ref_start = ref_start + 1

                    alt_base_dict = {}

                    # Note, all variants will be forced into the format delete insert
                    # Deleted bases in the ALT will be substituted for X
                    for int in range(transcript_gap_alt_n.posedit.pos.start.base,
                                     transcript_gap_alt_n.posedit.pos.end.base + 1, 1):
                        if int == alt_start:
                            alt_base_dict[int] = str(''.join(alternate_bases))
                        else:
                            alt_base_dict[int] = 'X'

                            # Generate the alt sequence
                    alternate_sequence_bases = []
                    for int in range(transcript_gap_n.posedit.pos.start.base, transcript_gap_n.posedit.pos.end.base + 1,
                                     1):
                        if int in alt_base_dict.keys():
                            alternate_sequence_bases.append(alt_base_dict[int])
                        else:
                            alternate_sequence_bases.append(ref_base_dict[int])
                    alternate_sequence = ''.join(alternate_sequence_bases)
                    alternate_sequence = alternate_sequence.replace('X', '')

                    # Update variant, map to genome using vm and normalize
                    transcript_gap_n.posedit.edit.alt = alternate_sequence

                    try:
                        transcript_gap_variant = vm.n_to_c(transcript_gap_n)
                    except:
                        transcript_gap_variant = transcript_gap_n

                    try:
                        hgvs_genomic = vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                        hgvs_genomic = hn.normalize(hgvs_genomic)
                    except Exception as e:
                        if str(e) == "base start position must be <= end position":
                            # Expansion out is required to map back to the genomic position
                            pre_base = sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.start.base - 2,
                                                    transcript_gap_n.posedit.pos.start.base - 1)
                            post_base = sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.end.base,
                                                     transcript_gap_n.posedit.pos.end.base + 1)
                            transcript_gap_n.posedit.pos.start.base = transcript_gap_n.posedit.pos.start.base - 1
                            transcript_gap_n.posedit.pos.end.base = transcript_gap_n.posedit.pos.end.base + 1
                            transcript_gap_n.posedit.edit.ref = pre_base + transcript_gap_n.posedit.edit.ref + post_base
                            transcript_gap_n.posedit.edit.alt = pre_base + transcript_gap_n.posedit.edit.alt + post_base
                            try:
                                transcript_gap_variant = vm.n_to_c(transcript_gap_n)
                            except:
                                transcript_gap_variant = transcript_gap_n
                            hgvs_genomic = vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                            hgvs_genomic = hn.normalize(hgvs_genomic)

                    # Bypass the next bit of gap code
                    expand_out = 'false'

            else:
                pass
        # No map to the flank of a gap or within the gap
        else:
            pass

    # CASCADING STATEMENTS WHICH CAPTURE t to g MAPPING OPTIONS
    # Remove identity bases
    if hgvs_c == stored_hgvs_c:
        expand_out = 'false'
    elif expand_out == 'false' or utilise_gap_code is False:
        pass
    # Correct expansion ref + 2
    elif expand_out == 'true' and (
            len(hgvs_genomic.posedit.edit.ref) == (len(stored_hgvs_c.posedit.edit.ref) + 2)):  # >= 3:
        hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base + 1
        hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base - 1
        hgvs_genomic.posedit.edit.ref = hgvs_genomic.posedit.edit.ref[1:-1]
        if hgvs_genomic.posedit.edit.alt is not None:
            hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt[1:-1]
    elif expand_out == 'true' and (
            len(hgvs_genomic.posedit.edit.ref) != (len(stored_hgvs_c.posedit.edit.ref) + 2)):  # >= 3:
        if expand_out == 'true' and len(hgvs_genomic.posedit.edit.ref) == 2:
            gn = hn.normalize(hgvs_genomic)
            pass

        # Likely if the start or end position aligns to a gap in the genomic sequence
        # Logic
        # We have checked that the variant does not cross boundaries, or is intronic
        # So is likely mapping to a genomic gap
        elif expand_out == 'true' and len(hgvs_genomic.posedit.edit.ref) <= 1:
            # Incorrect expansion, likely < ref + 2
            genomic_gap_variant = vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)
            try:
                hn.normalize(genomic_gap_variant)
            except Exception as e:
                if str(e) == 'base start position must be <= end position':
                    gap_start = genomic_gap_variant.posedit.pos.end.base
                    gap_end = genomic_gap_variant.posedit.pos.start.base
                    genomic_gap_variant.posedit.pos.start.base = gap_start
                    genomic_gap_variant.posedit.pos.end.base = gap_end
                # Remove alt
                try:
                    genomic_gap_variant.posedit.edit.alt = ''
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        pass
                # Should be a delins so will normalize statically and replace the reference bases
                genomic_gap_variant = hn.normalize(genomic_gap_variant)
                # Static map to c. and static normalize
                transcript_gap_variant = vm.g_to_t(genomic_gap_variant, hgvs_c.ac)
                stored_transcript_gap_variant = transcript_gap_variant
                transcript_gap_variant = hn.normalize(transcript_gap_variant)
                # if NM_ need the n. position
                if re.match('NM_', str(hgvs_c.ac)):
                    transcript_gap_n = no_norm_evm.c_to_n(transcript_gap_variant)
                    transcript_gap_alt_n = no_norm_evm.c_to_n(stored_hgvs_c)
                else:
                    transcript_gap_n = transcript_gap_variant
                    transcript_gap_alt_n = stored_hgvs_c

                # Ensure an ALT exists
                try:
                    if transcript_gap_alt_n.posedit.edit.alt is None:
                        transcript_gap_alt_n.posedit.edit.alt = 'X'
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        transcript_gap_n_delins_from_dup = transcript_gap_n.ac + ':' + transcript_gap_n.type + '.' + str(
                            transcript_gap_n.posedit.pos.start.base) + '_' + str(
                            transcript_gap_n.posedit.pos.end.base) + 'del' + transcript_gap_n.posedit.edit.ref + 'ins' + transcript_gap_n.posedit.edit.ref + transcript_gap_n.posedit.edit.ref
                        transcript_gap_n = hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)
                        transcript_gap_alt_n_delins_from_dup = transcript_gap_alt_n.ac + ':' + transcript_gap_alt_n.type + '.' + str(
                            transcript_gap_alt_n.posedit.pos.start.base) + '_' + str(
                            transcript_gap_alt_n.posedit.pos.end.base) + 'del' + transcript_gap_alt_n.posedit.edit.ref + 'ins' + transcript_gap_alt_n.posedit.edit.ref + transcript_gap_alt_n.posedit.edit.ref
                        transcript_gap_alt_n = hp.parse_hgvs_variant(transcript_gap_alt_n_delins_from_dup)

                # Split the reference and replacing alt sequence into a dictionary
                reference_bases = list(transcript_gap_n.posedit.edit.ref)
                if transcript_gap_alt_n.posedit.edit.alt is not None:
                    alternate_bases = list(transcript_gap_alt_n.posedit.edit.alt)
                else:
                    # Deletions with no ins
                    pre_alternate_bases = list(transcript_gap_alt_n.posedit.edit.ref)
                    alternate_bases = []
                    for base in pre_alternate_bases:
                        alternate_bases.append('X')

                # Create the dictionaries
                ref_start = transcript_gap_n.posedit.pos.start.base
                alt_start = transcript_gap_alt_n.posedit.pos.start.base
                ref_base_dict = {}
                for base in reference_bases:
                    ref_base_dict[ref_start] = str(base)
                    ref_start = ref_start + 1

                alt_base_dict = {}

                # Note, all variants will be forced into the format delete insert
                # Deleted bases in the ALT will be substituted for X
                for int in range(transcript_gap_alt_n.posedit.pos.start.base,
                                 transcript_gap_alt_n.posedit.pos.end.base + 1, 1):
                    if int == alt_start:
                        alt_base_dict[int] = str(''.join(alternate_bases))
                    else:
                        alt_base_dict[int] = 'X'

                # Generate the alt sequence
                alternate_sequence_bases = []
                for int in range(transcript_gap_n.posedit.pos.start.base, transcript_gap_n.posedit.pos.end.base + 1, 1):
                    if int in alt_base_dict.keys():
                        alternate_sequence_bases.append(alt_base_dict[int])
                    else:
                        alternate_sequence_bases.append(ref_base_dict[int])
                alternate_sequence = ''.join(alternate_sequence_bases)
                alternate_sequence = alternate_sequence.replace('X', '')

                # Update variant, map to genome using vm and normalize
                transcript_gap_n.posedit.edit.alt = alternate_sequence

                try:
                    transcript_gap_variant = vm.n_to_c(transcript_gap_n)
                except:
                    transcript_gap_variant = transcript_gap_n

                try:
                    hgvs_genomic = vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                    hgvs_genomic = hn.normalize(hgvs_genomic)
                except Exception as e:
                    if str(e) == "base start position must be <= end position":
                        # Expansion out is required to map back to the genomic position
                        pre_base = sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.start.base - 2,
                                                transcript_gap_n.posedit.pos.start.base - 1)
                        post_base = sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.end.base,
                                                 transcript_gap_n.posedit.pos.end.base + 1)
                        transcript_gap_n.posedit.pos.start.base = transcript_gap_n.posedit.pos.start.base - 1
                        transcript_gap_n.posedit.pos.end.base = transcript_gap_n.posedit.pos.end.base + 1
                        transcript_gap_n.posedit.edit.ref = pre_base + transcript_gap_n.posedit.edit.ref + post_base
                        transcript_gap_n.posedit.edit.alt = pre_base + transcript_gap_n.posedit.edit.alt + post_base
                        try:
                            transcript_gap_variant = vm.n_to_c(transcript_gap_n)
                        except:
                            transcript_gap_variant = transcript_gap_n
                        hgvs_genomic = vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                        hgvs_genomic = hn.normalize(hgvs_genomic)

    # Ins variants map badly - Especially between c. exon/exon boundary
    if hgvs_c.posedit.edit.type == 'ins' and hgvs_c.posedit.pos.start.offset == 0 and hgvs_c.posedit.pos.end.offset == 0:
        try:
            hn.normalize(hgvs_genomic)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if error == 'insertion length must be 1':
                if hgvs_c.type == 'c':
                    hgvs_t = vm.c_to_n(hgvs_c)
                else:
                    hgvs_t = copy.copy(hgvs_c)
                ins_ref = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 1, hgvs_t.posedit.pos.end.base)
                ins_alt = ins_ref[:1] + hgvs_t.posedit.edit.alt + ins_ref[-1:]
                ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(
                    hgvs_t.posedit.pos.end.base) + 'del' + ins_ref + 'ins' + ins_alt
                hgvs_t = hp.parse_hgvs_variant(ins_to_delins)
                try:
                    hgvs_c = vm.n_to_c(hgvs_t)
                except Exception:
                    hgvs_c = copy.copy(hgvs_t)
                try:
                    hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)
                except Exception as e:
                    error = str(e)
                    logger.warning('Ins mapping error in myt_to_g ' + error)

    return hgvs_genomic


"""
USE WITH MAPPER THAT DOES NOT REPLACE THE REFERENCE GENOMIC BASES AND DOED NOT NORMALIZE

Enhanced transcript to genome position mapping function using evm
Trys to ensure that a genomic position is always returned even if the c. or n. transcript
will not map to the specified genome build primary assembly.
Deals with transcript mapping to several genomic assemblies
Order 
Map to a single NC_ (or ALT) for the specified genome build
returns parsed hgvs g. object
"""


def noreplace_myevm_t_to_g(hgvs_c, evm, hdp, primary_assembly, vm, hn, hp, sf, no_norm_evm):
    try:
        hgvs_genomic = evm.t_to_g(hgvs_c)
        hn.normalize(hgvs_genomic)
    # This will fail on multiple refs for NC_
    except hgvs.exceptions.HGVSError as e:
        # Recover all available mapping options from UTA
        mapping_options = hdp.get_tx_mapping_options(hgvs_c.ac)

        if mapping_options == []:
            raise HGVSDataNotAvailableError("no g. mapping options available")

        for option in mapping_options:
            if re.match('blat', option[2]):
                continue
            if re.match('NC_', option[1]):
                chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
                if chr_num != 'false':
                    try:
                        hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                        break
                    except Exception as e:
                        attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + option[
                            1] + '~'
                        print e
                        continue

        # If not mapped, raise error
        try:
            hn.normalize(hgvs_genomic)
        except:
            for option in mapping_options:
                if re.match('blat', option[2]):
                    continue
                if re.match('NC_', option[1]):
                    chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
                    if chr_num != 'false':
                        try:
                            hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                            break
                        except Exception as e:
                            attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + option[
                                1] + '~'
                            print e
                            continue

            # If not mapped, raise error
            try:
                hn.normalize(hgvs_genomic)
            except:
                for option in mapping_options:
                    if re.match('blat', option[2]):
                        continue
                    if re.match('NC_', option[1]):
                        chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
                        if chr_num == 'false':
                            try:
                                hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                break
                            except Exception as e:
                                if re.search(option[1], attempted_mapping_error):
                                    pass
                                else:
                                    attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + \
                                                              option[
                                                                  1] + '~'
                                print e
                                continue
                try:
                    hn.normalize(hgvs_genomic)
                except:
                    for option in mapping_options:
                        if re.match('blat', option[2]):
                            continue
                        if re.match('NT_', option[1]):
                            chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]),
                                                                                        primary_assembly)
                            if chr_num != 'false':
                                try:
                                    hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                    break
                                except Exception as e:
                                    attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + \
                                                              option[
                                                                  1] + '~'
                                    print e
                                    continue
                    try:
                        hn.normalize(hgvs_genomic)
                    except:
                        for option in mapping_options:
                            if re.match('blat', option[2]):
                                continue
                            if re.match('NT_', option[1]):
                                chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]),
                                                                                            primary_assembly)
                                if chr_num == 'false':
                                    try:
                                        hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                        break
                                    except Exception as e:
                                        attempted_mapping_error = attempted_mapping_error + str(
                                            e) + "/" + hgvs_c.ac + "/" + \
                                                                  option[
                                                                      1] + '~'
                                        print e
                                        continue
                        try:
                            hn.normalize(hgvs_genomic)
                        except:
                            for option in mapping_options:
                                if re.match('blat', option[2]):
                                    continue
                                if re.match('NW_', option[1]):
                                    chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]),
                                                                                                primary_assembly)
                                    if chr_num != 'false':
                                        try:
                                            hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                            break
                                        except Exception as e:
                                            attempted_mapping_error = attempted_mapping_error + str(
                                                e) + "/" + hgvs_c.ac + "/" + \
                                                                      option[1] + '~'
                                            print e
                                            continue
                            try:
                                hn.normalize(hgvs_genomic)
                            except:
                                for option in mapping_options:
                                    if re.match('blat', option[2]):
                                        continue
                                    if re.match('NW_', option[1]):
                                        chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]),
                                                                                                    primary_assembly)
                                        if chr_num == 'false':
                                            try:
                                                hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                                break
                                            except Exception as e:
                                                attempted_mapping_error = attempted_mapping_error + str(
                                                    e) + "/" + hgvs_c.ac + "/" + \
                                                                          option[1] + '~'
                                                print e
                                                continue

                                # Only a RefSeqGene available
                                try:
                                    hn.normalize(hgvs_genomic)
                                except:
                                    for option in mapping_options:
                                        if re.match('blat', option[2]):
                                            continue
                                        if re.match('NG_', option[1]):
                                            try:
                                                hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                                break
                                            except Exception as e:
                                                attempted_mapping_error = attempted_mapping_error + str(
                                                    e) + "/" + hgvs_c.ac + "/" + \
                                                                          option[1] + '~'
                                                print e
                                                continue
    try:
        hgvs_genomic
    except Exception:
        raise HGVSDataNotAvailableError('No available t_to_g liftover')

    # Ins variants map badly - Especially between c. exon/exon boundary
    if hgvs_c.posedit.edit.type == 'ins' and hgvs_c.posedit.pos.start.offset == 0 and hgvs_c.posedit.pos.end.offset == 0:
        try:
            hn.normalize(hgvs_genomic)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if error == 'insertion length must be 1':
                if hgvs_c.type == 'c':
                    hgvs_t = vm.c_to_n(hgvs_c)
                else:
                    hgvs_t = copy.copy(hgvs_c)
                ins_ref = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 1, hgvs_t.posedit.pos.end.base)
                ins_alt = ins_ref[:1] + hgvs_t.posedit.edit.alt + ins_ref[-1:]
                ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(
                    hgvs_t.posedit.pos.end.base) + 'del' + ins_ref + 'ins' + ins_alt
                hgvs_t = hp.parse_hgvs_variant(ins_to_delins)
                try:
                    hgvs_c = vm.n_to_c(hgvs_t)
                except Exception:
                    hgvs_c = copy.copy(hgvs_t)
                try:
                    hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)
                except Exception as e:
                    error = str(e)
                    logger.warning('Ins mapping error in myt_to_g ' + error)

    return hgvs_genomic


"""
Enhanced transcript to genome position on a specified genomic reference using vm
Deals with mapping from transcript positions that do not exist in the genomic sequence
i.e. the stated position aligns to a genomic gap!
returns parsed hgvs g. object
"""


def myvm_t_to_g(hgvs_c, alt_chr, no_norm_evm, vm, hp, hn, sf, nr_vm):
    # store the input
    stored_hgvs_c = copy.deepcopy(hgvs_c)
    expand_out = 'false'
    utilise_gap_code = True

    # Gap gene black list
    try:
        gene_symbol = dbControls.data.get_gene_symbol_from_transcriptID(hgvs_c.ac)
    except Exception:
        utilise_gap_code = False
    else:
        # If the gene symbol is not in the list, the value False will be returned
        utilise_gap_code = gap_genes.gap_black_list(gene_symbol)
    # Warn gap code in use
    logger.warning("gap_compensation_mvm = " + str(utilise_gap_code))

    if utilise_gap_code is True and (
            hgvs_c.posedit.edit.type == 'identity' or hgvs_c.posedit.edit.type == 'del' or hgvs_c.posedit.edit.type == 'delins' or hgvs_c.posedit.edit.type == 'dup' or hgvs_c.posedit.edit.type == 'sub' or hgvs_c.posedit.edit.type == 'ins' or hgvs_c.posedit.edit.type == 'inv'):

        # if NM_ need the n. position
        if re.match('NM_', str(hgvs_c.ac)):
            hgvs_c = no_norm_evm.c_to_n(hgvs_c)

        # Check for intronic
        try:
            hn.normalize(hgvs_c)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if re.search('intronic variant', error):
                pass
            elif re.search('Length implied by coordinates must equal sequence deletion length', error) and re.match(
                    'NR_', hgvs_c.ac):
                hgvs_c.posedit.pos.end.base = hgvs_c.posedit.pos.start.base + len(hgvs_c.posedit.edit.ref) - 1

        # Check again before continuing
        if re.search('\d+\+', str(hgvs_c.posedit.pos)) or re.search('\d+\-', str(hgvs_c.posedit.pos)) or re.search(
                '\*\d+\+', str(hgvs_c.posedit.pos)) or re.search('\*\d+\-', str(hgvs_c.posedit.pos)):
            pass

        else:
            try:
                # For non-intronic sequence
                hgvs_t = copy.deepcopy(hgvs_c)
                # handle inversions
                if hgvs_t.posedit.edit.type == 'inv':
                    inv_alt = revcomp(hgvs_t.posedit.edit.ref)
                    t_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(
                        hgvs_t.posedit.pos.end.base) + 'del' + hgvs_t.posedit.edit.ref + 'ins' + inv_alt
                    hgvs_t_delins = hp.parse_hgvs_variant(t_delins)
                    pre_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                            hgvs_t.posedit.pos.start.base - 1)
                    post_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
                                             hgvs_t.posedit.pos.end.base + 1)
                    hgvs_t.posedit.edit.ref = pre_base + hgvs_t.posedit.edit.ref + post_base
                    inv_alt = pre_base + inv_alt + post_base
                    hgvs_t.posedit.pos.start.base = hgvs_t.posedit.pos.start.base - 1
                    start = hgvs_t.posedit.pos.start.base
                    hgvs_t.posedit.pos.start.base = start + 1
                    hgvs_t.posedit.pos.end.base = hgvs_t.posedit.pos.end.base + 1
                    end = hgvs_t.posedit.pos.end.base
                    hgvs_t.posedit.pos.start.base = start
                    hgvs_t.posedit.pos.end.base = end
                    hgvs_str = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(start) + '_' + str(
                        end) + 'del' + hgvs_t.posedit.edit.ref + 'ins' + inv_alt
                    hgvs_t = hp.parse_hgvs_variant(hgvs_str)
                if hgvs_c.posedit.edit.type == 'dup':
                    # hgvs_t = reverse_normalize.normalize(hgvs_t)
                    pre_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                            hgvs_t.posedit.pos.start.base - 1)
                    post_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
                                             hgvs_t.posedit.pos.end.base + 1)
                    alt = pre_base + hgvs_t.posedit.edit.ref + hgvs_t.posedit.edit.ref + post_base
                    ref = pre_base + hgvs_t.posedit.edit.ref + post_base
                    dup_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(
                        hgvs_t.posedit.pos.start.base - 1) + '_' + str(
                        (hgvs_t.posedit.pos.start.base + len(ref)) - 2) + 'del' + ref + 'ins' + alt
                    hgvs_t = hp.parse_hgvs_variant(dup_to_delins)
                elif hgvs_c.posedit.edit.type == 'ins':
                    ins_ref = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                           hgvs_t.posedit.pos.end.base + 1)
                    ins_alt = ins_ref[:2] + hgvs_t.posedit.edit.alt + ins_ref[-2:]
                    ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(
                        hgvs_t.posedit.pos.start.base - 1) + '_' + str(
                        hgvs_t.posedit.pos.end.base + 1) + 'del' + ins_ref + 'ins' + ins_alt
                    hgvs_t = hp.parse_hgvs_variant(ins_to_delins)
                else:
                    if str(hgvs_t.posedit.edit.alt) == 'None':
                        hgvs_t.posedit.edit.alt = ''
                    pre_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                            hgvs_t.posedit.pos.start.base - 1)
                    post_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
                                             hgvs_t.posedit.pos.end.base + 1)
                    hgvs_t.posedit.edit.ref = pre_base + hgvs_t.posedit.edit.ref + post_base
                    hgvs_t.posedit.edit.alt = pre_base + hgvs_t.posedit.edit.alt + post_base
                    hgvs_t.posedit.pos.start.base = hgvs_t.posedit.pos.start.base - 1
                    start = hgvs_t.posedit.pos.start.base
                    hgvs_t.posedit.pos.start.base = start + 1
                    hgvs_t.posedit.pos.end.base = hgvs_t.posedit.pos.end.base + 1
                    end = hgvs_t.posedit.pos.end.base
                    hgvs_t.posedit.pos.start.base = start
                    hgvs_t.posedit.pos.end.base = end
                    hgvs_str = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(start) + '_' + str(end) + str(
                        hgvs_t.posedit.edit)
                    hgvs_t = hp.parse_hgvs_variant(hgvs_str)
                hgvs_c = copy.deepcopy(hgvs_t)

                # Set expanded out test to true
                expand_out = 'true'

            except Exception:
                hgvs_c = hgvs_c

        if re.match('NM_', str(hgvs_c.ac)):
            try:
                hgvs_c = no_norm_evm.n_to_c(hgvs_c)
            except hgvs.exceptions.HGVSError as e:
                hgvs_c = copy.deepcopy(stored_hgvs_c)

        # Ensure the altered c. variant has not crossed intro exon boundaries
        hgvs_check_boundaries = copy.deepcopy(hgvs_c)
        try:
            h_variant = hn.normalize(hgvs_check_boundaries)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if re.search('spanning the exon-intron boundary', error):
                hgvs_c = copy.deepcopy(stored_hgvs_c)
        # Catch identity at the exon/intron boundary by trying to normalize ref only
        if hgvs_check_boundaries.posedit.edit.type == 'identity':
            reform_ident = str(hgvs_c).split(':')[0]
            reform_ident = reform_ident + ':' + stored_hgvs_c.type + '.' + str(hgvs_c.posedit.pos) + 'del' + str(
                hgvs_c.posedit.edit.ref)  # + 'ins' + str(hgvs_c.posedit.edit.alt)
            hgvs_reform_ident = hp.parse_hgvs_variant(reform_ident)
            try:
                hn.normalize(hgvs_reform_ident)
            except hgvs.exceptions.HGVSError as e:
                error = str(e)
                if re.search('spanning the exon-intron boundary', error) or re.search(
                        'Normalization of intronic variants', error):
                    hgvs_c = copy.deepcopy(stored_hgvs_c)

    hgvs_genomic = vm.t_to_g(hgvs_c, alt_chr)
    if hgvs_c.posedit.edit.type == 'identity' and hgvs_genomic.posedit.edit.type == 'delins' and hgvs_genomic.posedit.edit.alt == '' and expand_out != 'true':
        hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.ref
    if hgvs_genomic.posedit.edit.type == 'ins' and utilise_gap_code is True:
        try:
            hgvs_genomic = hn.normalize(hgvs_genomic)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if error == 'insertion length must be 1':
                ref = sf.fetch_seq(str(hgvs_genomic.ac), hgvs_genomic.posedit.pos.start.base - 1,
                                   hgvs_genomic.posedit.pos.end.base)
                hgvs_genomic.posedit.edit.ref = ref
                hgvs_genomic.posedit.edit.alt = ref[0:1] + hgvs_genomic.posedit.edit.alt + ref[-1:]
                hgvs_genomic = hn.normalize(hgvs_genomic)
            if error == 'base start position must be <= end position':
                start = hgvs_genomic.posedit.pos.start.base
                end = hgvs_genomic.posedit.pos.end.base
                hgvs_genomic.posedit.pos.start.base = end
                hgvs_genomic.posedit.pos.end.base = start
                hgvs_genomic = hn.normalize(hgvs_genomic)

    # Statements required to reformat the stored_hgvs_c into a useable synonym
    if (stored_hgvs_c.posedit.edit.ref == '' or stored_hgvs_c.posedit.edit.ref is None) and expand_out != 'false':
        if stored_hgvs_c.type == 'c':
            stored_hgvs_n = vm.c_to_n(stored_hgvs_c)
        else:
            stored_hgvs_n = stored_hgvs_c
        stored_ref = sf.fetch_seq(str(stored_hgvs_n.ac), stored_hgvs_n.posedit.pos.start.base - 1,
                                  stored_hgvs_n.posedit.pos.end.base)
        stored_hgvs_c.posedit.edit.ref = stored_ref

    if (hgvs_genomic.posedit.edit.ref == '' or hgvs_genomic.posedit.edit.ref is None) and expand_out != 'false':
        if hgvs_genomic.posedit.edit.type == 'ins':
            stored_ref = sf.fetch_seq(str(hgvs_genomic.ac), hgvs_genomic.posedit.pos.start.base - 1,
                                      hgvs_genomic.posedit.pos.end.base)
            stored_alt = stored_ref[:1] + hgvs_genomic.posedit.edit.alt + stored_ref[-1:]
            hgvs_genomic.posedit.edit.ref = stored_ref
            hgvs_genomic.posedit.edit.alt = stored_alt

    # First look for variants mapping to the flanks of gaps
    # either in the gap or on the flank but not fully within the gap
    if expand_out == 'true':
        nr_genomic = nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)
        try:
            hn.normalize(nr_genomic)
        except hgvs.exceptions.HGVSInvalidVariantError as e:
            error_type_1 = str(e)
            if re.match('Length implied by coordinates must equal sequence deletion length', str(e)) or str(
                    e) == 'base start position must be <= end position':
                # Effectively, this code is designed to handle variants that are directly proximal to
                # gap BOUNDARIES, but in some cases the replace reference function of hgvs mapping has removed bases
                # due to the deletion length being > the specified range.

                # Warn of variant location wrt the gap
                if re.match('Length implied by coordinates must equal sequence deletion length', str(e)):
                    logger.warning('Variant is proximal to the flank of a genomic gap')
                    genomic_gap_variant = vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)
                    try:
                        hn.normalize(genomic_gap_variant)
                    # Still a problem
                    except hgvs.exceptions.HGVSInvalidVariantError as e:
                        if 'base start position must be <= end position' in str(e) and \
                                'Length implied by coordinates must equal' in error_type_1:
                            make_gen_var = copy.copy(nr_genomic)
                            make_gen_var.posedit.edit.ref = sf.fetch_seq(nr_genomic.ac,
                                                                         nr_genomic.posedit.pos.start.base - 1,
                                                                         nr_genomic.posedit.pos.end.base)
                            genomic_gap_variant = make_gen_var
                            error_type_1 = None
                    else:
                        genomic_gap_variant = nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)

                if error_type_1 == 'base start position must be <= end position':
                    logger.warning('Variant is fully within a genomic gap')
                    genomic_gap_variant = vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)

                # Logic
                # We have checked that the variant does not cross boundaries, or is intronic
                # So is likely mapping to a genomic gap
                try:
                    hn.normalize(genomic_gap_variant)
                except Exception as e:
                    if str(e) == 'base start position must be <= end position':
                        # This will only happen when the variant is fully within the gap
                        gap_start = genomic_gap_variant.posedit.pos.end.base
                        gap_end = genomic_gap_variant.posedit.pos.start.base
                        genomic_gap_variant.posedit.pos.start.base = gap_start
                        genomic_gap_variant.posedit.pos.end.base = gap_end
                    if re.match('Length implied by coordinates must equal sequence deletion length', str(e)):
                        # This will only happen if the variant is flanking the gap but is
                        # not inside the gap
                        logger.warning('Variant is on the flank of a genomic gap but not within the gap')
                        gap_start = genomic_gap_variant.posedit.pos.start.base - 1
                        gap_end = genomic_gap_variant.posedit.pos.end.base + 1
                        genomic_gap_variant.posedit.pos.start.base = gap_start
                        genomic_gap_variant.posedit.pos.end.base = gap_end
                        genomic_gap_variant.posedit.edit.ref = ''
                        stored_hgvs_c = copy.deepcopy(hgvs_c)

                        # Remove alt
                    try:
                        genomic_gap_variant.posedit.edit.alt = ''
                    except Exception as e:
                        if str(e) == "'Dup' object has no attribute 'alt'":
                            pass

                            # Should be a delins so will normalize statically and replace the reference bases
                    genomic_gap_variant = hn.normalize(genomic_gap_variant)
                    # Static map to c. and static normalize
                    transcript_gap_variant = vm.g_to_t(genomic_gap_variant, hgvs_c.ac)
                    stored_transcript_gap_variant = transcript_gap_variant
                    if not re.match('Length implied by coordinates must equal sequence deletion length', str(e)):
                        try:
                            transcript_gap_variant = hn.normalize(transcript_gap_variant)
                        except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                            if ' Unsupported normalization of variants spanning the UTR-exon boundary' in str(e):
                                pass

                    # if NM_ need the n. position
                    if re.match('NM_', str(hgvs_c.ac)):
                        transcript_gap_n = no_norm_evm.c_to_n(transcript_gap_variant)
                        transcript_gap_alt_n = no_norm_evm.c_to_n(stored_hgvs_c)
                    else:
                        transcript_gap_n = transcript_gap_variant
                        transcript_gap_alt_n = stored_hgvs_c

                    # Ensure an ALT exists
                    try:
                        if transcript_gap_alt_n.posedit.edit.alt is None:
                            transcript_gap_alt_n.posedit.edit.alt = 'X'
                    except Exception as e:
                        if str(e) == "'Dup' object has no attribute 'alt'":
                            transcript_gap_n_delins_from_dup = transcript_gap_n.ac + ':' + transcript_gap_n.type + '.' + str(
                                transcript_gap_n.posedit.pos.start.base) + '_' + str(
                                transcript_gap_n.posedit.pos.end.base) + 'del' + transcript_gap_n.posedit.edit.ref + 'ins' + transcript_gap_n.posedit.edit.ref + transcript_gap_n.posedit.edit.ref
                            transcript_gap_n = hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)
                            transcript_gap_alt_n_delins_from_dup = transcript_gap_alt_n.ac + ':' + transcript_gap_alt_n.type + '.' + str(
                                transcript_gap_alt_n.posedit.pos.start.base) + '_' + str(
                                transcript_gap_alt_n.posedit.pos.end.base) + 'del' + transcript_gap_alt_n.posedit.edit.ref + 'ins' + transcript_gap_alt_n.posedit.edit.ref + transcript_gap_alt_n.posedit.edit.ref
                            transcript_gap_alt_n = hp.parse_hgvs_variant(transcript_gap_alt_n_delins_from_dup)

                    # Split the reference and replacing alt sequence into a dictionary
                    reference_bases = list(transcript_gap_n.posedit.edit.ref)
                    if transcript_gap_alt_n.posedit.edit.alt is not None:
                        alternate_bases = list(transcript_gap_alt_n.posedit.edit.alt)
                    else:
                        # Deletions with no ins
                        pre_alternate_bases = list(transcript_gap_alt_n.posedit.edit.ref)
                        alternate_bases = []
                        for base in pre_alternate_bases:
                            alternate_bases.append('X')

                    # Create the dictionaries
                    ref_start = transcript_gap_n.posedit.pos.start.base
                    alt_start = transcript_gap_alt_n.posedit.pos.start.base
                    ref_base_dict = {}
                    for base in reference_bases:
                        ref_base_dict[ref_start] = str(base)
                        ref_start = ref_start + 1

                    alt_base_dict = {}

                    # Note, all variants will be forced into the format delete insert
                    # Deleted bases in the ALT will be substituted for X
                    for int in range(transcript_gap_alt_n.posedit.pos.start.base,
                                     transcript_gap_alt_n.posedit.pos.end.base + 1, 1):
                        if int == alt_start:
                            alt_base_dict[int] = str(''.join(alternate_bases))
                        else:
                            alt_base_dict[int] = 'X'

                    # Generate the alt sequence
                    alternate_sequence_bases = []
                    for int in range(transcript_gap_n.posedit.pos.start.base, transcript_gap_n.posedit.pos.end.base + 1,
                                     1):
                        if int in alt_base_dict.keys():
                            alternate_sequence_bases.append(alt_base_dict[int])
                        else:
                            alternate_sequence_bases.append(ref_base_dict[int])
                    alternate_sequence = ''.join(alternate_sequence_bases)
                    alternate_sequence = alternate_sequence.replace('X', '')

                    # Update variant, map to genome using vm and normalize
                    transcript_gap_n.posedit.edit.alt = alternate_sequence

                    try:
                        transcript_gap_variant = vm.n_to_c(transcript_gap_n)
                    except:
                        transcript_gap_variant = transcript_gap_n

                    try:
                        hgvs_genomic = vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                        hgvs_genomic = hn.normalize(hgvs_genomic)
                    except Exception as e:
                        if str(e) == "base start position must be <= end position":
                            # Expansion out is required to map back to the genomic position
                            pre_base = sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.start.base - 2,
                                                    transcript_gap_n.posedit.pos.start.base - 1)
                            post_base = sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.end.base,
                                                     transcript_gap_n.posedit.pos.end.base + 1)
                            transcript_gap_n.posedit.pos.start.base = transcript_gap_n.posedit.pos.start.base - 1
                            transcript_gap_n.posedit.pos.end.base = transcript_gap_n.posedit.pos.end.base + 1
                            transcript_gap_n.posedit.edit.ref = pre_base + transcript_gap_n.posedit.edit.ref + post_base
                            transcript_gap_n.posedit.edit.alt = pre_base + transcript_gap_n.posedit.edit.alt + post_base
                            try:
                                transcript_gap_variant = vm.n_to_c(transcript_gap_n)
                            except:
                                transcript_gap_variant = transcript_gap_n
                            hgvs_genomic = vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                            hgvs_genomic = hn.normalize(hgvs_genomic)

                    # Bypass the next bit of gap code
                    expand_out = 'false'

            else:
                pass
        # No map to the flank of a gap or within the gap
        else:
            pass

            # CASCADING STATEMENTS WHICH CAPTURE t to g MAPPING OPTIONS
    # Remove identity bases
    if hgvs_c == stored_hgvs_c:
        expand_out = 'false'
    elif expand_out == 'false' or utilise_gap_code is False:
        pass
    # Correct expansion ref + 2
    elif expand_out == 'true' and (
            len(hgvs_genomic.posedit.edit.ref) == (len(stored_hgvs_c.posedit.edit.ref) + 2)):  # >= 3:
        hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base + 1
        hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base - 1
        hgvs_genomic.posedit.edit.ref = hgvs_genomic.posedit.edit.ref[1:-1]
        if hgvs_genomic.posedit.edit.alt is not None:
            hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt[1:-1]
    elif expand_out == 'true' and (
            len(hgvs_genomic.posedit.edit.ref) != (len(stored_hgvs_c.posedit.edit.ref) + 2)):  # >= 3:
        if expand_out == 'true' and len(hgvs_genomic.posedit.edit.ref) == 2:
            gn = hn.normalize(hgvs_genomic)
            pass

        # Likely if the start or end position aligns to a gap in the genomic sequence
        # Logic
        # We have checked that the variant does not cross boundaries, or is intronic
        # So is likely mapping to a genomic gap
        elif expand_out == 'true' and len(hgvs_genomic.posedit.edit.ref) <= 1:
            # Incorrect expansion, likely < ref + 2
            genomic_gap_variant = vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)
            try:
                hn.normalize(genomic_gap_variant)
            except Exception as e:
                if str(e) == 'base start position must be <= end position':
                    gap_start = genomic_gap_variant.posedit.pos.end.base
                    gap_end = genomic_gap_variant.posedit.pos.start.base
                    genomic_gap_variant.posedit.pos.start.base = gap_start
                    genomic_gap_variant.posedit.pos.end.base = gap_end
                # Remove alt
                try:
                    genomic_gap_variant.posedit.edit.alt = ''
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        pass
                # Should be a delins so will normalize statically and replace the reference bases
                genomic_gap_variant = hn.normalize(genomic_gap_variant)
                # Static map to c. and static normalize
                transcript_gap_variant = vm.g_to_t(genomic_gap_variant, hgvs_c.ac)
                stored_transcript_gap_variant = transcript_gap_variant
                transcript_gap_variant = hn.normalize(transcript_gap_variant)
                # if NM_ need the n. position
                if re.match('NM_', str(hgvs_c.ac)):
                    transcript_gap_n = no_norm_evm.c_to_n(transcript_gap_variant)
                    transcript_gap_alt_n = no_norm_evm.c_to_n(stored_hgvs_c)
                else:
                    transcript_gap_n = transcript_gap_variant
                    transcript_gap_alt_n = stored_hgvs_c

                # Ensure an ALT exists
                try:
                    if transcript_gap_alt_n.posedit.edit.alt is None:
                        transcript_gap_alt_n.posedit.edit.alt = 'X'
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        transcript_gap_n_delins_from_dup = transcript_gap_n.ac + ':' + transcript_gap_n.type + '.' + str(
                            transcript_gap_n.posedit.pos.start.base) + '_' + str(
                            transcript_gap_n.posedit.pos.end.base) + 'del' + transcript_gap_n.posedit.edit.ref + 'ins' + transcript_gap_n.posedit.edit.ref + transcript_gap_n.posedit.edit.ref
                        transcript_gap_n = hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)
                        transcript_gap_alt_n_delins_from_dup = transcript_gap_alt_n.ac + ':' + transcript_gap_alt_n.type + '.' + str(
                            transcript_gap_alt_n.posedit.pos.start.base) + '_' + str(
                            transcript_gap_alt_n.posedit.pos.end.base) + 'del' + transcript_gap_alt_n.posedit.edit.ref + 'ins' + transcript_gap_alt_n.posedit.edit.ref + transcript_gap_alt_n.posedit.edit.ref
                        transcript_gap_alt_n = hp.parse_hgvs_variant(transcript_gap_alt_n_delins_from_dup)

                # Split the reference and replacing alt sequence into a dictionary
                reference_bases = list(transcript_gap_n.posedit.edit.ref)
                if transcript_gap_alt_n.posedit.edit.alt is not None:
                    alternate_bases = list(transcript_gap_alt_n.posedit.edit.alt)
                else:
                    # Deletions with no ins
                    pre_alternate_bases = list(transcript_gap_alt_n.posedit.edit.ref)
                    alternate_bases = []
                    for base in pre_alternate_bases:
                        alternate_bases.append('X')

                # Create the dictionaries
                ref_start = transcript_gap_n.posedit.pos.start.base
                alt_start = transcript_gap_alt_n.posedit.pos.start.base
                ref_base_dict = {}
                for base in reference_bases:
                    ref_base_dict[ref_start] = str(base)
                    ref_start = ref_start + 1

                alt_base_dict = {}

                # Note, all variants will be forced into the format delete insert
                # Deleted bases in the ALT will be substituted for X
                for int in range(transcript_gap_alt_n.posedit.pos.start.base,
                                 transcript_gap_alt_n.posedit.pos.end.base + 1, 1):
                    if int == alt_start:
                        alt_base_dict[int] = str(''.join(alternate_bases))
                    else:
                        alt_base_dict[int] = 'X'

                # Generate the alt sequence
                alternate_sequence_bases = []
                for int in range(transcript_gap_n.posedit.pos.start.base, transcript_gap_n.posedit.pos.end.base + 1, 1):
                    if int in alt_base_dict.keys():
                        alternate_sequence_bases.append(alt_base_dict[int])
                    else:
                        alternate_sequence_bases.append(ref_base_dict[int])
                alternate_sequence = ''.join(alternate_sequence_bases)
                alternate_sequence = alternate_sequence.replace('X', '')

                # Update variant, map to genome using vm and normalize
                transcript_gap_n.posedit.edit.alt = alternate_sequence

                try:
                    transcript_gap_variant = vm.n_to_c(transcript_gap_n)
                except:
                    transcript_gap_variant = transcript_gap_n

                try:
                    hgvs_genomic = vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                    hgvs_genomic = hn.normalize(hgvs_genomic)
                except Exception as e:
                    if str(e) == "base start position must be <= end position":
                        # Expansion out is required to map back to the genomic position
                        pre_base = sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.start.base - 2,
                                                transcript_gap_n.posedit.pos.start.base - 1)
                        post_base = sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.end.base,
                                                 transcript_gap_n.posedit.pos.end.base + 1)
                        transcript_gap_n.posedit.pos.start.base = transcript_gap_n.posedit.pos.start.base - 1
                        transcript_gap_n.posedit.pos.end.base = transcript_gap_n.posedit.pos.end.base + 1
                        transcript_gap_n.posedit.edit.ref = pre_base + transcript_gap_n.posedit.edit.ref + post_base
                        transcript_gap_n.posedit.edit.alt = pre_base + transcript_gap_n.posedit.edit.alt + post_base
                        try:
                            transcript_gap_variant = vm.n_to_c(transcript_gap_n)
                        except:
                            transcript_gap_variant = transcript_gap_n
                        hgvs_genomic = vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                        hgvs_genomic = hn.normalize(hgvs_genomic)

    # Ins variants map badly - Especially between c. exon/exon boundary
    if hgvs_c.posedit.edit.type == 'ins' and hgvs_c.posedit.pos.start.offset == 0 and hgvs_c.posedit.pos.end.offset == 0:
        try:
            hn.normalize(hgvs_genomic)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if error == 'insertion length must be 1':
                if hgvs_c.type == 'c':
                    hgvs_t = vm.c_to_n(hgvs_c)
                else:
                    hgvs_t = copy.copy(hgvs_c)
                ins_ref = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 1, hgvs_t.posedit.pos.end.base)
                ins_alt = ins_ref[:1] + hgvs_t.posedit.edit.alt + ins_ref[-1:]
                ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(
                    hgvs_t.posedit.pos.end.base) + 'del' + ins_ref + 'ins' + ins_alt
                hgvs_t = hp.parse_hgvs_variant(ins_to_delins)
                try:
                    hgvs_c = vm.n_to_c(hgvs_t)
                except Exception:
                    hgvs_c = copy.copy(hgvs_t)
                try:
                    hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)
                except Exception as e:
                    error = str(e)
                    logger.warning('Ins mapping error in myt_to_g ' + error)

    return hgvs_genomic


"""
Simple hgvs g. to c. or n. mapping
returns parsed hgvs c. or n. object
"""


def myevm_g_to_t(evm, hgvs_genomic, alt_ac):
    hgvs_t = evm.g_to_t(hgvs_genomic, alt_ac)
    return hgvs_t


"""
parse p. strings into hgvs p. objects

MARKED FOR REMOVAL
"""

# def hgvs_protein(variant, hp):
#     # Set regular expressions for if statements
#     pat_p = re.compile(":p\.")  # Pattern looks for :g. Note (gene) has been removed
#     # If the :p. pattern is present in the input variant
#     if pat_p.search(variant):
#         # convert the input string into a hgvs object
#         var_p = hp.parse_hgvs_variant(variant)
#         return var_p


"""
Convert r. into c.
"""


def hgvs_r_to_c(hgvs_object):
    # check for LRG_t with r.
    if re.match('LRG', hgvs_object.ac):
        transcript_ac = dbControls.data.get_RefSeqTranscriptID_from_lrgTranscriptID(hgvs_object.ac)
        if transcript_ac == 'none':
            raise HGVSDataNotAvailableError('Unable to identify a relevant transcript for ' + hgvs_object.ac)
        else:
            hgvs_object.ac = transcript_ac
    hgvs_object.type = 'c'
    edit = str(hgvs_object.posedit.edit)
    edit = edit.upper()
    # lowercase the supported variant types
    edit = edit.replace('DEL', 'del')
    edit = edit.replace('INS', 'ins')
    edit = edit.replace('INV', 'inv')
    edit = edit.replace('DUP', 'dup')
    # edit = edit.replace('CON', 'con')
    # edit = edit.replace('TRA', 'tra')
    edit = edit.replace('U', 'T')
    hgvs_object.posedit.edit = edit
    return hgvs_object


"""
Convert c. into r.
"""


def hgvs_c_to_r(hgvs_object):
    hgvs_object.type = 'r'
    edit = str(hgvs_object.posedit.edit)
    edit = edit.lower()
    edit = edit.replace('t', 'u')
    hgvs_object.posedit.edit = edit
    return hgvs_object


"""
Input c. r. n. variant string
Use uta.py (hdp) to return the identity information for the transcript variant 
see hgvs.dataproviders.uta.py for details

MARKED FOR REMOVAL
"""

# def tx_identity_info(variant, hdp):
#     # Set regular expressions for if statements
#     pat_c = re.compile(":c\.")  # Pattern looks for :c. Note (gene) has been removed
#     pat_n = re.compile(":n\.")  # Pattern looks for :c. Note (gene) has been removed
#     pat_r = re.compile(":r\.")  # Pattern looks for :c. Note (gene) has been removed
#
#     # If the :c. pattern is present in the input variant
#     if pat_c.search(variant):
#         # Remove all text to the right and including pat_c
#         tx_ac = variant[:variant.index(':c.') + len(':c.')]
#         tx_ac = pat_c.sub('', tx_ac)
#         # Interface with the UTA database via get_tx_identity in uta.py
#         tx_id_info = hdp.get_tx_identity_info(tx_ac)
#         # NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
#         return tx_id_info
#
#     # If the :n. pattern is present in the input variant
#     if pat_n.search(variant):
#         # Remove all text to the right and including pat_c
#         tx_ac = variant[:variant.index(':n.') + len(':n.')]
#         tx_ac = pat_n.sub('', tx_ac)
#         # Interface with the UTA database via get_tx_identity in uta.py
#         tx_id_info = hdp.get_tx_identity_info(tx_ac)
#         # NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
#         return tx_id_info
#
#     # If the :r. pattern is present in the input variant
#     if pat_r.search(variant):
#         # Remove all text to the right and including pat_c
#         tx_ac = variant[:variant.index(':r.') + len(':r.')]
#         tx_ac = pat_r.sub('', tx_ac)
#         # Interface with the UTA database via get_tx_identity in uta.py
#         tx_id_info = hdp.get_tx_identity_info(tx_ac)
#         # NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
#         return tx_id_info


"""
Input c. r. nd accession string
Use uta.py (hdp) to return the identity information for the transcript variant 
see hgvs.dataproviders.uta.py for details

MARKED FOR REMOVAL
"""

# def tx_id_info(alt_ac, hdp):
#     tx_id_info = hdp.get_tx_identity_info(alt_ac)
#     # NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
#     return tx_id_info


"""
Use uta.py (hdp) to return the transcript information for a specified gene (HGNC SYMBOL)
see hgvs.dataproviders.uta.py for details


marked for removal
"""


def tx_for_gene(hgnc, hdp):
    # Interface with the UTA database via get_tx_for_gene in uta.py
    tx_for_gene = hdp.get_tx_for_gene(hgnc)
    return tx_for_gene


"""
Extract RefSeqGene Accession from transcript information
see hgvs.dataproviders.uta.py for details
"""


def ng_extract(tx_for_gene):
    # Set regular expressions for if statements
    pat_NG = re.compile("^NG_")  # Pattern looks for NG_ at beginning of a string
    # For each list in the list of lists tx_for_gene
    for list in tx_for_gene:
        # If the pattern NG_ is found in element 4
        if pat_NG.search(list[4]):
            # The gene accession is set to list element 4
            gene_ac = list[4]
            return gene_ac


"""
marked for removal
"""
# def int_start(var_g):
#   start = var_g.posedit.pos.start
#   # Stringify to get start co-ords
#   start = str(start)
#   # Make into an integer
#   int_start = int(start)
#   return int_start

"""
marked for removal
"""
# def int_end(var_g):
#   end = var_g.posedit.pos.end
#   # Stringify to get start co-ords
#   end = str(end)
#   # Make into an integer
#   int_end = int(end)
#   return int_end

"""
Returns exon information for a given transcript
e.g. how the exons align to the genomic reference
see hgvs.dataproviders.uta.py for details
"""


def tx_exons(tx_ac, alt_ac, alt_aln_method, hdp):
    # Interface with the UTA database via get_tx_exons in uta.py
    try:
        tx_exons = hdp.get_tx_exons(tx_ac, alt_ac, alt_aln_method)
    except hgvs.exceptions.HGVSError as e:
        tx_exons = 'hgvs Exception: ' + str(e)
        return tx_exons
    try:
        tx_exons[0]['alt_strand']
    except TypeError:
        tx_exons = 'error'
        return tx_exons
    # If on the reverse strand, reverse the order of elements
    if tx_exons[0]['alt_strand'] == -1:
        tx_exons = tx_exons[::-1]
        return tx_exons
    else:
        return tx_exons


"""
Automatically maps genomic positions onto all overlapping transcripts
"""


def relevant_transcripts(hgvs_genomic, evm, hdp, alt_aln_method, reverse_normalizer, hp):
    reverse_hn = reverse_normalizer
    # Pass relevant transcripts for the input variant to rts
    # Note, the evm method misses one end, the hdp. method misses the other. Combine both
    rts_list = hdp.get_tx_for_region(hgvs_genomic.ac, alt_aln_method, hgvs_genomic.posedit.pos.start.base - 1,
                                     hgvs_genomic.posedit.pos.end.base - 1)
    rts_dict = {}
    for tx_dat in rts_list:
        rts_dict[tx_dat[0]] = True
    rts_list_2 = evm.relevant_transcripts(hgvs_genomic)
    for tx_dat_2 in rts_list_2:
        rts_dict[tx_dat_2] = True
    rts = rts_dict.keys()

    # Project genomic variants to new transcripts
    # and  populate a code_var list
    #############################################
    # Open a list to store relevant transcripts
    code_var = []
    # Populate transcripts - The keys become the list elements from rel_trs
    for x in rts:
        y = x.rstrip()  # Chomp any whitespace from the right of x ($_) - Assign to y
        # Easy variant mapper used to map the input variant to the relevant transcripts
        # Check for coding transcripts
        try:
            variant = evm.g_to_t(hgvs_genomic, y)
        except hgvs.exceptions.HGVSError as e:
            # Check for non-coding transcripts
            try:
                variant = evm.g_to_t(hgvs_genomic, y)
            except hgvs.exceptions.HGVSError as e:
                continue
        except:
            continue

        # Corrective Normalisation of intronic descriptions in the antisense oriemtation
        pl = re.compile('\+')
        mi = re.compile('\-')
        ast = re.compile('\*')
        if pl.search(str(variant)) or mi.search(str(variant)) or ast.search(str(variant)):
            tx_ac = variant.ac
            alt_ac = hgvs_genomic.ac

            # Interface with the UTA database via get_tx_exons in uta.py
            try:
                tx_exons = hdp.get_tx_exons(tx_ac, alt_ac, alt_aln_method)
            except hgvs.exceptions.HGVSError as e:
                e
                tx_exons = 'hgvs Exception: ' + str(e)
                return tx_exons
            try:
                completion = tx_exons[0]['alt_strand']
            except TypeError:
                tx_exons = 'error'
                return tx_exons
            # If on the reverse strand, reverse the order of elements
            if tx_exons[0]['alt_strand'] == -1:
                tx_exons = tx_exons[::-1]
            else:
                pass

            # Gene orientation
            if tx_exons[0]['alt_strand'] == -1:
                antisense = 'true'
            else:
                antisense = 'false'

            # Pass if antisense = 'false'
            if antisense == 'false':
                pass
            else:
                # Reverse normalize hgvs_genomic
                rev_hgvs_genomic = reverse_hn.normalize(hgvs_genomic)
                # map back to coding
                variant = evm.g_to_t(rev_hgvs_genomic, tx_ac)

        strung = str(variant)
        try:
            hp.parse_hgvs_variant(strung)
        except hgvs.exceptions.HGVSError:
            continue
        except TypeError:
            continue
        else:
            code_var.append(str(variant))
    return code_var


"""
Take HGVS string, parse into hgvs object and validate
"""


def validate(input, hp, vr):
    hgvs_input = hp.parse_hgvs_variant(input)
    g = re.compile(":g.")
    p = re.compile(":p.")
    if p.search(input):
        if hasattr(hgvs_input.posedit.pos.start, 'offset'):
            pass
        else:
            hgvs_input.posedit.pos.start.offset = 0
        if hasattr(hgvs_input.posedit.pos.end, 'offset'):
            pass
        else:
            hgvs_input.posedit.pos.end.offset = 0
        if hasattr(hgvs_input.posedit.pos.start, 'datum'):
            pass
        else:
            hgvs_input.posedit.pos.start.datum = 0
        if hasattr(hgvs_input.posedit.pos.end, 'datum'):
            pass
        else:
            hgvs_input.posedit.pos.end.datum = 0
        if hasattr(hgvs_input.posedit.edit, 'ref_n'):
            pass
        else:
            hgvs_input.posedit.edit.ref_n = hgvs_input.posedit.pos.end.base - hgvs_input.posedit.pos.start.base + 1

    try:
        vr.validate(hgvs_input)
    except hgvs.exceptions.HGVSError as e:

        error = e
        return error

    else:
        error = 'false'
        return error


"""
marked for removal
"""
# def sequence_extractor(ac, hdp):
#   ac_seq = hdp.get_tx_seq(ac)
#   return ac_seq

"""
marked for removal
"""
# def ref_replace(e, hgvs_variant):
#   error = str(e)
#   match = re.findall('\(([GATC]+)\)', error)
#   new_ref = match[1]
#   hgvs_variant.posedit.edit.ref = new_ref
#   return hgvs_variant

"""
Search HGNC rest
"""


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


"""
Search Entrez databases with efetch and SeqIO
"""


def entrez_efetch(db, id, rettype, retmode):
    # IMPORT Bio modules
    # from Bio import Entrez
    Entrez.email = ENTREZ_ID
    # from Bio import SeqIO
    handle = Entrez.efetch(db=db, id=id, rettype=rettype, retmode=retmode)
    # Get record
    record = SeqIO.read(handle, "gb")
    # Place into text
    # text = handle.read()
    handle.close()
    return record


"""
search Entrez databases with efetch and read
"""


def entrez_read(db, id, retmode):
    # IMPORT Bio modules
    # from Bio import Entrez
    Entrez.email = ENTREZ_ID
    # from Bio import SeqIO
    handle = Entrez.efetch(db=db, id=id, retmode=retmode)
    # Get record
    record = Entrez.read(handle)
    # Place into text
    # text = handle.read()
    handle.close()
    return record


"""
Simple reverse complement function for nucleotide sequences
"""


def revcomp(bases):
    l2 = []
    l = list(bases)
    element = 0
    for base in l:
        element = element + 1
        if base == 'G':
            l2.append('C')
        if base == 'C':
            l2.append('G')
        if base == 'A':
            l2.append('T')
        if base == 'T':
            l2.append('A')
    revcomp = ''.join(l2)
    revcomp = revcomp[::-1]
    return revcomp


"""
Function designed to merge multiple HGVS variants (hgvs objects) into a single delins 
using 3 prime normalization
"""


def merge_hgvs_3pr(hgvs_variant_list, hp, vr, hn, vm, sf):
    # Ensure c. is mapped to the
    h_list = []

    # Sanity check and format the submitted variants
    for hgvs_v in hgvs_variant_list:
        # For testing include parser
        try:
            hgvs_v = hp.parse_hgvs_variant(hgvs_v)
        except Exception as e:
            print e
            pass

        # Validate
        vr.validate(hgvs_v)  # Let hgvs errors deal with invalid variants and not hgvs objects
        if hgvs_v.type == 'c':
            try:
                hgvs_v = vm.c_to_n(hgvs_v)
                h_list.append(hgvs_v)
            except:
                raise mergeHGVSerror("Unable to map from c. position to absolute position")
        elif hgvs_v.type == 'g':
            h_list.append(hgvs_v)
    if h_list != []:
        hgvs_variant_list = copy.deepcopy(h_list)

    # Define accession and start/end positions
    accession = None
    merge_start_pos = None
    merge_end_pos = None
    type = None
    full_list = []

    # Loop through the submitted variants and gather the required info
    for hgvs_v in hgvs_variant_list:
        # No intronic positions
        try:
            if hgvs_v.posedit.pos.start.offset != 0:
                raise mergeHGVSerror("Base-offset position submitted")
            if hgvs_v.posedit.pos.end.offset != 0:
                raise mergeHGVSerror("Base-offset position submitted")
        except AttributeError:
            pass

        # Normalize the variant (allow cross intron) which also adds the reference sequence (?)
        hgvs_v = hn.normalize(hgvs_v)

        # Set the accession and ensure that multiple reference sequences have not been queried
        if accession is None:
            accession = hgvs_v.ac
            type = hgvs_v.type
        else:
            if hgvs_v.ac != accession:
                raise mergeHGVSerror("More than one reference sequence submitted")
            else:
                pass

        # Set initial start and end positions
        if merge_start_pos is None:
            merge_start_pos = hgvs_v.posedit.pos.start.base
            merge_end_pos = hgvs_v.posedit.pos.end.base
            # Append to the final list of variants
            full_list.append(hgvs_v)
            continue
        # Ensure variants are in the correct order and not overlapping
        else:
            # ! hgvs_v.posedit.pos.start.base !>
            if hgvs_v.posedit.pos.start.base <= merge_end_pos:
                raise mergeHGVSerror("Submitted variants are out of order or their ranges overlap")
            else:
                # Create a fake variant to handle the missing sequence
                ins_seq = sf.fetch_seq(hgvs_v.ac, merge_end_pos, hgvs_v.posedit.pos.start.base - 1)
                gapping = hgvs_v.ac + ':' + hgvs_v.type + '.' + str(merge_end_pos + 1) + '_' + str(
                    hgvs_v.posedit.pos.start.base - 1) + 'delins' + ins_seq
                hgvs_gapping = hp.parse_hgvs_variant(gapping)
                full_list.append(hgvs_gapping)
                # update end_pos
                merge_end_pos = hgvs_v.posedit.pos.end.base
                # Append to the final list of variants
                full_list.append(hgvs_v)

    # Generate the alt sequence
    alt_sequence = ''
    for hgvs_v in full_list:
        ref_alt = hgvs2vcf.hgvs_ref_alt(hgvs_v, sf)
        alt_sequence = alt_sequence + ref_alt['alt']

    # Fetch the reference sequence and copy it for the basis of the alt sequence
    reference_sequence = sf.fetch_seq(accession, merge_start_pos - 1, merge_end_pos)
    # Generate an hgvs_delins
    if alt_sequence == '':
        delins = accession + ':' + type + '.' + str(merge_start_pos) + '_' + str(
            merge_end_pos) + 'del' + reference_sequence
    else:
        delins = accession + ':' + type + '.' + str(merge_start_pos) + '_' + str(
            merge_end_pos) + 'del' + reference_sequence + 'ins' + alt_sequence
    hgvs_delins = hp.parse_hgvs_variant(delins)
    try:
        hgvs_delins = vm.n_to_c(hgvs_delins)
    except:
        pass
    # Normalize (allow variants crossing into different exons)
    try:
        hgvs_delins = hn.normalize(hgvs_delins)
    except HGVSUnsupportedOperationError:
        pass
    return hgvs_delins


"""
Function designed to merge multiple HGVS variants (hgvs objects) into a single delins 
using 5 prime normalization
"""


def merge_hgvs_5pr(hgvs_variant_list, hp, vr, reverse_normalizer, vm, sf):
    reverse_hn = reverse_normalizer

    # Ensure c. is mapped to the
    h_list = []

    # Sanity check and format the submitted variants
    for hgvs_v in hgvs_variant_list:
        # For testing include parser
        try:
            hgvs_v = hp.parse_hgvs_variant(hgvs_v)
        except:
            pass

        # Validate
        vr.validate(hgvs_v)  # Let hgvs errors deal with invalid variants and not hgvs objects
        if hgvs_v.type == 'c':
            try:
                hgvs_v = vm.c_to_n(hgvs_v)
                h_list.append(hgvs_v)
            except:
                raise mergeHGVSerror("Unable to map from c. position to absolute position")
    if h_list != []:
        hgvs_variant_list = copy.deepcopy(h_list)

    # Define accession and start/end positions
    accession = None
    merge_start_pos = None
    merge_end_pos = None
    type = None
    full_list = []

    # Loop through the submitted variants and gather the required info
    for hgvs_v in hgvs_variant_list:
        try:
            # No intronic positions
            if hgvs_v.posedit.pos.start.offset != 0:
                raise mergeHGVSerror("Base-offset position submitted")
            if hgvs_v.posedit.pos.end.offset != 0:
                raise mergeHGVSerror("Base-offset position submitted")
        except AttributeError:
            pass

        # Normalize the variant (allow cross intron) which also adds the reference sequence (?)
        hgvs_v = reverse_hn.normalize(hgvs_v)

        # Set the accession and ensure that multiple reference sequences have not been queried
        if accession is None:
            accession = hgvs_v.ac
            type = hgvs_v.type
        else:
            if hgvs_v.ac != accession:
                raise mergeHGVSerror("More than one reference sequence submitted")
            else:
                pass

        # Set initial start and end positions
        if merge_start_pos is None:
            merge_start_pos = hgvs_v.posedit.pos.start.base
            merge_end_pos = hgvs_v.posedit.pos.end.base
            # Append to the final list of variants
            full_list.append(hgvs_v)
            continue
        # Ensure variants are in the correct order and not overlapping
        else:
            # ! hgvs_v.posedit.pos.start.base !>
            if hgvs_v.posedit.pos.start.base <= merge_end_pos:
                raise mergeHGVSerror("Submitted variants are out of order or their ranges overlap")
            else:
                # Create a fake variant to handle the missing sequence
                ins_seq = sf.fetch_seq(hgvs_v.ac, merge_end_pos, hgvs_v.posedit.pos.start.base - 1)
                gapping = hgvs_v.ac + ':' + hgvs_v.type + '.' + str(merge_end_pos + 1) + '_' + str(
                    hgvs_v.posedit.pos.start.base - 1) + 'delins' + ins_seq
                hgvs_gapping = hp.parse_hgvs_variant(gapping)
                full_list.append(hgvs_gapping)
                # update end_pos
                merge_end_pos = hgvs_v.posedit.pos.end.base
                # Append to the final list of variants
                full_list.append(hgvs_v)

    # Generate the alt sequence
    alt_sequence = ''
    for hgvs_v in full_list:
        ref_alt = hgvs2vcf.hgvs_ref_alt(hgvs_v, sf)
        alt_sequence = alt_sequence + ref_alt['alt']

    # Fetch the reference sequence and copy it for the basis of the alt sequence
    reference_sequence = sf.fetch_seq(accession, merge_start_pos - 1, merge_end_pos)

    # Generate an hgvs_delins
    if alt_sequence == '':
        delins = accession + ':' + type + '.' + str(merge_start_pos) + '_' + str(
            merge_end_pos) + 'del' + reference_sequence
    else:
        delins = accession + ':' + type + '.' + str(merge_start_pos) + '_' + str(
            merge_end_pos) + 'del' + reference_sequence + 'ins' + alt_sequence
    hgvs_delins = hp.parse_hgvs_variant(delins)
    try:
        hgvs_delins = vm.n_to_c(hgvs_delins)
    except:
        pass
    # Normalize (allow variants crossing into different exons)
    try:
        hgvs_delins = reverse_hn.normalize(hgvs_delins)
    except HGVSUnsupportedOperationError:
        pass
    return hgvs_delins


"""
Function designed to merge multiple pseudo VCF variants (strings) into a single HGVS delins 
using 5 prime normalization then return a 3 prime normalized final HGVS object
"""


def merge_pseudo_vcf(vcf_list, genome_build, reverse_normalizer, hn, hp):
    hgvs_list = []
    # Convert pseudo_vcf list into a HGVS list
    normalization_direction = 5
    for call in vcf_list:
        hgvs = pseudo_vcf2hgvs.pvcf_to_hgvs(call, genome_build, normalization_direction, reverse_normalizer, hn, hp)
        hgvs_list.append(hgvs)
    # Merge
    hgvs_delins = merge_hgvs_5pr(hgvs_list)
    # normalize 3 prime
    hgvs_delins = hn.normalize(hgvs_delins)
    # return
    return hgvs_delins


"""
HGVS allele handling function which takes a single HGVS allele description and 
separates each allele into a list of HGVS variants
"""


def hgvs_alleles(variant_description, hp, vr, hn, vm, sf):
    try:
        # Split up the description
        accession, remainder = variant_description.split(':')
        # Branch
        if re.search('[gcn]\.\d+\[', remainder):
            # NM_004006.2:c.2376[G>C];[(G>C)]
            # if re.search('\(', remainder):
            #   raise alleleVariantError('Unsupported format ' + remainder)
            # NM_004006.2:c.2376[G>C];[G>C]
            type, remainder = remainder.split('.')
            pos = re.match('\d+', remainder)
            pos = pos.group(0)
            remainder = remainder.replace(pos, '')
            remainder = remainder[1:-1]
            alleles = remainder.split('];[')
            my_alleles = []
            for posedit in alleles:
                if re.search('\(', posedit):
                    # NM_004006.2:c.2376[G>C];[(G>C)]
                    continue
                posedit_list = [posedit]
                current_allele = []
                for pe in posedit_list:
                    vrt = accession + ':' + type + '.' + str(pos) + pe
                    current_allele.append(vrt)
                my_alleles.append(current_allele)
        else:
            type, remainder = remainder.split('.')
            if re.search('\(;\)', remainder) and re.search('\];', remainder):
                # NM_004006.2:c.[296T>G];[476T>C](;)1083A>C(;)1406del
                pre_alleles = remainder.split('(;)')
                pre_merges = []
                alleles = []
                for allele in pre_alleles:
                    if re.match('\[', allele):
                        pre_merges.append(allele)
                    else:
                        alleles.append(allele)
                # Extract descriptions
                my_alleles = []
                # First alleles
                for posedits in alleles:
                    posedit_list = posedits.split(';')
                    current_allele = []
                    for pe in posedit_list:
                        vrt = accession + ':' + type + '.' + pe
                        current_allele.append(vrt)
                    my_alleles.append(current_allele)

                # Then Merges
                alleles = []
                remainder = ';'.join(pre_merges)
                remainder = remainder[1:-1]  # removes the first [ and the last ]
                alleles = remainder.split('];[')
                # now separate out the variants in each allele§
                for posedits in alleles:
                    posedit_list = posedits.split(';')
                    current_allele = []
                    for pe in posedit_list:
                        vrt = accession + ':' + type + '.' + pe
                        current_allele.append(vrt)
                    my_alleles.append(current_allele)
                # Now merge the alleles into a single variant
                merged_alleles = []
                for each_allele in my_alleles:
                    if re.search('\?', str(each_allele)):
                        # NM_004006.2:c.[2376G>C];[?]
                        continue
                    merge = []
                    allele = str(merge_hgvs_3pr(each_allele, hp, vr, hn, vm, sf))
                    merge.append(allele)
                    for variant in each_allele:
                        merged_alleles.append([variant])
                    # merged_alleles.append(merge)
                my_alleles = merged_alleles

            elif re.search('\(;\)', remainder):
                # If statement for uncertainties
                # NM_004006.2:c.[296T>G;476C>T];[476C>T](;)1083A>C
                if re.search('\[', remainder):
                    raise alleleVariantError('Unsupported format ' + type + '.' + remainder)
                # NM_004006.2:c.2376G>C(;)3103del
                # NM_000548.3:c.3623_3647del(;)3745_3756dup
                alleles = remainder.split('(;)')
                # now separate out the variants in each allele§
                my_alleles = []
                for posedits in alleles:
                    posedit_list = posedits.split(';')
                    current_allele = []
                    for pe in posedit_list:
                        vrt = accession + ':' + type + '.' + pe
                        current_allele.append(vrt)
                    my_alleles.append(current_allele)
            else:
                # If statement for uncertainties
                if re.search('\(', remainder):
                    raise alleleVariantError('Unsupported format ' + type + '.' + remainder)
                # NM_004006.2:c.[2376G>C];[3103del]
                # NM_004006.2:c.[2376G>C];[3103del]
                # NM_004006.2:c.[296T>G;476C>T;1083A>C];[296T>G;1083A>C]
                # NM_000548.3:c.[4358_4359del;4361_4372del]
                remainder = remainder[1:-1]  # removes the first [ and the last ]
                alleles = remainder.split('];[')
                # now separate out the variants in each allele§
                my_alleles = []
                for posedits in alleles:
                    posedit_list = posedits.split(';')
                    current_allele = []
                    for pe in posedit_list:
                        vrt = accession + ':' + type + '.' + pe
                        current_allele.append(vrt)
                    my_alleles.append(current_allele)
                # Now merge the alleles into a single variant
                merged_alleles = []

                for each_allele in my_alleles:
                    if re.search('\?', str(each_allele)):
                        # NM_004006.2:c.[2376G>C];[?]
                        continue
                    merge = []
                    allele = str(merge_hgvs_3pr(each_allele, hp, vr, hn, vm, sf))
                    merge.append(allele)
                    for variant in each_allele:
                        merged_alleles.append([variant])
                   # merged_alleles.append(merge)

                my_alleles = merged_alleles

        # Extract alleles into strings
        allele_strings = []
        for alleles_l in my_alleles:
            for allele in alleles_l:
                allele_strings.append(allele)
        my_alleles = allele_strings

        # return
        return my_alleles
    except Exception as e:
        import traceback
        exc_type, exc_value, last_traceback = sys.exc_info()
        te = traceback.format_exc()
        raise alleleVariantError(str(e))

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
