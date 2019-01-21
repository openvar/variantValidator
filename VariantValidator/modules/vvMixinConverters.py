import re
import os
import sys
import copy
from vvLogging import logger
import hgvs
import hgvs.exceptions
from hgvs.dataproviders import uta
from hgvs.dataproviders import seqfetcher
import hgvs.normalizer
import hgvs.validator
import hgvs.parser
import hgvs.variantmapper
import hgvs.sequencevariant
import vvMixinInit
import vvChromosomes
import vvHGVS
from urlparse import urlparse
import httplib2 as http
import json
from Bio import Entrez,SeqIO



#Error setup
from hgvs.exceptions import HGVSError, HGVSDataNotAvailableError, HGVSUnsupportedOperationError
class mergeHGVSerror(Exception):
    pass
class alleleVariantError(Exception):
    pass



class Mixin(vvMixinInit.Mixin):
    """
    r_to_c
    parses r. variant strings into hgvs object and maps to the c. equivalent.
    """
    def r_to_c(self, variant, evm):
        # convert the input string into a hgvs object by parsing
        var_r = self.hp.parse_hgvs_variant(variant)
        # map to the coding sequence
        var_c = evm.r_to_c(var_r)  # coding level variant
        variant = str(var_c)
        c_from_r = {'variant': variant, 'type': ':c.'}
        return c_from_r

    """ 
    Maps transcript variant descriptions onto specified RefSeqGene reference sequences
    Return an hgvs object containing the genomic sequence variant relative to the RefSeqGene 
    acession
    refseq_ac = RefSeqGene ac
    """


    def refseq(self, variant, vmOld, refseq_ac, hpOld, evm, hdpOld, primary_assembly):
        vr = hgvs.validator.Validator(self.hdp)
        # parse the variant into hgvs object
        var_c = self.hp.parse_hgvs_variant(variant)
        # map to the genomic co-ordinates using the easy variant mapper set to alt_aln_method = alt_aln_method
        var_g = self.myevm_t_to_g(var_c, evm, self.hdp, primary_assembly)
        # Get overlapping transcripts - forcing a splign alignment
        start_i = var_g.posedit.pos.start.base
        end_i = var_g.posedit.pos.end.base
        alt_ac = var_g.ac
        alt_aln_method = 'splign'
        transcripts = self.hdp.get_tx_for_region(alt_ac, alt_aln_method, start_i - 1, end_i)
        # Take the first transcript
        for trans in transcripts:
            tx_ac = trans[0]
            try:
                ref_c = self.vm.g_to_t(var_g, tx_ac, alt_aln_method='splign')
            except:
                continue
            else:
                # map the variant co-ordinates to the refseq Gene accession using vm
                ref_g_dict = {
                    'ref_g': '',
                    'error': 'false'
                }
                try:
                    ref_g_dict['ref_g'] = self.vm.t_to_g(ref_c, alt_ac=refseq_ac, alt_aln_method='splign')
                except:
                    e = sys.exc_info()[0]
                    ref_g_dict['error'] = e
                try:
                    vr.validate(ref_g_dict['ref_g'])
                except:
                    e = sys.exc_info()[0]
                    ref_g_dict['error'] = e
                if ref_g_dict['error'] == 'false':
                    return ref_g_dict
                else:
                    continue
        # Return as an error if all fail
        return ref_g_dict


    """
    Parses genomic variant strings into hgvs objects
    Maps genomic hgvs object into a coding hgvs object if the c accession string is provided
    returns a c. variant description string
    """


    def g_to_c(self, var_g, tx_ac, hpOld, evm):
        pat_g = re.compile("\:g\.")  # Pattern looks for :g.
        # If the :g. pattern is present in the input variant
        if pat_g.search(var_g):
            # convert the input string into a hgvs object by parsing
            var_g = self.hp.parse_hgvs_variant(var_g)
            # Map to coding variant
            var_c = str(evm.g_to_c(var_g, tx_ac))
            return var_c


    """
    Parses genomic variant strings into hgvs objects
    Maps genomic hgvs object into a non-coding hgvs object if the n accession string is provided
    returns a n. variant description string
    """


    def g_to_n(self, var_g, tx_ac, hpOld, evm):
        pat_g = re.compile("\:g\.")  # Pattern looks for :g.
        # If the :g. pattern is present in the input variant
        if pat_g.search(var_g):
            # convert the input string into a hgvs object by parsing
            var_g = self.hp.parse_hgvs_variant(var_g)
            # Map to coding variant
            var_n = str(evm.g_to_n(var_g, tx_ac))
            return var_n


    """
    Ensures variant strings are transcript c. or n.
    returns parsed hgvs c. or n. object
    """


    def coding(self, variant, hpOld):
        # If the :c. pattern is present in the input variant
        if re.search(':c.', variant) or re.search(':n.', variant):
            # convert the input string into a hgvs object
            var_c = self.hp.parse_hgvs_variant(variant)
            return var_c


    """
    Mapping transcript to genomic position
    Ensures variant strings are transcript c. or n.
    returns parsed hgvs g. object
    """


    def genomic(self, variant, evm, primary_assembly,hn):
        # Set regular expressions for if statements
        pat_g = re.compile("\:g\.")  # Pattern looks for :g.
        pat_n = re.compile("\:n\.")
        pat_c = re.compile("\:c\.")  # Pattern looks for :c.

        # If the :c. pattern is present in the input variant
        if pat_c.search(variant) or pat_n.search(variant):
            error = 'false'
            hgvs_var = self.hp.parse_hgvs_variant(variant)
            try:
                var_g = self.myevm_t_to_g(hgvs_var, evm, primary_assembly,hn)  # genomic level variant
            except hgvs.exceptions.HGVSError as e:
                error = e
            if error != 'false':
                var_g = 'error ' + str(e)
            return var_g

        # If the :g. pattern is present in the input variant
        elif (pat_g.search(variant)):  # or (pat_n.search(variant)):
            # convert the input string into a hgvs object
            var_g = self.hp.parse_hgvs_variant(variant)
            return var_g


    """
    Mapping transcript to protein prediction
    Ensures variant strings are transcript c.
    returns parsed hgvs p. object
    """




    """
    Function which takes a NORMALIZED hgvs Python transcript variant and maps to a specified protein reference sequence. A protein
    level hgvs python object is returned.
    
    Note the function currently assumes that the transcript description is correctly normalized having come from the 
    previous g_to_t function
    """





    """
    Ensures variant strings are g.
    returns parsed hgvs g. object
    """


    def hgvs_genomic(self, variant, hpOld):
        # Set regular expressions for if statements
        pat_g = re.compile("\:g\.")  # Pattern looks for :g. Note (gene) has been removed
        # If the :g. pattern is present in the input variant
        if pat_g.search(variant):
            # convert the input string into a hgvs object
            var_g = self.hp.parse_hgvs_variant(variant)
            return var_g


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


    def myevm_t_to_g(self,hgvs_c, no_norm_evm, primary_assembly, hn):

        # store the input
        stored_hgvs_c = copy.deepcopy(hgvs_c)
        expand_out = 'false'
        utilise_gap_code = True

        # Gap gene black list
        try:
            gene_symbol = self.db.get.get_gene_symbol_from_transcriptID(hgvs_c.ac)
        except Exception:
            utilise_gap_code = False
        else:
            # If the gene symbol is not in the list, the value False will be returned
            utilise_gap_code = vvChromosomes.gap_black_list(gene_symbol)
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
                        inv_alt = self.revcomp(hgvs_t.posedit.edit.ref)
                        t_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(
                            hgvs_t.posedit.pos.end.base) + 'del' + hgvs_t.posedit.edit.ref + 'ins' + inv_alt
                        hgvs_t_delins = self.hp.parse_hgvs_variant(t_delins)
                        pre_base = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                                hgvs_t.posedit.pos.start.base - 1)
                        post_base = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
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
                        hgvs_t = self.hp.parse_hgvs_variant(hgvs_str)
                    elif hgvs_c.posedit.edit.type == 'dup':
                        pre_base = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                                hgvs_t.posedit.pos.start.base - 1)
                        post_base = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
                                                 hgvs_t.posedit.pos.end.base + 1)
                        alt = pre_base + hgvs_t.posedit.edit.ref + hgvs_t.posedit.edit.ref + post_base
                        ref = pre_base + hgvs_t.posedit.edit.ref + post_base
                        dup_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(
                            hgvs_t.posedit.pos.start.base - 1) + '_' + str(
                            (hgvs_t.posedit.pos.start.base + len(ref)) - 2) + 'del' + ref + 'ins' + alt
                        hgvs_t = self.hp.parse_hgvs_variant(dup_to_delins)
                    elif hgvs_c.posedit.edit.type == 'ins':
                        ins_ref = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                               hgvs_t.posedit.pos.end.base + 1)
                        ins_alt = ins_ref[:2] + hgvs_t.posedit.edit.alt + ins_ref[-2:]
                        ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(
                            hgvs_t.posedit.pos.start.base - 1) + '_' + str(
                            hgvs_t.posedit.pos.end.base + 1) + 'del' + ins_ref + 'ins' + ins_alt
                        hgvs_t = self.hp.parse_hgvs_variant(ins_to_delins)
                    else:
                        if str(hgvs_t.posedit.edit.alt) == 'None':
                            hgvs_t.posedit.edit.alt = ''
                        pre_base = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                                hgvs_t.posedit.pos.start.base - 1)
                        post_base = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
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
                        hgvs_t = self.hp.parse_hgvs_variant(hgvs_str)
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
                hgvs_reform_ident = self.hp.parse_hgvs_variant(reform_ident)
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
            mapping_options = self.hdp.get_tx_mapping_options(hgvs_c.ac)

            if mapping_options == []:
                raise HGVSDataNotAvailableError(
                    "No alignment data between the specified transcript reference sequence and any GRCh37 and GRCh38 genomic reference sequences (including alternate chromosome assemblies, patches and RefSeqGenes) are available.")

            # Capture errors from attempted mappings
            attempted_mapping_error = ''

            for option in mapping_options:
                if re.match('blat', option[2]):
                    continue
                if re.match('NC_', option[1]):
                    chr_num = vvChromosomes.supported_for_mapping(str(option[1]), primary_assembly)
                    if chr_num != 'false':
                        try:
                            hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                        chr_num = vvChromosomes.supported_for_mapping(str(option[1]), primary_assembly)
                        if chr_num == 'false':
                            try:
                                hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                            chr_num = vvChromosomes.supported_for_mapping(str(option[1]), primary_assembly)
                            if chr_num != 'false':
                                try:
                                    hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                                chr_num = vvChromosomes.supported_for_mapping(str(option[1]),
                                                                                            primary_assembly)
                                if chr_num == 'false':
                                    try:
                                        hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                                    chr_num = vvChromosomes.supported_for_mapping(str(option[1]),
                                                                                                primary_assembly)
                                    if chr_num != 'false':
                                        try:
                                            hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                                        chr_num = vvChromosomes.supported_for_mapping(str(option[1]),
                                                                                                    primary_assembly)
                                        if chr_num == 'false':
                                            try:
                                                hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                                                hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                    ref = self.sf.fetch_seq(str(hgvs_genomic.ac), hgvs_genomic.posedit.pos.start.base - 1,
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
                stored_hgvs_n = self.vm.c_to_n(stored_hgvs_c)
            else:
                stored_hgvs_n = stored_hgvs_c
            stored_ref = self.sf.fetch_seq(str(stored_hgvs_n.ac), stored_hgvs_n.posedit.pos.start.base - 1,
                                      stored_hgvs_n.posedit.pos.end.base)
            stored_hgvs_c.posedit.edit.ref = stored_ref

        if (hgvs_genomic.posedit.edit.ref == '' or hgvs_genomic.posedit.edit.ref is None) and expand_out != 'false':
            if hgvs_genomic.posedit.edit.type == 'ins':
                stored_ref = self.sf.fetch_seq(str(hgvs_genomic.ac), hgvs_genomic.posedit.pos.start.base - 1,
                                          hgvs_genomic.posedit.pos.end.base)
                stored_alt = stored_ref[:1] + hgvs_genomic.posedit.edit.alt + stored_ref[-1:]
                hgvs_genomic.posedit.edit.ref = stored_ref
                hgvs_genomic.posedit.edit.alt = stored_alt

        # First look for variants mapping to the flanks of gaps
        # either in the gap or on the flank but not fully within the gap
        if expand_out == 'true':

            nr_genomic = self.nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)

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
                        genomic_gap_variant = self.vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)
                        try:
                            hn.normalize(genomic_gap_variant)
                        # Still a problem
                        except hgvs.exceptions.HGVSInvalidVariantError as e:
                            if 'base start position must be <= end position' in str(e) and \
                                    'Length implied by coordinates must equal' in error_type_1:
                                make_gen_var = copy.copy(nr_genomic)
                                make_gen_var.posedit.edit.ref = self.sf.fetch_seq(nr_genomic.ac,
                                                                             nr_genomic.posedit.pos.start.base - 1,
                                                                             nr_genomic.posedit.pos.end.base)
                                genomic_gap_variant = make_gen_var

                                error_type_1 = None
                        else:
                            genomic_gap_variant = self.nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)

                    if error_type_1 == 'base start position must be <= end position':
                        logger.warning('Variant is fully within a genomic gap')
                        genomic_gap_variant = self.vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)

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
                        transcript_gap_variant = self.vm.g_to_t(genomic_gap_variant, hgvs_c.ac)
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
                                transcript_gap_n = self.hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)
                                transcript_gap_alt_n_delins_from_dup = transcript_gap_alt_n.ac + ':' + transcript_gap_alt_n.type + '.' + str(
                                    transcript_gap_alt_n.posedit.pos.start.base) + '_' + str(
                                    transcript_gap_alt_n.posedit.pos.end.base) + 'del' + transcript_gap_alt_n.posedit.edit.ref + 'ins' + transcript_gap_alt_n.posedit.edit.ref + transcript_gap_alt_n.posedit.edit.ref
                                transcript_gap_alt_n = self.hp.parse_hgvs_variant(transcript_gap_alt_n_delins_from_dup)

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
                            transcript_gap_variant = self.vm.n_to_c(transcript_gap_n)
                        except:
                            transcript_gap_variant = transcript_gap_n

                        try:
                            hgvs_genomic = self.vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                            hgvs_genomic = hn.normalize(hgvs_genomic)
                        except Exception as e:
                            if str(e) == "base start position must be <= end position":
                                # Expansion out is required to map back to the genomic position
                                pre_base = self.sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.start.base - 2,
                                                        transcript_gap_n.posedit.pos.start.base - 1)
                                post_base = self.sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.end.base,
                                                         transcript_gap_n.posedit.pos.end.base + 1)
                                transcript_gap_n.posedit.pos.start.base = transcript_gap_n.posedit.pos.start.base - 1
                                transcript_gap_n.posedit.pos.end.base = transcript_gap_n.posedit.pos.end.base + 1
                                transcript_gap_n.posedit.edit.ref = pre_base + transcript_gap_n.posedit.edit.ref + post_base
                                transcript_gap_n.posedit.edit.alt = pre_base + transcript_gap_n.posedit.edit.alt + post_base
                                try:
                                    transcript_gap_variant = self.vm.n_to_c(transcript_gap_n)
                                except:
                                    transcript_gap_variant = transcript_gap_n
                                hgvs_genomic = self.vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
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
                genomic_gap_variant = self.vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)
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
                    transcript_gap_variant = self.vm.g_to_t(genomic_gap_variant, hgvs_c.ac)
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
                            transcript_gap_n = self.hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)
                            transcript_gap_alt_n_delins_from_dup = transcript_gap_alt_n.ac + ':' + transcript_gap_alt_n.type + '.' + str(
                                transcript_gap_alt_n.posedit.pos.start.base) + '_' + str(
                                transcript_gap_alt_n.posedit.pos.end.base) + 'del' + transcript_gap_alt_n.posedit.edit.ref + 'ins' + transcript_gap_alt_n.posedit.edit.ref + transcript_gap_alt_n.posedit.edit.ref
                            transcript_gap_alt_n = self.hp.parse_hgvs_variant(transcript_gap_alt_n_delins_from_dup)

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
                        transcript_gap_variant = self.vm.n_to_c(transcript_gap_n)
                    except:
                        transcript_gap_variant = transcript_gap_n

                    try:
                        hgvs_genomic = self.vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                        hgvs_genomic = hn.normalize(hgvs_genomic)
                    except Exception as e:
                        if str(e) == "base start position must be <= end position":
                            # Expansion out is required to map back to the genomic position
                            pre_base = self.sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.start.base - 2,
                                                    transcript_gap_n.posedit.pos.start.base - 1)
                            post_base = self.sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.end.base,
                                                     transcript_gap_n.posedit.pos.end.base + 1)
                            transcript_gap_n.posedit.pos.start.base = transcript_gap_n.posedit.pos.start.base - 1
                            transcript_gap_n.posedit.pos.end.base = transcript_gap_n.posedit.pos.end.base + 1
                            transcript_gap_n.posedit.edit.ref = pre_base + transcript_gap_n.posedit.edit.ref + post_base
                            transcript_gap_n.posedit.edit.alt = pre_base + transcript_gap_n.posedit.edit.alt + post_base
                            try:
                                transcript_gap_variant = self.vm.n_to_c(transcript_gap_n)
                            except:
                                transcript_gap_variant = transcript_gap_n
                            hgvs_genomic = self.vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                            hgvs_genomic = hn.normalize(hgvs_genomic)

        # Ins variants map badly - Especially between c. exon/exon boundary
        if hgvs_c.posedit.edit.type == 'ins' and hgvs_c.posedit.pos.start.offset == 0 and hgvs_c.posedit.pos.end.offset == 0:
            try:
                hn.normalize(hgvs_genomic)
            except hgvs.exceptions.HGVSError as e:
                error = str(e)
                if error == 'insertion length must be 1':
                    if hgvs_c.type == 'c':
                        hgvs_t = self.vm.c_to_n(hgvs_c)
                    else:
                        hgvs_t = copy.copy(hgvs_c)
                    ins_ref = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 1, hgvs_t.posedit.pos.end.base)
                    ins_alt = ins_ref[:1] + hgvs_t.posedit.edit.alt + ins_ref[-1:]
                    ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(
                        hgvs_t.posedit.pos.end.base) + 'del' + ins_ref + 'ins' + ins_alt
                    hgvs_t = self.hp.parse_hgvs_variant(ins_to_delins)
                    try:
                        hgvs_c = self.vm.n_to_c(hgvs_t)
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

    def noreplace_myevm_t_to_g(self,hgvs_c, evm, hdpOld, primary_assembly, vmOld, hn, hpOld, sfOld, no_norm_evm):
        try:
            hgvs_genomic = evm.t_to_g(hgvs_c)
            hn.normalize(hgvs_genomic)
        # This will fail on multiple refs for NC_
        except hgvs.exceptions.HGVSError as e:
            # Recover all available mapping options from UTA
            mapping_options = self.hdp.get_tx_mapping_options(hgvs_c.ac)

            if mapping_options == []:
                raise HGVSDataNotAvailableError("no g. mapping options available")

            for option in mapping_options:
                if re.match('blat', option[2]):
                    continue
                if re.match('NC_', option[1]):
                    chr_num = vvChromosomes.supported_for_mapping(str(option[1]), primary_assembly)
                    if chr_num != 'false':
                        try:
                            hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                        chr_num = vvChromosomes.supported_for_mapping(str(option[1]), primary_assembly)
                        if chr_num != 'false':
                            try:
                                hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                            chr_num = vvChromosomes.supported_for_mapping(str(option[1]), primary_assembly)
                            if chr_num == 'false':
                                try:
                                    hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                                chr_num = vvChromosomes.supported_for_mapping(str(option[1]),
                                                                                            primary_assembly)
                                if chr_num != 'false':
                                    try:
                                        hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                                    chr_num = vvChromosomes.supported_for_mapping(str(option[1]),
                                                                                                primary_assembly)
                                    if chr_num == 'false':
                                        try:
                                            hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                                        chr_num = vvChromosomes.supported_for_mapping(str(option[1]),
                                                                                                    primary_assembly)
                                        if chr_num != 'false':
                                            try:
                                                hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                                            chr_num = vvChromosomes.supported_for_mapping(str(option[1]),
                                                                                                        primary_assembly)
                                            if chr_num == 'false':
                                                try:
                                                    hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                                                    hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
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
                        hgvs_t = self.vm.c_to_n(hgvs_c)
                    else:
                        hgvs_t = copy.copy(hgvs_c)
                    ins_ref = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 1, hgvs_t.posedit.pos.end.base)
                    ins_alt = ins_ref[:1] + hgvs_t.posedit.edit.alt + ins_ref[-1:]
                    ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(
                        hgvs_t.posedit.pos.end.base) + 'del' + ins_ref + 'ins' + ins_alt
                    hgvs_t = self.hp.parse_hgvs_variant(ins_to_delins)
                    try:
                        hgvs_c = self.vm.n_to_c(hgvs_t)
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

    def myevm_g_to_t(self,evm, hgvs_genomic, alt_ac):
        hgvs_t = evm.g_to_t(hgvs_genomic, alt_ac)
        return hgvs_t
    def myvm_t_to_g(self, hgvs_c, alt_chr, no_norm_evm, hn):
        # store the input
        stored_hgvs_c = copy.deepcopy(hgvs_c)
        expand_out = 'false'
        utilise_gap_code = True

        # Gap gene black list
        try:
            gene_symbol = self.db.get.get_gene_symbol_from_transcriptID(hgvs_c.ac)
        except Exception:
            utilise_gap_code = False
        else:
            # If the gene symbol is not in the list, the value False will be returned
            utilise_gap_code = vvChromosomes.gap_black_list(gene_symbol)
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
                        inv_alt = self.revcomp(hgvs_t.posedit.edit.ref)
                        t_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(
                            hgvs_t.posedit.pos.end.base) + 'del' + hgvs_t.posedit.edit.ref + 'ins' + inv_alt
                        hgvs_t_delins = self.hp.parse_hgvs_variant(t_delins)
                        pre_base = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                                hgvs_t.posedit.pos.start.base - 1)
                        post_base = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
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
                        hgvs_t = self.hp.parse_hgvs_variant(hgvs_str)
                    if hgvs_c.posedit.edit.type == 'dup':
                        # hgvs_t = reverse_normalize.normalize(hgvs_t)
                        pre_base = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                                hgvs_t.posedit.pos.start.base - 1)
                        post_base = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
                                                 hgvs_t.posedit.pos.end.base + 1)
                        alt = pre_base + hgvs_t.posedit.edit.ref + hgvs_t.posedit.edit.ref + post_base
                        ref = pre_base + hgvs_t.posedit.edit.ref + post_base
                        dup_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(
                            hgvs_t.posedit.pos.start.base - 1) + '_' + str(
                            (hgvs_t.posedit.pos.start.base + len(ref)) - 2) + 'del' + ref + 'ins' + alt
                        hgvs_t = self.hp.parse_hgvs_variant(dup_to_delins)
                    elif hgvs_c.posedit.edit.type == 'ins':
                        ins_ref = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                               hgvs_t.posedit.pos.end.base + 1)
                        ins_alt = ins_ref[:2] + hgvs_t.posedit.edit.alt + ins_ref[-2:]
                        ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(
                            hgvs_t.posedit.pos.start.base - 1) + '_' + str(
                            hgvs_t.posedit.pos.end.base + 1) + 'del' + ins_ref + 'ins' + ins_alt
                        hgvs_t = self.hp.parse_hgvs_variant(ins_to_delins)
                    else:
                        if str(hgvs_t.posedit.edit.alt) == 'None':
                            hgvs_t.posedit.edit.alt = ''
                        pre_base = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                                hgvs_t.posedit.pos.start.base - 1)
                        post_base = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
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
                        hgvs_t = self.hp.parse_hgvs_variant(hgvs_str)
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
                hgvs_reform_ident = self.hp.parse_hgvs_variant(reform_ident)
                try:
                    hn.normalize(hgvs_reform_ident)
                except hgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if re.search('spanning the exon-intron boundary', error) or re.search(
                            'Normalization of intronic variants', error):
                        hgvs_c = copy.deepcopy(stored_hgvs_c)

        hgvs_genomic = self.vm.t_to_g(hgvs_c, alt_chr)
        if hgvs_c.posedit.edit.type == 'identity' and hgvs_genomic.posedit.edit.type == 'delins' and hgvs_genomic.posedit.edit.alt == '' and expand_out != 'true':
            hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.ref
        if hgvs_genomic.posedit.edit.type == 'ins' and utilise_gap_code is True:
            try:
                hgvs_genomic = hn.normalize(hgvs_genomic)
            except hgvs.exceptions.HGVSError as e:
                error = str(e)
                if error == 'insertion length must be 1':
                    ref = self.sf.fetch_seq(str(hgvs_genomic.ac), hgvs_genomic.posedit.pos.start.base - 1,
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
                stored_hgvs_n = self.vm.c_to_n(stored_hgvs_c)
            else:
                stored_hgvs_n = stored_hgvs_c
            stored_ref = self.sf.fetch_seq(str(stored_hgvs_n.ac), stored_hgvs_n.posedit.pos.start.base - 1,
                                      stored_hgvs_n.posedit.pos.end.base)
            stored_hgvs_c.posedit.edit.ref = stored_ref

        if (hgvs_genomic.posedit.edit.ref == '' or hgvs_genomic.posedit.edit.ref is None) and expand_out != 'false':
            if hgvs_genomic.posedit.edit.type == 'ins':
                stored_ref = self.sf.fetch_seq(str(hgvs_genomic.ac), hgvs_genomic.posedit.pos.start.base - 1,
                                          hgvs_genomic.posedit.pos.end.base)
                stored_alt = stored_ref[:1] + hgvs_genomic.posedit.edit.alt + stored_ref[-1:]
                hgvs_genomic.posedit.edit.ref = stored_ref
                hgvs_genomic.posedit.edit.alt = stored_alt

        # First look for variants mapping to the flanks of gaps
        # either in the gap or on the flank but not fully within the gap
        if expand_out == 'true':
            nr_genomic = self.nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)
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
                        genomic_gap_variant = self.vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)
                        try:
                            hn.normalize(genomic_gap_variant)
                        # Still a problem
                        except hgvs.exceptions.HGVSInvalidVariantError as e:
                            if 'base start position must be <= end position' in str(e) and \
                                    'Length implied by coordinates must equal' in error_type_1:
                                make_gen_var = copy.copy(nr_genomic)
                                make_gen_var.posedit.edit.ref = self.sf.fetch_seq(nr_genomic.ac,
                                                                             nr_genomic.posedit.pos.start.base - 1,
                                                                             nr_genomic.posedit.pos.end.base)
                                genomic_gap_variant = make_gen_var
                                error_type_1 = None
                        else:
                            genomic_gap_variant = self.nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)

                    if error_type_1 == 'base start position must be <= end position':
                        logger.warning('Variant is fully within a genomic gap')
                        genomic_gap_variant = self.vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)

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
                        transcript_gap_variant = self.vm.g_to_t(genomic_gap_variant, hgvs_c.ac)
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
                                transcript_gap_n = self.hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)
                                transcript_gap_alt_n_delins_from_dup = transcript_gap_alt_n.ac + ':' + transcript_gap_alt_n.type + '.' + str(
                                    transcript_gap_alt_n.posedit.pos.start.base) + '_' + str(
                                    transcript_gap_alt_n.posedit.pos.end.base) + 'del' + transcript_gap_alt_n.posedit.edit.ref + 'ins' + transcript_gap_alt_n.posedit.edit.ref + transcript_gap_alt_n.posedit.edit.ref
                                transcript_gap_alt_n = self.hp.parse_hgvs_variant(transcript_gap_alt_n_delins_from_dup)

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
                            transcript_gap_variant = self.vm.n_to_c(transcript_gap_n)
                        except:
                            transcript_gap_variant = transcript_gap_n

                        try:
                            hgvs_genomic = self.vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                            hgvs_genomic = hn.normalize(hgvs_genomic)
                        except Exception as e:
                            if str(e) == "base start position must be <= end position":
                                # Expansion out is required to map back to the genomic position
                                pre_base = self.sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.start.base - 2,
                                                        transcript_gap_n.posedit.pos.start.base - 1)
                                post_base = self.sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.end.base,
                                                         transcript_gap_n.posedit.pos.end.base + 1)
                                transcript_gap_n.posedit.pos.start.base = transcript_gap_n.posedit.pos.start.base - 1
                                transcript_gap_n.posedit.pos.end.base = transcript_gap_n.posedit.pos.end.base + 1
                                transcript_gap_n.posedit.edit.ref = pre_base + transcript_gap_n.posedit.edit.ref + post_base
                                transcript_gap_n.posedit.edit.alt = pre_base + transcript_gap_n.posedit.edit.alt + post_base
                                try:
                                    transcript_gap_variant = self.vm.n_to_c(transcript_gap_n)
                                except:
                                    transcript_gap_variant = transcript_gap_n
                                hgvs_genomic = self.vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
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
                genomic_gap_variant = self.self.vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)
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
                    transcript_gap_variant = self.vm.g_to_t(genomic_gap_variant, hgvs_c.ac)
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
                            transcript_gap_n = self.hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)
                            transcript_gap_alt_n_delins_from_dup = transcript_gap_alt_n.ac + ':' + transcript_gap_alt_n.type + '.' + str(
                                transcript_gap_alt_n.posedit.pos.start.base) + '_' + str(
                                transcript_gap_alt_n.posedit.pos.end.base) + 'del' + transcript_gap_alt_n.posedit.edit.ref + 'ins' + transcript_gap_alt_n.posedit.edit.ref + transcript_gap_alt_n.posedit.edit.ref
                            transcript_gap_alt_n = self.hp.parse_hgvs_variant(transcript_gap_alt_n_delins_from_dup)

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
                        transcript_gap_variant = self.vm.n_to_c(transcript_gap_n)
                    except:
                        transcript_gap_variant = transcript_gap_n

                    try:
                        hgvs_genomic = self.vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                        hgvs_genomic = hn.normalize(hgvs_genomic)
                    except Exception as e:
                        if str(e) == "base start position must be <= end position":
                            # Expansion out is required to map back to the genomic position
                            pre_base = self.sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.start.base - 2,
                                                    transcript_gap_n.posedit.pos.start.base - 1)
                            post_base = self.sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.end.base,
                                                     transcript_gap_n.posedit.pos.end.base + 1)
                            transcript_gap_n.posedit.pos.start.base = transcript_gap_n.posedit.pos.start.base - 1
                            transcript_gap_n.posedit.pos.end.base = transcript_gap_n.posedit.pos.end.base + 1
                            transcript_gap_n.posedit.edit.ref = pre_base + transcript_gap_n.posedit.edit.ref + post_base
                            transcript_gap_n.posedit.edit.alt = pre_base + transcript_gap_n.posedit.edit.alt + post_base
                            try:
                                transcript_gap_variant = self.vm.n_to_c(transcript_gap_n)
                            except:
                                transcript_gap_variant = transcript_gap_n
                            hgvs_genomic = self.vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                            hgvs_genomic = hn.normalize(hgvs_genomic)

        # Ins variants map badly - Especially between c. exon/exon boundary
        if hgvs_c.posedit.edit.type == 'ins' and hgvs_c.posedit.pos.start.offset == 0 and hgvs_c.posedit.pos.end.offset == 0:
            try:
                hn.normalize(hgvs_genomic)
            except hgvs.exceptions.HGVSError as e:
                error = str(e)
                if error == 'insertion length must be 1':
                    if hgvs_c.type == 'c':
                        hgvs_t = self.vm.c_to_n(hgvs_c)
                    else:
                        hgvs_t = copy.copy(hgvs_c)
                    ins_ref = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 1, hgvs_t.posedit.pos.end.base)
                    ins_alt = ins_ref[:1] + hgvs_t.posedit.edit.alt + ins_ref[-1:]
                    ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(
                        hgvs_t.posedit.pos.end.base) + 'del' + ins_ref + 'ins' + ins_alt
                    hgvs_t = self.hp.parse_hgvs_variant(ins_to_delins)
                    try:
                        hgvs_c = self.vm.n_to_c(hgvs_t)
                    except Exception:
                        hgvs_c = copy.copy(hgvs_t)
                    try:
                        hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)
                    except Exception as e:
                        error = str(e)
                        logger.warning('Ins mapping error in myt_to_g ' + error)

        return hgvs_genomic

    """
    parse p. strings into hgvs p. objects
    """


    def hgvs_protein(self, variant, hpOld):
        # Set regular expressions for if statements
        pat_p = re.compile("\:p\.")  # Pattern looks for :g. Note (gene) has been removed
        # If the :p. pattern is present in the input variant
        if pat_p.search(variant):
            # convert the input string into a hgvs object
            var_p = self.hp.parse_hgvs_variant(variant)
            return var_p


    """
    Convert r. into c.
    """


    def hgvs_r_to_c(self, hgvs_object):
        # check for LRG_t with r.
        if re.match('LRG', hgvs_object.ac):
            transcript_ac = self.db.get.get_RefSeqTranscriptID_from_lrgTranscriptID(hgvs_object.ac)
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


    def hgvs_c_to_r(self, hgvs_object):
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
    """


    def tx_identity_info(self, variant, hdpOld):
        # Set regular expressions for if statements
        pat_c = re.compile("\:c\.")  # Pattern looks for :c. Note (gene) has been removed
        pat_n = re.compile("\:n\.")  # Pattern looks for :c. Note (gene) has been removed
        pat_r = re.compile("\:r\.")  # Pattern looks for :c. Note (gene) has been removed

        # If the :c. pattern is present in the input variant
        if pat_c.search(variant):
            # Remove all text to the right and including pat_c
            tx_ac = variant[:variant.index(':c.') + len(':c.')]
            tx_ac = pat_c.sub('', tx_ac)
            # Interface with the UTA database via get_tx_identity in uta.py
            tx_id_info = self.hdp.get_tx_identity_info(tx_ac)
            # NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
            return tx_id_info

        # If the :n. pattern is present in the input variant
        if pat_n.search(variant):
            # Remove all text to the right and including pat_c
            tx_ac = variant[:variant.index(':n.') + len(':n.')]
            tx_ac = pat_n.sub('', tx_ac)
            # Interface with the UTA database via get_tx_identity in uta.py
            tx_id_info = self.hdp.get_tx_identity_info(tx_ac)
            # NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
            return tx_id_info

        # If the :r. pattern is present in the input variant
        if pat_r.search(variant):
            # Remove all text to the right and including pat_c
            tx_ac = variant[:variant.index(':r.') + len(':r.')]
            tx_ac = pat_r.sub('', tx_ac)
            # Interface with the UTA database via get_tx_identity in uta.py
            tx_id_info = self.hdp.get_tx_identity_info(tx_ac)
            # NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
            return tx_id_info


    """
    Input c. r. nd accession string
    Use uta.py (hdp) to return the identity information for the transcript variant 
    see hgvs.dataproviders.uta.py for details
    """


    def tx_id_info(self, alt_ac, hdpOld):
        tx_id_info = self.hdp.get_tx_identity_info(alt_ac)
        # NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
        return tx_id_info


    """
    Use uta.py (hdp) to return the transcript information for a specified gene (HGNC SYMBOL)
    see hgvs.dataproviders.uta.py for details
    """


    def tx_for_gene(self, hgnc, hdpOld):
        # Interface with the UTA database via get_tx_for_gene in uta.py
        tx_for_gene = self.hdp.get_tx_for_gene(hgnc)
        return tx_for_gene


    """
    Extract RefSeqGene Accession from transcript information
    see hgvs.dataproviders.uta.py for details
    """


    def ng_extract(self, tx_for_gene):
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
    Returns exon information for a given transcript
    e.g. how the exons align to the genomic reference
    see hgvs.dataproviders.uta.py for details
    """


    def tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        # Interface with the UTA database via get_tx_exons in uta.py
        try:
            tx_exons = self.hdp.get_tx_exons(tx_ac, alt_ac, alt_aln_method)
        except hgvs.exceptions.HGVSError as e:
            #e
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
            return tx_exons
        else:
            return tx_exons


    """
    Automatically maps genomic positions onto all overlapping transcripts
    """


    def relevant_transcripts(self, hgvs_genomic, evm, alt_aln_method,reverse_normalizer):
        # Pass relevant transcripts for the input variant to rts
        # Note, the evm method misses one end, the hdp. method misses the other. Combine both
        rts_list = self.hdp.get_tx_for_region(hgvs_genomic.ac, alt_aln_method, hgvs_genomic.posedit.pos.start.base-1, hgvs_genomic.posedit.pos.end.base-1)
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
                    tx_exons = self.hdp.get_tx_exons(tx_ac, alt_ac, alt_aln_method)
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
                    rev_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic)
                    # map back to coding
                    variant = evm.g_to_t(rev_hgvs_genomic, tx_ac)
            code_var.append(str(variant))
        return code_var


    """
    Take HGVS string, parse into hgvs object and validate
    """


    def validateHGVS(self, input):
        hgvs_input = self.hp.parse_hgvs_variant(input)
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
            self.vr.validate(hgvs_input)
        except hgvs.exceptions.HGVSError as e:

            error = e
            return error

        else:
            error = 'false'
            return error

    """
    Search HGNC rest
    """


    def hgnc_rest(self, path):
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


    def entrez_efetch(self, db, id, rettype, retmode):
        # IMPORT Bio modules
        # from Bio import Entrez
        Entrez.email = self.entrezID
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


    def entrez_read(self,db, id, retmode):
        # IMPORT Bio modules
        # from Bio import Entrez
        Entrez.email = self.entrezID
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


    def revcomp(self, bases):
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


    def merge_hgvs_3pr(self, hgvs_variant_list,hn):
        # Ensure c. is mapped to the
        h_list = []

        # Sanity check and format the submitted variants
        for hgvs_v in hgvs_variant_list:
            # For testing include parser
            try:
                hgvs_v = self.hp.parse_hgvs_variant(hgvs_v)
            except Exception as e:
                print e
                pass

            # Validate
            self.vr.validate(hgvs_v)  # Let hgvs errors deal with invalid variants and not hgvs objects
            if hgvs_v.type == 'c':
                try:
                    hgvs_v = self.vm.c_to_n(hgvs_v)
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
                    ins_seq = self.sf.fetch_seq(hgvs_v.ac, merge_end_pos, hgvs_v.posedit.pos.start.base - 1)
                    gapping = hgvs_v.ac + ':' + hgvs_v.type + '.' + str(merge_end_pos + 1) + '_' + str(
                        hgvs_v.posedit.pos.start.base - 1) + 'delins' + ins_seq
                    hgvs_gapping = self.hp.parse_hgvs_variant(gapping)
                    full_list.append(hgvs_gapping)
                    # update end_pos
                    merge_end_pos = hgvs_v.posedit.pos.end.base
                    # Append to the final list of variants
                    full_list.append(hgvs_v)

        # Generate the alt sequence
        alt_sequence = ''
        for hgvs_v in full_list:
            ref_alt = vvHGVS.hgvs_ref_alt(hgvs_v)
            alt_sequence = alt_sequence + ref_alt['alt']

        # Fetch the reference sequence and copy it for the basis of the alt sequence
        reference_sequence = self.sf.fetch_seq(accession, merge_start_pos - 1, merge_end_pos)
        # Generate an hgvs_delins
        if alt_sequence == '':
            delins = accession + ':' + type + '.' + str(merge_start_pos) + '_' + str(
                merge_end_pos) + 'del' + reference_sequence
        else:
            delins = accession + ':' + type + '.' + str(merge_start_pos) + '_' + str(
                merge_end_pos) + 'del' + reference_sequence + 'ins' + alt_sequence
        hgvs_delins = self.hp.parse_hgvs_variant(delins)
        try:
            hgvs_delins = self.vm.n_to_c(hgvs_delins)
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


    def merge_hgvs_5pr(self, hgvs_variant_list):
        # Ensure c. is mapped to the
        h_list = []

        # Sanity check and format the submitted variants
        for hgvs_v in hgvs_variant_list:
            # For testing include parser
            try:
                hgvs_v = self.hp.parse_hgvs_variant(hgvs_v)
            except:
                pass

            # Validate
            self.vr.validate(hgvs_v)  # Let hgvs errors deal with invalid variants and not hgvs objects
            if hgvs_v.type == 'c':
                try:
                    hgvs_v = self.vm.c_to_n(hgvs_v)
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
            hgvs_v = self.reverse_hn.normalize(hgvs_v)

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
                    ins_seq = self.sf.fetch_seq(hgvs_v.ac, merge_end_pos, hgvs_v.posedit.pos.start.base - 1)
                    gapping = hgvs_v.ac + ':' + hgvs_v.type + '.' + str(merge_end_pos + 1) + '_' + str(
                        hgvs_v.posedit.pos.start.base - 1) + 'delins' + ins_seq
                    hgvs_gapping = self.hp.parse_hgvs_variant(gapping)
                    full_list.append(hgvs_gapping)
                    # update end_pos
                    merge_end_pos = hgvs_v.posedit.pos.end.base
                    # Append to the final list of variants
                    full_list.append(hgvs_v)

        # Generate the alt sequence
        alt_sequence = ''
        for hgvs_v in full_list:
            ref_alt = vvHGVS.hgvs_ref_alt(hgvs_v)
            alt_sequence = alt_sequence + ref_alt['alt']

        # Fetch the reference sequence and copy it for the basis of the alt sequence
        reference_sequence = self.sf.fetch_seq(accession, merge_start_pos - 1, merge_end_pos)

        # Generate an hgvs_delins
        if alt_sequence == '':
            delins = accession + ':' + type + '.' + str(merge_start_pos) + '_' + str(
                merge_end_pos) + 'del' + reference_sequence
        else:
            delins = accession + ':' + type + '.' + str(merge_start_pos) + '_' + str(
                merge_end_pos) + 'del' + reference_sequence + 'ins' + alt_sequence
        hgvs_delins = self.hp.parse_hgvs_variant(delins)
        try:
            hgvs_delins = self.vm.n_to_c(hgvs_delins)
        except:
            pass
        # Normalize (allow variants crossing into different exons)
        try:
            hgvs_delins = self.reverse_hn.normalize(hgvs_delins)
        except HGVSUnsupportedOperationError:
            pass
        return hgvs_delins


    """
    Function designed to merge multiple pseudo VCF variants (strings) into a single HGVS delins 
    using 5 prime normalization then return a 3 prime normalized final HGVS object
    """


    def merge_pseudo_vcf(self, vcf_list, genome_build, hn):
        hgvs_list = []
        # Convert pseudo_vcf list into a HGVS list
        for call in vcf_list:
            x55hgvs = vvHGVS.pvcf_to_hgvs(call, genome_build, normalization_direction=5, validator=self)
            hgvs_list.append(x55hgvs)
        # Merge
        hgvs_delins = self.merge_hgvs_5pr(hgvs_list)
        # normalize 3 prime
        hgvs_delins = hn.normalize(hgvs_delins)
        # return
        return hgvs_delins


    """
    HGVS allele handling function which takes a single HGVS allele description and 
    separates each allele into a list of HGVS variants
    """


    def hgvs_alleles(self, variant_description):
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
                    # now separate out the variants in each allele |
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
                        allele = str(self.merge_hgvs_3pr(each_allele))
                        merge.append(allele)
                        merged_alleles.append(merge)
                    my_alleles = merged_alleles

                elif re.search('\(;\)', remainder):
                    # If statement for uncertainties
                    # NM_004006.2:c.[296T>G;476C>T];[476C>T](;)1083A>C
                    if re.search('\[', remainder):
                        raise alleleVariantError('Unsupported format ' + type + '.' + remainder)
                    # NM_004006.2:c.2376G>C(;)3103del
                    # NM_000548.3:c.3623_3647del(;)3745_3756dup
                    alleles = remainder.split('(;)')
                    # now separate out the variants in each allele |
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
                    # now separate out the variants in each allele |
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
                        print each_allele
                        if re.search('\?', str(each_allele)):
                            # NM_004006.2:c.[2376G>C];[?]
                            continue
                        merge = []
                        allele = str(self.merge_hgvs_3pr(each_allele))
                        merge.append(allele)
                        merged_alleles.append(merge)
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

    # Covert chromosomal HGVS description to RefSeqGene
    def chr_to_rsg(self, hgvs_genomic, hn, vrOld):
        # print 'chr_to_rsg triggered'
        hgvs_genomic = hn.normalize(hgvs_genomic)
        # split the description
        # Accessions
        chr_ac = hgvs_genomic.ac
        # Positions
        chr_start_pos = int(hgvs_genomic.posedit.pos.start.base)
        chr_end_pos = int(hgvs_genomic.posedit.pos.end.base)
        # edit
        chr_edit = hgvs_genomic.posedit.edit

        # Pre set variable, note there could be several
        rsg_data_set = []

        # Recover table from MySql
        all_info = self.db.get.get_g_to_g_info()
        for line in all_info:
            # Logic to identify the correct RefSeqGene
            rsg_data = {}
            if chr_ac == line[1] and chr_start_pos >= int(line[2]) and chr_end_pos <= int(line[3]):
                # query = "SELECT refSeqGeneID, refSeqChromosomeID, startPos, endPos, orientation, hgncSymbol FROM refSeqGene_loci"
                # (u'NG_034189.1', u'NC_000004.12', 190173122, 190177845, u'+', u'DUX4L1')
                # Set the values of the data dictionary
                rsg_data['rsg_ac'] = line[0]
                rsg_data['chr_ac'] = line[1]
                rsg_data['rsg_start'] = line[2]
                rsg_data['rsg_end'] = line[3]
                rsg_data['ori'] = line[4]
                rsg_data['gene'] = line[5]
                rsg_data_set.append(rsg_data)
            else:
                continue

        # Compile descriptions and validate
        descriptions = []
        for rsg_data in rsg_data_set:
            rsg_ac = rsg_data['rsg_ac']
            rsg_start = rsg_data['rsg_start']
            rsg_end = rsg_data['rsg_end']
            ori = rsg_data['ori']
            gene = rsg_data['gene']
            # String the description
            if ori == '+':
                rsg_description = rsg_ac + ':g.' + str(chr_start_pos - int(rsg_start) + 1) + '_' + str(
                    chr_end_pos - int(rsg_start) + 1) + str(chr_edit)
                hgvs_refseqgene = self.hp.parse_hgvs_variant(rsg_description)
                try:
                    hgvs_refseqgene = hn.normalize(hgvs_refseqgene)
                except:
                    error = 'Not in SeqRepo'
                    data = {'hgvs_refseqgene': str(hgvs_refseqgene), 'gene': gene, 'valid': str(error)}
                    descriptions.append(data)
                    continue
                try:
                    self.vr.validate(hgvs_refseqgene)
                except hgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if re.search('does not agree with reference sequence', error):
                        match = re.findall('\(([GATC]+)\)', error)
                        new_ref = match[1]
                        hgvs_refseqgene.posedit.edit.ref = new_ref
                        error = 'true'
                    else:
                        pass
                    data = {'hgvs_refseqgene': str(hgvs_refseqgene), 'gene': gene, 'valid': str(error)}
                else:
                    data = {'hgvs_refseqgene': str(hgvs_refseqgene), 'gene': gene, 'valid': 'true'}
                descriptions.append(data)
            if ori == '-':
                # Reverse complement of bases may be required. Let normalizer do the lifting for strings of bases
                # Look for scenarios with RC needed bases and extract the bases from the edit
                if re.search(r"((del[GATCUgatcu]+))", str(chr_edit)):
                    bases = re.search(r"((del[GATCUgatcu]+))", str(chr_edit))
                    bases = bases.group(1)
                    chr_edit = 'del' + str(chr_edit).replace(bases, '')
                if re.search(r"((ins[GATCUgatcu]+))", str(chr_edit)):
                    bases = re.search(r"((ins[GATCUgatcu]+))", str(chr_edit))
                    bases = bases.group(1)
                    ins_revcomp = self.revcomp(bases)
                    chr_edit = str(chr_edit).replace(bases, '') + 'ins' + ins_revcomp
                if re.search(r"((dup[GATCUgatcu]+))", str(chr_edit)):
                    bases = re.search(r"((dup[GATCUgatcu]+))", str(chr_edit))
                    bases = bases.group(1)
                    chr_edit = 'dup' + str(chr_edit).replace(bases, '')
                if re.search(r"((inv[GATCUgatcu]+))", str(chr_edit)):
                    bases = re.search(r"((inv[GATCUgatcu]+))", str(chr_edit))
                    bases = bases.group(1)
                    chr_edit = 'inv' + str(chr_edit).replace(bases, '')
                if re.search('>', str(chr_edit)) or re.search('=', str(chr_edit)):
                    chr_edit = str(chr_edit)
                    chr_edit = chr_edit.replace('A>', 't>')
                    chr_edit = chr_edit.replace('T>', 'a>')
                    chr_edit = chr_edit.replace('G>', 'c>')
                    chr_edit = chr_edit.replace('C>', 'g>')
                    chr_edit = chr_edit.replace('>A', '>t')
                    chr_edit = chr_edit.replace('>T', '>a')
                    chr_edit = chr_edit.replace('>G', '>c')
                    chr_edit = chr_edit.replace('>C', '>g')
                    chr_edit = chr_edit.replace('C=', 'g=')
                    chr_edit = chr_edit.replace('G=', 'c=')
                    chr_edit = chr_edit.replace('A=', 't=')
                    chr_edit = chr_edit.replace('T=', 'a=')
                    chr_edit = chr_edit.upper()

                rsg_description = rsg_ac + ':g.' + str(
                    (int(rsg_end) - int(rsg_start)) - (chr_end_pos - int(rsg_start)) + 1) + '_' + str(
                    (int(rsg_end) - int(rsg_start)) - (chr_start_pos - int(rsg_start)) + 1) + str(chr_edit)
                hgvs_refseqgene = self.hp.parse_hgvs_variant(rsg_description)
                try:
                    hgvs_refseqgene = hn.normalize(hgvs_refseqgene)
                except:
                    error = 'Not in SeqRepo'
                    data = {'hgvs_refseqgene': str(hgvs_refseqgene), 'gene': gene, 'valid': str(error)}
                    descriptions.append(data)
                    continue
                try:
                    self.vr.validate(hgvs_refseqgene)
                except hgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if re.search('does not agree with reference sequence', error):
                        match = re.findall('\(([GATC]+)\)', error)
                        new_ref = match[1]
                        hgvs_refseqgene.posedit.edit.ref = new_ref
                        error = 'true'
                    else:
                        pass
                    data = {'hgvs_refseqgene': str(hgvs_refseqgene), 'gene': gene, 'valid': str(error)}
                else:
                    data = {'hgvs_refseqgene': str(hgvs_refseqgene), 'gene': gene, 'valid': 'true'}
                descriptions.append(data)

        # Return the required data. This is a dictionary containing the rsg description, validation status and gene ID
        return descriptions


    # Covert RefSeqGene HGVS description to Chromosomal
    def rsg_to_chr(self, hgvs_refseqgene, primary_assembly, hn, vr):
        # normalize
        try:
            hgvs_refseqgene = hn.normalize(hgvs_refseqgene)
        except:
            pass
        # split the description
        # Accessions
        rsg_ac = hgvs_refseqgene.ac
        # Positions
        rsg_start_pos = int(hgvs_refseqgene.posedit.pos.start.base)
        rsg_end_pos = int(hgvs_refseqgene.posedit.pos.end.base)
        # edit
        rsg_edit = hgvs_refseqgene.posedit.edit

        # Pre set variable, note there could be several
        chr_data_set = []

        # Recover table from MySql
        all_info = self.db.get.get_g_to_g_info()
        for line in all_info:
            # Logic to identify the correct RefSeqGene
            chr_data = {}
            if rsg_ac == line[0] and primary_assembly == line[6]:
                # query = "SELECT refSeqGeneID, refSeqChromosomeID, startPos, endPos, orientation, hgncSymbol FROM refSeqGene_loci"
                # (u'NG_034189.1', u'NC_000004.12', 190173122, 190177845, u'+', u'DUX4L1')
                # Set the values of the data dictionary
                chr_data['rsg_ac'] = line[0]
                chr_data['chr_ac'] = line[1]
                chr_data['rsg_start'] = line[2]
                chr_data['rsg_end'] = line[3]
                chr_data['ori'] = line[4]
                chr_data['gene'] = line[5]
                chr_data_set.append(chr_data)
            else:
                continue

        # Compile descriptions and validate
        descriptions = []
        for chr_data in chr_data_set:
            chr_ac = chr_data['chr_ac']
            rsg_ac = chr_data['rsg_ac']
            chr_start = int(chr_data['rsg_start'])
            chr_end = int(chr_data['rsg_end'])
            ori = chr_data['ori']
            gene = chr_data['gene']
            # String the description
            if ori == '+':
                chr_description = chr_ac + ':g.' + str(chr_start + rsg_start_pos - 1) + '_' + str(
                    chr_start + rsg_end_pos - 1) + str(rsg_edit)
                hgvs_genomic = self.hp.parse_hgvs_variant(chr_description)
                hgvs_genomic = hn.normalize(hgvs_genomic)
                try:
                    vr.validate(hgvs_genomic)
                except hgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if re.search('does not agree with reference sequence', error):
                        match = re.findall('\(([GATC]+)\)', error)
                        new_ref = match[1]
                        hgvs_genomic.posedit.edit.ref = new_ref
                        error = 'true'
                    else:
                        pass
                    # # print str(e) + '\n3.'
                    data = {'hgvs_genomic': str(hgvs_genomic), 'gene': gene, 'valid': str(error)}
                else:
                    data = {'hgvs_genomic': str(hgvs_genomic), 'gene': gene, 'valid': 'true'}
                descriptions.append(data)
            if ori == '-':
                # Reverse complement of bases may be required. Let normalizer do the lifting for strings of bases
                # Look for scenarios with RC needed bases and extract the bases from the edit
                if re.search(r"((del[GATCUgatcu]+))", str(rsg_edit)):
                    bases = re.search(r"((del[GATCUgatcu]+))", str(rsg_edit))
                    bases = bases.group(1)
                    rsg_edit = 'del' + str(rsg_edit).replace(bases, '')
                if re.search(r"((ins[GATCUgatcu]+))", str(rsg_edit)):
                    bases = re.search(r"((ins[GATCUgatcu]+))", str(rsg_edit))
                    bases = bases.group(1)
                    ins_revcomp = self.revcomp(bases)
                    rsg_edit = str(rsg_edit).replace(bases, '') + 'ins' + ins_revcomp
                if re.search(r"((dup[GATCUgatcu]+))", str(rsg_edit)):
                    bases = re.search(r"((dup[GATCUgatcu]+))", str(rsg_edit))
                    bases = bases.group(1)
                    rsg_edit = 'dup' + str(rsg_edit).replace(bases, '')
                if re.search(r"((inv[GATCUgatcu]+))", str(rsg_edit)):
                    bases = re.search(r"((inv[GATCUgatcu]+))", str(rsg_edit))
                    bases = bases.group(1)
                    rsg_edit = 'inv' + str(rsg_edit).replace(bases, '')
                if re.search('>', str(rsg_edit)) or re.search('=', str(rsg_edit)):
                    rsg_edit = str(rsg_edit)
                    rsg_edit = rsg_edit.replace('A>', 't>')
                    rsg_edit = rsg_edit.replace('T>', 'a>')
                    rsg_edit = rsg_edit.replace('G>', 'c>')
                    rsg_edit = rsg_edit.replace('C>', 'g>')
                    rsg_edit = rsg_edit.replace('>A', '>t')
                    rsg_edit = rsg_edit.replace('>T', '>a')
                    rsg_edit = rsg_edit.replace('>G', '>c')
                    rsg_edit = rsg_edit.replace('>C', '>g')
                    rsg_edit = rsg_edit.replace('C=', 'g=')
                    rsg_edit = rsg_edit.replace('G=', 'c=')
                    rsg_edit = rsg_edit.replace('A=', 't=')
                    rsg_edit = rsg_edit.replace('T=', 'a=')
                    rsg_edit = rsg_edit.upper()

                chr_description = chr_ac + ':g.' + str(
                    int(chr_start) + (int(chr_end) - int(chr_start)) - rsg_end_pos + 1) + '_' + str(
                    int(chr_start) + (int(chr_end) - int(chr_start)) - rsg_start_pos + 1) + str(rsg_edit)

                hgvs_genomic = self.hp.parse_hgvs_variant(chr_description)
                hgvs_genomic = hn.normalize(hgvs_genomic)
                try:
                    vr.validate(hgvs_genomic)
                except hgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if re.search('does not agree with reference sequence', error):
                        match = re.findall('\(([GATC]+)\)', error)
                        new_ref = match[1]
                        hgvs_genomic.posedit.edit.ref = new_ref
                        error = 'true'
                    data = {'hgvs_genomic': str(hgvs_genomic), 'gene': gene, 'valid': str(error)}
                else:
                    data = {'hgvs_genomic': str(hgvs_genomic), 'gene': gene, 'valid': 'true'}
                descriptions.append(data)

        # Return the required data. This is a dictionary containing the rsg description, validation status and gene ID
        return descriptions
