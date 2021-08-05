import re
import copy
import logging
import vvhgvs
import vvhgvs.validator
from . import vvMixinInit
from . import seq_data
from . import hgvs_utils
from Bio import Entrez, SeqIO
from . import utils as fn
import sys
import traceback

from vvhgvs.exceptions import HGVSError, HGVSDataNotAvailableError, HGVSUnsupportedOperationError, \
     HGVSInvalidVariantError

logger = logging.getLogger(__name__)


class Mixin(vvMixinInit.Mixin):
    """
    This mixin contains converters that use the validator's configuration information.
    It inherits the Init mixin
    """
    # def r_to_c(self, variant, evm):
    #     """
    #     r_to_c
    #     parses r. variant strings into hgvs object and maps to the c. equivalent.
    #     """
    #     # convert the input string into a hgvs object by parsing
    #     var_r = self.hp.parse_hgvs_variant(variant)
    #     # map to the coding sequence
    #     var_c = evm.r_to_c(var_r)  # coding level variant
    #     variant = str(var_c)
    #     c_from_r = {'variant': variant, 'type': ':c.'}
    #     return c_from_r
    #
    # def refseq(self, variant, refseq_ac, evm, primary_assembly):
    #     """
    #     Maps transcript variant descriptions onto specified RefSeqGene reference sequences
    #     Return an hgvs object containing the genomic sequence variant relative to the RefSeqGene
    #     acession
    #     refseq_ac = RefSeqGene ac
    #     """
    #     vr = vvhgvs.validator.Validator(self.hdp)
    #     # parse the variant into hgvs object
    #     var_c = self.hp.parse_hgvs_variant(variant)
    #     # map to the genomic co-ordinates using the easy variant mapper set to alt_aln_method = alt_aln_method
    #     var_g = self.myevm_t_to_g(var_c, evm, self.hdp, primary_assembly)
    #     # Get overlapping transcripts - forcing a splign alignment
    #     start_i = var_g.posedit.pos.start.base
    #     end_i = var_g.posedit.pos.end.base
    #     alt_ac = var_g.ac
    #     alt_aln_method = 'splign'
    #     transcripts = self.hdp.get_tx_for_region(alt_ac, alt_aln_method, start_i - 1, end_i)
    #     # Take the first transcript
    #     ref_g_dict = {
    #         'ref_g': '',
    #         'error': 'false'
    #     }
    #     for trans in transcripts:
    #         tx_ac = trans[0]
    #         try:
    #             ref_c = self.vm.g_to_t(var_g, tx_ac, alt_aln_method='splign')
    #         except:
    #             continue
    #         else:
    #             try:
    #                 ref_g_dict['ref_g'] = self.vm.t_to_g(ref_c, alt_ac=refseq_ac, alt_aln_method='splign')
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
    #
    # def g_to_c(self, var_g, tx_ac, evm):
    #     """
    #     Parses genomic variant strings into hgvs objects
    #     Maps genomic hgvs object into a coding hgvs object if the c accession string is provided
    #     returns a c. variant description string
    #     """
    #     # If the :g. pattern is present in the input variant
    #     if ':g.' in var_g:
    #         # convert the input string into a hgvs object by parsing
    #         var_g = self.hp.parse_hgvs_variant(var_g)
    #         # Map to coding variant
    #         var_c = str(evm.g_to_c(var_g, tx_ac))
    #         return var_c
    #
    # def g_to_n(self, var_g, tx_ac, evm):
    #     """
    #     Parses genomic variant strings into hgvs objects
    #     Maps genomic hgvs object into a non-coding hgvs object if the n accession string is provided
    #     returns a n. variant description string
    #     """
    #     # If the :g. pattern is present in the input variant
    #     if ':g.' in var_g:
    #         # convert the input string into a hgvs object by parsing
    #         var_g = self.hp.parse_hgvs_variant(var_g)
    #         # Map to coding variant
    #         var_n = str(evm.g_to_n(var_g, tx_ac))
    #         return var_n

    def coding(self, variant):
        """
        Ensures variant strings are transcript c. or n.
        returns parsed hgvs c. or n. object
        """
        # If the :c. pattern is present in the input variant
        if ':c.' in variant or ':n.' in variant:
            # convert the input string into a hgvs object
            var_c = self.hp.parse_hgvs_variant(variant)
            return var_c

    def genomic(self, variant, evm, primary_assembly, hn):
        """
        Mapping transcript to genomic position
        Ensures variant strings are transcript c. or n.
        returns parsed hgvs g. object
        """
        # If the :c. pattern is present in the input variant
        if ':c.' in variant or ':n.' in variant:
            hgvs_var = self.hp.parse_hgvs_variant(variant)
            try:
                var_g = self.myevm_t_to_g(hgvs_var, evm, primary_assembly, hn)  # genomic level variant
            except vvhgvs.exceptions.HGVSError as e:
                return 'error ' + str(e)
            return var_g

        # If the :g. pattern is present in the input variant
        elif ':g.' in variant:  # or (pat_n.search(variant)):
            # convert the input string into a hgvs object
            var_g = self.hp.parse_hgvs_variant(variant)
            return var_g

    # def hgvs_genomic(self, variant):
    #     """
    #     Ensures variant strings are g.
    #     returns parsed hgvs g. object
    #     """
    #     # If the :g. pattern is present in the input variant
    #     if ':g.' in variant:
    #         # convert the input string into a hgvs object
    #         var_g = self.hp.parse_hgvs_variant(variant)
    #         return var_g

    def myevm_t_to_g(self, hgvs_c, no_norm_evm, primary_assembly, hn):
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
        # store the input
        stored_hgvs_c = copy.deepcopy(hgvs_c)
        expand_out = False

        # Gap gene black list
        try:
            gene_symbol = self.db.get_gene_symbol_from_transcript_id(hgvs_c.ac)
        except Exception:
            utilise_gap_code = False
        else:
            # If the gene symbol is not in the list, the value False will be returned
            utilise_gap_code = seq_data.gap_black_list(gene_symbol)
        # Warn gap code in use
        logger.debug("gap_compensation_myevm = " + str(utilise_gap_code))

        if utilise_gap_code is True and (hgvs_c.posedit.edit.type == 'identity' or hgvs_c.posedit.edit.type == 'del'
                                         or hgvs_c.posedit.edit.type == 'delins' or hgvs_c.posedit.edit.type == 'dup'
                                         or hgvs_c.posedit.edit.type == 'sub' or hgvs_c.posedit.edit.type == 'ins'
                                         or hgvs_c.posedit.edit.type == 'inv'):

            # if NM_ need the n. position
            if str(hgvs_c.ac).startswith('NM_'):
                hgvs_c = no_norm_evm.c_to_n(hgvs_c)

            # Check for intronic
            try:
                hn.normalize(hgvs_c)
            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)
                if 'intronic variant' not in error and \
                        'Length implied by coordinates must equal sequence deletion length' in error and \
                        hgvs_c.ac.startswith('NR_'):
                    hgvs_c.posedit.pos.end.base = hgvs_c.posedit.pos.start.base + len(hgvs_c.posedit.edit.ref) - 1

            # Check again before continuing
            if re.search(r'\d+\+', str(hgvs_c.posedit.pos)) or re.search(r'\d+-', str(hgvs_c.posedit.pos)) or \
                    re.search(r'\*\d+\+', str(hgvs_c.posedit.pos)) or re.search(r'\*\d+-', str(hgvs_c.posedit.pos)):
                pass

            else:
                try:
                    # For non-intronic sequence
                    hgvs_t = copy.deepcopy(hgvs_c)
                    if hgvs_t.posedit.edit.type == 'inv':
                        inv_alt = self.revcomp(hgvs_t.posedit.edit.ref)
                        t_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(
                            hgvs_t.posedit.pos.end.base) + 'del' + hgvs_t.posedit.edit.ref + 'ins' + inv_alt
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
                    expand_out = True

                except Exception:
                    hgvs_c = hgvs_c

            if str(hgvs_c.ac).startswith('NM_'):
                try:
                    hgvs_c = no_norm_evm.n_to_c(hgvs_c)
                except vvhgvs.exceptions.HGVSError:
                    hgvs_c = copy.deepcopy(stored_hgvs_c)

            # Ensure the altered c. variant has not crossed intro exon boundaries
            hgvs_check_boundaries = copy.deepcopy(hgvs_c)
            try:
                hn.normalize(hgvs_check_boundaries)
            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)
                if 'spanning the exon-intron boundary' in error:
                    hgvs_c = copy.deepcopy(stored_hgvs_c)
            # Catch identity at the exon/intron boundary by trying to normalize ref only
            if hgvs_check_boundaries.posedit.edit.type == 'identity':
                reform_ident = str(hgvs_c).split(':')[0]
                reform_ident = reform_ident + ':' + stored_hgvs_c.type + '.' + str(hgvs_c.posedit.pos) + 'del' + str(
                    hgvs_c.posedit.edit.ref)  # + 'ins' + str(hgvs_c.posedit.edit.alt)
                hgvs_reform_ident = self.hp.parse_hgvs_variant(reform_ident)
                try:
                    hn.normalize(hgvs_reform_ident)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'spanning the exon-intron boundary' in error or 'Normalization of intronic variants' in error:
                        hgvs_c = copy.deepcopy(stored_hgvs_c)

        # Capture errors from attempted mappings
        attempted_mapping_error = ''
        hgvs_genomic = None

        try:
            hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)
            hn.normalize(hgvs_genomic)  # Check the validity of the mapping
            # This will fail on multiple refs for NC_
        except vvhgvs.exceptions.HGVSError:
            # Recover all available mapping options from UTA
            mapping_options = self.hdp.get_tx_mapping_options(hgvs_c.ac)

            if not mapping_options:
                raise HGVSDataNotAvailableError(
                    "No alignment data between the specified transcript reference sequence and any GRCh37 and GRCh38 "
                    "genomic reference sequences (including alternate chromosome assemblies, patches and RefSeqGenes) "
                    "are available.")

            def search_through_options(hgvs_genomic, seqtype, chr_num_val, final=False):
                err = ''
                for option in mapping_options:
                    if option[2].startswith('blat'):
                        continue
                    if option[1].startswith(seqtype):
                        chr_num = seq_data.supported_for_mapping(str(option[1]), primary_assembly)
                        if final:
                            try:
                                hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
                                break
                            except Exception as e:
                                err += str(e) + "/" + hgvs_c.ac + "/" + option[1] + '~'
                                continue
                        if chr_num_val and chr_num != 'false':
                            try:
                                hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
                                break
                            except Exception as e:
                                err += str(e) + "/" + hgvs_c.ac + "/" + option[1] + '~'
                                continue
                        elif chr_num_val is False and chr_num == 'false':
                            try:
                                hgvs_genomic = self.vm.t_to_g(hgvs_c, str(option[1]))
                                break
                            except Exception as e:
                                err += str(e) + "/" + hgvs_c.ac + "/" + option[1] + '~'
                                continue

                return hgvs_genomic, err

            hgvs_genomic, new_error = search_through_options(hgvs_genomic, 'NC_', True)
            attempted_mapping_error += new_error

            # If not mapped, raise error
            try:
                hn.normalize(hgvs_genomic)
            except:
                hgvs_genomic, new_error = search_through_options(hgvs_genomic, 'NC_', False)
                attempted_mapping_error += new_error

                try:
                    hn.normalize(hgvs_genomic)
                except:
                    hgvs_genomic, new_error = search_through_options(hgvs_genomic, 'NT_', True)
                    attempted_mapping_error += new_error

                    try:
                        hn.normalize(hgvs_genomic)
                    except:
                        hgvs_genomic, new_error = search_through_options(hgvs_genomic, 'NT_', False)
                        attempted_mapping_error += new_error

                        try:
                            hn.normalize(hgvs_genomic)
                        except:
                            hgvs_genomic, new_error = search_through_options(hgvs_genomic, 'NW_', True)
                            attempted_mapping_error += new_error

                            try:
                                hn.normalize(hgvs_genomic)
                            except:
                                hgvs_genomic, new_error = search_through_options(hgvs_genomic, 'NW_', False)
                                attempted_mapping_error += new_error

                                # Only a RefSeqGene available
                                try:
                                    hn.normalize(hgvs_genomic)
                                except:
                                    hgvs_genomic, new_error = search_through_options(hgvs_genomic, 'NG_', True,
                                                                                     final=True)
                                    attempted_mapping_error += new_error

        # If not mapped, raise error
        if hgvs_genomic is None:
            raise HGVSDataNotAvailableError(attempted_mapping_error)

        if hgvs_c.posedit.edit.type == 'identity' and hgvs_genomic.posedit.edit.type == 'delins' and \
                hgvs_genomic.posedit.edit.alt == '' and not expand_out:
            hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.ref
        if hgvs_genomic.posedit.edit.type == 'ins' and utilise_gap_code is True:
            try:
                hgvs_genomic = hn.normalize(hgvs_genomic)
            except vvhgvs.exceptions.HGVSError as e:
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
        if (stored_hgvs_c.posedit.edit.ref == '' or stored_hgvs_c.posedit.edit.ref is None) and expand_out:
            if stored_hgvs_c.type == 'c':
                stored_hgvs_n = self.vm.c_to_n(stored_hgvs_c)
            else:
                stored_hgvs_n = stored_hgvs_c
            stored_ref = self.sf.fetch_seq(str(stored_hgvs_n.ac), stored_hgvs_n.posedit.pos.start.base - 1,
                                           stored_hgvs_n.posedit.pos.end.base)
            stored_hgvs_c.posedit.edit.ref = stored_ref

        if (hgvs_genomic.posedit.edit.ref == '' or hgvs_genomic.posedit.edit.ref is None) and expand_out:
            if hgvs_genomic.posedit.edit.type == 'ins':
                stored_ref = self.sf.fetch_seq(str(hgvs_genomic.ac), hgvs_genomic.posedit.pos.start.base - 1,
                                               hgvs_genomic.posedit.pos.end.base)
                stored_alt = stored_ref[:1] + hgvs_genomic.posedit.edit.alt + stored_ref[-1:]
                hgvs_genomic.posedit.edit.ref = stored_ref
                hgvs_genomic.posedit.edit.alt = stored_alt

        # First look for variants mapping to the flanks of gaps
        # either in the gap or on the flank but not fully within the gap
        if expand_out:
            nr_genomic = self.nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)

            try:
                hn.normalize(nr_genomic)
            except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                error_type_1 = str(e)
                if 'Length implied by coordinates must equal sequence deletion length' in str(e) or str(
                        e) == 'base start position must be <= end position':
                    # Effectively, this code is designed to handle variants that are directly proximal to
                    # gap BOUNDARIES, but in some cases the replace reference function of hgvs mapping has removed
                    # bases due to the deletion length being > the specified range.
                    genomic_gap_variant = None
                    # Warn of variant location wrt the gap
                    if 'Length implied by coordinates must equal sequence deletion length' in str(e):
                        logger.info('Variant is proximal to the flank of a genomic gap')
                        genomic_gap_variant = self.vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)
                        try:
                            hn.normalize(genomic_gap_variant)
                        # Still a problem
                        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                            if 'base start position must be <= end position' in str(e) and \
                                    'Length implied by coordinates must equal' in error_type_1:
                                make_gen_var = copy.copy(nr_genomic)
                                make_gen_var.posedit.edit.ref = self.sf.fetch_seq(
                                    nr_genomic.ac,
                                    nr_genomic.posedit.pos.start.base - 1,
                                    nr_genomic.posedit.pos.end.base
                                )
                                genomic_gap_variant = make_gen_var

                                error_type_1 = None
                        else:
                            genomic_gap_variant = self.nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)

                    if error_type_1 == 'base start position must be <= end position':
                        logger.info('Variant is fully within a genomic gap')
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
                        if 'Length implied by coordinates must equal sequence deletion length' in str(e):
                            # This will only happen if the variant is flanking the gap but is
                            # not inside the gap
                            logger.info('Variant is on the flank of a genomic gap but not within the gap')

                            # Test on the flank and if so, return

                            # Logic, normalize the c. variant and if a substitution (cannot normalize) then direct map
                            # Currently believe that sub.n is the only variant type which fits. ins can normalize
                            # and may also be a dup! Added identity also
                            try:
                                try:
                                    norm_stored_c = hn.normalize(stored_hgvs_c)
                                except HGVSUnsupportedOperationError:
                                    norm_stored_c = stored_hgvs_c
                                if norm_stored_c.posedit.edit.type == 'sub' or \
                                        norm_stored_c.posedit.edit.type == 'identity':
                                    flank_hgvs_genomic = self.vm.t_to_g(norm_stored_c, genomic_gap_variant.ac)
                                    self.vr.validate(flank_hgvs_genomic)
                                    return flank_hgvs_genomic

                            # Will occur if the variant still overlaps the gap / is in the gap
                            except HGVSInvalidVariantError:
                                pass

                            # If test fails, continue old processing
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
                            logger.debug("Except passed, %s", e)

                        # Should be a delins so will normalize statically and replace the reference bases
                        genomic_gap_variant = hn.normalize(genomic_gap_variant)
                        # Static map to c. and static normalize
                        transcript_gap_variant = self.vm.g_to_t(genomic_gap_variant, hgvs_c.ac)

                        if 'Length implied by coordinates must equal sequence deletion length' not in str(e):
                            try:
                                transcript_gap_variant = hn.normalize(transcript_gap_variant)
                            except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                                logger.debug("Except passed, %s", e)

                        # if NM_ need the n. position
                        if str(hgvs_c.ac).startswith('NM_'):
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
                                transcript_gap_n_delins_from_dup = fn.hgvs_dup2indel(transcript_gap_n)
                                transcript_gap_n = self.hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)
                                transcript_gap_alt_n_delins_from_dup = fn.hgvs_dup2indel(transcript_gap_alt_n)
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
                        for i in range(transcript_gap_alt_n.posedit.pos.start.base,
                                       transcript_gap_alt_n.posedit.pos.end.base + 1, 1):
                            if i == alt_start:
                                alt_base_dict[i] = str(''.join(alternate_bases))
                            else:
                                alt_base_dict[i] = 'X'

                                # Generate the alt sequence
                        alternate_sequence_bases = []
                        for i in range(transcript_gap_n.posedit.pos.start.base,
                                       transcript_gap_n.posedit.pos.end.base + 1,
                                       1):
                            if i in list(alt_base_dict.keys()):
                                alternate_sequence_bases.append(alt_base_dict[i])
                            elif i in list(ref_base_dict.keys()):
                                alternate_sequence_bases.append(ref_base_dict[i])
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
                                pre_base = self.sf.fetch_seq(
                                    transcript_gap_n.ac,
                                    transcript_gap_n.posedit.pos.start.base - 2,
                                    transcript_gap_n.posedit.pos.start.base - 1)
                                post_base = self.sf.fetch_seq(
                                    transcript_gap_n.ac,
                                    transcript_gap_n.posedit.pos.end.base,
                                    transcript_gap_n.posedit.pos.end.base + 1)
                                transcript_gap_n.posedit.pos.start.base = transcript_gap_n.posedit.pos.start.base - 1
                                transcript_gap_n.posedit.pos.end.base = transcript_gap_n.posedit.pos.end.base + 1
                                transcript_gap_n.posedit.edit.ref = pre_base + transcript_gap_n.posedit.edit.ref + \
                                    post_base
                                transcript_gap_n.posedit.edit.alt = pre_base + transcript_gap_n.posedit.edit.alt + \
                                    post_base
                                try:
                                    transcript_gap_variant = self.vm.n_to_c(transcript_gap_n)
                                except:
                                    transcript_gap_variant = transcript_gap_n
                                hgvs_genomic = self.vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                                hgvs_genomic = hn.normalize(hgvs_genomic)

                        # Bypass the next bit of gap code
                        expand_out = False

        # CASCADING STATEMENTS WHICH CAPTURE t to g MAPPING OPTIONS
        # Remove identity bases
        if hgvs_c == stored_hgvs_c:
            pass
        elif expand_out is False or utilise_gap_code is False:
            pass
        # Correct expansion ref + 2
        elif expand_out and (
                len(hgvs_genomic.posedit.edit.ref) == (len(stored_hgvs_c.posedit.edit.ref) + 2)):  # >= 3:
            hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base + 1
            hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base - 1
            hgvs_genomic.posedit.edit.ref = hgvs_genomic.posedit.edit.ref[1:-1]
            if hgvs_genomic.posedit.edit.alt is not None:
                hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt[1:-1]
        elif expand_out and (
                len(hgvs_genomic.posedit.edit.ref) != (len(stored_hgvs_c.posedit.edit.ref) + 2)):  # >= 3:
            if len(hgvs_genomic.posedit.edit.ref) == 2:
                hn.normalize(hgvs_genomic)

            # Likely if the start or end position aligns to a gap in the genomic sequence
            # Logic
            # We have checked that the variant does not cross boundaries, or is intronic
            # So is likely mapping to a genomic gap
            elif len(hgvs_genomic.posedit.edit.ref) <= 1:
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
                        logger.debug("Except passed, %s", e)
                    # Should be a delins so will normalize statically and replace the reference bases
                    genomic_gap_variant = hn.normalize(genomic_gap_variant)
                    # Static map to c. and static normalize
                    transcript_gap_variant = self.vm.g_to_t(genomic_gap_variant, hgvs_c.ac)
                    transcript_gap_variant = hn.normalize(transcript_gap_variant)
                    # if NM_ need the n. position
                    if str(hgvs_c.ac).startswith('NM_'):
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
                            transcript_gap_n_delins_from_dup = fn.hgvs_dup2indel(transcript_gap_n)
                            transcript_gap_n = self.hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)
                            transcript_gap_alt_n_delins_from_dup = fn.hgvs_dup2indel(transcript_gap_alt_n)
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
                    for i in range(transcript_gap_alt_n.posedit.pos.start.base,
                                   transcript_gap_alt_n.posedit.pos.end.base + 1, 1):
                        if i == alt_start:
                            alt_base_dict[i] = str(''.join(alternate_bases))
                        else:
                            alt_base_dict[i] = 'X'

                    # Generate the alt sequence
                    alternate_sequence_bases = []
                    for i in range(transcript_gap_n.posedit.pos.start.base,
                                   transcript_gap_n.posedit.pos.end.base + 1, 1):
                        if i in list(alt_base_dict.keys()):
                            alternate_sequence_bases.append(alt_base_dict[i])
                        elif i in list(ref_base_dict.keys()):
                            alternate_sequence_bases.append(ref_base_dict[i])
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
                            pre_base = self.sf.fetch_seq(transcript_gap_n.ac,
                                                         transcript_gap_n.posedit.pos.start.base - 2,
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
        if hgvs_c.posedit.edit.type == 'ins' and hgvs_c.posedit.pos.start.offset == 0 and \
                hgvs_c.posedit.pos.end.offset == 0:
            try:
                hn.normalize(hgvs_genomic)
            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)
                if error == 'insertion length must be 1':
                    if hgvs_c.type == 'c':
                        hgvs_t = self.vm.c_to_n(hgvs_c)
                    else:
                        hgvs_t = copy.copy(hgvs_c)
                    ins_ref = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 1,
                                                hgvs_t.posedit.pos.end.base)
                    ins_alt = ins_ref[:1] + hgvs_t.posedit.edit.alt + ins_ref[-1:]
                    ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(
                        hgvs_t.posedit.pos.start.base) + '_' + str(
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

    def noreplace_myevm_t_to_g(self, hgvs_c, variant):
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
        hgvs_genomic = None
        attempted_mapping_error = ''
        try:
            hgvs_genomic = variant.evm.t_to_g(hgvs_c)
            variant.hn.normalize(hgvs_genomic)
        # This will fail on multiple refs for NC_
        except vvhgvs.exceptions.HGVSError:
            # Recover all available mapping options from UTA
            mapping_options = self.hdp.get_tx_mapping_options(hgvs_c.ac)

            if not mapping_options:
                raise HGVSDataNotAvailableError("no g. mapping options available")

            def search_in_options(hgvs_genomic, seqtype, chr_num_val, final=False):
                err = ''
                for op in mapping_options:
                    if op[2].startswith('blat'):
                        continue
                    if op[1].startswith(seqtype):
                        if final:
                            try:
                                hgvs_genomic = self.vm.t_to_g(hgvs_c, str(op[1]))
                                break
                            except Exception as e:
                                err += str(e) + "/" + hgvs_c.ac + "/" + op[1] + '~'
                                continue
                        chr_num = seq_data.supported_for_mapping(str(op[1]), variant.primary_assembly)
                        if chr_num_val and chr_num != 'false':
                            try:
                                hgvs_genomic = self.vm.t_to_g(hgvs_c, str(op[1]))
                                break
                            except Exception as e:
                                err += str(e) + "/" + hgvs_c.ac + "/" + op[1] + '~'
                                continue
                        elif not chr_num_val and chr_num == 'false':
                            try:
                                hgvs_genomic = self.vm.t_to_g(hgvs_c, str(op[1]))
                                break
                            except Exception as e:
                                err += str(e) + "/" + hgvs_c.ac + "/" + op[1] + '~'
                                continue
                return hgvs_genomic, err

            hgvs_genomic, new_errors = search_in_options(hgvs_genomic, 'NC_', True)
            attempted_mapping_error += new_errors

            # If not mapped, raise error
            try:
                variant.hn.normalize(hgvs_genomic)
            except:
                hgvs_genomic, new_errors = search_in_options(hgvs_genomic, 'NC_', True)
                attempted_mapping_error += new_errors

                # If not mapped, raise error
                try:
                    variant.hn.normalize(hgvs_genomic)
                except:
                    hgvs_genomic, new_errors = search_in_options(hgvs_genomic, 'NC_', False)
                    attempted_mapping_error += new_errors
                    try:
                        variant.hn.normalize(hgvs_genomic)
                    except:
                        hgvs_genomic, new_errors = search_in_options(hgvs_genomic, 'NT_', True)
                        attempted_mapping_error += new_errors
                        try:
                            variant.hn.normalize(hgvs_genomic)
                        except:
                            hgvs_genomic, new_errors = search_in_options(hgvs_genomic, 'NT_', False)
                            attempted_mapping_error += new_errors
                            try:
                                variant.hn.normalize(hgvs_genomic)
                            except:
                                hgvs_genomic, new_errors = search_in_options(hgvs_genomic, 'NW_', True)
                                attempted_mapping_error += new_errors
                                try:
                                    variant.hn.normalize(hgvs_genomic)
                                except:
                                    hgvs_genomic, new_errors = search_in_options(hgvs_genomic, 'NW_', False)
                                    attempted_mapping_error += new_errors
                                    # Only a RefSeqGene available
                                    try:
                                        variant.hn.normalize(hgvs_genomic)
                                    except:
                                        hgvs_genomic, new_errors = search_in_options(hgvs_genomic, 'NG_', True,
                                                                                     final=True)
                                        attempted_mapping_error += new_errors
        if hgvs_genomic is None:
            raise HGVSDataNotAvailableError('No available t_to_g liftover')

        # Ins variants map badly - Especially between c. exon/exon boundary
        if hgvs_c.posedit.edit.type == 'ins' and hgvs_c.posedit.pos.start.offset == 0 and \
                hgvs_c.posedit.pos.end.offset == 0:
            try:
                variant.hn.normalize(hgvs_genomic)
            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)
                if error == 'insertion length must be 1':
                    if hgvs_c.type == 'c':
                        hgvs_t = self.vm.c_to_n(hgvs_c)
                    else:
                        hgvs_t = copy.copy(hgvs_c)
                    ins_ref = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 1,
                                                hgvs_t.posedit.pos.end.base)
                    ins_alt = ins_ref[:1] + hgvs_t.posedit.edit.alt + ins_ref[-1:]
                    ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + \
                        str(hgvs_t.posedit.pos.end.base) + 'del' + ins_ref + 'ins' + ins_alt
                    hgvs_t = self.hp.parse_hgvs_variant(ins_to_delins)
                    try:
                        hgvs_c = self.vm.n_to_c(hgvs_t)
                    except Exception:
                        hgvs_c = copy.copy(hgvs_t)
                    try:
                        hgvs_genomic = variant.no_norm_evm.t_to_g(hgvs_c)
                    except Exception as e:
                        error = str(e)
                        logger.warning('Ins mapping error in myt_to_g ' + error)

        return hgvs_genomic

    def myevm_g_to_t(self, evm, hgvs_genomic, alt_ac):
        """
        Enhanced transcript to genome position on a specified genomic reference using vm
        Deals with mapping from transcript positions that do not exist in the genomic sequence
        i.e. the stated position aligns to a genomic gap!
        returns parsed hgvs g. object
        """
        hgvs_t = evm.g_to_t(hgvs_genomic, alt_ac)
        return hgvs_t

    def myvm_t_to_g(self, hgvs_c, alt_chr, no_norm_evm, hn):
        # store the input
        stored_hgvs_c = copy.deepcopy(hgvs_c)
        expand_out = False

        # Gap gene black list
        try:
            gene_symbol = self.db.get_gene_symbol_from_transcript_id(hgvs_c.ac)
        except Exception:
            utilise_gap_code = False
        else:
            # If the gene symbol is not in the list, the value False will be returned
            utilise_gap_code = seq_data.gap_black_list(gene_symbol)
        # Warn gap code in use
        logger.debug("gap_compensation_mvm = " + str(utilise_gap_code))

        if utilise_gap_code and (hgvs_c.posedit.edit.type == 'identity' or hgvs_c.posedit.edit.type == 'del'
                                 or hgvs_c.posedit.edit.type == 'delins' or hgvs_c.posedit.edit.type == 'dup'
                                 or hgvs_c.posedit.edit.type == 'sub' or hgvs_c.posedit.edit.type == 'ins'
                                 or hgvs_c.posedit.edit.type == 'inv'):

            # if NM_ need the n. position
            if str(hgvs_c.ac).startswith("NM_"):
                hgvs_c = no_norm_evm.c_to_n(hgvs_c)

            # Check for intronic
            try:
                hn.normalize(hgvs_c)
            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)
                if 'intronic variant' in error:
                    logger.debug("Except passed, %s", e)
                elif 'Length implied by coordinates must equal sequence deletion length' in error and \
                        hgvs_c.ac.startswith('NR_'):
                    hgvs_c.posedit.pos.end.base = hgvs_c.posedit.pos.start.base + len(hgvs_c.posedit.edit.ref) - 1

            # Check again before continuing
            if re.search(r'\d+\+', str(hgvs_c.posedit.pos)) or re.search(r'\d+-', str(hgvs_c.posedit.pos)) or re.search(
                    r'\*\d+\+', str(hgvs_c.posedit.pos)) or re.search(r'\*\d+-', str(hgvs_c.posedit.pos)):
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
                    expand_out = True

                except Exception:
                    hgvs_c = hgvs_c

            if str(hgvs_c.ac).startswith('NM_'):
                try:
                    hgvs_c = no_norm_evm.n_to_c(hgvs_c)
                except vvhgvs.exceptions.HGVSError:
                    hgvs_c = copy.deepcopy(stored_hgvs_c)

            # Ensure the altered c. variant has not crossed intro exon boundaries
            hgvs_check_boundaries = copy.deepcopy(hgvs_c)
            try:
                hn.normalize(hgvs_check_boundaries)
            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)
                if 'spanning the exon-intron boundary' in error:
                    hgvs_c = copy.deepcopy(stored_hgvs_c)
            # Catch identity at the exon/intron boundary by trying to normalize ref only
            if hgvs_check_boundaries.posedit.edit.type == 'identity':
                reform_ident = str(hgvs_c).split(':')[0]
                reform_ident = reform_ident + ':' + stored_hgvs_c.type + '.' + str(hgvs_c.posedit.pos) + 'del' + str(
                    hgvs_c.posedit.edit.ref)  # + 'ins' + str(hgvs_c.posedit.edit.alt)
                hgvs_reform_ident = self.hp.parse_hgvs_variant(reform_ident)
                try:
                    hn.normalize(hgvs_reform_ident)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'spanning the exon-intron boundary' in error or 'Normalization of intronic variants' in error:
                        hgvs_c = copy.deepcopy(stored_hgvs_c)

        hgvs_genomic = self.vm.t_to_g(hgvs_c, alt_chr)
        if hgvs_c.posedit.edit.type == 'identity' and hgvs_genomic.posedit.edit.type == 'delins' and \
                hgvs_genomic.posedit.edit.alt == '' and expand_out is False:
            hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.ref
        if hgvs_genomic.posedit.edit.type == 'ins' and utilise_gap_code is True:
            try:
                hgvs_genomic = hn.normalize(hgvs_genomic)
            except vvhgvs.exceptions.HGVSError as e:
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
        if (stored_hgvs_c.posedit.edit.ref == '' or stored_hgvs_c.posedit.edit.ref is None) and expand_out:
            if stored_hgvs_c.type == 'c':
                stored_hgvs_n = self.vm.c_to_n(stored_hgvs_c)
            else:
                stored_hgvs_n = stored_hgvs_c
            stored_ref = self.sf.fetch_seq(str(stored_hgvs_n.ac), stored_hgvs_n.posedit.pos.start.base - 1,
                                           stored_hgvs_n.posedit.pos.end.base)
            stored_hgvs_c.posedit.edit.ref = stored_ref

        if (hgvs_genomic.posedit.edit.ref == '' or hgvs_genomic.posedit.edit.ref is None) and expand_out:
            if hgvs_genomic.posedit.edit.type == 'ins':
                stored_ref = self.sf.fetch_seq(str(hgvs_genomic.ac), hgvs_genomic.posedit.pos.start.base - 1,
                                               hgvs_genomic.posedit.pos.end.base)
                stored_alt = stored_ref[:1] + hgvs_genomic.posedit.edit.alt + stored_ref[-1:]
                hgvs_genomic.posedit.edit.ref = stored_ref
                hgvs_genomic.posedit.edit.alt = stored_alt

        # First look for variants mapping to the flanks of gaps
        # either in the gap or on the flank but not fully within the gap
        if expand_out:
            nr_genomic = self.nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)
            try:
                hn.normalize(nr_genomic)
            except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                error_type_1 = str(e)
                if 'Length implied by coordinates must equal sequence deletion length' in str(e) or str(
                        e) == 'base start position must be <= end position':
                    # Effectively, this code is designed to handle variants that are directly proximal to
                    # gap BOUNDARIES, but in some cases the replace reference function of hgvs mapping has removed bases
                    # due to the deletion length being > the specified range.
                    genomic_gap_variant = None
                    # Warn of variant location wrt the gap
                    if 'Length implied by coordinates must equal sequence deletion length' in str(e):
                        logger.info('Variant is proximal to the flank of a genomic gap')
                        genomic_gap_variant = self.vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)
                        try:
                            hn.normalize(genomic_gap_variant)
                        # Still a problem
                        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                            if 'base start position must be <= end position' in str(e) and \
                                    'Length implied by coordinates must equal' in error_type_1:
                                make_gen_var = copy.copy(nr_genomic)
                                make_gen_var.posedit.edit.ref = self.sf.fetch_seq(
                                    nr_genomic.ac,
                                    nr_genomic.posedit.pos.start.base - 1,
                                    nr_genomic.posedit.pos.end.base
                                )
                                genomic_gap_variant = make_gen_var
                                error_type_1 = None
                        else:
                            genomic_gap_variant = self.nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)

                    if error_type_1 == 'base start position must be <= end position':
                        logger.info('Variant is fully within a genomic gap')
                        genomic_gap_variant = self.vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)

                    # Logic
                    # We have checked that the variant does not cross boundaries, or is intronic
                    # So is likely mapping to a genomic gap
                    try:
                        hn.normalize(genomic_gap_variant)
                    except Exception as ea1:
                        if str(ea1) == 'base start position must be <= end position':
                            # This will only happen when the variant is fully within the gap
                            gap_start = genomic_gap_variant.posedit.pos.end.base
                            gap_end = genomic_gap_variant.posedit.pos.start.base
                            genomic_gap_variant.posedit.pos.start.base = gap_start
                            genomic_gap_variant.posedit.pos.end.base = gap_end
                        if 'Length implied by coordinates must equal sequence deletion length' in str(ea1):
                            # This will only happen if the variant is flanking the gap but is
                            # not inside the gap
                            logger.info('Variant is on the flank of a genomic gap but not within the gap')

                            # Test definately on the flank and if so, return
                            # Logic, normalize the c. variant and if a substitution (cannot normalize) then direct map
                            # Currently believe that sub.n is the only variant type which fits. ins can normalize
                            # and may also be a dup!
                            try:
                                norm_stored_c = hn.normalize(stored_hgvs_c)
                                if norm_stored_c.posedit.edit.type == 'sub' or \
                                        norm_stored_c.posedit.edit.type == 'identity':
                                    flank_hgvs_genomic = self.vm.t_to_g(norm_stored_c, genomic_gap_variant.ac)
                                    self.vr.validate(flank_hgvs_genomic)
                                    return flank_hgvs_genomic

                            # Will occur if the variant still overlaps the gap / is in the gap
                            except HGVSInvalidVariantError:
                                pass

                            # If test fails, continue old processing
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
                            logger.debug("Except passed, %s", e)

                        # Should be a delins so will normalize statically and replace the reference bases
                        genomic_gap_variant = hn.normalize(genomic_gap_variant)
                        # Static map to c. and static normalize
                        transcript_gap_variant = self.vm.g_to_t(genomic_gap_variant, hgvs_c.ac)
                        if 'Length implied by coordinates must equal sequence deletion length' not in str(ea1):
                            try:
                                transcript_gap_variant = hn.normalize(transcript_gap_variant)
                            except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                                logger.debug("Except passed, %s", e)

                        # if NM_ need the n. position
                        if str(hgvs_c.ac).startswith('NM_'):
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
                                transcript_gap_n_delins_from_dup = fn.hgvs_dup2indel(transcript_gap_n)
                                transcript_gap_n = self.hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)
                                transcript_gap_alt_n_delins_from_dup = fn.hgvs_dup2indel(transcript_gap_alt_n)
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
                        for i in range(transcript_gap_alt_n.posedit.pos.start.base,
                                       transcript_gap_alt_n.posedit.pos.end.base + 1, 1):
                            if i == alt_start:
                                alt_base_dict[i] = str(''.join(alternate_bases))
                            else:
                                alt_base_dict[i] = 'X'

                        # Generate the alt sequence
                        alternate_sequence_bases = []
                        for i in range(transcript_gap_n.posedit.pos.start.base,
                                       transcript_gap_n.posedit.pos.end.base + 1,
                                       1):
                            if i in list(alt_base_dict.keys()):
                                alternate_sequence_bases.append(alt_base_dict[i])
                            elif i in list(ref_base_dict.keys()):
                                alternate_sequence_bases.append(ref_base_dict[i])
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
                                pre_base = self.sf.fetch_seq(transcript_gap_n.ac,
                                                             transcript_gap_n.posedit.pos.start.base - 2,
                                                             transcript_gap_n.posedit.pos.start.base - 1)
                                post_base = self.sf.fetch_seq(transcript_gap_n.ac,
                                                              transcript_gap_n.posedit.pos.end.base,
                                                              transcript_gap_n.posedit.pos.end.base + 1)
                                transcript_gap_n.posedit.pos.start.base = transcript_gap_n.posedit.pos.start.base - 1
                                transcript_gap_n.posedit.pos.end.base = transcript_gap_n.posedit.pos.end.base + 1
                                transcript_gap_n.posedit.edit.ref = pre_base + transcript_gap_n.posedit.edit.ref + \
                                                                    post_base
                                transcript_gap_n.posedit.edit.alt = pre_base + transcript_gap_n.posedit.edit.alt + \
                                                                    post_base
                                try:
                                    transcript_gap_variant = self.vm.n_to_c(transcript_gap_n)
                                except:
                                    transcript_gap_variant = transcript_gap_n
                                hgvs_genomic = self.vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                                hgvs_genomic = hn.normalize(hgvs_genomic)

                        # Bypass the next bit of gap code
                        expand_out = False

        # CASCADING STATEMENTS WHICH CAPTURE t to g MAPPING OPTIONS
        # Remove identity bases
        if hgvs_c == stored_hgvs_c:
            expand_out = False
        elif expand_out is False or utilise_gap_code is False:
            pass
        # Correct expansion ref + 2
        elif expand_out and (
                len(hgvs_genomic.posedit.edit.ref) == (len(stored_hgvs_c.posedit.edit.ref) + 2)):  # >= 3:
            hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base + 1
            hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base - 1
            hgvs_genomic.posedit.edit.ref = hgvs_genomic.posedit.edit.ref[1:-1]
            if hgvs_genomic.posedit.edit.alt is not None:
                hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt[1:-1]
        elif expand_out and (
                len(hgvs_genomic.posedit.edit.ref) != (len(stored_hgvs_c.posedit.edit.ref) + 2)):  # >= 3:
            if len(hgvs_genomic.posedit.edit.ref) == 2:
                hn.normalize(hgvs_genomic)

            # Likely if the start or end position aligns to a gap in the genomic sequence
            # Logic
            # We have checked that the variant does not cross boundaries, or is intronic
            # So is likely mapping to a genomic gap
            elif len(hgvs_genomic.posedit.edit.ref) <= 1:
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
                        logger.debug("Except passed, %s", e)
                    # Should be a delins so will normalize statically and replace the reference bases
                    genomic_gap_variant = hn.normalize(genomic_gap_variant)
                    # Static map to c. and static normalize
                    transcript_gap_variant = self.vm.g_to_t(genomic_gap_variant, hgvs_c.ac)
                    transcript_gap_variant = hn.normalize(transcript_gap_variant)
                    # if NM_ need the n. position
                    if str(hgvs_c.ac).startswith('NM_'):
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
                            transcript_gap_n_delins_from_dup = fn.hgvs_dup2indel(transcript_gap_n)
                            transcript_gap_n = self.hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)
                            transcript_gap_alt_n_delins_from_dup = fn.hgvs_dup2indel(transcript_gap_alt_n)
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
                    for i in range(transcript_gap_alt_n.posedit.pos.start.base,
                                   transcript_gap_alt_n.posedit.pos.end.base + 1, 1):
                        if i == alt_start:
                            alt_base_dict[i] = str(''.join(alternate_bases))
                        else:
                            alt_base_dict[i] = 'X'

                    # Generate the alt sequence
                    alternate_sequence_bases = []
                    for i in range(transcript_gap_n.posedit.pos.start.base,
                                   transcript_gap_n.posedit.pos.end.base + 1, 1):
                        if i in list(alt_base_dict.keys()):
                            alternate_sequence_bases.append(alt_base_dict[i])
                        elif i in list(ref_base_dict.keys()):
                            alternate_sequence_bases.append(ref_base_dict[i])
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
                            pre_base = self.sf.fetch_seq(transcript_gap_n.ac,
                                                         transcript_gap_n.posedit.pos.start.base - 2,
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
        if hgvs_c.posedit.edit.type == 'ins' and hgvs_c.posedit.pos.start.offset == 0 and \
                hgvs_c.posedit.pos.end.offset == 0:
            try:
                hn.normalize(hgvs_genomic)
            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)
                if error == 'insertion length must be 1':
                    if hgvs_c.type == 'c':
                        hgvs_t = self.vm.c_to_n(hgvs_c)
                    else:
                        hgvs_t = copy.copy(hgvs_c)
                    ins_ref = self.sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 1,
                                                hgvs_t.posedit.pos.end.base)
                    ins_alt = ins_ref[:1] + hgvs_t.posedit.edit.alt + ins_ref[-1:]
                    ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + \
                        str(hgvs_t.posedit.pos.end.base) + 'del' + ins_ref + 'ins' + ins_alt
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

    # def hgvs_protein(self, variant, hpOld):
    #     """
    #     parse p. strings into hgvs p. objects
    #     """
    #     # If the :p. pattern is present in the input variant
    #     if ':p.' in variant:
    #         # convert the input string into a hgvs object
    #         var_p = self.hp.parse_hgvs_variant(variant)
    #         return var_p

    def hgvs_r_to_c(self, hgvs_object):
        """
        Convert r. into c.
        """
        # check for LRG_t with r.
        if 'LRG' in hgvs_object.ac:
            transcript_ac = self.db.get_refseq_transcript_id_from_lrg_transcript_id(hgvs_object.ac)
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

    # def hgvs_c_to_r(self, hgvs_object):
    #     """
    #     Convert c. into r.
    #     """
    #     hgvs_object.type = 'r'
    #     edit = str(hgvs_object.posedit.edit)
    #     edit = edit.lower()
    #     edit = edit.replace('t', 'u')
    #     hgvs_object.posedit.edit = edit
    #     return hgvs_object

    # def tx_identity_info(self, variant):
    #     """
    #     Input c. r. n. variant string
    #     Use uta.py (hdp) to return the identity information for the transcript variant
    #     see vvhgvs.dataproviders.uta.py for details
    #     """
    #     # If the :c. pattern is present in the input variant
    #     if ':c.' in variant:
    #         # Remove all text to the right and including pat_c
    #         tx_ac = variant[:variant.index(':c.') + len(':c.')]
    #         tx_ac = tx_ac.replace(':c.', '')
    #         # Interface with the UTA database via get_tx_identity in uta.py
    #         tx_id_info = self.hdp.get_tx_identity_info(tx_ac)
    #         # NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
    #         return tx_id_info
    #
    #     # If the :n. pattern is present in the input variant
    #     if ':n.' in variant:
    #         # Remove all text to the right and including pat_c
    #         tx_ac = variant[:variant.index(':n.') + len(':n.')]
    #         tx_ac = tx_ac.replace(':n.', '')
    #         # Interface with the UTA database via get_tx_identity in uta.py
    #         tx_id_info = self.hdp.get_tx_identity_info(tx_ac)
    #         # NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
    #         return tx_id_info
    #
    #     # If the :r. pattern is present in the input variant
    #     if ':r.' in variant:
    #         # Remove all text to the right and including pat_c
    #         tx_ac = variant[:variant.index(':r.') + len(':r.')]
    #         tx_ac = tx_ac.replace(':r.', '')
    #         # Interface with the UTA database via get_tx_identity in uta.py
    #         tx_id_info = self.hdp.get_tx_identity_info(tx_ac)
    #         # NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
    #         return tx_id_info

    # def tx_id_info(self, alt_ac):
    #     """
    #     Input c. r. nd accession string
    #     Use uta.py (hdp) to return the identity information for the transcript variant
    #     see vvhgvs.dataproviders.uta.py for details
    #     """
    #     tx_id_info = self.hdp.get_tx_identity_info(alt_ac)
    #     # NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
    #     return tx_id_info

    # def tx_for_gene(self, hgnc):
    #     """
    #     Use uta.py (hdp) to return the transcript information for a specified gene (HGNC SYMBOL)
    #     see vvhgvs.dataproviders.uta.py for details
    #     """
    #     # Interface with the UTA database via get_tx_for_gene in uta.py
    #     tx_for_gene = self.hdp.get_tx_for_gene(hgnc)
    #     return tx_for_gene

    # def ng_extract(self, tx_for_gene):
    #     """
    #     Extract RefSeqGene Accession from transcript information
    #     see vvhgvs.dataproviders.uta.py for details
    #     """
    #     # For each list in the list of lists tx_for_gene
    #     for item in tx_for_gene:
    #         # If the pattern NG_ is found in element 4
    #         if 'NG_' in item[4]:
    #             # The gene accession is set to list element 4
    #             gene_ac = item[4]
    #             return gene_ac

    def tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        """
        Returns exon information for a given transcript
        e.g. how the exons align to the genomic reference
        see vvhgvs.dataproviders.uta.py for details
        """
        # Interface with the UTA database via get_tx_exons in uta.py
        try:
            tx_exons = self.hdp.get_tx_exons(tx_ac, alt_ac, alt_aln_method)
        except vvhgvs.exceptions.HGVSError as e:
            #e
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

    def relevant_transcripts(self, hgvs_genomic, evm, alt_aln_method, reverse_normalizer):
        """
        Automatically maps genomic positions onto all overlapping transcripts
        """
        # Pass relevant transcripts for the input variant to rts
        # Note, the evm method misses one end, the hdp. method misses the other. Combine both
        rts_list = self.hdp.get_tx_for_region(hgvs_genomic.ac, alt_aln_method,
                                              hgvs_genomic.posedit.pos.start.base-1,
                                              hgvs_genomic.posedit.pos.end.base-1)
        rts_dict = {}
        for tx_dat in rts_list:
            rts_dict[tx_dat[0]] = True
        rts_list_2 = evm.relevant_transcripts(hgvs_genomic)
        for tx_dat_2 in rts_list_2:
            rts_dict[tx_dat_2] = True
        rts = list(rts_dict.keys())

        # First if we have a ins prepare for hgvs "ins" mishandling, which
        # causes failures on any ins->non ins case, start by making a forced
        # "delins", equivilent to the vcf format ins requirements, then use
        # if needed. This is similar to the 'Triple check' code in mappers.py,
        # but more limited.
        hgvs_genomic_forced_delins = None
        if hgvs_genomic.posedit.edit.type == 'ins':
            start = hgvs_genomic.posedit.pos.start.base
            base = self.sf.fetch_seq(
                    str(hgvs_genomic.ac),start_i=start - 1, end_i=start)
            alt = base + hgvs_genomic.posedit.edit.alt
            hgvs_genomic_forced_delins = vvhgvs.sequencevariant.SequenceVariant(
                    ac=hgvs_genomic.ac,
                    type="g",
                    posedit=vvhgvs.posedit.PosEdit(
                        vvhgvs.location.Interval(
                            start=vvhgvs.location.SimplePosition(base=start),
                            end=vvhgvs.location.SimplePosition(base=start),
                            uncertain=hgvs_genomic.posedit.pos.uncertain
                            ),
                        vvhgvs.edit.NARefAlt(ref=base, alt=alt)
                        )
                    )
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
            except vvhgvs.exceptions.HGVSError:
                curr_genomic = hgvs_genomic
                if hgvs_genomic_forced_delins:
                    curr_genomic = hgvs_genomic_forced_delins
                    try:
                        variant = evm.g_to_t(hgvs_genomic_forced_delins, y)
                    except vvhgvs.exceptions.HGVSError:
                        pass
                # Check for non-coding transcripts
                try:
                    variant = evm.g_to_t(curr_genomic, y)
                except vvhgvs.exceptions.HGVSError:
                    continue
            except Exception as err:
                logger.warning('non expected err type', str(err))
                continue
            # Corrective Normalisation of intronic descriptions in the antisense oriemtation
            if '+' in str(variant) or '-' in str(variant) or '*' in str(variant):
                tx_ac = variant.ac
                alt_ac = hgvs_genomic.ac

                # Interface with the UTA database via get_tx_exons in uta.py
                try:
                    tx_exons = self.hdp.get_tx_exons(tx_ac, alt_ac, alt_aln_method)
                except vvhgvs.exceptions.HGVSError as e:
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

                # Gene orientation
                if tx_exons[0]['alt_strand'] == -1:
                    antisense = True
                else:
                    antisense = False

                # Pass if antisense = 'false'
                try:
                    if antisense:
                        # Reverse normalize hgvs_genomic
                        rev_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic)
                        # map back to coding
                        variant = evm.g_to_t(rev_hgvs_genomic, tx_ac)
                except vvhgvs.exceptions.HGVSInvalidIntervalError:
                    variant = evm.g_to_t(hgvs_genomic, tx_ac)
            try:
                self.hp.parse_hgvs_variant(str(variant))
            except vvhgvs.exceptions.HGVSError:
                continue
            except TypeError:
                continue
            else:
                code_var.append(variant)
        return code_var

    def validateHGVS(self, query):
        """
        Take HGVS string, parse into hgvs object and validate
        """
        hgvs_input = self.hp.parse_hgvs_variant(query)

        if ':p.' in query:
            if not hasattr(hgvs_input.posedit.pos.start, 'offset'):
                hgvs_input.posedit.pos.start.offset = 0
            if not hasattr(hgvs_input.posedit.pos.end, 'offset'):
                hgvs_input.posedit.pos.end.offset = 0
            if not hasattr(hgvs_input.posedit.pos.start, 'datum'):
                hgvs_input.posedit.pos.start.datum = 0
            if not hasattr(hgvs_input.posedit.pos.end, 'datum'):
                hgvs_input.posedit.pos.end.datum = 0
            if not hasattr(hgvs_input.posedit.edit, 'ref_n'):
                hgvs_input.posedit.edit.ref_n = hgvs_input.posedit.pos.end.base - hgvs_input.posedit.pos.start.base + 1

        try:
            self.vr.validate(hgvs_input)
        except vvhgvs.exceptions.HGVSError as e:
            return e
        else:
            return 'false'

    # def hgnc_rest(self, path):
    #     """
    #     Search HGNC rest
    #     """
    #     data = {
    #         'record': '',
    #         'error': 'false'
    #     }
    #     # HGNC server
    #     headers = {
    #         'Accept': 'application/json',
    #     }
    #     uri = 'http://rest.genenames.org'
    #     target = urlparse(uri + path)
    #     method = 'GET'
    #     body = ''
    #     h = http.Http()
    #     # collect the response
    #     response, content = h.request(
    #         target.geturl(),
    #         method,
    #         body,
    #         headers)
    #     if response['status'] == '200':
    #         # assume that content is a json reply
    #         # parse content with the json module
    #         data['record'] = json.loads(content)
    #     else:
    #         data['error'] = "Unable to contact the HGNC database: Please try again later"
    #     return data

    def entrez_efetch(self, db, id, rettype, retmode):
        """
        Search Entrez databases with efetch and SeqIO
        """
        # from Bio import Entrez
        Entrez.email = self.entrez_email
        Entrez.tool = 'VariantValidator'
        if self.entrez_api_key:
            Entrez.api_key = self.entrez_api_key
        # from Bio import SeqIO
        handle = Entrez.efetch(db=db, id=id, rettype=rettype, retmode=retmode)
        # Get record
        record = SeqIO.read(handle, "gb")
        # Place into text
        # text = handle.read()
        handle.close()
        return record

    # def entrez_read(self,db, id, retmode):
    #     """
    #     search Entrez databases with efetch and read
    #     """
    #     # IMPORT Bio modules
    #     # from Bio import Entrez
    #     Entrez.email = self.entrezID
    #     # from Bio import SeqIO
    #     handle = Entrez.efetch(db=db, id=id, retmode=retmode)
    #     # Get record
    #     record = Entrez.read(handle)
    #     # Place into text
    #     # text = handle.read()
    #     handle.close()
    #     return record

    def revcomp(self, bases):
        """
        Simple reverse complement function for nucleotide sequences
        """
        l2 = []
        listbases = list(bases)
        element = 0
        for base in listbases:
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

    def merge_hgvs_3pr(self, hgvs_variant_list, hn, genomic_reference=False, final_norm=True):
        """
        Function designed to merge multiple HGVS variants (hgvs objects) into a single delins
        using 3 prime normalization
        """
        # Ensure c. is mapped to the
        h_list = []
        # Sanity check and format the submitted variants
        for hgvs_v in hgvs_variant_list:
            # For testing include parser
            try:
                hgvs_v = self.hp.parse_hgvs_variant(hgvs_v)
            except Exception as e:
                logger.debug("Except passed, %s" % e)

            # Validate
            try:
                self.vr.validate(hgvs_v)  # Let hgvs errors deal with invalid variants and not hgvs objects
            except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                if 'Cannot validate sequence of an intronic variant' in str(e):
                    if genomic_reference is not False:
                        pass
                    else:
                        raise fn.mergeHGVSerror("Intronic variants can only be validated if a genomic/gene reference "
                                                "sequence is also provided e.g. NC_000017.11(NM_000088.3):c.589-1G>T")
                else:
                    self.vr.validate(hgvs_v)  # Let hgvs errors deal with invalid variants and not hgvs objects

            if hgvs_v.type == 'c':
                try:
                    hgvs_v = self.vm.c_to_n(hgvs_v)
                    h_list.append(hgvs_v)
                except:
                    raise fn.mergeHGVSerror("Unable to map from c. position to absolute position")
            elif hgvs_v.type == 'g':
                h_list.append(hgvs_v)

        if h_list:
            hgvs_variant_list = copy.deepcopy(h_list)

        # Define accession and start/end positions
        accession = None
        merge_start_pos = None
        merge_end_pos = None
        seqtype = None
        full_list = []

        # Loop through the submitted variants to remove any identity variants, these will be re-created as required
        # Except if it is the forst variant in which case we need the start position. We cannot assume non-gap
        elec = 0
        cp_hgvs_v = []
        for hgvs_v in hgvs_variant_list:
            if hgvs_v.posedit.edit.type == "identity":
                if elec == 0:
                    hgvs_v.posedit.pos.end.base = hgvs_v.posedit.pos.start.base
                    hgvs_v.posedit.edit.ref = hgvs_v.posedit.edit.ref[0]
                    hgvs_v.posedit.edit.alt = hgvs_v.posedit.edit.alt[0]
                    cp_hgvs_v.append(hgvs_v)
                continue
            else:
                cp_hgvs_v.append(hgvs_v)

        # Loop through the submitted variants and gather the required info
        hgvs_variant_list = cp_hgvs_v
        for hgvs_v in hgvs_variant_list:
            # No intronic positions
            try:
                if hgvs_v.posedit.pos.start.offset != 0 and genomic_reference is False:
                    raise fn.mergeHGVSerror("Base-offset position submitted")
                if hgvs_v.posedit.pos.end.offset != 0 and genomic_reference is False:
                    raise fn.mergeHGVSerror("Base-offset position submitted")
            except AttributeError as e:
                logger.debug("Except passed, %s", e)

            # Normalize the variant (allow cross intron) which also adds the reference sequence (?)
            try:
                hgvs_v = hn.normalize(hgvs_v)
            except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                pass

            # Set the accession and ensure that multiple reference sequences have not been queried
            if accession is None:
                accession = hgvs_v.ac
                seqtype = hgvs_v.type
            else:
                if hgvs_v.ac != accession:
                    raise fn.mergeHGVSerror("More than one reference sequence submitted")

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
                if hgvs_v.posedit.pos.start.base > merge_end_pos:
                    if hgvs_v.posedit.pos.start.base > merge_end_pos + 1:
                        # Create a fake variant to handle the missing sequence
                        ins_seq = self.sf.fetch_seq(hgvs_v.ac, merge_end_pos, hgvs_v.posedit.pos.start.base - 1)
                        gapping = hgvs_v.ac + ':' + hgvs_v.type + '.' + str(merge_end_pos + 1) + '_' + str(
                            hgvs_v.posedit.pos.start.base - 1) + 'del' + ins_seq + 'ins' + ins_seq
                        hgvs_gapping = self.hp.parse_hgvs_variant(gapping)
                        full_list.append(hgvs_gapping)
                    # update end_pos
                    merge_end_pos = hgvs_v.posedit.pos.end.base
                    # Append to the final list of variants
                    full_list.append(hgvs_v)
                else:
                    raise fn.mergeHGVSerror("Submitted variants are out of order or their ranges overlap")

        # Generate the alt sequence
        alt_sequence = ''
        for hgvs_v in full_list:
            ref_alt = hgvs_utils.hgvs_ref_alt(hgvs_v, self.sf)
            alt_sequence = alt_sequence + ref_alt['alt']

        # Fetch the reference sequence and copy it for the basis of the alt sequence
        reference_sequence = self.sf.fetch_seq(accession, merge_start_pos - 1, merge_end_pos)
        # Generate an hgvs_delins
        if alt_sequence == '':
            delins = accession + ':' + seqtype + '.' + str(merge_start_pos) + '_' + str(
                merge_end_pos) + 'del' + reference_sequence
        else:
            delins = accession + ':' + seqtype + '.' + str(merge_start_pos) + '_' + str(
                merge_end_pos) + 'del' + reference_sequence + 'ins' + alt_sequence
        hgvs_delins = self.hp.parse_hgvs_variant(delins)
        try:
            hgvs_delins = self.vm.n_to_c(hgvs_delins)
        except Exception as e:
            logger.debug("Except passed, %s", e)
        # Normalize (allow variants crossing into different exons)
        if final_norm is True:
            try:
                hgvs_delins = hn.normalize(hgvs_delins)
            except HGVSUnsupportedOperationError as e:
                logger.debug("Except passed, %s", e)
        return hgvs_delins


    # Code is being saved as it may be used in the future
    ######################################################
    # def merge_hgvs_5pr(self, hgvs_variant_list, genomic_reference=False):
    #     """
    #     Function designed to merge multiple HGVS variants (hgvs objects) into a single delins
    #     using 5 prime normalization
    #     """
    #     # Ensure c. is mapped to the
    #     h_list = []
    #
    #     # Sanity check and format the submitted variants
    #     for hgvs_v in hgvs_variant_list:
    #         # For testing include parser
    #         try:
    #             hgvs_v = self.hp.parse_hgvs_variant(hgvs_v)
    #         except Exception as e:
    #             logger.debug("Except passed, %s", e)
    #
    #         # Validate
    #         self.vr.validate(hgvs_v)  # Let hgvs errors deal with invalid variants and not hgvs objects
    #         if hgvs_v.type == 'c':
    #             try:
    #                 hgvs_v = self.vm.c_to_n(hgvs_v)
    #                 h_list.append(hgvs_v)
    #             except:
    #                 raise fn.mergeHGVSerror("Unable to map from c. position to absolute position")
    #     if h_list:
    #         hgvs_variant_list = copy.deepcopy(h_list)
    #
    #     # Define accession and start/end positions
    #     accession = None
    #     merge_start_pos = None
    #     merge_end_pos = None
    #     seqtype = None
    #     full_list = []
    #
    #     # Loop through the submitted variants and gather the required info
    #     for hgvs_v in hgvs_variant_list:
    #         try:
    #             # No intronic positions
    #             if hgvs_v.posedit.pos.start.offset != 0:
    #                 raise fn.mergeHGVSerror("Base-offset position submitted")
    #             if hgvs_v.posedit.pos.end.offset != 0:
    #                 raise fn.mergeHGVSerror("Base-offset position submitted")
    #         except AttributeError as e:
    #             logger.debug("Except passed, %s", e)
    #
    #         # Normalize the variant (allow cross intron) which also adds the reference sequence (?)
    #         hgvs_v = self.reverse_hn.normalize(hgvs_v)
    #
    #         # Set the accession and ensure that multiple reference sequences have not been queried
    #         if accession is None:
    #             accession = hgvs_v.ac
    #             seqtype = hgvs_v.type
    #         else:
    #             if hgvs_v.ac != accession:
    #                 raise fn.mergeHGVSerror("More than one reference sequence submitted")
    #
    #         # Set initial start and end positions
    #         if merge_start_pos is None:
    #             merge_start_pos = hgvs_v.posedit.pos.start.base
    #             merge_end_pos = hgvs_v.posedit.pos.end.base
    #             # Append to the final list of variants
    #             full_list.append(hgvs_v)
    #             continue
    #         # Ensure variants are in the correct order and not overlapping
    #         else:
    #             # ! hgvs_v.posedit.pos.start.base !>
    #             if hgvs_v.posedit.pos.start.base <= merge_end_pos:
    #                 raise fn.mergeHGVSerror("Submitted variants are out of order or their ranges overlap")
    #             else:
    #                 # Create a fake variant to handle the missing sequence
    #                 ins_seq = self.sf.fetch_seq(hgvs_v.ac, merge_end_pos, hgvs_v.posedit.pos.start.base - 1)
    #                 gapping = hgvs_v.ac + ':' + hgvs_v.type + '.' + str(merge_end_pos + 1) + '_' + str(
    #                     hgvs_v.posedit.pos.start.base - 1) + 'delins' + ins_seq
    #                 hgvs_gapping = self.hp.parse_hgvs_variant(gapping)
    #                 full_list.append(hgvs_gapping)
    #                 # update end_pos
    #                 merge_end_pos = hgvs_v.posedit.pos.end.base
    #                 # Append to the final list of variants
    #                 full_list.append(hgvs_v)
    #
    #     # Generate the alt sequence
    #     alt_sequence = ''
    #     for hgvs_v in full_list:
    #         ref_alt = hgvs_utils.hgvs_ref_alt(hgvs_v, self.sf)
    #         alt_sequence = alt_sequence + ref_alt['alt']
    #
    #     # Fetch the reference sequence and copy it for the basis of the alt sequence
    #     reference_sequence = self.sf.fetch_seq(accession, merge_start_pos - 1, merge_end_pos)
    #
    #     # Generate an hgvs_delins
    #     if alt_sequence == '':
    #         delins = accession + ':' + seqtype + '.' + str(merge_start_pos) + '_' + str(
    #             merge_end_pos) + 'del' + reference_sequence
    #     else:
    #         delins = accession + ':' + seqtype + '.' + str(merge_start_pos) + '_' + str(
    #             merge_end_pos) + 'del' + reference_sequence + 'ins' + alt_sequence
    #     hgvs_delins = self.hp.parse_hgvs_variant(delins)
    #     try:
    #         hgvs_delins = self.vm.n_to_c(hgvs_delins)
    #     except Exception as e:
    #         logger.debug("Except passed, %s", e)
    #     # Normalize (allow variants crossing into different exons)
    #     try:
    #         hgvs_delins = self.reverse_hn.normalize(hgvs_delins)
    #     except HGVSUnsupportedOperationError as e:
    #         logger.debug("Except passed, %s", e)
    #    return hgvs_delins

    # def merge_pseudo_vcf(self, vcf_list, genome_build, hn):
    #     """
    #     Function designed to merge multiple pseudo VCF variants (strings) into a single HGVS delins
    #     using 5 prime normalization then return a 3 prime normalized final HGVS object
    #     """
    #     hgvs_list = []
    #     # Convert pseudo_vcf list into a HGVS list
    #     for call in vcf_list:
    #         x55hgvs = hgvs_utils.pvcf_to_hgvs(call, genome_build, normalization_direction=5, validator=self)
    #         hgvs_list.append(x55hgvs)
    #     # Merge
    #     hgvs_delins = self.merge_hgvs_5pr(hgvs_list)
    #     # normalize 3 prime
    #     hgvs_delins = hn.normalize(hgvs_delins)
    #     # return
    #     return hgvs_delins

    def hgvs_alleles(self, variant_description, hn, genomic_reference=False):
        """
        HGVS allele handling function which takes a single HGVS allele description and
        separates each allele into a list of HGVS variants
        """
        try:
            # Split up the description
            accession, remainder = variant_description.split(':')
            # Branch
            if re.search(r'[gcn]\.\d+\[', remainder):
                # NM_004006.2:c.2376[G>C];[(G>C)]
                # if re.search('\(', remainder):
                #   raise fn.alleleVariantError('Unsupported format ' + remainder)
                # NM_004006.2:c.2376[G>C];[G>C]
                type, remainder = remainder.split('.')
                pos = re.match(r'\d+', remainder)
                pos = pos.group(0)
                remainder = remainder.replace(pos, '')
                remainder = remainder[1:-1]
                alleles = remainder.split('];[')
                my_alleles = []
                for posedit in alleles:
                    if '(' in posedit:
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
                if '(;)' in remainder and '];' in remainder:
                    # NM_004006.2:c.[296T>G];[476T>C](;)1083A>C(;)1406del
                    pre_alleles = remainder.split('(;)')
                    pre_merges = []
                    alleles = []
                    for allele in pre_alleles:
                        if '[' in allele:
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
                        if '?' in str(each_allele):
                            # NM_004006.2:c.[2376G>C];[?]
                            continue
                        merge = []
                        allele = str(self.merge_hgvs_3pr(each_allele, hn, genomic_reference))
                        merge.append(allele)
                        for variant in each_allele:
                            merged_alleles.append([variant])
                    my_alleles = merged_alleles

                elif '(;)' in remainder:
                    # If statement for uncertainties
                    # NM_004006.2:c.[296T>G;476C>T];[476C>T](;)1083A>C
                    if '[' in remainder:
                        raise fn.alleleVariantError('Unsupported format ' + type + '.' + remainder)
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
                    if '(' in remainder:
                        raise fn.alleleVariantError('Unsupported format ' + type + '.' + remainder)
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
                        if '?' in str(each_allele):
                            # NM_004006.2:c.[2376G>C];[?]
                            continue
                        merge = []
                        allele = str(self.merge_hgvs_3pr(each_allele, hn, genomic_reference))
                        merge.append(allele)
                        for variant in each_allele:
                            merged_alleles.append([variant])
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
            exc_type, exc_value, last_traceback = sys.exc_info()
            logger.error(str(exc_type) + " " + str(exc_value))
            traceback.print_tb(last_traceback, file=sys.stdout)
            raise fn.alleleVariantError(str(e))

    def chr_to_rsg(self, hgvs_genomic, hn):
        """
        # Covert chromosomal HGVS description to RefSeqGene
        """
        # 'chr_to_rsg triggered'
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
        all_info = self.db.get_g_to_g_info()
        for line in all_info:
            # Logic to identify the correct RefSeqGene
            rsg_data = {}
            if chr_ac == line[1] and chr_start_pos >= int(line[2]) and chr_end_pos <= int(line[3]):
                # query = "SELECT refSeqGeneID, refSeqChromosomeID, startPos, endPos, orientation, hgncSymbol FROM
                # refSeqGene_loci"
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
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'does not agree with reference sequence' in error:
                        match = re.findall(r'\(([GATC]+)\)', error)
                        new_ref = match[1]
                        hgvs_refseqgene.posedit.edit.ref = new_ref
                        error = 'true'

                    data = {'hgvs_refseqgene': str(hgvs_refseqgene), 'gene': gene, 'valid': str(error)}
                else:
                    data = {'hgvs_refseqgene': str(hgvs_refseqgene), 'gene': gene, 'valid': 'true'}
                descriptions.append(data)
            if ori == '-':
                # Reverse complement of bases may be required. Let normalizer do the lifting for strings of bases
                # Look for scenarios with RC needed bases and extract the bases from the edit
                if re.search(r"(del[GATCUgatcu]+)", str(chr_edit)):
                    bases = re.search(r"(del[GATCUgatcu]+)", str(chr_edit))
                    bases = bases.group(1)
                    chr_edit = 'del' + str(chr_edit).replace(bases, '')
                if re.search(r"(ins[GATCUgatcu]+)", str(chr_edit)):
                    bases = re.search(r"(ins[GATCUgatcu]+)", str(chr_edit))
                    bases = bases.group(1)
                    ins_revcomp = self.revcomp(bases)
                    chr_edit = str(chr_edit).replace(bases, '') + 'ins' + ins_revcomp
                if re.search(r"(dup[GATCUgatcu]+)", str(chr_edit)):
                    bases = re.search(r"(dup[GATCUgatcu]+)", str(chr_edit))
                    bases = bases.group(1)
                    chr_edit = 'dup' + str(chr_edit).replace(bases, '')
                if re.search(r"(inv[GATCUgatcu]+)", str(chr_edit)):
                    bases = re.search(r"(inv[GATCUgatcu]+)", str(chr_edit))
                    bases = bases.group(1)
                    chr_edit = 'inv' + str(chr_edit).replace(bases, '')
                if '>' in str(chr_edit) or '=' in str(chr_edit):
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
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'does not agree with reference sequence' in error:
                        match = re.findall(r'\(([GATC]+)\)', error)
                        new_ref = match[1]
                        hgvs_refseqgene.posedit.edit.ref = new_ref
                        error = 'true'

                    data = {'hgvs_refseqgene': str(hgvs_refseqgene), 'gene': gene, 'valid': str(error)}
                else:
                    data = {'hgvs_refseqgene': str(hgvs_refseqgene), 'gene': gene, 'valid': 'true'}
                descriptions.append(data)

        # Return the required data. This is a dictionary containing the rsg description, validation status and gene ID
        return descriptions

    def rsg_to_chr(self, hgvs_refseqgene, primary_assembly, hn):
        """
        # Covert RefSeqGene HGVS description to Chromosomal

        :param hgvs_refseqgene:
        :param primary_assembly:
        :param hn: HGVS Normalizer
        :param vr:
        :return:
        """
        # normalize
        try:
            hgvs_refseqgene = hn.normalize(hgvs_refseqgene)
        except Exception as e:
            logger.debug("Except passed, %s", e)
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
        all_info = self.db.get_g_to_g_info()
        for line in all_info:
            # Logic to identify the correct RefSeqGene
            chr_data = {}
            if rsg_ac == line[0] and primary_assembly == line[6]:
                # query = "SELECT refSeqGeneID, refSeqChromosomeID, startPos, endPos, orientation,
                # hgncSymbol FROM refSeqGene_loci"
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
                    self.vr.validate(hgvs_genomic)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'does not agree with reference sequence' in error:
                        match = re.findall(r'\(([GATC]+)\)', error)
                        new_ref = match[1]
                        hgvs_genomic.posedit.edit.ref = new_ref
                        error = 'true'
                    data = {'hgvs_genomic': str(hgvs_genomic), 'gene': gene, 'valid': str(error)}
                else:
                    data = {'hgvs_genomic': str(hgvs_genomic), 'gene': gene, 'valid': 'true'}
                descriptions.append(data)
            if ori == '-':
                # Reverse complement of bases may be required. Let normalizer do the lifting for strings of bases
                # Look for scenarios with RC needed bases and extract the bases from the edit
                if re.search(r'(del[GATCUgatcu]+)', str(rsg_edit)):
                    bases = re.search(r"(del[GATCUgatcu]+)", str(rsg_edit))
                    bases = bases.group(1)
                    rsg_edit = 'del' + str(rsg_edit).replace(bases, '')
                if re.search(r"(ins[GATCUgatcu]+)", str(rsg_edit)):
                    bases = re.search(r"(ins[GATCUgatcu]+)", str(rsg_edit))
                    bases = bases.group(1)
                    ins_revcomp = self.revcomp(bases)
                    rsg_edit = str(rsg_edit).replace(bases, '') + 'ins' + ins_revcomp
                if re.search(r"(dup[GATCUgatcu]+)", str(rsg_edit)):
                    bases = re.search(r"(dup[GATCUgatcu]+)", str(rsg_edit))
                    bases = bases.group(1)
                    rsg_edit = 'dup' + str(rsg_edit).replace(bases, '')
                if re.search(r"(inv[GATCUgatcu]+)", str(rsg_edit)):
                    bases = re.search(r"(inv[GATCUgatcu]+)", str(rsg_edit))
                    bases = bases.group(1)
                    rsg_edit = 'inv' + str(rsg_edit).replace(bases, '')
                if '>' in str(rsg_edit) or '=' in str(rsg_edit):
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
                    self.vr.validate(hgvs_genomic)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'does not agree with reference sequence' in error:
                        match = re.findall(r'\(([GATC]+)\)', error)
                        new_ref = match[1]
                        hgvs_genomic.posedit.edit.ref = new_ref
                        error = 'true'
                    data = {'hgvs_genomic': str(hgvs_genomic), 'gene': gene, 'valid': str(error)}
                else:
                    data = {'hgvs_genomic': str(hgvs_genomic), 'gene': gene, 'valid': 'true'}
                descriptions.append(data)

        # Return the required data. This is a dictionary containing the rsg description, validation status and gene ID
        return descriptions

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
