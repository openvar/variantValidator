import hgvs
import re
import copy
from . import vvHGVS
from . import vvFunctions as fn
from .vvLogging import logger


def gapped_g_to_c(variant, validator, rel_var):
    """
    Gap aware projection from g. to c.
    """

    # Set variables for problem specific warnings
    gapped_alignment_warning = ''
    corrective_action_taken = ''
    gapped_transcripts = ''
    auto_info = ''
    disparity_deletion_in = []

    # Create a pseudo VCF so that normalization can be applied and a delins can be generated
    hgvs_genomic_variant = variant.hgvs_genomic
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = variant.reverse_normalizer.normalize(hgvs_genomic_variant)
    hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

    # VCF
    vcf_dict = vvHGVS.hgvs2vcf(reverse_normalized_hgvs_genomic, variant.primary_assembly,
                               variant.reverse_normalizer, validator.sf)
    chr = vcf_dict['chr']
    pos = vcf_dict['pos']
    ref = vcf_dict['ref']
    alt = vcf_dict['alt']

    # Generate an end position
    end = str(int(pos) + len(ref) - 1)
    pos = str(pos)

    # take a look at the input genomic variant for potential base salvage
    stash_ac = vcf_dict['chr']
    stash_pos = int(vcf_dict['pos'])
    stash_ref = vcf_dict['ref']
    stash_alt = vcf_dict['alt']
    stash_end = end
    stash_input = str(variant.stashed)
    # Re-Analyse genomic positions
    if 'NG_' in str(variant.hgvs_formatted):
        c = validator.hp.parse_hgvs_variant(rel_var[0])
        if hasattr(c.posedit.edit, 'ref') and c.posedit.edit.ref is not None:
            c.posedit.edit.ref = c.posedit.edit.ref.upper()
        if hasattr(c.posedit.edit, 'alt') and c.posedit.edit.alt is not None:
            c.posedit.edit.alt = c.posedit.edit.alt.upper()
        stash_input = validator.myevm_t_to_g(c, variant.no_norm_evm, variant.primary_assembly, variant.hn)
    if re.match('NC_', str(stash_input)) or re.match('NT_', str(stash_input)) or re.match('NW_',
                                                                                          str(
                                                                                              stash_input)):
        try:
            hgvs_stash = validator.hp.parse_hgvs_variant(stash_input)
        except:
            hgvs_stash = stash_input
        if hasattr(hgvs_stash.posedit.edit, 'ref') and hgvs_stash.posedit.edit.ref is not None:
            hgvs_stash.posedit.edit.ref = hgvs_stash.posedit.edit.ref.upper()
        if hasattr(hgvs_stash.posedit.edit, 'alt') and hgvs_stash.posedit.edit.alt is not None:
            hgvs_stash.posedit.edit.alt = hgvs_stash.posedit.edit.alt.upper()

        stash_ac = hgvs_stash.ac
        # MAKE A NO NORM HGVS2VCF
        stash_dict = vvHGVS.pos_lock_hgvs2vcf(hgvs_stash, variant.primary_assembly, variant.reverse_normalizer, validator.sf)
        stash_ac = hgvs_stash.ac
        stash_pos = int(stash_dict['pos'])
        stash_ref = stash_dict['ref']
        stash_alt = stash_dict['alt']
        # Generate an end position
        stash_end = str(stash_pos + len(stash_ref) - 1)

    # Store a not real deletion insertion
    stored_hgvs_not_delins = validator.hp.parse_hgvs_variant(str(
        hgvs_genomic_5pr.ac) + ':' + hgvs_genomic_5pr.type + '.' + pos + '_' + end + 'del' + ref + 'ins' + alt)
    stash_hgvs_not_delins = validator.hp.parse_hgvs_variant(
        stash_ac + ':' + hgvs_genomic_5pr.type + '.' + str(
            stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)

    # Set non-valid caution to false
    non_valid_caution = 'false'

    # make an empty rel_var
    nw_rel_var = []

    # loop through rel_var and amend where required
    for var in rel_var:
        # Store the current hgvs:c. description
        saved_hgvs_coding = validator.hp.parse_hgvs_variant(var)

        # Remove un-selected transcripts
        if validator.select_transcripts != 'all':
            tx_ac = saved_hgvs_coding.ac
            # If it's in the selected tx dict, keep it
            if tx_ac.split('.')[0] in list(validator.select_transcripts_dict.keys()):
                pass
            # If not get rid of it!
            else:
                continue

        # Get orientation of the gene wrt genome and a list of exons mapped to the genome
        ori = validator.tx_exons(tx_ac=saved_hgvs_coding.ac, alt_ac=hgvs_genomic_5pr.ac,
                            alt_aln_method=validator.alt_aln_method)
        orientation = int(ori[0]['alt_strand'])
        intronic_variant = 'false'

        if orientation == -1:
            # position genomic at its most 5 prime position
            try:
                query_genomic = variant.reverse_normalizer.normalize(variant.hgvs_genomic)
            except:
                query_genomic = variant.hgvs_genomic
            # Map to the transcript ant test for movement
            try:
                hgvs_seek_var = variant.evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
            except hgvs.exceptions.HGVSError as e:
                hgvs_seek_var = saved_hgvs_coding
            else:
                seek_var = fn.valstr(hgvs_seek_var)
                seek_ac = str(hgvs_seek_var.ac)
            if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                    saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                    hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                    saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset):
                pass
            else:
                hgvs_seek_var = saved_hgvs_coding

        elif orientation != -1:
            # position genomic at its most 3 prime position
            try:
                query_genomic = variant.hn.normalize(variant.hgvs_genomic)
            except:
                query_genomic = variant.hgvs_genomic
        # Map to the transcript and test for movement
        try:
            hgvs_seek_var = variant.evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
        except hgvs.exceptions.HGVSError as e:
            hgvs_seek_var = saved_hgvs_coding
        else:
            seek_var = fn.valstr(hgvs_seek_var)
            seek_ac = str(hgvs_seek_var.ac)
            if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                    saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                    hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                    saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset):
                pass
            else:
                hgvs_seek_var = saved_hgvs_coding

        try:
            intron_test = variant.hn.normalize(hgvs_seek_var)
        except hgvs.exceptions.HGVSUnsupportedOperationError as e:
            error = str(e)
            if re.match('Normalization of intronic variants is not supported', error) or re.match(
                    'Unsupported normalization of variants spanning the exon-intron boundary',
                    error):
                if re.match(
                        'Unsupported normalization of variants spanning the exon-intron boundary',
                        error):
                    intronic_variant = 'hard_fail'
                else:
                    # Double check to see whether the variant is actually intronic?
                    for exon in ori:
                        genomic_start = int(exon['alt_start_i'])
                        genomic_end = int(exon['alt_end_i'])
                        if (
                                hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                                hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                            intronic_variant = 'false'
                            break
                        else:
                            intronic_variant = 'true'

        if intronic_variant != 'hard_fail':
            if re.search(r'\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search(r'\d+\-', str(
                    hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+\+', str(
                hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+\-',
                                                         str(hgvs_seek_var.posedit.pos)):
                # Double check to see whether the variant is actually intronic?
                for exon in ori:
                    genomic_start = int(exon['alt_start_i'])
                    genomic_end = int(exon['alt_end_i'])
                    if (
                            hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                            hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                        intronic_variant = 'false'
                        break
                    else:
                        intronic_variant = 'true'

        if re.search(r'\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search(r'\d+\-', str(
                hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+\+', str(
            hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+\-', str(hgvs_seek_var.posedit.pos)):
            # Double check to see whether the variant is actually intronic?
            for exon in ori:
                genomic_start = int(exon['alt_start_i'])
                genomic_end = int(exon['alt_end_i'])
                if (
                        hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                        hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                    intronic_variant = 'false'
                    break
                else:
                    intronic_variant = 'true'

        # If exonic, process
        if intronic_variant != 'true':
            # map form reverse normalized g. to c.
            hgvs_from_5n_g = variant.no_norm_evm.g_to_t(hgvs_genomic_5pr, saved_hgvs_coding.ac)

            # Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping transcript lengths
            disparity_deletion_in = ['false', 'false']
            if stored_hgvs_not_delins != '':
                # Refresh hgvs_not_delins from stored_hgvs_not_delins
                hgvs_not_delins = copy.deepcopy(stored_hgvs_not_delins)
                # This test will only occur in dup of single base, insertion or substitution
                if not re.search('_', str(hgvs_not_delins.posedit.pos)):
                    if re.search('dup', hgvs_genomic_5pr.posedit.edit.type) or re.search('ins',
                                                                                         hgvs_genomic_5pr.posedit.edit.type):
                        # For gap in chr, map to t. - but becaouse we have pushed to 5 prime by norm, add 1 to end pos
                        plussed_hgvs_not_delins = copy.deepcopy(hgvs_not_delins)
                        plussed_hgvs_not_delins.posedit.pos.end.base = plussed_hgvs_not_delins.posedit.pos.end.base + 1
                        plussed_hgvs_not_delins.posedit.edit.ref = ''
                        transcript_variant = variant.no_norm_evm.g_to_t(plussed_hgvs_not_delins,
                                                                str(saved_hgvs_coding.ac))
                        if ((
                                transcript_variant.posedit.pos.end.base - transcript_variant.posedit.pos.start.base) > (
                                hgvs_genomic_5pr.posedit.pos.end.base - hgvs_genomic_5pr.posedit.pos.start.base)):
                            if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                                hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                start = hgvs_not_delins.posedit.pos.start.base - 1
                                end = hgvs_not_delins.posedit.pos.end.base
                                ref_bases = validator.sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                hgvs_not_delins.posedit.edit.ref = ref_bases
                                hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                   :1] + hgvs_not_delins.posedit.edit.alt[
                                                                         1:] + ref_bases[1:]
                            elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search(
                                    'del', str(hgvs_genomic_5pr.posedit.edit)):
                                hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                            elif re.search('ins',
                                           str(hgvs_genomic_5pr.posedit.edit)) and not re.search(
                                'del', str(hgvs_genomic_5pr.posedit.edit)):
                                hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                start = hgvs_not_delins.posedit.pos.start.base - 1
                                end = hgvs_not_delins.posedit.pos.end.base
                                ref_bases = validator.sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                hgvs_not_delins.posedit.edit.ref = ref_bases
                                hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                   :1] + hgvs_not_delins.posedit.edit.alt[
                                                                         1:] + ref_bases[1:]
                        else:
                            if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                                hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                start = hgvs_not_delins.posedit.pos.start.base - 1
                                end = hgvs_not_delins.posedit.pos.end.base
                                ref_bases = validator.sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                hgvs_not_delins.posedit.edit.ref = ref_bases
                                hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                   :1] + hgvs_not_delins.posedit.edit.alt[
                                                                         1:] + ref_bases[1:]
                            elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search(
                                    'del', str(hgvs_genomic_5pr.posedit.edit)):
                                hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                            elif re.search('ins',
                                           str(hgvs_genomic_5pr.posedit.edit)) and not re.search(
                                'del', str(hgvs_genomic_5pr.posedit.edit)):
                                hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                start = hgvs_not_delins.posedit.pos.start.base - 1
                                end = hgvs_not_delins.posedit.pos.end.base
                                ref_bases = validator.sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                hgvs_not_delins.posedit.edit.ref = ref_bases
                                hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                   :1] + hgvs_not_delins.posedit.edit.alt[
                                                                         1:] + ref_bases[1:]
                    else:
                        pass
                else:
                    pass

                try:
                    tx_hgvs_not_delins = variant.no_norm_evm.g_to_n(hgvs_not_delins, saved_hgvs_coding.ac)
                except hgvs.exceptions.HGVSInvalidIntervalError:
                    tx_hgvs_not_delins = variant.no_norm_evm.g_to_n(hgvs_genomic_5pr, saved_hgvs_coding.ac)
                except hgvs.exceptions.HGVSError:
                    if str(e) == 'start or end or both are beyond the bounds of transcript record':
                        tx_hgvs_not_delins = saved_hgvs_coding

                # Create normalized version of tx_hgvs_not_delins
                rn_tx_hgvs_not_delins = copy.deepcopy(tx_hgvs_not_delins)
                # Check for +ve base and adjust
                if (re.search(r'\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)) or re.search(r'\-',
                                                                                                str(
                                                                                                    rn_tx_hgvs_not_delins.posedit.pos.start))) and (
                        re.search(r'\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)) or re.search(
                    r'\-', str(rn_tx_hgvs_not_delins.posedit.pos.end))):
                    # Remove offsetting to span the gap
                    rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                    rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                    rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
                    rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                    try:
                        rn_tx_hgvs_not_delins.posedit.edit.alt = ''
                    except:
                        fn.exceptPass()
                elif re.search(r'\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                    # move tx end base to next available non-offset base
                    rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
                    rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                    rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                    if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                        test_tx_var = variant.no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                    else:
                        test_tx_var = rn_tx_hgvs_not_delins
                    # re-make genomic and tx
                    hgvs_not_delins = validator.myevm_t_to_g(test_tx_var, variant.no_norm_evm,
                                                        variant.primary_assembly, variant.hn)
                    rn_tx_hgvs_not_delins = variant.no_norm_evm.g_to_n(hgvs_not_delins,
                                                               str(saved_hgvs_coding.ac))
                elif re.search(r'\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                    # move tx start base to previous available non-offset base
                    rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                    rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                    if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                        test_tx_var = variant.no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                    else:
                        test_tx_var = rn_tx_hgvs_not_delins
                    # re-make genomic and tx
                    hgvs_not_delins = validator.myevm_t_to_g(test_tx_var, variant.no_norm_evm,
                                                        variant.primary_assembly, variant.hn)
                    rn_tx_hgvs_not_delins = variant.no_norm_evm.g_to_n(hgvs_not_delins,
                                                               str(saved_hgvs_coding.ac))
                    rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                #                                     else:
                #                                         pass

                # Check for -ve base and adjust
                elif re.search(r'\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search(r'\-',
                                                                                                str(
                                                                                                    rn_tx_hgvs_not_delins.posedit.pos.start)):
                    # Remove offsetting to span the gap
                    rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                    rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                    rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
                    rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                    try:
                        rn_tx_hgvs_not_delins.posedit.edit.alt = ''
                    except:
                        fn.exceptPass()
                elif re.search(r'\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                    # move tx end base back to next available non-offset base
                    rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                    # Delete the ref
                    rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                    # Add the additional base to the ALT
                    start = rn_tx_hgvs_not_delins.posedit.pos.end.base - 1
                    end = rn_tx_hgvs_not_delins.posedit.pos.end.base
                    ref_bases = validator.sf.fetch_seq(str(tx_hgvs_not_delins.ac), start, end)
                    rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + ref_bases
                    if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                        test_tx_var = variant.no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                    else:
                        test_tx_var = rn_tx_hgvs_not_delins
                    # re-make genomic and tx
                    hgvs_not_delins = validator.myevm_t_to_g(test_tx_var, variant.no_norm_evm,
                                                        variant.primary_assembly, variant.hn)
                    rn_tx_hgvs_not_delins = variant.no_norm_evm.g_to_n(hgvs_not_delins,
                                                               str(saved_hgvs_coding.ac))
                elif re.search(r'\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                    # move tx start base to previous available non-offset base
                    rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                    rn_tx_hgvs_not_delins.posedit.pos.start.base = rn_tx_hgvs_not_delins.posedit.pos.start.base - 1
                    rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                    if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                        test_tx_var = variant.no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                    else:
                        test_tx_var = rn_tx_hgvs_not_delins
                    # re-make genomic and tx
                    hgvs_not_delins = validator.myevm_t_to_g(test_tx_var, variant.no_norm_evm,
                                                        variant.primary_assembly, variant.hn)
                    rn_tx_hgvs_not_delins = variant.no_norm_evm.g_to_n(hgvs_not_delins,
                                                               str(saved_hgvs_coding.ac))
                    rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                else:
                    pass

                # Logic
                if len(hgvs_not_delins.posedit.edit.ref) < len(
                        rn_tx_hgvs_not_delins.posedit.edit.ref):
                    gap_length = len(rn_tx_hgvs_not_delins.posedit.edit.ref) - len(
                        hgvs_not_delins.posedit.edit.ref)
                    disparity_deletion_in = ['chromosome', gap_length]
                elif len(hgvs_not_delins.posedit.edit.ref) > len(
                        rn_tx_hgvs_not_delins.posedit.edit.ref):
                    gap_length = len(hgvs_not_delins.posedit.edit.ref) - len(
                        rn_tx_hgvs_not_delins.posedit.edit.ref)
                    disparity_deletion_in = ['transcript', gap_length]
                else:
                    hgvs_stash_t = validator.vm.g_to_t(stash_hgvs_not_delins, saved_hgvs_coding.ac)
                    if len(stash_hgvs_not_delins.posedit.edit.ref) > len(
                            hgvs_stash_t.posedit.edit.ref):
                        try:
                            variant.hn.normalize(hgvs_stash_t)
                        except:
                            fn.exceptPass()
                        else:
                            gap_length = len(stash_hgvs_not_delins.posedit.edit.ref) - len(
                                hgvs_stash_t.posedit.edit.ref)
                            disparity_deletion_in = ['transcript', gap_length]
                            try:
                                tx_hgvs_not_delins = validator.vm.c_to_n(hgvs_stash_t)
                            except:
                                tx_hgvs_not_delins = hgvs_stash_t
                            hgvs_not_delins = stash_hgvs_not_delins
                    elif hgvs_stash_t.posedit.pos.start.offset != 0 or hgvs_stash_t.posedit.pos.end.offset != 0:
                        disparity_deletion_in = ['transcript', 'Requires Analysis']
                        try:
                            tx_hgvs_not_delins = validator.vm.c_to_n(hgvs_stash_t)
                        except:
                            tx_hgvs_not_delins = hgvs_stash_t
                        hgvs_not_delins = stash_hgvs_not_delins
                        hgvs_genomic_5pr = stash_hgvs_not_delins
                    else:
                        pass

            # Final sanity checks
            try:
                validator.vm.g_to_t(hgvs_not_delins, tx_hgvs_not_delins.ac)
            except Exception as e:
                if str(e) == 'start or end or both are beyond the bounds of transcript record':
                    hgvs_not_delins = saved_hgvs_coding
                    disparity_deletion_in = ['false', 'false']
                logger.warning(str(e))
            try:
                variant.hn.normalize(tx_hgvs_not_delins)
            except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                error = str(e)
                if re.match('Normalization of intronic variants is not supported',
                            error) or re.match(
                    'Unsupported normalization of variants spanning the exon-intron boundary',
                    error):
                    if re.match(
                            'Unsupported normalization of variants spanning the exon-intron boundary',
                            error):
                        hgvs_not_delins = saved_hgvs_coding
                        disparity_deletion_in = ['false', 'false']
                    elif re.match('Normalization of intronic variants is not supported', error):
                        # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                        disparity_deletion_in = ['transcript', 'Requires Analysis']
                logger.warning(error)
            # Pre-processing of tx_hgvs_not_delins
            try:
                if tx_hgvs_not_delins.posedit.edit.alt is None:
                    tx_hgvs_not_delins.posedit.edit.alt = ''
            except Exception as e:
                if str(e) == "'Dup' object has no attribute 'alt'":
                    tx_hgvs_not_delins_delins_from_dup = tx_hgvs_not_delins.ac + ':' + tx_hgvs_not_delins.type + '.' + str(
                        tx_hgvs_not_delins.posedit.pos.start) + '_' + str(
                        tx_hgvs_not_delins.posedit.pos.end) + 'del' + tx_hgvs_not_delins.posedit.edit.ref + 'ins' + tx_hgvs_not_delins.posedit.edit.ref + tx_hgvs_not_delins.posedit.edit.ref
                    tx_hgvs_not_delins = validator.hp.parse_hgvs_variant(tx_hgvs_not_delins_delins_from_dup)

            # GAP IN THE TRANSCRIPT DISPARITY DETECTED
            if disparity_deletion_in[0] == 'transcript':
                gap_position = ''
                gapped_alignment_warning = str(
                    hgvs_genomic_5pr) + ' does not represent a true variant because it is an artefact of aligning the transcripts listed below with genome build ' + variant.primary_assembly

                # ANY VARIANT WHOLLY WITHIN THE GAP
                if (re.search(r'\+', str(tx_hgvs_not_delins.posedit.pos.start)) or re.search(r'\-',
                                                                                             str(
                                                                                                 tx_hgvs_not_delins.posedit.pos.start))) and (
                        re.search(r'\+', str(tx_hgvs_not_delins.posedit.pos.end)) or re.search(r'\-',
                                                                                               str(
                                                                                                   tx_hgvs_not_delins.posedit.pos.end))):
                    gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
                    # Copy the current variant
                    tx_gap_fill_variant = copy.deepcopy(tx_hgvs_not_delins)
                    try:
                        if tx_gap_fill_variant.posedit.edit.alt is None:
                            tx_gap_fill_variant.posedit.edit.alt = ''
                    except Exception as e:
                        if str(e) == "'Dup' object has no attribute 'alt'":
                            tx_gap_fill_variant_delins_from_dup = tx_gap_fill_variant.ac + ':' + tx_gap_fill_variant.type + '.' + str(
                                tx_gap_fill_variant.posedit.pos.start) + '_' + str(
                                tx_gap_fill_variant.posedit.pos.end) + 'del' + tx_gap_fill_variant.posedit.edit.ref + 'ins' + tx_gap_fill_variant.posedit.edit.ref + tx_gap_fill_variant.posedit.edit.ref
                            tx_gap_fill_variant = validator.hp.parse_hgvs_variant(
                                tx_gap_fill_variant_delins_from_dup)

                    # Identify which half of the NOT-intron the start position of the variant is in
                    if re.search(r'\-', str(tx_gap_fill_variant.posedit.pos.start)):
                        tx_gap_fill_variant.posedit.pos.start.base = tx_gap_fill_variant.posedit.pos.start.base - 1
                        tx_gap_fill_variant.posedit.pos.start.offset = int('0')  # int('+1')
                        tx_gap_fill_variant.posedit.pos.end.offset = int('0')  # int('-1')
                        tx_gap_fill_variant.posedit.edit.alt = ''
                        tx_gap_fill_variant.posedit.edit.ref = ''
                    elif re.search(r'\+', str(tx_gap_fill_variant.posedit.pos.start)):
                        tx_gap_fill_variant.posedit.pos.start.offset = int('0')  # int('+1')
                        tx_gap_fill_variant.posedit.pos.end.base = tx_gap_fill_variant.posedit.pos.end.base + 1
                        tx_gap_fill_variant.posedit.pos.end.offset = int('0')  # int('-1')
                        tx_gap_fill_variant.posedit.edit.alt = ''
                        tx_gap_fill_variant.posedit.edit.ref = ''

                    try:
                        tx_gap_fill_variant = validator.vm.n_to_c(tx_gap_fill_variant)
                    except:
                        fn.exceptPass()
                    genomic_gap_fill_variant = validator.vm.t_to_g(tx_gap_fill_variant,
                                                              reverse_normalized_hgvs_genomic.ac)
                    genomic_gap_fill_variant.posedit.edit.alt = genomic_gap_fill_variant.posedit.edit.ref

                    try:
                        c_tx_hgvs_not_delins = validator.vm.n_to_c(tx_hgvs_not_delins)
                    except Exception:
                        c_tx_hgvs_not_delins = copy.copy(tx_hgvs_not_delins)
                    genomic_gap_fill_variant_alt = validator.vm.t_to_g(c_tx_hgvs_not_delins,
                                                                  hgvs_genomic_5pr.ac)

                    # Ensure an ALT exists
                    try:
                        if genomic_gap_fill_variant_alt.posedit.edit.alt is None:
                            genomic_gap_fill_variant_alt.posedit.edit.alt = 'X'
                    except Exception as e:
                        if str(e) == "'Dup' object has no attribute 'alt'":
                            genomic_gap_fill_variant_delins_from_dup = genomic_gap_fill_variant.ac + ':' + genomic_gap_fill_variant.type + '.' + str(
                                genomic_gap_fill_variant.posedit.pos.start.base) + '_' + str(
                                genomic_gap_fill_variant.posedit.pos.end.base) + 'del' + genomic_gap_fill_variant.posedit.edit.ref + 'ins' + genomic_gap_fill_variant.posedit.edit.ref + genomic_gap_fill_variant.posedit.edit.ref
                            genomic_gap_fill_variant = validator.hp.parse_hgvs_variant(
                                genomic_gap_fill_variant_delins_from_dup)
                            genomic_gap_fill_variant_alt_delins_from_dup = genomic_gap_fill_variant_alt.ac + ':' + genomic_gap_fill_variant_alt.type + '.' + str(
                                genomic_gap_fill_variant_alt.posedit.pos.start.base) + '_' + str(
                                genomic_gap_fill_variant_alt.posedit.pos.end.base) + 'del' + genomic_gap_fill_variant_alt.posedit.edit.ref + 'ins' + genomic_gap_fill_variant_alt.posedit.edit.ref + genomic_gap_fill_variant_alt.posedit.edit.ref
                            genomic_gap_fill_variant_alt = validator.hp.parse_hgvs_variant(
                                genomic_gap_fill_variant_alt_delins_from_dup)

                    # Correct insertion alts
                    if genomic_gap_fill_variant_alt.posedit.edit.type == 'ins':
                        append_ref = validator.sf.fetch_seq(genomic_gap_fill_variant_alt.ac,
                                                       genomic_gap_fill_variant_alt.posedit.pos.start.base - 1,
                                                       genomic_gap_fill_variant_alt.posedit.pos.end.base)
                        genomic_gap_fill_variant_alt.posedit.edit.alt = append_ref[
                                                                            0] + genomic_gap_fill_variant_alt.posedit.edit.alt + \
                                                                        append_ref[1]

                    # Split the reference and replacing alt sequence into a dictionary
                    reference_bases = list(genomic_gap_fill_variant.posedit.edit.ref)
                    if genomic_gap_fill_variant_alt.posedit.edit.alt is not None:
                        alternate_bases = list(genomic_gap_fill_variant_alt.posedit.edit.alt)
                    else:
                        # Deletions with no ins
                        pre_alternate_bases = list(genomic_gap_fill_variant_alt.posedit.edit.ref)
                        alternate_bases = []
                        for base in pre_alternate_bases:
                            alternate_bases.append('X')

                    # Create the dictionaries
                    ref_start = genomic_gap_fill_variant.posedit.pos.start.base
                    alt_start = genomic_gap_fill_variant_alt.posedit.pos.start.base
                    ref_base_dict = {}
                    for base in reference_bases:
                        ref_base_dict[ref_start] = str(base)
                        ref_start = ref_start + 1

                    alt_base_dict = {}

                    # NEED TO SEARCH FOR RANGE = and replace with interval_range
                    # Need to search for int and replace with integer

                    # Note, all variants will be forced into the format delete insert
                    # Deleted bases in the ALT will be substituted for X
                    for integer in range(genomic_gap_fill_variant_alt.posedit.pos.start.base,
                                         genomic_gap_fill_variant_alt.posedit.pos.end.base + 1, 1):
                        if integer == alt_start:
                            alt_base_dict[integer] = str(''.join(alternate_bases))
                        else:
                            alt_base_dict[integer] = 'X'

                    # Generate the alt sequence
                    alternate_sequence_bases = []
                    for integer in range(genomic_gap_fill_variant.posedit.pos.start.base,
                                         genomic_gap_fill_variant.posedit.pos.end.base + 1, 1):
                        if integer in list(alt_base_dict.keys()):
                            alternate_sequence_bases.append(alt_base_dict[integer])
                        else:
                            alternate_sequence_bases.append(ref_base_dict[integer])
                    alternate_sequence = ''.join(alternate_sequence_bases)
                    alternate_sequence = alternate_sequence.replace('X', '')

                    # Add the new alt to the gap fill variant and generate transcript variant
                    genomic_gap_fill_variant.posedit.edit.alt = alternate_sequence
                    hgvs_refreshed_variant = validator.vm.g_to_t(genomic_gap_fill_variant,
                                                            tx_gap_fill_variant.ac)

                    # Set warning
                    gap_size = str(len(genomic_gap_fill_variant.posedit.edit.ref) - 2)
                    disparity_deletion_in[1] = [gap_size]
                    auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                        stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + gap_size + ' genomic base(s) that fail to align to transcript ' + str(
                        tx_hgvs_not_delins.ac)
                    non_valid_caution = 'true'

                    # Alignment position
                    for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                    if re.match('NM_', str(for_location_c)):
                        for_location_c = variant.no_norm_evm.n_to_c(tx_hgvs_not_delins)
                    if re.match(r'\-', str(for_location_c.posedit.pos.start.offset)):
                        gps = for_location_c.posedit.pos.start.base - 1
                        gpe = for_location_c.posedit.pos.start.base
                    else:
                        gps = for_location_c.posedit.pos.start.base
                        gpe = for_location_c.posedit.pos.start.base + 1
                    gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                    auto_info = auto_info + '%s' % (gap_position)
                else:
                    if tx_hgvs_not_delins.posedit.pos.start.offset == 0 and tx_hgvs_not_delins.posedit.pos.end.offset == 0:
                        # In this instance, we have identified a transcript gap but the n. version of
                        # the transcript variant but do not have a position which actually hits the gap,
                        # so the variant likely spans the gap, and is not picked up by an offset.
                        try:
                            c1 = validator.vm.n_to_c(tx_hgvs_not_delins)
                        except:
                            c1 = tx_hgvs_not_delins
                        g1 = validator.nr_vm.t_to_g(c1, variant.hgvs_genomic.ac)
                        g3 = validator.nr_vm.t_to_g(c1, variant.hgvs_genomic.ac)
                        g2 = validator.vm.t_to_g(c1, variant.hgvs_genomic.ac)
                        ng2 = variant.hn.normalize(g2)
                        g3.posedit.pos.end.base = g3.posedit.pos.start.base + (
                                len(g3.posedit.edit.ref) - 1)
                        try:
                            c2 = validator.vm.g_to_t(g3, c1.ac)
                            if c2.posedit.pos.start.offset == 0 and c2.posedit.pos.end.offset == 0:
                                pass
                            else:
                                tx_hgvs_not_delins = c2
                                try:
                                    tx_hgvs_not_delins = validator.vm.c_to_n(tx_hgvs_not_delins)
                                except hgvs.exceptions.HGVSError:
                                    fn.exceptPass()
                        except hgvs.exceptions.HGVSInvalidVariantError:
                            fn.exceptPass()

                    if re.search(r'\+', str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                            r'\+', str(tx_hgvs_not_delins.posedit.pos.end)):
                        auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                            stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                            disparity_deletion_in[
                                1]) + ' genomic base(s) that fail to align to transcript ' + str(
                            tx_hgvs_not_delins.ac)
                        non_valid_caution = 'true'
                        try:
                            c2 = validator.vm.n_to_c(tx_hgvs_not_delins)
                        except:
                            c2 = tx_hgvs_not_delins
                        c1 = copy.deepcopy(c2)
                        c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                        c1.posedit.pos.start.offset = 0
                        c1.posedit.pos.end = c2.posedit.pos.start
                        c1.posedit.edit.ref = ''
                        c1.posedit.edit.alt = ''
                        if orientation != -1:
                            g1 = validator.vm.t_to_g(c1, variant.hgvs_genomic.ac)
                            g2 = validator.vm.t_to_g(c2, variant.hgvs_genomic.ac)
                            g1.posedit.edit.alt = g1.posedit.edit.ref
                        else:
                            g1 = validator.vm.t_to_g(c2, variant.hgvs_genomic.ac)
                            g2 = validator.vm.t_to_g(c1, variant.hgvs_genomic.ac)
                            g2.posedit.edit.alt = g2.posedit.edit.ref
                        reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                        alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                        g3 = copy.deepcopy(g1)
                        g3.posedit.pos.end.base = g2.posedit.pos.end.base
                        g3.posedit.edit.ref = reference
                        g3.posedit.edit.alt = alternate
                        c3 = validator.vm.g_to_t(g3, c1.ac)
                        hgvs_refreshed_variant = c3
                        # Alignment position
                        for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                        if re.match('NM_', str(for_location_c)):
                            for_location_c = variant.no_norm_evm.n_to_c(tx_hgvs_not_delins)
                            gps = for_location_c.posedit.pos.start.base
                            gpe = for_location_c.posedit.pos.start.base + 1
                        gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                        # Warn update
                        auto_info = auto_info + '%s' % (gap_position)
                    elif re.search(r'\+', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search(
                            r'\+', str(tx_hgvs_not_delins.posedit.pos.start)):
                        auto_info = auto_info + 'Genome position ' + str(
                            stored_hgvs_not_delins.ac) + ':g.' + str(
                            stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(
                            disparity_deletion_in[1]) + '-bp gap in transcript ' + str(
                            tx_hgvs_not_delins.ac)
                        gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
                        non_valid_caution = 'true'
                        try:
                            c1 = validator.vm.n_to_c(tx_hgvs_not_delins)
                        except:
                            c1 = tx_hgvs_not_delins
                        c2 = copy.deepcopy(c1)
                        c2.posedit.pos.start = c1.posedit.pos.end
                        c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                        c2.posedit.pos.end.offset = 0
                        c2.posedit.edit.ref = ''
                        c2.posedit.edit.alt = ''
                        if orientation != -1:
                            g1 = validator.vm.t_to_g(c1, variant.hgvs_genomic.ac)
                            g2 = validator.vm.t_to_g(c2, variant.hgvs_genomic.ac)
                            g2.posedit.edit.alt = g2.posedit.edit.ref
                        else:
                            g1 = validator.vm.t_to_g(c2, variant.hgvs_genomic.ac)
                            g2 = validator.vm.t_to_g(c1, variant.hgvs_genomic.ac)
                            g1.posedit.edit.alt = g1.posedit.edit.ref
                        reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                        alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                        g3 = copy.deepcopy(g1)
                        g3.posedit.pos.end.base = g2.posedit.pos.end.base
                        g3.posedit.edit.ref = reference
                        g3.posedit.edit.alt = alternate
                        c3 = validator.vm.g_to_t(g3, c1.ac)
                        hgvs_refreshed_variant = c3
                        # Alignment position
                        for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                        if re.match('NM_', str(for_location_c)):
                            for_location_c = variant.no_norm_evm.n_to_c(tx_hgvs_not_delins)
                        gps = for_location_c.posedit.pos.end.base
                        gpe = for_location_c.posedit.pos.end.base + 1
                        gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                        # Warn update
                        auto_info = auto_info + '%s' % (gap_position)
                    elif re.search(r'\-',
                                   str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                        r'\-', str(tx_hgvs_not_delins.posedit.pos.end)):
                        auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                            stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                            disparity_deletion_in[
                                1]) + ' genomic base(s) that fail to align to transcript ' + str(
                            tx_hgvs_not_delins.ac)
                        non_valid_caution = 'true'
                        try:
                            c2 = validator.vm.n_to_c(tx_hgvs_not_delins)
                        except:
                            c2 = tx_hgvs_not_delins
                        c1 = copy.deepcopy(c2)
                        c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                        c1.posedit.pos.start.offset = 0
                        c1.posedit.pos.end = c2.posedit.pos.start
                        c1.posedit.edit.ref = ''
                        c1.posedit.edit.alt = ''
                        if orientation != -1:
                            g1 = validator.vm.t_to_g(c1, variant.hgvs_genomic.ac)
                            g2 = validator.vm.t_to_g(c2, variant.hgvs_genomic.ac)
                            g1.posedit.edit.alt = g1.posedit.edit.ref
                        else:
                            g1 = validator.vm.t_to_g(c2, variant.hgvs_genomic.ac)
                            g2 = validator.vm.t_to_g(c1, variant.hgvs_genomic.ac)
                            g2.posedit.edit.alt = g2.posedit.edit.ref
                        reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                        alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                        g3 = copy.deepcopy(g1)
                        g3.posedit.pos.end.base = g2.posedit.pos.end.base
                        g3.posedit.edit.ref = reference
                        g3.posedit.edit.alt = alternate
                        c3 = validator.vm.g_to_t(g3, c1.ac)
                        hgvs_refreshed_variant = c3
                        # Alignment position
                        for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                        if re.match('NM_', str(for_location_c)):
                            for_location_c = variant.no_norm_evm.n_to_c(tx_hgvs_not_delins)
                        gps = for_location_c.posedit.pos.start.base - 1
                        gpe = for_location_c.posedit.pos.start.base
                        gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                        # Warn update
                        auto_info = auto_info + '%s' % (gap_position)
                    elif re.search(r'\-', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search(
                            r'\-', str(tx_hgvs_not_delins.posedit.pos.start)):
                        auto_info = auto_info + 'Genome position ' + str(
                            stored_hgvs_not_delins.ac) + ':g.' + str(
                            stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(
                            disparity_deletion_in[1]) + '-bp gap in transcript ' + str(
                            tx_hgvs_not_delins.ac)
                        gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
                        non_valid_caution = 'true'
                        try:
                            c1 = validator.vm.n_to_c(tx_hgvs_not_delins)
                        except:
                            c1 = tx_hgvs_not_delins
                        c2 = copy.deepcopy(c1)
                        c2.posedit.pos.start = c1.posedit.pos.end
                        c2.posedit.pos.end.base = c1.posedit.pos.end.base
                        c2.posedit.pos.end.offset = 0
                        c2.posedit.edit.ref = ''
                        c2.posedit.edit.alt = ''
                        g2 = validator.vm.t_to_g(c2, variant.hgvs_genomic.ac)
                        c2 = validator.vm.g_to_t(g2, c2.ac)
                        reference = c1.posedit.edit.ref + c2.posedit.edit.ref[1:]
                        alternate = c1.posedit.edit.alt + c2.posedit.edit.ref[1:]
                        c3 = copy.deepcopy(c1)
                        c3.posedit.pos.end = c2.posedit.pos.end
                        c3.posedit.edit.ref = ''  # reference
                        c3.posedit.edit.alt = alternate
                        hgvs_refreshed_variant = c3
                        # Alignment position
                        for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                        if re.match('NM_', str(for_location_c)):
                            for_location_c = variant.no_norm_evm.n_to_c(tx_hgvs_not_delins)
                        gps = for_location_c.posedit.pos.end.base - 1
                        gpe = for_location_c.posedit.pos.end.base
                        gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                        # Warn update
                        auto_info = auto_info + '%s' % (gap_position)
                    else:
                        auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                            stored_hgvs_not_delins.posedit.pos) + ' contains ' + str(
                            disparity_deletion_in[
                                1]) + ' genomic base(s) that fail to align to transcript ' + str(
                            tx_hgvs_not_delins.ac) + '\n'
                        hgvs_refreshed_variant = tx_hgvs_not_delins
                        gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)

            # GAP IN THE CHROMOSOME
            elif disparity_deletion_in[0] == 'chromosome':
                # Set warning variables
                gap_position = ''
                gapped_alignment_warning = str(
                    hgvs_genomic_5pr) + ' does not represent a true variant because it is an artefact of aligning the transcripts listed below with genome build ' + variant.primary_assembly
                hgvs_refreshed_variant = tx_hgvs_not_delins
                # Warn
                auto_info = auto_info + str(hgvs_refreshed_variant.ac) + ':c.' + str(
                    hgvs_refreshed_variant.posedit.pos) + ' contains ' + str(disparity_deletion_in[
                                                                                 1]) + ' transcript base(s) that fail to align to chromosome ' + str(
                    variant.hgvs_genomic.ac) + '\n'
                gapped_transcripts = gapped_transcripts + str(hgvs_refreshed_variant.ac) + ' '
            else:
                # Try the push
                hgvs_stash = copy.deepcopy(stash_hgvs_not_delins)
                stash_ac = hgvs_stash.ac
                # Make a hard left and hard right not delins g.
                stash_dict_right = vvHGVS.hard_right_hgvs2vcf(hgvs_stash, variant.primary_assembly, variant.hn, validator.sf)
                stash_pos_right = int(stash_dict_right['pos'])
                stash_ref_right = stash_dict_right['ref']
                stash_alt_right = stash_dict_right['alt']
                stash_end_right = str(stash_pos_right + len(stash_ref_right) - 1)
                stash_hgvs_not_delins_right = validator.hp.parse_hgvs_variant(
                    stash_ac + ':' + hgvs_stash.type + '.' + str(
                        stash_pos_right) + '_' + stash_end_right + 'del' + stash_ref_right + 'ins' + stash_alt_right)
                stash_dict_left = vvHGVS.hard_left_hgvs2vcf(hgvs_stash, variant.primary_assembly,
                                                            variant.reverse_normalizer, validator.sf)
                stash_pos_left = int(stash_dict_left['pos'])
                stash_ref_left = stash_dict_left['ref']
                stash_alt_left = stash_dict_left['alt']
                stash_end_left = str(stash_pos_left + len(stash_ref_left) - 1)
                stash_hgvs_not_delins_left = validator.hp.parse_hgvs_variant(
                    stash_ac + ':' + hgvs_stash.type + '.' + str(
                        stash_pos_left) + '_' + stash_end_left + 'del' + stash_ref_left + 'ins' + stash_alt_left)
                # Map in-situ to the transcript left and right
                try:
                    tx_hard_right = validator.vm.g_to_t(stash_hgvs_not_delins_right, saved_hgvs_coding.ac)
                except Exception as e:
                    tx_hard_right = saved_hgvs_coding
                else:
                    normalize_stash_right = variant.hn.normalize(stash_hgvs_not_delins_right)
                    if str(normalize_stash_right.posedit) == str(stash_hgvs_not_delins.posedit):
                        tx_hard_right = saved_hgvs_coding
                try:
                    tx_hard_left = validator.vm.g_to_t(stash_hgvs_not_delins_left, saved_hgvs_coding.ac)
                except Exception as e:
                    tx_hard_left = saved_hgvs_coding
                else:
                    normalize_stash_left = variant.hn.normalize(stash_hgvs_not_delins_left)
                    if str(normalize_stash_left.posedit) == str(stash_hgvs_not_delins.posedit):
                        tx_hard_left = saved_hgvs_coding
                # The Logic - Currently limited to genome gaps
                if len(stash_hgvs_not_delins_right.posedit.edit.ref) < len(
                        tx_hard_right.posedit.edit.ref):
                    tx_hard_right = variant.hn.normalize(tx_hard_right)
                    gap_position = ''
                    gapped_alignment_warning = str(
                        hgvs_genomic_5pr) + ' may be an artefact of aligning the transcripts listed below with genome build ' + variant.primary_assembly
                    hgvs_refreshed_variant = tx_hard_right
                    gapped_transcripts = gapped_transcripts + str(tx_hard_right.ac) + ' '
                elif len(stash_hgvs_not_delins_left.posedit.edit.ref) < len(
                        tx_hard_left.posedit.edit.ref):
                    tx_hard_left = variant.hn.normalize(tx_hard_left)
                    gap_position = ''
                    gapped_alignment_warning = str(
                        hgvs_genomic_5pr) + ' may be an artefact of aligning the transcripts listed below with genome build ' + variant.primary_assembly
                    hgvs_refreshed_variant = tx_hard_left
                    gapped_transcripts = gapped_transcripts + str(tx_hard_left.ac) + ' '
                else:
                    # Keep the same by re-setting rel_var
                    hgvs_refreshed_variant = saved_hgvs_coding

            # Edit the output
            if re.match('NM_', str(hgvs_refreshed_variant.ac)) and not re.search('c', str(
                    hgvs_refreshed_variant.type)):
                hgvs_refreshed_variant = variant.evm.n_to_c(hgvs_refreshed_variant)
            else:
                pass
            try:
                hgvs_refreshed_variant = variant.hn.normalize(hgvs_refreshed_variant)
                if hgvs_refreshed_variant.posedit.edit.type == 'delins' and \
                        hgvs_refreshed_variant.posedit.edit.ref[-1] == \
                        hgvs_refreshed_variant.posedit.edit.alt[-1]:
                    hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[
                                                              0:-1]
                    hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[
                                                              0:-1]
                    hgvs_refreshed_variant.posedit.pos.end.base = hgvs_refreshed_variant.posedit.pos.end.base - 1
                    hgvs_refreshed_variant = variant.hn.normalize(hgvs_refreshed_variant)
                elif hgvs_refreshed_variant.posedit.edit.type == 'delins' and \
                        hgvs_refreshed_variant.posedit.edit.ref[0] == \
                        hgvs_refreshed_variant.posedit.edit.alt[0]:
                    hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[
                                                              1:]
                    hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[
                                                              1:]
                    hgvs_refreshed_variant.posedit.pos.start.base = hgvs_refreshed_variant.posedit.pos.start.base + 1
                    hgvs_refreshed_variant = variant.hn.normalize(hgvs_refreshed_variant)
            except Exception as e:
                error = str(e)
                # Ensure the final variant is not intronic nor does it cross exon boundaries
                if re.match('Normalization of intronic variants is not supported',
                            error) or re.match(
                    'Unsupported normalization of variants spanning the exon-intron boundary',
                    error):
                    hgvs_refreshed_variant = saved_hgvs_coding
                else:
                    pass
                fn.exceptPass()
            # Send to empty nw_rel_var
            nw_rel_var.append(hgvs_refreshed_variant)

        # Otherwise these variants need to be set
        else:
            corrective_action_taken = ''
            gapped_alignment_warning = ''
            # Send to empty nw_rel_var
            nw_rel_var.append(saved_hgvs_coding)

    data = {
        'gapped_alignment_warning': gapped_alignment_warning,
        'corrective_action_taken': corrective_action_taken,
        'auto_info': auto_info,
        'disparity_deletion_in': disparity_deletion_in,
        'gapped_transcripts': gapped_transcripts
    }
    return data, nw_rel_var
