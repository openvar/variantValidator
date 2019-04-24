import hgvs
import re
import copy
import time
import sys
from .vvLogging import logger
from . import vvHGVS
from .variant import Variant
from . import vvChromosomes
from . import vvFunctions as fn
from . import gapped_mapping


def gene_to_transcripts(variant, validator):
    g_query = validator.hp.parse_hgvs_variant(str(variant.hgvs_formatted))

    # Genomic coordinates can be validated immediately
    error = 'false'
    try:
        validator.vr.validate(g_query)
    except hgvs.exceptions.HGVSError as e:
        error = str(e)
    except KeyError:
        error = 'Reference sequence ' + variant.hgvs_genomic.ac + ' is either not supported or does not exist'
    if error != 'false':
        reason = 'Invalid variant description'
        variant.warnings += ': ' + str(error)
        logger.warning(str(error))
        return True

    # Set test to see if Norm alters the coords
    g_test = variant.hn.normalize(g_query)

    # Perform test
    if g_query.posedit.pos != g_test.posedit.pos:
        # my_variant.warnings += ': ' + 'Input variant description normalized to ' + str(g_test)
        variant.hgvs_genomic = g_test
    else:
        variant.hgvs_genomic = g_query

    # Collect rel_var
    # rel_var is a keyworded list of relevant transcripts with associated coding variants
    """
    Initial simple projection from the provided g. position all overlapping
    transcripts
    """
    rel_var = validator.relevant_transcripts(variant.hgvs_genomic, variant.evm, validator.alt_aln_method, variant.reverse_normalizer)

    # Double check rel_vars have not been missed when mapping from a RefSeqGene
    if len(rel_var) != 0 and 'NG_' in variant.hgvs_genomic.ac:
        for var in rel_var:
            hgvs_coding_variant = validator.hp.parse_hgvs_variant(var)
            try:
                variant.hgvs_genomic = validator.myevm_t_to_g(hgvs_coding_variant, variant.no_norm_evm,
                                                 variant.primary_assembly, variant.hn)
            except hgvs.exceptions.HGVSError as e:
                try_rel_var = []
            else:
                try_rel_var = validator.relevant_transcripts(variant.hgvs_genomic, variant.evm, validator.alt_aln_method,
                                                        variant.reverse_normalizer)
            if len(try_rel_var) > len(rel_var):
                rel_var = try_rel_var
                break
            else:
                continue

    #  Tripple check this assumption by querying the gene position database
    if len(rel_var) == 0:
        vcf_dict = vvHGVS.hgvs2vcf(variant.hgvs_genomic, variant.primary_assembly, variant.reverse_normalizer, validator.sf)
        not_di = str(variant.hgvs_genomic.ac) + ':g.' + str(vcf_dict['pos']) + '_' + str(
            int(vcf_dict['pos']) + (len(vcf_dict['ref']) - 1)) + 'del' + vcf_dict['ref'] + 'ins' + \
                 vcf_dict['alt']
        hgvs_not_di = validator.hp.parse_hgvs_variant(not_di)
        rel_var = validator.relevant_transcripts(hgvs_not_di, variant.evm, validator.alt_aln_method,
                                            variant.reverse_normalizer)

    # list return statements
    """
    If mapping to transcripts has been unsuccessful, provide relevant details
    """
    if len(rel_var) == 0:

        # Check for NG_
        if str(variant.hgvs_formatted).startswith('NG_'):
            # parse
            hgvs_refseqgene = validator.hp.parse_hgvs_variant(str(variant.hgvs_formatted))
            # Convert to chromosomal position
            refseqgene_data = validator.rsg_to_chr(hgvs_refseqgene, variant.primary_assembly, variant.hn, validator.vr)
            # There should only ever be one description returned
            refseqgene_data = refseqgene_data[0]

            # Extract data
            if refseqgene_data['valid'] == 'true':
                input = refseqgene_data['hgvs_genomic']
                # re_submit
                # Tag the line so that it is not written out
                variant.warnings += ': ' + str(variant.hgvs_formatted) + ' automapped to genome position ' + str(input)
                query = Variant(variant.original, quibble=input, warnings=variant.warnings,
                                        primary_assembly=variant.primary_assembly, order=variant.order)

                coding = 'intergenic'
                validator.batch_list.append(query)
            else:
                error = 'Mapping unavailable for RefSeqGene ' + str(variant.hgvs_formatted) + ' using alignment method = ' + validator.alt_aln_method
                variant.warnings += ': ' + str(error)
                logger.warning(str(error))
                return True

        # Chromosome build is not supported or intergenic???
        else:
            sfm = vvChromosomes.supported_for_mapping(variant.hgvs_genomic.ac, variant.primary_assembly)
            if sfm == 'true':
                try:
                    validator.vr.validate(variant.hgvs_genomic)
                except hgvs.exceptions.HGVSError as e:
                    error = str(e)
                    variant.warnings += ': ' + str(error)
                    logger.warning(str(error))
                    return True
                else:
                    # Map to RefSeqGene if available
                    refseqgene_data = validator.chr_to_rsg(variant.hgvs_genomic, variant.hn, validator.vr)
                    rsg_data = ''
                    # Example {'gene': 'NTHL1', 'hgvs_refseqgene': 'NG_008412.1:g.3455_3464delCAAACACACA', 'valid': 'true'}
                    for data in refseqgene_data:
                        if data['valid'] == 'true':
                            data['hgvs_refseqgene'] = validator.hp.parse_hgvs_variant(data['hgvs_refseqgene'])
                            data['hgvs_refseqgene'] = fn.valstr(data['hgvs_refseqgene'])
                            rsg_data = rsg_data + data['hgvs_refseqgene'] + ' (' + data['gene'] + '), '

                    error = 'No transcripts found that fully overlap the described variation in the genomic sequence'
                    # set output type flag
                    variant.output_type_flag = 'intergenic'
                    # set genomic and where available RefSeqGene outputs
                    variant.warnings += ': ' + str(error)
                    variant.genomic_g = fn.valstr(variant.hgvs_genomic)
                    variant.genomic_r = str(rsg_data.split('(')[0])
                    logger.warning(str(error))
                    return True
            else:
                error = 'Please ensure the requested chromosome version relates to a supported genome build. Supported genome builds are: GRCh37, GRCh38, hg19 and hg38'
                variant.warnings += ': ' + str(error)
                logger.warning(str(error))
                return True

    else:
        # Tag the line so that it is not written out
        variant.write = False

        data, nw_rel_var = gapped_mapping.gapped_g_to_c(variant, validator, rel_var)

        # Warn the user that the g. description is not valid
        if data['gapped_alignment_warning'] != '':
            if data['disparity_deletion_in'][0] == 'transcript':
                corrective_action_taken = 'Automap has deleted  ' + str(
                    data['disparity_deletion_in'][1]) + ' bp from chromosomal reference sequence ' + str(
                    variant.hgvs_genomic.ac) + ' to ensure perfect alignment with transcript reference sequence(s)' + data['gapped_transcripts']
            if data['disparity_deletion_in'][0] == 'chromosome':
                corrective_action_taken = 'Automap has added  ' + str(
                    data['disparity_deletion_in'][1]) + ' bp to chromosomal reference sequence ' + str(
                    variant.hgvs_genomic.ac) + ' to ensure perfect alignment with transcript reference sequence(s) ' + data['gapped_transcripts']

        # Add additional data to the front of automap
        if data['auto_info'] != '':
            automap = data['auto_info'] + '\n' + 'false'

        rel_var = nw_rel_var

        # Set the values and append to batch_list
        for c_description in rel_var:
            query = Variant(variant.original, quibble=str(c_description), warnings=variant.warnings,
                                    primary_assembly=variant.primary_assembly, order=variant.order)
            validator.batch_list.append(query)
        logger.warning("Continue reached when mapping transcript types to variants")
        # Call next description
        return True
    return False


def transcripts_to_gene(variant, validator):
    """This seems to use the quibble and not the HGVS formatted variant format."""

    # Flag for validation
    valid = False
    boundary = 'false'
    warning = ''
    caution = ''
    error = ''
    # Collect information for genomic level validation
    obj = validator.hp.parse_hgvs_variant(str(variant.hgvs_formatted))

    tx_ac = obj.ac

    input = str(variant.quibble)
    formatted_variant = str(variant.hgvs_formatted)

    # Do we keep it?
    if validator.select_transcripts != 'all':
        if tx_ac in list(validator.select_transcripts_dict_plus_version.keys()):
            pass
        # If not get rid of it!
        else:
            # By marking it as Do Not Write and continuing through the validation loop
            variant.write = False
            return True
    else:
        pass

    print(variant.hgvs_formatted)
    print(variant.quibble)
    # Set a cross_variant object
    cross_variant = 'false'
    # Se rec_var to '' so it can be updated later
    rec_var = ''

    # First task is to get the genomic equivalent, and print useful error messages if it can't be found.
    try:
        to_g = validator.myevm_t_to_g(obj, variant.no_norm_evm, variant.primary_assembly, variant.hn)
        print('Genomic:', to_g)
        genomic_ac = to_g.ac
    except hgvs.exceptions.HGVSDataNotAvailableError as e:
        if ('~' in str(e) and 'Alignment is incomplete' in str(e)) or "No relevant genomic mapping options" in str(e):
            # Unable to map the input variant onto a genomic position
            if '~' in str(e) and 'Alignment is incomplete' in str(e):
                error_list = str(e).split('~')[:-1]
                combos = [
                    'Full alignment data between the specified transcript reference sequence and all GRCh37 and GRCh38 '
                    'genomic reference sequences (including alternate chromosome assemblies, patches and RefSeqGenes) '
                    'are not available: Consequently the input variant description cannot be fully validated and is '
                    'not supported: Use the Gene to Transcripts function to determine whether an updated transcript '
                    'reference sequence is available']
                # Partial alignment data is available for the following genomic reference sequences: ']
                error = '; '.join(combos)
                error = error.replace(': ;', ': ')
            else:
                error = str(e)
                error = error + ': Consequently the input variant description cannot be fully validated and is not ' \
                                'supported: Use the Gene to Transcripts function to determine whether an updated ' \
                                'transcript reference sequence is available'
            variant.warnings += ': ' + error
            logger.warning(error)
            return True
        try:
            gene_symbol = validator.db.get_gene_symbol_from_transcriptID(tx_ac)
        except:
            gene_symbol = None
        if gene_symbol is None:
            error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, ' \
                    'please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' to  ' \
                    'https://variantvalidator.org/ref_finder/, or select an alternative genome build'
        else:
            error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, ' \
                    'please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' or ' + \
                    gene_symbol + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'

        variant.warnings += ': ' + error
        logger.warning(error)
        return True
    except TypeError:
        try:
            gene_symbol = validator.db.get_gene_symbol_from_transcriptID(tx_ac)
        except:
            gene_symbol = 'none'
        if gene_symbol == 'none':
            error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, ' \
                    'please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' to  ' \
                    'https://variantvalidator.org/ref_finder/, or select an alternative genome build'
        else:
            error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, ' \
                    'please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' or ' + \
                    gene_symbol + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
        variant.warnings += ': ' + error
        logger.warning(error)
        return True

    # Get orientation of the gene wrt genome and a list of exons mapped to the genome
    ori = validator.tx_exons(tx_ac=tx_ac, alt_ac=genomic_ac, alt_aln_method=validator.alt_aln_method)
    #print('exons:', ori)
    orientation = int(ori[0]['alt_strand'])
    intronic_variant = 'false'

    # Collect variant sequence information via normalisation (normalizer) or if intronic via mapping
    # INTRONIC OFFSETS - Required for Exon table
    # Variable to collect offset to exon boundary
    ex_offset = 0
    plus = re.compile(r"\d\+\d")  # finds digit + digit
    minus = re.compile(r"\d\-\d")  # finds digit - digit

    geno = re.compile(r':g.')
    if plus.search(input) or minus.search(input):
        if 'error' in str(to_g):
            if validator.alt_aln_method != 'genebuild':
                error = "If the following error message does not address the issue and the problem persists please " \
                        "contact admin: " + str(to_g)
                variant.warnings += ': ' + error
                logger.warning(error)
                return True

            else:
                error = "If the following error message does not address the issue and the problem persists please " \
                        "contact admin: " + str(to_g)
                variant.warnings += ': ' + error
                logger.warning(error)
                return True

        else:
            # Insertions at exon boundaries are miss-handled by vm.g_to_t
            if (obj.posedit.edit.type == 'ins' and
                obj.posedit.pos.start.offset == 0 and
                obj.posedit.pos.end.offset != 0) or (obj.posedit.edit.type == 'ins' and
                                                     obj.posedit.pos.start.offset != 0 and
                                                     obj.posedit.pos.end.offset == 0):
                formatted_variant = str(obj)
            else:
                # Normalize was I believe to replace ref. Mapping does this anyway
                # to_g = variant.hn.normalize(to_g)
                formatted_variant = str(validator.myevm_g_to_t(variant.evm, to_g, tx_ac))
                tx_ac = ''

    elif ':g.' in input:
        if plus.search(formatted_variant) or minus.search(formatted_variant):
            to_g = validator.genomic(formatted_variant, variant.no_norm_evm, variant.primary_assembly, variant.hn)
            if 'error' in str(to_g):
                if validator.alt_aln_method != 'genebuild':
                    error = "If the following error message does not address the issue and the problem persists " \
                            "please contact admin: " + str(to_g)
                    variant.warnings += ': ' + error
                    logger.warning(error)
                    return True

                else:
                    error = "If the following error message does not address the issue and the problem persists " \
                            "please contact admin: " + str(to_g)
                    variant.warnings += ': ' + error
                    logger.warning(error)
                    return True
        else:
            # Insertions at exon boundaries are miss-handled by vm.g_to_t
            if (obj.posedit.edit.type == 'ins' and
                obj.posedit.pos.start.offset == 0 and
                obj.posedit.pos.end.offset != 0) or (obj.posedit.edit.type == 'ins' and
                                                     obj.posedit.pos.start.offset != 0 and
                                                     obj.posedit.pos.end.offset == 0):
                formatted_variant = str(obj)
            else:
                # Normalize was I believe to replace ref. Mapping does this anyway
                # to_g = hn.normalize(to_g)
                formatted_variant = str(validator.myevm_g_to_t(variant.evm, to_g, tx_ac))
                tx_ac = ''

    else:
        # Normalize the variant
        error = 'false'
        try:
            h_variant = variant.hn.normalize(obj)
        except hgvs.exceptions.HGVSUnsupportedOperationError as e:
            error = str(e)
            if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                h_variant = obj
                formatted_variant = formatted_variant
                caution = 'This coding sequence variant description spans at least one intron'
                automap = 'Use of the corresponding genomic sequence variant descriptions may be invalid. ' \
                          'Please refer to https://www35.lamp.le.ac.uk/recommendations/'
                variant.warnings += ': ' + caution + ': ' + automap
                logger.warning(caution + ": " + automap)
        else:
            formatted_variant = str(h_variant)

        # tx_ac = ''
        # # Create a crosser (exon boundary crossed) variant
        # crossed_variant = str(variant.evm._maybe_normalize(obj))
        # if formatted_variant == crossed_variant:
        #     cross_variant = 'false'
        # else:
        #     hgvs_crossed_variant = variant.evm._maybe_normalize(obj)
        #     cross_variant = [
        #         "Coding sequence allowing for exon boundary crossing (default = no crossing)",
        #         crossed_variant, hgvs_crossed_variant.ac]
        #     cr_available = 'true'
        #
        # # control of cross_variant
        # if boundary == 'false':
        #     cross_variant = 'false'

        # Moved this forwards and removed the previous section as it doesn't seem to be used anywhere

        error = validator.validateHGVS(formatted_variant)
        if error == 'false':
            valid = True
        else:
            variant.warnings += ': ' + str(error)
            logger.warning(str(error))
            return True

    # Tackle the plus intronic offset
    cck = False
    if plus.search(input):
        # Regular expression catches the start of the interval only based on .00+00 pattern
        inv_start = re.compile(r"\.\d+\+\d")
        if inv_start.search(input):
            cck = True
    if minus.search(input):
        # Regular expression catches the start of the interval only based on .00-00 pattern
        inv_start = re.compile(r"\.\d+\-\d")
        if inv_start.search(input):
            cck = True

    # COORDINATE CHECKER
    # hgvs will handle incorrect coordinates so need to automap errors
    # Make sure any input intronic coordinates are correct
    # Get the desired transcript
    pat_r = re.compile(':r.')
    pat_g = re.compile(':g.')
    if cck:
        # This should only ever hit coding and RNA variants
        if 'del' in formatted_variant:
            # RNA - looking at trapped variant which was saved before RNA converted to cDNA
            #TODO: rename variant.trapped to variant.pre_RNA_conversion or something similar so it makes sense.
            if ':r.' in variant.trapped:
                coding = validator.coding(formatted_variant, validator.hp)
                trans_acc = coding.ac
                # c to Genome coordinates - Map the variant to the genome
                pre_var = validator.genomic(formatted_variant, variant.no_norm_evm, variant.primary_assembly,variant.hn)
                # genome back to C coordinates
                post_var = validator.myevm_g_to_t(variant.evm, pre_var, trans_acc)

                test = validator.hp.parse_hgvs_variant(input)
                if post_var.posedit.pos.start.base != test.posedit.pos.start.base or \
                        post_var.posedit.pos.end.base != test.posedit.pos.end.base:
                    caution = 'The entered coordinates do not agree with the intron/exon boundaries for the selected ' \
                              'transcript:'
                    automap = 'Automap has corrected the coordinates to match the intron/exon boundaries for the ' \
                              'selected transcript'
                    # automapping of variant completed
                    # Change to rna variant
                    # TODO: Need to look this section over. Doesn't make any sense.
                    # THERE IS NO SUCH THING AS QUERY. THIS WOULDN'T HAVE WORKED AND ISN'T RUN IN ANY TESTS
                    query = variant  # Deliberately won't work so I can fix this once I have an appropriate test.
                    posedit = query.posedit
                    posedit = posedit.lower()
                    query.posedit = posedit
                    query.type = 'r'
                    post_var = str(query)
                    automap = variant.trapped + ' automapped to ' + str(post_var)
                    variant.warnings += ': ' + str(caution) + ': ' + str(automap)

                    # Kill current line and append for re-submission
                    # Tag the line so that it is not written out
                    variant.write = False
                    # Set the values and append to batch_list
                    hgvs_vt = validator.hp.parse_hgvs_variant(str(post_var))
                    assert str(hgvs_vt) == str(post_var)
                    query = Variant(variant.original, quibble=fn.valstr(hgvs_vt), warnings=automap,
                                    primary_assembly=variant.primary_assembly, order=variant.order)
                    validator.batch_list.append(query)

            # Coding
            else:
                coding = validator.coding(formatted_variant, validator.hp)
                trans_acc = coding.ac
                # c to Genome coordinates - Map the variant to the genome
                pre_var = validator.hp.parse_hgvs_variant(formatted_variant)
                try:
                    pre_var = validator.myevm_t_to_g(pre_var, variant.no_norm_evm, variant.primary_assembly,
                                                variant.hn)
                except Exception as e:
                    error = str(e)
                    if error == 'expected from_start_i <= from_end_i':
                        error = 'Automap is unable to correct the input exon/intron boundary coordinates, ' \
                                'please check your variant description'
                        variant.warnings += ': ' + error
                        return True
                    else:
                        fn.exceptPass()
                # genome back to C coordinates
                try:
                    post_var = validator.myevm_g_to_t(variant.evm, pre_var, trans_acc)
                except hgvs.exceptions.HGVSError as error:
                    variant.warnings += ': ' + str(error)
                    logger.warning(str(error))
                    return True
                query = post_var
                test = validator.hp.parse_hgvs_variant(input)

                if post_var.posedit.pos.start.base != test.posedit.pos.start.base or \
                        post_var.posedit.pos.end.base != test.posedit.pos.end.base:
                    caution = 'The entered coordinates do not agree with the intron/exon boundaries for the ' \
                              'selected transcript:'
                    automap = 'Automap has corrected the coordinates to match the intron/exon boundaries for the ' \
                              'selected transcript'
                    # automapping of variant completed
                    automap = variant.trapped + ' automapped to ' + str(post_var)
                    variant.warnings += str(caution) + ': ' + str(automap)
                    relevant = "Select the automapped transcript and click Submit to analyse"

                    # Kill current line and append for re-submission
                    # Tag the line so that it is not written out
                    variant.write = False
                    # Set the values and append to batch_list
                    hgvs_vt = validator.hp.parse_hgvs_variant(str(post_var))
                    assert str(hgvs_vt) == str(post_var)
                    query = Variant(variant.original, quibble=fn.valstr(hgvs_vt), warnings=automap, primary_assembly=variant.primary_assembly, order=variant.order)
                    validator.batch_list.append(query)

        else:  # del not in formatted_variant
            if pat_r.search(variant.trapped):
                coding = validator.coding(formatted_variant, validator.hp)
                trans_acc = coding.ac
                # c to Genome coordinates - Map the variant to the genome
                pre_var = validator.genomic(formatted_variant, variant.no_norm_evm, variant.primary_assembly,variant.hn)
                # genome back to C coordinates
                post_var = validator.myevm_g_to_t(variant.evm, pre_var, trans_acc)

                test = validator.hp.parse_hgvs_variant(input)
                if post_var.posedit.pos.start.base != test.posedit.pos.start.base or post_var.posedit.pos.end.base != test.posedit.pos.end.base:
                    caution = 'The entered coordinates do not agree with the intron/exon boundaries for the selected transcript:'
                    automap = 'Automap has corrected the coordinates to match the intron/exon boundaries for the selected transcript'
                    # automapping of variant completed
                    # Change to rna variant
                    # TODO: As before this section needs fixing
                    # THERE IS NO SUCH THING AS QUERY. THIS WOULDN'T HAVE WORKED AND ISN'T RUN IN ANY TESTS
                    query = variant
                    posedit = query.posedit
                    posedit = posedit.lower()
                    query.posedit = posedit
                    query.type = 'r'
                    post_var = str(query)
                    automap = input + ' automapped to ' + post_var
                    variant.warnings += ': ' + str(caution) + ': ' + str(
                        automap)
                    relevant = "Select the automapped transcript and click Submit to analyse"

                    # Kill current line and append for re-submission
                    # Tag the line so that it is not written out
                    variant.write = False
                    # Set the values and append to batch_list
                    hgvs_vt = validator.hp.parse_hgvs_variant(str(post_var))
                    assert str(hgvs_vt) == str(post_var)
                    query = Variant(variant.original, quibble=fn.valstr(hgvs_vt), warnings=automap, primary_assembly=variant.primary_assembly, order=variant.order)
                    validator.batch_list.append(query)

            else:
                coding = validator.coding(formatted_variant, validator.hp)
                trans_acc = coding.ac
                # c to Genome coordinates - Map the variant to the genome
                pre_var = validator.genomic(formatted_variant, variant.no_norm_evm, variant.primary_assembly,variant.hn)

                # genome back to C coordinates
                post_var = validator.myevm_g_to_t(variant.evm, pre_var, trans_acc)

                test = validator.hp.parse_hgvs_variant(input)
                if post_var.posedit.pos.start.base != test.posedit.pos.start.base or post_var.posedit.pos.end.base != test.posedit.pos.end.base:
                    caution = 'The entered coordinates do not agree with the intron/exon boundaries for the selected transcript:'
                    automap = 'Automap has corrected the coordinates to match the intron/exon boundaries for the selected transcript'
                    # automapping of variant completed
                    automap = str(variant.trapped) + ' automapped to ' + str(post_var)
                    variant.warnings += ': ' + str(caution) + ': ' + str(
                        automap)
                    relevant = "Select the automapped transcript and click Submit to analyse"

                    # Kill current line and append for re-submission
                    # Tag the line so that it is not written out
                    variant.write = False
                    # Set the values and append to batch_list
                    hgvs_vt = validator.hp.parse_hgvs_variant(str(post_var))
                    assert str(hgvs_vt) == str(post_var)
                    query = Variant(variant.original, quibble=fn.valstr(hgvs_vt), warnings=automap, primary_assembly=variant.primary_assembly, order=variant.order)
                    validator.batch_list.append(query)

    # If cck not true
    elif pat_r.search(variant.trapped):
        # set input hgvs object
        hgvs_rna_input = validator.hp.parse_hgvs_variant(variant.trapped)  # Traps the hgvs variant of r. for further use
        inp = str(validator.hgvs_r_to_c(hgvs_rna_input))
        # Regex
        if plus.search(input) or minus.search(input):
            to_g = validator.genomic(inp, variant.no_norm_evm, variant.primary_assembly, variant.hn)
            if 'error' in str(to_g):
                error = "If the following error message does not address the issue and the problem persists " \
                        "please contact admin: " + to_g
                variant.warnings += ': ' + error
                logger.warning(error)
                return True

            else:
                # Set variants pre and post genomic norm
                hgvs_inp = validator.myevm_g_to_t(variant.evm, to_g, tx_ac=obj.ac)
                to_g = variant.hn.normalize(to_g)
                hgvs_otp = validator.myevm_g_to_t(variant.evm, to_g, tx_ac=obj.ac)
        else:
            # Set variants pre and post RNA norm
            hgvs_inp = validator.hp.parse_hgvs_variant(inp)
            try:
                hgvs_otp = variant.hn.normalize(hgvs_inp)
            except hgvs.exceptions.HGVSError as e:
                hgvs_otp = hgvs_inp

        # Set remaining variables
        hgvs_otp.posedit.edit = str(hgvs_otp.posedit.edit).lower()
        otp = str(hgvs_otp)
        query = str(hgvs_otp.posedit.pos)
        test = str(hgvs_inp.posedit.pos)
        query = query.replace('T', 'U')
        query = query.replace('ENSU', 'ENST')
        test = test.replace('T', 'U')
        test = test.replace('ENSU', 'ENST')
        output = otp.replace(':c.', ':r.')
        # Apply coordinates test
        if query != test:
            caution = 'The variant description ' + input + ' requires alteration to comply with HGVS variant ' \
                                                           'nomenclature:'
            # automapping of variant completed
            automap = variant.trapped + ' automapped to ' + output
            variant.warnings += ': ' + caution + ': ' + automap

            # Kill current line and append for re-submission
            # Tag the line so that it is not written out
            variant.write = False
            # Set the values and append to batch_list
            hgvs_vt = validator.hp.parse_hgvs_variant(str(query))
            assert str(hgvs_vt) == str(query)
            query = Variant(variant.original, quibble=fn.valstr(hgvs_vt), warnings=automap, primary_assembly=variant.primary_assembly, order=variant.order)
            validator.batch_list.append(query)

    elif pat_g.search(input):
        pass

    else:
        query = validator.hp.parse_hgvs_variant(formatted_variant)
        test = validator.hp.parse_hgvs_variant(input)
        if query.posedit.pos != test.posedit.pos:
            caution = 'The variant description ' + input + ' requires alteration to comply with HGVS variant ' \
                                                           'nomenclature:'
            # automapping of variant completed
            automap = str(test) + ' automapped to ' + str(query)
            variant.warnings += ': ' + caution + ': ' + automap

            # Kill current line and append for re-submission
            # Tag the line so that it is not written out
            variant.write = False
            # Set the values and append to batch_list
            hgvs_vt = validator.hp.parse_hgvs_variant(str(query))
            assert str(hgvs_vt) == str(query)
            query = Variant(variant.original, quibble=fn.valstr(hgvs_vt), warnings=automap, primary_assembly=variant.primary_assembly, order=variant.order)
            validator.batch_list.append(query)

    # VALIDATION of intronic variants
    pre_valid = validator.hp.parse_hgvs_variant(input)
    post_valid = validator.hp.parse_hgvs_variant(formatted_variant)

    # valid is false if the input contains a \d+\d, \d-\d or :g.
    if not valid:
        genomic_validation = str(validator.genomic(input, variant.no_norm_evm, variant.primary_assembly, variant.hn))
        if fn.valstr(pre_valid) != fn.valstr(post_valid):
            if variant.reftype != ':g.':
                if caution == '':
                    caution = fn.valstr(pre_valid) + ' automapped to ' + fn.valstr(post_valid)
                variant.warnings += ': ' + caution
                logger.warning(caution)

        # Apply validation to intronic variant descriptions (should be valid but make sure)
        error = validator.validateHGVS(genomic_validation)
        if error == 'false':
            valid = True
        else:
            variant.warnings += ': ' + error
            return True

    assert valid is True
    # If valid is False we won't reach this part, so I can remove the if condition

    var_tab = 'true'
    cores = "HGVS-compliant variant descriptions" + warning

    # v0.1a1 edit
    if fn.valstr(pre_valid) != fn.valstr(post_valid):
        if variant.reftype == ':g.':
            if caution == '':
                caution = fn.valstr(pre_valid) + ' automapped to ' + fn.valstr(post_valid)
            variant.warnings += ': ' + str(caution)

    # COLLECT VARIANT DESCRIPTIONS
    ##############################

    # Coding sequence - BASED ON NORMALIZED VARIANT IF EXONIC
    hgvs_coding = validator.coding(formatted_variant, validator.hp)
    boundary = re.compile('exon-intron boundary')
    spanning = re.compile('exon/intron')

    try:
        hgvs_coding = variant.hn.normalize(hgvs_coding)
    except hgvs.exceptions.HGVSError as e:
        error = str(e)

    # Gap compensating code status
    gap_compensation = True

    # Gap gene black list
    try:
        gene_symbol = validator.db.get_gene_symbol_from_transcriptID(hgvs_coding.ac)
    except Exception:
        fn.exceptPass()
    else:
        # If the gene symbol is not in the list, the value False will be returned
        gap_compensation = vvChromosomes.gap_black_list(gene_symbol)

    # Intron spanning variants
    if 'boundary' in str(error) or 'spanning' in str(error):
        try:
            hgvs_coding = variant.evm._maybe_normalize(hgvs_coding)
            gap_compensation = False
        except hgvs.exceptions.HGVSError as error:
            variant.warnings += ': ' + str(error)
            logger.warning(str(error))
            return True

    # Warn status
    logger.warning("gap_compensation_1 = " + str(gap_compensation))
    coding = fn.valstr(hgvs_coding)

    # RNA sequence
    hgvs_rna = copy.deepcopy(hgvs_coding)
    hgvs_rna = validator.hgvs_c_to_r(hgvs_rna)
    rna = str(hgvs_rna)

    # Genomic sequence
    hgvs_genomic = validator.myevm_t_to_g(hgvs_coding, variant.no_norm_evm, variant.primary_assembly, variant.hn)
    final_hgvs_genomic = hgvs_genomic

    # genomic_possibilities
    # 1. take the simple 3 pr normalized hgvs_genomic
    # 2. Lock in hgvs_genomic at its most 5 prime position wrt genome
    hgvs_genomic_possibilities = []

    # Loop out gap finding code under these circumstances!
    if gap_compensation is True:

        hgvs_genomic, gapped_transcripts, auto_info, suppress_c_normalization, hgvs_coding = gapped_mapping.g_to_t_compensation(variant, validator, ori, hgvs_coding, rec_var)

    else:
        stored_hgvs_genomic_variant = hgvs_genomic
        suppress_c_normalization = 'false'
        gapped_alignment_warning = ''
        auto_info = ''
        genomic = fn.valstr(hgvs_genomic)

    # Create pseudo VCF based on amended hgvs_genomic
    hgvs_genomic_variant = hgvs_genomic
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = variant.reverse_normalizer.normalize(hgvs_genomic_variant)

    hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

    # Create vcf
    vcf_dict = vvHGVS.hgvs2vcf(reverse_normalized_hgvs_genomic, variant.primary_assembly,
                               variant.reverse_normalizer, validator.sf)
    chr = vcf_dict['chr']
    pos = vcf_dict['pos']
    ref = vcf_dict['ref']
    alt = vcf_dict['alt']

    # Create a VCF call
    vcf_component_list = [str(chr), str(pos), str(ref), (alt)]
    vcf_genomic = '-'.join(vcf_component_list)

    # DO NOT DELETE
    # Generate an end position
    end = str(int(pos) + len(ref) - 1)
    pos = str(pos)

    # DO NOT DELETE
    stored_hgvs_not_delins = validator.hp.parse_hgvs_variant(str(
        hgvs_genomic_5pr.ac) + ':' + hgvs_genomic_5pr.type + '.' + pos + '_' + end + 'del' + ref + 'ins' + alt)

    # Apply gap code to re-format hgvs_coding
    # Store the current hgvs:c. description
    saved_hgvs_coding = copy.deepcopy(hgvs_coding)

    # Get orientation of the gene wrt genome and a list of exons mapped to the genome
    ori = validator.tx_exons(tx_ac=saved_hgvs_coding.ac, alt_ac=hgvs_genomic_5pr.ac,
                        alt_aln_method=validator.alt_aln_method)
    orientation = int(ori[0]['alt_strand'])

    # Look for normalized variant options that do not match hgvs_coding
    hgvs_genomic = copy.deepcopy(hgvs_genomic_variant)
    if orientation == -1:
        # position genomic at its most 5 prime position
        try:
            query_genomic = variant.reverse_normalizer.normalize(hgvs_genomic)
        except:
            query_genomic = hgvs_genomic
        # Map to the transcript ant test for movement
        try:
            hgvs_seek_var = variant.evm.g_to_t(query_genomic, hgvs_coding.ac)
        except hgvs.exceptions.HGVSError as e:
            hgvs_seek_var = saved_hgvs_coding
        else:
            seek_var = fn.valstr(hgvs_seek_var)
            seek_ac = str(hgvs_seek_var.ac)
        if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                hgvs_coding.posedit.pos.start.base + hgvs_coding.posedit.pos.start.offset) and (
                hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                hgvs_coding.posedit.pos.end.base + hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
            pass
        else:
            hgvs_seek_var = saved_hgvs_coding

    elif orientation != -1:
        # position genomic at its most 3 prime position
        try:
            query_genomic = variant.hn.normalize(hgvs_genomic)
        except:
            query_genomic = hgvs_genomic
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

    # Loop out gap finding code under these circumstances!
    logger.warning("gap_compensation_2 = " + str(gap_compensation))
    if gap_compensation is True:
        logger.warning('g_to_t gap code 2 active')
        # is it in an exon?
        is_it_in_an_exon = 'no'
        for exon in ori:
            genomic_start = int(exon['alt_start_i'])
            genomic_end = int(exon['alt_end_i'])
            # Take from stored copy
            # hgvs_genomic_5pr = copy.deepcopy(stored_hgvs_genomic_5pr)
            if (
                    hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                    hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                is_it_in_an_exon = 'yes'
        if is_it_in_an_exon == 'yes':
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

                hard_fail = 'false'
                try:
                    tx_hgvs_not_delins = variant.no_norm_evm.g_to_n(hgvs_not_delins, saved_hgvs_coding.ac)
                except Exception as e:
                    if str(e) == 'start or end or both are beyond the bounds of transcript record':
                        tx_hgvs_not_delins = hgvs_coding
                        hard_fail = 'true'

                # Create normalized version of tx_hgvs_not_delins
                rn_tx_hgvs_not_delins = copy.deepcopy(tx_hgvs_not_delins)
                # Check for +ve base and adjust
                if re.search(r'\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search(r'\+',
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
                    rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                    # move tx start base to previous available non-offset base
                    rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
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
                    # rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base - 1
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
                    re_capture_tx_variant = []
                    for internal_possibility in hgvs_genomic_possibilities:

                        if internal_possibility == '':
                            continue

                        hgvs_t_possibility = validator.vm.g_to_t(internal_possibility, hgvs_coding.ac)
                        if hgvs_t_possibility.posedit.edit.type == 'ins':
                            try:
                                hgvs_t_possibility = validator.vm.c_to_n(hgvs_t_possibility)
                            except:
                                fn.exceptPass()
                            ins_ref = validator.sf.fetch_seq(hgvs_t_possibility.ac,
                                                        hgvs_t_possibility.posedit.pos.start.base - 1,
                                                        hgvs_t_possibility.posedit.pos.start.base + 1)
                            try:
                                hgvs_t_possibility = validator.vm.n_to_c(hgvs_t_possibility)
                            except:
                                fn.exceptPass()
                            hgvs_t_possibility.posedit.edit.ref = ins_ref
                            hgvs_t_possibility.posedit.edit.alt = ins_ref[
                                                                      0] + hgvs_t_possibility.posedit.edit.alt + \
                                                                  ins_ref[1]
                        if internal_possibility.posedit.edit.type == 'ins':
                            ins_ref = validator.sf.fetch_seq(internal_possibility.ac,
                                                        internal_possibility.posedit.pos.start.base - 1,
                                                        internal_possibility.posedit.pos.end.base)
                            internal_possibility.posedit.edit.ref = ins_ref
                            internal_possibility.posedit.edit.alt = ins_ref[
                                                                        0] + internal_possibility.posedit.edit.alt + \
                                                                    ins_ref[1]

                        if len(hgvs_t_possibility.posedit.edit.ref) < len(
                                internal_possibility.posedit.edit.ref):
                            gap_length = len(internal_possibility.posedit.edit.ref) - len(
                                hgvs_t_possibility.posedit.edit.ref)
                            re_capture_tx_variant = ['transcript', gap_length, hgvs_t_possibility]
                            hgvs_not_delins = internal_possibility
                            hgvs_genomic_5pr = internal_possibility
                            break

                    if re_capture_tx_variant != []:
                        try:
                            tx_hgvs_not_delins = validator.vm.c_to_n(re_capture_tx_variant[2])
                        except:
                            tx_hgvs_not_delins = re_capture_tx_variant[2]
                        disparity_deletion_in = re_capture_tx_variant[0:-1]
                    else:
                        pass

            # Final sanity checks
            try:
                validator.vm.g_to_t(hgvs_not_delins, tx_hgvs_not_delins.ac)
            except Exception as e:
                if str(e) == 'start or end or both are beyond the bounds of transcript record':
                    logger.warning(str(e))
                    return True
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
                        logger.warning(error)
                        return True
                    elif re.match('Normalization of intronic variants is not supported', error):
                        # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                        disparity_deletion_in = ['transcript', 'Requires Analysis']

            if hard_fail == 'true':
                disparity_deletion_in = ['false', 'false']

            # Recreate hgvs_genomic
            if disparity_deletion_in[0] == 'transcript':
                hgvs_genomic = hgvs_not_delins

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
                        g1 = validator.nr_vm.t_to_g(c1, hgvs_genomic.ac)
                        g3 = validator.nr_vm.t_to_g(c1, hgvs_genomic.ac)
                        g2 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
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
                            g1 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
                            g2 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                            g1.posedit.edit.alt = g1.posedit.edit.ref
                        else:
                            g1 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                            g2 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
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
                            g1 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
                            g2 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                            g2.posedit.edit.alt = g2.posedit.edit.ref
                        else:
                            g1 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                            g2 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
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
                            g1 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
                            g2 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                            g1.posedit.edit.alt = g1.posedit.edit.ref
                        else:
                            g1 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                            g2 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
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
                        c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                        c2.posedit.pos.end.offset = 0
                        c2.posedit.edit.ref = ''
                        c2.posedit.edit.alt = ''
                        if orientation != -1:
                            g1 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
                            g2 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                            g2.posedit.edit.alt = g2.posedit.edit.ref
                        else:
                            g1 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                            g2 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
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
                    hgvs_genomic.ac) + '\n'
                gapped_transcripts = gapped_transcripts + str(hgvs_refreshed_variant.ac) + ' '
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

            # Sort out equality to equality c. events where the code will add 2 additional bases
            if hgvs_coding.posedit.edit.type == 'identity' and hgvs_refreshed_variant.posedit.edit.type == 'identity':  # and len(hgvs_refreshed_variant.posedit.edit.ref) ==  (len(hgvs_coding.posedit.edit.ref) + 2):
                pass
            else:
                hgvs_coding = copy.deepcopy(hgvs_refreshed_variant)
            coding = fn.valstr(hgvs_coding)
            formatted_variant = coding

    # OBTAIN THE RefSeqGene coordinates
    # Attempt 1 = UTA
    sequences_for_tx = validator.hdp.get_tx_mapping_options(hgvs_coding.ac)
    recovered_rsg = []

    for sequence in sequences_for_tx:
        if re.match('^NG_', sequence[1]):
            recovered_rsg.append(sequence[1])
    recovered_rsg.sort()
    recovered_rsg.reverse()

    if len(recovered_rsg) > 0 and 'NG_' in recovered_rsg[0]:
        refseqgene_ac = recovered_rsg[0]
    else:
        refseqgene_ac = ''

    # Given the difficulties with mapping to and from RefSeqGenes, we now solely rely on UTA
    if refseqgene_ac != '':
        hgvs_refseq = validator.vm.t_to_g(hgvs_coding, refseqgene_ac)
        # Normalize the RefSeqGene Variant to the correct position
        try:
            hgvs_refseq = variant.hn.normalize(hgvs_refseq)
        except Exception as e:
            # if re.search('insertion length must be 1', error):
            hgvs_refseq = 'RefSeqGene record not available'
            refseq = 'RefSeqGene record not available'
            hgvs_refseq_ac = 'RefSeqGene record not available'
            pass
        else:
            refseq = fn.valstr(hgvs_refseq)
            hgvs_refseq_ac = hgvs_refseq.ac
    else:
        hgvs_refseq = 'RefSeqGene record not available'
        refseq = 'RefSeqGene record not available'
        hgvs_refseq_ac = 'RefSeqGene record not available'

    # Predicted effect on protein
    protein_dict = validator.myc_to_p(hgvs_coding, variant.evm, re_to_p=False)
    if protein_dict['error'] == '':
        hgvs_protein = protein_dict['hgvs_protein']
        protein = str(hgvs_protein)
    else:
        error = protein_dict['error']
        variant.warnings += ': ' + str(error)
        if error == 'Cannot identify an in-frame Termination codon in the variant mRNA sequence':
            hgvs_protein = protein_dict['hgvs_protein']
            protein = str(hgvs_protein)
        else:
            logger.error(error)
            return True

    # Gene orientation wrt genome
    ori = validator.tx_exons(tx_ac=hgvs_coding.ac, alt_ac=hgvs_genomic.ac,
                        alt_aln_method=validator.alt_aln_method)
    ori = int(ori[0]['alt_strand'])

    # Look for normalized variant options that do not match hgvs_coding
    # boundary crossing normalization
    # Re-Save the required variants
    hgvs_seek_var = copy.deepcopy(hgvs_coding)
    saved_hgvs_coding = copy.deepcopy(hgvs_coding)

    if ori == -1:
        # position genomic at its most 5 prime position
        try:
            query_genomic = variant.reverse_normalizer.normalize(hgvs_genomic)
        except:
            query_genomic = hgvs_genomic
        # Map to the transcript and test for movement
        try:
            hgvs_seek_var = variant.evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
        except hgvs.exceptions.HGVSError as e:
            hgvs_seek_var = saved_hgvs_coding
        else:
            seek_var = fn.valstr(hgvs_seek_var)
            seek_ac = str(hgvs_seek_var.ac)
        if saved_hgvs_coding.posedit.edit.type != hgvs_seek_var.posedit.edit.type:
            rec_var = 'false'
            hgvs_seek_var = saved_hgvs_coding
            seek_var = fn.valstr(hgvs_seek_var)
            seek_ac = str(hgvs_seek_var.ac)
        elif suppress_c_normalization == 'true':
            rec_var = 'false'
            hgvs_seek_var = saved_hgvs_coding
            seek_var = fn.valstr(hgvs_seek_var)
            seek_ac = str(hgvs_seek_var.ac)
        elif (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
            try:
                automap = fn.valstr(saved_hgvs_coding) + ' normalized to ' + fn.valstr(hgvs_seek_var)
                hgvs_coding = hgvs_seek_var
                coding = fn.valstr(hgvs_coding)
                variant.warnings += ': ' + automap
                rng = variant.hn.normalize(query_genomic)
            except NotImplementedError:
                pass
            try:
                c_for_p = validator.vm.g_to_t(rng, hgvs_coding.ac)
            except hgvs.exceptions.HGVSInvalidIntervalError as e:
                c_for_p = seek_var
            try:
                # Predicted effect on protein
                protein_dict = validator.myc_to_p(c_for_p, variant.evm, re_to_p=False)
                if protein_dict['error'] == '':
                    hgvs_protein = protein_dict['hgvs_protein']
                    protein = str(hgvs_protein)
                else:
                    error = protein_dict['error']
                    if error == 'Cannot identify an in-frame Termination codon in the variant mRNA sequence':
                        hgvs_protein = protein_dict['hgvs_protein']
                        variant.warnings += ': ' + str(error)
                # Replace protein description in vars table
                protein = str(hgvs_protein)
            except NotImplementedError:
                fn.exceptPass()
        else:
            # Double check protein position by normalize genomic, and normalize back to c. for normalize or not to normalize issue
            coding = fn.valstr(hgvs_coding)

    elif ori != -1:
        # position genomic at its most 3 prime position
        try:
            query_genomic = variant.hn.normalize(hgvs_genomic)
        except:
            query_genomic = hgvs_genomic
        # Map to the transcript and test for movement
        try:
            hgvs_seek_var = variant.evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
        except hgvs.exceptions.HGVSError as e:
            hgvs_seek_var = saved_hgvs_coding
        else:
            seek_var = fn.valstr(hgvs_seek_var)
            seek_ac = str(hgvs_seek_var.ac)
        if saved_hgvs_coding.posedit.edit.type != hgvs_seek_var.posedit.edit.type:
            rec_var = 'false'
            hgvs_seek_var = saved_hgvs_coding
            seek_var = fn.valstr(hgvs_seek_var)
            seek_ac = str(hgvs_seek_var.ac)
        elif suppress_c_normalization == 'true':
            rec_var = 'false'
            hgvs_seek_var = saved_hgvs_coding
            seek_var = fn.valstr(hgvs_seek_var)
            seek_ac = str(hgvs_seek_var.ac)
        elif (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
            try:
                automap = fn.valstr(saved_hgvs_coding) + ' normalized to ' + fn.valstr(hgvs_seek_var)
                hgvs_coding = hgvs_seek_var
                coding = fn.valstr(hgvs_coding)
                variant.warnings += ': ' + automap
            except NotImplementedError:
                fn.exceptPass()
        else:
            # Double check protein position by reverse_norm genomic, and normalize back to c. for normalize or not to normalize issue
            coding = fn.valstr(hgvs_coding)
            rng = variant.reverse_normalizer.normalize(query_genomic)
            try:
                # Diagram where - = intron and E = Exon

                # 3 prime
                # ---------EEEEEEEEEEEEEEEEE-----------
                #          <
                # Result, normalize of new variant will baulk at intronic
                # 5 prime
                #                          <
                # Result, normalize of new variant will be happy
                c_for_p = validator.vm.g_to_t(rng, hgvs_coding.ac)
                try:
                    variant.hn.normalize(c_for_p)
                except hgvs.exceptions.HGVSError as e:
                    fn.exceptPass()
                else:
                    # hgvs_protein = va_func.protein(str(c_for_p), variant.evm, hp)
                    protein_dict = validator.myc_to_p(c_for_p, variant.evm, re_to_p=False)
                    if protein_dict['error'] == '':
                        hgvs_protein = protein_dict['hgvs_protein']
                        protein = str(hgvs_protein)
                    else:
                        error = protein_dict['error']
                        if error == 'Cannot identify an in-frame Termination codon in the variant mRNA sequence':
                            hgvs_protein = protein_dict['hgvs_protein']
                            variant.warnings += ': ' + str(error)
                    # Replace protein description in vars table
                    protein = str(hgvs_protein)
            except Exception:
                fn.exceptPass()

    # Check for up-to-date transcript version
    updated_transcript_variant = 'None'
    tx_id_info = validator.hdp.get_tx_identity_info(hgvs_coding.ac)
    uta_gene_symbol = tx_id_info[6]
    tx_for_gene = validator.hdp.get_tx_for_gene(uta_gene_symbol)
    ac_root, ac_version = hgvs_coding.ac.split('.')
    version_tracking = '0'
    update = ''
    for accession in tx_for_gene:
        try:
            if re.match(ac_root, accession[3]):
                query_version = accession[3].split('.')[1]
                if int(query_version) > int(ac_version) and int(query_version) > int(
                        version_tracking):
                    version_tracking = query_version
                    update = accession[3]
        except ValueError:
            fn.exceptPass()

    if update != '':
        hgvs_updated = copy.deepcopy(hgvs_coding)
        hgvs_updated.ac = update
        try:
            validator.vr.validate(hgvs_updated)
        # Updated reference sequence
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if re.search('does not agree with reference sequence', str(error)):
                match = re.findall(r'\(([GATC]+)\)', error)
                new_ref = match[1]
                hgvs_updated.posedit.edit.ref = new_ref
                validator.vr.validate(hgvs_updated)
                updated_transcript_variant = hgvs_updated
            else:
                pass
        updated_transcript_variant = hgvs_updated
        variant.warnings += ': ' + 'A more recent version of the selected reference sequence ' + hgvs_coding.ac + ' is available (' + updated_transcript_variant.ac + ')' + ': ' + str(
            updated_transcript_variant) + ' MUST be fully validated prior to use in reports: select_variants=' + fn.valstr(
            updated_transcript_variant)

    variant.coding = str(hgvs_coding)
    variant.genomic_r = str(hgvs_refseq)
    variant.genomic_g = str(hgvs_genomic)
    variant.protein = str(hgvs_protein)

    # if gap_compensation is True:
    #     variant.test_stash_tx_left = test_stash_tx_left
    #     variant.test_stash_tx_right = test_stash_tx_right

    return False