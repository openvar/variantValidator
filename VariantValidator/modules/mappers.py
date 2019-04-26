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

        gap_mapper = gapped_mapping.GapMapper(variant, validator)

        data, nw_rel_var = gap_mapper.gapped_g_to_c(rel_var)

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

    plus = re.compile(r"\d\+\d")  # finds digit + digit
    minus = re.compile(r"\d\-\d")  # finds digit - digit

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

                    # Kill current line and append for re-submission
                    # Tag the line so that it is not written out
                    variant.write = False
                    # Set the values and append to batch_list
                    hgvs_vt = validator.hp.parse_hgvs_variant(str(post_var))
                    assert str(hgvs_vt) == str(post_var)
                    query = Variant(variant.original, quibble=fn.valstr(hgvs_vt), warnings=automap, primary_assembly=variant.primary_assembly, order=variant.order)
                    validator.batch_list.append(query)

        else:  # del not in formatted_variant
            if ':r.' in variant.trapped:
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

                    # Kill current line and append for re-submission
                    # Tag the line so that it is not written out
                    variant.write = False
                    # Set the values and append to batch_list
                    hgvs_vt = validator.hp.parse_hgvs_variant(str(post_var))
                    assert str(hgvs_vt) == str(post_var)
                    query = Variant(variant.original, quibble=fn.valstr(hgvs_vt), warnings=automap, primary_assembly=variant.primary_assembly, order=variant.order)
                    validator.batch_list.append(query)

    # If cck not true
    elif ':r.' in variant.trapped:
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

    elif ':g.' in input:
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

    # Genomic sequence
    hgvs_genomic = validator.myevm_t_to_g(hgvs_coding, variant.no_norm_evm, variant.primary_assembly, variant.hn)

    # Create gap_mapper object instance
    gap_mapper = gapped_mapping.GapMapper(variant, validator)

    # --- GAP MAPPING 1 ---
    # Loop out gap finding code under these circumstances!
    if gap_compensation is True:
        hgvs_genomic, suppress_c_normalization, hgvs_coding = gap_mapper.g_to_t_compensation(ori, hgvs_coding, rec_var)

    else:
        suppress_c_normalization = 'false'

    # Create pseudo VCF based on amended hgvs_genomic
    # hgvs_genomic_variant = hgvs_genomic
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = variant.reverse_normalizer.normalize(hgvs_genomic)

    # Get orientation of the gene wrt genome and a list of exons mapped to the genome
    ori = validator.tx_exons(tx_ac=hgvs_coding.ac, alt_ac=reverse_normalized_hgvs_genomic.ac,
                        alt_aln_method=validator.alt_aln_method)

    # --- GAP MAPPING 2 ---
    # Loop out gap finding code under these circumstances!
    logger.warning("gap_compensation_2 = " + str(gap_compensation))
    if gap_compensation is True:
        hgvs_coding = gap_mapper.g_to_t_gapped_mapping_stage2(ori, hgvs_coding, hgvs_genomic)

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
    else:
        hgvs_refseq = 'RefSeqGene record not available'

    # Predicted effect on protein
    protein_dict = validator.myc_to_p(hgvs_coding, variant.evm, re_to_p=False)
    if protein_dict['error'] == '':
        hgvs_protein = protein_dict['hgvs_protein']
    else:
        error = protein_dict['error']
        variant.warnings += ': ' + str(error)
        if error == 'Cannot identify an in-frame Termination codon in the variant mRNA sequence':
            hgvs_protein = protein_dict['hgvs_protein']
        else:
            logger.error(error)
            return True

    # Gene orientation wrt genome
    ori = validator.tx_exons(tx_ac=hgvs_coding.ac, alt_ac=hgvs_genomic.ac,
                        alt_aln_method=validator.alt_aln_method)
    ori = int(ori[0]['alt_strand'])

    # Look for normalized variant options that do not match hgvs_coding
    # boundary crossing normalization
    hgvs_seek_var, query_genomic = gap_mapper.get_hgvs_seek_var(hgvs_genomic, hgvs_coding,
                                                                ori=ori, with_query_genomic=True)

    if hgvs_coding.posedit.edit.type != hgvs_seek_var.posedit.edit.type:
        pass
    elif suppress_c_normalization == 'true':
        pass
    elif (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
            hgvs_coding.posedit.pos.start.base + hgvs_coding.posedit.pos.start.offset) and (
            hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
            hgvs_coding.posedit.pos.end.base + hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
        try:
            automap = fn.valstr(hgvs_coding) + ' normalized to ' + fn.valstr(hgvs_seek_var)
            hgvs_coding = hgvs_seek_var
            variant.warnings += ': ' + automap
        except NotImplementedError:
            fn.exceptPass()
        if ori == -1:
            rng = variant.hn.normalize(query_genomic)
            try:
                c_for_p = validator.vm.g_to_t(rng, hgvs_coding.ac)
            except hgvs.exceptions.HGVSInvalidIntervalError as e:
                c_for_p = fn.valstr(hgvs_seek_var)
            try:
                # Predicted effect on protein
                protein_dict = validator.myc_to_p(c_for_p, variant.evm, re_to_p=False)
                if protein_dict['error'] == '':
                    hgvs_protein = protein_dict['hgvs_protein']
                else:
                    error = protein_dict['error']
                    if error == 'Cannot identify an in-frame Termination codon in the variant mRNA sequence':
                        hgvs_protein = protein_dict['hgvs_protein']
                        variant.warnings += ': ' + str(error)
            except NotImplementedError:
                fn.exceptPass()
    elif ori == 1:
        # Double check protein position by reverse_norm genomic, and normalize back to c. for normalize or not to normalize issue
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
                else:
                    error = protein_dict['error']
                    if error == 'Cannot identify an in-frame Termination codon in the variant mRNA sequence':
                        hgvs_protein = protein_dict['hgvs_protein']
                        variant.warnings += ': ' + str(error)
                # Replace protein description in vars table
        except Exception:
            fn.exceptPass()

    # Check for up-to-date transcript version
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
        variant.warnings += ': ' + 'A more recent version of the selected reference sequence ' + hgvs_coding.ac + ' is available (' + updated_transcript_variant.ac + ')' + ': ' + str(
            updated_transcript_variant) + ' MUST be fully validated prior to use in reports: select_variants=' + fn.valstr(
            updated_transcript_variant)

    variant.coding = str(hgvs_coding)
    variant.genomic_r = str(hgvs_refseq)
    variant.genomic_g = str(hgvs_genomic)
    variant.protein = str(hgvs_protein)

    return False


def final_tx_to_multiple_genomic(variant, validator, tx_variant, rec_var):

    warnings = ''

    # Multiple genomic variants
    # multi_gen_vars = []
    hgvs_coding = validator.hp.parse_hgvs_variant(str(tx_variant))
    # Gap gene black list
    try:
        gene_symbol = validator.db.get_gene_symbol_from_transcriptID(hgvs_coding.ac)
    except Exception:
        fn.exceptPass()
    else:
        # If the gene symbol is not in the list, the value False will be returned
        gap_compensation = vvChromosomes.gap_black_list(gene_symbol)

    # Look for variants spanning introns
    try:
        hgvs_coding = variant.hn.normalize(hgvs_coding)
    except hgvs.exceptions.HGVSUnsupportedOperationError as e:
        error = str(e)
        if re.search('boundary', str(error)) or re.search('spanning', str(error)):
            gap_compensation = False
        else:
            pass
    except hgvs.exceptions.HGVSError:
        fn.exceptPass()

    # Warn gap code status
    logger.warning("gap_compensation_3 = " + str(gap_compensation))
    multi_g = []
    multi_list = []
    mapping_options = validator.hdp.get_tx_mapping_options(hgvs_coding.ac)
    for alt_chr in mapping_options:
        if (re.match('NC_', alt_chr[1]) or re.match('NT_', alt_chr[1]) or re.match('NW_',
                                                                                   alt_chr[1])) and \
                alt_chr[2] == validator.alt_aln_method:
            multi_list.append(alt_chr[1])

    for alt_chr in multi_list:
        try:
            # Re set ori
            ori = validator.tx_exons(tx_ac=hgvs_coding.ac, alt_ac=alt_chr,
                                alt_aln_method=validator.alt_aln_method)
            orientation = int(ori[0]['alt_strand'])
            hgvs_alt_genomic = validator.myvm_t_to_g(hgvs_coding, alt_chr, variant.no_norm_evm, variant.hn)
            # Set hgvs_genomic accordingly
            hgvs_genomic = copy.deepcopy(hgvs_alt_genomic)

            # genomic_possibilities
            # 1. take the simple 3 pr normalized hgvs_genomic
            # 2. Lock in hgvs_genomic at its most 5 prime position wrt genome
            hgvs_genomic_possibilities = []

            gap_mapper = gapped_mapping.GapMapper(variant, validator)

            # Loop out gap code under these circumstances!
            if gap_compensation is True:

                hgvs_alt_genomic, hgvs_coding = gap_mapper.g_to_t_gap_compensation_version3(hgvs_alt_genomic, hgvs_coding, ori, alt_chr, rec_var)

                # logger.warning('g_to_t gap code 3 active')
                # rn_hgvs_genomic = variant.reverse_normalizer.normalize(hgvs_alt_genomic)
                # hgvs_genomic_possibilities.append(rn_hgvs_genomic)
                # if orientation != -1:
                #     try:
                #         chromosome_normalized_hgvs_coding = variant.reverse_normalizer.normalize(
                #             hgvs_coding)
                #     except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                #         chromosome_normalized_hgvs_coding = hgvs_coding
                # else:
                #     try:
                #         chromosome_normalized_hgvs_coding = variant.hn.normalize(hgvs_coding)
                #     except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                #         error = str(e)
                #         chromosome_normalized_hgvs_coding = hgvs_coding
                #
                # most_3pr_hgvs_genomic = validator.myvm_t_to_g(chromosome_normalized_hgvs_coding,
                #                                          alt_chr,
                #                                          variant.no_norm_evm, variant.hn)
                # hgvs_genomic_possibilities.append(most_3pr_hgvs_genomic)
                #
                # # First to the right
                # hgvs_stash = copy.deepcopy(hgvs_coding)
                # try:
                #     hgvs_stash = variant.no_norm_evm.c_to_n(hgvs_stash)
                # except:
                #     fn.exceptPass()
                # try:
                #     stash_ac = hgvs_stash.ac
                #     stash_dict = vvHGVS.hard_right_hgvs2vcf(hgvs_stash, variant.primary_assembly, variant.hn, validator.sf)
                #     stash_pos = int(stash_dict['pos'])
                #     stash_ref = stash_dict['ref']
                #     stash_alt = stash_dict['alt']
                #     # Generate an end position
                #     stash_end = str(stash_pos + len(stash_ref) - 1)
                #     # make a not real deletion insertion
                #     stash_hgvs_not_delins = validator.hp.parse_hgvs_variant(
                #         stash_ac + ':' + hgvs_stash.type + '.' + str(
                #             stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)
                #     try:
                #         stash_hgvs_not_delins = variant.no_norm_evm.n_to_c(stash_hgvs_not_delins)
                #     except:
                #         fn.exceptPass()
                #         # Store a tx copy for later use
                #     test_stash_tx_right = copy.deepcopy(stash_hgvs_not_delins)
                #     stash_genomic = validator.myvm_t_to_g(test_stash_tx_right, hgvs_alt_genomic.ac,
                #                                      variant.no_norm_evm, variant.hn)
                #     # Stash the outputs if required
                #     # test variants = NC_000006.11:g.90403795G= (causes double identity)
                #     #                 NC_000002.11:g.73675227_73675228insCTC (? incorrect assumed insertion position)
                #     #                 NC_000003.11:g.14561629_14561630GC= NC_000003.11:g.14561629_14561630insG (Odd gap position)
                #     # if test_stash_tx_right.posedit.edit.type == 'identity' and stash_genomic.posedit.edit.type == 'identity':
                #     # pass
                #     if len(test_stash_tx_right.posedit.edit.ref) == ((
                #                                                              stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                #         stash_tx_right = test_stash_tx_right
                #         if hasattr(test_stash_tx_right.posedit.edit,
                #                    'alt') and test_stash_tx_right.posedit.edit.alt is not None:
                #             alt = test_stash_tx_right.posedit.edit.alt
                #         else:
                #             alt = ''
                #         if hasattr(stash_genomic.posedit.edit,
                #                    'alt') and stash_genomic.posedit.edit.alt is not None:
                #             g_alt = stash_genomic.posedit.edit.alt
                #         else:
                #             g_alt = ''
                #         if (len(alt) - (
                #                 test_stash_tx_right.posedit.pos.end.base - test_stash_tx_right.posedit.pos.start.base) + 1) != (
                #                 len(g_alt) - (
                #                 stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                #             hgvs_genomic_possibilities.append(stash_genomic)
                #         else:
                #             hgvs_genomic_possibilities.append('')
                #     elif test_stash_tx_right.posedit.edit.type == 'identity':
                #         reform_ident = str(test_stash_tx_right).split(':')[0]
                #         reform_ident = reform_ident + ':c.' + str(test_stash_tx_right.posedit.pos) + 'del' + str(
                #             test_stash_tx_right.posedit.edit.ref)  # + 'ins' + str(test_stash_tx_right.posedit.edit.alt)
                #         hgvs_reform_ident = validator.hp.parse_hgvs_variant(reform_ident)
                #         try:
                #             variant.hn.normalize(hgvs_reform_ident)
                #         except hgvs.exceptions.HGVSError as e:
                #             error = str(e)
                #             if re.search('spanning the exon-intron boundary', error):
                #                 stash_tx_right = test_stash_tx_right
                #                 hgvs_genomic_possibilities.append('')
                #         else:
                #             stash_tx_right = test_stash_tx_right
                #             hgvs_genomic_possibilities.append(stash_genomic)
                #     else:
                #         try:
                #             variant.hn.normalize(test_stash_tx_right)
                #         except hgvs.exceptions.HGVSUnsupportedOperationError:
                #             hgvs_genomic_possibilities.append('')
                #         else:
                #             stash_tx_right = test_stash_tx_right
                #             hgvs_genomic_possibilities.append(stash_genomic)
                # except hgvs.exceptions.HGVSError as e:
                #     fn.exceptPass()
                # except ValueError:
                #     fn.exceptPass()
                #
                # # Then to the left
                # hgvs_stash = copy.deepcopy(hgvs_coding)
                # try:
                #     hgvs_stash = variant.no_norm_evm.c_to_n(hgvs_stash)
                # except:
                #     fn.exceptPass()
                # try:
                #     stash_ac = hgvs_stash.ac
                #     stash_dict = vvHGVS.hard_left_hgvs2vcf(hgvs_stash, variant.primary_assembly,
                #                                            variant.reverse_normalizer, validator.sf)
                #     stash_pos = int(stash_dict['pos'])
                #     stash_ref = stash_dict['ref']
                #     stash_alt = stash_dict['alt']
                #     # Generate an end position
                #     stash_end = str(stash_pos + len(stash_ref) - 1)
                #     # make a not real deletion insertion
                #     stash_hgvs_not_delins = validator.hp.parse_hgvs_variant(
                #         stash_ac + ':' + hgvs_stash.type + '.' + str(
                #             stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)
                #     try:
                #         stash_hgvs_not_delins = variant.no_norm_evm.n_to_c(stash_hgvs_not_delins)
                #     except:
                #         fn.exceptPass()
                #         # Store a tx copy for later use
                #     test_stash_tx_left = copy.deepcopy(stash_hgvs_not_delins)
                #     stash_genomic = validator.myvm_t_to_g(test_stash_tx_left, hgvs_alt_genomic.ac,
                #                                      variant.no_norm_evm, variant.hn)
                #     # Stash the outputs if required
                #     # test variants = NC_000006.11:g.90403795G= (causes double identity)
                #     #                 NC_000002.11:g.73675227_73675228insCTC
                #     #                 NC_000003.11:g.14561629_14561630GC= NC_000003.11:g.14561629_14561630insG (Odd gap position)
                #     # if test_stash_tx_left.posedit.edit.type == 'identity' and stash_genomic.posedit.edit.type == 'identity':
                #     # pass
                #     if len(test_stash_tx_left.posedit.edit.ref) == ((
                #                                                             stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):  # len(stash_genomic.posedit.edit.ref):
                #         stash_tx_left = test_stash_tx_left
                #         if hasattr(test_stash_tx_left.posedit.edit,
                #                    'alt') and test_stash_tx_left.posedit.edit.alt is not None:
                #             alt = test_stash_tx_left.posedit.edit.alt
                #         else:
                #             alt = ''
                #         if hasattr(stash_genomic.posedit.edit,
                #                    'alt') and stash_genomic.posedit.edit.alt is not None:
                #             g_alt = stash_genomic.posedit.edit.alt
                #         else:
                #             g_alt = ''
                #         if (len(alt) - (
                #                 test_stash_tx_left.posedit.pos.end.base - test_stash_tx_left.posedit.pos.start.base) + 1) != (
                #                 len(g_alt) - (
                #                 stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                #             hgvs_genomic_possibilities.append(stash_genomic)
                #         else:
                #             hgvs_genomic_possibilities.append('')
                #     elif test_stash_tx_left.posedit.edit.type == 'identity':
                #         reform_ident = str(test_stash_tx_left).split(':')[0]
                #         reform_ident = reform_ident + ':c.' + str(test_stash_tx_left.posedit.pos) + 'del' + str(
                #             test_stash_tx_left.posedit.edit.ref)  # + 'ins' + str(test_stash_tx_left.posedit.edit.alt)
                #         hgvs_reform_ident = validator.hp.parse_hgvs_variant(reform_ident)
                #         try:
                #             variant.hn.normalize(hgvs_reform_ident)
                #         except hgvs.exceptions.HGVSError as e:
                #             error = str(e)
                #             if re.search('spanning the exon-intron boundary', error):
                #                 stash_tx_left = test_stash_tx_left
                #                 hgvs_genomic_possibilities.append('')
                #         else:
                #             stash_tx_left = test_stash_tx_left
                #             hgvs_genomic_possibilities.append(stash_genomic)
                #     else:
                #         try:
                #             variant.hn.normalize(test_stash_tx_left)
                #         except hgvs.exceptions.HGVSUnsupportedOperationError:
                #             hgvs_genomic_possibilities.append('')
                #         else:
                #             stash_tx_left = test_stash_tx_left
                #             hgvs_genomic_possibilities.append(stash_genomic)
                # except hgvs.exceptions.HGVSError as e:
                #     fn.exceptPass()
                # except ValueError:
                #     fn.exceptPass()
                #
                # # direct mapping from reverse_normalized transcript insertions in the delins format
                # try:
                #     if hgvs_coding.posedit.edit.type == 'ins':
                #         most_5pr_hgvs_transcript_variant = copy.deepcopy(hgvs_coding)
                #         most_3pr_hgvs_transcript_variant = variant.reverse_normalizer.normalize(hgvs_coding)
                #         try:
                #             n_3pr = validator.vm.c_to_n(most_3pr_hgvs_transcript_variant)
                #             n_5pr = validator.vm.c_to_n(most_5pr_hgvs_transcript_variant)
                #         except:
                #             n_3pr = most_3pr_hgvs_transcript_variant
                #             n_5pr = most_5pr_hgvs_transcript_variant
                #         # Make into a delins by adding the ref bases to the variant ref and alt
                #         pr3_ref = validator.sf.fetch_seq(hgvs_coding.ac, n_3pr.posedit.pos.start.base - 1,
                #                                     n_3pr.posedit.pos.end.base)
                #         pr5_ref = validator.sf.fetch_seq(hgvs_coding.ac, n_5pr.posedit.pos.start.base - 1,
                #                                     n_5pr.posedit.pos.end.base)
                #         most_3pr_hgvs_transcript_variant.posedit.edit.ref = pr3_ref
                #         most_5pr_hgvs_transcript_variant.posedit.edit.ref = pr5_ref
                #         most_3pr_hgvs_transcript_variant.posedit.edit.alt = pr3_ref[
                #                                                                 0] + most_3pr_hgvs_transcript_variant.posedit.edit.alt + \
                #                                                             pr3_ref[1]
                #         most_5pr_hgvs_transcript_variant.posedit.edit.alt = pr5_ref[
                #                                                                 0] + most_5pr_hgvs_transcript_variant.posedit.edit.alt + \
                #                                                             pr5_ref[1]
                #         # Map to the genome
                #         genomic_from_most_3pr_hgvs_transcript_variant = validator.vm.t_to_g(
                #             most_3pr_hgvs_transcript_variant, hgvs_genomic.ac)
                #         genomic_from_most_5pr_hgvs_transcript_variant = validator.vm.t_to_g(
                #             most_5pr_hgvs_transcript_variant, hgvs_genomic.ac)
                #
                #         # Normalize - If the variant spans a gap it should then form a static genomic variant
                #         try:
                #             genomic_from_most_3pr_hgvs_transcript_variant = variant.hn.normalize(
                #                 genomic_from_most_3pr_hgvs_transcript_variant)
                #         except hgvs.exceptions.HGVSInvalidVariantError as e:
                #             error = str(e)
                #             if error == 'base start position must be <= end position':
                #                 start = genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.start.base
                #                 end = genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.end.base
                #                 genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.start.base = end
                #                 genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.end.base = start
                #                 genomic_from_most_3pr_hgvs_transcript_variant = variant.hn.normalize(
                #                     genomic_from_most_3pr_hgvs_transcript_variant)
                #         try:
                #             genomic_from_most_5pr_hgvs_transcript_variant = variant.hn.normalize(
                #                 genomic_from_most_5pr_hgvs_transcript_variant)
                #         except hgvs.exceptions.HGVSInvalidVariantError as e:
                #             error = str(e)
                #             if error == 'base start position must be <= end position':
                #                 start = genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.start.base
                #                 end = genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.end.base
                #                 genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.start.base = end
                #                 genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.end.base = start
                #                 genomic_from_most_5pr_hgvs_transcript_variant = variant.hn.normalize(
                #                     genomic_from_most_5pr_hgvs_transcript_variant)
                #
                #         try:
                #             if genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.alt is None:
                #                 genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.alt = ''
                #         except Exception as e:
                #             if str(e) == "'Dup' object has no attribute 'alt'":
                #                 genomic_from_most_3pr_hgvs_transcript_variant_delins_from_dup = genomic_from_most_3pr_hgvs_transcript_variant.ac + ':' + genomic_from_most_3pr_hgvs_transcript_variant.type + '.' + str(
                #                     genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.start.base) + '_' + str(
                #                     genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.end.base) + 'del' + genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.ref + 'ins' + genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.ref + genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.ref
                #                 genomic_from_most_3pr_hgvs_transcript_variant = validator.hp.parse_hgvs_variant(
                #                     genomic_from_most_3pr_hgvs_transcript_variant_delins_from_dup)
                #
                #         try:
                #             if most_3pr_hgvs_transcript_variant.posedit.edit.alt is None:
                #                 most_3pr_hgvs_transcript_variant.posedit.edit.alt = ''
                #         except Exception as e:
                #             if str(e) == "'Dup' object has no attribute 'alt'":
                #                 most_3pr_hgvs_transcript_variant_delins_from_dup = most_3pr_hgvs_transcript_variant.ac + ':' + most_3pr_hgvs_transcript_variant.type + '.' + str(
                #                     most_3pr_hgvs_transcript_variant.posedit.pos.start.base) + '_' + str(
                #                     most_3pr_hgvs_transcript_variant.posedit.pos.end.base) + 'del' + most_3pr_hgvs_transcript_variant.posedit.edit.ref + 'ins' + most_3pr_hgvs_transcript_variant.posedit.edit.ref + most_3pr_hgvs_transcript_variant.posedit.edit.ref
                #                 most_3pr_hgvs_transcript_variant = validator.hp.parse_hgvs_variant(
                #                     most_3pr_hgvs_transcript_variant_delins_from_dup)
                #
                #         try:
                #             if genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.alt is None:
                #                 genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.alt = ''
                #         except Exception as e:
                #             if str(e) == "'Dup' object has no attribute 'alt'":
                #                 genomic_from_most_5pr_hgvs_transcript_variant_delins_from_dup = genomic_from_most_5pr_hgvs_transcript_variant.ac + ':' + genomic_from_most_5pr_hgvs_transcript_variant.type + '.' + str(
                #                     genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.start.base) + '_' + str(
                #                     genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.end.base) + 'del' + genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.ref + 'ins' + genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.ref + genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.ref
                #                 genomic_from_most_5pr_hgvs_transcript_variant = validator.hp.parse_hgvs_variant(
                #                     genomic_from_most_5pr_hgvs_transcript_variant_delins_from_dup)
                #
                #         try:
                #             if most_5pr_hgvs_transcript_variant.posedit.edit.alt is None:
                #                 most_5pr_hgvs_transcript_variant.posedit.edit.alt = ''
                #         except Exception as e:
                #             if str(e) == "'Dup' object has no attribute 'alt'":
                #                 most_5pr_hgvs_transcript_variant_delins_from_dup = most_5pr_hgvs_transcript_variant.ac + ':' + most_5pr_hgvs_transcript_variant.type + '.' + str(
                #                     most_5pr_hgvs_transcript_variant.posedit.pos.start.base) + '_' + str(
                #                     most_5pr_hgvs_transcript_variant.posedit.pos.end.base) + 'del' + most_5pr_hgvs_transcript_variant.posedit.edit.ref + 'ins' + most_5pr_hgvs_transcript_variant.posedit.edit.ref + most_5pr_hgvs_transcript_variant.posedit.edit.ref
                #                 most_5pr_hgvs_transcript_variant = validator.hp.parse_hgvs_variant(
                #                     most_5pr_hgvs_transcript_variant_delins_from_dup)
                #
                #         if len(
                #                 genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.alt) < len(
                #             most_3pr_hgvs_transcript_variant.posedit.edit.alt):
                #             hgvs_genomic_possibilities.append(
                #                 genomic_from_most_3pr_hgvs_transcript_variant)
                #         if len(
                #                 genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.alt) < len(
                #             most_5pr_hgvs_transcript_variant.posedit.edit.alt):
                #             hgvs_genomic_possibilities.append(
                #                 genomic_from_most_5pr_hgvs_transcript_variant)
                #
                # except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                #     error = str(e)
                #     if re.match('Normalization of intronic variants is not supported',
                #                 error) or re.match(
                #         'Unsupported normalization of variants spanning the exon-intron boundary',
                #         error):
                #         pass
                #     fn.exceptPass()
                #
                # # Set variables for problem specific warnings
                # gapped_alignment_warning = ''
                # corrective_action_taken = ''
                # gapped_transcripts = ''
                # auto_info = ''
                #
                # # Mark as not disparity detected
                # disparity_deletion_in = ['false', 'false']
                # # Loop through to see if a gap can be located
                # possibility_counter = 0
                # for possibility in hgvs_genomic_possibilities:
                #     possibility_counter = possibility_counter + 1
                #     # Loop out stash possibilities which will not spot gaps so are empty
                #     if possibility == '':
                #         continue
                #
                #     # Use VCF generation code to push hgvs_genomic as for 5 prime as possible to uncover gaps
                #     hgvs_genomic_variant = possibility
                #     stored_hgvs_genomic_variant = copy.deepcopy(hgvs_genomic_variant)
                #
                #     # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
                #     try:
                #         reverse_normalized_hgvs_genomic = variant.reverse_normalizer.normalize(
                #             hgvs_genomic_variant)
                #     except hgvs.exceptions.HGVSError as e:
                #         # Strange error caused by gap in genomic
                #         error = str(e)
                #         if re.search('base start position must be <= end position', error):
                #             if hgvs_genomic.posedit.edit.type == 'delins':
                #                 start = hgvs_genomic.posedit.pos.start.base
                #                 end = hgvs_genomic.posedit.pos.end.base
                #                 lhb = validator.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                #                 rhb = validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                #                 hgvs_genomic.posedit.edit.ref = lhb + rhb
                #                 hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                #                 hgvs_genomic.posedit.pos.start.base = end
                #                 hgvs_genomic.posedit.pos.end.base = start
                #                 reverse_normalized_hgvs_genomic = variant.reverse_normalizer.normalize(
                #                     hgvs_genomic)
                #             if hgvs_genomic.posedit.edit.type == 'del':
                #                 start = hgvs_genomic.posedit.pos.start.base
                #                 end = hgvs_genomic.posedit.pos.end.base
                #                 lhb = validator.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                #                 rhb = validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                #                 hgvs_genomic.posedit.edit.ref = lhb + rhb
                #                 hgvs_genomic.posedit.edit.alt = lhb + rhb
                #                 hgvs_genomic.posedit.pos.start.base = end
                #                 hgvs_genomic.posedit.pos.end.base = start
                #                 reverse_normalized_hgvs_genomic = variant.reverse_normalizer.normalize(
                #                     hgvs_genomic)
                #         if re.search('insertion length must be 1', error):
                #             if hgvs_genomic.posedit.edit.type == 'ins':
                #                 start = hgvs_genomic.posedit.pos.start.base
                #                 end = hgvs_genomic.posedit.pos.end.base
                #                 ref_bases = validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, end)
                #                 lhb = validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                #                 rhb = validator.sf.fetch_seq(str(hgvs_genomic.ac), start, end)
                #                 hgvs_genomic.posedit.edit.ref = lhb + rhb
                #                 hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                #                 reverse_normalized_hgvs_genomic = variant.reverse_normalizer.normalize(
                #                     hgvs_genomic)
                #
                #     hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)
                #     # Store a copy for later use
                #     stored_hgvs_genomic_5pr = copy.deepcopy(hgvs_genomic_5pr)
                #
                #     # Make VCF
                #     vcf_dict = vvHGVS.hgvs2vcf(reverse_normalized_hgvs_genomic, variant.primary_assembly,
                #                                variant.reverse_normalizer, validator.sf)
                #     chr = vcf_dict['chr']
                #     pos = vcf_dict['pos']
                #     ref = vcf_dict['ref']
                #     alt = vcf_dict['alt']
                #
                #     # Look for exonic gaps within transcript or chromosome
                #     no_normalized_c = 'false'  # Mark true to not produce an additional normalization of c.
                #
                #     # Generate an end position
                #     end = str(int(pos) + len(ref) - 1)
                #     pos = str(pos)
                #
                #     # Store a not real deletion insertion to test for gapping
                #     stored_hgvs_not_delins = validator.hp.parse_hgvs_variant(str(
                #         hgvs_genomic_5pr.ac) + ':' + hgvs_genomic_5pr.type + '.' + pos + '_' + end + 'del' + ref + 'ins' + alt)
                #     v = [chr, pos, ref, alt]
                #
                #     # Save a copy of current hgvs_coding
                #     try:
                #         saved_hgvs_coding = variant.no_norm_evm.g_to_t(stored_hgvs_not_delins,
                #                                                hgvs_coding.ac)
                #     except Exception as e:
                #         if str(
                #                 e) == 'start or end or both are beyond the bounds of transcript record':
                #             saved_hgvs_coding = hgvs_coding
                #             continue
                #
                #     # Detect intronic variation using normalization
                #     intronic_variant = 'false'
                #     # Look for normalized variant options that do not match hgvs_coding
                #     if orientation == -1:
                #         # position genomic at its most 5 prime position
                #         try:
                #             query_genomic = variant.reverse_normalizer.normalize(hgvs_genomic)
                #         except:
                #             query_genomic = hgvs_genomic
                #         # Map to the transcript ant test for movement
                #         try:
                #             hgvs_seek_var = variant.evm.g_to_t(query_genomic, hgvs_coding.ac)
                #         except hgvs.exceptions.HGVSError as e:
                #             hgvs_seek_var = saved_hgvs_coding
                #         else:
                #             seek_var = fn.valstr(hgvs_seek_var)
                #             seek_ac = str(hgvs_seek_var.ac)
                #         if (
                #                 hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                #                 hgvs_coding.posedit.pos.start.base + hgvs_coding.posedit.pos.start.offset) and (
                #                 hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                #                 hgvs_coding.posedit.pos.end.base + hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
                #             pass
                #         else:
                #             hgvs_seek_var = saved_hgvs_coding
                #
                #     elif orientation != -1:
                #         # position genomic at its most 3 prime position
                #         try:
                #             query_genomic = variant.hn.normalize(hgvs_genomic)
                #         except:
                #             query_genomic = hgvs_genomic
                #         # Map to the transcript and test for movement
                #         try:
                #             hgvs_seek_var = variant.evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
                #         except hgvs.exceptions.HGVSError as e:
                #             hgvs_seek_var = saved_hgvs_coding
                #         seek_var = fn.valstr(hgvs_seek_var)
                #         seek_ac = str(hgvs_seek_var.ac)
                #         if (
                #                 hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                #                 hgvs_coding.posedit.pos.start.base + hgvs_coding.posedit.pos.start.offset) and (
                #                 hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                #                 hgvs_coding.posedit.pos.end.base + hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
                #             pass
                #         else:
                #             hgvs_seek_var = saved_hgvs_coding
                #
                #     try:
                #         intron_test = variant.hn.normalize(hgvs_seek_var)
                #     except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                #         error = str(e)
                #         if re.match('Normalization of intronic variants is not supported',
                #                     error) or re.match(
                #             'Unsupported normalization of variants spanning the exon-intron boundary',
                #             error):
                #             if re.match(
                #                     'Unsupported normalization of variants spanning the exon-intron boundary',
                #                     error):
                #                 intronic_variant = 'hard_fail'
                #             else:
                #                 # Double check to see whether the variant is actually intronic?
                #                 for exon in ori:
                #                     genomic_start = int(exon['alt_start_i'])
                #                     genomic_end = int(exon['alt_end_i'])
                #                     if (
                #                             hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                #                             hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                #                         intronic_variant = 'false'
                #                         break
                #                     else:
                #                         intronic_variant = 'true'
                #
                #     if intronic_variant != 'hard_fail':
                #         if re.search(r'\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search(r'\d+\-',
                #                                                                             str(
                #                                                                                 hgvs_seek_var.posedit.pos)) or re.search(
                #             r'\*\d+\+', str(
                #                 hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+\-', str(
                #             hgvs_seek_var.posedit.pos)):
                #             # Double check to see whether the variant is actually intronic?
                #             for exon in ori:
                #                 genomic_start = int(exon['alt_start_i'])
                #                 genomic_end = int(exon['alt_end_i'])
                #                 if (
                #                         hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                #                         hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                #                     intronic_variant = 'false'
                #                     break
                #                 else:
                #                     intronic_variant = 'true'
                #
                #     if intronic_variant != 'true':
                #         # Flag RefSeqGene for ammendment
                #         # amend_RefSeqGene = 'false'
                #         # Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping transcript lengths
                #         if stored_hgvs_not_delins != '':
                #             # Refresh hgvs_not_delins from stored_hgvs_not_delins
                #             hgvs_not_delins = copy.deepcopy(stored_hgvs_not_delins)
                #             # This test will only occur in dup of single base, insertion or substitution
                #             if not re.search('_', str(hgvs_not_delins.posedit.pos)):
                #                 if re.search('dup',
                #                              hgvs_genomic_5pr.posedit.edit.type) or re.search(
                #                     'ins', hgvs_genomic_5pr.posedit.edit.type):
                #                     # For gap in chr, map to t. - but becaouse we have pushed to 5 prime by norm, add 1 to end pos
                #                     plussed_hgvs_not_delins = copy.deepcopy(hgvs_not_delins)
                #                     plussed_hgvs_not_delins.posedit.pos.end.base = plussed_hgvs_not_delins.posedit.pos.end.base + 1
                #                     plussed_hgvs_not_delins.posedit.edit.ref = ''
                #                     transcript_variant = variant.no_norm_evm.g_to_t(plussed_hgvs_not_delins,
                #                                                             str(
                #                                                                 saved_hgvs_coding.ac))
                #                     if ((
                #                             transcript_variant.posedit.pos.end.base - transcript_variant.posedit.pos.start.base) > (
                #                             hgvs_genomic_5pr.posedit.pos.end.base - hgvs_genomic_5pr.posedit.pos.start.base)):
                #                         if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                #                             hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                #                             start = hgvs_not_delins.posedit.pos.start.base - 1
                #                             end = hgvs_not_delins.posedit.pos.end.base
                #                             ref_bases = validator.sf.fetch_seq(str(hgvs_not_delins.ac), start,
                #                                                           end)
                #                             hgvs_not_delins.posedit.edit.ref = ref_bases
                #                             hgvs_not_delins.posedit.edit.alt = ref_bases[
                #                                                                :1] + hgvs_not_delins.posedit.edit.alt[
                #                                                                      1:] + ref_bases[
                #                                                                            1:]
                #                         elif re.search('ins', str(
                #                                 hgvs_genomic_5pr.posedit.edit)) and re.search('del',
                #                                                                               str(
                #                                                                                   hgvs_genomic_5pr.posedit.edit)):
                #                             hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                #                         elif re.search('ins', str(
                #                                 hgvs_genomic_5pr.posedit.edit)) and not re.search(
                #                             'del',
                #                             str(
                #                                 hgvs_genomic_5pr.posedit.edit)):
                #                             hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                #                             start = hgvs_not_delins.posedit.pos.start.base - 1
                #                             end = hgvs_not_delins.posedit.pos.end.base
                #                             ref_bases = validator.sf.fetch_seq(str(hgvs_not_delins.ac), start,
                #                                                           end)
                #                             hgvs_not_delins.posedit.edit.ref = ref_bases
                #                             hgvs_not_delins.posedit.edit.alt = ref_bases[
                #                                                                :1] + hgvs_not_delins.posedit.edit.alt[
                #                                                                      1:] + ref_bases[
                #                                                                            1:]
                #                     else:
                #                         if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                #                             hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                #                             start = hgvs_not_delins.posedit.pos.start.base - 1
                #                             end = hgvs_not_delins.posedit.pos.end.base
                #                             ref_bases = validator.sf.fetch_seq(str(hgvs_not_delins.ac), start,
                #                                                           end)
                #                             hgvs_not_delins.posedit.edit.ref = ref_bases
                #                             hgvs_not_delins.posedit.edit.alt = ref_bases[
                #                                                                :1] + hgvs_not_delins.posedit.edit.alt[
                #                                                                      1:] + ref_bases[
                #                                                                            1:]
                #                         elif re.search('ins', str(
                #                                 hgvs_genomic_5pr.posedit.edit)) and re.search('del',
                #                                                                               str(
                #                                                                                   hgvs_genomic_5pr.posedit.edit)):
                #                             hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                #                         elif re.search('ins', str(
                #                                 hgvs_genomic_5pr.posedit.edit)) and not re.search(
                #                             'del',
                #                             str(
                #                                 hgvs_genomic_5pr.posedit.edit)):
                #                             hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                #                             start = hgvs_not_delins.posedit.pos.start.base - 1
                #                             end = hgvs_not_delins.posedit.pos.end.base
                #                             ref_bases = validator.sf.fetch_seq(str(hgvs_not_delins.ac), start,
                #                                                           end)
                #                             hgvs_not_delins.posedit.edit.ref = ref_bases
                #                             hgvs_not_delins.posedit.edit.alt = ref_bases[
                #                                                                :1] + hgvs_not_delins.posedit.edit.alt[
                #                                                                      1:] + ref_bases[
                #                                                                            1:]
                #                 else:
                #                     pass
                #             else:
                #                 pass
                #             tx_hgvs_not_delins = variant.no_norm_evm.g_to_n(hgvs_not_delins,
                #                                                     saved_hgvs_coding.ac)
                #             # Create normalized version of tx_hgvs_not_delins
                #             rn_tx_hgvs_not_delins = copy.deepcopy(tx_hgvs_not_delins)
                #             # Check for +1 base and adjust
                #             if re.search(r'\+',
                #                          str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search(
                #                 r'\+',
                #                 str(
                #                     rn_tx_hgvs_not_delins.posedit.pos.start)):
                #                 # Remove offsetting to span the gap
                #                 rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                #                 rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                #                 rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
                #                 rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                #                 try:
                #                     rn_tx_hgvs_not_delins.posedit.edit.alt = ''
                #                 except:
                #                     fn.exceptPass()
                #
                #             elif re.search(r'\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                #                 # move tx end base to next available non-offset base
                #                 rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
                #                 rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                #                 rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                #                 if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                #                     test_tx_var = variant.no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                #                 else:
                #                     test_tx_var = rn_tx_hgvs_not_delins
                #                 # re-make genomic and tx
                #                 hgvs_not_delins = validator.myvm_t_to_g(test_tx_var, alt_chr,
                #                                                    variant.no_norm_evm, variant.hn)
                #                 rn_tx_hgvs_not_delins = variant.no_norm_evm.g_to_n(hgvs_not_delins,
                #                                                            str(
                #                                                                saved_hgvs_coding.ac))
                #             elif re.search(r'\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                #                 # move tx start base to previous available non-offset base
                #                 rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                #                 rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                #                 if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                #                     test_tx_var = variant.no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                #                 else:
                #                     test_tx_var = rn_tx_hgvs_not_delins
                #                 # re-make genomic and tx
                #                 hgvs_not_delins = validator.myvm_t_to_g(test_tx_var, alt_chr,
                #                                                    variant.no_norm_evm, variant.hn)
                #                 rn_tx_hgvs_not_delins = variant.no_norm_evm.g_to_n(hgvs_not_delins,
                #                                                            str(
                #                                                                saved_hgvs_coding.ac))
                #                 rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                #             #                                                 else:
                #             #                                                     pass
                #
                #             # Check for -ve base and adjust
                #             elif re.search(r'\-',
                #                            str(
                #                                rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search(
                #                 r'\-',
                #                 str(
                #                     rn_tx_hgvs_not_delins.posedit.pos.start)):
                #                 # Remove offsetting to span the gap
                #                 rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                #                 rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                #                 rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
                #                 rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                #                 try:
                #                     rn_tx_hgvs_not_delins.posedit.edit.alt = ''
                #                 except:
                #                     fn.exceptPass()
                #             elif re.search(r'\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                #                 # move tx end base back to next available non-offset base
                #                 rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                #                 # Delete the ref
                #                 rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                #                 # Add the additional base to the ALT
                #                 start = rn_tx_hgvs_not_delins.posedit.pos.end.base - 1
                #                 end = rn_tx_hgvs_not_delins.posedit.pos.end.base
                #                 ref_bases = validator.sf.fetch_seq(str(tx_hgvs_not_delins.ac), start, end)
                #                 rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + ref_bases
                #                 if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                #                     test_tx_var = variant.no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                #                 else:
                #                     test_tx_var = rn_tx_hgvs_not_delins
                #                 # re-make genomic and tx
                #                 hgvs_not_delins = validator.myvm_t_to_g(test_tx_var, alt_chr,
                #                                                    variant.no_norm_evm, variant.hn)
                #                 rn_tx_hgvs_not_delins = variant.no_norm_evm.g_to_n(hgvs_not_delins,
                #                                                            str(
                #                                                                saved_hgvs_coding.ac))
                #             elif re.search(r'\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                #                 # move tx start base to previous available non-offset base
                #                 rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                #                 rn_tx_hgvs_not_delins.posedit.pos.start.base = rn_tx_hgvs_not_delins.posedit.pos.start.base - 1
                #                 rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                #                 if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                #                     test_tx_var = variant.no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                #                 else:
                #                     test_tx_var = rn_tx_hgvs_not_delins
                #                 # re-make genomic and tx
                #                 hgvs_not_delins = validator.myvm_t_to_g(test_tx_var, alt_chr,
                #                                                    variant.no_norm_evm, variant.hn)
                #                 rn_tx_hgvs_not_delins = variant.no_norm_evm.g_to_n(hgvs_not_delins,
                #                                                            str(
                #                                                                saved_hgvs_coding.ac))
                #                 rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                #             else:
                #                 fn.exceptPass()
                #
                #             # Logic
                #             if len(hgvs_not_delins.posedit.edit.ref) < len(
                #                     rn_tx_hgvs_not_delins.posedit.edit.ref):
                #                 gap_length = len(rn_tx_hgvs_not_delins.posedit.edit.ref) - len(
                #                     hgvs_not_delins.posedit.edit.ref)
                #                 disparity_deletion_in = ['chromosome', gap_length]
                #             elif len(hgvs_not_delins.posedit.edit.ref) > len(
                #                     rn_tx_hgvs_not_delins.posedit.edit.ref):
                #                 gap_length = len(hgvs_not_delins.posedit.edit.ref) - len(
                #                     rn_tx_hgvs_not_delins.posedit.edit.ref)
                #                 disparity_deletion_in = ['transcript', gap_length]
                #             else:
                #                 re_capture_tx_variant = []
                #                 for possibility in hgvs_genomic_possibilities:
                #                     if possibility == '':
                #                         continue
                #                     hgvs_t_possibility = validator.vm.g_to_t(possibility, hgvs_coding.ac)
                #                     if hgvs_t_possibility.posedit.edit.type == 'ins':
                #                         try:
                #                             hgvs_t_possibility = validator.vm.c_to_n(hgvs_t_possibility)
                #                         except:
                #                             continue
                #                         if hgvs_t_possibility.posedit.pos.start.offset != 0 or hgvs_t_possibility.posedit.pos.end.offset != 0:
                #                             continue
                #                         ins_ref = validator.sf.fetch_seq(hgvs_t_possibility.ac,
                #                                                     hgvs_t_possibility.posedit.pos.start.base - 1,
                #                                                     hgvs_t_possibility.posedit.pos.start.base + 1)
                #                         try:
                #                             hgvs_t_possibility = validator.vm.n_to_c(hgvs_t_possibility)
                #                         except:
                #                             continue
                #                         hgvs_t_possibility.posedit.edit.ref = ins_ref
                #                         hgvs_t_possibility.posedit.edit.alt = ins_ref[
                #                                                                   0] + hgvs_t_possibility.posedit.edit.alt + \
                #                                                               ins_ref[1]
                #                     if possibility.posedit.edit.type == 'ins':
                #                         ins_ref = validator.sf.fetch_seq(possibility.ac,
                #                                                     possibility.posedit.pos.start.base - 1,
                #                                                     possibility.posedit.pos.end.base)
                #                         possibility.posedit.edit.ref = ins_ref
                #                         possibility.posedit.edit.alt = ins_ref[
                #                                                            0] + possibility.posedit.edit.alt + \
                #                                                        ins_ref[1]
                #                     if len(hgvs_t_possibility.posedit.edit.ref) < len(
                #                             possibility.posedit.edit.ref):
                #                         gap_length = len(possibility.posedit.edit.ref) - len(
                #                             hgvs_t_possibility.posedit.edit.ref)
                #                         re_capture_tx_variant = ['transcript', gap_length,
                #                                                  hgvs_t_possibility]
                #                         hgvs_not_delins = possibility
                #                         hgvs_genomic_5pr = possibility
                #                         break
                #
                #                 if re_capture_tx_variant != []:
                #                     try:
                #                         tx_hgvs_not_delins = validator.vm.c_to_n(re_capture_tx_variant[2])
                #                     except:
                #                         tx_hgvs_not_delins = re_capture_tx_variant[2]
                #                     disparity_deletion_in = re_capture_tx_variant[0:-1]
                #                 else:
                #                     pass
                #
                #         # Final sanity checks
                #         try:
                #             validator.vm.g_to_t(hgvs_not_delins, tx_hgvs_not_delins.ac)
                #         except Exception as e:
                #             if str(
                #                     e) == 'start or end or both are beyond the bounds of transcript record':
                #                 continue
                #         try:
                #             variant.hn.normalize(tx_hgvs_not_delins)
                #         except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                #             error = str(e)
                #
                #             if re.match('Normalization of intronic variants is not supported',
                #                         error) or re.match(
                #                 'Unsupported normalization of variants spanning the exon-intron boundary',
                #                 error):
                #                 if re.match(
                #                         'Unsupported normalization of variants spanning the exon-intron boundary',
                #                         error):
                #                     continue
                #                 elif re.match('Normalization of intronic variants is not supported',
                #                               error):
                #                     # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                #                     disparity_deletion_in = ['transcript', 'Requires Analysis']
                #
                #         # Recreate hgvs_genomic
                #         if disparity_deletion_in[0] == 'transcript':
                #             hgvs_genomic = hgvs_not_delins
                #
                #         # Find oddly placed gaps where the tx variant is encompassed in the gap
                #         if disparity_deletion_in[0] == 'false' and (
                #                 possibility_counter == 3 or possibility_counter == 4):
                #             rg = variant.reverse_normalizer.normalize(hgvs_not_delins)
                #             rtx = validator.vm.g_to_t(rg, tx_hgvs_not_delins.ac)
                #             fg = variant.hn.normalize(hgvs_not_delins)
                #             ftx = validator.vm.g_to_t(fg, tx_hgvs_not_delins.ac)
                #             if (
                #                     rtx.posedit.pos.start.offset == 0 and rtx.posedit.pos.end.offset == 0) and (
                #                     ftx.posedit.pos.start.offset != 0 and ftx.posedit.pos.end.offset != 0):
                #                 exons = validator.hdp.get_tx_exons(ftx.ac, hgvs_not_delins.ac, validator.alt_aln_method)
                #                 exonic = False
                #                 for ex_test in exons:
                #                     if ftx.posedit.pos.start.base in range(ex_test[6], ex_test[
                #                         7]) and ftx.posedit.pos.end.base in range(ex_test[6],
                #                                                                   ex_test[7]):
                #                         exonic = True
                #                 if exonic is True:
                #                     hgvs_not_delins = fg
                #                     hgvs_genomic = fg
                #                     hgvs_genomic_5pr = fg
                #                     try:
                #                         tx_hgvs_not_delins = validator.vm.c_to_n(ftx)
                #                     except Exception:
                #                         tx_hgvs_not_delins = ftx
                #                     disparity_deletion_in = ['transcript', 'Requires Analysis']
                #
                #         # Pre-processing of tx_hgvs_not_delins
                #         try:
                #             if tx_hgvs_not_delins.posedit.edit.alt is None:
                #                 tx_hgvs_not_delins.posedit.edit.alt = ''
                #         except Exception as e:
                #             if str(e) == "'Dup' object has no attribute 'alt'":
                #                 tx_hgvs_not_delins_delins_from_dup = tx_hgvs_not_delins.ac + ':' + tx_hgvs_not_delins.type + '.' + str(
                #                     tx_hgvs_not_delins.posedit.pos.start) + '_' + str(
                #                     tx_hgvs_not_delins.posedit.pos.end) + 'del' + tx_hgvs_not_delins.posedit.edit.ref + 'ins' + tx_hgvs_not_delins.posedit.edit.ref + tx_hgvs_not_delins.posedit.edit.ref
                #                 tx_hgvs_not_delins = validator.hp.parse_hgvs_variant(
                #                     tx_hgvs_not_delins_delins_from_dup)
                #
                #         if disparity_deletion_in[0] == 'transcript':
                #             # amend_RefSeqGene = 'true'
                #             # ANY VARIANT WHOLLY WITHIN THE GAP
                #             if (re.search(r'\+',
                #                           str(tx_hgvs_not_delins.posedit.pos.start)) or re.search(
                #                 r'\-', str(tx_hgvs_not_delins.posedit.pos.start))) and (
                #                     re.search(r'\+',
                #                               str(tx_hgvs_not_delins.posedit.pos.end)) or re.search(
                #                 r'\-', str(tx_hgvs_not_delins.posedit.pos.end))):
                #                 gapped_transcripts = gapped_transcripts + ' ' + str(
                #                     tx_hgvs_not_delins.ac)
                #
                #                 # Copy the current variant
                #                 tx_gap_fill_variant = copy.deepcopy(tx_hgvs_not_delins)
                #                 try:
                #                     if tx_gap_fill_variant.posedit.edit.alt is None:
                #                         tx_gap_fill_variant.posedit.edit.alt = ''
                #                 except Exception as e:
                #                     if str(e) == "'Dup' object has no attribute 'alt'":
                #                         tx_gap_fill_variant_delins_from_dup = tx_gap_fill_variant.ac + ':' + tx_gap_fill_variant.type + '.' + str(
                #                             tx_gap_fill_variant.posedit.pos.start) + '_' + str(
                #                             tx_gap_fill_variant.posedit.pos.end) + 'del' + tx_gap_fill_variant.posedit.edit.ref + 'ins' + tx_gap_fill_variant.posedit.edit.ref + tx_gap_fill_variant.posedit.edit.ref
                #                         tx_gap_fill_variant = validator.hp.parse_hgvs_variant(
                #                             tx_gap_fill_variant_delins_from_dup)
                #
                #                 # Identify which half of the NOT-intron the start position of the variant is in
                #                 if re.search(r'\-', str(tx_gap_fill_variant.posedit.pos.start)):
                #                     tx_gap_fill_variant.posedit.pos.start.base = tx_gap_fill_variant.posedit.pos.start.base - 1
                #                     tx_gap_fill_variant.posedit.pos.start.offset = int(
                #                         '0')  # int('+1')
                #                     tx_gap_fill_variant.posedit.pos.end.offset = int(
                #                         '0')  # int('-1')
                #                     tx_gap_fill_variant.posedit.edit.alt = ''
                #                     tx_gap_fill_variant.posedit.edit.ref = ''
                #                 elif re.search(r'\+', str(tx_gap_fill_variant.posedit.pos.start)):
                #                     tx_gap_fill_variant.posedit.pos.start.offset = int(
                #                         '0')  # int('+1')
                #                     tx_gap_fill_variant.posedit.pos.end.base = tx_gap_fill_variant.posedit.pos.end.base + 1
                #                     tx_gap_fill_variant.posedit.pos.end.offset = int(
                #                         '0')  # int('-1')
                #                     tx_gap_fill_variant.posedit.edit.alt = ''
                #                     tx_gap_fill_variant.posedit.edit.ref = ''
                #
                #                 try:
                #                     tx_gap_fill_variant = validator.vm.n_to_c(tx_gap_fill_variant)
                #                 except:
                #                     fn.exceptPass()
                #                 genomic_gap_fill_variant = validator.vm.t_to_g(tx_gap_fill_variant,
                #                                                           reverse_normalized_hgvs_genomic.ac)
                #                 genomic_gap_fill_variant.posedit.edit.alt = genomic_gap_fill_variant.posedit.edit.ref
                #
                #                 try:
                #                     c_tx_hgvs_not_delins = validator.vm.n_to_c(tx_hgvs_not_delins)
                #                 except Exception:
                #                     c_tx_hgvs_not_delins = copy.copy(tx_hgvs_not_delins)
                #                 genomic_gap_fill_variant_alt = validator.vm.t_to_g(c_tx_hgvs_not_delins,
                #                                                               hgvs_genomic_5pr.ac)
                #
                #                 # Ensure an ALT exists
                #                 try:
                #                     if genomic_gap_fill_variant_alt.posedit.edit.alt is None:
                #                         genomic_gap_fill_variant_alt.posedit.edit.alt = 'X'
                #                 except Exception as e:
                #                     if str(e) == "'Dup' object has no attribute 'alt'":
                #                         genomic_gap_fill_variant_delins_from_dup = genomic_gap_fill_variant.ac + ':' + genomic_gap_fill_variant.type + '.' + str(
                #                             genomic_gap_fill_variant.posedit.pos.start.base) + '_' + str(
                #                             genomic_gap_fill_variant.posedit.pos.end.base) + 'del' + genomic_gap_fill_variant.posedit.edit.ref + 'ins' + genomic_gap_fill_variant.posedit.edit.ref + genomic_gap_fill_variant.posedit.edit.ref
                #                         genomic_gap_fill_variant = validator.hp.parse_hgvs_variant(
                #                             genomic_gap_fill_variant_delins_from_dup)
                #                         genomic_gap_fill_variant_alt_delins_from_dup = genomic_gap_fill_variant_alt.ac + ':' + genomic_gap_fill_variant_alt.type + '.' + str(
                #                             genomic_gap_fill_variant_alt.posedit.pos.start.base) + '_' + str(
                #                             genomic_gap_fill_variant_alt.posedit.pos.end.base) + 'del' + genomic_gap_fill_variant_alt.posedit.edit.ref + 'ins' + genomic_gap_fill_variant_alt.posedit.edit.ref + genomic_gap_fill_variant_alt.posedit.edit.ref
                #                         genomic_gap_fill_variant_alt = validator.hp.parse_hgvs_variant(
                #                             genomic_gap_fill_variant_alt_delins_from_dup)
                #
                #                 # Correct insertion alts
                #                 if genomic_gap_fill_variant_alt.posedit.edit.type == 'ins':
                #                     append_ref = validator.sf.fetch_seq(genomic_gap_fill_variant_alt.ac,
                #                                                    genomic_gap_fill_variant_alt.posedit.pos.start.base - 1,
                #                                                    genomic_gap_fill_variant_alt.posedit.pos.end.base)
                #                     genomic_gap_fill_variant_alt.posedit.edit.alt = append_ref[
                #                                                                         0] + genomic_gap_fill_variant_alt.posedit.edit.alt + \
                #                                                                     append_ref[1]
                #
                #                 # Split the reference and replacing alt sequence into a dictionary
                #                 reference_bases = list(genomic_gap_fill_variant.posedit.edit.ref)
                #                 if genomic_gap_fill_variant_alt.posedit.edit.alt is not None:
                #                     alternate_bases = list(
                #                         genomic_gap_fill_variant_alt.posedit.edit.alt)
                #                 else:
                #                     # Deletions with no ins
                #                     pre_alternate_bases = list(
                #                         genomic_gap_fill_variant_alt.posedit.edit.ref)
                #                     alternate_bases = []
                #                     for base in pre_alternate_bases:
                #                         alternate_bases.append('X')
                #
                #                 # Create the dictionaries
                #                 ref_start = genomic_gap_fill_variant.posedit.pos.start.base
                #                 alt_start = genomic_gap_fill_variant_alt.posedit.pos.start.base
                #                 ref_base_dict = {}
                #                 for base in reference_bases:
                #                     ref_base_dict[ref_start] = str(base)
                #                     ref_start = ref_start + 1
                #
                #                 alt_base_dict = {}
                #
                #                 # NEED TO SEARCH FOR RANGE = and replace with interval_range
                #                 # Need to search for int and replace with integer
                #
                #                 # Note, all variants will be forced into the format delete insert
                #                 # Deleted bases in the ALT will be substituted for X
                #                 for integer in range(
                #                         genomic_gap_fill_variant_alt.posedit.pos.start.base,
                #                         genomic_gap_fill_variant_alt.posedit.pos.end.base + 1, 1):
                #                     if integer == alt_start:
                #                         alt_base_dict[integer] = str(''.join(alternate_bases))
                #                     else:
                #                         alt_base_dict[integer] = 'X'
                #
                #                 # Generate the alt sequence
                #                 alternate_sequence_bases = []
                #                 for integer in range(
                #                         genomic_gap_fill_variant.posedit.pos.start.base,
                #                         genomic_gap_fill_variant.posedit.pos.end.base + 1,
                #                         1):
                #                     if integer in list(alt_base_dict.keys()):
                #                         alternate_sequence_bases.append(alt_base_dict[integer])
                #                     else:
                #                         alternate_sequence_bases.append(ref_base_dict[integer])
                #                 alternate_sequence = ''.join(alternate_sequence_bases)
                #                 alternate_sequence = alternate_sequence.replace('X', '')
                #
                #                 # Add the new alt to the gap fill variant and generate transcript variant
                #                 genomic_gap_fill_variant.posedit.edit.alt = alternate_sequence
                #                 hgvs_refreshed_variant = validator.vm.g_to_t(genomic_gap_fill_variant,
                #                                                         tx_gap_fill_variant.ac)
                #
                #                 # Set warning
                #                 gap_size = str(len(genomic_gap_fill_variant.posedit.edit.ref) - 2)
                #                 disparity_deletion_in[1] = [gap_size]
                #                 auto_info = auto_info + str(
                #                     stored_hgvs_not_delins.ac) + ':g.' + str(
                #                     stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + gap_size + ' genomic base(s) that fail to align to transcript ' + str(
                #                     tx_hgvs_not_delins.ac)
                #                 non_valid_caution = 'true'
                #
                #                 # Alignment position
                #                 for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                #                 if re.match('NM_', str(for_location_c)):
                #                     for_location_c = variant.no_norm_evm.n_to_c(tx_hgvs_not_delins)
                #                 if re.match(r'\-', str(for_location_c.posedit.pos.start.offset)):
                #                     gps = for_location_c.posedit.pos.start.base - 1
                #                     gpe = for_location_c.posedit.pos.start.base
                #                 else:
                #                     gps = for_location_c.posedit.pos.start.base
                #                     gpe = for_location_c.posedit.pos.start.base + 1
                #                 gap_position = ' between positions c.' + str(gps) + '_' + str(
                #                     gpe) + '\n'
                #                 auto_info = auto_info + '%s' % (gap_position)
                #
                #             else:
                #                 if tx_hgvs_not_delins.posedit.pos.start.offset == 0 and tx_hgvs_not_delins.posedit.pos.end.offset == 0:
                #                     # In this instance, we have identified a transcript gap but the n. version of
                #                     # the transcript variant but do not have a position which actually hits the gap,
                #                     # so the variant likely spans the gap, and is not picked up by an offset.
                #                     try:
                #                         c1 = validator.vm.n_to_c(tx_hgvs_not_delins)
                #                     except:
                #                         c1 = tx_hgvs_not_delins
                #                     g1 = validator.nr_vm.t_to_g(c1, hgvs_genomic.ac)
                #                     g3 = validator.nr_vm.t_to_g(c1, hgvs_genomic.ac)
                #                     g2 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
                #                     ng2 = variant.hn.normalize(g2)
                #                     g3.posedit.pos.end.base = g3.posedit.pos.start.base + (
                #                             len(g3.posedit.edit.ref) - 1)
                #                     try:
                #                         c2 = validator.vm.g_to_t(g3, c1.ac)
                #                         if c2.posedit.pos.start.offset == 0 and c2.posedit.pos.end.offset == 0:
                #                             pass
                #                         else:
                #                             tx_hgvs_not_delins = c2
                #                             try:
                #                                 tx_hgvs_not_delins = validator.vm.c_to_n(tx_hgvs_not_delins)
                #                             except hgvs.exceptions.HGVSError:
                #                                 fn.exceptPass()
                #                     except hgvs.exceptions.HGVSInvalidVariantError:
                #                         fn.exceptPass()
                #
                #                 if re.search(r'\+', str(
                #                         tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                #                     r'\+',
                #                     str(
                #                         tx_hgvs_not_delins.posedit.pos.end)):
                #                     auto_info = auto_info + str(
                #                         stored_hgvs_not_delins.ac) + ':g.' + str(
                #                         stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                #                         disparity_deletion_in[
                #                             1]) + ' genomic base(s) that fail to align to transcript ' + str(
                #                         tx_hgvs_not_delins.ac)
                #                     non_valid_caution = 'true'
                #                     try:
                #                         c2 = validator.vm.n_to_c(tx_hgvs_not_delins)
                #                     except:
                #                         c2 = tx_hgvs_not_delins
                #                     c1 = copy.deepcopy(c2)
                #                     c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                #                     c1.posedit.pos.start.offset = 0
                #                     c1.posedit.pos.end = c2.posedit.pos.start
                #                     c1.posedit.edit.ref = ''
                #                     c1.posedit.edit.alt = ''
                #                     if orientation != -1:
                #                         g1 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
                #                         g2 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                #                         g1.posedit.edit.alt = g1.posedit.edit.ref
                #                     else:
                #                         g1 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                #                         g2 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
                #                         g2.posedit.edit.alt = g2.posedit.edit.ref
                #                     reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                #                     alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                #                     g3 = copy.deepcopy(g1)
                #                     g3.posedit.pos.end.base = g2.posedit.pos.end.base
                #                     g3.posedit.edit.ref = reference
                #                     g3.posedit.edit.alt = alternate
                #                     c3 = validator.vm.g_to_t(g3, c1.ac)
                #                     hgvs_refreshed_variant = c3
                #                     # Alignment position
                #                     for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                #                     if re.match('NM_', str(for_location_c)):
                #                         for_location_c = variant.no_norm_evm.n_to_c(tx_hgvs_not_delins)
                #                         gps = for_location_c.posedit.pos.start.base
                #                         gpe = for_location_c.posedit.pos.start.base + 1
                #                     gap_position = ' between positions c.' + str(gps) + '_' + str(
                #                         gpe) + '\n'
                #                     # Warn update
                #                     auto_info = auto_info + '%s' % (gap_position)
                #                 elif re.search(r'\+', str(
                #                         tx_hgvs_not_delins.posedit.pos.end)) and not re.search(r'\+',
                #                                                                                str(
                #                                                                                    tx_hgvs_not_delins.posedit.pos.start)):
                #                     auto_info = auto_info + 'Genome position ' + str(
                #                         stored_hgvs_not_delins.ac) + ':g.' + str(
                #                         stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(
                #                         disparity_deletion_in[1]) + '-bp gap in transcript ' + str(
                #                         tx_hgvs_not_delins.ac)
                #                     gapped_transcripts = gapped_transcripts + ' ' + str(
                #                         tx_hgvs_not_delins.ac)
                #                     non_valid_caution = 'true'
                #                     try:
                #                         c1 = validator.vm.n_to_c(tx_hgvs_not_delins)
                #                     except:
                #                         c1 = tx_hgvs_not_delins
                #                     c2 = copy.deepcopy(c1)
                #                     c2.posedit.pos.start = c1.posedit.pos.end
                #                     c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                #                     c2.posedit.pos.end.offset = 0
                #                     c2.posedit.edit.ref = ''
                #                     c2.posedit.edit.alt = ''
                #                     if orientation != -1:
                #                         g1 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
                #                         g2 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                #                         g2.posedit.edit.alt = g2.posedit.edit.ref
                #                     else:
                #                         g1 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                #                         g2 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
                #                         g1.posedit.edit.alt = g1.posedit.edit.ref
                #                     reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                #                     alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                #                     g3 = copy.deepcopy(g1)
                #                     g3.posedit.pos.end.base = g2.posedit.pos.end.base
                #                     g3.posedit.edit.ref = reference
                #                     g3.posedit.edit.alt = alternate
                #                     c3 = validator.vm.g_to_t(g3, c1.ac)
                #                     hgvs_refreshed_variant = c3
                #                     # Alignment position
                #                     for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                #                     if re.match('NM_', str(for_location_c)):
                #                         for_location_c = variant.no_norm_evm.n_to_c(tx_hgvs_not_delins)
                #                     gps = for_location_c.posedit.pos.end.base
                #                     gpe = for_location_c.posedit.pos.end.base + 1
                #                     gap_position = ' between positions c.' + str(gps) + '_' + str(
                #                         gpe) + '\n'
                #                     # Warn update
                #                     auto_info = auto_info + '%s' % (gap_position)
                #                 elif re.search(r'\-', str(
                #                         tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                #                     r'\-',
                #                     str(
                #                         tx_hgvs_not_delins.posedit.pos.end)):
                #                     auto_info = auto_info + str(
                #                         stored_hgvs_not_delins.ac) + ':g.' + str(
                #                         stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                #                         disparity_deletion_in[
                #                             1]) + ' genomic base(s) that fail to align to transcript ' + str(
                #                         tx_hgvs_not_delins.ac)
                #                     non_valid_caution = 'true'
                #                     try:
                #                         c2 = validator.vm.n_to_c(tx_hgvs_not_delins)
                #                     except:
                #                         c2 = tx_hgvs_not_delins
                #                     c1 = copy.deepcopy(c2)
                #                     c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                #                     c1.posedit.pos.start.offset = 0
                #                     c1.posedit.pos.end = c2.posedit.pos.start
                #                     c1.posedit.edit.ref = ''
                #                     c1.posedit.edit.alt = ''
                #                     if orientation != -1:
                #                         g1 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
                #                         g2 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                #                         g1.posedit.edit.alt = g1.posedit.edit.ref
                #                     else:
                #                         g1 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                #                         g2 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
                #                         g2.posedit.edit.alt = g2.posedit.edit.ref
                #                     reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                #                     alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                #                     g3 = copy.deepcopy(g1)
                #                     g3.posedit.pos.end.base = g2.posedit.pos.end.base
                #                     g3.posedit.edit.ref = reference
                #                     g3.posedit.edit.alt = alternate
                #                     c3 = validator.vm.g_to_t(g3, c1.ac)
                #                     hgvs_refreshed_variant = c3
                #                     # Alignment position
                #                     for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                #                     if re.match('NM_', str(for_location_c)):
                #                         for_location_c = variant.no_norm_evm.n_to_c(tx_hgvs_not_delins)
                #                     gps = for_location_c.posedit.pos.start.base - 1
                #                     gpe = for_location_c.posedit.pos.start.base
                #                     gap_position = ' between positions c.' + str(gps) + '_' + str(
                #                         gpe) + '\n'
                #                     # Warn update
                #                     auto_info = auto_info + '%s' % (gap_position)
                #                 elif re.search(r'\-', str(
                #                         tx_hgvs_not_delins.posedit.pos.end)) and not re.search(r'\-',
                #                                                                                str(
                #                                                                                    tx_hgvs_not_delins.posedit.pos.start)):
                #                     auto_info = auto_info + 'Genome position ' + str(
                #                         stored_hgvs_not_delins.ac) + ':g.' + str(
                #                         stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(
                #                         disparity_deletion_in[1]) + '-bp gap in transcript ' + str(
                #                         tx_hgvs_not_delins.ac)
                #                     gapped_transcripts = gapped_transcripts + ' ' + str(
                #                         tx_hgvs_not_delins.ac)
                #                     non_valid_caution = 'true'
                #                     try:
                #                         c1 = validator.vm.n_to_c(tx_hgvs_not_delins)
                #                     except:
                #                         c1 = tx_hgvs_not_delins
                #                     c2 = copy.deepcopy(c1)
                #                     c2.posedit.pos.start = c1.posedit.pos.end
                #                     c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                #                     c2.posedit.pos.end.offset = 0
                #                     c2.posedit.edit.ref = ''
                #                     c2.posedit.edit.alt = ''
                #                     if orientation != -1:
                #                         g1 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
                #                         g2 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                #                         g2.posedit.edit.alt = g2.posedit.edit.ref
                #                     else:
                #                         g1 = validator.vm.t_to_g(c2, hgvs_genomic.ac)
                #                         g2 = validator.vm.t_to_g(c1, hgvs_genomic.ac)
                #                         g1.posedit.edit.alt = g1.posedit.edit.ref
                #                     reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                #                     alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                #                     g3 = copy.deepcopy(g1)
                #                     g3.posedit.pos.end.base = g2.posedit.pos.end.base
                #                     g3.posedit.edit.ref = reference
                #                     g3.posedit.edit.alt = alternate
                #                     c3 = validator.vm.g_to_t(g3, c1.ac)
                #                     hgvs_refreshed_variant = c3
                #                     # Alignment position
                #                     for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                #                     if re.match('NM_', str(for_location_c)):
                #                         for_location_c = variant.no_norm_evm.n_to_c(tx_hgvs_not_delins)
                #                     gps = for_location_c.posedit.pos.end.base - 1
                #                     gpe = for_location_c.posedit.pos.end.base
                #                     gap_position = ' between positions c.' + str(gps) + '_' + str(
                #                         gpe) + '\n'
                #                     # Warn update
                #                     auto_info = auto_info + '%s' % (gap_position)
                #                 else:
                #                     auto_info = auto_info + str(
                #                         stored_hgvs_not_delins.ac) + ':g.' + str(
                #                         stored_hgvs_not_delins.posedit.pos) + ' contains ' + str(
                #                         disparity_deletion_in[
                #                             1]) + ' genomic base(s) that fail to align to transcript ' + str(
                #                         tx_hgvs_not_delins.ac) + '\n'
                #                     hgvs_refreshed_variant = tx_hgvs_not_delins
                #
                #         # GAP IN THE CHROMOSOME
                #         elif disparity_deletion_in[0] == 'chromosome':
                #             # amend_RefSeqGene = 'true'
                #             if possibility_counter == 3:
                #                 hgvs_refreshed_variant = stash_tx_right
                #             elif possibility_counter == 4:
                #                 hgvs_refreshed_variant = stash_tx_left
                #             else:
                #                 hgvs_refreshed_variant = chromosome_normalized_hgvs_coding
                #             # Warn
                #             auto_info = auto_info + str(hgvs_refreshed_variant.ac) + ':c.' + str(
                #                 hgvs_refreshed_variant.posedit.pos) + ' contains ' + str(
                #                 disparity_deletion_in[
                #                     1]) + ' transcript base(s) that fail to align to chromosome ' + str(
                #                 hgvs_genomic.ac) + '\n'
                #         else:
                #             # Keep the same by re-setting rel_var
                #             hgvs_refreshed_variant = hgvs_coding
                #         # amend_RefSeqGene = 'false'
                #
                #         # Edit the output
                #         if re.match('NM_', str(hgvs_refreshed_variant.ac)) and not re.search('c',
                #                                                                              str(
                #                                                                                  hgvs_refreshed_variant.type)):
                #             hgvs_refreshed_variant = variant.no_norm_evm.n_to_c(hgvs_refreshed_variant)
                #         else:
                #             pass
                #
                #         try:
                #             variant.hn.normalize(hgvs_refreshed_variant)
                #         except Exception as e:
                #             error = str(e)
                #             # Ensure the final variant is not intronic nor does it cross exon boundaries
                #             if re.match('Normalization of intronic variants is not supported',
                #                         error) or re.match(
                #                 'Unsupported normalization of variants spanning the exon-intron boundary',
                #                 error):
                #                 hgvs_refreshed_variant = saved_hgvs_coding
                #             else:
                #                 continue
                #
                #         # Quick check to make sure the coding variant has not changed
                #         try:
                #             to_test = variant.hn.normalize(hgvs_refreshed_variant)
                #         except:
                #             to_test = hgvs_refreshed_variant
                #         if str(to_test.posedit.edit) != str(hgvs_coding.posedit.edit):
                #             # Try the next available genomic option
                #             if hgvs_coding.posedit.edit.type == 'identity' and to_test.posedit.edit.type == 'identity':
                #                 hgvs_coding = to_test
                #             else:
                #                 continue
                #
                #         # Update hgvs_genomic
                #         hgvs_alt_genomic = validator.myvm_t_to_g(hgvs_refreshed_variant, alt_chr,
                #                                             variant.no_norm_evm, variant.hn)
                #         if hgvs_alt_genomic.posedit.edit.type == 'identity':
                #             re_c = validator.vm.g_to_t(hgvs_alt_genomic, hgvs_refreshed_variant.ac)
                #             if (variant.hn.normalize(re_c)) != (variant.hn.normalize(hgvs_refreshed_variant)):
                #                 shuffle_left_g = copy.copy(hgvs_alt_genomic)
                #                 shuffle_left_g.posedit.edit.ref = ''
                #                 shuffle_left_g.posedit.edit.alt = ''
                #                 shuffle_left_g.posedit.pos.start.base = shuffle_left_g.posedit.pos.start.base - 1
                #                 shuffle_left_g.posedit.pos.end.base = shuffle_left_g.posedit.pos.end.base - 1
                #                 shuffle_left_g = variant.reverse_normalizer.normalize(shuffle_left_g)
                #                 re_c = validator.vm.g_to_t(shuffle_left_g, hgvs_refreshed_variant.ac)
                #                 if (variant.hn.normalize(re_c)) != (variant.hn.normalize(hgvs_refreshed_variant)):
                #                     hgvs_alt_genomic = shuffle_left_g
                #
                #                     # If it is intronic, these vairables will not have been set
                #     else:
                #         # amend_RefSeqGene = 'false'
                #         no_normalized_c = 'false'
                #
                #     # Break if gap has been detected
                #     if disparity_deletion_in[0] != 'false':
                #         break
                #
                # # Normailse hgvs_genomic
                # try:
                #     hgvs_alt_genomic = variant.hn.normalize(hgvs_alt_genomic)
                # except hgvs.exceptions.HGVSError as e:
                #     # Strange error caused by gap in genomic
                #     error = str(e)
                #     if re.search('base start position must be <= end position', error) and \
                #             disparity_deletion_in[0] == 'chromosome':
                #         if hgvs_alt_genomic.posedit.edit.type == 'delins':
                #             start = hgvs_alt_genomic.posedit.pos.start.base
                #             end = hgvs_alt_genomic.posedit.pos.end.base
                #             lhb = validator.sf.fetch_seq(str(hgvs_alt_genomic.ac), end - 1, end)
                #             rhb = validator.sf.fetch_seq(str(hgvs_alt_genomic.ac), start - 1, start)
                #             hgvs_alt_genomic.posedit.edit.ref = lhb + rhb
                #             hgvs_alt_genomic.posedit.edit.alt = lhb + hgvs_alt_genomic.posedit.edit.alt + rhb
                #             hgvs_alt_genomic.posedit.pos.start.base = end
                #             hgvs_alt_genomic.posedit.pos.end.base = start
                #             hgvs_alt_genomic = variant.hn.normalize(hgvs_alt_genomic)
                #         if hgvs_alt_genomic.posedit.edit.type == 'del':
                #             start = hgvs_alt_genomic.posedit.pos.start.base
                #             end = hgvs_alt_genomic.posedit.pos.end.base
                #             lhb = validator.sf.fetch_seq(str(hgvs_alt_genomic.ac), end - 1, end)
                #             rhb = validator.sf.fetch_seq(str(hgvs_alt_genomic.ac), start - 1, start)
                #             hgvs_alt_genomic.posedit.edit.ref = lhb + rhb
                #             hgvs_alt_genomic.posedit.edit.alt = lhb + rhb
                #             hgvs_alt_genomic.posedit.pos.start.base = end
                #             hgvs_alt_genomic.posedit.pos.end.base = start
                #             hgvs_alt_genomic = variant.hn.normalize(hgvs_alt_genomic)

                # Refresh the :g. variant
                multi_g.append(hgvs_alt_genomic)
            else:
                multi_g.append(hgvs_alt_genomic)
                corrective_action_taken = 'false'

        # In this instance, the gap code has generally found an incomplete-alignment rather than a
        # truly gapped alignment.
        except KeyError:
            warnings = warnings + ': Suspected incomplete alignment between transcript %s and ' \
                                  'genomic reference sequence %s' % (hgvs_coding.ac,
                                                                     alt_chr)
            continue
        except hgvs.exceptions.HGVSError as e:
            logger.error(str(e))
            logger.debug(str(e))
            continue

    if multi_g != []:

        multi_gen_vars = multi_g  # '|'.join(multi_g)
    else:
        multi_gen_vars = []

    return multi_gen_vars, hgvs_coding
