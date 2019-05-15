import hgvs
import re
import copy
import hgvs.exceptions
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
    # rel_var is a key-worded list of relevant transcripts with associated coding variants
    """
    Initial simple projection from the provided g. position all overlapping
    transcripts
    """
    rel_var = validator.relevant_transcripts(variant.hgvs_genomic, variant.evm, validator.alt_aln_method,
                                             variant.reverse_normalizer)

    # Double check rel_vars have not been missed when mapping from a RefSeqGene
    if len(rel_var) != 0 and 'NG_' in variant.hgvs_genomic.ac:
        for var in rel_var:
            hgvs_coding_variant = validator.hp.parse_hgvs_variant(var)
            try:
                variant.hgvs_genomic = validator.myevm_t_to_g(hgvs_coding_variant, variant.no_norm_evm,
                                                              variant.primary_assembly, variant.hn)
            except hgvs.exceptions.HGVSError:
                try_rel_var = []
            else:
                try_rel_var = validator.relevant_transcripts(variant.hgvs_genomic, variant.evm,
                                                             validator.alt_aln_method, variant.reverse_normalizer)
            if len(try_rel_var) > len(rel_var):
                rel_var = try_rel_var
                break
            else:
                continue

    #  Tripple check this assumption by querying the gene position database
    if len(rel_var) == 0:
        vcf_dict = vvHGVS.hgvs2vcf(variant.hgvs_genomic, variant.primary_assembly, variant.reverse_normalizer,
                                   validator.sf)
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
                genomic_input = refseqgene_data['hgvs_genomic']
                # re_submit
                # Tag the line so that it is not written out
                variant.warnings += ': ' + str(variant.hgvs_formatted) + ' automapped to genome position ' + \
                                    str(genomic_input)
                query = Variant(variant.original, quibble=genomic_input, warnings=variant.warnings,
                                primary_assembly=variant.primary_assembly, order=variant.order)

                validator.batch_list.append(query)
            else:
                error = 'Mapping unavailable for RefSeqGene ' + str(variant.hgvs_formatted) + \
                        ' using alignment method = ' + validator.alt_aln_method
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
                    # Example {'gene': 'NTHL1', 'hgvs_refseqgene': 'NG_008412.1:g.3455_3464delCAAACACACA',
                    # 'valid': 'true'}
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
                error = 'Please ensure the requested chromosome version relates to a supported genome build. ' \
                        'Supported genome builds are: GRCh37, GRCh38, hg19 and hg38'
                variant.warnings += ': ' + str(error)
                logger.warning(str(error))
                return True

    else:
        # Tag the line so that it is not written out
        variant.write = False

        gap_mapper = gapped_mapping.GapMapper(variant, validator)

        data, nw_rel_var = gap_mapper.gapped_g_to_c(rel_var)

        # # Warn the user that the g. description is not valid
        # if data['gapped_alignment_warning'] != '':
        #     if data['disparity_deletion_in'][0] == 'transcript':
        #         corrective_action_taken = 'Automap has deleted  ' + str(
        #             data['disparity_deletion_in'][1]) + ' bp from chromosomal reference sequence ' + str(
        #             variant.hgvs_genomic.ac) + ' to ensure perfect alignment with transcript reference sequence(s)'\
        #                                   + data['gapped_transcripts']
        #     if data['disparity_deletion_in'][0] == 'chromosome':
        #         corrective_action_taken = 'Automap has added  ' + str(
        #             data['disparity_deletion_in'][1]) + ' bp to chromosomal reference sequence ' + str(
        #             variant.hgvs_genomic.ac) + ' to ensure perfect alignment with transcript reference sequence(s) '\
        #                                   + data['gapped_transcripts']
        #
        # # Add additional data to the front of automap
        # if data['auto_info'] != '':
        #     automap = data['auto_info'] + '\n' + 'false'

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
    caution = ''
    error = ''
    # Collect information for genomic level validation
    obj = validator.hp.parse_hgvs_variant(str(variant.hgvs_formatted))

    tx_ac = obj.ac

    quibble_input = str(variant.quibble)
    formatted_variant = str(variant.hgvs_formatted)

    # Do we keep it?
    if validator.select_transcripts != 'all':
        if tx_ac not in list(validator.select_transcripts_dict_plus_version.keys()):
            # By marking it as Do Not Write and continuing through the validation loop
            variant.write = False
            return True

    # Se rec_var to '' so it can be updated later
    rec_var = ''

    # First task is to get the genomic equivalent, and print useful error messages if it can't be found.
    try:
        to_g = validator.myevm_t_to_g(obj, variant.no_norm_evm, variant.primary_assembly, variant.hn)
        genomic_ac = to_g.ac
    except hgvs.exceptions.HGVSDataNotAvailableError as e:
        if ('~' in str(e) and 'Alignment is incomplete' in str(e)) or "No relevant genomic mapping options" in str(e):
            # Unable to map the input variant onto a genomic position
            if '~' in str(e) and 'Alignment is incomplete' in str(e):
                error = 'Full alignment data between the specified transcript reference sequence and all GRCh37 ' \
                        'and GRCh38 genomic reference sequences (including alternate chromosome assemblies, ' \
                        'patches and RefSeqGenes) are not available: Consequently the input variant description ' \
                        'cannot be fully validated and is not supported: Use the Gene to Transcripts function to ' \
                        'determine whether an updated transcript reference sequence is available'
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
    minus = re.compile(r"\d-\d")  # finds digit - digit

    if plus.search(quibble_input) or minus.search(quibble_input):
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

    elif ':g.' in quibble_input:
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

    else:
        # Normalize the variant
        try:
            h_variant = variant.hn.normalize(obj)
        except hgvs.exceptions.HGVSUnsupportedOperationError as e:
            error = str(e)
            if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                formatted_variant = formatted_variant
                caution = 'This coding sequence variant description spans at least one intron'
                automap = 'Use of the corresponding genomic sequence variant descriptions may be invalid. ' \
                          'Please refer to https://www35.lamp.le.ac.uk/recommendations/'
                variant.warnings += ': ' + caution + ': ' + automap
                logger.warning(caution + ": " + automap)
        else:
            formatted_variant = str(h_variant)

        error = validator.validateHGVS(formatted_variant)
        if error == 'false':
            valid = True
        else:
            variant.warnings += ': ' + str(error)
            logger.warning(str(error))
            return True

    # Tackle the plus intronic offset
    cck = False
    if plus.search(quibble_input):
        # Regular expression catches the start of the interval only based on .00+00 pattern
        inv_start = re.compile(r"\.\d+\+\d")
        if inv_start.search(quibble_input):
            cck = True
    if minus.search(quibble_input):
        # Regular expression catches the start of the interval only based on .00-00 pattern
        inv_start = re.compile(r"\.\d+-\d")
        if inv_start.search(quibble_input):
            cck = True

    # COORDINATE CHECKER
    # hgvs will handle incorrect coordinates so need to automap errors
    # Make sure any input intronic coordinates are correct
    # Get the desired transcript

    if cck:
        # This should only ever hit coding and RNA variants
        if 'del' in formatted_variant:
            # RNA - looking at trapped variant which was saved before RNA converted to cDNA
            if ':r.' in variant.pre_RNA_conversion:
                coding = validator.coding(formatted_variant, validator.hp)
                trans_acc = coding.ac
                # c to Genome coordinates - Map the variant to the genome
                pre_var = validator.genomic(formatted_variant, variant.no_norm_evm, variant.primary_assembly,
                                            variant.hn)
                # genome back to C coordinates
                post_var = validator.myevm_g_to_t(variant.evm, pre_var, trans_acc)

                test = validator.hp.parse_hgvs_variant(quibble_input)
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
                    automap = variant.pre_RNA_conversion + ' automapped to ' + str(post_var)
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
                test = validator.hp.parse_hgvs_variant(quibble_input)

                if post_var.posedit.pos.start.base != test.posedit.pos.start.base or \
                        post_var.posedit.pos.end.base != test.posedit.pos.end.base:
                    caution = 'The entered coordinates do not agree with the intron/exon boundaries for the ' \
                              'selected transcript:'
                    # automapping of variant completed
                    automap = variant.pre_RNA_conversion + ' automapped to ' + str(post_var)
                    variant.warnings += str(caution) + ': ' + str(automap)

                    # Kill current line and append for re-submission
                    # Tag the line so that it is not written out
                    variant.write = False
                    # Set the values and append to batch_list
                    hgvs_vt = validator.hp.parse_hgvs_variant(str(post_var))
                    assert str(hgvs_vt) == str(post_var)
                    query = Variant(variant.original, quibble=fn.valstr(hgvs_vt), warnings=automap,
                                    primary_assembly=variant.primary_assembly, order=variant.order)
                    validator.batch_list.append(query)

        else:  # del not in formatted_variant
            if ':r.' in variant.pre_RNA_conversion:
                coding = validator.coding(formatted_variant, validator.hp)
                trans_acc = coding.ac
                # c to Genome coordinates - Map the variant to the genome
                pre_var = validator.genomic(formatted_variant, variant.no_norm_evm, variant.primary_assembly,
                                            variant.hn)
                # genome back to C coordinates
                post_var = validator.myevm_g_to_t(variant.evm, pre_var, trans_acc)

                test = validator.hp.parse_hgvs_variant(quibble_input)
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
                    automap = quibble_input + ' automapped to ' + post_var
                    variant.warnings += ': ' + str(caution) + ': ' + str(
                        automap)

                    # Kill current line and append for re-submission
                    # Tag the line so that it is not written out
                    variant.write = False
                    # Set the values and append to batch_list
                    hgvs_vt = validator.hp.parse_hgvs_variant(str(post_var))
                    assert str(hgvs_vt) == str(post_var)
                    query = Variant(variant.original, quibble=fn.valstr(hgvs_vt), warnings=automap,
                                    primary_assembly=variant.primary_assembly, order=variant.order)
                    validator.batch_list.append(query)

            else:
                coding = validator.coding(formatted_variant, validator.hp)
                trans_acc = coding.ac
                # c to Genome coordinates - Map the variant to the genome
                pre_var = validator.genomic(formatted_variant, variant.no_norm_evm, variant.primary_assembly,
                                            variant.hn)

                # genome back to C coordinates
                post_var = validator.myevm_g_to_t(variant.evm, pre_var, trans_acc)

                test = validator.hp.parse_hgvs_variant(quibble_input)
                if post_var.posedit.pos.start.base != test.posedit.pos.start.base or \
                        post_var.posedit.pos.end.base != test.posedit.pos.end.base:
                    caution = 'The entered coordinates do not agree with the intron/exon boundaries for the ' \
                              'selected transcript:'
                    # automapping of variant completed
                    automap = str(variant.pre_RNA_conversion) + ' automapped to ' + str(post_var)
                    variant.warnings += ': ' + str(caution) + ': ' + str(
                        automap)

                    # Kill current line and append for re-submission
                    # Tag the line so that it is not written out
                    variant.write = False
                    # Set the values and append to batch_list
                    hgvs_vt = validator.hp.parse_hgvs_variant(str(post_var))
                    assert str(hgvs_vt) == str(post_var)
                    query = Variant(variant.original, quibble=fn.valstr(hgvs_vt), warnings=automap,
                                    primary_assembly=variant.primary_assembly, order=variant.order)
                    validator.batch_list.append(query)

    # If cck not true
    elif ':r.' in variant.pre_RNA_conversion:
        # set input hgvs object
        # Traps the hgvs variant of r. for further use
        hgvs_rna_input = validator.hp.parse_hgvs_variant(variant.pre_RNA_conversion)
        inp = str(validator.hgvs_r_to_c(hgvs_rna_input))
        # Regex
        if plus.search(quibble_input) or minus.search(quibble_input):
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
            except hgvs.exceptions.HGVSError:
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
            caution = 'The variant description ' + quibble_input + ' requires alteration to comply with HGVS variant ' \
                                                           'nomenclature:'
            # automapping of variant completed
            automap = variant.pre_RNA_conversion + ' automapped to ' + output
            variant.warnings += ': ' + caution + ': ' + automap

            # Kill current line and append for re-submission
            # Tag the line so that it is not written out
            variant.write = False
            # Set the values and append to batch_list
            hgvs_vt = validator.hp.parse_hgvs_variant(str(query))
            assert str(hgvs_vt) == str(query)
            query = Variant(variant.original, quibble=fn.valstr(hgvs_vt), warnings=automap,
                            primary_assembly=variant.primary_assembly, order=variant.order)
            validator.batch_list.append(query)

    elif ':g.' not in quibble_input:
        query = validator.hp.parse_hgvs_variant(formatted_variant)
        test = validator.hp.parse_hgvs_variant(quibble_input)
        if query.posedit.pos != test.posedit.pos:
            caution = 'The variant description ' + quibble_input + ' requires alteration to comply with HGVS variant ' \
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
            query = Variant(variant.original, quibble=fn.valstr(hgvs_vt), warnings=automap,
                            primary_assembly=variant.primary_assembly, order=variant.order)
            validator.batch_list.append(query)

    # VALIDATION of intronic variants
    pre_valid = validator.hp.parse_hgvs_variant(quibble_input)
    post_valid = validator.hp.parse_hgvs_variant(formatted_variant)

    # valid is false if the input contains a \d+\d, \d-\d or :g.
    if not valid:
        genomic_validation = str(validator.genomic(quibble_input, variant.no_norm_evm, variant.primary_assembly,
                                                   variant.hn))
        if fn.valstr(pre_valid) != fn.valstr(post_valid):
            if variant.reftype != ':g.':
                if caution == '':
                    caution = fn.valstr(pre_valid) + ' automapped to ' + fn.valstr(post_valid)
                variant.warnings += ': ' + caution
                logger.warning(caution)

        # Apply validation to intronic variant descriptions (should be valid but make sure)
        error = validator.validateHGVS(genomic_validation)
        if error != 'false':
            variant.warnings += ': ' + error
            return True

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
        if sequence[1].startswith('NG_'):
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
        except Exception:
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
            except hgvs.exceptions.HGVSInvalidIntervalError:
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
        # Double check protein position by reverse_norm genomic, and normalize back to c.
        # for normalize or not to normalize issue
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
            except hgvs.exceptions.HGVSError:
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
            if 'does not agree with reference sequence' in error:
                match = re.findall(r'\(([GATC]+)\)', error)
                new_ref = match[1]
                hgvs_updated.posedit.edit.ref = new_ref
                validator.vr.validate(hgvs_updated)

        updated_transcript_variant = hgvs_updated
        variant.warnings += ': ' + 'A more recent version of the selected reference sequence ' + hgvs_coding.ac + \
                            ' is available (' + updated_transcript_variant.ac + ')' + ': ' + \
                            str(updated_transcript_variant) + ' MUST be fully validated prior to use in reports: ' \
                            'select_variants=' + fn.valstr(updated_transcript_variant)

    variant.coding = str(hgvs_coding)
    variant.genomic_r = str(hgvs_refseq)
    variant.genomic_g = str(hgvs_genomic)
    variant.protein = str(hgvs_protein)

    return False


def final_tx_to_multiple_genomic(variant, validator, tx_variant):

    warnings = ''
    rec_var = ''
    gap_compensation = True

    # Multiple genomic variants
    # multi_gen_vars = []
    variant.hgvs_coding = validator.hp.parse_hgvs_variant(str(tx_variant))
    # Gap gene black list
    try:
        gene_symbol = validator.db.get_gene_symbol_from_transcriptID(variant.hgvs_coding.ac)
    except Exception:
        fn.exceptPass()
    else:
        # If the gene symbol is not in the list, the value False will be returned
        gap_compensation = vvChromosomes.gap_black_list(gene_symbol)

    # Look for variants spanning introns
    try:
        variant.hn.normalize(variant.hgvs_coding)
    except hgvs.exceptions.HGVSUnsupportedOperationError as e:
        error = str(e)
        if 'boundary' in error or 'spanning' in error:
            gap_compensation = False

    except hgvs.exceptions.HGVSError:
        fn.exceptPass()

    # Warn gap code status
    logger.warning("gap_compensation_3 = " + str(gap_compensation))
    multi_g = []
    multi_list = []
    mapping_options = validator.hdp.get_tx_mapping_options(variant.hgvs_coding.ac)
    for alt_chr in mapping_options:
        if ('NC_' in alt_chr[1] or 'NT_' in alt_chr[1] or 'NW_' in alt_chr[1]) and \
                alt_chr[2] == validator.alt_aln_method:
            multi_list.append(alt_chr[1])

    for alt_chr in multi_list:
        try:
            # Re set ori
            ori = validator.tx_exons(tx_ac=variant.hgvs_coding.ac, alt_ac=alt_chr,
                                     alt_aln_method=validator.alt_aln_method)

            hgvs_alt_genomic = validator.myvm_t_to_g(variant.hgvs_coding, alt_chr, variant.no_norm_evm, variant.hn)

            gap_mapper = gapped_mapping.GapMapper(variant, validator)

            # Loop out gap code under these circumstances!
            if gap_compensation:
                hgvs_alt_genomic, hgvs_coding = gap_mapper.g_to_t_gap_compensation_version3(
                    hgvs_alt_genomic, variant.hgvs_coding, ori, alt_chr, rec_var)
                variant.hgvs_coding = hgvs_coding

                # Refresh the :g. variant
                multi_g.append(hgvs_alt_genomic)
            else:
                multi_g.append(hgvs_alt_genomic)

        # In this instance, the gap code has generally found an incomplete-alignment rather than a
        # truly gapped alignment.
        except KeyError:
            warnings = warnings + ': Suspected incomplete alignment between transcript %s and ' \
                                  'genomic reference sequence %s' % (variant.hgvs_coding.ac, alt_chr)
        except hgvs.exceptions.HGVSError as e:
            logger.error(str(e))
            logger.debug(str(e))

    return multi_g
