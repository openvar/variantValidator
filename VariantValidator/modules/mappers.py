import hgvs
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

