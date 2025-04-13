import re
import vvhgvs.exceptions
import copy
import logging
from VariantValidator.modules.variant import Variant
from VariantValidator.modules import seq_data, initial_formatting
from VariantValidator.modules import utils as fn
import VariantValidator.modules.rna_formatter
from VariantValidator.modules import complex_descriptions, use_checking, \
        expanded_repeats, methyl_syntax
from VariantValidator.modules.vvMixinConverters import AlleleSyntaxError
from VariantValidator.modules.hgvs_utils import hgvs_delins_parts_to_hgvs_obj,\
        unset_hgvs_obj_ref

logger = logging.getLogger(__name__)


def initial_format_conversions(variant, validator, select_transcripts_dict_plus_version):

    # VCF type 1
    toskip = vcf2hgvs_stage1(variant, validator)
    if toskip:
        return True

    # API type non-HGVS
    # e.g. Chr16:2099572TC>T
    toskip = vcf2hgvs_stage2(variant, validator)
    if toskip:
        return True

    toskip = vcf2hgvs_stage3(variant, validator)
    if toskip:
        return True

    toskip = gene_symbol_catch(variant, validator, select_transcripts_dict_plus_version)
    if toskip:
        return True

    # NG_:c. or NC_:c.
    toskip = refseq_catch(variant, validator, select_transcripts_dict_plus_version)
    if toskip:
        return True

    # Find not_sub type in input e.g. GGGG>G
    toskip = vcf2hgvs_stage4(variant, validator)
    if toskip:
        return True

    # Extract variants from HGVS allele descriptions
    # http://varnomen.hgvs.org/recommendations/DNA/variant/alleles/
    toskip = allele_parser(variant, validator, validator)
    if toskip:
        return True

    # Conversions
    # are not currently supported. The HGVS format for conversions
    # is rarely seen wrt genomic sequencing data and needs to be re-evaluated
    # so abort before hgvs object conversion to avoid errors & parsing overhead
    if 'con' in str(variant.quibble):
        variant.warnings.append('Conversions are no longer valid HGVS Sequence Variant Descriptions')
        logger.warning('Conversions are no longer valid HGVS Sequence Variant Descriptions')
        return True

    toskip = use_checking.pre_parsing_global_common_mistakes(variant)
    if toskip:
        return True

    # Remove & store Methylation Syntax suffix before hgvs object parsing
    methyl_syntax.methyl_syntax(variant)

    # Uncertain positions (converts to hgvs object so must be post split/fix)
    toskip = uncertain_pos(variant, validator)
    if toskip:
        return True

    # Expanded repeat->delins code can not handle Uncertain positions yet
    # also does hgvs object conversion if it triggers
    if type(variant.quibble) is str: # not Uncertain
        toskip = convert_expanded_repeat(variant, validator)
        if toskip:
            return True

    # Catches del12/ins21 type variants, can usfully trigger on the outupt of the expanded repeat conversions
    # and can not currently handle uncertain positions
    if type(variant.quibble) is str: # not Uncertain
        toskip = indel_catching(variant, validator)
        if toskip:
            return True

    # Quibble should now be correctly formatted hgvs & work for object parsing
    # or else already have been parsed already
    if type(variant.quibble) is str:
        try:
            toskip = final_hgvs_convert(variant, validator)
        except:
            # Check for common mistakes
            toskip = use_checking.refseq_common_mistakes(variant)
            if toskip:
                return True
            # try again if corrected
            try:
                toskip = final_hgvs_convert(variant, validator)
            except vvhgvs.exceptions.HGVSParseError as err:
                variant.warnings.append("HgvsSyntaxError: " + str(err))
                return True
            except vvhgvs.exceptions.HGVSError as err:
                variant.warnings.append(f"HgvsParserError: Unknown error during"
                                        "reading of variant {variant.quibble}")
                return True
        # fail if un-corrected errors persist (warning should already have been generated)
        if toskip:
            return True

    # Tackle compound variant descriptions NG or NC (NM_) i.e. correctly input NG/NC_(NM_):c.
    intronic_converter(variant, validator)
    return False

def final_hgvs_convert(variant,validator):
    """
    For use in the final hgvs str ->hgvs obj conversion.
    Requires a fully checked out text variant quibble
    Avoids issues with XX_000XX(XX_000XX): type variants by parsing more
    directly
    Returns skipvar i.e true if something went wrong
    """
    seq_ac, _sep, type_posedit = variant.quibble.partition(':')
    var_type, _sep, posedit = type_posedit.partition('.')
    if var_type == 'c':
        posedit = validator.hp.parse_c_posedit(posedit)
    elif var_type == 'g':
        posedit = validator.hp.parse_g_posedit(posedit)
    elif var_type == 'm':
        posedit = validator.hp.parse_m_posedit(posedit)
    elif var_type == 'n':
        posedit = validator.hp.parse_n_posedit(posedit)
    elif var_type == 'p':
        posedit = validator.hp.parse_p_posedit(posedit)
    elif var_type == 'r':
        if 'T' in posedit:
            e = 'The IUPAC RNA alphabet dictates that RNA variants must use '+\
                    'the character u in place of t'
            variant.warnings.append(e)
            return True
        posedit = validator.hp.parse_r_posedit(posedit)
    else:
        e = "VariantSyntaxError: The detected variant sequence type of "+\
                f"{var_type} ' was not one of the allowed HGVS type "+\
                "characters of c, g, m, n, p, or r"
        variant.warnings.append(e)
        return True

    variant.quibble = vvhgvs.sequencevariant.SequenceVariant(
            ac = seq_ac,
            type = var_type,
            posedit = posedit
            )
    return False


def vcf2hgvs_stage1(variant, validator):
    """
    VCF2HGVS stage 1. converts chr-pos-ref-alt into chr:posRef>Alt
    The output format is a common mistake caused by inaccurate conversion of
    VCF variants into HGVS - hence the need for conversion step 2
    """
    skipvar = False
    vcf_data = re.split(r'[-:]',variant.quibble)
    if len(vcf_data) < 4:
        logger.debug("Completed VCF-HVGS step 1 for %s", variant.quibble)
        return False
    poss_genome = vcf_data[0].lower()
    if ('grch3' in poss_genome or 'hg' in poss_genome) and poss_genome[-1].isdigit():
        vcf_data = vcf_data[1:]
        variant.quibble = '-'.join(vcf_data)
        # TODO test assembly given against settings
        if len(vcf_data) < 4:
            variant.warnings.append("Insufficient or incorrect  VCF elements provided. "
                                    "Elements required are chr-pos-ref-alt")
            return True
    # no coordinate found
    if not vcf_data[1].isdigit():
        logger.debug("Completed VCF-HVGS step 1 for %s", variant.quibble)
        return False
    # ref is present and not . or DNA, or alt is present and not DNA or a list of DNA possibilities
    # for now we leave ',' in, this is used to do vcf type multi-alt variants. We correct this as a
    # later step since we may get psudo-hgvs input with the same pattern too.
    if vcf_data[2] == '.':
        vcf_data[2] = ''
    if  vcf_data[3] == '.':
        vcf_data[3] = ''

    if re.search(r'[^CGAT]',vcf_data[2]) or re.search(r'[^CGAT,]',vcf_data[3]):
        logger.debug("Completed VCF-HVGS step 1 for %s", variant.quibble)
        return False
    if vcf_data[2] and vcf_data[3]:
        variant.quibble = f'{vcf_data[0]}:{vcf_data[1]}{vcf_data[2]}>{vcf_data[3]}'
    elif vcf_data[2]:
        variant.warnings = ['Not stating ALT bases is ambiguous because VCF specification 4.0 would treat ' +
                            variant.quibble + ' as a deletion whereas VCF specification 4.1 onwards would treat ' +
                            variant.quibble + ' as ALT = REF']
        variant.warnings.append('VariantValidator has output both alternatives')
        logger.info('Not stating ALT bases is ambiguous because VCF specification 4.0 would treat %s as a deletion '
                    'whereas VCF specification 4.1 onwards would treat %s as ALT = REF. Validator will output '
                    'both alternatives.', variant.quibble, variant.quibble)
        variant.write = False
        input_a = '%s:%s%s>%s' % (vcf_data[0], vcf_data[1], vcf_data[2], 'del')
        input_b = '%s:%s%s>%s' % (vcf_data[0], vcf_data[1], vcf_data[2], vcf_data[2])
        query_a = Variant(variant.original, quibble=input_a, warnings=variant.warnings,
                          primary_assembly=variant.primary_assembly, order=variant.order)
        query_b = Variant(variant.original, quibble=input_b, warnings=variant.warnings,
                          primary_assembly=variant.primary_assembly, order=variant.order)
        validator.batch_list.append(query_a)
        validator.batch_list.append(query_b)
        logger.info("Submitting new variant with format %s", input_a)
        logger.info("Submitting new variant with format %s", input_b)
        skipvar = True
    elif vcf_data[3]:
        variant.quibble = f'{vcf_data[0]}:{vcf_data[1]}ins{vcf_data[3]}'

    logger.debug("Completed VCF-HVGS step 1 for %s", variant.quibble)
    return skipvar


def vcf2hgvs_stage2(variant, validator):
    """
    VCF2HGVS conversion step 2 identifies the correct chromosomal reference
    sequence based upon the non compliant identifier e.g. <Chr16>:2099572TC>T.
    The data is currently stored in variantanalyser.supported_chromosome_builds.
    Anticipated future builds will be transferred to MySQL which can be more
    easily updated and maintained.
    LRGs and LRG_ts also need to be assigned the correct reference sequence identifier.
    The LRG ID data ia stored in the VariantValidator MySQL database.
    The reference sequence type is also assigned.
    """
    skipvar = False
    if (re.search(r'\w+:', variant.quibble) or re.search(r'\w+\(\w+\):', variant.quibble)) and not \
            (re.search(r'\w+:[gcnmrpoGCNMRPO]\.', variant.quibble) or re.search(r'\w+\(\w+\):[gcnmrpoGCNMRPO]\.',
                                                                              variant.quibble)):
        if re.search(r'\w+:[gcnmrpo]', variant.quibble) and not re.search(r'\w+:[gcnmrpo]\.', variant.quibble):
            # Missing dot
            pass
        else:
            try:
                if 'GRCh37' in variant.quibble or 'hg19' in variant.quibble:
                    variant.primary_assembly = 'GRCh37'
                    validator.selected_assembly = 'GRCh37'
                    variant.format_quibble()
                elif 'GRCh38' in variant.quibble or 'hg38' in variant.quibble:
                    variant.primary_assembly = 'GRCh38'
                    validator.selected_assembly = 'GRCh38'
                    variant.format_quibble()
                # Remove all content in brackets
                input_list = variant.quibble.split(':')
                pos_ref_alt = str(input_list[1])
                position_and_edit = input_list[1]
                if not re.match(r'N[CGTWMRP]_', variant.quibble) and not re.match(r'LRG_', variant.quibble):
                    chr_num = str(input_list[0])
                    chr_num = chr_num.upper().strip()
                    if re.match('CHR', chr_num):
                        chr_num = chr_num.replace('CHR', '')
                    # Use selected assembly
                    accession = seq_data.to_accession(chr_num, validator.selected_assembly)
                    if accession is None:
                        variant.warnings.append(chr_num + ' is not part of genome build ' + validator.selected_assembly)
                        logger.warning(chr_num + ' is not part of genome build ' + validator.selected_assembly)
                        skipvar = True
                else:
                    accession = input_list[0]
                if '>' in variant.quibble:
                    if 'del' in variant.quibble:
                        pos = re.match(r'\d+', pos_ref_alt)
                        position = pos.group(0)
                        old_ref, old_alt = pos_ref_alt.split('>')
                        old_ref = old_ref.replace(position, '')
                        position = int(position) - 1
                        required_base = validator.sf.fetch_seq(accession, start_i=position - 1, end_i=position)
                        ref = required_base + old_ref
                        alt = required_base
                        position_and_edit = str(position) + ref + '>' + alt
                    elif 'ins' in variant.quibble:
                        pos = re.match(r'\d+', pos_ref_alt)
                        position = pos.group(0)
                        old_ref, old_alt = pos_ref_alt.split('>')
                        # old_ref = old_ref.replace(position, '')
                        position = int(position) - 1
                        required_base = validator.sf.fetch_seq(accession, start_i=position - 1, end_i=position)
                        ref = required_base
                        alt = required_base + old_alt
                        position_and_edit = str(position) + ref + '>' + alt

                # Assign reference sequence type
                if accession in ["NC_012920.1", "NC_001807.4"]:
                    ref_type = ":m."
                else:
                    ref_type = validator.db.ref_type_assign(accession)

                # Sort LRG formatting
                if re.match('LRG_', accession):
                    if ref_type == ':g.':
                        accession = validator.db.get_refseq_id_from_lrg_id(accession)
                    else:
                        accession = validator.db.get_refseq_transcript_id_from_lrg_transcript_id(accession)
                else:
                    accession = accession
                variant.quibble = str(accession) + ref_type + str(position_and_edit)

            except Exception as e:
                logger.debug("Except passed, %s", e)

    # Descriptions lacking the colon : or the dot .
    if re.search(r'[gcnmrp]\.', variant.quibble) and not re.search(r':[gcnmrp]\.', variant.quibble):
        error = 'Unable to identify a colon (:) in the variant description %s. A colon is required in HGVS variant ' \
                'descriptions to separate the reference accession from the reference type i.e. <accession>:<type>. ' \
                'e.g. :c.' % variant.quibble
        variant.warnings.append(error)
        logger.warning(error)
        skipvar = True
    elif re.search(r':[gcnmrp]', variant.quibble) and not re.search(r':[gcnmrp]\.', variant.quibble):
        error = 'Unable to identify a dot (.) in the variant description %s following the reference sequence ' \
                'type (g,c,n,r, or p). A dot is required in HGVS variant ' \
                'descriptions to separate the reference type from the variant position i.e. <accession>:<type>. ' \
                'e.g. :g.' % variant.quibble
        variant.warnings.append(error)
        logger.warning(error)
        skipvar = True
    elif re.search(r':[GCNMRPO]\.', variant.quibble):
        error = 'Reference type incorrectly stated in the variant description %s ' \
                'Valid types are g,c,n,r, or p' % variant.quibble
        variant.warnings.append(error)
        logger.warning(error)
        match = re.search(r':[GCNMRPO]\.', variant.quibble)[0]
        variant.quibble = variant.quibble.replace(match, match.lower())

    # Ambiguous chr reference
    logger.debug("Completed VCF-HVGS step 2 for %s", variant.quibble)
    return skipvar


def vcf2hgvs_stage3(variant, validator):
    """
    VCF2HGVS conversion step 3 is similar to step 2 but handles
    formats like Chr16:g.2099572TC>T which are provided by Alamut and other
    software
    """
    skipvar = False
    if (re.search(r'\w+:[gcnmrpGCMNRP]\.', variant.quibble) or re.search(r'\w+\(\w+\):[gcnmrpGCMNRP]\.',
                                                                         variant.quibble)) \
            and not re.match(r'N[CGTWMRP]_', variant.quibble):

        # Take out lowercase Accession characters
        lower_cased_list = variant.quibble.split(':')
        if re.search('LRG', lower_cased_list[0], re.IGNORECASE):
            lower_case_accession = lower_cased_list[0]
            lower_case_accession = lower_case_accession.replace('l', 'L')
            lower_case_accession = lower_case_accession.replace('r', 'R')
            lower_case_accession = lower_case_accession.replace('g', 'G')
        else:
            lower_case_accession = lower_cased_list[0]
            lower_case_accession = lower_case_accession.upper()
        variant.quibble = ''.join(lower_cased_list[1:])
        variant.quibble = lower_case_accession + ':' + variant.quibble
        if 'LRG_' not in variant.quibble and 'ENS' not in variant.quibble and not re.match('N[MRPC]_', variant.quibble):
            try:
                if re.search('GRCh37', variant.quibble, re.IGNORECASE) or \
                        re.search('hg19', variant.quibble, re.IGNORECASE):
                    variant.primary_assembly = 'GRCh37'
                    validator.selected_assembly = 'GRCh37'
                    variant.format_quibble()
                if re.search('GRCh38', variant.quibble, re.IGNORECASE) or \
                        re.search('hg38', variant.quibble, re.IGNORECASE):
                    variant.primary_assembly = 'GRCh38'
                    validator.selected_assembly = 'GRCh38'
                    variant.format_quibble()
                input_list = variant.quibble.split(':')
                query_a_symbol = input_list[0]
                is_it_a_gene = validator.db.get_hgnc_symbol(query_a_symbol)
                if is_it_a_gene == 'none':
                    position_and_edit = input_list[1]
                    chr_num = str(input_list[0])
                    chr_num = chr_num.upper().strip()
                    if re.match('CHR', chr_num):
                        chr_num = chr_num.replace('CHR', '')  # Use selected assembly
                    accession = seq_data.to_accession(chr_num, validator.selected_assembly)
                    if accession is None:
                        variant.warnings.append(chr_num + ' is not part of genome build ' + validator.selected_assembly)
                        logger.warning(chr_num + ' is not part of genome build ' + validator.selected_assembly)
                        skipvar = True
                    variant.quibble = str(accession) + ':' + str(position_and_edit)
            except Exception as e:
                logger.debug("Except passed, %s", e)

    logger.debug("Completed VCF-HGVS step 3 for %s", variant.quibble)
    return skipvar


def gene_symbol_catch(variant, validator, select_transcripts_dict_plus_version):
    """
    Searches for gene symbols that have been used as reference sequence
    identifiers. Provides a sufficiently repremanding warning, but also provides
    correctly formatted variant descriptions with appropriate transcript
    reference sequence identifiers i.e. NM_ ....
    Note: the output from the function must be validated because VV has no way
    of knowing which the users intended reference sequence was, and the exon
    boundaries etc of the alternative transcript variants may not be equivalent
    """
    skipvar = False
    query_a_symbol, _sep, tx_edit = variant.quibble.partition(':')
    if not (query_a_symbol[:3] in ['NC','NM_','NR_','ENST'] or
            query_a_symbol.startswith('ENST')
            ) and re.search(r'\w+:[cn]\.', variant.quibble):
        try:
            is_it_a_gene = validator.db.get_hgnc_symbol(query_a_symbol)
            if is_it_a_gene != 'none':
                uta_symbol = validator.db.get_uta_symbol(is_it_a_gene)
                available_transcripts = validator.hdp.get_tx_for_gene(uta_symbol)
                select_from_these_transcripts = []
                for tx in available_transcripts:
                    if validator.alt_aln_method == 'splign' and ('NM_' in tx[3] or 'NR_' in tx[3]):
                        if tx[3] not in select_from_these_transcripts:
                            select_from_these_transcripts.append(tx[3])
                    elif validator.alt_aln_method == 'genebuild' and 'ENST' in tx[3]:
                        if tx[3] not in select_from_these_transcripts:
                            select_from_these_transcripts.append(tx[3])
                select_from_these_transcripts = '|'.join(select_from_these_transcripts)
                if validator.select_transcripts != 'all' and validator.select_transcripts != 'raw':
                    variant.write = False
                    for transcript in list(select_transcripts_dict_plus_version.keys()):
                        if transcript == "mane":
                            for tx in select_from_these_transcripts.split('|'):
                                annotation = validator.db.get_transcript_annotation(tx)
                                if '"mane_select": true' in annotation or '"mane_plus_clinical": true' in annotation:
                                    transcript = tx
                                else:
                                    continue
                        elif transcript == "mane_select":
                            for tx in select_from_these_transcripts.split('|'):
                                annotation = validator.db.get_transcript_annotation(tx)
                                if '"mane_select": true' in annotation:
                                    transcript = tx
                                else:
                                    continue

                        variant.warnings.append('InvalidReferenceError: HGVS variant nomenclature does not '
                                                'allow the use of a gene symbol (' +
                                            query_a_symbol + ') in place of a valid reference sequence')
                        refreshed_description = transcript + ':' + tx_edit
                        query = Variant(variant.original, quibble=refreshed_description,
                                        warnings=variant.warnings, primary_assembly=variant.primary_assembly,
                                        order=variant.order)
                        validator.batch_list.append(query)
                        logger.info('HGVS variant nomenclature does not allow the use of a gene symbol (' +
                                    query_a_symbol + ') in place of a valid reference sequence')
                        logger.info("Submitting new variant with format %s", refreshed_description)
                else:
                    variant.warnings.append('InvalidReferenceError: HGVS variant nomenclature does not allow '
                                            'the use of a gene symbol ('
                                            + query_a_symbol + ') in place of a valid reference sequence: Re-submit ' +
                                            str(variant.quibble) + ' and specify transcripts from the following: ' +
                                            'select_transcripts=' + select_from_these_transcripts)
                    logger.warning('HGVS variant nomenclature does not allow the use of a gene symbol (' +
                                   query_a_symbol + ') in place of a valid reference sequence: Re-submit ' +
                                   str(variant.quibble) + ' and specify transcripts from the following: ' +
                                   'select_transcripts=' + select_from_these_transcripts)
                skipvar = True
        except Exception as e:
            logger.debug("Except passed, %s", e)
    logger.debug("Gene symbol reference catching complete")
    return skipvar


def refseq_catch(variant, validator, select_transcripts_dict_plus_version):
    """
    Similar to the GENE_SYMBOL:c. n. types function, but spots RefSeqGene or
    Chromosomal reference sequence identifiers used in the context of c. variant
    descriptions
    """
    skipvar = False
    if re.search(r'\w+:[cn]', variant.quibble):
        try:
            if variant.quibble.startswith('NG_'):
                ref_seq_gene_id = variant.quibble.split(':')[0]
                tx_edit = variant.quibble.split(':')[1]
                gene_symbol = validator.db.get_gene_symbol_from_refseq_id(ref_seq_gene_id)
                if gene_symbol != 'none':
                    uta_symbol = validator.db.get_uta_symbol(gene_symbol)
                    available_transcripts = validator.hdp.get_tx_for_gene(uta_symbol)
                    select_from_these_transcripts = []
                    for tx in available_transcripts:
                        if 'NM_' in tx[3] or 'NR_' in tx[3] or 'ENST' in tx[3]:
                            if tx[3] not in select_from_these_transcripts:
                                select_from_these_transcripts.append(tx[3])
                    select_from_these_transcripts = '|'.join(select_from_these_transcripts)
                    if validator.select_transcripts != 'all' and validator.select_transcripts != 'raw':
                        variant.write = False
                        for transcript in list(select_transcripts_dict_plus_version.keys()):
                            variant.warnings = ['NG_:c.PositionVariation descriptions should not be used unless a '
                                                'transcript reference sequence has also been provided e.g. '
                                                'NG_(NM_):c.PositionVariation']
                            refreshed_description = ref_seq_gene_id + '(' + transcript + ')' + ':' + tx_edit
                            query = Variant(variant.original, quibble=refreshed_description,
                                            warnings=variant.warnings, primary_assembly=variant.primary_assembly,
                                            order=variant.order)

                            logger.info('NG_:c.PositionVariation descriptions should not be used unless a transcript '
                                        'reference sequence has also been provided e.g. NG_(NM_):c.PositionVariation. '
                                        'Resubmitting corrected version.')
                            validator.batch_list.append(query)
                            logger.info("Submitting new variant with format %s", refreshed_description)
                    else:
                        variant.warnings.append('A transcript reference sequence has not been provided e.g. '
                                                'NG_(NM_):c.PositionVariation. Re-submit ' + variant.quibble +
                                                ' but also specify transcripts from the following: ' +
                                                'select_transcripts=' + select_from_these_transcripts)
                        logger.warning('A transcript reference sequence has not been provided e.g. '
                                       'NG_(NM_):c.PositionVariation. Re-submit ' + variant.quibble + ' but also '
                                       'specify transcripts from the following: select_transcripts=' +
                                       select_from_these_transcripts)
                    skipvar = True
                else:
                    variant.warnings.append('A transcript reference sequence has not been provided e.g. '
                                            'NG_(NM_):c.PositionVariation')
                    logger.warning(
                        'A transcript reference sequence has not been provided e.g. NG_(NM_):c.PositionVariation')
                skipvar = True
            elif variant.quibble.startswith('NC_'):
                variant.warnings.append('A transcript reference sequence has not been provided e.g. '
                                        'NC_(NM_):c.PositionVariation. Unable to predict available transcripts '
                                        'because chromosomal position is not specified')
                logger.warning(
                    'A transcript reference sequence has not been provided e.g. NC_(NM_):c.PositionVariation. '
                    'Unable to predict available transcripts because chromosomal position is not specified')
                skipvar = True
        except Exception as e:
            logger.debug("Except passed, %s", e)

    logger.debug("Chromosomal/RefSeqGene reference catching complete")
    return skipvar


def vcf2hgvs_stage4(variant, validator):
    """
    VCF2HGVS conversion step 4 has two purposes
    1. VCF is frequently inappropriately converted into HGVS like descriptions
    such as GGGG>G which is actually a delins, del or ins. The function assigns
    the correct edit type
    2. Detects and extracts multiple ALT sequences into HGVS descriptions and
    automatically submits them for validation
    """
    skipvar = False
    not_sub = variant.quibble
    not_sub_find = re.compile(r"([GATCgatc]+)>([GATCgatc]+)")
    if '>' in not_sub and '(' not in not_sub and ')' not in not_sub:
        try:
            # If the length of either side of the substitution delimer (>) is >1
            matches = not_sub_find.search(not_sub)
            if len(matches.group(1)) > 1 or len(matches.group(2)) > 1 or \
                    ('>' in not_sub and ',' in not_sub):
                # Search for and remove range
                interval_range = re.compile(r"([0-9]+)_([0-9]+)")
                if interval_range.search(not_sub):
                    m = not_sub_find.search(not_sub)
                    start = m.group(1)
                    delete = m.group(2)
                    beginning_string, middle_string = not_sub.split(':')
                    middle_string = middle_string.split('_')[0]
                    end_string = start + '>' + delete
                    not_sub = beginning_string + ':' + middle_string + end_string
                # Split description
                ref_ac, _sep, remainder1 = not_sub.partition(':')
                ref_type, _sep, posedit = remainder1.partition('.')
                pos_ref, _sep, insert = posedit.partition('>')
                # If we have a list on inserts rather than 1 we need to resubmit
                # and abort! No need to continue with known over-loaded input
                if ',' in insert:
                    header = ref_ac + ':' + ref_type + '.' + pos_ref + '>'
                    alt_list = insert.split(',')
                    # Assemble and re-submit
                    for alt in alt_list:
                        variant.warnings = ['Multiple ALT sequences detected: '
                                            'auto-submitting all possible combinations']
                        variant.write = False
                        refreshed_description = header + alt
                        query = Variant(variant.original, quibble=refreshed_description,
                                        warnings=variant.warnings, primary_assembly=variant.primary_assembly,
                                        order=variant.order)
                        validator.batch_list.append(query)
                        logger.info('Multiple ALT sequences detected. Auto-submitting all possible combinations.')
                        logger.info("Submitting new variant with format %s", refreshed_description)
                    skipvar = True
                    return skipvar

                # Split ref from position using matches
                r = re.compile(r"([0-9]+)([GATCgatc]+)")
                m = r.search(pos_ref)
                delete = m.group(2)
                starts = posedit.split(delete)[0]
                hgvs_re_try = hgvs_delins_parts_to_hgvs_obj(
                        ref_ac,
                        ref_type,
                        starts,delete[0],insert,
                        offset_pos=True)
                hgvs_re_try.posedit.edit.ref = delete
                start_pos = str(hgvs_re_try.posedit.pos.start)
                if '-' in start_pos:
                    base, offset = start_pos.split('-')
                    new_offset = 0 - int(offset) + (len(delete))
                    hgvs_re_try.posedit.pos.end.base = int(base)
                    hgvs_re_try.posedit.pos.end.offset = int(new_offset) - 1
                    hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                            ref_ac,
                            ref_type,
                            hgvs_re_try.posedit.pos,delete,insert,
                            offset_pos=True)
                elif '+' in start_pos:
                    base, offset = start_pos.split('+')
                    end_pos = int(base) + (len(delete) - int(offset) - 1)
                    new_offset = 0 + int(offset) + (len(delete) - 1)
                    hgvs_re_try.posedit.pos.end.base = int(end_pos)
                    hgvs_re_try.posedit.pos.end.offset = int(new_offset)
                    hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                            ref_ac,
                            ref_type,
                            hgvs_re_try.posedit.pos,delete,insert,
                            offset_pos=True)
                else:
                    end_pos = int(start_pos) + (len(delete) - 1)
                    hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                            ref_ac,
                            ref_type,
                            start_pos,delete,insert,
                            end=end_pos,
                            offset_pos=True)

                # attempt to normalise output
                try:
                    not_delins = str(variant.hn.normalize(hgvs_not_delins))
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'Normalization of intronic variants is not supported' in error:
                        not_delins = str(hgvs_not_delins)
                    else:
                        variant.warnings.append(error)
                        logger.warning(str(e))
                        skipvar = True
                # Create warning
                automap = variant.quibble + ' automapped to ' + not_delins
                variant.warnings.append(automap)
                # Change input to normalized variant
                variant.quibble = not_delins
        except Exception as e:
            logger.debug("Except passed, %s", e)
    logger.debug("Completed VCF-HVGS step 4 for %s", variant.quibble)

    return skipvar

def convert_expanded_repeat(my_variant, validator):

    # Remove gene symbols from reference sequences
    if "(" in my_variant.quibble and ")" in my_variant.quibble:
        initial_formatting.remove_gene_symbol_from_ref(my_variant, validator)

    # Format expanded repeat syntax into a usable hgvs variant
    """
    Waiting for HGVS nomenclature changes
    """
    try:
        has_ex_repeat = expanded_repeats.convert_tandem(my_variant, validator, my_variant.primary_assembly,
                                                 "all")
    except expanded_repeats.RepeatSyntaxError as e:
        my_variant.warnings = [str(e)]
        return True
    except vvhgvs.exceptions.HGVSInvalidVariantError as e:
        my_variant.warnings = ["HgvsSyntaxError: " + str(e)]
        return True
    except vvhgvs.exceptions.HGVSDataNotAvailableError as e:
        if "invalid coordinates:" in str(e):
            my_variant.warnings = [(f"ExonBoundaryError: Stated position "
                                    f"does not correspond with an exon boundary for "
                                    f"transcript {my_variant.quibble.split(':')[0]}")]
            return True
    except Exception as e:
        # import traceback
        # traceback.print_exc()
        my_variant.warnings = ["ExpandedRepeatError: " + str(e)]
        return True

    try:
        has_ex_repeat
    except UnboundLocalError:
        return False
    else:
        if not has_ex_repeat:
            return False

    if my_variant.quibble != my_variant.expanded_repeat["variant"]:
        if re.search("\d+_", my_variant.quibble):
            my_variant.warnings.append(f"ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly"
                                       f" in the submitted description {my_variant.quibble}. The corrected format"
                                       f" would be {my_variant.expanded_repeat['variant'].split('[')[0]}"
                                       f"[int], where int requires you to update the number of repeats")
            return True
        else:
            my_variant.warnings.append(f"ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly"
                                       f" in the submitted description {my_variant.quibble}. The corrected description is "
                                       f"{my_variant.expanded_repeat['variant']}")
    ins_bases = (my_variant.expanded_repeat["repeat_sequence"] *
                 int(my_variant.expanded_repeat["copy_number"]))
    start_pos, _sep, end_pos = my_variant.expanded_repeat['position'].partition('_')
    repeat_to_delins = hgvs_delins_parts_to_hgvs_obj(
            my_variant.expanded_repeat['reference'],
            my_variant.expanded_repeat['prefix'],
            start_pos,
            '',
            ins_bases,
            end=end_pos)

    try:
        repeat_to_delins = my_variant.hn.normalize(repeat_to_delins)
    except vvhgvs.exceptions.HGVSUnsupportedOperationError:
        pass
    my_variant.quibble = repeat_to_delins #fn.valstr(repeat_to_delins)
    my_variant.warnings.append(f"ExpandedRepeatWarning: {my_variant.expanded_repeat['variant']} "
                               f"should only be used as an annotation for the core "
                               f"HGVS descriptions provided")
    return False

def indel_catching(variant, validator):
    """
    Warns that descriptions such as c.ins12 or g.del69 are not HGVS compliant
    Note - Updated to include dup28 examples
    Strips the trailing numbers and tries to parse the description into an
    hgvs object.
    If parses, provides a warning including links to the VarNomen web page, but
    continues validation
    If not, an error message is generated and the loop continues
    """
    edit_pass = re.compile(r'_\d+$')
    edit_fail = re.compile(r'\d+$')
    if edit_fail.search(variant.quibble):
        if not edit_pass.search(variant.quibble) and 'fs' not in variant.quibble:
            # Log 'dup in variant'
            dup_in_quibble = False
            error = 'Trailing digits are not permitted in HGVS variant descriptions'
            issue_link = 'http://varnomen.hgvs.org/recommendations/DNA/variant/'
            if 'dup' in variant.quibble:
                dup_in_quibble = True
                variant.quibble = variant.quibble.replace('dup', 'del')
            if '(' in variant.quibble:
                # Tackle compound variant descriptions NG or NC (NM_) i.e. correctly input NG/NC_(NM_):c.
                final_hgvs_convert(variant, validator)
                intronic_converter(variant, validator)
            else:
                final_hgvs_convert(variant, validator)
            try:
                validator.vr.validate(variant.quibble)
            except vvhgvs.exceptions.HGVSError as e:
                if 'Length implied by coordinates must equal ' \
                   'sequence deletion length' in str(e) and dup_in_quibble is True:
                    variant.warnings.append('Length implied by coordinates must equal sequence duplication length')
                else:
                    variant.warnings.append(str(e))
                variant.warnings.append(error)
                variant.warnings.append('Refer to ' + issue_link)
                logger.info(e)
                return True

            # Remove them so that the string SHOULD parse
            if dup_in_quibble is True:
                variant.quibble.posedit.edit.alt = variant.quibble.posedit.edit.ref + variant.quibble.posedit.edit.ref
            variant.warnings.append(error)
            variant.warnings.append('Refer to ' + issue_link)
            logger.info(error)
            return False

    logger.debug("Ins/Del reference catching complete for %s", variant.quibble)
    return False


def intronic_converter(variant, validator, skip_check=False, uncertain=False):
    """
    Fully HGVS compliant intronic variant descriptions take the format e.g
    NG_007400.1(NM_000088.3):c.589-1G>T. However, hgvs cannot parse and map
    these variant strings.
    This function:
    Removes the g. reference sequence
    NG_007400.1(NM_000088.3):c.589-1G>T ---> (NM_000088.3):c.589-1G>T
    Removes the parintheses
    (NM_000088.3):c.589-1G>T ---> NM_000088.3:c.589-1G>T
    hgvs can now parse the string into an hgvs variant object and manipulate it
    We now can parse in such variants but they still need fixing before mapping
    """
    compounder = re.compile(r'\((NM_|NR_|ENST)')
    compounder2 = re.compile(r'\(LRG_\d+t')
    if type(variant.quibble) is str:
        parsed = False
        acc_section, _sep, remainder = variant.quibble.partition(":")
    else:
        parsed = True
        acc_section = variant.quibble.ac

    if compounder.search(acc_section) or compounder2.search(acc_section):
        # Convert LRG transcript
        if compounder2.search(acc_section):
            lrg_transcript = acc_section.split("(")[1].replace(")", "")
            refseq_transcript = validator.db.get_refseq_transcript_id_from_lrg_transcript_id(lrg_transcript)
            if parsed:
                variant.quibble.ac = variant.quibble.ac.replace(lrg_transcript, refseq_transcript)
                acc_section = variant.quibble.ac
            else:
                variant.quibble = variant.quibble.replace(lrg_transcript, refseq_transcript)
                acc_section, _sep, _remain = variant.quibble.partition(":")
            variant.warnings.append(f"Reference sequence {lrg_transcript} updated to {refseq_transcript}")

        # Find pattern e.g. +0000 and assign to a variable
        if uncertain is True:
            genomic, _sep ,transcript = acc_section.partition('(')
            transcript = transcript.replace(')', '')
            if type(variant.quibble) is str:
                variant.quibble = f"{transcript}:{remainder}"
            else:
                variant.quibble.ac = transcript
            return variant
        else:
            genomic_ref = acc_section.split('(')[0]
        transy = re.search(r"((NM_|ENST|NR_).+)", acc_section)
        transy = transy.group(1)
        transy = transy.replace(')', '')

        # Add the edited variant for next stage error processing e.g. exon boundaries.

        if parsed:
            variant.quibble.ac = transy
            hgvs_transy = variant.quibble
        else:
            variant.quibble = variant.quibble.replace(acc_section,transy)
            try:
                hgvs_transy = validator.hp.parse_hgvs_variant(variant.quibble)
            except vvhgvs.exceptions.HGVSError:
                # Allele syntax caught here
                if "[" in variant.quibble and "]" in variant.quibble:
                    return genomic_ref
                else:
                    raise
        if skip_check is True:
            return genomic_ref
        else:
            # Check the specified base is correct
            hgvs_genomic = validator.nr_vm.c_to_g(variant.quibble, genomic_ref,
                                                  alt_aln_method=validator.alt_aln_method)
        try:
            validator.vr.validate(hgvs_genomic)
        except vvhgvs.exceptions.HGVSError as e:
            if 'Length implied by coordinates must equal sequence deletion length' in str(e) \
                    and not re.search(r'\d+$', variant.original):
                pass
            elif "does not agree with reference sequence" in str(e):
                previous_exception = e
                try:
                    check_g = validator.vm.t_to_g(hgvs_transy, hgvs_genomic.ac,
                                                  alt_aln_method=validator.alt_aln_method)
                    validator.vm.g_to_t(check_g, hgvs_transy.ac, alt_aln_method=validator.alt_aln_method)
                except vvhgvs.exceptions.HGVSError as e:
                    if "start or end or both are beyond the bounds of transcript record" in str(e):
                        if hgvs_transy.posedit.pos.start.offset != 0:
                            hgvs_transy.posedit.pos.start.offset = 1
                        if hgvs_transy.posedit.pos.end.offset != 0:
                            hgvs_transy.posedit.pos.end.offset = 1
                        remap_intronic(hgvs_transy, hgvs_genomic, variant, validator)
                        variant.warnings.append(e)
                        raise
                try:
                    variant.hn.normalize(hgvs_transy)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                    if "Unsupported normalization of variants spanning the exon-intron boundary" in str(e):
                        pass
                    else:
                        raise previous_exception
                else:
                    raise
            else:
                raise e

        # Check re-mapping of intronic variants
        remap_intronic(hgvs_genomic, hgvs_transy, variant, validator)

    logger.debug("HVGS typesetting complete")

def remap_intronic(hgvs_transy, hgvs_genomic, variant, validator):
    # Check re-mapping of intronic variants
    try:
        if hgvs_transy.posedit.pos.start.offset != 0 or hgvs_transy.posedit.pos.end.offset != 0:
            try:
                check_intronic_mapping = validator.nr_vm.g_to_t(hgvs_genomic, hgvs_transy.ac,
                                                                alt_aln_method=validator.alt_aln_method)
                if check_intronic_mapping.posedit.pos == hgvs_transy.posedit.pos:
                    pass
                else:
                    variant.warnings.append(f'ExonBoundaryError: {hgvs_transy.posedit.pos} does not match the exon '
                                            f'boundaries for the alignment of {hgvs_transy.ac} to {hgvs_genomic.ac}')
            except vvhgvs.exceptions.HGVSError as e:
                pass
    except AttributeError:
        pass

def allele_parser(variant, validation, validator):
    """
    HGVS allele string parsing function Occurance #1
    Takes a single HGVS allele description and separates each allele into a
    list of HGVS variants. The variants are then automatically submitted for
    validation.
    Note: In this context, it is inappropriate to validate descriptions
    containing intronic variant descriptions. In such instances, allele
    descriptions should be re-submitted by the user at the gene or genome level
    """
    caution = ''

    if (re.search(r':[gcnr].\[', variant.quibble) and ';' in variant.quibble) or (
            re.search(r':[gcrn].\d+\[', variant.quibble) and ';' in variant.quibble) or ('(;)'
                                                                                         in variant.quibble):
        # Edit compound descriptions
        genomic_ref = intronic_converter(variant, validator, skip_check=True)
        if genomic_ref is None:
            if re.match(r'NC_', variant.quibble):
                genomic_reference = variant.quibble.split(':')[0]
            else:
                genomic_reference = False
        elif 'NC_' in genomic_ref or 'NG_' in genomic_ref:
            genomic_reference = genomic_ref
        else:
            genomic_reference = False

        # handle LRG inputs
        if re.match(r'^LRG', variant.quibble):
            if re.match(r'^LRG\d+', variant.quibble):
                string, remainder = variant.quibble.split(':')
                reference = string.replace('LRG', 'LRG_')
                variant.quibble = reference + ':' + remainder
                caution = string + ' updated to ' + reference
            if not re.match(r'^LRG_\d+', variant.quibble):
                pass
            elif re.match(r'^LRG_\d+:g.', variant.quibble) or re.match(r'^LRG_\d+:p.',
                                                                       variant.quibble) \
                    or re.match(r'^LRG_\d+:c.', variant.quibble) or re.match(r'^LRG_\d+:n.',
                                                                             variant.quibble):
                lrg_reference, variation = variant.quibble.split(':')
                refseqgene_reference = validation.db.get_refseq_id_from_lrg_id(lrg_reference)
                if refseqgene_reference != 'none':
                    variant.quibble = refseqgene_reference + ':' + variation
                    if caution == '':
                        caution = lrg_reference + ':' + variation + ' automapped to ' + \
                                  refseqgene_reference + ':' + variation
                    else:
                        caution = caution + ': ' + lrg_reference + ':' + variation + ' automapped to equivalent ' \
                                                                                     'RefSeq record' + \
                                  refseqgene_reference + ':' + variation
                    variant.warnings.append(caution)
                    logger.info(caution)
            elif re.match(r'^LRG_\d+t\d+:c.', variant.quibble) or re.match(r'^LRG_\d+t\d+:n.',
                                                                           variant.quibble) or \
                    re.match(r'^LRG_\d+t\d+:p.', variant.quibble) or re.match(r'^LRG_\d+t\d+:g.',
                                                                              variant.quibble):
                lrg_reference, variation = variant.quibble.split(':')
                refseqtranscript_reference = validation.db.get_refseq_transcript_id_from_lrg_transcript_id(
                    lrg_reference)
                if refseqtranscript_reference != 'none':
                    variant.quibble = refseqtranscript_reference + ':' + variation
                    if caution == '':
                        caution = lrg_reference + ':' + variation + ' automapped to equivalent RefSeq record ' + \
                                  refseqtranscript_reference + ':' + variation
                    else:
                        caution = caution + ': ' + lrg_reference + ':' + variation + ' automapped to equivalent ' \
                                                                                     'RefSeq record' + \
                                  refseqtranscript_reference + ':' + variation
                    variant.warnings.append(caution)
                    logger.info(caution)
            else:
                pass
        try:
            # Submit to allele extraction function
            try:
                alleles = validation.hgvs_alleles(variant, genomic_reference)
            except fn.alleleVariantError as e:
                # import traceback
                # traceback.print_exc()
                variant.warnings.append(str(e))
                logger.warning(str(e))
                return True
            except AlleleSyntaxError as e:
                variant.warnings.append(str(e))
                return True

            variant.warnings.append('The alleleic description is in the correct syntax and all possible variant '
                                    'descriptions have been extracted')
            variant.warnings.append('Each variant is validated independently and users must update the original '
                                    'description accordingly based on these validations')

            logger.info('Automap has extracted possible variant descriptions, resubmitting')
            for allele in alleles:
                query = Variant(variant.original, quibble=allele, warnings=variant.warnings, write=True,
                                primary_assembly=variant.primary_assembly, order=variant.order)
                validation.batch_list.append(query)
                logger.info("Submitting new variant with format %s", allele)
            variant.write = False
            return True
        except fn.alleleVariantError as e:
            if "Cannot validate sequence of an intronic variant" in str(e):
                variant.warnings.append('Intronic positions not supported for HGVS Allele descriptions')
                logger.warning('Intronic positions not supported for HGVS Allele descriptions')
                return True
            elif "No transcript definition for " in str(e):
                variant.warnings.append(str(e))
                logger.warning(str(e))
                return True
            else:
                raise fn.VariantValidatorError(str(e))
    logger.debug("HVGS String allele parsing pass 1 complete")
    return False


def lrg_to_refseq(variant, validator):
    """
    LRG and LRG_t reference sequence identifiers need to be replaced with
    equivalent RefSeq identifiers. The lookup data is stored in the
    VariantValidator  MySQL database
    Currently only used post obj conversion
    """
    caution = ''
    if variant.refsource == 'LRG':
        if re.match(r'^LRG\d+', variant.hgvs_formatted.ac):
            reference = variant.hgvs_formatted.ac.replace('LRG', 'LRG_')
            caution = variant.hgvs_formatted.ac + ' updated to ' + reference + ': '
            variant.hgvs_formatted.ac = reference
            variant.set_quibble(variant.hgvs_formatted)

        if re.match(r'^LRG_\d+t\d+$', variant.quibble.ac):
            lrg_reference = variant.quibble.ac
            refseqtrans_reference = validator.db.get_refseq_transcript_id_from_lrg_transcript_id(lrg_reference)
            if refseqtrans_reference != 'none':
                old_var_str = str(variant.hgvs_formatted)
                variant.hgvs_formatted.ac = refseqtrans_reference
                variant.set_quibble(variant.hgvs_formatted)
                caution += old_var_str + ' automapped to equivalent RefSeq record ' \
                                                             '' + str(variant.hgvs_formatted)
                variant.warnings.append(caution)
                logger.info(caution)

        elif re.match(r'^LRG_\d+p\d+$', variant.quibble.ac):
            lrg_reference = variant.quibble.ac
            refseqprot_reference = validator.db.get_refseq_protein_id_from_lrg_protein_id(lrg_reference)
            if refseqprot_reference != 'none':
                old_var_str = str(variant.hgvs_formatted)
                variant.hgvs_formatted.ac = refseqprot_reference
                variant.set_quibble(variant.hgvs_formatted)
                caution +=  old_var_str + ' automapped to equivalent RefSeq record ' \
                                                             '' + str(variant.hgvs_formatted)
                variant.warnings.append(caution)
                logger.info(caution)

        elif re.match(r'^LRG_\d+$', variant.quibble.ac):
            lrg_reference = variant.quibble.ac
            refseqgene_reference = validator.db.get_refseq_id_from_lrg_id(lrg_reference)
            if refseqgene_reference != 'none':
                old_var_str = str(variant.hgvs_formatted)
                variant.hgvs_formatted.ac = refseqgene_reference
                variant.set_quibble(variant.hgvs_formatted)
                caution +=  old_var_str + ' automapped to equivalent RefSeq record ' \
                                                             '' + str(variant.hgvs_formatted)
                variant.warnings.append(caution)
                logger.info(caution)


def mitochondrial(variant, validator):
    """Will check if variant is mitochondrial and if so it will reformat the type to 'm' and save a value to the variant
    hgvs_genomic attribute"""

    if variant.reftype == ':m.' or variant.hgvs_formatted.ac == 'NC_012920.1' or \
            variant.hgvs_formatted.ac == 'NC_001807.4':

        # set flag
        variant.output_type_flag = 'mitochondrial'

        # Ensure the correct reference sequence type is used, if not, warn the user
        hgvs_mito = copy.deepcopy(variant.hgvs_formatted)
        if hgvs_mito.type == 'g' and (hgvs_mito.ac == 'NC_012920.1' or hgvs_mito.ac == 'NC_001807.4'):
            hgvs_mito.type = 'm'
            if "NC_012920.1" in hgvs_mito.ac and "hg19" in variant.selected_assembly:
                wrn = "NC_012920.1 is not associated with genome build hg19, instead use genome build GRCh37"
                variant.warnings.append(wrn)
                return True
            elif "NC_001807.4" in hgvs_mito.ac and "GRCh37" in variant.selected_assembly:
                wrn = "NC_001807.4 is not associated with genome build GRCh37, instead use genome build hg19"
                variant.warnings.append(wrn)
                return True
            else:
                wrn = "The given reference sequence (%s) does not match the DNA type (g). For %s, please use (m). " \
                  "For g. variants, please use a linear genomic reference sequence" % (hgvs_mito.ac, hgvs_mito.ac)
            variant.warnings.append(wrn)

        # Validate the variant description
        try:
            validator.vr.validate(hgvs_mito)
        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.warning(error)
            return True
        except KeyError:
            error = 'Currently unable to validate ' + hgvs_mito.ac + ' sequence variation'
            variant.warnings.append(error)
            logger.warning(error)
            return True
        else:

            # Check for movement during normalization
            try:
                norm_check = variant.hn.normalize(variant.hgvs_formatted)
                if hgvs_mito.posedit.pos != norm_check.posedit.pos:
                    norm_check.type = "m"
                    error = "%s updated to %s" % (fn.valstr(hgvs_mito), fn.valstr(norm_check))
                    variant.warnings.append(error)
                    logger.warning(error)

            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)
                variant.warnings.append(error)
                logger.warning(error)
                return True

            # Are there any transcripts?
            rel_var = validator.relevant_transcripts(hgvs_mito, variant.evm, validator.alt_aln_method,
                                                     variant.reverse_normalizer, validator.select_transcripts)

            # Add a description of the reference sequence type and continue
            variant.hgvs_genomic = hgvs_mito
            if len(rel_var) == 0:
                variant.genomic_g = unset_hgvs_obj_ref(hgvs_mito)
                variant.description = 'Homo sapiens mitochondrion, complete genome'
                logger.info('Homo sapiens mitochondrion, complete genome')
                return True

    return False


def proteins(variant, validator):
    """Handle protein sequences"""
    if variant.reftype == ':p.':
        error = None
        hgvs_object = None
        # Try to validate the variant
        hgvs_object = variant.hgvs_formatted
        try:
            validator.vr.validate(hgvs_object)

        except AttributeError as e:

            if "AARefAlt' object has no attribute 'ref_n'" in str(e):
                start_pos = hgvs_object.posedit.pos.start.pos
                end_pos = hgvs_object.posedit.pos.end.pos
                posedit = hgvs_object.posedit
                posedit = str(posedit).split(str(hgvs_object.posedit.edit))[0]
                posedit = posedit.replace("(", "").replace(")", "")

                if "_" in posedit:
                    start_edit, end_edit = posedit.split("_")
                    if "(" in start_edit:
                        start_edit = start_edit.replace("(", "")
                    start_aa = str(start_edit).replace(str(start_pos), "")
                    end_aa = str(end_edit).replace(str(end_pos), "")
                    start_aa_sl = fn.three_to_one(start_aa)
                    end_aa_sl = fn.three_to_one(end_aa)
                    check_start_aa = validator.sf.fetch_seq(hgvs_object.ac,
                                                            start_i=start_pos - 1,
                                                            end_i=start_pos)
                    check_end_aa = validator.sf.fetch_seq(hgvs_object.ac,
                                                          start_i=end_pos - 1,
                                                          end_i=end_pos)
                else:
                    start_edit = posedit
                    end_edit = posedit
                    start_aa = str(start_edit).replace(str(start_pos), "")
                    end_aa = str(end_edit).replace(str(end_pos), "")
                    start_aa_sl = fn.three_to_one(start_aa)
                    end_aa_sl = fn.three_to_one(end_aa)
                    check_start_aa = validator.sf.fetch_seq(hgvs_object.ac,
                                                            start_i=start_pos - 1,
                                                            end_i=start_pos)
                    check_end_aa = validator.sf.fetch_seq(hgvs_object.ac,
                                                          start_i=end_pos - 1,
                                                          end_i=end_pos)

                if start_aa_sl != check_start_aa and end_aa_sl != check_end_aa:

                    e1 = "The amino acid at position %s of %s is %s not %s" % (start_pos,
                                                                               hgvs_object.ac,
                                                                               check_start_aa,
                                                                               start_aa_sl)
                    e2 = "The amino acid at position %s of %s is %s not %s" % (end_pos,
                                                                               hgvs_object.ac,
                                                                               check_end_aa,
                                                                               end_aa_sl)

                    variant.warnings.extend([e1, e2])

                elif start_aa_sl != check_start_aa:
                    e1 = "The amino acid at position %s of %s is %s not %s" % (start_pos,
                                                                               hgvs_object.ac,
                                                                               check_start_aa,
                                                                               start_aa_sl)
                    variant.warnings.extend([e1])

                elif end_aa_sl != check_end_aa:
                    e1 = "The amino acid at position %s of %s is %s not %s" % (end_pos,
                                                                               hgvs_object.ac,
                                                                               check_end_aa,
                                                                               end_aa_sl)
                    variant.warnings.extend([e1])

                else:
                    error = str(
                        hgvs_object) + ' is HGVS compliant and contains a valid reference amino acid description'
                    reason = 'Protein level variant descriptions are not fully supported due to redundancy' \
                             ' in the genetic code'
                    variant.warnings.extend([reason, error])
                    logger.warning(reason + ": " + error)
                    variant.protein = hgvs_object
                    return True

        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
        if error:
            variant.warnings.append(error)
            logger.warning(error)
            return True
        else:
            try:
                validator.vr.validate(hgvs_object)
            except AttributeError as e:

                if "AARefAlt' object has no attribute 'ref_n'" in str(e):
                    start_pos = hgvs_object.posedit.pos.start.pos
                    end_pos = hgvs_object.posedit.pos.end.pos
                    posedit = hgvs_object.posedit
                    posedit = str(posedit).split(str(hgvs_object.posedit.edit))[0]
                    if "_" in posedit:
                        start_edit, end_edit = posedit.split("_")
                        start_aa = str(start_edit).replace(str(start_pos), "")
                        end_aa = str(end_edit).replace(str(end_pos), "")
                        start_aa_sl = fn.three_to_one(start_aa)
                        end_aa_sl = fn.three_to_one(end_aa)

                        check_start_aa = validator.sf.fetch_seq(hgvs_object.ac,
                                                                start_i=start_pos - 1,
                                                                end_i=start_pos)
                        check_end_aa = validator.sf.fetch_seq(hgvs_object.ac,
                                                              start_i=end_pos - 1,
                                                              end_i=end_pos)
                    else:
                        start_edit = posedit
                        end_edit = posedit
                        start_aa = str(start_edit).replace(str(start_pos), "")
                        end_aa = str(end_edit).replace(str(end_pos), "")
                        start_aa_sl = fn.three_to_one(start_aa)
                        end_aa_sl = fn.three_to_one(end_aa)
                        check_start_aa = validator.sf.fetch_seq(hgvs_object.ac,
                                                                start_i=start_pos - 1,
                                                                end_i=start_pos)
                        check_end_aa = validator.sf.fetch_seq(hgvs_object.ac,
                                                              start_i=end_pos - 1,
                                                              end_i=end_pos)

                    if start_aa_sl != check_start_aa and end_aa_sl != check_end_aa:
                        e1 = "The amino acid at position %s of %s is %s not %s" % (start_pos,
                                                                                   hgvs_object.ac,
                                                                                   check_start_aa,
                                                                                   start_aa_sl)
                        e2 = "The amino acid at position %s of %s is %s not %s" % (end_pos,
                                                                                   hgvs_object.ac,
                                                                                   check_end_aa,
                                                                                   end_aa_sl)
                        variant.warnings.extend([e1, e2])

                    elif start_aa_sl != check_start_aa:
                        e1 = "The amino acid at position %s of %s is %s not %s" % (start_pos,
                                                                                   hgvs_object.ac,
                                                                                   check_start_aa,
                                                                                   start_aa_sl)
                        variant.warnings.extend([e1])

                    elif end_aa_sl != check_end_aa:
                        e1 = "The amino acid at position %s of %s is %s not %s" % (end_pos,
                                                                                   hgvs_object.ac,
                                                                                   check_end_aa,
                                                                                   end_aa_sl)
                        variant.warnings.extend([e1])

                    else:
                        error = str(
                            hgvs_object) + ' is HGVS compliant and contains a valid reference amino acid description'
                        reason = 'Protein level variant descriptions are not fully supported due to redundancy' \
                                 ' in the genetic code'
                        variant.warnings.extend([reason, error])
                        logger.warning(reason + ": " + error)
                        variant.protein = hgvs_object
                        return True

                    return True

            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)
            else:
                error = str(
                    hgvs_object) + ' is HGVS compliant and contains a valid reference amino acid description'
            reason = 'Protein level variant descriptions are not fully supported due to redundancy' \
                     ' in the genetic code'
            variant.warnings.extend([reason, error])
            variant.protein = hgvs_object
            logger.warning(reason + ": " + error)
            return True
    return False


def rna(variant, validator):
    """
    convert r, into c.
    """
    if variant.reftype == ':r.' or ":r." in variant.original:
        if ":r.(" in str(variant.hgvs_formatted):
            if type(variant.hgvs_formatted) is str:
                strip_prediction = str(variant.hgvs_formatted).replace("(", "")
                strip_prediction = strip_prediction[:-1]
                hgvs_input = validator.hp.parse_hgvs_variant(strip_prediction)
            else:
                hgvs_input = variant.hgvs_formatted
                hgvs_input.posedit.pos.uncertain = False
                #hgvs_input.posedit.uncertain = False
        else:
            hgvs_input = variant.hgvs_formatted

        tx_info = validator.hdp.get_tx_identity_info(hgvs_input.ac)
        if tx_info[3] is None:
            error = "Invalid variant type for non-coding transcript. Instead use n."
            variant.warnings.append(error)
            logger.warning(str(error))
            return True
        # Change to coding variant
        variant.reftype = ':c.'
        # Change input to reflect!
        try:
            hgvs_c = validator.hgvs_r_to_c(hgvs_input)
        except vvhgvs.exceptions.HGVSDataNotAvailableError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.warning(str(error))
            return True
        variant.hgvs_formatted = hgvs_c

        # Create variant.rna_data dictionary
        rnd = VariantValidator.modules.rna_formatter.RnaDescriptions(validator.alt_aln_method,
                                                                     variant.primary_assembly,
                                                                     validator)
        try:
            rnd.check_syntax(str(variant.original), variant)
        except VariantValidator.modules.rna_formatter.RnaVariantSyntaxError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.warning(str(error))
            return True
        except vvhgvs.exceptions.HGVSParseError:
            try:
                rnd.check_syntax(str(variant.quibble), variant)
            except VariantValidator.modules.rna_formatter.RnaVariantSyntaxError as e:
                error = str(e)
                variant.warnings.append(error)
                logger.warning(str(error))
                return True

        # Add data to the variant object
        variant.rna_data = rnd.output_dict()

    return False


def uncertain_pos(variant, validator):
    """
    check for uncertain positions in the variant and return unsupported warning
    """
    try:
        to_check = variant.quibble
        posedit = to_check.split(':')[1]
        if '(' in posedit or ')' in posedit:
            if 'p.' in posedit or '(;)' in posedit or '(:)' in posedit:
                return False
            elif re.search("ins\(\d+\)$", posedit) or re.search("ins\(\d+_\d+\)$", posedit):
                return False
            elif ":r.(" in to_check:
                return False
            else:
                if ("[" in posedit or "]" in posedit) and not re.search("\[\d+\]", posedit):
                    return False
                try:
                    complex_descriptions.uncertain_positions(variant, validator)
                except complex_descriptions.IncompatibleTypeError:
                    use_checking.refseq_common_mistakes(variant)
                    # import traceback
                    # traceback.print_exc()
                    return True
                except complex_descriptions.InvalidRangeError as e:
                    variant.warnings.append(str(e))
                    # import traceback
                    # traceback.print_exc()
                    return True
                except Exception as e:
                    variant.warnings.append(str(e))
                    # import traceback
                    # traceback.print_exc()
                    return True

            return False
        else:
            return False
    except Exception:
        return False


# <LICENSE>
# Copyright (C) 2016-2025 VariantValidator Contributors
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
