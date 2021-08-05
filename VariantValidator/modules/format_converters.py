import re
import vvhgvs.exceptions
import copy
import logging
from .variant import Variant
from . import seq_data
from . import utils as fn

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

    # Uncertain positions
    toskip = uncertain_pos(variant)
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

    toskip = indel_catching(variant, validator)
    if toskip:
        return True

    # Tackle compound variant descriptions NG or NC (NM_) i.e. correctly input NG/NC_(NM_):c.
    intronic_converter(variant, validator)

    return False


def vcf2hgvs_stage1(variant, validator):
    """
    VCF2HGVS stage 1. converts chr-pos-ref-alt into chr:posRef>Alt
    The output format is a common mistake caused by inaccurate conversion of
    VCF variants into HGVS - hence the need for conversion step 2
    """
    skipvar = False

    if re.search(r'[-:]\d+[-:][GATC]+[-:][GATC]+', variant.quibble):
        variant.quibble = variant.quibble.replace(':', '-')
        # Extract primary_assembly if provided
        if re.match(r'GRCh3\d+-', variant.quibble) or re.match(r'hg\d+-', variant.quibble):
            in_list = variant.quibble.split('-')
            validator.selected_assembly = in_list[0]
            variant.quibble = '-'.join(in_list[1:])
        pre_input = variant.quibble
        vcf_elements = pre_input.split('-')
        variant.quibble = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], vcf_elements[3])
    elif re.search(r'[-:]\d+[-:][GATC]+[-:]', variant.quibble):
        variant.quibble = variant.quibble.replace(':', '-')
        # Extract primary_assembly if provided
        if re.match(r'GRCh3\d+-', variant.quibble) or re.match(r'hg\d+-', variant.quibble):
            in_list = variant.quibble.split('-')
            validator.selected_assembly = in_list[0]
            variant.quibble = '-'.join(in_list[1:])
        pre_input = variant.quibble
        vcf_elements = pre_input.split('-')
        variant.warnings = ['Not stating ALT bases is ambiguous because VCF specification 4.0 would treat ' +
                            pre_input + ' as a deletion whereas VCF specification 4.1 onwards would treat ' +
                            pre_input + ' as ALT = REF']
        variant.warnings.append('VariantValidator has output both alternatives')
        logger.info('Not stating ALT bases is ambiguous because VCF specification 4.0 would treat %s as a deletion '
                    'whereas VCF specification 4.1 onwards would treat %s as ALT = REF. Validator will output '
                    'both alternatives.', pre_input, pre_input)
        variant.write = False
        input_a = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], 'del')
        input_b = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], vcf_elements[2])
        query_a = Variant(variant.original, quibble=input_a, warnings=variant.warnings,
                          primary_assembly=variant.primary_assembly, order=variant.order)
        query_b = Variant(variant.original, quibble=input_b, warnings=variant.warnings,
                          primary_assembly=variant.primary_assembly, order=variant.order)
        validator.batch_list.append(query_a)
        validator.batch_list.append(query_b)
        logger.info("Submitting new variant with format %s", input_a)
        logger.info("Submitting new variant with format %s", input_b)
        skipvar = True
    elif re.search(r'[-:]\d+[-:][-:][GATC]+', variant.quibble) or \
            re.search(r'[-:]\d+[-:][.][-:][GATC]+', variant.quibble):
        variant.quibble = variant.quibble.replace(':', '-')
        if re.search('-.-', variant.quibble):
            variant.quibble = variant.quibble.replace('-.-', '-ins-')
        if re.search('--', variant.quibble):
            variant.quibble = variant.quibble.replace('--', '-ins-')
        # Extract primary_assembly if provided
        if re.match(r'GRCh3\d+-', variant.quibble) or re.match(r'hg\d+-', variant.quibble):
            in_list = variant.quibble.split('-')
            variant.quibble = '-'.join(in_list[1:])
        pre_input = variant.quibble
        vcf_elements = pre_input.split('-')
        variant.quibble = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], vcf_elements[3])

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
            (re.search(r'\w+:[gcnmrp]\.', variant.quibble) or re.search(r'\w+\(\w+\):[gcnmrp]\.', variant.quibble)):
        if re.search(r'\w+:[gcnmrp]', variant.quibble) and not re.search(r'\w+:[gcnmrp]\.', variant.quibble):
            # Missing dot
            pass
        else:
            try:
                if 'GRCh37' in variant.quibble or 'hg19' in variant.quibble:
                    variant.primary_assembly = 'GRCh37'
                    validator.selected_assembly = 'GRCh37'
                    variant.quibble.format_quibble()
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
                ref_type = validator.db.ref_type_assign(accession)
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

    # Descriptions lacking the colon :
    if re.search(r'[gcnmrp]\.', variant.quibble) and not re.search(r':[gcnmrp]\.', variant.quibble):
        error = 'Unable to identify a colon (:) in the variant description %s. A colon is required in HGVS variant ' \
                'descriptions to separate the reference accession from the reference type i.e. <accession>:<type>. ' \
                'e.g. :c.' % variant.quibble
        variant.warnings.append(error)
        logger.warning(error)
        skipvar = True

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
    if (re.search(r'\w+:[gcnmrp]\.', variant.quibble) or re.search(r'\w+\(\w+\):[gcnmrp]\.', variant.quibble)) \
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
    if re.search(r'\w+:[cn]\.', variant.quibble):
        try:
            pre_input = variant.quibble.split(':')
            query_a_symbol = pre_input[0]
            tx_edit = pre_input[1]
            is_it_a_gene = validator.db.get_hgnc_symbol(query_a_symbol)
            if is_it_a_gene != 'none':
                uta_symbol = validator.db.get_uta_symbol(is_it_a_gene)
                available_transcripts = validator.hdp.get_tx_for_gene(uta_symbol)
                select_from_these_transcripts = []
                for tx in available_transcripts:
                    if 'NM_' in tx[3] or 'NR_' in tx[3]:
                        if tx[3] not in select_from_these_transcripts:
                            select_from_these_transcripts.append(tx[3])
                select_from_these_transcripts = '|'.join(select_from_these_transcripts)
                if validator.select_transcripts != 'all':
                    variant.write = False
                    for transcript in list(select_transcripts_dict_plus_version.keys()):
                        variant.warnings = ['HGVS variant nomenclature does not allow the use of a gene symbol (' +
                                            query_a_symbol + ') in place of a valid reference sequence']
                        refreshed_description = transcript + ':' + tx_edit
                        query = Variant(variant.original, quibble=refreshed_description,
                                        warnings=variant.warnings, primary_assembly=variant.primary_assembly,
                                        order=variant.order)
                        validator.batch_list.append(query)
                        logger.info('HGVS variant nomenclature does not allow the use of a gene symbol (' +
                                    query_a_symbol + ') in place of a valid reference sequence')
                        logger.info("Submitting new variant with format %s", refreshed_description)
                else:
                    variant.warnings.append('HGVS variant nomenclature does not allow the use of a gene symbol ('
                                            + query_a_symbol + ') in place of a valid reference sequence: Re-submit ' +
                                            variant.quibble + ' and specify transcripts from the following: ' +
                                            'select_transcripts=' + select_from_these_transcripts)
                    logger.warning('HGVS variant nomenclature does not allow the use of a gene symbol (' +
                                   query_a_symbol + ') in place of a valid reference sequence: Re-submit ' +
                                   variant.quibble + ' and specify transcripts from the following: ' +
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
                        if 'NM_' in tx[3] or 'NR_' in tx[3]:
                            if tx[3] not in select_from_these_transcripts:
                                select_from_these_transcripts.append(tx[3])
                    select_from_these_transcripts = '|'.join(select_from_these_transcripts)
                    if validator.select_transcripts != 'all':
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
                    ('>' in variant.quibble and ',' in variant.quibble):
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
                split_colon = not_sub.split(':')
                ref_ac = split_colon[0]
                remainder1 = split_colon[1]
                split_dot = remainder1.split('.')
                ref_type = split_dot[0]
                remainder = split_dot[1]
                posedit = remainder
                split_greater = remainder.split('>')
                insert = split_greater[1]
                remainder = split_greater[0]
                # Split remainder using matches
                r = re.compile(r"([0-9]+)([GATCgatc]+)")
                try:
                    m = r.search(remainder)
                    delete = m.group(2)
                    starts = posedit.split(delete)[0]
                    re_try = ref_ac + ':' + ref_type + '.' + starts + 'del' + delete[0] + 'ins' + insert
                    hgvs_re_try = validator.hp.parse_hgvs_variant(re_try)
                    hgvs_re_try.posedit.edit.ref = delete
                    start_pos = str(hgvs_re_try.posedit.pos.start)
                    if '-' in start_pos:
                        base, offset = start_pos.split('-')
                        new_offset = 0 - int(offset) + (len(delete))
                        end_pos = int(base)
                        hgvs_re_try.posedit.pos.end.base = int(end_pos)
                        hgvs_re_try.posedit.pos.end.offset = int(new_offset) - 1
                        not_delins = ref_ac + ':' + ref_type + '.' + start_pos + '_' + str(
                            hgvs_re_try.posedit.pos.end) + 'del' + delete + 'ins' + insert
                    elif '+' in start_pos:
                        base, offset = start_pos.split('+')
                        end_pos = int(base) + (len(delete) - int(offset) - 1)
                        new_offset = 0 + int(offset) + (len(delete) - 1)
                        hgvs_re_try.posedit.pos.end.base = int(end_pos)
                        hgvs_re_try.posedit.pos.end.offset = int(new_offset)
                        not_delins = ref_ac + ':' + ref_type + '.' + start_pos + '_' + str(
                            hgvs_re_try.posedit.pos.end) + 'del' + delete + 'ins' + insert
                    else:
                        end_pos = int(start_pos) + (len(delete) - 1)
                        not_delins = ref_ac + ':' + ref_type + '.' + start_pos + '_' + str(
                            end_pos) + 'del' + delete + 'ins' + insert
                except Exception as e:
                    logger.debug("Except passed, %s", e)
                    not_delins = not_sub
                # Parse into hgvs object
                hgvs_not_delins = None
                try:
                    hgvs_not_delins = validator.hp.parse_hgvs_variant(not_delins)
                except vvhgvs.exceptions.HGVSError as e:
                    # Sort out multiple ALTS from VCF inputs
                    if '>' in not_delins and ',' in not_delins:
                        header, alts = not_delins.split('>')
                        # Split up the alts into a list
                        alt_list = alts.split(',')
                        # Assemble and re-submit
                        for alt in alt_list:
                            variant.warnings = ['Multiple ALT sequences detected: '
                                                'auto-submitting all possible combinations']
                            variant.write = False
                            refreshed_description = header + '>' + alt
                            query = Variant(variant.original, quibble=refreshed_description,
                                            warnings=variant.warnings, primary_assembly=variant.primary_assembly,
                                            order=variant.order)

                            validator.batch_list.append(query)
                            logger.info('Multiple ALT sequences detected. Auto-submitting all possible combinations.')
                            logger.info("Submitting new variant with format %s", refreshed_description)
                        skipvar = True
                    else:
                        error = str(e)
                        variant.warnings.append(error)
                        logger.warning(str(e))
                        skipvar = True

                try:
                    not_delins = str(variant.hn.normalize(hgvs_not_delins))
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'Normalization of intronic variants is not supported' in error:
                        not_delins = not_delins
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
            try:
                hgvs_quibble = validator.hp.parse_hgvs_variant(variant.quibble)
            except vvhgvs.exceptions.HGVSError:
                # Tackle compound variant descriptions NG or NC (NM_) i.e. correctly input NG/NC_(NM_):c.
                intronic_converter(variant, validator)
                hgvs_quibble = validator.hp.parse_hgvs_variant(variant.quibble)
            try:
                validator.vr.validate(hgvs_quibble)
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
                variant.quibble = str(hgvs_quibble).replace('del', 'dup')
            variant.warnings.append(error)
            variant.warnings.append('Refer to ' + issue_link)
            logger.info(error)
            return False

    logger.debug("Ins/Del reference catching complete for %s", variant.quibble)
    return False


def intronic_converter(variant, validator, skip_check=False):
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
    """
    compounder = re.compile(r'\(NM_')
    if compounder.search(variant.quibble):
        # Find pattern e.g. +0000 and assign to a variable
        genomic_ref = variant.quibble.split('(')[0]
        transy = re.search(r"(NM_.+)", variant.quibble)
        transy = transy.group(1)
        transy = transy.replace(')', '')
        # Add the edited variant for next stage error processing e.g. exon boundaries.
        variant.quibble = transy
        if skip_check is True:
            return genomic_ref
        else:
            # Check the specified base is correct
            hgvs_genomic = validator.nr_vm.c_to_g(validator.hp.parse_hgvs_variant(transy), genomic_ref)
        try:
            validator.vr.validate(hgvs_genomic)
        except vvhgvs.exceptions.HGVSError as e:
            if 'Length implied by coordinates must equal sequence deletion length' in str(e) \
                    and not re.search(r'\d+$', variant.quibble):
                pass
            else:
                validator.vr.validate(hgvs_genomic)

    logger.debug("HVGS typesetting complete")


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
            re.search(r':[gcrn].\d+\[', variant.quibble) and ';' in variant.quibble) or ('(;)' in variant.quibble):

        # Edit compound descriptions
        genomic_ref = intronic_converter(variant, validator, skip_check=True)
        if genomic_ref is None:
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
            elif re.match(r'^LRG_\d+:g.', variant.quibble) or re.match(r'^LRG_\d+:p.', variant.quibble) \
                    or re.match(r'^LRG_\d+:c.', variant.quibble) or re.match(r'^LRG_\d+:n.', variant.quibble):
                lrg_reference, variation = variant.quibble.split(':')
                refseqgene_reference = validation.db.get_refseq_id_from_lrg_id(lrg_reference)
                if refseqgene_reference != 'none':
                    variant.quibble = refseqgene_reference + ':' + variation
                    if caution == '':
                        caution = lrg_reference + ':' + variation + ' automapped to ' + \
                                  refseqgene_reference + ':' + variation
                    else:
                        caution = caution + ': ' + lrg_reference + ':' + variation + ' automapped to ' + \
                                  refseqgene_reference + ':' + variation
                    variant.warnings.append(caution)
                    logger.info(caution)
            elif re.match(r'^LRG_\d+t\d+:c.', variant.quibble) or re.match(r'^LRG_\d+t\d+:n.', variant.quibble) or \
                    re.match(r'^LRG_\d+t\d+:p.', variant.quibble) or re.match(r'^LRG_\d+t\d+:g.', variant.quibble):
                lrg_reference, variation = variant.quibble.split(':')
                refseqtranscript_reference = validation.db.get_refseq_transcript_id_from_lrg_transcript_id(
                    lrg_reference)
                if refseqtranscript_reference != 'none':
                    variant.quibble = refseqtranscript_reference + ':' + variation
                    if caution == '':
                        caution = lrg_reference + ':' + variation + ' automapped to ' + \
                                  refseqtranscript_reference + ':' + variation
                    else:
                        caution = caution + ': ' + lrg_reference + ':' + variation + ' automapped to ' + \
                                  refseqtranscript_reference + ':' + variation
                    variant.warnings.append(caution)
                    logger.info(caution)
            else:
                pass
        try:
            # Submit to allele extraction function
            try:
                alleles = validation.hgvs_alleles(variant.quibble, variant.hn, genomic_reference)
            except fn.alleleVariantError as e:
                variant.warnings.append(str(e))
                logger.warning(str(e))
                return True
            variant.warnings.append('Automap has extracted possible variant descriptions')
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
    """
    caution = ''
    if variant.refsource == 'LRG':
        if re.match(r'^LRG\d+', variant.hgvs_formatted.ac):
            reference = variant.hgvs_formatted.ac.replace('LRG', 'LRG_')
            caution = variant.hgvs_formatted.ac + ' updated to ' + reference + ': '
            variant.hgvs_formatted.ac = reference
            variant.set_quibble(str(variant.hgvs_formatted))

        if re.match(r'^LRG_\d+t\d+:', variant.quibble):
            lrg_reference, variation = variant.quibble.split(':')
            refseqtrans_reference = validator.db.get_refseq_transcript_id_from_lrg_transcript_id(lrg_reference)
            if refseqtrans_reference != 'none':
                variant.hgvs_formatted.ac = refseqtrans_reference
                variant.set_quibble(str(variant.hgvs_formatted))
                caution += lrg_reference + ':' + variation + ' automapped to ' + refseqtrans_reference + ':' + variation
                variant.warnings.append(caution)
                logger.info(caution)
        elif re.match(r'^LRG_\d+:', variant.quibble):
            lrg_reference, variation = variant.quibble.split(':')
            refseqgene_reference = validator.db.get_refseq_id_from_lrg_id(lrg_reference)
            if refseqgene_reference != 'none':
                variant.hgvs_formatted.ac = refseqgene_reference
                variant.set_quibble(str(variant.hgvs_formatted))
                caution += lrg_reference + ':' + variation + ' automapped to ' + refseqgene_reference + ':' + variation
                variant.warnings.append(caution)
                logger.info(caution)


def mitochondrial(variant, validator):
    """Will check if variant is mitochondrial and if so it will reformat the type to 'm' and save a value to the variant
    hgvs_genomic attribute"""

    if variant.reftype == ':m.' or variant.hgvs_formatted.ac == 'NC_012920.1' or \
            variant.hgvs_formatted.ac == 'NC_001807.4':

        # set flag
        variant.output_type_flag = 'mitochondrial'

        hgvs_mito = copy.deepcopy(variant.hgvs_formatted)
        if hgvs_mito.type == 'g' and (hgvs_mito.ac == 'NC_012920.1' or hgvs_mito.ac == 'NC_001807.4'):
            hgvs_mito.type = 'm'
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
            # Any transcripts?
            rel_var = validator.relevant_transcripts(hgvs_mito, variant.evm, validator.alt_aln_method,
                                                     variant.reverse_normalizer)
            variant.hgvs_genomic = hgvs_mito
            if len(rel_var) == 0:
                variant.genomic_g = fn.valstr(hgvs_mito)
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
        try:
            hgvs_object = validator.hp.parse_hgvs_variant(variant.hgvs_formatted)
        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
        try:
            validator.vr.validate(hgvs_object)
        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
        if error:
            variant.warnings.append(error)
            logger.warning(error)
            return True
        else:
            # Get accurate descriptions from the relevant databases
            # RefSeq databases
            if validator.alt_aln_method != 'genebuild':
                # Gene description  - requires GenBank search to get all the required info, i.e. transcript variant ID
                # Look for the accession in our database
                # Connect to database and send request
                # record = validator.entrez_efetch(db="nuccore", id=accession, rettype="gb", retmode="text")

                try:
                    validator.vr.validate(hgvs_object)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                else:
                    error = str(
                        hgvs_object) + ' is HGVS compliant and contains a valid reference amino acid description'
                reason = 'Protein level variant descriptions are not fully supported due to redundancy' \
                         ' in the genetic code'
                variant.warnings.extend([reason, error])
                variant.protein = str(hgvs_object)
                logger.warning(reason + ": " + error)
                return True
    return False


def rna(variant, validator):
    """
    convert r, into c.
    """
    if variant.reftype == ':r.':
        hgvs_input = validator.hp.parse_hgvs_variant(str(variant.hgvs_formatted))
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

    return False


def uncertain_pos(variant):
    """
    check for uncertain positions in the variant and return unsupported warning
    """
    try:
        to_check = variant.quibble
        posedit = to_check.split(':')[1]
        if '(' in posedit or ')' in posedit:
            if 'p.' in posedit or '[' in posedit or ']' in posedit or '(;)' in posedit or '(:)' in posedit:
                return False
            error = 'Uncertain positions are not currently supported'
            variant.warnings.append(error)
            logger.warning(str(error))
            return True
        else:
            return False
    except Exception:
        return False


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
