import re
import hgvs
from .vvLogging import logger
from .variant import Variant
from . import vvChromosomes
from . import vvFunctions as fn


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
        variant.warnings = 'Not stating ALT bases is ambiguous because VCF specification 4.0 would treat ' + \
                           pre_input + ' as a deletion whereas VCF specification 4.1 onwards would treat ' + \
                           pre_input + ' as ALT = REF'
        variant.warnings += ': VariantValidator has output both alternatives'
        logger.resub('Not stating ALT bases is ambiguous because VCF specification 4.0 would treat ' +
                     pre_input + ' as a deletion whereas VCF specification 4.1 onwards would treat ' + pre_input +
                     ' as ALT = REF. Validator will output both alternatives.')
        variant.write = False
        input_A = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], 'del')
        input_B = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], vcf_elements[2])
        queryA = Variant(variant.original, quibble=input_A, warnings=variant.warnings,
                         primary_assembly=variant.primary_assembly, order=variant.order)
        queryB = Variant(variant.original, quibble=input_B, warnings=variant.warnings,
                         primary_assembly=variant.primary_assembly, order=variant.order)
        validator.batch_list.append(queryA)
        validator.batch_list.append(queryB)
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
            selected_assembly = in_list[0]
            variant.quibble = '-'.join(in_list[1:])
        pre_input = variant.quibble
        vcf_elements = pre_input.split('-')
        variant.quibble = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], vcf_elements[3])

    logger.trace("Completed VCF-HVGS step 1", variant)

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
    if re.search(r'\w+\:', variant.quibble) and not re.search(r'\w+\:[gcnmrp]\.', variant.quibble):
        if re.search(r'\w+\:[gcnmrp]', variant.quibble) and not re.search(r'\w+\:[gcnmrp]\.', variant.quibble):
            # Missing dot
            pass
        else:
            try:
                if 'GRCh37' in variant.quibble or 'hg19' in variant.quibble:
                    variant.primary_assembly = 'GRCh37'
                elif 'GRCh38' in variant.quibble or 'hg38' in variant.quibble:
                    variant.primary_assembly = 'GRCh38'
                input_list = variant.quibble.split(':')
                pos_ref_alt = str(input_list[1])
                positionAndEdit = input_list[1]
                if not re.match(r'N[CGTWMRP]_', variant.quibble) and not re.match(r'LRG_', variant.quibble):
                    chr_num = str(input_list[0])
                    chr_num = chr_num.upper().strip()
                    if re.match('CHR', chr_num):
                        chr_num = chr_num.replace('CHR', '')
                    # Use selected assembly
                    accession = vvChromosomes.to_accession(chr_num, validator.selected_assembly)
                    if accession is None:
                        variant.warnings += ': ' + chr_num + \
                                                 ' is not part of genome build ' + validator.selected_assembly
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
                        positionAndEdit = str(position) + ref + '>' + alt
                    elif 'ins' in variant.quibble:
                        pos = re.match(r'\d+', pos_ref_alt)
                        position = pos.group(0)
                        old_ref, old_alt = pos_ref_alt.split('>')
                        # old_ref = old_ref.replace(position, '')
                        position = int(position) - 1
                        required_base = validator.sf.fetch_seq(accession, start_i=position - 1, end_i=position)
                        ref = required_base
                        alt = required_base + old_alt
                        positionAndEdit = str(position) + ref + '>' + alt
                # Assign reference sequence type
                ref_type = validator.db.ref_type_assign(accession)
                if re.match('LRG_', accession):
                    if ref_type == ':g.':
                        accession = validator.db.get_RefSeqGeneID_from_lrgID(accession)
                    else:
                        accession = validator.db.get_RefSeqTranscriptID_from_lrgTranscriptID(accession)
                else:
                    accession = accession
                variant.quibble = str(accession) + ref_type + str(positionAndEdit)

            except:
                fn.exceptPass(variant)

    # Descriptions lacking the colon :
    if re.search(r'[gcnmrp]\.', variant.quibble) and not re.search(r':[gcnmrp]\.', variant.quibble):
        error = 'Unable to identify a colon (:) in the variant description %s. A colon is required in HGVS variant ' \
                'descriptions to separate the reference accession from the reference type i.e. <accession>:<type>. ' \
                'e.g. :c.' % (variant.quibble)
        variant.warnings += ': ' + error
        logger.warning(error)
        skipvar = True

    # Ambiguous chr reference
    logger.trace("Completed VCF-HVGS step 2", variant)

    return skipvar


def vcf2hgvs_stage3(variant, validator):
    """
    VCF2HGVS conversion step 3 is similar to step 2 but handles
    formats like Chr16:g.2099572TC>T which are provided by Alamut and other
    software
    """
    skipvar = False
    if re.search(r'\w+:[gcnmrp]\.', variant.quibble) and not re.match(r'N[CGTWMRP]_', variant.quibble):
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
                if 'GRCh37' in variant.quibble or 'hg19' in variant.quibble:
                    variant.primary_assembly = 'GRCh37'
                elif 'GRCh38' in variant.quibble or 'hg38' in variant.quibble:
                    variant.primary_assembly = 'GRCh38'
                input_list = variant.quibble.split(':')
                query_a_symbol = input_list[0]
                is_it_a_gene = validator.db.get_hgnc_symbol(query_a_symbol)
                if is_it_a_gene == 'none':
                    positionAndEdit = input_list[1]
                    chr_num = str(input_list[0])
                    chr_num = chr_num.upper().strip()
                    if re.match('CHR', chr_num):
                        chr_num = chr_num.replace('CHR', '')  # Use selected assembly
                    accession = vvChromosomes.to_accession(chr_num, validator.selected_assembly)
                    if accession is None:
                        variant.warnings += ': ' + chr_num + \
                                               ' is not part of genome build ' + validator.selected_assembly
                        skipvar = True
                    variant.quibble = str(accession) + ':' + str(positionAndEdit)
            except Exception:
                fn.exceptPass(variant)

    logger.trace("Completed VCF-HGVS step 3", variant)
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
    if re.search(r'\w+\:[cn]\.', variant.quibble):
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
                        variant.warnings = 'HGVS variant nomenclature does not allow the use of a gene symbol (' + \
                                              query_a_symbol + ') in place of a valid reference sequence'
                        refreshed_description = transcript + ':' + tx_edit
                        query = Variant(variant.original, quibble=refreshed_description,
                                        warnings=variant.warnings, primary_assembly=variant.primary_assembly,
                                        order=variant.order)
                        validator.batch_list.append(query)
                        logger.resub('HGVS variant nomenclature does not allow the use of a gene symbol (' + \
                                     query_a_symbol + ') in place of a valid reference sequence')
                else:
                    variant.warnings += ': ' + 'HGVS variant nomenclature does not allow the use of a gene symbol ('\
                                        + query_a_symbol + ') in place of a valid reference sequence: Re-submit ' + \
                                        variant.quibble + ' and specify transcripts from the following: ' + \
                                        'select_transcripts=' + select_from_these_transcripts
                    logger.warning('HGVS variant nomenclature does not allow the use of a gene symbol (' + \
                                   query_a_symbol + ') in place of a valid reference sequence: Re-submit ' +
                                   variant.quibble + ' and specify transcripts from the following: ' +
                                   'select_transcripts=' + select_from_these_transcripts)
                skipvar = True
        except:
            fn.exceptPass()
    logger.trace("Gene symbol reference catching complete", variant)
    return skipvar


def refseq_catch(variant, validator, select_transcripts_dict_plus_version):
    """
    Similar to the GENE_SYMBOL:c. n. types function, but spots RefSeqGene or
    Chromosomal reference sequence identifiers used in the context of c. variant
    descriptions
    """
    skipvar = False
    if re.search(r'\w+\:[cn]', variant.quibble):
        try:
            if variant.quibble.startswith('NG_'):
                refSeqGeneID = variant.quibble.split(':')[0]
                tx_edit = variant.quibble.split(':')[1]
                gene_symbol = validator.db.get_gene_symbol_from_refSeqGeneID(refSeqGeneID)
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
                            variant.warnings = 'NG_:c.PositionVariation descriptions should not be used unless a ' \
                                               'transcript reference sequence has also been provided e.g. ' \
                                               'NG_(NM_):c.PositionVariation'
                            refreshed_description = refSeqGeneID + '(' + transcript + ')' + ':' + tx_edit
                            query = Variant(variant.original, quibble=refreshed_description,
                                            warnings=variant.warnings, primary_assembly=variant.primary_assembly,
                                            order=variant.order)

                            logger.resub('NG_:c.PositionVariation descriptions should not be used unless a transcript '
                                         'reference sequence has also been provided e.g. NG_(NM_):c.PositionVariation. '
                                         'Resubmitting corrected version.')
                            validator.batch_list.append(query)
                    else:
                        variant.warnings += ': ' + 'A transcript reference sequence has not been provided e.g. ' \
                                                   'NG_(NM_):c.PositionVariation. Re-submit ' + variant.quibble + \
                                            ' but also specify transcripts from the following: ' + 'select_transcripts='\
                                            + select_from_these_transcripts
                        logger.warning('A transcript reference sequence has not been provided e.g. '
                                       'NG_(NM_):c.PositionVariation. Re-submit ' + variant.quibble + ' but also '
                                       'specify transcripts from the following: select_transcripts=' +
                                       select_from_these_transcripts)
                    skipvar = True
                else:
                    variant.warnings += ': ' + 'A transcript reference sequence has not been provided e.g. ' \
                                               'NG_(NM_):c.PositionVariation'
                    logger.warning(
                        'A transcript reference sequence has not been provided e.g. NG_(NM_):c.PositionVariation')
                skipvar = True
            elif variant.quibble.startswith('NC_'):
                variant.warnings += ': ' + 'A transcript reference sequence has not been provided e.g. ' \
                                           'NC_(NM_):c.PositionVariation. Unable to predict available transcripts ' \
                                           'because chromosomal position is not specified'
                logger.warning(
                    'A transcript reference sequence has not been provided e.g. NC_(NM_):c.PositionVariation. '
                    'Unable to predict available transcripts because chromosomal position is not specified')
                skipvar = True
        except:
            fn.exceptPass()

    logger.trace("Chromosomal/RefSeqGene reference catching complete", variant)

    return skipvar


def vcf2hgvs_stage4(variant, validator, hn):
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
    if not_sub_find.search(not_sub):
        try:
            # If the length of either side of the substitution delimer (>) is >1
            matches = not_sub_find.search(not_sub)
            if len(matches.group(1)) > 1 or len(matches.group(2)) > 1 or re.search(
                    r"([GATCgatc]+)>([GATCgatc]+),([GATCgatc]+)", variant.quibble):
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
                remainder = split_colon[1]
                split_dot = remainder.split('.')
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
                    start = m.group(1)
                    delete = m.group(2)
                    starts = posedit.split(delete)[0]
                    re_try = ref_ac + ':' + ref_type + '.' + starts + 'del' + delete[0] + 'ins' + insert
                    hgvs_re_try = validator.hp.parse_hgvs_variant(re_try)
                    hgvs_re_try.posedit.edit.ref = delete
                    start_pos = str(hgvs_re_try.posedit.pos.start)
                    if re.search(r'\-', start_pos):
                        base, offset = start_pos.split('-')
                        new_offset = 0 - int(offset) + (len(delete))
                        end_pos = int(base)
                        hgvs_re_try.posedit.pos.end.base = int(end_pos)
                        hgvs_re_try.posedit.pos.end.offset = int(new_offset) - 1
                        not_delins = ref_ac + ':' + ref_type + '.' + start_pos + '_' + str(
                            hgvs_re_try.posedit.pos.end) + 'del' + delete + 'ins' + insert
                    elif re.search(r'\+', start_pos):
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
                except:
                    fn.exceptPass()
                    not_delins = not_sub
                # Parse into hgvs object
                try:
                    hgvs_not_delins = validator.hp.parse_hgvs_variant(not_delins)
                except hgvs.exceptions.HGVSError as e:
                    # Sort out multiple ALTS from VCF inputs
                    if re.search(r"([GATCgatc]+)>([GATCgatc]+),([GATCgatc]+)", not_delins):
                        header, alts = not_delins.split('>')
                        # Split up the alts into a list
                        alt_list = alts.split(',')
                        # Assemble and re-submit
                        for alt in alt_list:
                            variant.warnings = 'Multiple ALT sequences detected: auto-submitting all possible combinations'
                            variant.write = False
                            refreshed_description = header + '>' + alt
                            query = Variant(variant.original, quibble=refreshed_description,
                                            warnings=variant.warnings, primary_assembly=variant.primary_assembly,
                                            order=variant.order)

                            validator.batch_list.append(query)
                            logger.resub(
                                'Multiple ALT sequences detected. Auto-submitting all possible combinations.')
                        skipvar = True
                    else:
                        error = str(e)
                        variant.warnings += ': ' + error
                        logger.warning(str(e))
                        skipvar = True

                try:
                    not_delins = str(hn.normalize(hgvs_not_delins))
                except hgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if re.search('Normalization of intronic variants is not supported', error):
                        not_delins = not_delins
                    else:
                        issue_link = ''
                        variant.warnings += ': ' + str(error)
                        logger.warning(str(e))
                        skipvar = True
                # Create warning
                caution = 'Variant description ' + variant.quibble + ' is not HGVS compliant'
                automap = variant.quibble + ' automapped to ' + not_delins
                variant.warnings += ': ' + automap
                # Change input to normalized variant
                variant.quibble = not_delins
        except:
            fn.exceptPass()
    logger.trace("Completed VCF-HVGS step 4", variant)
    return skipvar


def indel_catching(variant, validator):
    """
    Warns that descriptions such as c.ins12 or g.del69 are not HGVS compliant
    Strips the trailing numbers and tries to parse the description into an
    hgvs object.
    If parses, provides a warning including links to the VarNomen web page, but
    continues validation
    If not, an error message is generated and the loop continues
    """
    edit_pass = re.compile(r'_\d+$')
    edit_fail = re.compile(r'\d+$')
    if edit_fail.search(variant.quibble):
        if not edit_pass.search(variant.quibble):
            failed = variant.quibble
            # Catch the trailing digits
            digits = re.search(r"(\d+$)", failed)
            digits = digits.group(1)
            remove = str(digits) + 'end_anchor'
            failed = failed + 'end_anchor'
            failed = failed.replace(remove, '')

            # Remove them so that the string SHOULD parse
            try:
                hgvs_failed = validator.hp.parse_hgvs_variant(failed)
            except hgvs.exceptions.HGVSError as e:
                error = 'The syntax of the input variant description is invalid '
                if failed.endswith('ins'):
                    issue_link = 'http://varnomen.hgvs.org/recommendations/DNA/variant/insertion/'
                    error = error + ' please refer to ' + issue_link
                variant.warnings += error
                logger.warning(str(error) + " " + str(e))
                return True

            hgvs_failed.posedit.edit = str(hgvs_failed.posedit.edit).replace(digits, '')
            failed = str(hgvs_failed)
            automap = 'Non HGVS compliant variant description ' + variant.quibble + ' automapped to ' + failed
            variant.warnings += ': ' + automap
            logger.warning(automap)
            variant.quibble = failed

    logger.trace("Ins/Del reference catching complete", variant)
    return False


def intronic_converter(variant):
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
        transy = re.search(r"(NM_.+)", variant.quibble)
        transy = transy.group(1)
        transy = transy.replace(')', '')
        variant.quibble = transy
    logger.trace("HVGS typesetting complete", variant)

