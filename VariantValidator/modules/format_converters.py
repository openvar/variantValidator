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
from . import hgvs_position_utils

logger = logging.getLogger(__name__)

class GeneSymbolFormatError(Exception):
    pass

def initial_format_conversions(variant, validator, select_transcripts_dict_plus_version, batch_list):

    # VCF type 1
    toskip = vcf2hgvs_stage1(variant, batch_list)
    if toskip:
        return True

    # API type non-HGVS
    # e.g. Chr16:2099572TC>T
    toskip = vcf2hgvs_stage2(variant, validator)
    if toskip:
        return True

    toskip = gene_symbol_catch(variant, validator, select_transcripts_dict_plus_version, batch_list)
    if toskip:
        return True

    # NG_:c. or NC_:c.
    toskip = refseq_catch(variant, validator, select_transcripts_dict_plus_version, batch_list)
    if toskip:
        return True

    # Find not_sub type in input e.g. GGGG>G
    toskip = vcf2hgvs_stage4(variant, batch_list)
    if toskip:
        return True

    # Extract variants from HGVS allele descriptions
    # http://varnomen.hgvs.org/recommendations/DNA/variant/alleles/
    logger.info("Try find allele syntax in %s", variant.quibble)
    toskip = allele_parser(variant, validator, validator, batch_list)
    if toskip:
        return True

    # Conversions
    # are not currently supported. The HGVS format for conversions
    # is rarely seen wrt genomic sequencing data and needs to be re-evaluated
    # so abort before hgvs object conversion to avoid errors & parsing overhead
    if 'con' in variant.quibble:
        variant.warnings.append('Conversions are no longer valid HGVS Sequence Variant Descriptions')
        logger.info('Conversions are no longer valid HGVS Sequence Variant Descriptions')
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
    logger.info("Try find Expanded Repeat syntax in %s", variant.quibble)
    if isinstance(variant.quibble, str): # not Uncertain
        toskip = convert_expanded_repeat(variant, validator)
        if toskip:
            return True

    # Catches del12/ins21 type variants, can usfully trigger on the outupt of the expanded repeat conversions
    # and can not currently handle uncertain positions
    if isinstance(variant.quibble, str): # not Uncertain
        toskip = indel_catching(variant, validator)
        if toskip:
            return True

    # Quibble should now be correctly formatted hgvs & work for object parsing
    # or else already have been parsed already
    if isinstance(variant.quibble, str):
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


def vcf2hgvs_stage1(variant, batch_list):
    """
    VCF2HGVS stage 1. converts chr-pos-ref-alt into chr:posRef>Alt
    The output format is a common mistake caused by inaccurate conversion of
    VCF variants into HGVS - hence the need for conversion step 2
    """
    skipvar = False
    vcf_data = re.split(r'[-:]', variant.quibble)

    if len(vcf_data) < 4:
        logger.debug("Completed VCF-HVGS step 1 for %s", variant.quibble)
        return False

    # VCF CNV descriptions can be directly converted to a HGVS like description
    alt_lower = vcf_data[3].lower()
    if vcf_data[2].isdigit() and ('del' in alt_lower or 'inv' in alt_lower):

        if not any(
                type_prefix in part.lower()
                for part in vcf_data
                for type_prefix in ('g.', 'a.', 't.', 'c.', 'n.', 'm.', 'o.')
        ):
            logger.info("CNV format identified in %s", variant.quibble)
            cnv_var = f"{vcf_data[0]}:{vcf_data[1]}_{vcf_data[2]}{alt_lower}"
            logger.info("CNV identified, and mapped to %s", cnv_var)
            variant.warnings.append(
                f"VcfConversionWarning: CNV identified, and mapped to {cnv_var}"
            )
            variant.quibble = cnv_var

    poss_genome = vcf_data[0].lower()
    if ('grch3' in poss_genome or 'hg' in poss_genome) and poss_genome[-1].isdigit():
        vcf_data = vcf_data[1:]
        variant.quibble = '-'.join(vcf_data)

        # TODO test assembly given against settings
        if len(vcf_data) < 4:
            variant.warnings.append(
                "Insufficient or incorrect VCF elements provided. "
                "Elements required are chr-pos-ref-alt"
            )
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
    if vcf_data[3] == '.':
        vcf_data[3] = ''

    if re.search(r'[^CGAT]', vcf_data[2]) or re.search(r'[^CGAT,]', vcf_data[3]):
        logger.debug("Completed VCF-HVGS step 1 for %s", variant.quibble)
        return False

    if vcf_data[2] and vcf_data[3]:
        variant.quibble = f'{vcf_data[0]}:{vcf_data[1]}{vcf_data[2]}>{vcf_data[3]}'

    elif vcf_data[2]:
        variant.warnings = [
            'Not stating ALT bases is ambiguous because VCF specification 4.0 would treat ' +
            variant.quibble + ' as a deletion whereas VCF specification 4.1 onwards would treat ' +
            variant.quibble + ' as ALT = REF'
        ]
        variant.warnings.append('VariantValidator has output both alternatives')

        logger.info(
            'Not stating ALT bases is ambiguous because VCF specification 4.0 would treat %s as a deletion '
            'whereas VCF specification 4.1 onwards would treat %s as ALT = REF. Validator will output '
            'both alternatives.',
            variant.quibble,
            variant.quibble
        )

        variant.write = False

        input_a = f'{vcf_data[0]}:{vcf_data[1]}{vcf_data[2]}>del'
        input_b = f'{vcf_data[0]}:{vcf_data[1]}{vcf_data[2]}>{vcf_data[2]}'

        query_a = Variant(
            variant.original,
            quibble=input_a,
            warnings=variant.warnings,
            primary_assembly=variant.primary_assembly,
            order=variant.order
        )
        query_b = Variant(
            variant.original,
            quibble=input_b,
            warnings=variant.warnings,
            primary_assembly=variant.primary_assembly,
            order=variant.order
        )

        batch_list.append(query_a)
        batch_list.append(query_b)
        logger.info("Submitting new variant with format %s", input_a)
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
    Now updated to also handle <Chr16>(hg38):g.2099572TC>T, or NM_(GRCh3*) (or
    worse NC_(GRCh3*), type variant descriptions.
    """
    skipvar = False
    seq_id_part, _sep, type_posedit = variant.quibble.partition(':')
    var_type, _sep, posedit = type_posedit.partition('.')

    # We used to target variant descriptions without types here, and then move on to user input like
    # <Chr16>:g. in the next stage, but we need to start by fixing the reference identifier first
    # either way.
    # first abort on un-handleable types i.e no ':' and bad 'type' with no '.'
    if not type_posedit:
        if re.search(r'[gcnmrp]\.', variant.quibble):
            error = (
                'Unable to identify a colon (:) in the variant description %s. A colon is required in HGVS variant '
                'descriptions to separate the reference accession from the reference type i.e. <accession>:<type>. '
                'e.g. :c.' % variant.quibble
            )
            variant.warnings.append(error)
            logger.info(error)
            return True
        else:
            return False  # do we want to full abort on this yet?

    if type_posedit and not posedit and type_posedit[:1].lower() in ['g', 'c', 'n', 'm', 'r', 'p']:
        error = (
            'Unable to identify a dot (.) in the variant description %s following the reference sequence '
            'type (g,c,n,r, or p). A dot is required in HGVS variant '
            'descriptions to separate the reference type from the variant position i.e. <accession>:<type>. '
            'e.g. :g.' % variant.quibble
        )
        variant.warnings.append(error)
        logger.info(error)
        skipvar = True

    if posedit and var_type in ['G', 'C', 'N', 'M', 'R', 'P', 'O']:
        error = (
            'Reference type incorrectly stated in the variant description %s '
            'Valid types are g,c,n,r, or p' % variant.quibble
        )
        variant.warnings.append(error)
        logger.info(error)
        var_type = var_type.lower()

    # now check for and remove GRC/hg genome builds
    specifed_ref = False
    upper_seq_id_part = seq_id_part.upper()

    if 'GRCH37' in upper_seq_id_part or 'HG19' in upper_seq_id_part:
        specifed_ref = 'GRCh37'
    elif 'GRCH38' in upper_seq_id_part or 'HG38' in upper_seq_id_part:
        specifed_ref = 'GRCh38'

    if specifed_ref:
        if validator.selected_assembly and validator.selected_assembly != specifed_ref:
            variant.warnings.append(
                specifed_ref + ' from ' + seq_id_part +
                ' does not match the selected genome build of ' +
                validator.selected_assembly
            )
            return True

        if variant.primary_assembly and variant.primary_assembly != specifed_ref:
            variant.warnings.append(
                specifed_ref + ' from ' + seq_id_part +
                ' does not match the selected genome build of ' +
                validator.selected_assembly
            )
            return True

        if not validator.selected_assembly:
            variant.selected_assembly = specifed_ref

        if not variant.primary_assembly:
            variant.primary_assembly = specifed_ref

        upper_seq_id_part = seq_id_part.upper()
        build_start = upper_seq_id_part.find(specifed_ref.upper())

        if build_start != -1:
            build_end = build_start + len(specifed_ref)

            if build_start > 0 and seq_id_part[build_start - 1] == '(':
                build_start -= 1
            if build_end < len(seq_id_part) and seq_id_part[build_end] == ')':
                build_end += 1

            seq_id_part = seq_id_part[:build_start] + seq_id_part[build_end:]

        upper_seq_id_part = seq_id_part.upper()

    # now fix case, all non LRG Ref name characters should be upper case
    # LRG_ has t as a special case but gene symbols are fine in upper case too
    if 'LRG' in upper_seq_id_part:
        if 'l' in seq_id_part or 'r' in seq_id_part or 'g' in seq_id_part:
            seq_id_part = seq_id_part.replace('l', 'L')
            seq_id_part = seq_id_part.replace('r', 'R')
            seq_id_part = seq_id_part.replace('g', 'G')
    else:
        seq_id_part = upper_seq_id_part

    # Handle LRG references separately because untyped LRG transcripts and
    # genomic LRG references need to be mapped to different RefSeq accessions.
    if seq_id_part.startswith('LRG_'):
        # For now preserve LRG_ in non VCF type inputs, this may want changing
        # later but that means test updates too.
        if var_type:
            accession = seq_id_part
        elif 't' in seq_id_part:
            accession = validator.db.get_refseq_transcript_id_from_lrg_transcript_id(
                seq_id_part
            )
        else:
            accession = validator.db.get_refseq_id_from_lrg_id(seq_id_part)

    # Recognised RefSeq and Ensembl accessions can be retained directly.
    elif seq_id_part.startswith(
            ('NC_', 'NG_', 'NT_', 'NW_', 'NM_', 'NR_', 'NP_', 'ENST')):
        accession = seq_id_part

    else:
        chr_in = False

        if seq_id_part.startswith('CHR'):
            seq_id_part = seq_id_part[3:]
            chr_in = True

        accession = seq_data.to_accession(
            seq_id_part,
            variant.primary_assembly
        )

        # we limit full erroring out to variants without specified type as they are presumed to be VCF type
        # and so not valid as gene symbols, but just skip otherwise this preserves the old behaviour.
        # Also hard fail if the var started with chr but none was found, may also want to hard fail for
        # non c/t types or at least for g
        if accession is None and (chr_in or not posedit):
            error = (
                seq_id_part + ' is not part of genome build ' +
                validator.selected_assembly
            )
            variant.warnings.append(error)
            logger.info(error)
            return True

        elif accession is None:  # non vcf type failure
            accession = seq_id_part

    # finally fix posedit if we have a VCF style input (i.e. no '.' split type)
    if not posedit:
        posedit = var_type

        if '>' in posedit and ('del' in posedit or 'ins' in posedit):
            pos = re.match(r'\d+', posedit)
            position_string = pos.group(0)
            old_ref, old_alt = posedit.split('>')

            position = int(position_string) - 1
            required_base = validator.sf.fetch_seq(
                accession,
                start_i=position - 1,
                end_i=position
            )

            if 'del' in posedit:
                old_ref = old_ref[len(position_string):]
                ref = required_base + old_ref
                alt = required_base
            else:
                ref = required_base
                alt = required_base + old_alt

            posedit = f'{position}{ref}>{alt}'

        # Assign reference sequence type if it is missing
        if accession in ["NC_012920.1", "NC_001807.4"]:
            var_type = "m"
        else:
            var_type = validator.db.ref_type_assign(accession)[1]

    # this among other things will force gene type variants to at least specify c/n etc.
    if not posedit:
        error = (
            'Unable to identify a dot (.) in the variant description %s following the reference sequence '
            'type (g,c,n,r, or p). A dot is required in HGVS variant '
            'descriptions to separate the reference type from the variant position i.e. <accession>:<type>. '
            'e.g. :g.' % variant.quibble
        )
        variant.warnings.append(error)
        logger.info(error)
        skipvar = True

    variant.quibble = accession + ':' + var_type + '.' + posedit
    logger.debug("Completed VCF-HVGS step 2 for %s", variant.quibble)

    return skipvar


def gene_symbol_catch(variant, validator, select_transcripts_dict_plus_version, batch_list):
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
    logger.info("Extracted query_a_symbol: %s", query_a_symbol)

    if (
            not query_a_symbol.startswith(
                ('NC_', 'NW_', 'NM_', 'NR_', 'NG_', 'LRG', 'ENST')
            )
            and tx_edit.startswith(('c.', 'n.'))
    ):
        try:
            if '(' in query_a_symbol and ')' in query_a_symbol:
                variant.warnings.append(
                    f"ReferenceSequenceError: {query_a_symbol} is an invalid "
                    f"reference sequence identifier"
                )
                raise GeneSymbolFormatError(
                    f"{query_a_symbol} contains parentheses, which is not allowed "
                    f"in gene symbols."
                )

            is_it_a_gene = validator.db.get_hgnc_symbol(query_a_symbol)

            if is_it_a_gene == 'none':
                warning = (
                    f"{query_a_symbol} identified as the reference sequence for "
                    f"{variant.quibble} and is not a valid, and is also not a valid "
                    f"gene symbol"
                )
                variant.warnings.append(warning)
                logger.info(warning)
                skipvar = True

            else:
                uta_symbol = validator.db.get_uta_symbol(is_it_a_gene)
                available_transcripts = validator.hdp.get_tx_for_gene(uta_symbol)
                select_from_these_transcripts = []

                for tx in available_transcripts:
                    transcript = tx[3]

                    if (
                            validator.alt_aln_method == 'splign'
                            and transcript.startswith(('NM_', 'NR_'))
                            and transcript not in select_from_these_transcripts
                    ):
                        select_from_these_transcripts.append(transcript)

                    elif (
                            validator.alt_aln_method == 'genebuild'
                            and transcript.startswith('ENST')
                            and transcript not in select_from_these_transcripts
                    ):
                        select_from_these_transcripts.append(transcript)

                if validator.select_transcripts not in ('all', 'raw'):
                    variant.write = False

                    for transcript in select_transcripts_dict_plus_version:
                        if transcript == 'mane':
                            for tx in select_from_these_transcripts:
                                annotation = validator.db.get_transcript_annotation(tx)
                                if (
                                        '"mane_select": true' in annotation
                                        or '"mane_plus_clinical": true' in annotation
                                ):
                                    transcript = tx
                                else:
                                    continue

                        elif transcript == 'mane_select':
                            for tx in select_from_these_transcripts:
                                annotation = validator.db.get_transcript_annotation(tx)
                                if '"mane_select": true' in annotation:
                                    transcript = tx
                                else:
                                    continue

                        variant.warnings.append(
                            'InvalidReferenceError: HGVS variant nomenclature does not '
                            'allow the use of a gene symbol (' +
                            query_a_symbol +
                            ') in place of a valid reference sequence'
                        )

                        refreshed_description = transcript + ':' + tx_edit

                        query = Variant(
                            variant.original,
                            quibble=refreshed_description,
                            warnings=variant.warnings,
                            primary_assembly=variant.primary_assembly,
                            order=variant.order
                        )
                        batch_list.append(query)

                        logger.info(
                            'HGVS variant nomenclature does not allow the use of a '
                            'gene symbol (%s) in place of a valid reference sequence',
                            query_a_symbol
                        )
                        logger.info(
                            "Submitting new variant with format %s",
                            refreshed_description
                        )

                else:
                    select_from_these_transcripts = '|'.join(
                        select_from_these_transcripts
                    )

                    variant.warnings.append(
                        'InvalidReferenceError: HGVS variant nomenclature does not allow '
                        'the use of a gene symbol (' +
                        query_a_symbol +
                        ') in place of a valid reference sequence: Re-submit ' +
                        str(variant.quibble) +
                        ' and specify transcripts from the following: ' +
                        'select_transcripts=' +
                        select_from_these_transcripts
                    )

                    logger.info(
                        'HGVS variant nomenclature does not allow the use of a gene '
                        'symbol (' +
                        query_a_symbol +
                        ') in place of a valid reference sequence: Re-submit ' +
                        str(variant.quibble) +
                        ' and specify transcripts from the following: '
                        'select_transcripts=' +
                        select_from_these_transcripts
                    )

                skipvar = True

        except Exception as e:
            logger.debug("Except passed, %s", e)

    logger.debug("Gene symbol reference catching complete")
    return skipvar


def refseq_catch(variant, validator, select_transcripts_dict_plus_version, batch_list):
    """
    Similar to the GENE_SYMBOL:c. n. types function, but spots RefSeqGene or
    Chromosomal reference sequence identifiers used in the context of c. variant
    descriptions
    """
    skipvar = False
    query_a_seq, _sep, tx_edit = variant.quibble.partition(':')

    if tx_edit.startswith(('c.', 'n.')):
        # remove, handle, and sometimes fix, some complex broken input types, e.g.
        # 'NC_000017.11(NC_000017.11(ENST00000357654.9)' warn and abort otherwise.
        # Since genes like ENST(GEN_ID):c. may be found tolerate junk, but
        # implicitly insist on GenomicReferenceID(TranscriptReferenceID) *not*
        # GenomicReferenceID(GeneID)(TranscriptReferenceID)
        curr_query_a_ref_seq = ''
        tx_seq_found = False

        while '(' in query_a_seq and not tx_seq_found:
            query_a_seq_test, _sep, query_a_tx_seq = query_a_seq.partition('(')

            if query_a_seq_test.startswith(('NM_', 'NR_', 'LRG', 'ENST')):
                # presume we have a gene symbol as the next component if we have one and abort
                # without complaining
                tx_seq_found = True

            elif query_a_seq_test.startswith(('NG_', 'NC_', 'NW_', 'NT_')):
                query_a_seq_test = query_a_seq_test.replace(')', '')

                if not curr_query_a_ref_seq or query_a_seq_test == curr_query_a_ref_seq:
                    if curr_query_a_ref_seq:
                        variant.quibble = query_a_seq + ':' + tx_edit

                    curr_query_a_ref_seq = query_a_seq_test
                    query_a_seq = query_a_tx_seq

                    if query_a_tx_seq.startswith(('NM_', 'NR_', 'LRG', 'ENST')):
                        tx_seq_found = True

                else:
                    variant.warnings.append(
                        'HgvsSyntaxError: '
                        'Multiple genomic reference sequences have been '
                        'provided in the same transcript variant description'
                        + variant.quibble + ' should at most have one genomic'
                        ' reference, starting with NG_, NC_, NW_, or NT_, '
                        'paired with the relevant transcript, in the form '
                        'GenomicReferenceID(TranscriptReferenceID). Please '
                        're-submit with your favored genomic reference or '
                        'submit each pair separately.'
                    )
                    return True

            elif query_a_tx_seq.startswith(('NM_', 'NR_', 'LRG', 'ENST')):
                tx_seq_found = True
                query_a_seq = query_a_tx_seq

            else:
                variant.warnings.append(
                    'InvalidReferenceError: '
                    'A transcript type variant description ( ' + variant.quibble + ' ) '
                    ' has been submitted with the reference specified in a manner '
                    'recognised as GenomicReferenceID(TranscriptReferenceID), but the '
                    'apparent Transcript Reference ID was not recognised as an expected'
                    ' type, i.e. either a RefSeq transcript, which should start with NM_'
                    ', or NR_ or an ENSEMBL transcript which should start with  ENST. '
                )
                return True

        try:
            if query_a_seq.startswith('NG_') and not tx_seq_found:
                ref_seq_gene_id = query_a_seq
                gene_symbol = validator.db.get_gene_symbol_from_refseq_id(ref_seq_gene_id)

                if gene_symbol != 'none':
                    uta_symbol = validator.db.get_uta_symbol(gene_symbol)
                    available_transcripts = validator.hdp.get_tx_for_gene(uta_symbol)
                    select_from_these_transcripts = []

                    for tx in available_transcripts:
                        transcript = tx[3]

                        if (
                                transcript.startswith(('NM_', 'NR_', 'ENST'))
                                and transcript not in select_from_these_transcripts
                        ):
                            select_from_these_transcripts.append(transcript)

                    if validator.select_transcripts not in ('all', 'raw'):
                        variant.write = False

                        for transcript in select_transcripts_dict_plus_version:
                            variant.warnings = [
                                'NG_:c.PositionVariation descriptions should not be used unless a '
                                'transcript reference sequence has also been provided e.g. '
                                'NG_(NM_):c.PositionVariation'
                            ]

                            refreshed_description = (
                                ref_seq_gene_id + '(' + transcript + ')' + ':' + tx_edit
                            )

                            query = Variant(
                                variant.original,
                                quibble=refreshed_description,
                                warnings=variant.warnings,
                                primary_assembly=variant.primary_assembly,
                                order=variant.order
                            )

                            logger.info(
                                'NG_:c.PositionVariation descriptions should not be used unless a transcript '
                                'reference sequence has also been provided e.g. NG_(NM_):c.PositionVariation. '
                                'Resubmitting corrected version.'
                            )
                            batch_list.append(query)
                            logger.info(
                                "Submitting new variant with format %s",
                                refreshed_description
                            )

                    else:
                        select_from_these_transcripts = '|'.join(
                            select_from_these_transcripts
                        )

                        variant.warnings.append(
                            'A transcript reference sequence has not been provided e.g. '
                            'NG_(NM_):c.PositionVariation. Re-submit ' + variant.quibble +
                            ' but also specify transcripts from the following: '
                            'select_transcripts=' + select_from_these_transcripts
                        )
                        logger.info(
                            'A transcript reference sequence has not been provided e.g. '
                            'NG_(NM_):c.PositionVariation. Re-submit ' + variant.quibble +
                            ' but also specify transcripts from the following: '
                            'select_transcripts=' + select_from_these_transcripts
                        )

                    skipvar = True

                else:
                    variant.warnings.append(
                        'A transcript reference sequence has not been provided e.g. '
                        'NG_(NM_):c.PositionVariation'
                    )
                    logger.info(
                        'A transcript reference sequence has not been provided e.g. '
                        'NG_(NM_):c.PositionVariation'
                    )
                    skipvar = True

            elif query_a_seq.startswith(('NC_', 'NW_', 'NT_')) and not tx_seq_found:
                variant.warnings.append(
                    'A transcript reference sequence has not been provided e.g. '
                    'NC_(NM_):c.PositionVariation. Unable to predict available transcripts '
                    'because chromosomal position is not specified'
                )
                logger.info(
                    'A transcript reference sequence has not been provided e.g. '
                    'NC_(NM_):c.PositionVariation. Unable to predict available transcripts '
                    'because chromosomal position is not specified'
                )
                skipvar = True

            elif curr_query_a_ref_seq and not tx_seq_found:
                genomic_prefix = curr_query_a_ref_seq.split('_', 1)[0]

                variant.warnings.append(
                    f'A transcript reference sequence has not been provided e.g. '
                    f'{genomic_prefix}_(NM_):c.PositionVariation'
                )
                logger.info(
                    f'A transcript reference sequence has not been provided e.g. '
                    f'{genomic_prefix}_(NM_):c.PositionVariation'
                )
                skipvar = True

        except Exception as e:
            logger.debug("Except passed, %s", e)

    logger.debug("Chromosomal/RefSeqGene reference catching complete")
    return skipvar

def vcf2hgvs_stage4(variant, batch_list):
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

    if '>' in not_sub and '(' not in not_sub and ')' not in not_sub:
        try:
            # Split description
            ref_ac, _sep, remainder = not_sub.partition(':')
            ref_type, _sep, posedit = remainder.partition('.')
            pos_ref, _sep, insert = posedit.partition('>')

            # Separate the position from the reference sequence
            ref_start = next(
                i for i, char in enumerate(pos_ref)
                if char.upper() in 'GATC'
            )
            position = pos_ref[:ref_start]
            delete = pos_ref[ref_start:]

            # If the length of either side of the substitution delimiter (>) is >1
            if len(delete) > 1 or len(insert) > 1:

                # Remove supplied range; end position is derived from REF length
                position = position.partition('_')[0]

                # If we have a list of inserts rather than 1 we need to resubmit
                # and abort! No need to continue with known over-loaded input
                if ',' in insert:
                    header = ref_ac + ':' + ref_type + '.' + pos_ref + '>'

                    # Assemble and re-submit
                    for alt in insert.split(','):
                        variant.warnings = [
                            'Multiple ALT sequences detected: '
                            'auto-submitting all possible combinations'
                        ]
                        variant.write = False
                        refreshed_description = header + alt

                        query = Variant(
                            variant.original,
                            quibble=refreshed_description,
                            warnings=variant.warnings,
                            primary_assembly=variant.primary_assembly,
                            order=variant.order
                        )
                        batch_list.append(query)

                        logger.info(
                            'Multiple ALT sequences detected. '
                            'Auto-submitting all possible combinations.'
                        )
                        logger.info(
                            "Submitting new variant with format %s",
                            refreshed_description
                        )

                    return True

                hgvs_re_try = hgvs_delins_parts_to_hgvs_obj(
                    ref_ac,
                    ref_type,
                    position,
                    delete[0],
                    insert,
                    offset_pos=True
                )
                hgvs_re_try.posedit.edit.ref = delete

                start = hgvs_re_try.posedit.pos.start
                offset = start.offset

                if hgvs_position_utils.offset_is_negative(
                        hgvs_re_try, check_start=True):
                    new_offset = offset + len(delete)
                    hgvs_re_try.posedit.pos.end.base = start.base
                    hgvs_re_try.posedit.pos.end.offset = new_offset - 1

                elif hgvs_position_utils.offset_is_positive(
                        hgvs_re_try, check_start=True):
                    end_pos = start.base + (len(delete) - offset - 1)
                    new_offset = offset + (len(delete) - 1)
                    hgvs_re_try.posedit.pos.end.base = end_pos
                    hgvs_re_try.posedit.pos.end.offset = new_offset

                else:
                    hgvs_re_try.posedit.pos.end.base = (
                        start.base + len(delete) - 1
                    )

                hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                    ref_ac,
                    ref_type,
                    hgvs_re_try.posedit.pos,
                    delete,
                    insert,
                    offset_pos=True
                )

                # attempt to normalise output
                try:
                    not_delins = str(variant.hn.normalize(hgvs_not_delins))
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'Normalization of intronic variants is not supported' in error:
                        not_delins = str(hgvs_not_delins)
                    else:
                        variant.warnings.append(error)
                        logger.info(error)
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
    if '(' in my_variant.quibble and ')' in my_variant.quibble:
        initial_formatting.remove_gene_symbol_from_ref(my_variant, validator)

    # Format expanded repeat syntax into a usable hgvs variant
    """
    Waiting for HGVS nomenclature changes
    """
    logger.info(
        "Checking my variant: %s for expanded repeats in format_converters",
        my_variant.quibble
    )

    try:
        has_ex_repeat = expanded_repeats.convert_tandem(
            my_variant,
            validator,
            my_variant.primary_assembly,
            "all"
        )
    except expanded_repeats.RepeatSyntaxError as e:
        my_variant.warnings = [str(e)]
        return True
    except vvhgvs.exceptions.HGVSInvalidVariantError as e:
        my_variant.warnings = ["HgvsSyntaxError: " + str(e)]
        return True
    except vvhgvs.exceptions.HGVSDataNotAvailableError as e:
        error = str(e)
        if "invalid coordinates:" in error:
            transcript = my_variant.quibble.partition(':')[0]
            my_variant.warnings = [
                f"ExonBoundaryError: Stated position "
                f"does not correspond with an exon boundary for "
                f"transcript {transcript}"
            ]
            return True
    except Exception as e:
        my_variant.warnings = ["ExpandedRepeatError: " + str(e)]
        return True

    if not has_ex_repeat:
        return False

    expanded_repeat_variant = my_variant.expanded_repeat["variant"]
    expanded_repeat_str = str(expanded_repeat_variant)

    if my_variant.quibble != expanded_repeat_str:
        if re.search(r"\d+_", my_variant.quibble):
            corrected_format = expanded_repeat_str.partition('[')[0]
            my_variant.warnings.append(
                f"ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly"
                f" in the submitted description {my_variant.quibble}. The corrected format"
                f" would be {corrected_format}"
                f"[int], where int requires you to update the number of repeats"
            )
            return True
        else:
            my_variant.warnings.append(
                f"ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly"
                f" in the submitted description {my_variant.quibble}. The corrected description is "
                f"{expanded_repeat_str}"
            )

    repeat_to_delins = copy.deepcopy(expanded_repeat_variant)
    repeat_to_delins.posedit.expanded_rep = False

    try:
        repeat_to_delins = my_variant.hn.normalize(repeat_to_delins)
    except vvhgvs.exceptions.HGVSUnsupportedOperationError:
        pass

    logger.info(
        "Expanded repeat to delins normalised: %s",
        repeat_to_delins
    )

    my_variant.quibble = repeat_to_delins
    my_variant.warnings.append(
        f"ExpandedRepeatWarning: {expanded_repeat_variant} "
        f"should only be used as an annotation for the core "
        f"HGVS descriptions provided"
    )

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
    if isinstance(variant.quibble, str):
        parsed = False
        acc_section, _sep, remainder = variant.quibble.partition(':')
    else:
        parsed = True
        acc_section = variant.quibble.ac

    pos_genomic, _sep, pos_transcript = acc_section.partition('(')

    if pos_transcript and (
            pos_transcript.startswith(('NM_', 'NR_', 'ENST')) or
            (
                pos_transcript.startswith('LRG_') and
                pos_transcript[4:5].isdigit() and
                't' in pos_transcript
            )
    ):
        # Convert LRG transcript
        if pos_transcript.startswith('LRG_'):
            lrg_transcript = pos_transcript.replace(')', '')
            refseq_transcript = (
                validator.db.get_refseq_transcript_id_from_lrg_transcript_id(
                    lrg_transcript
                )
            )

            if parsed:
                variant.quibble.ac = variant.quibble.ac.replace(
                    lrg_transcript,
                    refseq_transcript
                )
                acc_section = variant.quibble.ac
            else:
                variant.quibble = variant.quibble.replace(
                    lrg_transcript,
                    refseq_transcript
                )
                acc_section, _sep, _remain = variant.quibble.partition(':')

            variant.warnings.append(
                f"Reference sequence {lrg_transcript} updated to "
                f"{refseq_transcript}"
            )

        if uncertain:
            genomic, _sep, transcript = acc_section.partition('(')
            transcript = transcript.replace(')', '')

            if isinstance(variant.quibble, str):
                variant.quibble = f"{transcript}:{remainder}"
            else:
                variant.quibble.ac = transcript
                variant.quibble.rel_ac = genomic

            return variant

        genomic_ref, _sep, transcript = acc_section.partition('(')
        transy = transcript.replace(')', '')

        # Add the edited variant for next stage error processing e.g. exon boundaries.
        if parsed:
            variant.quibble.ac = transy
            variant.quibble.rel_ac = genomic_ref
            hgvs_transy = variant.quibble

        else:
            variant.quibble = variant.quibble.replace(acc_section, transy)

            try:
                hgvs_transy = validator.hp.parse_hgvs_variant(
                    variant.quibble
                )
            except vvhgvs.exceptions.HGVSError:
                # Allele syntax caught here
                if '[' in variant.quibble and ']' in variant.quibble:
                    return genomic_ref
                raise

        if skip_check:
            return genomic_ref

        # Check the specified base is correct
        hgvs_genomic = validator.nr_vm.c_to_g(
            variant.quibble,
            genomic_ref,
            alt_aln_method=validator.alt_aln_method
        )

        try:
            validator.vr.validate(hgvs_genomic)

        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)

            if (
                    'Length implied by coordinates must equal sequence deletion length'
                    in error
                    and not re.search(r'\d+$', variant.original)
            ):
                pass

            elif 'does not agree with reference sequence' in error:
                previous_exception = e

                try:
                    check_g = validator.vm.t_to_g(
                        hgvs_transy,
                        hgvs_genomic.ac,
                        alt_aln_method=validator.alt_aln_method
                    )
                    validator.vm.g_to_t(
                        check_g,
                        hgvs_transy.ac,
                        alt_aln_method=validator.alt_aln_method
                    )

                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)

                    if (
                            'start or end or both are beyond the bounds of '
                            'transcript record' in error
                    ):
                        if hgvs_position_utils.position_is_intronic(
                                hgvs_transy, check_start=True):
                            hgvs_transy.posedit.pos.start.offset = 1

                        if hgvs_position_utils.position_is_intronic(
                                hgvs_transy, check_end=True):
                            hgvs_transy.posedit.pos.end.offset = 1

                        remap_intronic(
                            hgvs_transy,
                            hgvs_genomic,
                            variant,
                            validator
                        )
                        variant.warnings.append(e)
                        raise

                try:
                    variant.hn.normalize(hgvs_transy)

                except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                    if (
                            'Unsupported normalization of variants spanning '
                            'the exon-intron boundary' in str(e)
                    ):
                        pass
                    else:
                        raise previous_exception
                else:
                    raise

            else:
                raise

        # Check re-mapping of intronic variants
        remap_intronic(
            hgvs_genomic,
            hgvs_transy,
            variant,
            validator
        )

    logger.debug("HVGS typesetting complete")

def remap_intronic(hgvs_transy, hgvs_genomic, variant, validator):
    # Check re-mapping of intronic variants
    try:
        if not hgvs_position_utils.position_is_intronic(
                hgvs_transy, check_start=True, check_end=True):
            return

        try:
            check_intronic_mapping = validator.nr_vm.g_to_t(
                hgvs_genomic,
                hgvs_transy.ac,
                alt_aln_method=validator.alt_aln_method
            )

            if check_intronic_mapping.posedit.pos != hgvs_transy.posedit.pos:
                variant.warnings.append(
                    f'ExonBoundaryError: {hgvs_transy.posedit.pos} does not match the exon '
                    f'boundaries for the alignment of {hgvs_transy.ac} to {hgvs_genomic.ac}'
                )

        except vvhgvs.exceptions.HGVSError:
            pass

    except AttributeError:
        pass

def allele_parser(variant, validation, validator, batch_list):
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
    ac_part, _sep, var = variant.quibble.partition(':')
    is_digit, _sep, end = var.partition('[')

    if (
            (var.startswith(('g.[', 'c.[', 'n.[', 'm.[', 'r.[')) and ';' in var)
            or (
                is_digit.startswith(('g.', 'c.', 'n.', 'm.', 'r.'))
                and ';' in end
                and is_digit[2:].isdigit()
            )
            or '(;)' in var
    ):
        # Edit compound descriptions
        genomic_ref = intronic_converter(
            variant,
            validator,
            skip_check=True
        )

        # Determine whether a genomic reference has already been supplied.
        #
        # Path 1:
        # intronic_converter() could not derive a genomic reference, so use the
        # submitted accession if it is already a genomic RefSeq accession.
        #
        # Path 2:
        # intronic_converter() successfully derived a genomic reference from the
        # submitted transcript/gene accession. Use it if it is a supported genomic
        # RefSeq accession.
        if genomic_ref is None:
            if ac_part.startswith(('NC_', 'NG_', 'NT_', 'NW_')):
                genomic_reference = ac_part
            else:
                genomic_reference = False

        elif genomic_ref.startswith(('NC_', 'NG_', 'NT_', 'NW_')):
            genomic_reference = genomic_ref

        else:
            genomic_reference = False

        # Handle LRG inputs
        if variant.quibble.startswith('LRG'):
            lrg_reference, separator, variation = variant.quibble.partition(':')

            # Correct LRG123 -> LRG_123
            if (
                    lrg_reference.startswith('LRG')
                    and not lrg_reference.startswith('LRG_')
                    and lrg_reference[3:].isdigit()
            ):
                old_reference = lrg_reference
                lrg_reference = 'LRG_' + lrg_reference[3:]
                variant.quibble = lrg_reference + separator + variation
                caution = old_reference + ' updated to ' + lrg_reference

            # Split the LRG identifier into genomic/transcript components.
            lrg_id = lrg_reference.removeprefix('LRG_')
            lrg_number, transcript_sep, transcript_number = lrg_id.partition('t')

            is_lrg = lrg_number.isdigit()
            is_lrg_transcript = (
                is_lrg
                and bool(transcript_sep)
                and transcript_number.isdigit()
            )

            if is_lrg_transcript and variation.startswith(
                    ('c.', 'n.', 'p.', 'g.')):
                refseqtranscript_reference = (
                    validation.db.get_refseq_transcript_id_from_lrg_transcript_id(
                        lrg_reference
                    )
                )

                if refseqtranscript_reference != 'none':
                    original_description = lrg_reference + ':' + variation
                    mapped_description = (
                        refseqtranscript_reference + ':' + variation
                    )
                    variant.quibble = mapped_description

                    if not caution:
                        caution = (
                            original_description +
                            ' automapped to equivalent RefSeq record ' +
                            mapped_description
                        )
                    else:
                        caution += (
                            ': ' + original_description +
                            ' automapped to equivalent RefSeq record ' +
                            mapped_description
                        )

                    variant.warnings.append(caution)
                    logger.info(caution)

            elif is_lrg and not transcript_sep and variation.startswith(
                    ('g.', 'p.', 'c.', 'n.')):
                refseqgene_reference = (
                    validation.db.get_refseq_id_from_lrg_id(
                        lrg_reference
                    )
                )

                if refseqgene_reference != 'none':
                    original_description = lrg_reference + ':' + variation
                    mapped_description = (
                        refseqgene_reference + ':' + variation
                    )
                    variant.quibble = mapped_description

                    if not caution:
                        caution = (
                            original_description +
                            ' automapped to ' +
                            mapped_description
                        )
                    else:
                        caution += (
                            ': ' + original_description +
                            ' automapped to equivalent RefSeq record ' +
                            mapped_description
                        )

                    variant.warnings.append(caution)
                    logger.info(caution)

        logger.info("Try allele extraction")

        try:
            # Submit to allele extraction function
            try:
                alleles = validation.hgvs_alleles(
                    variant,
                    genomic_reference
                )

            except fn.alleleVariantError as e:
                error = str(e)
                variant.warnings.append(error)
                logger.info(error)
                return True

            except AlleleSyntaxError as e:
                variant.warnings.append(str(e))
                return True

            variant.warnings.append(
                'The alleleic description is in the correct syntax and all possible '
                'variant descriptions have been extracted'
            )
            variant.warnings.append(
                'Each variant is validated independently and users must update the '
                'original description accordingly based on these validations'
            )

            logger.info(
                'Automap has extracted possible variant descriptions, resubmitting'
            )

            for allele in alleles:
                query = Variant(
                    variant.original,
                    quibble=allele,
                    warnings=variant.warnings,
                    write=True,
                    primary_assembly=variant.primary_assembly,
                    order=variant.order
                )
                batch_list.append(query)
                logger.info(
                    "Submitting new variant with format %s",
                    allele
                )

            variant.write = False
            return True

        except fn.alleleVariantError as e:
            error = str(e)

            if "Cannot validate sequence of an intronic variant" in error:
                variant.warnings.append(
                    'Intronic positions not supported for HGVS Allele descriptions'
                )
                logger.info(
                    'Intronic positions not supported for HGVS Allele descriptions'
                )
                return True

            elif "No transcript definition for " in error:
                variant.warnings.append(error)
                logger.info(error)
                return True

            else:
                raise fn.VariantValidatorError(error)

    logger.debug("HVGS String allele parsing pass 1 complete")
    return False


def lrg_to_refseq(variant, validator):
    """
    LRG and LRG_t reference sequence identifiers need to be replaced with
    equivalent RefSeq identifiers. The lookup data is stored in the
    VariantValidator MySQL database.
    Currently only used post obj conversion.
    """
    caution = ''

    if variant.refsource != 'LRG':
        return

    hgvs_formatted = variant.hgvs_formatted

    # Correct LRG123 to LRG_123
    if hgvs_formatted.ac.startswith('LRG') and hgvs_formatted.ac[3:4].isdigit():
        old_reference = hgvs_formatted.ac
        reference = 'LRG_' + old_reference[3:]
        caution = old_reference + ' updated to ' + reference + ': '
        hgvs_formatted.ac = reference
        variant.set_quibble(hgvs_formatted)

    lrg_reference = variant.quibble.ac

    if not lrg_reference.startswith('LRG_'):
        return

    lrg_id = lrg_reference[4:]

    # LRG transcript: LRG_123t1
    lrg_number, separator, transcript_number = lrg_id.partition('t')

    if (
            separator
            and lrg_number.isdigit()
            and transcript_number.isdigit()
    ):
        refseq_reference = (
            validator.db.get_refseq_transcript_id_from_lrg_transcript_id(
                lrg_reference
            )
        )

    else:
        # LRG protein: LRG_123p1
        lrg_number, separator, protein_number = lrg_id.partition('p')

        if (
                separator
                and lrg_number.isdigit()
                and protein_number.isdigit()
        ):
            refseq_reference = (
                validator.db.get_refseq_protein_id_from_lrg_protein_id(
                    lrg_reference
                )
            )

        # Genomic LRG: LRG_123
        elif lrg_id.isdigit():
            refseq_reference = validator.db.get_refseq_id_from_lrg_id(
                lrg_reference
            )

        else:
            return

    if refseq_reference != 'none':
        old_var_str = str(hgvs_formatted)
        hgvs_formatted.ac = refseq_reference
        variant.set_quibble(hgvs_formatted)

        caution += (
            old_var_str +
            ' automapped to equivalent RefSeq record ' +
            str(hgvs_formatted)
        )
        variant.warnings.append(caution)
        logger.info(caution)


def mitochondrial(variant, validator):
    """
    Will check if variant is mitochondrial and if so it will reformat the type
    to 'm' and save a value to the variant hgvs_genomic attribute.
    """
    mitochondrial_accessions = ('NC_012920.1', 'NC_001807.4')

    if (
            variant.reftype == ':m.'
            or variant.hgvs_formatted.ac in mitochondrial_accessions
    ):
        # Set flag
        variant.output_type_flag = 'mitochondrial'

        # Ensure the correct reference sequence type is used, if not, warn the user
        hgvs_mito = copy.deepcopy(variant.hgvs_formatted)

        if hgvs_mito.type == 'g' and hgvs_mito.ac in mitochondrial_accessions:
            hgvs_mito.type = 'm'

            if hgvs_mito.ac == 'NC_012920.1' and 'hg19' in variant.selected_assembly:
                wrn = (
                    'NC_012920.1 is not associated with genome build hg19, '
                    'instead use genome build GRCh37'
                )
                variant.warnings.append(wrn)
                return True

            elif hgvs_mito.ac == 'NC_001807.4' and 'GRCh37' in variant.selected_assembly:
                wrn = (
                    'NC_001807.4 is not associated with genome build GRCh37, '
                    'instead use genome build hg19'
                )
                variant.warnings.append(wrn)
                return True

            wrn = (
                'The given reference sequence (%s) does not match the DNA type (g). '
                'For %s, please use (m). For g. variants, please use a linear '
                'genomic reference sequence'
                % (hgvs_mito.ac, hgvs_mito.ac)
            )
            variant.warnings.append(wrn)

        # Validate the variant description
        try:
            validator.vr.validate(hgvs_mito)

        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.info(error)
            return True

        except KeyError:
            error = (
                'Currently unable to validate ' +
                hgvs_mito.ac +
                ' sequence variation'
            )
            variant.warnings.append(error)
            logger.info(error)
            return True

        # Check for movement during normalization
        try:
            norm_check = variant.hn.normalize(variant.hgvs_formatted)

            if hgvs_mito.posedit.pos != norm_check.posedit.pos:
                norm_check.type = 'm'
                error = '%s updated to %s' % (
                    fn.valstr(hgvs_mito),
                    fn.valstr(norm_check)
                )
                variant.warnings.append(error)
                logger.info(error)

        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.info(error)
            return True

        # Are there any transcripts?
        rel_var = validator.relevant_transcripts(
            hgvs_mito,
            variant.evm,
            validator.alt_aln_method,
            variant.reverse_normalizer,
            validator.select_transcripts
        )

        # Add a description of the reference sequence type and continue
        variant.hgvs_genomic = hgvs_mito

        if not rel_var:
            variant.genomic_g = unset_hgvs_obj_ref(hgvs_mito)
            variant.description = 'Homo sapiens mitochondrion, complete genome'
            logger.info('Homo sapiens mitochondrion, complete genome')
            return True

    return False


def _check_protein_reference(hgvs_object, validator, variant):
    """Check stated protein reference amino acids against the reference sequence."""
    start = hgvs_object.posedit.pos.start
    end = hgvs_object.posedit.pos.end

    start_pos = start.base
    end_pos = end.base

    check_start_aa = validator.sf.fetch_seq(
        hgvs_object.ac,
        start_i=start_pos - 1,
        end_i=start_pos
    )

    if start_pos == end_pos:
        check_end_aa = check_start_aa
    else:
        check_end_aa = validator.sf.fetch_seq(
            hgvs_object.ac,
            start_i=end_pos - 1,
            end_i=end_pos
        )

    valid_reference = True

    if start.aa != check_start_aa:
        variant.warnings.append(
            "The amino acid at position %s of %s is %s not %s" % (
                start_pos,
                hgvs_object.ac,
                check_start_aa,
                start.aa
            )
        )
        valid_reference = False

    if end.aa != check_end_aa:
        variant.warnings.append(
            "The amino acid at position %s of %s is %s not %s" % (
                end_pos,
                hgvs_object.ac,
                check_end_aa,
                end.aa
            )
        )
        valid_reference = False

    return valid_reference

def proteins(variant, validator):
    """Handle protein sequences."""
    if variant.reftype != ':p.':
        return False

    hgvs_object = variant.hgvs_formatted

    try:
        validator.vr.validate(hgvs_object)

    except AttributeError as e:
        if "AARefAlt' object has no attribute 'ref_n'" not in str(e):
            raise

        if not _check_protein_reference(hgvs_object, validator, variant):
            return True

    except vvhgvs.exceptions.HGVSError as e:
        error = str(e)
        variant.warnings.append(error)
        logger.info(error)
        return True

    error = (
        str(hgvs_object)
        + ' is HGVS compliant and contains a valid reference amino acid description'
    )
    reason = (
        'Protein level variant descriptions are not fully supported due to redundancy'
        ' in the genetic code'
    )

    variant.warnings.extend([reason, error])
    variant.protein = hgvs_object
    logger.info(reason + ": " + error)

    return True


def rna(variant, validator):
    """
    convert r, into c.
    """
    if variant.reftype == ':r.' or ":r." in variant.original:
        hgvs_input = variant.hgvs_formatted

        if hgvs_input.posedit.uncertain:
            hgvs_input.posedit.uncertain = False

        tx_info = validator.hdp.get_tx_identity_info(hgvs_input.ac)
        if tx_info[3] is None:
            error = "Invalid variant type for non-coding transcript. Instead use n."
            variant.warnings.append(error)
            logger.info(error)
            return True

        # Change to coding variant
        variant.reftype = ':c.'

        # Change input to reflect!
        try:
            hgvs_c = validator.hgvs_r_to_c(hgvs_input)
        except vvhgvs.exceptions.HGVSDataNotAvailableError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.info(error)
            return True
        variant.hgvs_formatted = hgvs_c

        # Create variant.rna_data dictionary
        rnd = VariantValidator.modules.rna_formatter.RnaDescriptions(validator.alt_aln_method,
                                                                     variant.primary_assembly,
                                                                     validator, variant)
        try:
            rnd.check_syntax(variant.original, variant)
        except VariantValidator.modules.rna_formatter.RnaVariantSyntaxError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.info(error)
            return True
        except vvhgvs.exceptions.HGVSParseError:
            try:
                rnd.check_syntax(variant.quibble, variant)
            except VariantValidator.modules.rna_formatter.RnaVariantSyntaxError as e:
                error = str(e)
                variant.warnings.append(error)
                logger.info(error)
                return True

        # Add data to the variant object
        variant.rna_data = rnd.output_dict()

    return False


def uncertain_pos(variant, validator):
    """
    Check for uncertain positions in the variant and return unsupported warning.
    """
    try:
        to_check = variant.quibble
        _accession, _sep, posedit = to_check.partition(':')

        if '(' not in posedit and ')' not in posedit:
            return False

        if 'p.' in posedit or '(;)' in posedit or '(:)' in posedit:
            return False

        if re.search(r'ins\(\d+(?:_\d+)?\)$', posedit):
            return False

        if ':r.(' in to_check:
            return False

        if (
                ('[' in posedit or ']' in posedit)
                and not re.search(r'\[\d+\]', posedit)
        ):
            return False

        logger.info(
            "Checking for uncertain positions in %s",
            variant.original
        )

        try:
            complex_descriptions.uncertain_positions(
                variant,
                validator
            )

        except complex_descriptions.UncertainConversionError as e:
            variant.warnings.append(
                f"UncertainConversionError: {e}"
            )
            return True

        except complex_descriptions.IncompatibleTypeError:
            use_checking.refseq_common_mistakes(variant)
            return True

        except (
                complex_descriptions.InvalidRangeError,
                complex_descriptions.FuzzyPositionError,
                complex_descriptions.FuzzyRangeError
        ) as e:
            variant.warnings.append(str(e))
            return True

        except Exception as e:
            variant.warnings.append(str(e))
            return True

        return False

    except Exception:
        return False


# <LICENSE>
# Copyright (C) 2016-2026 VariantValidator Contributors
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