import re
import logging

logger = logging.getLogger(__name__)

def remove_gene_symbol_from_ref(my_variant, validator):
    try:
    # 1. Requested warnings from https://github.com/openvar/variantValidator/issues/195
        if re.search(r'\(.+?\)', my_variant.quibble.ac):  # Pattern looks for (....)
            gene_symbol_query = re.search(r'\(.+?\)', my_variant.quibble.ac).group(0)
            gene_symbol_query = gene_symbol_query.replace('(', '')
            gene_symbol_query = gene_symbol_query.replace(')', '')
            is_it_a_gene = validator.db.get_hgnc_symbol(gene_symbol_query)
            if is_it_a_gene != 'none':
                warning = ("VariantSyntaxError: Removing redundant gene symbol %s from variant "
                           "description") % is_it_a_gene
                my_variant.quibble.ac = my_variant.quibble.ac.replace(f'({gene_symbol_query})', '')
                my_variant.warnings.append(warning)
                logger.warning(warning)
    except AttributeError:
        if re.search(r'\(.+?\)', my_variant.quibble):  # Pattern looks for (....)
            gene_symbol_query = re.search(r'\(.+?\)', my_variant.quibble).group(0)
            gene_symbol_query = gene_symbol_query.replace('(', '')
            gene_symbol_query = gene_symbol_query.replace(')', '')
            is_it_a_gene = validator.db.get_hgnc_symbol(gene_symbol_query)
            if is_it_a_gene != 'none':
                warning = ("VariantSyntaxError: Removing redundant gene symbol %s from variant "
                           "description") % is_it_a_gene
                my_variant.quibble = my_variant.quibble.replace(f'({gene_symbol_query})', '')
                my_variant.warnings.append(warning)
                logger.warning(warning)


def initial_user_formattng(my_variant, validator):
    # INITIAL USER INPUT FORMATTING
    """
    In this section of the code we are compiling HGVS errors and providing improved warnings/error
    messages
    """
    remove_gene_symbol_from_ref(my_variant, validator)

    if re.search('del[GATC]+', my_variant.original) or re.search('inv[GATC]+', my_variant.original) \
            or \
            re.search('dup[GATC]+', my_variant.original) or re.search('ins[GATC]+', my_variant.original):

        if not re.search('ins[GATC]+', my_variant.original):
            warning = "VariantSyntaxError: Removing redundant reference bases from variant description"
            my_variant.warnings.append(warning)
            logger.warning(warning)

    # 2. expand options for issue https://github.com/openvar/variantValidator/issues/338
    # Basically, all reference sequences must be upper case, so we make an upper-case query accession
    # to test the input accession against and try to spot a discrepancy
    # The exception to the rule is LTG transcripts e.g. LRG_1t1 which we handle immediately below!
    upper_case_accession = my_variant.quibble.ac.upper()
    original_ac, _sep, _remain = my_variant.original.partition(':')
    uc_original_ac = original_ac.upper()
    if uc_original_ac[:3] == "LRG":
        if "LRG" != original_ac[:3]:
            e = "This not a valid HGVS description, due to characters being in the wrong case. " \
                "Please check the use of upper- and lowercase characters."
            my_variant.warnings.append(str(e))
            logger.warning(str(e))
        if "T" in original_ac:
            e = "This not a valid HGVS description, due to characters being in the wrong case. " \
                "Please check the use of upper- and lowercase characters."
            my_variant.warnings.append(str(e))
            logger.warning(str(e))
            my_variant.quibble.ac = my_variant.quibble.ac.replace("T", "t")

    # Reference sequence types other than LRG
    elif original_ac != uc_original_ac and uc_original_ac[:3] != "LRG":
        # See issue #357
        if (uc_original_ac[:3] == 'CHR' or
                uc_original_ac[:4] == "GRCH" or
                uc_original_ac[:2] == "HG"):  # M already handled
            e = "This is not a valid HGVS variant description, because no reference sequence ID " \
                "has been provided"
        else:
            e = "This not a valid HGVS description, due to characters being in the wrong case. " \
                "Please check the use of upper- and lowercase characters."
        my_variant.warnings.append(str(e))
        logger.warning(str(e))
    elif (uc_original_ac[:3] == 'CHR' or
          uc_original_ac[:4] == "GRCH" or
          uc_original_ac[:2] == "HG"):
        e = "This is not a valid HGVS variant description, because no reference sequence ID " \
            "has been provided"
        my_variant.warnings.append(e)
        logger.warning(e)
