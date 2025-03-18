import re
import vvhgvs
import vvhgvs.exceptions
import vvhgvs.variantmapper
import logging
from . import utils as fn
from . import format_converters
import copy
from . import hgvs_utils

logger = logging.getLogger(__name__)

def pre_parsing_global_common_mistakes(my_variant):
    """
    A set of common error types that need to be found/handled before parsing variants into objects
    to compile HGVS errors and provide improved warnings/error messages initially from
    # INITIAL USER INPUT FORMATTING, which will end up as# INITIAL POST-OBJECT USER INPUT FORMATTING
    This may in fact want to be merged into the later use checking functions in the long term,
    or else may grow to handle more if some of these are converted to post-obj parsing instead
    """
    # test that it is not just a number or a numeric ID
    # since numeric ids may contain a : reverse quibble substitutions if otherwise fully numeric
    # e.g 1:111111 2:435636 12:30 would be treated as appropriate NC_ otherwise
    if re.match(r'^[\d\s\.,:;\-\+]+$',my_variant.original.strip()):
        my_variant.quibble = my_variant.original.strip()
    if re.match(r'\d', my_variant.quibble) and re.match(r'^[\d\s\.,:;\-\+]+$', my_variant.quibble):
        warning = "InvalidVariantError: VariantValidator operates on variant descriptions, but " +\
            f'this variant "{my_variant.quibble}" only contains numeric characters (and ' +\
            "possibly numeric associated punctuation), so can not be analysed. Did you enter this"+\
            " incorrectly, for example entering the numeric ID of a variant, instead of it`s " +\
            "description, or else enter just a within-sequence location, without specifying the " +\
            "actual variation?"
        my_variant.warnings.append(warning)
        return True

    # Find concatenated descriptions

    concat_descriptions = ["p\\.", "c\\.", "r\\.", "g\\.", "n\\."]  # Escape dots for regex
    pattern = f"({'|'.join(concat_descriptions)})"

    # Find all matches
    matches = re.findall(pattern, my_variant.original.strip())

    # Check if two or more matches are found
    if len(matches) >= 2:
        warning = (f"InvalidVariantError: {my_variant.original.strip()} is a concatenation of "
                   f"{'& '.join(matches)} descriptions, which is not compliant with the HGVS nomenclature standard")
        my_variant.warnings.append(warning)
        return True

    # some additional formats we have seen repeatedly
    matches = re.findall(":\w+:[crnpg]\.", my_variant.original.strip())
    if len(matches) >= 1:
        my_variant.warnings.append("VariantSyntaxError: HGVS descriptions contain a single colon between the reference "
                                   "sequence ID and the reference sequence type in the format "
                                   "'reference_sequence_ID:type.'")
        bad_chars = [x[1:-3] for x in matches]
        warning = (f"VariantSyntaxError: Illegal addition of the invalid characters [{'& '.join(bad_chars)}] between "
                   f"the two colons")
        my_variant.warnings.append(warning)
        return True

    matches = re.findall("[crng]\.[GATCgatc]\d+[GATCgatc]", my_variant.original.strip())
    if len(matches) >= 1:
        reformat = [f"{x[0:3]}>{x[-1]}" for x in matches]
        warning = (f"VariantSyntaxError: The format(s) {'& '.join(matches)} is(are) not compliant with the HGVS "
                   f"nomenclature standard for nucleotide variant descriptions. Did you mean {'& '.join(reformat)}?")
        my_variant.warnings.append(warning)
        return True

    invalid = my_variant.format_quibble()
    if invalid:
        if re.search(r'\w+:[gcnmrp],', my_variant.quibble):
            error = 'Variant description ' + my_variant.quibble + ' contained the , character between '\
                    '<type> and <position> in the expected pattern <accession>:<type>.<position> and ' \
                    'has been auto-corrected'
            my_variant.quibble = my_variant.quibble.replace(',', '.')
            my_variant.warnings.append(error)
            logger.warning(error)

        # Upper case type see issue #338
        elif re.search(r":[GCNMR].", str(my_variant.quibble)):
            rs_type_upper = re.search(r":[GCNMR].", str(my_variant.quibble))
            e = "This not a valid HGVS description, due to characters being in the wrong case. " \
                "Please check the use of upper- and lowercase characters."
            my_variant.warnings.append(str(e))
            logger.warning(str(e))
            my_variant.quibble = my_variant.quibble.replace(rs_type_upper.group(0),
                                                            rs_type_upper.group(0).lower())

        elif re.search(r'\w+:[gcnmrp]', my_variant.quibble) and not \
                re.search(r'\w+:[gcnmrp]\.', my_variant.quibble):
            error = 'Variant description ' + my_variant.quibble + ' lacks the . character between ' \
                    '<type> and <position> in the expected pattern <accession>:<type>.<position>'
            my_variant.warnings.append(error)
            logger.warning(error)
            return True

        else:
            error = 'Variant description ' + my_variant.quibble + ' is not in an accepted format'
            my_variant.warnings.append(error)
            logger.warning(error)
            return True

    # Here we handle syntax errors in ins and delins variants
    # https://github.com/openvar/variantValidator/issues/359
    if re.search("ins$", my_variant.quibble):
        my_variant.warnings.append("The inserted sequence must be provided for insertions or "
                                   "deletion-insertions")
        try:
            if "_" not in my_variant.quibble.split(":")[1] and \
                    "del" not in my_variant.quibble.split(":")[1]:
                my_variant.warnings.append("An insertion must be provided with the two positions "
                                           "between which the insertion has taken place")
        except IndexError:
            pass
        return True
    if re.search("ins\(\d+\)$", my_variant.quibble):
        my_variant.warnings.append("The length of the variant is not formatted following the "
                                   "HGVS guidelines. Please rewrite e.g. '(10)' to 'N[10]'"
                                   "(where N is an unknown nucleotide)")
        try:
            if "_" not in my_variant.quibble.split(":")[1] and \
                    "del" not in my_variant.quibble.split(":")[1]:
                my_variant.warnings.append("An insertion must be provided with the two positions "
                                           "between which the insertion has taken place")
        except IndexError:
            pass
        return True

    if re.search("ins\d+$", my_variant.quibble):
        my_variant.warnings.append("The length of the variant is not formatted following the HGVS "
                                   "guidelines. Please rewrite e.g. '10' to 'N[10]'"
                                   "(where N is an unknown nucleotide)")
        try:
            if "_" not in my_variant.quibble.split(":")[1] and \
                    "del" not in my_variant.quibble.split(":")[1]:
                my_variant.warnings.append("An insertion must be provided with the two positions "
                                           "between which the insertion has taken place")
        except IndexError:
            pass
        return True

    if re.search("ins\(\d+_\d+\)$", my_variant.quibble):
        my_variant.warnings.append("The length of the variant is not formatted following the HGVS "
                                   "guidelines. Please rewrite e.g. '(10_20)' to 'N[(10_20)]'"
                                   "(where N is an unknown nucleotide and [(10_20)] is an uncertain"
                                   " number of N nucleotides ranging from 10 to 20)")
        return True

    if re.search("ins\[\(\d+_\d+\)\]$", my_variant.quibble):
        counts = re.findall("\d+", my_variant.quibble.split("ins")[1])

        if int(counts[1]) < int(counts[0]):
            wrn = "The length of the variant is not formatted following the HGVS guidelines. " \
                  "Please rewrite (%s_%s) to N[(%s_%s)]" % (counts[0], counts[1],
                                                            counts[1], counts[0])
            my_variant.warnings.append(wrn)
        elif int(counts[1]) == int(counts[0]):
            wrn = "The length of the variant is not formatted following the HGVS guidelines. " \
                  "Please rewrite (%s_%s) to N[(%s)]" % (counts[0], counts[1], counts[1])
            my_variant.warnings.append(wrn)

        try:
            if not re.search("\d_\d", my_variant.quibble.split("ins")[0]) and \
                    "del" not in my_variant.quibble.split(":")[1]:
                my_variant.warnings.append("An insertion must be provided with the two positions "
                                           "between which the insertion has taken place")
        except IndexError:
            pass

        if my_variant.warnings == []:
            wrn = "The variant description is syntactically correct " \
                  "but no further validation is possible because the description contains " \
                  "uncertainty"
            my_variant.warnings.append(wrn)
        return True

    if re.search("(?:delins|del|ins)[NGATC]\[\d+\]$", my_variant.quibble) or \
            re.search("(?:delins|del|ins)\[[NGATC]\[\d+\];", my_variant.quibble):

        match = re.search("(?:delins|del|ins)", my_variant.quibble)[0]

        if re.search(f"{match}\[[GATCN]+\[\d+\];", my_variant.quibble):
            sections = my_variant.quibble.split(match)[1]
            sections = sections[1:-1]
            sections_listed = sections.split(";")
            sections_edited = []
            for stn in sections_listed:
                if '[' in stn and ']' in stn:
                    sections_edited.append(stn)
                else:
                    sections_edited.append(stn + "[1]")

            ins_seq_in_full = []
            for each_stn in sections_edited:
                bases, count = each_stn.split("[")
                count = int(count.replace("]", ""))
                for i in range(count):
                    ins_seq_in_full.append(bases)

        else:
            bases, count = my_variant.quibble.split("[")
            bases = bases.split(match)[1]
            count = int(count.replace("]", ""))
            ins_seq_in_full = []
            for i in range(count):
                ins_seq_in_full.append(bases)
        ins_seq_in_full = "".join(ins_seq_in_full)
        vt_in_full = my_variant.quibble.split(match)[0] + match + ins_seq_in_full
        warn = "%s may also be written as %s" % (my_variant.quibble, vt_in_full)
        my_variant.warnings.append(warn)

        try:
            if "_" not in my_variant.quibble.split(":")[1] and \
                    "del" not in my_variant.quibble.split(":")[1]:
                my_variant.warnings.append("An insertion must be provided with the two positions "
                                           "between which the insertion has taken place")
        except IndexError:
            pass

        # overwrite the current quibble for now instead of re-submiting for validation
        my_variant.quibble=vt_in_full
    return False

def refseq_common_mistakes(variant):
    """
    Evolving list of common mistakes, see sections below
    This is used both pre and post hgvs text to object conversion
    """
    # NM_ .g
    if type(variant.quibble) is str:
        acc4 = variant.quibble[:4]
        acc3 = variant.quibble[:3]
    else:
        acc4 = variant.quibble.ac[:4]
        acc3 = variant.quibble.ac[:3]

    if (acc3 in ['NM_','NR_'] or acc4 == 'ENST') and variant.reftype == ':g.':
        if acc3 == 'NR_' or variant.transcript_type == 'n':
            suggestion = str(variant.quibble).replace(':g.', ':n.')
        else:
            suggestion = str(variant.quibble).replace(':g.', ':c.')
        error = 'Transcript reference sequence input as genomic (g.) reference sequence. ' \
                'Did you mean ' + suggestion + '?'
        variant.warnings.append(error)
        logger.warning(error)
        return True

    # NR_ c.
    if variant.transcript_type == "n" and variant.reftype == ':c.':
        suggestion = str(variant.quibble).replace(':c.', ':n.')
        error = 'Non-coding transcript reference sequence input as coding (c.) reference sequence. ' \
                'Did you mean ' + suggestion + '?'
        variant.warnings.append(error)
        logger.warning(error)
        return True

    # NP_ c.
    if (acc3 == "NP_" or acc4 == "ENSP") and variant.reftype in [':c.', ':n.', ':g.', ':r.']:
        error = f'Protein reference sequence input as Nucleotide ({variant.reftype}) variant.'
        variant.warnings.append(error)
        logger.warning(error)
        return True

    # NM_ n.
    if variant.transcript_type == "c" and variant.reftype == ':n.':
        suggestion = str(variant.quibble).replace(':n.', ':c.')
        error = 'Coding transcript reference sequence input as non-coding transcript (n.) reference sequence. ' \
                'Did you mean ' + suggestion + '?'
        variant.warnings.append(error)
        logger.warning(error)
        return True

    # NM_ NC_ NG_ NR_ p.
    if (acc3 in ['NM_', 'NR_', 'NC_', 'NG_',] or acc4 == 'ENST') and variant.reftype == ':p.':
        error = 'Using a nucleotide reference sequence (NM_ NR_ NG_ NC_) to specify protein-level (p.) variation is ' \
                'not HGVS compliant. Please select an appropriate protein reference sequence (NP_)'
        variant.warnings.append(error)
        logger.warning(error)
        return True

    # NG_ c or NC_c..
    if acc3 in ['NG_', 'NC_'] and variant.reftype == ':c.':
        suggestion = 'For additional assistance, submit ' + str(variant.quibble) + ' to VariantValidator'
        error = 'NG_:c.PositionVariation descriptions should not be used unless a transcript reference sequence has ' \
                'also been provided e.g. NG_(NM_):c.PositionVariation'
        variant.warnings.extend([error, suggestion])
        logger.warning(error)
        return True

    return False


def structure_checks(variant, validator):
    """
    An evolving set of variant structure and content searches which identify
    and warn users about inappropriate use of HGVS
    Primarily, this code filters out variants that cannot realistically be
    auto corrected and will cause the downstream functions to return errors
    """
    if type(variant.quibble) is not str:
        input_parses = copy.deepcopy(variant.quibble)
    else:
        input_parses = validator.hp.parse_hgvs_variant(variant.quibble)
    variant.input_parses = input_parses
    variant.gene_symbol = validator.db.get_gene_symbol_from_transcript_id(variant.input_parses.ac)

    if variant.gene_symbol == 'none':
        variant.gene_symbol = ''
    if input_parses.type == 'g' or input_parses.type == 'm':
        check = structure_checks_g(variant, validator)
        if check:
            return True

    elif input_parses.type == 'c':
        check = structure_checks_c(variant, validator)
        if check:
            # Also check intron boundaries to provide additional warnings
            if "beyond the bounds" in str(variant.warnings) and (variant.input_parses.posedit.pos.start.offset != 0 or
                                                                 variant.input_parses.posedit.pos.end.offset != 0):
                if variant.input_parses.posedit.pos.start.offset != 0:
                    variant.input_parses.posedit.pos.start.offset = 1
                if variant.input_parses.posedit.pos.end.offset != 0:
                    variant.input_parses.posedit.pos.end.offset =1
                hgvs_genomic_vt = validator.myevm_t_to_g(variant.input_parses, variant.no_norm_evm,
                                                         variant.primary_assembly, variant.hn)
                format_converters.remap_intronic(variant.input_parses, hgvs_genomic_vt, variant, validator)
            return True

    elif input_parses.type == 'n':
        check = structure_checks_n(variant, validator)
        if check:
            return True
    else:
        pass


def structure_checks_g(variant, validator):
    """
    Structure checks for when reftype is genomic
    """
    if variant.input_parses.ac[:3] not in ['NC_', 'NG_', 'NT_', 'NW_']:
        error = 'Invalid reference sequence identifier (' + variant.input_parses.ac + ')'
        variant.warnings.append(error)
        logger.warning(error)
        return True

    try:
        validator.vr.validate(variant.input_parses)
    except Exception as e:
        if "does not agree with reference sequence ()" in str(e):
            e = "The specified coordinate is outside the boundaries of reference sequence %s " % variant.input_parses.ac

        if "base start position must be <= end position" in str(e) and (
                "NC_012920.1" in variant.hgvs_formatted.ac or
                "NC_001807.4" in variant.hgvs_formatted.ac):

            if variant.hgvs_formatted.ac not in variant.original:
                err = "This is not a valid HGVS variant description, because no reference sequence ID has " \
                      "been provided, instead use %s" % str(variant.hgvs_formatted)
                variant.warnings.append(err)
            variant.warnings.append("The variant positions are valid but we cannot normalize variants spanning "
                                    "the origin of circular reference sequences")
            return True

        elif "insertion length must be 1" in str(e) and "(" in str(variant.input_parses.posedit.pos) and ")" in \
                str(variant.input_parses.posedit.pos):
            return True
        elif "insertion length must be 1" in str(e) and "(" not in str(variant.input_parses.posedit.pos) and ")" not in\
                str(variant.input_parses.posedit.pos):
            ins_warning = (f'Insertion length must be 1 e.g. '
                           f'{str(int(variant.input_parses.posedit.pos.start.base))}'
                           f'_{str(int(variant.input_parses.posedit.pos.start.base)+1)}'
                           f'ins{variant.input_parses.posedit.edit.alt}')
            variant.warnings.append(ins_warning)
            for warning in variant.warnings:
                if warning == "insertion length must be 1":
                    variant.warnings.remove(warning)
            return True

            return True
        elif "Length implied by coordinates must equal sequence deletion length" in str(e) and \
             "(" in str(variant.input_parses.posedit.pos) and ")" in str(variant.input_parses.posedit.pos):
            return True
        else:
            error = str(e)
            variant.warnings.append(error)
            logger.warning(error)
            return True

    # Additional test
    try:
        np = variant.hn.normalize(variant.input_parses)
    except vvhgvs.exceptions.HGVSError as e:
        error = str(e)
        variant.warnings.append(error)
        logger.warning(error)
        return True

    # Look for variants in runs of N bases
    try:
        if "N" in variant.input_parses.posedit.edit.ref:
            error = (f"UncertainSequenceError: The submitted variant description {variant.input_parses} refers to a "
                     f"genomic reference region with "
                     f"an uncertain base composition (N)")
            variant.warnings.append(error)
            logger.warning(error)
            return True
    except TypeError:
        pass

    return False


def structure_checks_c(variant, validator):
    """
    structure checks for when reftype is coding
    :param variant:
    :param validator:
    :return:
    """

    if '*' in str(variant.input_parses) or 'c.-' in str(variant.input_parses):
        # Catch variation in UTRs
        # These should be in the sequence so can be directly validated. Need to pass to n.
        try:
            validator.vr.validate(variant.input_parses)
        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
            if 'datums is ill-defined' in error:
                called_ref = variant.input_parses.posedit.edit.ref

                try:
                    to_n = variant.evm.c_to_n(variant.input_parses)
                except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                    error = str(e)
                    variant.warnings.append(error)
                    logger.warning(error)
                    return True
                actual_ref = to_n.posedit.edit.ref

                if called_ref != actual_ref:
                    error = 'Variant reference (' + called_ref + ') does not agree with reference sequence ' \
                                                                 '(' + actual_ref + ')'
                    variant.warnings.append(error)
                    logger.warning(error)
                    return True
                else:
                    if variant.input_parses.posedit.edit.type == "ins":
                        variant.input_parses.posedit.edit.ref = None
                    else:
                        variant.input_parses.posedit.edit.ref = ''
                    variant.hgvs_formatted = variant.input_parses

            else:
                if 'bounds' in error or 'intronic variant' in error:
                    try:
                        variant.hn.normalize(variant.input_parses)
                    except vvhgvs.exceptions.HGVSError as e:
                        logger.debug("Except passed, %s", e)

                    if 'bounds' in error:
                        try:
                            identity_info = validator.hdp.get_tx_identity_info(variant.input_parses.ac)
                            ref_start = identity_info[3]
                            ref_end = identity_info[4]
                            if '-' in str(variant.input_parses.posedit.pos.start) and \
                                    variant.input_parses.posedit.pos.start.offset == 0:
                                # upstream positions
                                boundary = -ref_start
                                remainder = variant.input_parses.posedit.pos.start.base - boundary
                                variant.input_parses.posedit.pos.start.base = boundary
                                variant.input_parses.posedit.pos.start.offset = remainder
                            if '-' in str(variant.input_parses.posedit.pos.end) and \
                                    variant.input_parses.posedit.pos.end.offset == 0:
                                boundary = -ref_start
                                remainder = variant.input_parses.posedit.pos.end.base - boundary
                                variant.input_parses.posedit.pos.end.base = boundary
                                variant.input_parses.posedit.pos.end.offset = remainder
                            if '*' in str(variant.input_parses.posedit.pos.start) and \
                                    variant.input_parses.posedit.pos.start.offset == 0:
                                # downstream positions
                                tot_end_pos = str(variant.input_parses.posedit.pos.start).replace('*', '')
                                ts_seq = validator.sf.fetch_seq(variant.input_parses.ac)
                                boundary = len(ts_seq) - ref_end
                                variant.input_parses.posedit.pos.start.base = boundary
                                offset = int(tot_end_pos) - boundary
                                variant.input_parses.posedit.pos.start.offset = offset
                            if '*' in str(variant.input_parses.posedit.pos.end) and \
                                    variant.input_parses.posedit.pos.end.offset == 0:
                                tot_end_pos = str(variant.input_parses.posedit.pos.end).replace('*', '')
                                ts_seq = validator.sf.fetch_seq(variant.input_parses.ac)
                                boundary = len(ts_seq) - ref_end
                                variant.input_parses.posedit.pos.end.base = boundary
                                offset = int(tot_end_pos) - boundary
                                variant.input_parses.posedit.pos.end.offset = offset

                            # Create a lose vm instance
                            variant.lose_vm = vvhgvs.variantmapper.VariantMapper(validator.hdp,
                                                                                 replace_reference=True,
                                                                                 prevalidation_level=None)

                            report_gen = validator.myevm_t_to_g(variant.input_parses, variant.no_norm_evm,
                                                                variant.primary_assembly, variant.hn)
                            report_gen = variant.hn.normalize(report_gen)
                            error = 'Using a transcript reference sequence to specify a variant position that lies ' \
                                    'outside of the reference sequence is not HGVS-compliant: ' \
                                    'Instead re-submit ' + fn.valstr(report_gen)
                        except Exception as e:
                            logger.debug("Except passed, %s", e)
                        variant.warnings.append(error)
                        logger.warning(error)
                        return True

        try:
            variant.input_parses = variant.evm.c_to_n(variant.input_parses)
        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.warning(e)
            return True

        if 'n.1-' in str(variant.input_parses):
            input_parses = variant.evm.n_to_c(variant.input_parses)
            error = 'Using a transcript reference sequence to specify a variant position that lies outside of the ' \
                    'reference sequence is not HGVS-compliant. Instead re-submit '
            genomic_position = validator.myevm_t_to_g(input_parses, variant.no_norm_evm, variant.primary_assembly,
                                                      variant.hn)
            genomic_position = variant.hn.normalize(genomic_position)
            error = error + fn.valstr(genomic_position)
            variant.warnings.append(error)
            logger.warning(error)
            return True

        # Re-map input_parses back to c. variant
        variant.input_parses = variant.evm.n_to_c(variant.input_parses)

        # Intronic positions in UTRs
        if re.search(r'\d-\d', str(variant.input_parses)) or re.search(r'\d\+\d', str(variant.input_parses)):
            # Can we go c-g-c
            try:
                to_genome = validator.myevm_t_to_g(variant.input_parses, variant.no_norm_evm,
                                                   variant.primary_assembly, variant.hn)
                to_tx = variant.evm.g_to_t(to_genome, variant.input_parses.ac)
            except vvhgvs.exceptions.HGVSInvalidIntervalError as e:
                error = str(e)
                if 'bounds' in error:
                    try:
                        identity_info = validator.hdp.get_tx_identity_info(variant.input_parses.ac)
                        ref_start = identity_info[3]
                        ref_end = identity_info[4]
                        if '-' in str(variant.input_parses.posedit.pos.start):
                            # upstream positions
                            boundary = -ref_start
                            remainder = variant.input_parses.posedit.pos.start.base - boundary
                            variant.input_parses.posedit.pos.start.base = boundary
                            variant.input_parses.posedit.pos.start.offset = remainder
                        if '-' in str(variant.input_parses.posedit.pos.end):
                            boundary = -ref_start
                            remainder = variant.input_parses.posedit.pos.end.base - boundary
                            variant.input_parses.posedit.pos.end.base = boundary
                            variant.input_parses.posedit.pos.end.offset = remainder
                        if '*' in str(variant.input_parses.posedit.pos.start):
                            # downstream positions
                            tot_end_pos = str(variant.input_parses.posedit.pos.start).replace('*', '')
                            ts_seq = validator.sf.fetch_seq(variant.input_parses.ac)
                            boundary = len(ts_seq) - ref_end
                            variant.input_parses.posedit.pos.start.base = boundary
                            te1, te2 = tot_end_pos.split('+')
                            tot_end_pos = int(te1) + int(te2)
                            offset = tot_end_pos - boundary
                            variant.input_parses.posedit.pos.start.offset = offset
                        if '*' in str(variant.input_parses.posedit.pos.end):
                            tot_end_pos = str(variant.input_parses.posedit.pos.end).replace('*', '')
                            ts_seq = validator.sf.fetch_seq(variant.input_parses.ac)
                            boundary = len(ts_seq) - ref_end
                            variant.input_parses.posedit.pos.end.base = boundary
                            te1, te2 = tot_end_pos.split('+')
                            tot_end_pos = int(te1) + int(te2)
                            offset = tot_end_pos - boundary
                            variant.input_parses.posedit.pos.end.offset = offset

                        report_gen = validator.myevm_t_to_g(variant.input_parses, variant.no_norm_evm,
                                                            variant.primary_assembly, variant.hn)
                        report_gen = variant.hn.normalize(report_gen)
                        error = 'Using a transcript reference sequence to specify a variant position that lies ' \
                                'outside of the reference sequence is not HGVS-compliant. Instead re-submit '\
                                + fn.valstr(report_gen)
                    except Exception as e:
                        logger.debug("Except passed, %s", e)
                variant.warnings.append(error)
                logger.warning(error)
                return True

            except vvhgvs.exceptions.HGVSDataNotAvailableError as e:
                error = str(e)
                if 'Alignment is incomplete' in error:
                    e_list = error.split('~')
                    gens = []
                    for el in e_list:
                        el_l = el.split('/')
                        if el_l[-1] == '':
                            continue
                        gens.append(el_l[-1])
                    acs = '; '.join(gens)
                    error = 'Cannot map ' + fn.valstr(variant.input_parses) + ' to a genomic position. '\
                            + variant.input_parses.ac + ' can only be partially aligned to genomic reference ' \
                            'sequences ' + acs
                    variant.warnings.append(error)
                    logger.warning(error)
                    return True

    elif re.search(r'\d-', str(variant.input_parses)) or re.search(r'\d\+', str(variant.input_parses)):
        # Create a no_replace vm instance
        variant.no_replace_vm = vvhgvs.variantmapper.VariantMapper(validator.hdp,
                                                                   replace_reference=False,
                                                                   prevalidation_level=None)

        # Quick look at syntax validation
        try:
            validator.vr.validate(variant.input_parses)
        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)
            if 'bounds' in error:
                try:
                    report_gen = validator.myevm_t_to_g(variant.input_parses, variant.no_norm_evm,
                                                        variant.primary_assembly, variant.hn)
                    report_gen = variant.hn.normalize(report_gen)
                except vvhgvs.exceptions.HGVSError as e:
                    logger.debug("Except passed, %s", e)
                else:
                    error = 'Using a transcript reference sequence to specify a variant position that lies outside of '\
                            'the reference sequence is not HGVS-compliant. Instead re-submit ' + fn.valstr(report_gen)
                variant.warnings.append(error)
                logger.warning(error)
                return True
            elif 'insertion length must be 1' in error:
                    if "-" in str(variant.input_parses.posedit.pos.start.offset):
                        start_offset = str(variant.input_parses.posedit.pos.start.offset - 1)
                        end_offset = str(variant.input_parses.posedit.pos.start.offset)
                    elif ("-" not in str(variant.input_parses.posedit.pos.start.offset) and
                          variant.input_parses.posedit.pos.start.offset != 0):
                        start_offset = f"+{str(variant.input_parses.posedit.pos.start.offset)}"
                        end_offset = f"+{str(variant.input_parses.posedit.pos.start.offset + 1)}"
                    if "(" not in str(variant.input_parses.posedit.pos):
                        ins_warning = (f'Insertion length must be 1 e.g. '
                                       f'{variant.input_parses.posedit.pos.start.base}{start_offset}'
                                       f'_{str(int(variant.input_parses.posedit.pos.start.base))}{end_offset}'
                                       f'ins{variant.input_parses.posedit.edit.alt}')
                        variant.warnings.append(ins_warning)
                        for warning in variant.warnings:
                            if warning == "insertion length must be 1":
                                variant.warnings.remove(warning)
                    return True

            elif 'base start position must be <= end position' in error:
                correction = copy.deepcopy(variant.input_parses)
                st = variant.input_parses.posedit.pos.start
                ed = variant.input_parses.posedit.pos.end
                correction.posedit.pos.start = ed
                correction.posedit.pos.end = st
                error = error + ': Did you mean ' + str(correction) + '?'
                variant.warnings.append(error)
                logger.warning(error)
                return True

        # Create a specific minimal evm with no normalizer and no replace_reference
        # Have to use this method due to potential multi chromosome error, note normalizes but does not replace sequence
        output = None
        try:
            output = validator.noreplace_myevm_t_to_g(variant.input_parses, variant)
        except vvhgvs.exceptions.HGVSDataNotAvailableError:
            tx_info = validator.hdp.get_tx_identity_info(variant.input_parses.ac)
            if (variant.input_parses.posedit.pos.end.base > int(tx_info[4]) or variant.input_parses.posedit.pos.end.base
                > int(tx_info[4])) and ("*" not in str(variant.input_parses.posedit.pos.end) or "*" not in
                                        str(variant.input_parses.posedit.pos.start)):
                errors = ["CDSError: Variant start position and/or end position are beyond the CDS end position "
                          "and likely also beyond the end of the selected reference sequence"]
            else:
                errors = ['Required information for ' + variant.input_parses.ac + ' is missing from the Universal '
                          'Transcript Archive', 'Query gene2transcripts with search term %s for '
                          'available transcripts' % variant.input_parses.ac.split('.')[0]]
            variant.warnings.extend(errors)
            logger.info(str(errors))
            return True
        except ValueError as e:
            error = str(e)
            if '> end' in error:
                error = 'Interval start position ' + str(variant.input_parses.posedit.pos.start) + ' > interval end '\
                        'position ' + str(variant.input_parses.posedit.pos.end)
                variant.warnings.append(error)
                logger.warning(error)
                return True
        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)
            if 'base start position must be <= end position' in error:
                # correction = copy.deepcopy(variant.input_parses)
                # st = variant.input_parses.posedit.pos.start
                # ed = variant.input_parses.posedit.pos.end
                # correction.posedit.pos.start = ed
                # correction.posedit.pos.end = st
                # error = error + ': Did you mean ' + str(correction) + '?'
                error = 'Interval start position ' + str(variant.input_parses.posedit.pos.start) + ' > interval end' \
                        ' position ' + str(variant.input_parses.posedit.pos.end)
                variant.warnings.append(error)
                logger.warning(error)
                return True
            else:
                variant.warnings.append(error)
                logger.warning(error)
                return True

        try:
            variant.evm.g_to_t(output, variant.input_parses.ac)
        except vvhgvs.exceptions.HGVSError as e:
            if "Alignment is incomplete" in str(e):
                output = hgvs_utils.incomplete_alignment_mapping_t_to_g(validator, variant)
                if output is None:
                    error = str(e)
                    variant.warnings.append(error)
                    logger.warning(error)
                    return True
            else:
                error = str(e)
                variant.warnings.append(error)
                logger.warning(error)
                return True

        # Check that the reference is correct by direct mapping without replacing reference
        check_ref_g = variant.no_replace_vm.t_to_g(variant.input_parses, output.ac,
                                                   alt_aln_method=validator.alt_aln_method)
        check_ref_t = variant.no_replace_vm.g_to_t(check_ref_g, variant.input_parses.ac,
                                                   alt_aln_method=validator.alt_aln_method)

        # Snapshot current variant error log
        if "*" in str(check_ref_t) and "*" not in str(variant.input_parses):
            convert = "%s auto-mapped to %s" % (variant.input_parses, check_ref_t)
            variant.warnings.append(convert)
        snap = copy.copy(variant.warnings)

        # Look for syntax errors
        try:
            validator.vr.validate(check_ref_t)
        except vvhgvs.exceptions.HGVSError as e:
            if "intron" not in str(e) and "bounds" not in str(e) and "insertion length must be 1" not in str(e) and \
                    "base start position must be <= end position" not in str(e):
                error = str(e)
                variant.warnings.append(error)
                logger.warning(error)
        try:
            validator.vr.validate(output)
        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.warning(error)

        # Check for additional warnings
        if len(variant.warnings) > len(snap):
            return True

    else:
        # All other variation
        try:
            validator.vr.validate(variant.input_parses)
        except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
            logger.debug("Except passed, %s", e)
        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)
            # This catches errors in introns
            if 'base start position must be <= end position' in error:
                error = 'Interval start position ' + str(variant.input_parses.posedit.pos.start) + ' > interval end '\
                        'position ' + str(variant.input_parses.posedit.pos.end)
            if "(" not in str(variant.input_parses.posedit.pos):
                variant.warnings.append(error)
                logger.warning(error)
            if 'insertion length must be 1' in error:
                if "(" not in str(variant.input_parses.posedit.pos):
                    ins_warning = (f'Insertion length must be 1 e.g. '
                                   f'{str(int(variant.input_parses.posedit.pos.start.base))}'
                                   f'_{str(int(variant.input_parses.posedit.pos.start.base)+1)}'
                                   f'ins{variant.input_parses.posedit.edit.alt}')
                    variant.warnings.append(ins_warning)
                    for warning in variant.warnings:
                        if warning == "insertion length must be 1":
                            variant.warnings.remove(warning)
            return True

        except vvhgvs.exceptions.HGVSDataNotAvailableError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.warning(error)
            return True
        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
            if 'bounds' in error:
                error += ' (' + variant.input_parses.ac + ')'
                variant.warnings.append(error)
                logger.warning(error)
                return True

    return False


def structure_checks_n(variant, validator):
    """
    structure checks for reftype nucleotide
    :param variant:
    :param validator:
    :return:
    """
    if '+' in str(variant.input_parses) or '-' in str(variant.input_parses):
        # Catch variation in UTRs
        # These should be in the sequence so can be directly validated. Need to pass to n.
        try:
            validator.vr.validate(variant.input_parses)
        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
            if 'intronic variant' in error:
                pass
            elif 'datums is ill-defined' in error:
                called_ref = variant.input_parses.posedit.edit.ref
                to_n = variant.evm.c_to_n(variant.input_parses)
                actual_ref = to_n.posedit.edit.ref
                if called_ref != actual_ref:
                    error = 'Variant reference (' + called_ref + ') does not agree with reference sequence (' + \
                            actual_ref + ')'
                    variant.warnings.append(error)
                    logger.warning(str(error))
                    return True
                else:
                    variant.input_parses.posedit.edit.ref = ''
                    variant.hgvs_formatted = variant.input_parses

            elif 'base must be >=1 for datum = SEQ_START or CDS_END' in error:
                error = 'The given coordinate is outside the bounds of the reference sequence.'

                try:
                    if '-' in str(variant.input_parses.posedit.pos.start):
                        # upstream positions
                        boundary = 1
                        remainder = variant.input_parses.posedit.pos.start.base - boundary
                        remainder = remainder + 1
                        variant.input_parses.posedit.pos.start.base = boundary
                        variant.input_parses.posedit.pos.start.offset = remainder
                    if '-' in str(variant.input_parses.posedit.pos.end):
                        boundary = 1
                        remainder = variant.input_parses.posedit.pos.end.base - boundary
                        remainder = remainder + 1
                        variant.input_parses.posedit.pos.end.base = boundary
                        variant.input_parses.posedit.pos.end.offset = remainder
                    report_gen = validator.myevm_t_to_g(variant.input_parses, variant.no_norm_evm,
                                                        variant.primary_assembly, variant.hn)
                    report_gen = variant.hn.normalize(report_gen)
                    error = 'Using a transcript reference sequence to specify a variant position that lies outside of' \
                            ' the reference sequence is not HGVS-compliant. Instead re-submit ' + fn.valstr(report_gen)
                except Exception as e:
                    logger.debug("Except passed, %s", e)
                variant.warnings.append(error)
                logger.warning(error)
                return True
            else:
                variant.warnings.append(error)
                logger.warning(error)
                return True

    if 'n.1-' in str(variant.input_parses):
        error = 'Using a transcript reference sequence to specify a variant position that lies outside of the ' \
                'reference sequence is not HGVS-compliant. Instead re-submit '
        genomic_position = validator.myevm_t_to_g(variant.input_parses, variant.no_norm_evm, variant.primary_assembly,
                                                  variant.hn)
        genomic_position = variant.hn.normalize(genomic_position)
        error = error + fn.valstr(genomic_position)
        variant.warnings.append(error)
        logger.warning(error)
        return True

    if re.search(r'\d-', str(variant.input_parses)) or re.search(r'\d\+', str(variant.input_parses)):
        # Quick look at syntax validation
        try:
            validator.vr.validate(variant.input_parses)
        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)
            if 'bounds' in error:
                try:
                    report_gen = validator.myevm_t_to_g(variant.input_parses, variant.no_norm_evm,
                                                        variant.primary_assembly, variant.hn)
                    report_gen = variant.hn.normalize(report_gen)
                except vvhgvs.exceptions.HGVSError as e:
                    logger.debug("Except passed, %s", e)
                else:
                    error = 'Using a transcript reference sequence to specify a variant position that lies outside of '\
                            'the reference sequence is not HGVS-compliant. Instead re-submit ' + fn.valstr(report_gen)
                variant.warnings.append(error)
                logger.warning(error)
                return True
            elif 'insertion length must be 1' in error:
                variant.warnings.append(error)
                logger.warning(error)
                return True
            elif 'base start position must be <= end position' in error:
                correction = copy.deepcopy(variant.input_parses)
                st = variant.input_parses.posedit.pos.start
                ed = variant.input_parses.posedit.pos.end
                correction.posedit.pos.start = ed
                correction.posedit.pos.end = st
                error = error + ': Did you mean ' + str(correction) + '?'
                # error = 'Interval start position ' + str(input_parses.posedit.pos.start) + ' > interval end
                # position ' + str(input_parses.posedit.pos.end)
                variant.warnings.append(error)
                logger.warning(error)
                return True
            elif 'Cannot validate sequence of an intronic variant' in error:
                try:
                    test_g = validator.myevm_t_to_g(variant.input_parses, variant.no_norm_evm, variant.primary_assembly,
                                                    variant.hn)
                    variant.evm.g_to_t(test_g, variant.input_parses.ac)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'bounds' in error:
                        report_gen = validator.myevm_t_to_g(variant.input_parses, variant.no_norm_evm,
                                                            variant.primary_assembly, variant.hn)
                        report_gen = variant.hn.normalize(report_gen)
                        error = 'Using a transcript reference sequence to specify a variant position that lies ' \
                                'outside of the reference sequence is not HGVS-compliant. Instead re-submit ' + \
                                fn.valstr(report_gen)
                        variant.warnings.append(error)
                        logger.warning(error)
                        return True

        # Create a specific minimal evm with no normalizer and no replace_reference
        # Have to use this method due to potential multi chromosome error, note, normalizes but does not replace
        # sequence
        output = None
        try:
            output = validator.noreplace_myevm_t_to_g(variant.input_parses, variant)
        except vvhgvs.exceptions.HGVSDataNotAvailableError:
            errors = ['Required information for ' + variant.input_parses.ac + ' is missing from the Universal '
                                                                              'Transcript Archive',
                      'Query gene2transcripts with search term %s for '
                      'available transcripts' % variant.input_parses.ac.split('.')[0]]
            variant.warnings.extend(errors)
            logger.info(str(errors))
            return True
        except ValueError as e:
            error = str(e)
            if '> end' in error:
                error = 'Interval start position ' + str(
                    variant.input_parses.posedit.pos.start) + ' > interval end position ' + str(
                    variant.input_parses.posedit.pos.end)
                variant.warnings.append(error)
                logger.warning(error)
                return True
        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)
            if 'base start position must be <= end position' in error:
                correction = copy.deepcopy(variant.input_parses)
                st = variant.input_parses.posedit.pos.start
                ed = variant.input_parses.posedit.pos.end
                correction.posedit.pos.start = ed
                correction.posedit.pos.end = st
                # error = error + ': Did you mean ' + str(correction) + '?'
                error = 'Interval start position ' + str(
                    variant.input_parses.posedit.pos.start) + ' > interval end position ' + str(
                    variant.input_parses.posedit.pos.end)
                variant.warnings.append(error)
                logger.warning(error)
                return True
        try:
            validator.vr.validate(output)
        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.warning(error)
            return True

    else:
        # All other variation
        try:
            validator.vr.validate(variant.input_parses)
        except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
            logger.debug("Except passed, %s", e)
        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)
            # if re.search('Length implied by coordinates', error):
            #     # Applies to del and inv
            #     # NOTE, there has been no normalization at all so this error is valid here
            #     my_variant.warnings += ': ' + str(error)
            #     continue
            # # Will apply to > del and inv
            # if re.search('does not agree with reference sequence', error):
            #     my_variant.warnings += ': ' + str(error)
            #     continue
            # # ensures x_y for insertions
            # if re.search('insertion length must be 1', error):
            #     my_variant.warnings += ': ' + str(error)
            #     continue
            # # Boundary issue
            # if re.search('Variant coordinate is out of the bound of CDS region', error):
            #     my_variant.warnings += ': ' + str(error)
            #     continue
            # This catches errors in introns
            if 'base start position must be <= end position' in error:
                # correction = copy.deepcopy(variant.input_parses)
                # st = variant.input_parses.posedit.pos.start
                # ed = variant.input_parses.posedit.pos.end
                # correction.posedit.pos.start = ed
                # correction.posedit.pos.end = st
                # error = error + ': Did you mean ' + str(correction) + '?'
                error = 'Interval start position ' + str(
                    variant.input_parses.posedit.pos.start) + ' > interval end position ' + str(
                    variant.input_parses.posedit.pos.end)
                logger.warning(error)
                variant.warnings.append(error)
                return True
            variant.warnings.append(error)
            logger.warning(error)
            return True
        except vvhgvs.exceptions.HGVSDataNotAvailableError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.warning(error)
            return True
        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
            if 'bounds' in error:
                error = error + ' (' + variant.input_parses.ac + ')'
                variant.warnings.append(error)
                logger.warning(error)
                return True
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
