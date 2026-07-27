import re
import vvhgvs
import vvhgvs.exceptions
import vvhgvs.variantmapper
import logging
from . import utils as fn
from . import format_converters
import copy
from . import hgvs_utils, hgvs_position_utils

logger = logging.getLogger(__name__)

class InvalidVariantError(Exception):
    pass


def pre_parsing_global_common_mistakes(my_variant):
    """
    A set of common error types that need to be found/handled before parsing
    variants into objects, to compile HGVS errors and provide improved
    warnings/error messages initially from INITIAL USER INPUT FORMATTING.

    This may in fact want to be merged into the later use checking functions
    in the long term, or else may grow to handle more if some of these are
    converted to post-object parsing instead.
    """

    original = my_variant.original
    original_stripped = original.strip()
    quibble = my_variant.quibble

    # Test that it is not just a number or a numeric ID.
    # Numeric IDs may contain punctuation such as:
    # 1:111111, 2:435636, 12:30.

    if re.fullmatch(r'[\d\s.,:;\-+]+', original_stripped):
        quibble = original_stripped
        my_variant.quibble = quibble

    if quibble and quibble[0].isdigit() and re.fullmatch(
            r'[\d\s.,:;\-+]+',
            quibble
    ):
        warning = (
            "InvalidVariantError: VariantValidator operates on variant "
            f'descriptions, but this variant "{original}" only contains '
            "numeric characters (and possibly numeric associated punctuation), "
            "so can not be analysed. Did you enter this incorrectly, for "
            "example entering a gene ID without specifying the actual "
            "variation? If so, try our genes to transcripts tool "
            "https://variantvalidator.org/service/gene2trans/"
        )
        my_variant.warnings.append(warning)
        return True

    if re.fullmatch(r'[\w+]+', quibble):
        warning = (
            "InvalidVariantError: VariantValidator operates on variant "
            f'descriptions, but this variant "{original}" only contains '
            "alphanumeric characters so can not be analysed. Did you enter "
            "this incorrectly, for example entering a gene symbol without "
            "specifying the actual variation? If so, try our genes to "
            "transcripts tool "
            "https://variantvalidator.org/service/gene2trans/"
        )
        my_variant.warnings.append(warning)
        return True

    # Find concatenated descriptions.
    matches = re.findall(r'[pcrgn]\.', original_stripped)

    if len(matches) >= 2:
        warning = (
            f"InvalidVariantError: {original_stripped} is a concatenation of "
            f"{'& '.join(matches)} descriptions, which is not compliant with "
            "the HGVS nomenclature standard"
        )
        my_variant.warnings.append(warning)
        return True

    # Additional malformed formats seen repeatedly.
    matches = re.findall(r':\w+:[crnpg]\.', original_stripped)

    if matches:
        my_variant.warnings.append(
            "VariantSyntaxError: HGVS descriptions contain a single colon "
            "between the reference sequence ID and the reference sequence "
            "type in the format 'reference_sequence_ID:type.'"
        )

        bad_chars = [match[1:-3] for match in matches]

        warning = (
            "VariantSyntaxError: Illegal addition of the invalid characters "
            f"[{'& '.join(bad_chars)}] between the two colons in "
            f"{original_stripped}"
        )
        my_variant.warnings.append(warning)
        return True

    matches = re.findall(
        r':\w+\.\d+:[crnpg]\.',
        original_stripped
    )

    if matches:
        my_variant.warnings.append(
            "VariantSyntaxError: HGVS descriptions contain a single colon "
            "between the reference sequence ID and the reference sequence "
            "type in the format 'reference_sequence_ID:type.'"
        )

        bad_chars = [match[1:-3] for match in matches]

        warning = (
            "VariantSyntaxError: Illegal addition of the invalid characters "
            f"[{'& '.join(bad_chars)}] between the two colons in "
            f"{original_stripped}"
        )
        my_variant.warnings.append(warning)
        return True

    matches = re.findall(
        r'[crng]\.[GATCgatc]\d+[GATCgatc]',
        original_stripped
    )

    if matches:
        reformat = [
            f"{match[:3]}>{match[-1]}"
            for match in matches
        ]

        warning = (
            f"VariantSyntaxError: The format(s) {'& '.join(matches)} is(are) "
            "not compliant with the HGVS nomenclature standard for nucleotide "
            f"variant descriptions. Did you mean {'& '.join(reformat)}?"
        )
        my_variant.warnings.append(warning)
        return True

    invalid = my_variant.format_quibble()

    # format_quibble() may modify quibble, so refresh our local reference.
    quibble = my_variant.quibble

    if invalid:
        if re.search(r'\w+:[gcnmrp],', quibble):
            error = (
                f"Variant description {quibble} contained the , character "
                "between <type> and <position> in the expected pattern "
                "<accession>:<type>.<position> and has been auto-corrected"
            )

            quibble = quibble.replace(',', '.')
            my_variant.quibble = quibble
            my_variant.warnings.append(error)
            logger.info(error)

        else:
            rs_type_upper = re.search(r':[GCNMR].', quibble)

            if rs_type_upper:
                error = (
                    "This not a valid HGVS description, due to characters "
                    "being in the wrong case. Please check the use of upper- "
                    "and lowercase characters."
                )
                my_variant.warnings.append(error)
                logger.info(error)

                matched_type = rs_type_upper.group(0)
                quibble = quibble.replace(
                    matched_type,
                    matched_type.lower()
                )
                my_variant.quibble = quibble

            elif (
                    re.search(r'\w+:[gcnmrp]', quibble)
                    and not re.search(r'\w+:[gcnmrp]\.', quibble)
            ):
                error = (
                    f"Variant description {quibble} lacks the . character "
                    "between <type> and <position> in the expected pattern "
                    "<accession>:<type>.<position>"
                )
                my_variant.warnings.append(error)
                logger.info(error)
                return True

            elif (
                    (
                        re.search(r'\(ENST\d+\.\d+\):', quibble)
                        or re.search(r'\(N[MRCG]_\d+\.\d+\):', quibble)
                        or re.search(r'\(LRG_\d+t\d+\):', quibble)
                    )
                    and not quibble.startswith('NC_')
            ):
                reference_region, variation = quibble.split(':', 1)
                reference = reference_region.split('(', 1)[1].replace(')', '')

                new_quibble = f"{reference}:{variation}"

                my_variant.warnings.append(
                    "VariantSyntaxError: Stripping unnecessary characters "
                    f"from {quibble} and updating to {new_quibble}"
                )

                quibble = new_quibble
                my_variant.quibble = quibble

            else:
                error = (
                    f"InvalidVariantError: Variant description {original} "
                    "is not in an accepted format"
                )
                my_variant.warnings.append(error)
                logger.info(error)
                raise InvalidVariantError(error)

    # Here we handle syntax errors in ins and delins variants.
    # https://github.com/openvar/variantValidator/issues/359

    if quibble.endswith('ins'):
        my_variant.warnings.append(
            "The inserted sequence must be provided for insertions or "
            "deletion-insertions"
        )

        try:
            variation = quibble.split(':', 1)[1]

            if '_' not in variation and 'del' not in variation:
                my_variant.warnings.append(
                    "An insertion must be provided with the two positions "
                    "between which the insertion has taken place"
                )
        except IndexError:
            pass

        return True

    if re.search(r'ins\(\d+\)$', quibble):
        my_variant.warnings.append(
            "The length of the variant is not formatted following the HGVS "
            "guidelines. Please rewrite e.g. '(10)' to 'N[10]'"
            "(where N is an unknown nucleotide)"
        )

        try:
            variation = quibble.split(':', 1)[1]

            if '_' not in variation and 'del' not in variation:
                my_variant.warnings.append(
                    "An insertion must be provided with the two positions "
                    "between which the insertion has taken place"
                )
        except IndexError:
            pass

        return True

    if re.search(r'ins\d+$', quibble):
        my_variant.warnings.append(
            "The length of the variant is not formatted following the HGVS "
            "guidelines. Please rewrite e.g. '10' to 'N[10]'"
            "(where N is an unknown nucleotide)"
        )

        try:
            variation = quibble.split(':', 1)[1]

            if '_' not in variation and 'del' not in variation:
                my_variant.warnings.append(
                    "An insertion must be provided with the two positions "
                    "between which the insertion has taken place"
                )
        except IndexError:
            pass

        return True

    if re.search(r'ins\(\d+_\d+\)$', quibble):
        my_variant.warnings.append(
            "The length of the variant is not formatted following the HGVS "
            "guidelines. Please rewrite e.g. '(10_20)' to 'N[(10_20)]'"
            "(where N is an unknown nucleotide and [(10_20)] is an uncertain "
            "number of N nucleotides ranging from 10 to 20)"
        )
        return True

    if re.search(r'ins\[\(\d+_\d+\)\]$', quibble):
        insertion = quibble.split('ins', 1)[1]
        counts = re.findall(r'\d+', insertion)

        lower = int(counts[0])
        upper = int(counts[1])

        if upper < lower:
            warning = (
                "The length of the variant is not formatted following the "
                "HGVS guidelines. Please rewrite "
                f"({lower}_{upper}) to N[({upper}_{lower})]"
            )
            my_variant.warnings.append(warning)

        elif upper == lower:
            warning = (
                "The length of the variant is not formatted following the "
                "HGVS guidelines. Please rewrite "
                f"({lower}_{upper}) to N[({upper})]"
            )
            my_variant.warnings.append(warning)

        try:
            before_ins = quibble.split('ins', 1)[0]
            variation = quibble.split(':', 1)[1]

            if (
                    not re.search(r'\d_\d', before_ins)
                    and 'del' not in variation
            ):
                my_variant.warnings.append(
                    "An insertion must be provided with the two positions "
                    "between which the insertion has taken place"
                )
        except IndexError:
            pass

        if not my_variant.warnings:
            warning = (
                "The variant description is syntactically correct but no "
                "further validation is possible because the description "
                "contains uncertainty"
            )
            my_variant.warnings.append(warning)

        return True

    repeat_match = re.search(
        r'(?:delins|del|ins)[NGATC]+\[\d+\]$',
        quibble
    )

    section_repeat_match = re.search(
        r'(?:delins|del|ins)\[[NGATC]+\[\d+\];',
        quibble
    )

    if repeat_match or section_repeat_match:
        edit_match = re.search(r'(?:delins|del|ins)', quibble)
        edit_type = edit_match.group(0)

        if re.search(
                rf'{edit_type}\[[GATCN]+\[\d+\];',
                quibble
        ):
            sections = quibble.split(edit_type, 1)[1][1:-1]
            sections_listed = sections.split(';')

            expanded_sections = []

            for section in sections_listed:
                if '[' in section and ']' in section:
                    bases, count = section.split('[', 1)
                    count = int(count.rstrip(']'))
                else:
                    bases = section
                    count = 1

                expanded_sections.append(bases * count)

            ins_seq_in_full = ''.join(expanded_sections)

        else:
            repeat_section = quibble.split(edit_type, 1)[1]
            bases, count = repeat_section.split('[', 1)
            count = int(count.rstrip(']'))
            ins_seq_in_full = bases * count

        vt_in_full = (
            quibble.split(edit_type, 1)[0]
            + edit_type
            + ins_seq_in_full
        )

        warning = f"{quibble} may also be written as {vt_in_full}"
        my_variant.warnings.append(warning)

        try:
            variation = quibble.split(':', 1)[1]

            if '_' not in variation and 'del' not in variation:
                my_variant.warnings.append(
                    "An insertion must be provided with the two positions "
                    "between which the insertion has taken place"
                )
        except IndexError:
            pass

        # Overwrite the current quibble for now instead of re-submitting
        # for validation.
        quibble = vt_in_full
        my_variant.quibble = quibble

    # Ranges in substitutions.
    if '>' in quibble and re.search(r'\d+_', quibble):
        my_variant.warnings.append(
            "VariantSyntaxError: Base substitution (>) submitted with a "
            f"reference sequence range in {quibble}"
        )
        return True

    return False

def refseq_common_mistakes(variant, validator):
    """
    Check common reference sequence/type mistakes in unparsed HGVS input.

    Reference sequence/type mismatches are reported as ReferenceTypeErrors
    and terminate processing. Suggested corrected types are provided for
    guidance only; the submitted variant is not modified.

    Compound genomic/transcript descriptions such as NG_(NM_):c.,
    NC_(NR_):n., NT_(ENST):c., and NW_(LRG_t):c. are checked against
    the type of the embedded transcript and, when valid, allowed to
    continue to intronic_converter().

    ENST transcript type is determined from transcript CDS limits.
    LRG transcript references are treated as coding transcripts.

    Requires variant.quibble to be a string.
    """
    if not isinstance(variant.quibble, str):
        logger.error(
            "refseq_common_mistakes called with HGVS object %s",
            variant.quibble,
            stack_info=True
        )
        return False

    quibble = variant.quibble
    accession, _sep, _posedit = quibble.partition(':')

    acc3 = accession[:3]
    acc4 = accession[:4]

    logger.info(
        "Checking common reference sequence/type mistakes in %s",
        quibble
    )

    def type_suggestion(new_type):
        """
        Generate a corrected HGVS type without modifying the variant.
        """
        ref, sep, type_posedit = quibble.partition(':')
        _old_type, dot, posedit = type_posedit.partition('.')
        return f"{ref}{sep}{new_type}{dot}{posedit}"

    def transcript_reference_type(transcript_ref):
        """
        Determine whether a transcript reference is coding or non-coding.

        Returns 'c', 'n', or None when the transcript type cannot be
        determined here.
        """
        if transcript_ref.startswith('NM_'):
            return 'c'

        if transcript_ref.startswith('NR_'):
            return 'n'

        if transcript_ref.startswith('LRG_'):
            if (
                    transcript_ref[4:5].isdigit()
                    and 't' in transcript_ref
            ):
                return 'c'

            return None

        if transcript_ref.startswith('ENST'):
            try:
                tx_limits = validator.hdp.get_tx_limits(transcript_ref)
            except vvhgvs.exceptions.HGVSDataNotAvailableError:
                return None

            if tx_limits['cds_start_i'] is None:
                return 'n'

            return 'c'

        return None

    def genomic_transcript_error(transcript_ref, expected_type):
        """
        Report a mismatch between an embedded transcript and HGVS type.
        """
        suggestion = type_suggestion(expected_type)

        error = (
            'ReferenceTypeError: Transcript reference sequence '
            f'{transcript_ref} is not compatible with '
            f'{variant.reftype} variation. '
            f'Did you mean {suggestion}?'
        )
        variant.warnings.append(error)
        logger.info(error)

    # Genomic reference sequence.
    if acc3 in ('NC_', 'NG_', 'NT_', 'NW_'):
        _genomic_ref, separator, transcript_ref = accession.partition('(')

        # Genomic reference submitted as p.
        if variant.reftype == ':p.':
            error = (
                'ReferenceTypeError: Using a nucleotide reference sequence '
                '(NM_ NR_ NC_ NG_ NT_ NW_) to specify protein-level (p.) '
                'variation is not HGVS compliant. Please select an '
                'appropriate protein reference sequence (NP_)'
            )
            variant.warnings.append(error)
            logger.info(error)
            return True

        # Genomic reference submitted as c. or n.
        if variant.reftype in (':c.', ':n.'):
            if not separator:
                error = (
                    'ReferenceTypeError: Genomic reference sequence '
                    f'{variant.reftype} descriptions should not be used '
                    'unless a transcript reference sequence has also been '
                    'provided e.g. NG_(NM_):c.PositionVariation'
                )
                suggestion = (
                    f'For additional assistance, submit {variant.quibble} '
                    'to VariantValidator'
                )

                variant.warnings.extend([error, suggestion])
                logger.info(error)
                return True

            transcript_ref = transcript_ref.rstrip(')')
            transcript_type = transcript_reference_type(transcript_ref)

            # Unsupported or unavailable transcript information is left for
            # later validation to report appropriately.
            if transcript_type is None:
                return False

            requested_type = variant.reftype[1]

            if transcript_type != requested_type:
                genomic_transcript_error(
                    transcript_ref,
                    transcript_type
                )
                return True

            # Valid compound genomic/transcript descriptions are processed
            # later by intronic_converter().
            return False

        # Genomic DNA reference submitted as RNA.
        if variant.reftype == ':r.':
            error = (
                'ReferenceTypeError: Genomic DNA reference sequence input '
                'as RNA (r.) reference sequence.'
            )
            variant.warnings.append(error)
            logger.info(error)
            return True

        # NC_/NG_/NT_/NW_:g. is valid.
        return False

    # Determine the actual transcript type for simple transcript references.
    transcript_type = None

    if acc3 in ('NM_', 'NR_') or acc4 == 'ENST':
        transcript_type = transcript_reference_type(accession)

    # Transcript submitted as g.
    if (
            transcript_type in ('c', 'n')
            and variant.reftype == ':g.'
    ):
        suggestion = type_suggestion(transcript_type)

        error = (
            'ReferenceTypeError: Transcript reference sequence input as '
            'genomic (g.) reference sequence. '
            f'Did you mean {suggestion}?'
        )
        variant.warnings.append(error)
        logger.info(error)
        return True

    # Non-coding transcript submitted as c.
    if (
            transcript_type == 'n'
            and variant.reftype == ':c.'
    ):
        suggestion = type_suggestion('n')

        error = (
            'ReferenceTypeError: Non-coding transcript reference sequence '
            'input as coding (c.) reference sequence. '
            f'Did you mean {suggestion}?'
        )
        variant.warnings.append(error)
        logger.info(error)
        return True

    # NP_/ENSP submitted as a nucleotide variant.
    if (
            (acc3 == 'NP_' or acc4 == 'ENSP')
            and variant.reftype in (':c.', ':n.', ':g.', ':r.')
    ):
        error = (
            'ReferenceTypeError: Protein reference sequence input as '
            f'Nucleotide ({variant.reftype}) variant.'
        )
        variant.warnings.append(error)
        logger.info(error)
        return True

    # Coding transcript submitted as n.
    if (
            transcript_type == 'c'
            and variant.reftype == ':n.'
    ):
        suggestion = type_suggestion('c')

        error = (
            'ReferenceTypeError: Coding transcript reference sequence input '
            'as non-coding transcript (n.) reference sequence. '
            f'Did you mean {suggestion}?'
        )
        variant.warnings.append(error)
        logger.info(error)
        return True

    # Nucleotide transcript reference submitted as p.
    if (
            (
                acc3 in ('NM_', 'NR_')
                or acc4 == 'ENST'
            )
            and variant.reftype == ':p.'
    ):
        error = (
            'ReferenceTypeError: Using a nucleotide reference sequence '
            '(NM_ NR_ NC_ NG_ NT_ NW_) to specify protein-level (p.) '
            'variation is not HGVS compliant. Please select an appropriate '
            'protein reference sequence (NP_)'
        )
        variant.warnings.append(error)
        logger.info(error)
        return True

    return False


def refseq_type_mismatch(variant, validator):
    """
    Check reference sequence/type mismatches on a parsed HGVS variant.

    Reference sequence/type mismatches are reported as ReferenceTypeErrors
    and terminate processing. Suggested corrected types are provided for
    guidance only; the parsed HGVS object is not modified.

    ENST transcript type is determined from transcript CDS limits rather
    than relying on variant.transcript_type.

    This provides an object-stage safety net for variants that reach this
    point without having undergone the corresponding string-stage check.

    Requires variant.quibble to be an HGVS object.
    """
    if isinstance(variant.quibble, str):
        return False

    hgvs_variant = variant.quibble
    accession = hgvs_variant.ac
    variant_type = hgvs_variant.type

    acc3 = accession[:3]
    acc4 = accession[:4]

    def type_suggestion(new_type):
        """
        Generate a corrected type suggestion from the parsed HGVS object.
        """
        return f"{accession}:{new_type}.{hgvs_variant.posedit}"

    def transcript_reference_type(transcript_ref):
        """
        Determine whether a transcript reference is coding or non-coding.

        Returns 'c', 'n', or None when the transcript type cannot be
        determined here.
        """
        if transcript_ref.startswith('NM_'):
            return 'c'

        if transcript_ref.startswith('NR_'):
            return 'n'

        if transcript_ref.startswith('LRG_'):
            if (
                    transcript_ref[4:5].isdigit()
                    and 't' in transcript_ref
            ):
                return 'c'

            return None

        if transcript_ref.startswith('ENST'):
            try:
                tx_limits = validator.hdp.get_tx_limits(transcript_ref)
            except vvhgvs.exceptions.HGVSDataNotAvailableError:
                return None

            if tx_limits['cds_start_i'] is None:
                return 'n'

            return 'c'

        return None

    # Determine transcript type independently of variant.transcript_type.
    transcript_type = None

    if acc3 in ('NM_', 'NR_') or acc4 == 'ENST':
        transcript_type = transcript_reference_type(accession)

    # Transcript submitted as g.
    if (
            transcript_type in ('c', 'n')
            and variant_type == 'g'
    ):
        suggestion = type_suggestion(transcript_type)

        error = (
            'ReferenceTypeError: Transcript reference sequence input as '
            'genomic (g.) reference sequence. '
            f'Did you mean {suggestion}?'
        )
        variant.warnings.append(error)
        logger.info(error)
        return True

    # Non-coding transcript submitted as c.
    if (
            transcript_type == 'n'
            and variant_type == 'c'
    ):
        suggestion = type_suggestion('n')

        error = (
            'ReferenceTypeError: Non-coding transcript reference sequence '
            'input as coding (c.) reference sequence. '
            f'Did you mean {suggestion}?'
        )
        variant.warnings.append(error)
        logger.info(error)
        return True

    # NP_/ENSP submitted as a nucleotide variant.
    if (
            (acc3 == 'NP_' or acc4 == 'ENSP')
            and variant_type in ('c', 'n', 'g', 'r')
    ):
        error = (
            'ReferenceTypeError: Protein reference sequence input as '
            f'Nucleotide (:{variant_type}.) variant.'
        )
        variant.warnings.append(error)
        logger.info(error)
        return True

    # Coding transcript submitted as n.
    if (
            transcript_type == 'c'
            and variant_type == 'n'
    ):
        suggestion = type_suggestion('c')

        error = (
            'ReferenceTypeError: Coding transcript reference sequence input '
            'as non-coding transcript (n.) reference sequence. '
            f'Did you mean {suggestion}?'
        )
        variant.warnings.append(error)
        logger.info(error)
        return True

    # Nucleotide reference submitted as p.
    if (
            (
                acc3 in (
                    'NM_', 'NR_', 'NC_', 'NG_', 'NT_', 'NW_'
                )
                or acc4 == 'ENST'
            )
            and variant_type == 'p'
    ):
        error = (
            'ReferenceTypeError: Using a nucleotide reference sequence '
            '(NM_ NR_ NC_ NG_ NT_ NW_) to specify protein-level (p.) '
            'variation is not HGVS compliant. Please select an appropriate '
            'protein reference sequence (NP_)'
        )
        variant.warnings.append(error)
        logger.info(error)
        return True

    # Genomic reference submitted as c. or n. without transcript context.
    if (
            acc3 in ('NC_', 'NG_', 'NT_', 'NW_')
            and variant_type in ('c', 'n')
    ):
        error = (
            'ReferenceTypeError: Genomic reference sequence '
            f':{variant_type}. descriptions should not be used unless a '
            'transcript reference sequence has also been provided e.g. '
            'NG_(NM_):c.PositionVariation'
        )
        suggestion = (
            f'For additional assistance, submit {hgvs_variant} '
            'to VariantValidator'
        )

        variant.warnings.extend([error, suggestion])
        logger.info(error)
        return True

    # Genomic DNA reference submitted as RNA.
    if (
            acc3 in ('NC_', 'NG_', 'NT_', 'NW_')
            and variant_type == 'r'
    ):
        error = (
            'ReferenceTypeError: Genomic DNA reference sequence input '
            'as RNA (r.) reference sequence.'
        )
        variant.warnings.append(error)
        logger.info(error)
        return True

    return False

def structure_checks(variant, validator):
    """
    An evolving set of variant structure and content searches which identify
    and warn users about inappropriate use of HGVS.

    Primarily, this code filters out variants that cannot realistically be
    auto-corrected and will cause the downstream functions to return errors.
    """
    if not isinstance(variant.quibble, str):
        input_parses = copy.deepcopy(variant.quibble)
    else:
        input_parses = validator.hp.parse_hgvs_variant(
            variant.quibble
        )

    variant.input_parses = input_parses

    variant.gene_symbol = (
        validator.db.get_gene_symbol_from_transcript_id(
            variant.input_parses.ac
        )
    )

    if variant.gene_symbol == 'none':
        variant.gene_symbol = ''

    if input_parses.type in ('g', 'm'):
        if structure_checks_g(variant, validator):
            return True

    elif input_parses.type == 'c':
        if structure_checks_c(variant, validator):

            # Also check intron boundaries to provide additional warnings.
            if (
                    "beyond the bounds" in str(variant.warnings)
                    and hgvs_position_utils.either_position_is_intronic(
                        variant.input_parses
                    )
            ):
                if hgvs_position_utils.start_position_is_intronic(
                        variant.input_parses
                ):
                    variant.input_parses.posedit.pos.start.offset = 1

                if hgvs_position_utils.end_position_is_intronic(
                        variant.input_parses
                ):
                    variant.input_parses.posedit.pos.end.offset = 1

                if variant.genomic_context_ac is not None:
                    hgvs_genomic_vt = validator.vm.t_to_g(
                        variant.input_parses,
                        variant.genomic_context_ac
                    )
                else:
                    hgvs_genomic_vt = validator.myevm_t_to_g(
                        variant.input_parses,
                        variant.no_norm_evm,
                        variant.primary_assembly,
                        variant.hn,
                        variant
                    )

                format_converters.remap_intronic(
                    variant.input_parses,
                    hgvs_genomic_vt,
                    variant,
                    validator
                )

            return True

    elif input_parses.type == 'n':
        if structure_checks_n(variant, validator):
            return True

    return False


def structure_checks_g(variant, validator):
    """
    Structure checks for when reftype is genomic.
    """
    if not variant.input_parses.ac.startswith(
            ('NC_', 'NG_', 'NT_', 'NW_')
    ):
        error = (
            f'Invalid reference sequence identifier '
            f'({variant.input_parses.ac})'
        )
        variant.warnings.append(error)
        logger.info(error)
        return True

    try:
        validator.vr.validate(variant.input_parses)

    except Exception as e:
        error = str(e)

        if "does not agree with reference sequence ()" in error:
            error = (
                f"The specified coordinate is outside the boundaries of "
                f"reference sequence {variant.input_parses.ac} "
            )

        if (
                "base start position must be <= end position" in error
                and variant.hgvs_formatted.ac in (
                    "NC_012920.1",
                    "NC_001807.4"
                )
        ):
            if variant.hgvs_formatted.ac not in variant.original:
                warning = (
                    "This is not a valid HGVS variant description, because "
                    "no reference sequence ID has been provided, instead use "
                    f"{variant.hgvs_formatted}"
                )
                variant.warnings.append(warning)

            variant.warnings.append(
                "The variant positions are valid but we cannot normalize "
                "variants spanning the origin of circular reference sequences"
            )
            return True

        elif (
                "insertion length must be 1" in error
                and "(" in str(variant.input_parses.posedit.pos)
                and ")" in str(variant.input_parses.posedit.pos)
        ):
            return True

        elif (
                "insertion length must be 1" in error
                and "(" not in str(variant.input_parses.posedit.pos)
                and ")" not in str(variant.input_parses.posedit.pos)
        ):
            start = variant.input_parses.posedit.pos.start.base

            ins_warning = (
                f'Insertion length must be 1 e.g. '
                f'{start}_{start + 1}'
                f'ins{variant.input_parses.posedit.edit.alt}'
            )

            variant.warnings.append(ins_warning)

            if "insertion length must be 1" in variant.warnings:
                variant.warnings.remove(
                    "insertion length must be 1"
                )

            return True

        elif (
                "Length implied by coordinates must equal sequence "
                "deletion length" in error
                and "(" in str(variant.input_parses.posedit.pos)
                and ")" in str(variant.input_parses.posedit.pos)
        ):
            return True

        else:
            variant.warnings.append(error)
            logger.info(error)
            return True

    # Additional normalization test.
    try:
        variant.hn.normalize(variant.input_parses)

    except vvhgvs.exceptions.HGVSError as e:
        error = str(e)
        variant.warnings.append(error)
        logger.info(error)
        return True

    # Look for variants in runs of N bases.
    ref = getattr(
        variant.input_parses.posedit.edit,
        'ref',
        None
    )

    if ref is not None and "N" in ref:
        error = (
            f"UncertainSequenceError: The submitted variant description "
            f"{variant.input_parses} refers to a genomic reference region "
            f"with an uncertain base composition (N)"
        )
        variant.warnings.append(error)
        logger.info(error)
        return True

    return False

def structure_checks_c(variant, validator):
    """
    Structure checks for when reftype is coding.

    :param variant:
    :param validator:
    :return:
    """

    # Catch variation in UTRs.
    # These should be in the sequence so can be directly validated.
    if (
            hgvs_position_utils.start_is_5_prime_utr(variant.input_parses)
            or hgvs_position_utils.end_is_3_prime_utr(variant.input_parses)
    ):
        logger.info(
            f"Check datum definition for variant {variant.input_parses}"
        )

        try:
            validator.vr.validate(variant.input_parses)

        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)

            if 'datums is ill-defined' in error:
                logger.info(
                    f"The specified datum is ill-defined for variant "
                    f"{variant.input_parses}"
                )

                called_ref = variant.input_parses.posedit.edit.ref

                try:
                    to_n = variant.evm.c_to_n(variant.input_parses)
                except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                    error = str(e)
                    variant.warnings.append(error)
                    logger.info(error)
                    return True

                actual_ref = to_n.posedit.edit.ref

                if called_ref != actual_ref:
                    error = (
                        f'Variant reference ({called_ref}) does not agree '
                        f'with reference sequence ({actual_ref})'
                    )
                    variant.warnings.append(error)
                    logger.info(error)
                    return True

                if variant.input_parses.posedit.edit.type == "ins":
                    variant.input_parses.posedit.edit.ref = None
                else:
                    variant.input_parses.posedit.edit.ref = ''

                variant.hgvs_formatted = variant.input_parses

            elif 'bounds' in error or 'intronic variant' in error:
                try:
                    variant.hn.normalize(variant.input_parses)
                except vvhgvs.exceptions.HGVSError as e:
                    logger.debug("Except passed, %s", e)

                if 'bounds' in error:
                    try:
                        identity_info = validator.hdp.get_tx_identity_info(
                            variant.input_parses.ac
                        )
                        ref_start = identity_info[3]
                        ref_end = identity_info[4]

                        if (
                                hgvs_position_utils.start_is_5_prime_utr(
                                    variant.input_parses
                                )
                                and not hgvs_position_utils.start_position_is_intronic(
                                    variant.input_parses
                                )
                        ):
                            boundary = -ref_start
                            remainder = (
                                variant.input_parses.posedit.pos.start.base
                                - boundary
                            )
                            variant.input_parses.posedit.pos.start.base = boundary
                            variant.input_parses.posedit.pos.start.offset = remainder

                        if (
                                hgvs_position_utils.end_is_5_prime_utr(
                                    variant.input_parses
                                )
                                and not hgvs_position_utils.end_position_is_intronic(
                                    variant.input_parses
                                )
                        ):
                            boundary = -ref_start
                            remainder = (
                                variant.input_parses.posedit.pos.end.base
                                - boundary
                            )
                            variant.input_parses.posedit.pos.end.base = boundary
                            variant.input_parses.posedit.pos.end.offset = remainder

                        if (
                                hgvs_position_utils.start_is_3_prime_utr(
                                    variant.input_parses
                                )
                                and not hgvs_position_utils.start_position_is_intronic(
                                    variant.input_parses
                                )
                        ):
                            tot_end_pos = (
                                variant.input_parses.posedit.pos.start.base
                            )
                            ts_seq = validator.sf.fetch_seq(
                                variant.input_parses.ac
                            )
                            boundary = len(ts_seq) - ref_end
                            variant.input_parses.posedit.pos.start.base = boundary
                            offset = tot_end_pos - boundary
                            variant.input_parses.posedit.pos.start.offset = offset

                        if (
                                hgvs_position_utils.end_is_3_prime_utr(
                                    variant.input_parses
                                )
                                and not hgvs_position_utils.end_position_is_intronic(
                                    variant.input_parses
                                )
                        ):
                            tot_end_pos = (
                                variant.input_parses.posedit.pos.end.base
                            )
                            ts_seq = validator.sf.fetch_seq(
                                variant.input_parses.ac
                            )
                            boundary = len(ts_seq) - ref_end
                            variant.input_parses.posedit.pos.end.base = boundary
                            offset = tot_end_pos - boundary
                            variant.input_parses.posedit.pos.end.offset = offset

                        # Create a loose vm instance.
                        variant.lose_vm = vvhgvs.variantmapper.VariantMapper(
                            validator.hdp,
                            replace_reference=True,
                            prevalidation_level=None
                        )

                        report_gen = validator.myevm_t_to_g(
                            variant.input_parses,
                            variant.no_norm_evm,
                            variant.primary_assembly,
                            variant.hn,
                            variant
                        )
                        report_gen = variant.hn.normalize(report_gen)

                        error = (
                            'Using a transcript reference sequence to specify '
                            'a variant position that lies outside of the '
                            'reference sequence is not HGVS-compliant: '
                            'Instead re-submit '
                            + fn.valstr(report_gen)
                        )

                    except Exception as e:
                        logger.debug("Except passed, %s", e)

                    variant.warnings.append(error)
                    logger.info(error)
                    return True

        try:
            variant.input_parses = variant.evm.c_to_n(
                variant.input_parses
            )
        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.info(e)
            return True

        if (
                variant.input_parses.posedit.pos.start.base == 1
                and hgvs_position_utils.start_offset_is_negative(
                    variant.input_parses
                )
        ):
            input_parses = variant.evm.n_to_c(
                variant.input_parses
            )

            error = (
                'Using a transcript reference sequence to specify a variant '
                'position that lies outside of the reference sequence is not '
                'HGVS-compliant. Instead re-submit '
            )

            genomic_position = validator.myevm_t_to_g(
                input_parses,
                variant.no_norm_evm,
                variant.primary_assembly,
                variant.hn,
                variant
            )
            genomic_position = variant.hn.normalize(genomic_position)

            error += fn.valstr(genomic_position)
            variant.warnings.append(error)
            logger.info(error)
            return True

        # Re-map input_parses back to c. variant.
        variant.input_parses = variant.evm.n_to_c(
            variant.input_parses
        )

        # Intronic positions in UTRs.
        if hgvs_position_utils.either_position_is_intronic(
                variant.input_parses
        ):
            logger.info(
                f"Check UTR boundaries for variant {variant.input_parses}"
            )

            try:
                to_genome = validator.myevm_t_to_g(
                    variant.input_parses,
                    variant.no_norm_evm,
                    variant.primary_assembly,
                    variant.hn,
                    variant
                )
                variant.evm.g_to_t(
                    to_genome,
                    variant.input_parses.ac
                )

            except vvhgvs.exceptions.HGVSInvalidIntervalError as e:
                error = str(e)

                logger.info(
                    f"Check UTR boundaries for variant "
                    f"{variant.input_parses}: error {error}"
                )

                if 'bounds' in error:
                    try:
                        identity_info = validator.hdp.get_tx_identity_info(
                            variant.input_parses.ac
                        )
                        ref_start = identity_info[3]
                        ref_end = identity_info[4]

                        if hgvs_position_utils.start_is_5_prime_utr(
                                variant.input_parses
                        ):
                            boundary = -ref_start
                            remainder = (
                                variant.input_parses.posedit.pos.start.base
                                - boundary
                            )
                            variant.input_parses.posedit.pos.start.base = boundary
                            variant.input_parses.posedit.pos.start.offset = remainder

                        if hgvs_position_utils.end_is_5_prime_utr(
                                variant.input_parses
                        ):
                            boundary = -ref_start
                            remainder = (
                                variant.input_parses.posedit.pos.end.base
                                - boundary
                            )
                            variant.input_parses.posedit.pos.end.base = boundary
                            variant.input_parses.posedit.pos.end.offset = remainder

                        if hgvs_position_utils.start_is_3_prime_utr(
                                variant.input_parses
                        ):
                            tot_end_pos = (
                                variant.input_parses.posedit.pos.start.base
                                + variant.input_parses.posedit.pos.start.offset
                            )
                            ts_seq = validator.sf.fetch_seq(
                                variant.input_parses.ac
                            )
                            boundary = len(ts_seq) - ref_end
                            variant.input_parses.posedit.pos.start.base = boundary
                            offset = tot_end_pos - boundary
                            variant.input_parses.posedit.pos.start.offset = offset

                        if hgvs_position_utils.end_is_3_prime_utr(
                                variant.input_parses
                        ):
                            tot_end_pos = (
                                variant.input_parses.posedit.pos.end.base
                                + variant.input_parses.posedit.pos.end.offset
                            )
                            ts_seq = validator.sf.fetch_seq(
                                variant.input_parses.ac
                            )
                            boundary = len(ts_seq) - ref_end
                            variant.input_parses.posedit.pos.end.base = boundary
                            offset = tot_end_pos - boundary
                            variant.input_parses.posedit.pos.end.offset = offset

                        report_gen = validator.myevm_t_to_g(
                            variant.input_parses,
                            variant.no_norm_evm,
                            variant.primary_assembly,
                            variant.hn,
                            variant
                        )
                        report_gen = variant.hn.normalize(report_gen)

                        error = (
                            'Using a transcript reference sequence to specify '
                            'a variant position that lies outside of the '
                            'reference sequence is not HGVS-compliant. '
                            'Instead re-submit '
                            + fn.valstr(report_gen)
                        )

                    except Exception as e:
                        logger.debug("Except passed, %s", e)

                variant.warnings.append(error)
                logger.info(error)
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

                    error = (
                        'Cannot map '
                        + fn.valstr(variant.input_parses)
                        + ' to a genomic position. '
                        + variant.input_parses.ac
                        + ' can only be partially aligned to genomic reference '
                        'sequences '
                        + acs
                    )

                    variant.warnings.append(error)
                    logger.info(error)
                    return True

    elif hgvs_position_utils.either_position_is_intronic(
            variant.input_parses
    ):
        # Create a no_replace vm instance.
        variant.no_replace_vm = vvhgvs.variantmapper.VariantMapper(
            validator.hdp,
            replace_reference=False,
            prevalidation_level=None
        )

        # Quick look at syntax validation.
        try:
            validator.vr.validate(variant.input_parses)

        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)

            if 'bounds' in error:
                try:
                    report_gen = validator.myevm_t_to_g(
                        variant.input_parses,
                        variant.no_norm_evm,
                        variant.primary_assembly,
                        variant.hn,
                        variant
                    )
                    report_gen = variant.hn.normalize(report_gen)

                except vvhgvs.exceptions.HGVSError as e:
                    logger.debug("Except passed, %s", e)

                else:
                    error = (
                        'Using a transcript reference sequence to specify a '
                        'variant position that lies outside of the reference '
                        'sequence is not HGVS-compliant. Instead re-submit '
                        + fn.valstr(report_gen)
                    )

                variant.warnings.append(error)
                logger.info(error)
                return True

            elif 'insertion length must be 1' in error:
                if hgvs_position_utils.start_offset_is_negative(
                        variant.input_parses
                ):
                    start_offset = str(
                        variant.input_parses.posedit.pos.start.offset - 1
                    )
                    end_offset = str(
                        variant.input_parses.posedit.pos.start.offset
                    )

                elif hgvs_position_utils.start_offset_is_positive(
                        variant.input_parses
                ):
                    start_offset = (
                        f"+{variant.input_parses.posedit.pos.start.offset}"
                    )
                    end_offset = (
                        f"+{variant.input_parses.posedit.pos.start.offset + 1}"
                    )

                if "(" not in str(variant.input_parses.posedit.pos):
                    ins_warning = (
                        f'Insertion length must be 1 e.g. '
                        f'{variant.input_parses.posedit.pos.start.base}'
                        f'{start_offset}_'
                        f'{int(variant.input_parses.posedit.pos.start.base)}'
                        f'{end_offset}'
                        f'ins{variant.input_parses.posedit.edit.alt}'
                    )

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

                error = (
                    error
                    + ': Did you mean '
                    + str(correction)
                    + '?'
                )

                variant.warnings.append(error)
                logger.info(error)
                return True

        # Create a specific minimal evm with no normalizer and no
        # replace_reference.
        output = None

        try:
            if variant.genomic_context_ac is not None:
                output = validator.vm.t_to_g(
                    variant.input_parses,
                    variant.genomic_context_ac
                )
            else:
                output = validator.noreplace_myevm_t_to_g(
                    variant.input_parses,
                    variant
                )

        except vvhgvs.exceptions.HGVSDataNotAvailableError:
            tx_info = validator.hdp.get_tx_identity_info(
                variant.input_parses.ac
            )

            if (
                    (
                        variant.input_parses.posedit.pos.start.base
                        > int(tx_info[4])
                        or variant.input_parses.posedit.pos.end.base
                        > int(tx_info[4])
                    )
                    and not (
                        hgvs_position_utils.start_is_3_prime_utr(
                            variant.input_parses
                        )
                        or hgvs_position_utils.end_is_3_prime_utr(
                            variant.input_parses
                        )
                    )
            ):
                errors = [
                    "CDSError: Variant start position and/or end position "
                    "are beyond the CDS end position and likely also beyond "
                    "the end of the selected reference sequence"
                ]
            else:
                errors = [
                    'Required information for '
                    + variant.input_parses.ac
                    + ' is missing from the Universal Transcript Archive',
                    'Query gene2transcripts with search term %s for '
                    'available transcripts'
                    % variant.input_parses.ac.split('.')[0]
                ]

            variant.warnings.extend(errors)
            logger.info(str(errors))
            return True

        except ValueError as e:
            error = str(e)

            if '> end' in error:
                error = (
                    'Interval start position '
                    + str(variant.input_parses.posedit.pos.start)
                    + ' > interval end position '
                    + str(variant.input_parses.posedit.pos.end)
                )

                variant.warnings.append(error)
                logger.info(error)
                return True

        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)

            if 'base start position must be <= end position' in error:
                error = (
                    'Interval start position '
                    + str(variant.input_parses.posedit.pos.start)
                    + ' > interval end position '
                    + str(variant.input_parses.posedit.pos.end)
                )

            variant.warnings.append(error)
            logger.info(error)
            return True

        try:
            variant.evm.g_to_t(
                output,
                variant.input_parses.ac
            )

        except vvhgvs.exceptions.HGVSError as e:
            if "Alignment is incomplete" in str(e):
                output = hgvs_utils.incomplete_alignment_mapping_t_to_g(
                    validator,
                    variant
                )

                if output is None:
                    error = str(e)
                    variant.warnings.append(error)
                    logger.info(error)
                    return True

            else:
                error = str(e)
                variant.warnings.append(error)
                logger.info(error)
                return True

        # Check that the reference is correct by direct mapping without
        # replacing reference.
        check_ref_g = variant.no_replace_vm.t_to_g(
            variant.input_parses,
            output.ac,
            alt_aln_method=validator.alt_aln_method
        )

        logger.info(
            f"Output {output} mapped to {check_ref_g}"
        )

        check_ref_t = variant.no_replace_vm.g_to_t(
            check_ref_g,
            variant.input_parses.ac,
            alt_aln_method=validator.alt_aln_method
        )

        logger.info(
            f"{check_ref_g} mapped to {check_ref_t}"
        )

        # Snapshot current variant error log.
        if (
                hgvs_position_utils.end_is_3_prime_utr(check_ref_t)
                and not hgvs_position_utils.end_is_3_prime_utr(
                    variant.input_parses
                )
        ):
            convert = (
                f"{variant.input_parses} auto-mapped to {check_ref_t}"
            )
            variant.warnings.append(convert)

        snap = copy.copy(variant.warnings)

        # Look for syntax errors.
        try:
            validator.vr.validate(check_ref_t)
        except vvhgvs.exceptions.HGVSError as e:
            if (
                    "intron" not in str(e)
                    and "bounds" not in str(e)
                    and "insertion length must be 1" not in str(e)
                    and "base start position must be <= end position"
                    not in str(e)
            ):
                error = str(e)
                variant.warnings.append(error)
                logger.info(error)

        try:
            validator.vr.validate(output)
        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.info(error)

        # Check for additional warnings.
        if len(variant.warnings) > len(snap):
            return True

    else:
        # All other variation.
        try:
            validator.vr.validate(variant.input_parses)

        except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
            logger.debug("Except passed, %s", e)

        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)

            if 'base start position must be <= end position' in error:
                error = (
                    'Interval start position '
                    + str(variant.input_parses.posedit.pos.start)
                    + ' > interval end position '
                    + str(variant.input_parses.posedit.pos.end)
                )

            if "(" not in str(variant.input_parses.posedit.pos):
                variant.warnings.append(error)
                logger.info(error)

            if 'insertion length must be 1' in error:
                if "(" not in str(variant.input_parses.posedit.pos):
                    ins_warning = (
                        f'Insertion length must be 1 e.g. '
                        f'{int(variant.input_parses.posedit.pos.start.base)}_'
                        f'{int(variant.input_parses.posedit.pos.start.base) + 1}'
                        f'ins{variant.input_parses.posedit.edit.alt}'
                    )

                    variant.warnings.append(ins_warning)

                    for warning in variant.warnings:
                        if warning == "insertion length must be 1":
                            variant.warnings.remove(warning)

            return True

        except vvhgvs.exceptions.HGVSDataNotAvailableError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.info(error)
            return True

        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)

            if 'bounds' in error:
                error += f' ({variant.input_parses.ac})'
                variant.warnings.append(error)
                logger.info(error)
                return True

    return False


def structure_checks_n(variant, validator):
    """
    Structure checks for reftype nucleotide.

    :param variant:
    :param validator:
    :return:
    """

    # Intronic n. variants.
    if hgvs_position_utils.either_position_is_intronic(
            variant.input_parses
    ):
        # Quick look at syntax validation.
        try:
            validator.vr.validate(variant.input_parses)

        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)

            if 'intronic variant' in error:
                pass

            elif 'datums is ill-defined' in error:
                called_ref = variant.input_parses.posedit.edit.ref

                # input_parses is already an n. variant, so do not convert
                # c. -> n. here. Validate the supplied reference against the
                # transcript sequence via the existing mapping path.
                try:
                    if variant.genomic_context_ac is not None:
                        to_genome = validator.vm.t_to_g(
                            variant.input_parses,
                            variant.genomic_context_ac
                        )
                    else:
                        to_genome = validator.noreplace_myevm_t_to_g(
                            variant.input_parses,
                            variant
                        )

                    to_n = validator.vm.g_to_t(
                        to_genome,
                        variant.input_parses.ac
                    )

                    actual_ref = to_n.posedit.edit.ref

                except vvhgvs.exceptions.HGVSError as mapping_error:
                    error = str(mapping_error)
                    variant.warnings.append(error)
                    logger.info(error)
                    return True

                if called_ref != actual_ref:
                    error = (
                        f'Variant reference ({called_ref}) does not agree '
                        f'with reference sequence ({actual_ref})'
                    )
                    variant.warnings.append(error)
                    logger.info(error)
                    return True

                variant.input_parses.posedit.edit.ref = ''
                variant.hgvs_formatted = variant.input_parses

            elif (
                    'base must be >=1 for datum = SEQ_START or CDS_END'
                    in error
            ):
                error = (
                    'The given coordinate is outside the bounds of the '
                    'reference sequence.'
                )

                try:
                    if (
                            variant.input_parses.posedit.pos.start.base < 1
                            and not hgvs_position_utils.start_position_is_intronic(
                                variant.input_parses
                            )
                    ):
                        boundary = 1
                        remainder = (
                            variant.input_parses.posedit.pos.start.base
                            - boundary
                        )
                        remainder += 1
                        variant.input_parses.posedit.pos.start.base = boundary
                        variant.input_parses.posedit.pos.start.offset = remainder

                    if (
                            variant.input_parses.posedit.pos.end.base < 1
                            and not hgvs_position_utils.end_position_is_intronic(
                                variant.input_parses
                            )
                    ):
                        boundary = 1
                        remainder = (
                            variant.input_parses.posedit.pos.end.base
                            - boundary
                        )
                        remainder += 1
                        variant.input_parses.posedit.pos.end.base = boundary
                        variant.input_parses.posedit.pos.end.offset = remainder

                    report_gen = validator.myevm_t_to_g(
                        variant.input_parses,
                        variant.no_norm_evm,
                        variant.primary_assembly,
                        variant.hn,
                        variant
                    )
                    report_gen = variant.hn.normalize(report_gen)

                    error = (
                        'Using a transcript reference sequence to specify a '
                        'variant position that lies outside of the reference '
                        'sequence is not HGVS-compliant. Instead re-submit '
                        + fn.valstr(report_gen)
                    )

                except Exception as e:
                    logger.debug("Except passed, %s", e)

                variant.warnings.append(error)
                logger.info(error)
                return True

            elif 'Cannot validate sequence of an intronic variant' in error:
                try:
                    if variant.genomic_context_ac is not None:
                        test_g = validator.vm.t_to_g(
                            variant.input_parses,
                            variant.genomic_context_ac
                        )
                    else:
                        test_g = validator.myevm_t_to_g(
                            variant.input_parses,
                            variant.no_norm_evm,
                            variant.primary_assembly,
                            variant.hn,
                            variant
                        )

                    variant.evm.g_to_t(
                        test_g,
                        variant.input_parses.ac
                    )

                except vvhgvs.exceptions.HGVSError as mapping_error:
                    error = str(mapping_error)

                    if 'bounds' in error:
                        if variant.genomic_context_ac is not None:
                            report_gen = validator.vm.t_to_g(
                                variant.input_parses,
                                variant.genomic_context_ac
                            )
                        else:
                            report_gen = validator.myevm_t_to_g(
                                variant.input_parses,
                                variant.no_norm_evm,
                                variant.primary_assembly,
                                variant.hn,
                                variant
                            )

                        report_gen = variant.hn.normalize(report_gen)

                        error = (
                            'Using a transcript reference sequence to specify '
                            'a variant position that lies outside of the '
                            'reference sequence is not HGVS-compliant. '
                            'Instead re-submit '
                            + fn.valstr(report_gen)
                        )

                        variant.warnings.append(error)
                        logger.info(error)
                        return True

            else:
                variant.warnings.append(error)
                logger.info(error)
                return True

        # n.1-x is outside the transcript reference sequence.
        if (
                variant.input_parses.posedit.pos.start.base == 1
                and hgvs_position_utils.start_offset_is_negative(
                    variant.input_parses
                )
        ):
            error = (
                'Using a transcript reference sequence to specify a variant '
                'position that lies outside of the reference sequence is not '
                'HGVS-compliant. Instead re-submit '
            )

            if variant.genomic_context_ac is not None:
                genomic_position = validator.vm.t_to_g(
                    variant.input_parses,
                    variant.genomic_context_ac
                )
            else:
                genomic_position = validator.myevm_t_to_g(
                    variant.input_parses,
                    variant.no_norm_evm,
                    variant.primary_assembly,
                    variant.hn,
                    variant
                )

            genomic_position = variant.hn.normalize(
                genomic_position
            )

            error += fn.valstr(genomic_position)
            variant.warnings.append(error)
            logger.info(error)
            return True

        # Create a specific minimal mapping with no normalizer and no
        # replace_reference.
        output = None

        try:
            if variant.genomic_context_ac is not None:
                output = validator.vm.t_to_g(
                    variant.input_parses,
                    variant.genomic_context_ac
                )
            else:
                output = validator.noreplace_myevm_t_to_g(
                    variant.input_parses,
                    variant
                )

        except vvhgvs.exceptions.HGVSDataNotAvailableError:
            errors = [
                'Required information for '
                + variant.input_parses.ac
                + ' is missing from the Universal Transcript Archive',
                'Query gene2transcripts with search term %s for '
                'available transcripts'
                % variant.input_parses.ac.split('.')[0]
            ]
            variant.warnings.extend(errors)
            logger.info(str(errors))
            return True

        except ValueError as e:
            error = str(e)

            if '> end' in error:
                error = (
                    'Interval start position '
                    + str(variant.input_parses.posedit.pos.start)
                    + ' > interval end position '
                    + str(variant.input_parses.posedit.pos.end)
                )
                variant.warnings.append(error)
                logger.info(error)
                return True

        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)

            if 'base start position must be <= end position' in error:
                error = (
                    'Interval start position '
                    + str(variant.input_parses.posedit.pos.start)
                    + ' > interval end position '
                    + str(variant.input_parses.posedit.pos.end)
                )

            variant.warnings.append(error)
            logger.info(error)
            return True

        try:
            validator.vr.validate(output)

        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.info(error)
            return True

    else:
        # All other variation.
        try:
            validator.vr.validate(variant.input_parses)

        except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
            logger.debug("Except passed, %s", e)

        except vvhgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)

            if 'base start position must be <= end position' in error:
                error = (
                    'Interval start position '
                    + str(variant.input_parses.posedit.pos.start)
                    + ' > interval end position '
                    + str(variant.input_parses.posedit.pos.end)
                )

            variant.warnings.append(error)
            logger.info(error)
            return True

        except vvhgvs.exceptions.HGVSDataNotAvailableError as e:
            error = str(e)
            variant.warnings.append(error)
            logger.info(error)
            return True

        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)

            if 'bounds' in error:
                error += f' ({variant.input_parses.ac})'
                variant.warnings.append(error)
                logger.info(error)
                return True

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
