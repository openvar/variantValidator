import re
import copy
from VariantValidator.modules import utils as vv_utils, format_converters, hgvs_position_utils
from VariantValidator.modules.hgvs_utils import hgvs_obj_from_existing_edit
import vvhgvs.exceptions
from vvhgvs.location import Interval
from vvhgvs.enums import ValidationLevel
from vvhgvs.exceptions import HGVSUnsupportedOperationError
import logging
logger = logging.getLogger(__name__)

class UncertainConversionError(Exception):
    pass

# fuzzy but specified ended intervals, take a normal interval for both ends
# can take a pair of Intervals or BaseOffsetIntervals
class FEInterval(Interval):
    uncertain = True
    pass
    def validate(self):
        if self.start:
            (res, msg) = self.start.validate()
            if res != ValidationLevel.VALID:
                logger.info("Validation failed for start interval {}: {}".format(self.start, msg))
                return (res, msg)
        if self.end:
            (res, msg) = self.end.validate()
            if res != ValidationLevel.VALID:
                logger.info("Validation failed for end interval {}: {}".format(self.end, msg))
                return (res, msg)
        # Check start less than or equal to end
        # for now overlap is allowed so long as some is distinct
        if not self.start or not self.end:
            logger.info("No start or end interval specified")
            return (ValidationLevel.VALID, None)
        try:
            if self.start.start <= self.end.end:
                return (ValidationLevel.VALID, None)
            else:
                return (ValidationLevel.ERROR, "base start position must be <= end position")
        except HGVSUnsupportedOperationError as err:
            logger.info("HGVSUnsupportedOperationError {}".format(err))
            return (ValidationLevel.WARNING, str(err))

    def format(self, conf=None):
        if self.start is None:
            return ""
        if self.end is None or self.start == self.end:
            return self.start.format(conf)
        # always uncertain
        iv = "(" +self.start.format(conf) + ")_(" + self.end.format(conf) + ")"
        return iv



# Errors
class FuzzyPositionError(Exception):
    pass


class FuzzyRangeError(Exception):
    pass


class InvalidRangeError(Exception):
    pass


class IncompatibleTypeError(Exception):
    pass


class HgvsParseError(Exception):
    pass


def fuzzy_ends(my_variant, validator):
    """
    Detect fuzzy/unknown HGVS positions.

    Fuzzy range syntax that cannot be represented directly by the HGVS parser
    is handled as a string. Where an HGVS object is already available, inspect
    the structured position object rather than converting it back to a string.

    Returns False if no fuzzy end is detected; otherwise raises the appropriate
    fuzzy-position exception.
    """
    quibble = my_variant.quibble

    # This route operates on unsupported fuzzy syntax, so quibble is expected
    # to be a string here.
    accession, _sep, type_posedit = quibble.partition(':')

    if accession.startswith(('NC_', 'NG_')):
        format_converters.intronic_converter(
            my_variant,
            validator,
            uncertain=True
        )
        quibble = my_variant.quibble
        accession, _sep, type_posedit = quibble.partition(':')

    # Cache the variant type for the exon-boundary checks below.
    var_type = type_posedit[:1]

    def parse_position(position, intronic_positions):
        """
        Convert a fuzzy-range position to its base coordinate.

        Intronic positions retain their exon-boundary base so it can be
        validated later.
        """
        try:
            return int(position)
        except ValueError:
            pass

        if '+' in position:
            base = position.partition('+')[0]
            try:
                base = int(base)
                intronic_positions.append(base)
            except ValueError:
                pass
            return base

        if '-' in position:
            base = position.partition('-')[0]
            try:
                base = int(base)
                intronic_positions.append(base)
            except ValueError:
                pass
            return base

        return position

    def check_intronic_boundaries(intronic_positions):
        if var_type not in ('c', 'n', 'r') or not intronic_positions:
            return

        exon_boundaries = vv_utils.get_exon_boundary_list(
            my_variant,
            validator
        )

        for pos in intronic_positions:
            if str(pos) not in exon_boundaries:
                raise FuzzyPositionError(
                    f"ExonBoundaryError: Position {pos} does not correspond "
                    f"to an exon boundary for transcript {accession} aligned "
                    f"to {my_variant.primary_assembly} genomic reference "
                    f"sequence {exon_boundaries[1]}"
                )

    tx_info = None

    def utr_to_int(position):
        nonlocal tx_info

        if isinstance(position, int):
            return position

        if '*' not in position:
            return int(position)

        if tx_info is None:
            tx_info = validator.hdp.get_tx_identity_info(accession)

        return int(position.replace('*', '')) + int(tx_info[4]) - 1

    # Unknown start and unknown end:
    #
    # (?_100)_(200_?)del
    if re.match(
            r"(NM_|NR_|NC_|NG_|ENST)\d+\.\d+:(g|c)\."
            r"\(\?\_(\+|-|\*)?\d+(\+|-)?\d*\)"
            r"\_\((\+|-|\*)?\d+(\+|-)?\d*_\?\)"
            r"(del|dup|inv)$",
            quibble
    ):
        parts = quibble.split(')_(', 1)
        num1 = parts[0].rsplit('_', 1)[1]
        num2 = parts[1].partition('_')[0]

        logger.info("fuzzy end detected %s %s", num1, num2)

        intronic_positions = []

        num1 = parse_position(num1, intronic_positions)
        num2 = parse_position(num2, intronic_positions)

        check_intronic_boundaries(intronic_positions)

        num1 = utr_to_int(num1)
        num2 = utr_to_int(num2)

        if num1 >= num2:
            my_variant.warnings.append(
                "Uncertain positions are not fully supported, however the "
                "start position is > the end position"
            )
        else:
            my_variant.warnings.append(
                "Uncertain positions are not fully supported, however the "
                "syntax is valid"
            )

        raise FuzzyRangeError(
            "Fuzzy/unknown variant start and end positions "
            "in submitted variant description"
        )

    # Known start range, unknown end:
    #
    # (100_200)_(300_?)del
    elif re.match(
            r"(NM_|NR_|NC_|NG_|ENST)\d+\.\d+:(g|c)\."
            r"\((\+|-|\*)?\d+(\+|-)?\d*_(\+|-|\*)?\d+(\+|-)?\d*\)"
            r"\_\((\+|-|\*)?\d+(\+|-)?\d*_\?\)"
            r"(del|dup|inv)$",
            quibble
    ):
        parts = quibble.split(')_(', 1)

        first_range = parts[0].partition('(')[2]
        num1, num2 = first_range.split('_', 1)
        num3 = parts[1].partition('_')[0]

        intronic_positions = []

        num1 = parse_position(num1, intronic_positions)
        num2 = parse_position(num2, intronic_positions)
        num3 = parse_position(num3, intronic_positions)

        check_intronic_boundaries(intronic_positions)

        num1 = utr_to_int(num1)
        num2 = utr_to_int(num2)
        num3 = utr_to_int(num3)

        if num1 <= num2 <= num3:
            my_variant.warnings.append(
                "Uncertain positions are not fully supported, however the "
                "syntax is valid"
            )
        else:
            my_variant.warnings.append(
                "Uncertain positions are not fully supported, however the "
                "provided positions are out of order"
            )

        raise FuzzyRangeError(
            "Fuzzy/unknown variant end position in submitted variant "
            "description"
        )

    # Unknown start, known end range:
    #
    # (?_100)_(200_300)del
    elif re.match(
            r"(NM_|NR_|NC_|NG_|ENST)\d+\.\d+:(g|c)\."
            r"\(\?\_(\+|-|\*)?\d+(\+|-)?\d*\)"
            r"\_\((\+|-|\*)?\d+(\+|-)?\d*_"
            r"(\+|-|\*)?\d+(\+|-)?\d*\)"
            r"(del|dup|inv)$",
            quibble
    ):
        parts = quibble.split(')_(', 1)

        num1 = parts[0].rsplit('_', 1)[1]
        second_range = parts[1].partition(')')[0]
        num2, num3 = second_range.split('_', 1)

        intronic_positions = []

        num1 = parse_position(num1, intronic_positions)
        num2 = parse_position(num2, intronic_positions)
        num3 = parse_position(num3, intronic_positions)

        check_intronic_boundaries(intronic_positions)

        num1 = utr_to_int(num1)
        num2 = utr_to_int(num2)
        num3 = utr_to_int(num3)

        if num1 <= num2 <= num3:
            my_variant.warnings.append(
                "Uncertain positions are not fully supported, however the "
                "syntax is valid"
            )
        else:
            my_variant.warnings.append(
                "Uncertain positions are not fully supported, however the "
                "provided positions are out of order"
            )

        raise FuzzyRangeError(
            "Fuzzy/unknown variant end position in submitted variant "
            "description"
        )

    # If parsing has already produced an HGVS object, inspect its structured
    # positions rather than serialising the object and searching for '?'.
    else:
        logger.info(
            "%s has a fuzzy position but is not a range",
            my_variant.quibble
        )

        try:
            position = my_variant.hgvs_formatted.posedit.pos
            start = position.start
            end = position.end

            start_fuzzy = getattr(start, 'base', None) is None
            end_fuzzy = getattr(end, 'base', None) is None

            if end_fuzzy and not start_fuzzy:
                logger.info(
                    "Fuzzy/unknown variant end position in submitted "
                    "variant description"
                )
                raise FuzzyPositionError(
                    "Fuzzy/unknown variant end position in submitted "
                    "variant description"
                )

            if start_fuzzy and not end_fuzzy:
                logger.info(
                    "Fuzzy/unknown variant start position in submitted "
                    "variant description"
                )
                raise FuzzyPositionError(
                    "Fuzzy/unknown variant start position in submitted "
                    "variant description"
                )

            if start_fuzzy and end_fuzzy:
                logger.info(
                    "Fuzzy/unknown variant start and end positions in "
                    "submitted variant description"
                )
                raise FuzzyPositionError(
                    "Fuzzy/unknown variant start and end positions in "
                    "submitted variant description"
                )

        except AttributeError as e:
            logger.info("%s: %s", my_variant.quibble, e)

            # We are back at the unsupported-input boundary here, so direct
            # string inspection is appropriate.
            if isinstance(my_variant.quibble, str) and '?' in my_variant.quibble:
                logger.info(
                    "Fuzzy/unknown variant end position in submitted "
                    "variant description"
                )
                raise FuzzyPositionError(
                    "Fuzzy/unknown variant end position in submitted "
                    "variant description"
                )

    return False

def uncertain_positions(my_variant, validator):
    """
    Handle uncertain HGVS positions that cannot be parsed directly.

    The submitted description remains a string only while the unsupported
    uncertain-position syntax is interpreted. Once valid HGVS objects have
    been constructed, quibble is kept as an HGVS object.
    """
    evm = validator.no_norm_evm

    hgvs_accession, _sep, variation = my_variant.quibble.partition(':')

    if not variation.startswith(('g.(', 'c.(', 'n.(', 'r.(')):
        return

    logger.info(
        "Found potential uncertain positions in %s",
        my_variant.quibble
    )

    # Resolve LRG identifiers locally. Do not write an intermediate string
    # back to quibble.
    if hgvs_accession.startswith('LRG_'):
        if 't' in hgvs_accession:
            hgvs_accession = (
                validator.db.get_refseq_transcript_id_from_lrg_transcript_id(
                    hgvs_accession
                )
            )
        else:
            hgvs_accession = validator.db.get_refseq_id_from_lrg_id(
                hgvs_accession
            )

    # ------------------------------------------------------------------
    # Two uncertain ranges, e.g.
    # g.(90136803_90144453)_(90159675_90261231)dup
    # c.(375+1_376-1)_(672+1_673-1)del
    # ------------------------------------------------------------------
    if ')_(' in variation and '?' not in variation:
        accession = hgvs_accession
        var_type, _sep, positions_edit = variation.partition('.(')
        position_1, _sep, position_2_edit = positions_edit.partition(')_(')
        position_2, _sep, edit_string = position_2_edit.partition(')')

        my_variant.reftype = f':{var_type}.'

        try:
            if var_type == 'p':
                edit = validator.hp.parse_pro_edit(edit_string)
            elif var_type in ('c', 't', 'r'):
                edit = validator.hp.parse_rna_edit(edit_string)
            else:
                edit = validator.hp.parse_dna_edit(edit_string)
        except vvhgvs.exceptions.HGVSError as e:
            raise HgvsParseError(str(e))

        start_1, _sep, end_1 = position_1.partition('_')
        start_2, _sep, end_2 = position_2.partition('_')

        if not end_1:
            end_1 = start_1

        if not end_2:
            end_2 = start_2

        try:
            parsed_v1 = hgvs_obj_from_existing_edit(
                accession,
                var_type,
                start_1,
                copy.copy(edit),
                end=end_1
            )

            parsed_v2 = hgvs_obj_from_existing_edit(
                accession,
                var_type,
                start_2,
                copy.copy(edit),
                end=end_2
            )

        except vvhgvs.exceptions.HGVSError as e:
            raise HgvsParseError(str(e))

        for parsed_variant, position in (
                (parsed_v1, position_1),
                (parsed_v2, position_2)
        ):
            try:
                validator.vr.validate(parsed_variant)

            except vvhgvs.exceptions.HGVSError as e:
                error = str(e)

                if 'is not known to be compatible with variant type' in error:
                    raise IncompatibleTypeError(error)

                if 'base start position must be <= end position' in error:
                    raise InvalidRangeError(
                        f'{error} in position {parsed_variant.posedit.pos}'
                    )

                if not hgvs_position_utils.either_position_is_intronic(parsed_variant):
                    raise InvalidRangeError(
                        f'{position} is an invalid range for '
                        f'accession {accession}'
                    )

        if parsed_v1.posedit.pos.end.base >= parsed_v2.posedit.pos.start.base:
            raise InvalidRangeError(
                f'Position {parsed_v1.posedit.pos} is > or overlaps '
                f'{parsed_v2.posedit.pos}'
            )

        my_variant.warnings.append(
            'Uncertain positions are not fully supported, however the '
            'syntax is valid'
        )
        my_variant.reformat_output = 'uncertain_pos'

        # Authoritative representation of the complete uncertain variant.
        uncertain_variant = hgvs_obj_from_existing_edit(
            accession,
            var_type,
            FEInterval(
                start=parsed_v1.posedit.pos,
                end=parsed_v2.posedit.pos
            ),
            copy.copy(edit)
        )
        uncertain_variant.posedit.pos.uncertain = True

        my_variant.quibble = uncertain_variant

        # --------------------------------------------------------------
        # Chromosomal genomic route
        # --------------------------------------------------------------
        if accession.startswith('NC_'):
            my_variant.hgvs_genomic = uncertain_variant

            if (
                    validator.select_transcripts != 'select'
                    and validator.select_transcripts.startswith(
                        ('NM_', 'ENST')
                    )
                    and '|' not in validator.select_transcripts
                    and '[' not in validator.select_transcripts
            ):
                pass

            elif validator.select_transcripts != 'select':
                validator.select_transcripts = 'select'
                my_variant.warnings.append(
                    'Only a single transcript can be processed, updating '
                    'to select. Where no select transcript is identified '
                    'a suitable transcript will be used'
                )

            # relevant_transcripts() expects an ordinary HGVS interval whose
            # start and end expose .base directly. The authoritative
            # uncertain_variant contains an FEInterval, so construct a
            # separate operational HGVS object spanning the complete region.
            #
            # This object is used only for transcript discovery and never
            # replaces the uncertain HGVS representation.
            mapping_variant = hgvs_obj_from_existing_edit(
                accession,
                var_type,
                parsed_v1.posedit.pos.start,
                copy.copy(edit),
                end=parsed_v2.posedit.pos.end
            )
            mapping_variant.posedit.pos.uncertain = True

            rel_var = validator.relevant_transcripts(
                mapping_variant,
                evm,
                validator.alt_aln_method,
                my_variant.reverse_normalizer,
                validator.select_transcripts
            )

            if rel_var:
                my_variant.output_type_flag = 'gene'

                if validator.select_transcripts == 'select':
                    for transcript_variant in rel_var:
                        annotation = (
                            validator.db.get_transcript_annotation(
                                transcript_variant.ac
                            )
                        )

                        if '"select": "MANE"' in annotation:
                            validator.select_transcripts = (
                                transcript_variant.ac
                            )
                            break

                        if (
                                '"select": "RefSeq"' in annotation
                                or '"select": "Ensembl"' in annotation
                        ):
                            validator.select_transcripts = (
                                transcript_variant.ac
                            )
                            continue

                        if (
                                validator.alt_aln_method == 'splign'
                                and transcript_variant.ac.startswith('NM_')
                        ):
                            validator.select_transcripts = (
                                transcript_variant.ac
                            )

                        elif (
                                validator.alt_aln_method == 'genebuild'
                                and transcript_variant.ac.startswith('ENST')
                        ):
                            validator.select_transcripts = (
                                transcript_variant.ac
                            )

                ptv1 = validator.relevant_transcripts(
                    parsed_v1,
                    evm,
                    validator.alt_aln_method,
                    my_variant.reverse_normalizer,
                    validator.select_transcripts
                )

                ptv2 = validator.relevant_transcripts(
                    parsed_v2,
                    evm,
                    validator.alt_aln_method,
                    my_variant.reverse_normalizer,
                    validator.select_transcripts
                )

                ptv1 = next(
                    transcript_variant
                    for transcript_variant in ptv1
                    if transcript_variant.ac == validator.select_transcripts
                )

                ptv2 = next(
                    transcript_variant
                    for transcript_variant in ptv2
                    if transcript_variant.ac == validator.select_transcripts
                )

                if (
                        ptv1.posedit.pos.start.base
                        < ptv2.posedit.pos.start.base
                ):
                    t_position_1 = ptv1.posedit.pos
                    t_position_2 = ptv2.posedit.pos
                else:
                    t_position_1 = ptv2.posedit.pos
                    t_position_2 = ptv1.posedit.pos

                tx_variant = hgvs_obj_from_existing_edit(
                    ptv1.ac,
                    ptv1.type,
                    FEInterval(
                        start=t_position_1,
                        end=t_position_2
                    ),
                    copy.copy(edit)
                )
                tx_variant.posedit.pos.uncertain = True

                my_variant.hgvs_coding = tx_variant
                my_variant.hgvs_transcript_variant = tx_variant
                my_variant.quibble = tx_variant

            else:
                my_variant.output_type_flag = 'intergenic'
                my_variant.warnings.append(
                    'Selected transcript does not span the entire range '
                    'of the genomic variation'
                )

        # --------------------------------------------------------------
        # Transcript route
        # --------------------------------------------------------------
        elif accession.startswith(('NM_', 'NR_', 'ENST')):
            try:
                pgv1 = evm.t_to_g(parsed_v1)
            except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                raise UncertainConversionError(
                    f'Validation of {parsed_v1} - {e}'
                )

            try:
                pgv2 = evm.t_to_g(parsed_v2)
            except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                raise UncertainConversionError(
                    f'Validation of {parsed_v2} - {e}'
                )

            if pgv1.posedit.pos.start.base < pgv2.posedit.pos.start.base:
                g_position_1 = pgv1.posedit.pos
                g_position_2 = pgv2.posedit.pos
            else:
                g_position_1 = pgv2.posedit.pos
                g_position_2 = pgv1.posedit.pos

            gen_variant = hgvs_obj_from_existing_edit(
                pgv1.ac,
                pgv1.type,
                FEInterval(
                    start=g_position_1,
                    end=g_position_2
                ),
                copy.copy(edit)
            )
            gen_variant.posedit.pos.uncertain = True

            my_variant.hgvs_genomic = gen_variant
            my_variant.hgvs_coding = uncertain_variant
            my_variant.hgvs_transcript_variant = uncertain_variant
            my_variant.quibble = uncertain_variant
            my_variant.output_type_flag = 'gene'

        return

    # ------------------------------------------------------------------
    # Single uncertain range, e.g. g.(100_200)del
    # ------------------------------------------------------------------
    if '?' not in variation:
        if ')(' in variation:
            raise InvalidRangeError(
                'Invalid range submitted, missing underscore between '
                'stated uncertain positions'
            )

        logger.info(
            'Uncertain Position route 2 for %s:%s',
            hgvs_accession,
            variation
        )

        accession = hgvs_accession
        var_type, _sep, position_edit = variation.partition('.(')
        position, _sep, edit_string = position_edit.partition(')')

        my_variant.reftype = f':{var_type}.'

        try:
            if var_type == 'p':
                edit = validator.hp.parse_pro_edit(edit_string)
                eq_edit = validator.hp.parse_pro_edit('=')

            elif var_type in ('c', 't', 'r'):
                edit = validator.hp.parse_rna_edit(edit_string)
                eq_edit = validator.hp.parse_rna_edit('=')

            else:
                edit = validator.hp.parse_dna_edit(edit_string)
                eq_edit = validator.hp.parse_dna_edit('=')

            start, _sep, end = position.partition('_')

            if not end:
                end = start

            parsed_v1 = hgvs_obj_from_existing_edit(
                accession,
                var_type,
                start,
                eq_edit,
                end=end
            )

        except vvhgvs.exceptions.HGVSError as e:
            raise HgvsParseError(str(e))

        try:
            validator.vr.validate(parsed_v1)

        except vvhgvs.exceptions.HGVSError as e:
            error = str(e)

            if 'is not known to be compatible with variant type' in error:
                raise IncompatibleTypeError(error)

            if 'base start position must be <= end position' in error:
                raise InvalidRangeError(
                    f'{error} in position {parsed_v1.posedit.pos}'
                )

            if not hgvs_position_utils.either_position_is_intronic(parsed_v1):
                raise InvalidRangeError(
                    f'{position} is an invalid range for '
                    f'accession {accession}'
                )

        my_variant.warnings.append(
            'Uncertain positions are not fully supported, however the '
            'syntax is valid'
        )
        my_variant.reformat_output = 'uncertain_pos'

        complete_variant = hgvs_obj_from_existing_edit(
            accession,
            var_type,
            parsed_v1.posedit.pos,
            copy.copy(edit)
        )
        complete_variant.posedit.pos.uncertain = True

        my_variant.quibble = complete_variant

        # --------------------------------------------------------------
        # Chromosomal genomic route
        # --------------------------------------------------------------
        if accession.startswith('NC_'):
            my_variant.hgvs_genomic = complete_variant

            if (
                    validator.select_transcripts != 'select'
                    and validator.select_transcripts.startswith(
                        ('NM_', 'ENST')
                    )
                    and '|' not in validator.select_transcripts
                    and '[' not in validator.select_transcripts
            ):
                pass

            elif validator.select_transcripts != 'select':
                validator.select_transcripts = 'select'
                my_variant.warnings.append(
                    'Only a single transcript can be processed, '
                    'updating to Select'
                )

            rel_var = validator.relevant_transcripts(
                parsed_v1,
                evm,
                validator.alt_aln_method,
                my_variant.reverse_normalizer,
                validator.select_transcripts
            )

            if rel_var:
                my_variant.output_type_flag = 'gene'

                if validator.select_transcripts == 'select':
                    for transcript_variant in rel_var:
                        annotation = (
                            validator.db.get_transcript_annotation(
                                transcript_variant.ac
                            )
                        )

                        if '"select": "MANE"' in annotation:
                            validator.select_transcripts = (
                                transcript_variant.ac
                            )
                            break

                        if (
                                '"select": "RefSeq"' in annotation
                                or '"select": "Ensembl"' in annotation
                        ):
                            validator.select_transcripts = (
                                transcript_variant.ac
                            )

                ptv1 = next(
                    transcript_variant
                    for transcript_variant in rel_var
                    if transcript_variant.ac == validator.select_transcripts
                )

                tx_variant = hgvs_obj_from_existing_edit(
                    ptv1.ac,
                    ptv1.type,
                    ptv1.posedit.pos,
                    copy.copy(edit)
                )
                tx_variant.posedit.pos.uncertain = True

                my_variant.hgvs_coding = tx_variant
                my_variant.hgvs_transcript_variant = tx_variant
                my_variant.quibble = tx_variant

            else:
                my_variant.output_type_flag = 'intergenic'
                my_variant.warnings.append(
                    'Selected transcript does not span the entire range '
                    'of the genomic variation'
                )

        # --------------------------------------------------------------
        # Transcript route
        # --------------------------------------------------------------
        elif accession.startswith(('NM_', 'NR_', 'ENST')):
            try:
                pgv1 = evm.t_to_g(parsed_v1)
            except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                raise UncertainConversionError(
                    f'Validation of {parsed_v1} - {e}'
                )

            gen_variant = hgvs_obj_from_existing_edit(
                pgv1.ac,
                pgv1.type,
                pgv1.posedit.pos,
                copy.copy(edit)
            )
            gen_variant.posedit.pos.uncertain = True

            my_variant.hgvs_genomic = gen_variant
            my_variant.hgvs_coding = complete_variant
            my_variant.hgvs_transcript_variant = complete_variant
            my_variant.quibble = complete_variant
            my_variant.output_type_flag = 'gene'

        return

    # ------------------------------------------------------------------
    # Fuzzy ends
    # ------------------------------------------------------------------
    logger.info(
        'Looking for fuzzy ends in %s',
        my_variant.quibble
    )
    fuzzy_ends(my_variant, validator)

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
