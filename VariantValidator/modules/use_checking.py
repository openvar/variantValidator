import re
import vvhgvs
import vvhgvs.exceptions
import vvhgvs.variantmapper
import logging
from . import utils as fn
import copy

logger = logging.getLogger(__name__)


def refseq_common_mistakes(variant):
    """
    Evolving list of common mistakes, see sections below
    """
    # NM_ .g
    if (variant.quibble.startswith('NM_') or variant.quibble.startswith('NR_')) and variant.reftype == ':g.':
        suggestion = variant.quibble.replace(':g.', ':c.')
        error = 'Transcript reference sequence input as genomic (g.) reference sequence. ' \
                'Did you mean ' + suggestion + '?'
        variant.warnings.append(error)
        logger.warning(error)
        return True
    # NR_ c.
    if variant.quibble.startswith('NR_') and variant.reftype == ':c.':
        suggestion = variant.quibble.replace(':c.', ':n.')
        error = 'Non-coding transcript reference sequence input as coding (c.) reference sequence. ' \
                'Did you mean ' + suggestion + '?'
        variant.warnings.append(error)
        logger.warning(error)
        return True
    # NM_ n.
    if variant.quibble.startswith('NM_') and variant.reftype == ':n.':
        suggestion = variant.quibble.replace(':n.', ':c.')
        error = 'Coding transcript reference sequence input as non-coding transcript (n.) reference sequence. ' \
                'Did you mean ' + suggestion + '?'
        variant.warnings.append(error)
        logger.warning(error)
        return True

    # NM_ NC_ NG_ NR_ p.
    if (variant.quibble.startswith('NM_') or variant.quibble.startswith('NR_') or variant.quibble.startswith('NC_') or
            variant.quibble.startswith('NG_')) and variant.reftype == ':p.':
        error = 'Using a nucleotide reference sequence (NM_ NR_ NG_ NC_) to specify protein-level (p.) variation is ' \
                'not HGVS compliant. Please select an appropriate protein reference sequence (NP_)'
        variant.warnings.append(error)
        logger.warning(error)
        return True

    # NG_ c or NC_c..
    if (variant.quibble.startswith('NG_') or variant.quibble.startswith('NC_')) and variant.reftype == ':c.':
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
    input_parses = validator.hp.parse_hgvs_variant(variant.quibble)
    variant.input_parses = input_parses
    variant.gene_symbol = validator.db.get_gene_symbol_from_transcript_id(variant.input_parses.ac)
    if variant.gene_symbol == 'none':
        variant.gene_symbol = ''
    if input_parses.type == 'g':
        check = structure_checks_g(variant, validator)
        if check:
            return True

    elif input_parses.type == 'c':
        check = structure_checks_c(variant, validator)
        if check:
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
    if not variant.quibble.startswith('NC_') and not variant.quibble.startswith('NG_') \
            and not variant.quibble.startswith('NT_') and not variant.quibble.startswith('NW_'):
        error = 'Invalid reference sequence identifier (' + variant.input_parses.ac + ')'
        variant.warnings.append(error)
        logger.warning(error)
        return True

    try:
        validator.vr.validate(variant.input_parses)
    except Exception as e:
        if "does not agree with reference sequence ()" in str(e):
            e = "The specified coordinate is outside the boundaries of reference sequence %s " % variant.input_parses.ac
        error = str(e)
        variant.warnings.append(error)
        logger.warning(error)
        return True

    # Additional test
    try:
        variant.hn.normalize(variant.input_parses)
    except vvhgvs.exceptions.HGVSError as e:
        error = str(e)
        variant.warnings.append(error)
        logger.warning(error)
        return True

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
                                                                               prevalidation_level=None
                                                                               )

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
                variant.warnings.append(error)
                logger.warning(error)
                return True

        # Create a specific minimal evm with no normalizer and no replace_reference
        # Have to use this method due to potential multi chromosome error, note normalizes but does not replace sequence
        output = None
        try:
            output = validator.noreplace_myevm_t_to_g(variant.input_parses, variant)
        except vvhgvs.exceptions.HGVSDataNotAvailableError:
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
            error = str(e)
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
            # This catches errors in introns
            if 'base start position must be <= end position' in error:
                # correction = variant.input_parses
                # st = variant.input_parses.posedit.pos.start
                # ed = variant.input_parses.posedit.pos.end
                # correction.posedit.pos.start = ed
                # correction.posedit.pos.end = st
                # error = error + ': Did you mean ' + str(correction) + '?'
                error = 'Interval start position ' + str(variant.input_parses.posedit.pos.start) + ' > interval end '\
                        'position ' + str(variant.input_parses.posedit.pos.end)
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
                    back_to_n = variant.evm.g_to_t(test_g, variant.input_parses.ac)
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
