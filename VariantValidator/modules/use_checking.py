import re
import hgvs
from . import vvFunctions as fn
from .vvLogging import logger
import copy


def refseq_common_mistakes(variant):
    """
    Evolving list of common mistakes, see sections below
    """
    # NM_ .g
    if (variant.quibble.startswith('NM_') or variant.quibble.startswith('NR_')) and variant.reftype == ':g.':
        suggestion = variant.quibble.replace(':g.', ':c.')
        error = 'Transcript reference sequence input as genomic (g.) reference sequence. ' \
                'Did you mean ' + suggestion + '?'
        variant.warnings += ': ' + error
        logger.warning(error)
        return True
    # NR_ c.
    if variant.quibble.startswith('NR_') and variant.reftype == ':c.':
        suggestion = variant.quibble.replace(':c.', ':n.')
        error = 'Non-coding transcript reference sequence input as coding (c.) reference sequence. ' \
                'Did you mean ' + suggestion + '?'
        variant.warnings += ': ' + error
        logger.warning(error)
        return True
    # NM_ n.
    if variant.quibble.startswith('NM_') and variant.reftype == ':n.':
        suggestion = variant.quibble.replace(':n.', ':c.')
        error = 'Coding transcript reference sequence input as non-coding transcript (n.) reference sequence. ' \
                'Did you mean ' + suggestion + '?'
        variant.warnings += ': ' + error
        logger.warning(error)
        return True

    # NM_ NC_ NG_ NR_ p.
    if (variant.quibble.startswith('NM_') or variant.quibble.startswith('NR_') or variant.quibble.startswith('NC_') or
        variant.quibble.startswith('NG_')) and variant.reftype == ':p.':
        issue_link = 'http://varnomen.hgvs.org/recommendations/protein/'
        error = 'Using a nucleotide reference sequence (NM_ NR_ NG_ NC_) to specify protein-level (p.) variation is ' \
                'not HGVS compliant. Please select an appropriate protein reference sequence (NP_)'
        variant.warnings += ': ' + error
        logger.warning(error)
        return True

    # NG_ c or NC_c..
    if (variant.quibble.startswith('NG_') or variant.quibble.startswith('NC_')) and variant.reftype == ':c.':
        suggestion = ': For additional assistance, submit ' + str(variant.quibble) + ' to VariantValidator'
        error = 'NG_:c.PositionVariation descriptions should not be used unless a transcript reference sequence has ' \
                'also been provided e.g. NG_(NM_):c.PositionVariation' + suggestion
        variant.warnings += ': ' + error
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
        error = 'Invalid reference sequence identifier (' + variant.hgvs_formatted.ac + ')'
        variant.warnings += ': ' + str(error)
        logger.warning(error)
        return True

    try:
        validator.vr.validate(variant.hgvs_formatted)
    except Exception as e:
        error = str(e)
        variant.warnings += ': ' + str(error)
        logger.warning(error)
        return True

    # Additional test
    try:
        variant.hn.normalize(variant.hgvs_formatted)
    except hgvs.exceptions.HGVSError as e:
        error = str(e)
        variant.warnings += ': ' + str(error)
        logger.warning(error)
        return True

    return False


def structure_checks_c(variant, validator):
    """
    structure checks for when reftype is coding
    :param variant:
    :param validator:
    :param hn:
    :return:
    """

    if '*' in str(variant.hgvs_formatted) or 'c.-' in str(variant.hgvs_formatted):
        # Catch variation in UTRs
        # These should be in the sequence so can be directly validated. Need to pass to n.
        try:
            validator.vr.validate(variant.hgvs_formatted)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if 'datums is ill-defined' in error:
                called_ref = variant.hgvs_formatted.posedit.edit.ref
                try:
                    to_n = variant.evm.c_to_n(variant.hgvs_formatted)
                except hgvs.exceptions.HGVSInvalidVariantError as e:
                    error = str(e)
                    variant.warnings += ': ' + error
                    logger.warning(error)
                    return True
                actual_ref = to_n.posedit.edit.ref
                if called_ref != actual_ref:
                    error = 'Variant reference (' + called_ref + ') does not agree with reference sequence ' \
                                                                 '(' + actual_ref + ')'
                    variant.warnings += ': ' + error
                    logger.warning(error)
                    return True
                else:
                    variant.hgvs_formatted.posedit.edit.ref = ''
            else:
                if 'bounds' in error or 'intronic variant' in error:
                    try:
                        variant.hn.normalize(variant.hgvs_formatted)
                    except hgvs.exceptions.HGVSError:
                        fn.exceptPass()

                    if 'bounds' in error:
                        try:
                            identity_info = validator.hdp.get_tx_identity_info(variant.hgvs_formatted.ac)
                            ref_start = identity_info[3]
                            ref_end = identity_info[4]
                            if '-' in str(variant.hgvs_formatted.posedit.pos.start) and variant.hgvs_formatted.posedit.pos.start.offset == 0:
                                # upstream positions
                                boundary = -ref_start
                                remainder = variant.hgvs_formatted.posedit.pos.start.base - boundary
                                variant.hgvs_formatted.posedit.pos.start.base = boundary
                                variant.hgvs_formatted.posedit.pos.start.offset = remainder
                            if '-' in str(variant.hgvs_formatted.posedit.pos.end) and variant.hgvs_formatted.posedit.pos.end.offset == 0:
                                boundary = -ref_start
                                remainder = variant.hgvs_formatted.posedit.pos.end.base - boundary
                                variant.hgvs_formatted.posedit.pos.end.base = boundary
                                variant.hgvs_formatted.posedit.pos.end.offset = remainder
                            if '*' in str(variant.hgvs_formatted.posedit.pos.start) and variant.hgvs_formatted.posedit.pos.start.offset == 0:
                                # downstream positions
                                tot_end_pos = str(variant.hgvs_formatted.posedit.pos.start).replace('*', '')
                                ts_seq = validator.sf.fetch_seq(variant.hgvs_formatted.ac)
                                boundary = len(ts_seq) - ref_end
                                variant.hgvs_formatted.posedit.pos.start.base = boundary
                                offset = int(tot_end_pos) - boundary
                                variant.hgvs_formatted.posedit.pos.start.offset = offset
                            if '*' in str(variant.hgvs_formatted.posedit.pos.end) and variant.hgvs_formatted.posedit.pos.end.offset == 0:
                                tot_end_pos = str(variant.hgvs_formatted.posedit.pos.end).replace('*', '')
                                ts_seq = validator.sf.fetch_seq(variant.hgvs_formatted.ac)
                                boundary = len(ts_seq) - ref_end
                                variant.hgvs_formatted.posedit.pos.end.base = boundary
                                offset = int(tot_end_pos) - boundary
                                variant.hgvs_formatted.posedit.pos.end.offset = offset

                            # Create a lose vm instance
                            variant.lose_vm = hgvs.variantmapper.VariantMapper(validator.hdp,
                                                                               replace_reference=True,
                                                                               prevalidation_level=None
                                                                               )

                            report_gen = validator.myevm_t_to_g(variant.hgvs_formatted, variant.no_norm_evm,
                                                                variant.primary_assembly, variant.hn)
                            error = 'Using a transcript reference sequence to specify a variant position that lies ' \
                                    'outside of the reference sequence is not HGVS-compliant: ' \
                                    'Instead use ' + fn.valstr(report_gen)
                        except Exception:
                            fn.exceptPass()
                        variant.warnings += ': ' + error
                        logger.warning(error)
                        return True

        try:
            variant.hgvs_formatted = variant.evm.c_to_n(variant.hgvs_formatted)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            variant.warnings += ': ' + error
            logger.warning(e)
            return True

        if 'n.1-' in str(variant.hgvs_formatted):
            input_parses = variant.evm.n_to_c(variant.hgvs_formatted)
            error = 'Using a transcript reference sequence to specify a variant position that lies outside of the ' \
                    'reference sequence is not HGVS-compliant. Instead use '
            genomic_position = validator.myevm_t_to_g(input_parses, variant.no_norm_evm, variant.primary_assembly,
                                                      variant.hn)
            error = error + fn.valstr(genomic_position)
            variant.warnings += ': ' + error
            logger.warning(error)
            return True

        # Re-map input_parses back to c. variant
        variant.hgvs_formatted = variant.evm.n_to_c(variant.hgvs_formatted)

        # Intronic positions in UTRs
        if re.search(r'\d\-\d', str(variant.hgvs_formatted)) or re.search(r'\d\+\d', str(variant.hgvs_formatted)):
            # Can we go c-g-c
            try:
                to_genome = validator.myevm_t_to_g(variant.hgvs_formatted, variant.no_norm_evm,
                                                   variant.primary_assembly, variant.hn)
                to_tx = variant.evm.g_to_t(to_genome, variant.hgvs_formatted.ac)
            except hgvs.exceptions.HGVSInvalidIntervalError as e:
                error = str(e)
                if 'bounds' in error:
                    try:
                        identity_info = validator.hdp.get_tx_identity_info(variant.hgvs_formatted.ac)
                        ref_start = identity_info[3]
                        ref_end = identity_info[4]
                        if '-' in str(variant.hgvs_formatted.posedit.pos.start):
                            # upstream positions
                            boundary = -ref_start
                            remainder = variant.hgvs_formatted.posedit.pos.start.base - boundary
                            variant.hgvs_formatted.posedit.pos.start.base = boundary
                            variant.hgvs_formatted.posedit.pos.start.offset = remainder
                        if '-' in str(variant.hgvs_formatted.posedit.pos.end):
                            boundary = -ref_start
                            remainder = variant.hgvs_formatted.posedit.pos.end.base - boundary
                            variant.hgvs_formatted.posedit.pos.end.base = boundary
                            variant.hgvs_formatted.posedit.pos.end.offset = remainder
                        if '*' in str(variant.hgvs_formatted.posedit.pos.start):
                            # downstream positions
                            tot_end_pos = str(variant.hgvs_formatted.posedit.pos.start).replace('*', '')
                            ts_seq = validator.sf.fetch_seq(variant.hgvs_formatted.ac)
                            boundary = len(ts_seq) - ref_end
                            variant.hgvs_formatted.posedit.pos.start.base = boundary
                            te1, te2 = tot_end_pos.split('+')
                            tot_end_pos = int(te1) + int(te2)
                            offset = tot_end_pos - boundary
                            variant.hgvs_formatted.posedit.pos.start.offset = offset
                        if '*' in str(variant.hgvs_formatted.posedit.pos.end):
                            tot_end_pos = str(variant.hgvs_formatted.posedit.pos.end).replace('*', '')
                            ts_seq = validator.sf.fetch_seq(variant.hgvs_formatted.ac)
                            boundary = len(ts_seq) - ref_end
                            variant.hgvs_formatted.posedit.pos.end.base = boundary
                            te1, te2 = tot_end_pos.split('+')
                            tot_end_pos = int(te1) + int(te2)
                            offset = tot_end_pos - boundary
                            variant.hgvs_formatted.posedit.pos.end.offset = offset

                        report_gen = validator.myevm_t_to_g(variant.hgvs_formatted, variant.no_norm_evm,
                                                       variant.primary_assembly, variant.hn)
                        error = 'Using a transcript reference sequence to specify a variant position that lies ' \
                                'outside of the reference sequence is not HGVS-compliant. Instead use '\
                                + fn.valstr(report_gen)
                    except Exception:
                        fn.exceptPass()
                variant.warnings += ': ' + error
                logger.warning(error)
                return True

            except hgvs.exceptions.HGVSDataNotAvailableError as e:
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
                    error = 'Cannot map ' + fn.valstr(variant.hgvs_formatted) + ' to a genomic position. '\
                            + variant.hgvs_formatted.ac + ' can only be partially aligned to genomic reference ' \
                                                          'sequences ' + acs
                    variant.warnings += ': ' + error
                    logger.warning(error)
                    return True

    elif re.search(r'\d-', str(variant.hgvs_formatted)) or re.search(r'\d\+', str(variant.hgvs_formatted)):
        # Quick look at syntax validation
        try:
            validator.vr.validate(variant.hgvs_formatted)
        except hgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)
            if 'bounds' in error:
                try:
                    report_gen = validator.myevm_t_to_g(variant.hgvs_formatted, variant.no_norm_evm,
                                                        variant.primary_assembly, variant.hn)
                except hgvs.exceptions.HGVSError:
                    fn.exceptPass()
                else:
                    error = 'Using a transcript reference sequence to specify a variant position that lies outside of '\
                            'the reference sequence is not HGVS-compliant. Instead use ' + fn.valstr(report_gen)
                variant.warnings += ': ' + error
                logger.warning(error)
                return True
            elif 'insertion length must be 1' in error:
                variant.warnings += ': ' + error
                logger.warning(error)
                return True
            elif 'base start position must be <= end position' in error:
                correction = copy.deepcopy(variant.hgvs_formatted)
                st = variant.hgvs_formatted.posedit.pos.start
                ed = variant.hgvs_formatted.posedit.pos.end
                correction.posedit.pos.start = ed
                correction.posedit.pos.end = st
                error = error + ': Did you mean ' + str(correction) + '?'
                variant.warnings += ': ' + error
                logger.warning(error)
                return True

        # Create a specific minimal evm with no normalizer and no replace_reference
        # Have to use this method due to potential multi chromosome error, note normalizes but does not replace sequence
        try:
            output = validator.noreplace_myevm_t_to_g(variant.hgvs_formatted, variant.evm, validator.hdp,
                                                      variant.primary_assembly, validator.vm, variant.hn, validator.hp,
                                                      validator.sf, variant.no_norm_evm)
        except hgvs.exceptions.HGVSDataNotAvailableError:
            tx_ac = variant.hgvs_formatted.ac
            try:
                gene_symbol = validator.db.get_gene_symbol_from_transcriptID(tx_ac)
            except:
                error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, ' \
                        'please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' to  ' \
                        'https://variantvalidator.org/ref_finder/, or select an alternative genome build'
            else:
                error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, ' \
                        'please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' or ' \
                        + gene_symbol + ' to  https://variantvalidator.org/ref_finder/, or select an alternative ' \
                                        'genome build'
            variant.warnings += ': ' + error
            logger.warning(error)
            return True
        except ValueError as e:
            error = str(e)
            if '> end' in error:
                error = 'Interval start position ' + str(variant.hgvs_formatted.posedit.pos.start) + ' > interval end '\
                        'position ' + str(variant.hgvs_formatted.posedit.pos.end)
                variant.warnings += ': ' + error
                logger.warning(error)
                return True
        except hgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)
            if 'base start position must be <= end position' in error:
                # correction = copy.deepcopy(variant.hgvs_formatted)
                # st = variant.hgvs_formatted.posedit.pos.start
                # ed = variant.hgvs_formatted.posedit.pos.end
                # correction.posedit.pos.start = ed
                # correction.posedit.pos.end = st
                # error = error + ': Did you mean ' + str(correction) + '?'
                error = 'Interval start position ' + str(variant.hgvs_formatted.posedit.pos.start) + ' > interval end' \
                        ' position ' + str(variant.hgvs_formatted.posedit.pos.end)
                variant.warnings += ': ' + error
                logger.warning(error)
                return True
            else:
                variant.warnings += ': ' + error
                logger.warning(error)
                return True

        try:
            variant.evm.g_to_t(output, variant.hgvs_formatted.ac)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            variant.warnings += ': ' + error
            logger.warning(error)
            return True

        try:
            validator.vr.validate(output)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            variant.warnings += ': ' + error
            logger.warning(error)
            return True

    else:
        # All other variation
        try:
            validator.vr.validate(variant.hgvs_formatted)
        except hgvs.exceptions.HGVSUnsupportedOperationError:
            fn.exceptPass()
        except hgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)
            # This catches errors in introns
            if 'base start position must be <= end position' in error:
                # correction = variant.hgvs_formatted
                # st = variant.hgvs_formatted.posedit.pos.start
                # ed = variant.hgvs_formatted.posedit.pos.end
                # correction.posedit.pos.start = ed
                # correction.posedit.pos.end = st
                # error = error + ': Did you mean ' + str(correction) + '?'
                error = 'Interval start position ' + str(variant.hgvs_formatted.posedit.pos.start) + ' > interval end '\
                        'position ' + str(variant.hgvs_formatted.posedit.pos.end)
            variant.warnings += ': ' + error
            logger.warning(error)
            return True

        except hgvs.exceptions.HGVSDataNotAvailableError as e:
            error = str(e)
            variant.warnings += ': ' + error
            logger.warning(error)
            return True
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if 'bounds' in error:
                error += ' (' + variant.hgvs_formatted.ac + ')'
                variant.warnings += ': ' + error
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
    if '+' in str(variant.hgvs_formatted) or '-' in str(variant.hgvs_formatted):
        # Catch variation in UTRs
        # These should be in the sequence so can be directly validated. Need to pass to n.
        try:
            validator.vr.validate(variant.hgvs_formatted)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if 'intronic variant' in error:
                pass
            elif 'datums is ill-defined' in error:
                called_ref = variant.hgvs_formatted.posedit.edit.ref
                to_n = variant.evm.c_to_n(variant.hgvs_formatted)
                actual_ref = to_n.posedit.edit.ref
                if called_ref != actual_ref:
                    error = 'Variant reference (' + called_ref + ') does not agree with reference sequence (' + actual_ref + ')'
                    variant.warnings += ': ' + str(error)
                    logger.warning(str(error))
                    return True
                else:
                    variant.hgvs_formatted.posedit.edit.ref = ''
                    formatted_variant = str(variant.hgvs_formatted)

            elif 'base must be >=1 for datum = SEQ_START or CDS_END' in error:
                error = 'The given coordinate is outside the bounds of the reference sequence.'

                try:
                    if '-' in str(variant.hgvs_formatted.posedit.pos.start):
                        # upstream positions
                        boundary = 1
                        remainder = variant.hgvs_formatted.posedit.pos.start.base - boundary
                        remainder = remainder + 1
                        variant.hgvs_formatted.posedit.pos.start.base = boundary
                        variant.hgvs_formatted.posedit.pos.start.offset = remainder
                    if '-' in str(variant.hgvs_formatted.posedit.pos.end):
                        boundary = 1
                        remainder = variant.hgvs_formatted.posedit.pos.end.base - boundary
                        remainder = remainder + 1
                        variant.hgvs_formatted.posedit.pos.end.base = boundary
                        variant.hgvs_formatted.posedit.pos.end.offset = remainder
                    report_gen = validator.myevm_t_to_g(variant.hgvs_formatted, variant.no_norm_evm, variant.primary_assembly,
                                                   variant.hn)
                    error = 'Using a transcript reference sequence to specify a variant position that lies outside of the reference sequence is not HGVS-compliant. Instead use ' + fn.valstr(
                        report_gen)
                except Exception:
                    fn.exceptPass()
                variant.warnings += ': ' + error
                logger.warning(error)
                return True
            else:
                variant.warnings += ': ' + error
                logger.warning(error)
                return True

    if 'n.1-' in str(variant.hgvs_formatted):
        error = 'Using a transcript reference sequence to specify a variant position that lies outside of the reference sequence is not HGVS-compliant. Instead use '
        genomic_position = validator.myevm_t_to_g(variant.hgvs_formatted, variant.no_norm_evm, variant.primary_assembly,
                                             variant.hn)
        error = error + fn.valstr(genomic_position)
        variant.warnings += ': ' + error
        logger.warning(error)
        return True

    if re.search(r'\d-', str(variant.hgvs_formatted)) or re.search(r'\d\+', str(variant.hgvs_formatted)):
        # Quick look at syntax validation
        try:
            validator.vr.validate(variant.hgvs_formatted)
        except hgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)
            if 'bounds' in error:
                try:
                    report_gen = validator.myevm_t_to_g(variant.hgvs_formatted, variant.no_norm_evm, variant.primary_assembly,
                                                   variant.hn)
                except hgvs.exceptions.HGVSError as e:
                    fn.exceptPass()
                else:
                    error = 'Using a transcript reference sequence to specify a variant position that lies outside of the reference sequence is not HGVS-compliant. Instead use ' + fn.valstr(
                        report_gen)
                variant.warnings += ': ' + error
                logger.warning(error)
                return True
            elif 'insertion length must be 1' in error:
                variant.warnings += ': ' + error
                logger.warning(error)
                return True
            elif 'base start position must be <= end position' in error:
                correction = copy.deepcopy(variant.hgvs_formatted)
                st = variant.hgvs_formatted.posedit.pos.start
                ed = variant.hgvs_formatted.posedit.pos.end
                correction.posedit.pos.start = ed
                correction.posedit.pos.end = st
                error = error + ': Did you mean ' + str(correction) + '?'
                # error = 'Interval start position ' + str(input_parses.posedit.pos.start) + ' > interval end position ' + str(input_parses.posedit.pos.end)
                variant.warnings += ': ' + error
                logger.warning(error)
                return True
            elif 'Cannot validate sequence of an intronic variant' in error:
                try:
                    test_g = validator.myevm_t_to_g(variant.hgvs_formatted, variant.no_norm_evm, variant.primary_assembly,
                                               variant.hn)
                    back_to_n = variant.evm.g_to_t(test_g, variant.hgvs_formatted.ac)
                except hgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'bounds' in error:
                        report_gen = validator.myevm_t_to_g(variant.hgvs_formatted, variant.no_norm_evm,
                                                       variant.primary_assembly, variant.hn)
                        error = 'Using a transcript reference sequence to specify a variant position that lies outside of the reference sequence is not HGVS-compliant. Instead use ' + fn.valstr(
                            report_gen)
                        variant.warnings += ': ' + error
                        logger.warning(error)
                        return True

        # Create a specific minimal evm with no normalizer and no replace_reference
        # Have to use this method due to potential multi chromosome error, note, normalizes but does not replace sequence
        try:
            output = validator.noreplace_myevm_t_to_g(variant.hgvs_formatted, variant.evm, validator.hdp, variant.primary_assembly, validator.vm, variant.hn,
                                                 validator.hp, validator.sf, variant.no_norm_evm)
        except hgvs.exceptions.HGVSDataNotAvailableError as e:
            tx_ac = variant.hgvs_formatted.ac
            try:
                gene_symbol = validator.db.get_gene_symbol_from_transcriptID(tx_ac)
            except:
                gene_symbol = None
            if gene_symbol is None:
                error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
            else:
                error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' or ' + gene_symbol + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
            variant.warnings += ': ' + error
            logger.warning(error)
            return True
        except ValueError as e:
            error = str(e)
            if '> end' in error:
                error = 'Interval start position ' + str(
                    variant.hgvs_formatted.posedit.pos.start) + ' > interval end position ' + str(
                    variant.hgvs_formatted.posedit.pos.end)
                variant.warnings += ': ' + error
                logger.warning(error)
                return True
        except hgvs.exceptions.HGVSInvalidVariantError as e:
            error = str(e)
            if 'base start position must be <= end position' in error:
                correction = copy.deepcopy(variant.hgvs_formatted)
                st = variant.hgvs_formatted.posedit.pos.start
                ed = variant.hgvs_formatted.posedit.pos.end
                correction.posedit.pos.start = ed
                correction.posedit.pos.end = st
                error = error + ': Did you mean ' + str(correction) + '?'
                error = 'Interval start position ' + str(
                    variant.hgvs_formatted.posedit.pos.start) + ' > interval end position ' + str(
                    variant.hgvs_formatted.posedit.pos.end)
                variant.warnings += ': ' + error
                logger.warning(error)
                return True
        try:
            validator.vr.validate(output)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            variant.warnings += ': ' + error
            logger.warning(error)
            return True

    else:
        # All other variation
        try:
            validator.vr.validate(variant.hgvs_formatted)
        except hgvs.exceptions.HGVSUnsupportedOperationError:
            fn.exceptPass()
        except hgvs.exceptions.HGVSInvalidVariantError as e:
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
                # correction = copy.deepcopy(variant.hgvs_formatted)
                # st = variant.hgvs_formatted.posedit.pos.start
                # ed = variant.hgvs_formatted.posedit.pos.end
                # correction.posedit.pos.start = ed
                # correction.posedit.pos.end = st
                # error = error + ': Did you mean ' + str(correction) + '?'
                error = 'Interval start position ' + str(
                    variant.hgvs_formatted.posedit.pos.start) + ' > interval end position ' + str(
                    variant.hgvs_formatted.posedit.pos.end)
                logger.warning(error)
                variant.warnings += ': ' + error
                return True
            variant.warnings += ': ' + error
            logger.warning(error)
            return True
        except hgvs.exceptions.HGVSDataNotAvailableError as e:
            error = str(e)
            variant.warnings += ': ' + error
            logger.warning(error)
            return True
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if 'bounds' in error:
                error = error + ' (' + variant.hgvs_formatted.ac + ')'
                variant.warnings += ': ' + error
                logger.warning(error)
                return True
    return False
