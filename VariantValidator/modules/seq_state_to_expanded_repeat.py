import logging
import re
from VariantValidator.modules import utils

# Custom exceptions for better error granularity
class RepeatedUnitError(Exception):
    """Raised when the repeated unit is empty or invalid."""
    pass

class StartPositionError(Exception):
    """Raised when the start position is missing or invalid."""
    pass

class VariantFormatError(Exception):
    """Raised when the variant type is not recognized or improperly formatted."""
    pass

class VariantDataError(Exception):
    """Raised when the variant data is incomplete or malformed."""
    pass

logger = logging.getLogger(__name__)


def reassemble_expanded_repeat_variant(reference, reference_start, reference_end, reference_type,
                                       repeated_unit, inserted_sequence, validator, remove_units=0,
                                       converted_from_intronic=False, alt_aln_method=None):
    """
    Assemble the final expanded repeat representation in HGVS-like format.

    Example: NG_012232.1:g.3_6T[20]

    Args explain how all parts are used to compute total repeat count.
    """
    if not repeated_unit:
        raise RepeatedUnitError("Repeated unit must be a non-empty string.")

    unit_len = len(repeated_unit)
    inserted_repeat_units = len(inserted_sequence) // unit_len
    reference_repeat_units = (reference_end - reference_start + 1) // unit_len
    total_units = reference_repeat_units + inserted_repeat_units - remove_units

    # Begin constructing the variant with positional metadata

    if converted_from_intronic is False:
        variant_identity_format = f"{reference}:{reference_type}{reference_start}_{reference_end}="
    else:
        variant_identity_format = f"{reference}:g.{reference_start}_{reference_end}="

    # If the variant is coding (c.), convert to nucleotide (n.) and back to ensure alignment
    if reference_type == "c." and "NC_" not in reference:
        reparse = validator.hp.parse(variant_identity_format.replace("c.", "n."))
        variant_identity_format = utils.valstr(validator.vm.n_to_c(reparse))

    elif converted_from_intronic is not False:
        variant_identity_format = utils.valstr(validator.vm.g_to_t(validator.hp.parse(variant_identity_format), converted_from_intronic,
                                                                   alt_aln_method=alt_aln_method))
        orientation = validator.hdp.get_tx_exons(converted_from_intronic, reference, alt_aln_method=alt_aln_method)[0][3]
        if orientation != 1:
            # If the transcript is on the reverse strand, we need reverse complement the repeat bases
            repeated_unit = validator.revcomp(repeated_unit)

    # Replace = with actual repeat structure like T[20]
    return variant_identity_format.replace("=", f"{repeated_unit}[{total_units}]")

def decipher_repeated_unit(sequence):
    """
    Detect the smallest unit that repeats to form the given DNA sequence.
    If no repeat unit found, return the original sequence.
    """
    n = len(sequence)
    for unit_len in range(1, n // 2 + 1):
        unit = sequence[:unit_len]
        multiplier = n // unit_len
        if unit * multiplier == sequence[:unit_len * multiplier]:
            return unit
    return sequence  # Default: no recognizable repeat pattern

def decipher_start_of_full_reference_repeated_sequence(reference, repeated_unit, start, validator, window_size=100):
    """
    Given a position where a repeat starts, fetch a larger window of sequence
    and walk backwards in-memory to find the full beginning of the repeat block.
    """
    if not repeated_unit:
        raise RepeatedUnitError("Repeated unit must be a non-empty string.")
    if start is None or not isinstance(start, int) or start <= 0:
        raise StartPositionError("Start position must be a positive integer.")

    unit_len = len(repeated_unit)
    l_start = start - 1  # convert to 0-based

    # Determine how far back we can fetch
    fetch_start = max(0, l_start - window_size)
    fetch_end = l_start + unit_len

    # Fetch once
    seq_window = validator.sf.fetch_seq(reference, start_i=fetch_start, end_i=fetch_end)

    # Position of the unit in the fetched window
    local_pos = l_start - fetch_start
    current_chunk = seq_window[local_pos:local_pos + unit_len]

    if current_chunk != repeated_unit:
        return None  # Repeat does not match expected unit

    # Walk backwards within the in-memory window
    while local_pos - unit_len >= 0:
        prev_chunk = seq_window[local_pos - unit_len:local_pos]
        if prev_chunk != repeated_unit:
            break
        local_pos -= unit_len

    # Convert back to reference coordinate (1-based)
    return fetch_start + local_pos + 1


def decipher_end_of_full_reference_repeated_sequence(reference, repeated_unit, start, validator, window_size=100):
    """
    Walk forward from a given repeat position to find the full extent of the repeated block.
    Fetch a larger window once to avoid multiple fetch_seq calls.
    """
    if not repeated_unit:
        raise RepeatedUnitError("Repeated unit must be a non-empty string.")
    if start is None or not isinstance(start, int) or start <= 0:
        raise StartPositionError("Start position must be a positive integer.")

    unit_len = len(repeated_unit)

    # Try two offsets for frame alignment
    for start_0 in [start - 1 - unit_len, start - 1]:
        if start_0 < 0:
            continue

        # Fetch a single forward window
        fetch_end = start_0 + window_size
        try:
            seq_window = validator.sf.fetch_seq(reference, start_i=start_0, end_i=fetch_end)
        except Exception:
            continue

        # Check if the repeated unit matches at this start
        local_pos = 0
        if seq_window[local_pos:local_pos + unit_len] != repeated_unit:
            continue

        # Walk forward in-memory
        while local_pos + unit_len <= len(seq_window) - unit_len:
            next_chunk = seq_window[local_pos + unit_len:local_pos + 2 * unit_len]
            if next_chunk != repeated_unit:
                break
            local_pos += unit_len

        return start_0 + local_pos + unit_len  # convert to 1-based end position

    return None

def convert_seq_state_to_expanded_repeat(variant, validator, genomic_reference=None, known_repeat_unit=None):
    """
    Main interface: convert a HGVS variant string to an expanded repeat representation.

    Supports: insertions (ins), deletions (del), duplications (dup), identity (=)
    """
    if any(key in variant for key in ['ins', 'del', 'dup', '=']):
        hgvs_variant = validator.hp.parse(variant)
        # Use ENST normalizer or splign depending on variant format
        hn = validator.genebuild_normalizer if "ENST" in variant else validator.splign_normalizer
        alt_aln_method = "genebuild" if "ENST" in variant else "splign"
        hgvs_variant = hn.normalize(hgvs_variant)

        # Determine reference type and extract sequence info
        converted_from_intronic = False
        if ":g." in variant:
            reference_type = "g."
            parts = str(hgvs_variant).split(':g.')
        elif ":c." in variant or ":n." in variant:
            if ":c." in variant:
                reference_type = "c."
                hgvs_n = validator.vm.c_to_n(hgvs_variant)
            elif ":n." in variant:
                hgvs_n = hgvs_variant
                reference_type = "n."
            logger.info(f"variant is a transcript variant {variant}")
            parts = str(hgvs_n).split(':n.')
            logger.info(f"Converted to n. coordinates: {hgvs_n}")

            # For intronic variants, we need to ensure need to do some extra processing
            if ((hgvs_n.posedit.pos.start.offset != 0 or hgvs_n.posedit.pos.end.offset != 0) and
                    genomic_reference is None):
                raise VariantFormatError(f"Intronic variants are currently not supported: {variant}")
            elif ((hgvs_n.posedit.pos.start.offset != 0 or hgvs_n.posedit.pos.end.offset != 0) and
                    genomic_reference is not None):
                try:
                    hgvs_genomic = validator.vm.t_to_g(hgvs_variant, genomic_reference, alt_aln_method=alt_aln_method)
                except Exception:
                    raise VariantFormatError(f"Unable to map intronic variant {variant} to "
                                             f"genomic reference {genomic_reference}. ")
                else:
                    converted_from_intronic = hgvs_variant.ac
                    logger.info(f"Converted to genomic coordinates: {hgvs_genomic}")
                    if re.search("ins$", str(hgvs_genomic)) and re.search("=$", str(hgvs_variant)):
                        hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.ref
                    parts = str(hgvs_genomic).split(':g.')
            logger.info(f"parts: {parts}")
        else:
            logger.error(f"Invalid variant format: {variant}")
            raise VariantFormatError("Variant must contain one of ':g.', ':c.', or ':n.'")

        if len(parts) != 2:
            logger.error(f"Variant splitting failed. Got parts: {parts}")
            raise VariantFormatError("Invalid variant format after split.")

        reference, seq_state_info = parts

        # Process insertions
        if "ins" in seq_state_info:
            position, sequence = seq_state_info.split('ins')
            variant_type = 'ins'
            start, end = map(int, position.split('_'))
            logger.info(f"Detected {variant_type} at position {start} to {end} with sequence {sequence}")
            if known_repeat_unit is None:
                repeated_unit = decipher_repeated_unit(sequence)
            else:
                if known_repeat_unit in sequence:
                    repeated_unit = known_repeat_unit
                else:
                    repeated_unit = validator.revcomp(known_repeat_unit)
            logger.info(f"Repeated unit: {repeated_unit}")
            reference_start = decipher_start_of_full_reference_repeated_sequence(reference, repeated_unit, start,
                                                                                 validator)
            reference_end = decipher_end_of_full_reference_repeated_sequence(reference, repeated_unit, start, validator)
            logger.info(f"Reference start position: {reference_start} and end position: {reference_end}")
            return reassemble_expanded_repeat_variant(reference, reference_start, reference_end, reference_type,
                                                      repeated_unit, sequence, validator,
                                                      converted_from_intronic=converted_from_intronic,
                                                      alt_aln_method=alt_aln_method)

        # Process identity cases (=)
        elif "=" in seq_state_info:
            seq_state_info = seq_state_info.replace('=', "")
            position, sequence = re.split(r'(?<=\d)(?=[A-Za-z])', seq_state_info)
            start, end = map(int, position.split('_'))
            if known_repeat_unit is None:
                repeated_unit = decipher_repeated_unit(sequence)
            else:
                if known_repeat_unit in sequence:
                    repeated_unit = known_repeat_unit
                else:
                    repeated_unit = validator.revcomp(known_repeat_unit)
            logger.info(f"Repeated unit: {repeated_unit}")
            reference_start = decipher_start_of_full_reference_repeated_sequence(reference, repeated_unit, start, validator)
            reference_end = decipher_end_of_full_reference_repeated_sequence(reference, repeated_unit, start, validator)

            logger.info(f"Converted from {converted_from_intronic} to expanded repeat variant: {reference_start}, {reference_end}, {reference_type}, {repeated_unit}")

            return reassemble_expanded_repeat_variant(reference, reference_start, reference_end, reference_type,
                                                      repeated_unit, "", validator,
                                                      converted_from_intronic=converted_from_intronic,
                                                      alt_aln_method=alt_aln_method)

        # Process duplications
        elif "dup" in seq_state_info:
            position, sequence = seq_state_info.split('dup')
            original_sequence = sequence
            sequence *= 2  # Expand the duplication so it resembles an insertion
            start, end = map(int, position.split('_'))
            if known_repeat_unit is None:
                repeated_unit = decipher_repeated_unit(sequence)
            else:
                if known_repeat_unit in sequence:
                    repeated_unit = known_repeat_unit
                else:
                    repeated_unit = validator.revcomp(known_repeat_unit)
            logger.info(f"Repeated unit: {repeated_unit}")
            reference_start = decipher_start_of_full_reference_repeated_sequence(reference, repeated_unit, start, validator)
            reference_end = decipher_end_of_full_reference_repeated_sequence(reference, repeated_unit, start, validator)
            return reassemble_expanded_repeat_variant(reference, reference_start, reference_end, reference_type,
                                                      repeated_unit, original_sequence, validator,
                                                      converted_from_intronic=converted_from_intronic,
                                                      alt_aln_method=alt_aln_method)

        # Process deletions
        elif "del" in seq_state_info:
            position, sequence = seq_state_info.split('del')
            start, end = map(int, position.split('_'))
            if known_repeat_unit is None:
                repeated_unit = decipher_repeated_unit(sequence)
            else:
                if known_repeat_unit in sequence:
                    repeated_unit = known_repeat_unit
                else:
                    repeated_unit = validator.revcomp(known_repeat_unit)
            logger.info(f"Repeated unit: {repeated_unit}")
            reference_start = decipher_start_of_full_reference_repeated_sequence(reference, repeated_unit, start, validator)
            reference_end = decipher_end_of_full_reference_repeated_sequence(reference, repeated_unit, start, validator)
            remove_units = len(sequence) // len(repeated_unit)
            return reassemble_expanded_repeat_variant(reference, reference_start, reference_end, reference_type,
                                                      repeated_unit, repeated_unit, validator,
                                                      remove_units=remove_units,
                                                      converted_from_intronic=converted_from_intronic,
                                                      alt_aln_method=alt_aln_method)

    raise VariantFormatError("Variant must be an HGVS-formatted string of type 'ins', 'del', 'dup', or '='.")


if __name__ == '__main__':
    import VariantValidator
    validator = VariantValidator.Validator()
    result = convert_seq_state_to_expanded_repeat("NM_002111.8:c.54_116=", validator, genomic_reference="NC_000004.11")
    print(result)


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
