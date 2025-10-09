import logging
import re
from VariantValidator.modules import utils
from VariantValidator.modules.hgvs_utils import hgvs_delins_parts_to_hgvs_obj, to_vv_hgvs
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

    # Begin constructing the variant with positional metadata, the vv_hgvs code treats all non dup
    # types as "delins" (deliberately not "indel" due to it's more official usage)
    variant = hgvs_delins_parts_to_hgvs_obj(reference,reference_type, reference_start, '', '',end=reference_end)

    # If the variant is coding normalise to ensure alignment
    if reference_type in ["c","n"] and "NC_" not in reference:
        variant = hgvs_delins_parts_to_hgvs_obj(reference,'n', reference_start, '', '',end=reference_end)
        if reference_type == "c":
            variant = validator.vm.n_to_c(variant)
        if reference.startswith('ENS'):
            variant_n = validator.genebuild_normalizer.normalize(variant)
        else:
            variant_n = validator.splign_normalizer.normalize(variant)
        assert variant_n.posedit.edit.ref == variant.posedit.edit.ref
    elif converted_from_intronic is not False:
        variant.type = 'g'
        variant = validator.vm.g_to_t(variant, converted_from_intronic,alt_aln_method=alt_aln_method)
        orientation = validator.hdp.get_tx_exons(
                converted_from_intronic, reference, alt_aln_method=alt_aln_method)[0][3]
        if orientation != 1:
            # If the transcript is on the reverse strand, we need reverse complement the repeat bases
            repeated_unit = validator.revcomp(repeated_unit)
    # Switch from = type use for mapping checks, final version only, the expanded repeat type output will
    # lose annotation on normalisation (this avoids re-parsing on VRS output production etc)
    variant = to_vv_hgvs(variant) # restore vv type extras if lost in mapping etc
    variant.posedit.edit.alt = repeated_unit * total_units
    variant.posedit.expanded_rep = repeated_unit
    return variant

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
    This function currently forces at least one block of a repeat to exist, or
    it returns None!
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
    Now assumes that the previous section is valid, so that it can start from the end of a block.
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

    return start

def convert_seq_state_to_expanded_repeat(variant, validator, genomic_reference=None, known_repeat_unit=None):
    """
    Main interface: convert a HGVS variant string to an expanded repeat representation.

    Supports: insertions (ins), deletions (del), duplications (dup), identity (=)
    """
    # We should only get variant data in as a hgvs object, but null variants may be None or ''
    # Return the same output as given input in this case
    if not variant:
        return variant
    if not variant.posedit.edit.type in ['ins', 'del', 'dup', 'identity']:
        raise VariantFormatError("Variant must be a HGVS format object of a equal, or length "
                                 "change type, (i.e. 'ins', 'del', 'dup', or '=').")

    # Use ENST normalizer or splign depending on variant format
    hn = validator.genebuild_normalizer if "ENST" in variant.ac else validator.splign_normalizer
    alt_aln_method = "genebuild" if "ENST" in variant.ac else "splign"
    variant = hn.normalize(variant)

    # Determine reference type and extract sequence info
    converted_from_intronic = False
    working_hgvs = False # this will be a G variant for intronic, if fixed
    hgvs_genomic = False # store the genomic map for mapping change detection
    if variant.type == 'g':
        reference_type = variant.type
        working_hgvs = variant
    elif variant.type in ["c", "n"]:
        if variant.type == 'c':
            reference_type = variant.type
            hgvs_n = validator.vm.c_to_n(variant)
            working_hgvs = hgvs_n
        else:
            hgvs_n = variant
            reference_type = variant.type
            working_hgvs = variant
        logger.info(f"variant is a transcript variant {str(variant)}")
        logger.info(f"Converted to n. coordinates: {hgvs_n}")
        # For intronic variants, we need to ensure need to do some extra processing
        if ((hgvs_n.posedit.pos.start.offset != 0 or hgvs_n.posedit.pos.end.offset != 0) and
                genomic_reference is None):
            raise VariantFormatError(f"Intronic variants are currently not supported: {variant}")
        elif ((hgvs_n.posedit.pos.start.offset != 0 or hgvs_n.posedit.pos.end.offset != 0) and
                genomic_reference is not None):
            try:
                hgvs_genomic = validator.vm.t_to_g(variant, genomic_reference, alt_aln_method=alt_aln_method)
            except Exception:
                raise VariantFormatError(f"Unable to map intronic variant {variant} to "
                                         f"genomic reference {genomic_reference}. ")
            else:
                converted_from_intronic = variant.ac
                logger.info(f"Converted to genomic coordinates: {hgvs_genomic}")
                if re.search("ins$", str(hgvs_genomic)) and variant.posedit.edit.type == 'identity':
                    hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.ref
                working_hgvs = hgvs_genomic
    else:
        logger.error(f"Invalid variant format: {variant}")
        raise VariantFormatError("Variant must contain one of ':g.', ':c.', or ':n.'")

    # set errors possibly triggered from multiple locations
    ins_not_repeat_err = "At least one repeat must exist in the genome for "\
             "HGVS Repeated Sequences"

    # derive repeat unit if none given (using ref sequence, except for del type)
    if known_repeat_unit is None:
        if working_hgvs.posedit.edit.ref:
            repeated_unit = decipher_repeated_unit(working_hgvs.posedit.edit.ref)
        else:
            repeated_unit = decipher_repeated_unit(working_hgvs.posedit.edit.alt)
    else:
        sequence = working_hgvs.posedit.edit.ref
        # ins normalises to dup, so long as at == number of repeats exist in the genomic
        # flank, so this should be rare
        if not sequence:
            sequence = working_hgvs.posedit.edit.alt
        if known_repeat_unit in sequence:
            repeated_unit = known_repeat_unit
        else:
            repeated_unit = validator.revcomp(known_repeat_unit)
    logger.info(f"Repeated unit: {repeated_unit}")

    # get ref and type (to avoid repeated lookup)
    reference = working_hgvs.ac
    edit_type = working_hgvs.posedit.edit.type

    # test for variants that have been converted into other sequence alterations,
    # or now contain non-repeat based sequence, rather than even length changes
    # This is sometimes an expected outcome with != alignments
    if hgvs_genomic and hgvs_genomic.posedit.edit.ref != variant.posedit.edit.ref:
        # There could be some corner cases here where the edge sequence makes up
        # for uneven deletions in the mapping but these should normalise out, so
        # we return on a simple condition for now
        if edit_type == 'delins' or \
            len(hgvs_genomic.posedit.edit.ref) % len(repeated_unit) or \
            int( len(hgvs_genomic.posedit.edit.ref)/len(repeated_unit)
            ) * repeated_unit != hgvs_genomic.posedit.edit.ref:
            err_str ="Variant format no longer valid for repeat after map to "+\
                f"{hgvs_genomic} (from {variant})"
            logger.info(err_str)
            raise VariantFormatError(err_str)
    elif  working_hgvs and working_hgvs.posedit.edit.ref and \
            working_hgvs.posedit.edit.ref[:len(repeated_unit)] != repeated_unit:
        err_str =f"Variant format not valid for repeat {repeated_unit} "\
                f"(from {variant}) this should only happen for mapped "\
                "consequences of expanded repeat input over regions of "\
                "alignment mismatch"
        logger.info(err_str)
        raise VariantFormatError(err_str)


    # Process insertions (normalisation should currently turn some of these into dups)
    # a non normalised ins can happen internally, or at one of the ends, end or start
    # ins locs just need -1 to the starting base, but the logic of the start finding
    # requires us to start -1 repeat unit of the end if is at the end
    if edit_type == "ins":
        logger.info(f"Detected ins at position {str(working_hgvs.posedit.pos)} "+\
                "with sequence {working_hgvs.posedit.edit.alt}")
        reference_start = decipher_start_of_full_reference_repeated_sequence(
                reference, repeated_unit, working_hgvs.posedit.pos.start.base + 1,validator)

        if reference_start is None:
            reference_start = decipher_start_of_full_reference_repeated_sequence(
                reference, repeated_unit,
                working_hgvs.posedit.pos.start.base - len(repeated_unit) + 1,
                validator)
            # if fall back fails then no repeats exist in genome at this loc, error out
            if reference_start is None:
                raise VariantFormatError(ins_not_repeat_err)
        reference_end = decipher_end_of_full_reference_repeated_sequence(
                reference, repeated_unit, working_hgvs.posedit.pos.end.base,validator)
        logger.info(f"Reference start position: {reference_start} and end position: {reference_end}")
        return reassemble_expanded_repeat_variant(reference, reference_start, reference_end, reference_type,
                                                  repeated_unit, working_hgvs.posedit.edit.alt, validator,
                                                  converted_from_intronic=converted_from_intronic,
                                                  alt_aln_method=alt_aln_method)

    # Process identity cases (=)
    elif edit_type == "identity":
        reference_start = decipher_start_of_full_reference_repeated_sequence(
                reference, repeated_unit, working_hgvs.posedit.pos.start.base, validator)
        reference_end = decipher_end_of_full_reference_repeated_sequence(
                reference, repeated_unit, working_hgvs.posedit.pos.end.base, validator)
        logger.info(f"Converted from {converted_from_intronic} to expanded repeat variant: {reference_start}, {reference_end}, {reference_type}, {repeated_unit} via {str(working_hgvs)}")

        return reassemble_expanded_repeat_variant(reference, reference_start, reference_end, reference_type,
                                                  repeated_unit, "", validator,
                                                  converted_from_intronic=converted_from_intronic,
                                                  alt_aln_method=alt_aln_method)

    # Process duplications
    elif edit_type == "dup":
        original_sequence = working_hgvs.posedit.edit.ref
        reference_start = decipher_start_of_full_reference_repeated_sequence(
                reference, repeated_unit, working_hgvs.posedit.pos.start.base, validator)
        reference_end = decipher_end_of_full_reference_repeated_sequence(
                reference, repeated_unit, working_hgvs.posedit.pos.end.base, validator)
        return reassemble_expanded_repeat_variant(reference, reference_start, reference_end, reference_type,
                                                  repeated_unit, original_sequence, validator,
                                                  converted_from_intronic=converted_from_intronic,
                                                  alt_aln_method=alt_aln_method)

    # Process deletions
    elif edit_type == "del":
        reference_start = decipher_start_of_full_reference_repeated_sequence(
                reference, repeated_unit, working_hgvs.posedit.pos.start.base, validator)
        reference_end = decipher_end_of_full_reference_repeated_sequence(
                reference, repeated_unit, working_hgvs.posedit.pos.end.base, validator)
        remove_units = len(working_hgvs.posedit.edit.ref) // len(repeated_unit)
        return reassemble_expanded_repeat_variant(reference, reference_start, reference_end, reference_type,
                                                  repeated_unit, repeated_unit, validator,
                                                  remove_units=remove_units,
                                                  converted_from_intronic=converted_from_intronic,
                                                  alt_aln_method=alt_aln_method)


def quick_testfunc():
    import VariantValidator
    validator = VariantValidator.Validator()
    variant = validator.hp.parse("NM_002111.8:c.54_116=")
    result = convert_seq_state_to_expanded_repeat(variant, validator, genomic_reference="NC_000004.11")
    print(result)
    return result

if __name__ == '__main__': quick_testfunc()


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
