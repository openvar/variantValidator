import logging

from VariantValidator.modules.hgvs_utils import (
    hgvs_delins_parts_to_hgvs_obj,
    to_vv_hgvs,
)


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


def reassemble_expanded_repeat_variant(
        reference,
        reference_start,
        reference_end,
        reference_type,
        repeated_unit,
        inserted_sequence,
        validator,
        remove_units=0,
        converted_from_intronic=False,
        alt_aln_method=None,
):
    """
    Assemble the final expanded repeat representation in HGVS-like format.

    Example: NG_012232.1:g.3_6T[20]

    Args explain how all parts are used to compute total repeat count.
    """
    if not repeated_unit:
        raise RepeatedUnitError(
            "Repeated unit must be a non-empty string."
        )

    unit_len = len(repeated_unit)
    inserted_repeat_units = len(inserted_sequence) // unit_len
    reference_repeat_units = (
        reference_end - reference_start + 1
    ) // unit_len
    total_units = (
        reference_repeat_units
        + inserted_repeat_units
        - remove_units
    )

    # Begin constructing the variant with positional metadata. The vv_hgvs
    # code treats all non-dup types as "delins" (deliberately not "indel"
    # due to its more official usage).
    variant = hgvs_delins_parts_to_hgvs_obj(
        reference,
        reference_type,
        reference_start,
        "",
        "",
        end=reference_end,
    )

    # If the variant is coding, normalise to ensure alignment.
    if (
            reference_type in ("c", "n")
            and not reference.startswith("NC_")
    ):
        variant = hgvs_delins_parts_to_hgvs_obj(
            reference,
            "n",
            reference_start,
            "",
            "",
            end=reference_end,
        )

        if reference_type == "c":
            variant = validator.vm.n_to_c(variant)

        if reference.startswith("ENS"):
            variant_n = validator.genebuild_normalizer.normalize(
                variant
            )
        else:
            variant_n = validator.splign_normalizer.normalize(
                variant
            )

        assert (
            variant_n.posedit.edit.ref
            == variant.posedit.edit.ref
        )

    elif converted_from_intronic:
        variant.type = "g"
        variant = validator.vm.g_to_t(
            variant,
            converted_from_intronic,
            alt_aln_method=alt_aln_method,
        )

        orientation = validator.hdp.get_tx_exons(
            converted_from_intronic,
            reference,
            alt_aln_method=alt_aln_method,
        )[0][3]

        if orientation != 1:
            # If the transcript is on the reverse strand, reverse
            # complement the repeat bases.
            repeated_unit = validator.revcomp(repeated_unit)

    # Switch from = type used for mapping checks to the final version.
    # Expanded repeat output loses annotation on normalisation, so restore
    # the VV type extras here rather than re-parsing during VRS output.
    variant = to_vv_hgvs(variant)
    variant.posedit.edit.alt = repeated_unit * total_units
    variant.posedit.expanded_rep = repeated_unit

    return variant


def decipher_repeated_unit(sequence):
    """
    Detect the smallest unit that repeats to form the given DNA sequence.

    If no repeat unit is found, return the original sequence.
    """
    sequence_len = len(sequence)

    for unit_len in range(1, sequence_len // 2 + 1):
        unit = sequence[:unit_len]
        multiplier = sequence_len // unit_len

        if (
                unit * multiplier
                == sequence[:unit_len * multiplier]
        ):
            return unit

    return sequence


def decipher_start_of_full_reference_repeated_sequence(
        reference,
        repeated_unit,
        start,
        validator,
        window_size=100,
):
    """
    Given a position where a repeat starts, fetch a larger window of
    sequence and walk backwards in-memory to find the full beginning of
    the repeat block.

    This function currently forces at least one block of a repeat to
    exist, or it returns None.
    """
    if not repeated_unit:
        raise RepeatedUnitError(
            "Repeated unit must be a non-empty string."
        )

    if not isinstance(start, int) or start <= 0:
        raise StartPositionError(
            "Start position must be a positive integer."
        )

    unit_len = len(repeated_unit)
    local_start = start - 1

    # Determine how far back we can fetch.
    fetch_start = max(0, local_start - window_size)
    fetch_end = local_start + unit_len

    # Fetch once.
    seq_window = validator.sf.fetch_seq(
        reference,
        start_i=fetch_start,
        end_i=fetch_end,
    )

    # Position of the unit in the fetched window.
    local_pos = local_start - fetch_start
    current_chunk = seq_window[
        local_pos:local_pos + unit_len
    ]

    if current_chunk != repeated_unit:
        return None

    # Walk backwards within the in-memory window.
    while local_pos - unit_len >= 0:
        prev_chunk = seq_window[
            local_pos - unit_len:local_pos
        ]

        if prev_chunk != repeated_unit:
            break

        local_pos -= unit_len

    # Convert back to reference coordinate (1-based).
    return fetch_start + local_pos + 1


def decipher_end_of_full_reference_repeated_sequence(
        reference,
        repeated_unit,
        start,
        validator,
        window_size=100,
):
    """
    Walk forward from a given repeat position to find the full extent of
    the repeated block.

    Fetch a larger window once to avoid multiple fetch_seq calls. This
    assumes that the previous section is valid, so it can start from the
    end of a block.
    """
    if not repeated_unit:
        raise RepeatedUnitError(
            "Repeated unit must be a non-empty string."
        )

    if not isinstance(start, int) or start <= 0:
        raise StartPositionError(
            "Start position must be a positive integer."
        )

    unit_len = len(repeated_unit)

    # Try two offsets for frame alignment.
    for start_0 in (start - 1 - unit_len, start - 1):
        if start_0 < 0:
            continue

        # Fetch a single forward window.
        fetch_end = start_0 + window_size

        try:
            seq_window = validator.sf.fetch_seq(
                reference,
                start_i=start_0,
                end_i=fetch_end,
            )
        except Exception:
            continue

        # Check if the repeated unit matches at this start.
        if seq_window[:unit_len] != repeated_unit:
            continue

        local_pos = 0

        # Walk forward in-memory.
        while local_pos + unit_len <= len(seq_window) - unit_len:
            next_chunk = seq_window[
                local_pos + unit_len:
                local_pos + (2 * unit_len)
            ]

            if next_chunk != repeated_unit:
                break

            local_pos += unit_len

        # Convert to 1-based end position.
        return start_0 + local_pos + unit_len

    return start


def convert_seq_state_to_expanded_repeat(
        variant,
        validator,
        genomic_reference=None,
        known_repeat_unit=None,
):
    """
    Main interface: convert a HGVS variant object to an expanded repeat
    representation.

    Supports: insertions (ins), deletions (del), duplications (dup),
    identity (=).
    """
    # We should only get variant data in as an HGVS object, but null
    # variants may be None or ''. Return the same output in this case.
    if not variant:
        return variant

    if variant.posedit.edit.type not in (
            "ins",
            "del",
            "dup",
            "identity",
    ):
        raise VariantFormatError(
            "Variant must be a HGVS format object of a equal, or length "
            "change type, (i.e. 'ins', 'del', 'dup', or '=')."
        )

    # Use ENST normalizer or splign depending on variant format.
    if variant.ac.startswith("ENST"):
        hn = validator.genebuild_normalizer
        alt_aln_method = "genebuild"
    else:
        hn = validator.splign_normalizer
        alt_aln_method = "splign"

    variant = hn.normalize(variant)

    # Determine reference type and extract sequence information.
    converted_from_intronic = False
    working_hgvs = False
    hgvs_genomic = False

    if variant.type == "g":
        reference_type = variant.type
        working_hgvs = variant

    elif variant.type in ("c", "n"):
        reference_type = variant.type

        if variant.type == "c":
            hgvs_n = validator.vm.c_to_n(variant)
            working_hgvs = hgvs_n
        else:
            hgvs_n = variant
            working_hgvs = variant

        logger.info(
            "variant is a transcript variant %s",
            variant,
        )
        logger.info(
            "Converted to n. coordinates: %s",
            hgvs_n,
        )

        # Intronic variants require genomic mapping for the additional
        # processing required to identify the repeat.
        is_intronic = (
            hgvs_n.posedit.pos.start.offset != 0
            or hgvs_n.posedit.pos.end.offset != 0
        )

        if is_intronic and genomic_reference is None:
            raise VariantFormatError(
                "Intronic variants are currently not supported: "
                f"{variant}"
            )

        if is_intronic:
            try:
                hgvs_genomic = validator.vm.t_to_g(
                    variant,
                    genomic_reference,
                    alt_aln_method=alt_aln_method,
                )
            except Exception:
                raise VariantFormatError(
                    f"Unable to map intronic variant {variant} to "
                    f"genomic reference {genomic_reference}. "
                )
            else:
                converted_from_intronic = variant.ac

                logger.info(
                    "Converted to genomic coordinates: %s",
                    hgvs_genomic,
                )

                # Identity variants can map across alignment differences
                # as a delins with reference sequence but an empty alt.
                # The previous implementation detected this by checking
                # whether the rendered HGVS description ended in "ins".
                # Inspect the HGVS edit object directly instead.
                genomic_edit = hgvs_genomic.posedit.edit
                if (
                        genomic_edit.type == "delins"
                        and not genomic_edit.alt
                        and variant.posedit.edit.type == "identity"
                ):
                    genomic_edit.alt = genomic_edit.ref

                working_hgvs = hgvs_genomic

    else:
        logger.error(
            "Invalid variant format: %s",
            variant,
        )
        raise VariantFormatError(
            "Variant must contain one of ':g.', ':c.', or ':n.'"
        )

    ins_not_repeat_err = (
        "At least one repeat must exist in the genome for "
        "HGVS Repeated Sequences"
    )

    # Derive the repeat unit if none is supplied.
    if known_repeat_unit is None:
        if working_hgvs.posedit.edit.ref:
            repeated_unit = decipher_repeated_unit(
                working_hgvs.posedit.edit.ref
            )
        else:
            repeated_unit = decipher_repeated_unit(
                working_hgvs.posedit.edit.alt
            )
    else:
        sequence = working_hgvs.posedit.edit.ref

        # ins normalises to dup so long as the appropriate number of
        # repeats exist in the genomic flank. This should therefore be
        # rare.
        if not sequence:
            sequence = working_hgvs.posedit.edit.alt

        if known_repeat_unit in sequence:
            repeated_unit = known_repeat_unit
        else:
            repeated_unit = validator.revcomp(
                known_repeat_unit
            )

    logger.info(
        "Repeated unit: %s",
        repeated_unit,
    )

    # Store frequently accessed values.
    reference = working_hgvs.ac
    edit = working_hgvs.posedit.edit
    edit_type = edit.type
    unit_len = len(repeated_unit)

    # Test for variants converted into other sequence alterations or
    # containing non-repeat sequence rather than even length changes.
    # This is sometimes an expected outcome with != alignments.
    if (
            hgvs_genomic
            and hgvs_genomic.posedit.edit.ref
            != variant.posedit.edit.ref
    ):
        genomic_ref = hgvs_genomic.posedit.edit.ref

        # There could be corner cases where edge sequence compensates for
        # uneven deletions during mapping, but these should normalise out.
        if (
                edit_type == "delins"
                or len(genomic_ref) % unit_len
                or (
                    (len(genomic_ref) // unit_len) * repeated_unit
                    != genomic_ref
                )
        ):
            err_str = (
                "Variant format no longer valid for repeat after map to "
                f"{hgvs_genomic} (from {variant})"
            )
            logger.info(err_str)
            raise VariantFormatError(err_str)

    elif (
            working_hgvs
            and edit.ref
            and edit.ref[:unit_len] != repeated_unit
    ):
        err_str = (
            f"Variant format not valid for repeat {repeated_unit} "
            f"(from {variant}) this should only happen for mapped "
            "consequences of expanded repeat input over regions of "
            "alignment mismatch"
        )
        logger.info(err_str)
        raise VariantFormatError(err_str)

    # Process insertions. Normalisation should currently turn some of
    # these into duplications. A non-normalised insertion can happen
    # internally or at one of the ends.
    if edit_type == "ins":
        logger.info(
            "Detected ins at position %s with sequence %s",
            working_hgvs.posedit.pos,
            edit.alt,
        )

        reference_start = (
            decipher_start_of_full_reference_repeated_sequence(
                reference,
                repeated_unit,
                working_hgvs.posedit.pos.start.base + 1,
                validator,
            )
        )

        if reference_start is None:
            reference_start = (
                decipher_start_of_full_reference_repeated_sequence(
                    reference,
                    repeated_unit,
                    (
                        working_hgvs.posedit.pos.start.base
                        - unit_len
                        + 1
                    ),
                    validator,
                )
            )

            # If fallback fails then no repeats exist in the genome at
            # this location.
            if reference_start is None:
                raise VariantFormatError(ins_not_repeat_err)

        reference_end = (
            decipher_end_of_full_reference_repeated_sequence(
                reference,
                repeated_unit,
                working_hgvs.posedit.pos.end.base,
                validator,
            )
        )

        logger.info(
            "Reference start position: %s and end position: %s",
            reference_start,
            reference_end,
        )

        return reassemble_expanded_repeat_variant(
            reference,
            reference_start,
            reference_end,
            reference_type,
            repeated_unit,
            edit.alt,
            validator,
            converted_from_intronic=converted_from_intronic,
            alt_aln_method=alt_aln_method,
        )

    # Process identity cases (=).
    if edit_type == "identity":
        reference_start = (
            decipher_start_of_full_reference_repeated_sequence(
                reference,
                repeated_unit,
                working_hgvs.posedit.pos.start.base,
                validator,
            )
        )
        reference_end = (
            decipher_end_of_full_reference_repeated_sequence(
                reference,
                repeated_unit,
                working_hgvs.posedit.pos.end.base,
                validator,
            )
        )

        logger.info(
            "Converted from %s to expanded repeat variant: "
            "%s, %s, %s, %s via %s",
            converted_from_intronic,
            reference_start,
            reference_end,
            reference_type,
            repeated_unit,
            working_hgvs,
        )

        return reassemble_expanded_repeat_variant(
            reference,
            reference_start,
            reference_end,
            reference_type,
            repeated_unit,
            "",
            validator,
            converted_from_intronic=converted_from_intronic,
            alt_aln_method=alt_aln_method,
        )

    # Process duplications.
    if edit_type == "dup":
        reference_start = (
            decipher_start_of_full_reference_repeated_sequence(
                reference,
                repeated_unit,
                working_hgvs.posedit.pos.start.base,
                validator,
            )
        )
        reference_end = (
            decipher_end_of_full_reference_repeated_sequence(
                reference,
                repeated_unit,
                working_hgvs.posedit.pos.end.base,
                validator,
            )
        )

        return reassemble_expanded_repeat_variant(
            reference,
            reference_start,
            reference_end,
            reference_type,
            repeated_unit,
            edit.ref,
            validator,
            converted_from_intronic=converted_from_intronic,
            alt_aln_method=alt_aln_method,
        )

    # Process deletions.
    if edit_type == "del":
        reference_start = (
            decipher_start_of_full_reference_repeated_sequence(
                reference,
                repeated_unit,
                working_hgvs.posedit.pos.start.base,
                validator,
            )
        )
        reference_end = (
            decipher_end_of_full_reference_repeated_sequence(
                reference,
                repeated_unit,
                working_hgvs.posedit.pos.end.base,
                validator,
            )
        )
        remove_units = len(edit.ref) // unit_len

        return reassemble_expanded_repeat_variant(
            reference,
            reference_start,
            reference_end,
            reference_type,
            repeated_unit,
            repeated_unit,
            validator,
            remove_units=remove_units,
            converted_from_intronic=converted_from_intronic,
            alt_aln_method=alt_aln_method,
        )


# def quick_testfunc():
#     import VariantValidator
#     validator = VariantValidator.Validator()
#     variant = validator.hp.parse("NM_002111.8:c.54_116=")
#     result = convert_seq_state_to_expanded_repeat(
#         variant,
#         validator,
#         genomic_reference="NC_000004.11"
#     )
#     print(result)
#     return result

# if __name__ == '__main__': quick_testfunc()


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
