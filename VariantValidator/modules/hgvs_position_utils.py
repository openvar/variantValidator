from vvhgvs.enums import Datum


# ---------------------------------------------------------------------------
# Exceptions and validation
# ---------------------------------------------------------------------------

class HgvsPositionException(Exception):
    """Raised when a position helper is given a non-HGVS object."""
    pass


def is_object(hgvs_object):
    """
    Check that the supplied value is a parsed HGVS object.

    Parameters
    ----------
    hgvs_object
        Object expected to be a parsed HGVS SequenceVariant.

    Raises
    ------
    HgvsPositionException
        If the supplied value does not have an HGVS accession attribute.
    """
    try:
        hgvs_object.ac
    except AttributeError:
        raise HgvsPositionException(
            f"Variant {hgvs_object} is not a parsed hgvs object"
        )


# ---------------------------------------------------------------------------
# Intronic position helpers
# ---------------------------------------------------------------------------

def start_position_is_intronic(hgvs_object):
    """
    Return whether the start position is intronic.
    """
    is_object(hgvs_object)
    return hgvs_object.posedit.pos.start.offset != 0


def end_position_is_intronic(hgvs_object):
    """
    Return whether the end position is intronic.
    """
    is_object(hgvs_object)
    return hgvs_object.posedit.pos.end.offset != 0


def both_positions_are_intronic(hgvs_object):
    """
    Return whether both the start and end positions are intronic.
    """
    return (
        start_position_is_intronic(hgvs_object)
        and end_position_is_intronic(hgvs_object)
    )


def either_position_is_intronic(hgvs_object):
    """
    Return whether either the start or end position is intronic.
    """
    return (
        start_position_is_intronic(hgvs_object)
        or end_position_is_intronic(hgvs_object)
    )


def only_one_position_is_intronic(hgvs_object):
    """
    Return whether exactly one of the start or end positions is intronic.
    """
    return (
        start_position_is_intronic(hgvs_object)
        != end_position_is_intronic(hgvs_object)
    )


# ---------------------------------------------------------------------------
# Positive intronic offset helpers
# ---------------------------------------------------------------------------

def start_offset_is_positive(hgvs_object):
    """
    Return whether the start position has a positive intronic offset.
    """
    is_object(hgvs_object)
    return hgvs_object.posedit.pos.start.offset > 0


def end_offset_is_positive(hgvs_object):
    """
    Return whether the end position has a positive intronic offset.
    """
    is_object(hgvs_object)
    return hgvs_object.posedit.pos.end.offset > 0


def both_offsets_are_positive(hgvs_object):
    """
    Return whether both start and end offsets are positive.
    """
    return (
        start_offset_is_positive(hgvs_object)
        and end_offset_is_positive(hgvs_object)
    )


def either_offset_is_positive(hgvs_object):
    """
    Return whether either the start or end offset is positive.
    """
    return (
        start_offset_is_positive(hgvs_object)
        or end_offset_is_positive(hgvs_object)
    )


# ---------------------------------------------------------------------------
# Negative intronic offset helpers
# ---------------------------------------------------------------------------

def start_offset_is_negative(hgvs_object):
    """
    Return whether the start position has a negative intronic offset.
    """
    is_object(hgvs_object)
    return hgvs_object.posedit.pos.start.offset < 0


def end_offset_is_negative(hgvs_object):
    """
    Return whether the end position has a negative intronic offset.
    """
    is_object(hgvs_object)
    return hgvs_object.posedit.pos.end.offset < 0


def both_offsets_are_negative(hgvs_object):
    """
    Return whether both start and end offsets are negative.
    """
    return (
        start_offset_is_negative(hgvs_object)
        and end_offset_is_negative(hgvs_object)
    )


def either_offset_is_negative(hgvs_object):
    """
    Return whether either the start or end offset is negative.
    """
    return (
        start_offset_is_negative(hgvs_object)
        or end_offset_is_negative(hgvs_object)
    )


# ---------------------------------------------------------------------------
# 5-prime UTR helpers
# ---------------------------------------------------------------------------

def start_is_5_prime_utr(hgvs_object):
    """
    Return whether the start position is in the 5-prime UTR.

    If the end position is in the 5-prime UTR, the start position will
    necessarily also be in the 5-prime UTR.
    """
    is_object(hgvs_object)
    return hgvs_object.posedit.pos.start.base < 0


def end_is_5_prime_utr(hgvs_object):
    """
    Return whether the end position is in the 5-prime UTR.

    If True, the complete interval is within the 5-prime UTR.
    """
    is_object(hgvs_object)
    return hgvs_object.posedit.pos.end.base < 0


# ---------------------------------------------------------------------------
# 3-prime UTR inspection helpers
# ---------------------------------------------------------------------------

def start_is_3_prime_utr(hgvs_object):
    """
    Return whether the start position is in the 3-prime UTR.

    If True, the complete interval is within the 3-prime UTR.
    """
    is_object(hgvs_object)
    return hgvs_object.posedit.pos.start.datum == Datum.CDS_END


def end_is_3_prime_utr(hgvs_object):
    """
    Return whether the end position is in the 3-prime UTR.

    If the start position is in the 3-prime UTR, the end position will
    necessarily also be in the 3-prime UTR.
    """
    is_object(hgvs_object)
    return hgvs_object.posedit.pos.end.datum == Datum.CDS_END


# ---------------------------------------------------------------------------
# 3-prime UTR modification helpers
# ---------------------------------------------------------------------------

def set_start_as_3_prime_utr(hgvs_object):
    """
    Set the start position datum to the CDS end.
    """
    is_object(hgvs_object)
    hgvs_object.posedit.pos.start.datum = Datum.CDS_END


def set_end_as_3_prime_utr(hgvs_object):
    """
    Set the end position datum to the CDS end.
    """
    is_object(hgvs_object)
    hgvs_object.posedit.pos.end.datum = Datum.CDS_END


def set_both_positions_as_3_prime_utr(hgvs_object):
    """
    Set both position datums to the CDS end.
    """
    set_start_as_3_prime_utr(hgvs_object)
    set_end_as_3_prime_utr(hgvs_object)


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
