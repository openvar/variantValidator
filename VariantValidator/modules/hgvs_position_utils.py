from vvhgvs.enums import Datum

class HgvsPositionException(Exception):
    pass

def is_object(hgvs_object):
    try:
        hgvs_object.ac
    except AttributeError:
        raise HgvsPositionException(
            f"Variant {hgvs_object} is not a parsed hgvs object"
        )


def position_is_intronic(hgvs_object, check_start=False, check_end=False, check_start_and_end=False
):

    is_object(hgvs_object)

    if check_start_and_end:
        return (
            hgvs_object.posedit.pos.start.offset != 0
            and hgvs_object.posedit.pos.end.offset != 0
        )
    elif check_start:
        return hgvs_object.posedit.pos.start.offset != 0
    elif check_end:
        return hgvs_object.posedit.pos.end.offset != 0

    return False


def offset_is_positive(hgvs_object, check_start=False, check_end=False):

    is_object(hgvs_object)

    if check_start:
        return hgvs_object.posedit.pos.start.offset > 0
    elif check_end:
        return hgvs_object.posedit.pos.end.offset > 0

    return False


def offset_is_negative(hgvs_object, check_start=False, check_end=False):

    is_object(hgvs_object)

    if check_start:
        return hgvs_object.posedit.pos.start.offset < 0
    elif check_end:
        return hgvs_object.posedit.pos.end.offset < 0

    return False


def is_5_prime_utr(hgvs_object, check_start=False, check_end=False):

    is_object(hgvs_object)

    if check_start:
        return hgvs_object.posedit.pos.start.base < 0
    elif check_end:
        return hgvs_object.posedit.pos.end.base < 0

    return False


def is_3_prime_utr(hgvs_object, check_start=False, check_end=False):

    is_object(hgvs_object)

    if check_start:
        return hgvs_object.posedit.pos.start.datum == Datum.CDS_END
    elif check_end:
        return hgvs_object.posedit.pos.end.datum == Datum.CDS_END

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
