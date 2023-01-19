# Errors
class FuzzyPositionError(Exception):
    pass


class FuzzyRangeError(Exception):
    pass


def fuzzy_ends(my_variant):
    """
    :param my_variant:
    :return: False if no fuzzy end detected otherwise raises Exception and provides information as to where the
    fuzzy end is located
    """

    if "?" in str(my_variant.hgvs_formatted.posedit.pos):
        if "?" in str(my_variant.hgvs_formatted.posedit.pos.end) and "?" not in str(
                my_variant.hgvs_formatted.posedit.pos.start):
            raise FuzzyPositionError("Fuzzy/unknown variant end position in submitted variant description")
        elif "?" in str(my_variant.hgvs_formatted.posedit.pos.start) and "?" not in str(
                my_variant.hgvs_formatted.posedit.pos.end):
            raise FuzzyPositionError("Fuzzy/unknown variant start position in submitted variant description")
        else:
            raise FuzzyPositionError("Fuzzy/unknown variant start and end positions "
                                     "in submitted variant description")

    return False

# <LICENSE>
# Copyright (C) 2016-2023 VariantValidator Contributors
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
