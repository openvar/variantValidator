
def methyl_syntax(my_variant):
    """
    Remove and store a methylation syntax suffix before HGVS object parsing.

    :param my_variant:
    :return: Updated variant if methylation syntax is detected, otherwise None.
    """
    quibble = my_variant.quibble

    if "|" not in quibble:
        return None

    if "|gom" in quibble:
        my_variant.reformat_output = "|gom"
    elif "|lom" in quibble:
        my_variant.reformat_output = "|lom"
    elif "|met=" in quibble:
        my_variant.reformat_output = "|met="
    else:
        return None

    my_variant.quibble = quibble.split("|", 1)[0] + "="
    return my_variant


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
