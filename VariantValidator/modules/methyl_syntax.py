
def methyl_syntax(my_variant):
    """
    :param my_variant:
    :return: False if no | is detected otherwise raises Exception and provides information as to where the
    | is located
    """

    if "|" in my_variant.quibble:
        if "gom" in my_variant.quibble or "lom" in my_variant.quibble or "met=" in my_variant.quibble:
            if "|gom" in my_variant.quibble:
                my_variant.reformat_output = "|gom"
            if "|lom" in my_variant.quibble:
                my_variant.reformat_output = "|lom"
            if "|met=" in my_variant.quibble:
                my_variant.reformat_output = "|met="
            met_var = str(my_variant.quibble.split("|")[0]) + "="
            my_variant.quibble = met_var
            return my_variant


# <LICENSE>
# Copyright (C) 2016-2024 VariantValidator Contributors
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
