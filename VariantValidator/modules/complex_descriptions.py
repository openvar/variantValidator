import re
import vvhgvs.exceptions

# Errors
class FuzzyPositionError(Exception):
    pass


class FuzzyRangeError(Exception):
    pass


class InvalidRangeError(Exception):
    pass


class IncompatibleTypeError(Exception):
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


def uncertain_positions(my_variant, validator):

    print("AWOOO")
    print(dir(my_variant))
    print(my_variant.quibble)

    print("Try it")

    # Formats like NC_000005.9:g.(90136803_90144453)_(90159675_90261231)dup
    if ")_(" in my_variant.quibble:
        accession_and_type, positions_and_edit = my_variant.quibble.split(".(")
        position_1, position_2 = positions_and_edit.split(")_(")
        position_2, variation = position_2.split(")")
        position_1 = position_1.replace(")", "")
        position_2 = position_2.replace("(", "")
        print(accession_and_type)
        print(position_1)
        print(position_2)
        print(variation)
        v1 = f"{accession_and_type}.{position_1}{variation}"
        v2 = f"{accession_and_type}.{position_2}{variation}"
        my_variant.reftype = f":{accession_and_type.split(':')[1]}."
        print(my_variant.reftype)
        try:
            parsed_v1 = validator.hp.parse(v1)
            validator.vr.validate(parsed_v1)
        except vvhgvs.exceptions.HGVSError as e:
            if "is not known to be compatible with variant type" in str(e):
                raise IncompatibleTypeError(str(e))
            import traceback
            traceback.print_exc()
            raise InvalidRangeError(f"{position_1} is an invlaid range for "
                                    f"accession {accession_and_type.split(':')[0]}")
        try:
            parsed_v2 = validator.hp.parse(v2)
            validator.vr.validate(parsed_v2)
        except vvhgvs.exceptions.HGVSError as e:
            if "is not known to be compatible with variant type" in str(e):
                raise IncompatibleTypeError(str(e))
            import traceback
            raise InvalidRangeError(f"{position_1} is an invlaid range for "
                                    f"accession {accession_and_type.split(':')[0]}")

        my_variant.warnings = ["Uncertain positions are not fully supported, however the syntax is valid"]
        if "NC_" in my_variant.quibble:
            my_variant.hgvs_genomic = my_variant.quibble
        if "NM_" in my_variant.quibble or "NR_" in my_variant.quibble or "ENST" in my_variant.quibble:
            my_variant.hgvs_coding = my_variant.quibble
            my_variant.hgvs_transcript_variant = my_variant.quibble


    return

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
