import pytest
from VariantValidator.modules.use_checking import InvalidVariantError

from VariantValidator.modules.use_checking import (
    pre_parsing_global_common_mistakes,
)


class MockVariant:
    def __init__(self, text):
        self.original = text
        self.quibble = text
        self.warnings = []

    def format_quibble(self):
        return False

class InvalidFormatVariant(MockVariant):
    def format_quibble(self):
        return True


def test_numeric_only_variant():
    variant = MockVariant("123456")

    assert pre_parsing_global_common_mistakes(variant) is True
    assert any("only contains numeric characters" in w for w in variant.warnings)


def test_numeric_with_punctuation():
    variant = MockVariant("1:12345")

    assert pre_parsing_global_common_mistakes(variant) is True
    assert any("only contains numeric characters" in w for w in variant.warnings)


def test_alphanumeric_only_variant():
    variant = MockVariant("BRCA1")

    assert pre_parsing_global_common_mistakes(variant) is True
    assert any("only contains alphanumeric characters" in w for w in variant.warnings)


def test_concatenated_c_and_p_description():
    variant = MockVariant("NM_000001.1:c.123A>G:p.Gly41Arg")

    assert pre_parsing_global_common_mistakes(variant) is True
    assert any("concatenation" in w for w in variant.warnings)


def test_double_colon_with_extra_identifier():
    variant = MockVariant("NM_000001.1:abc:c.123A>G")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any("single colon" in w for w in variant.warnings)
    assert any("Illegal addition" in w for w in variant.warnings)

def test_old_style_substitution_format():
    variant = MockVariant("NM_000001.1:c.A123T")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any("not compliant" in w for w in variant.warnings)
    assert any("Did you mean" in w for w in variant.warnings)


def test_upper_case_reference_type_corrected():
    variant = InvalidFormatVariant("NM_000001.1:C.123A>G")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble == "NM_000001.1:c.123A>G"
    assert any("wrong case" in w for w in variant.warnings)


def test_missing_dot_after_reference_type():
    variant = InvalidFormatVariant("NM_000001.1:c123A>G")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any("lacks the ." in w for w in variant.warnings)


def test_bracketed_reference_is_stripped():
    variant = InvalidFormatVariant("(NM_000001.1):c.123A>G")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble == "NM_000001.1:c.123A>G"
    assert any("Stripping unnecessary characters" in w for w in variant.warnings)


def test_invalid_variant_raises():
    class InvalidVariant(MockVariant):
        def format_quibble(self):
            return True

    variant = InvalidVariant("not@a_variant")
    variant.quibble = "not@a_variant"

    with pytest.raises(InvalidVariantError):
        pre_parsing_global_common_mistakes(variant)

    assert any(
        "not in an accepted format" in warning
        for warning in variant.warnings
    )


def test_ins_without_sequence():
    variant = MockVariant("NM_000001.1:c.1ins")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any("inserted sequence" in w for w in variant.warnings)


def test_ins_length_parentheses():
    variant = MockVariant("NM_000001.1:c.1ins(10)")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any("N[10]" in w for w in variant.warnings)


def test_ins_length_plain_number():
    variant = MockVariant("NM_000001.1:c.1ins10")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any("N[10]" in w for w in variant.warnings)


def test_ins_uncertain_range():
    variant = MockVariant("NM_000001.1:c.1_2ins(10_20)")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any("N[(10_20)]" in w for w in variant.warnings)


def test_substitution_with_range():
    variant = MockVariant("NM_000001.1:c.1_2A>G")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any("Base substitution" in w for w in variant.warnings)


def test_ins_uncertain_range_reversed():
    variant = MockVariant("NM_000001.1:c.1_2ins[(20_10)]")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any("Please rewrite (20_10) to N[(10_20)]" in w for w in variant.warnings)


def test_ins_uncertain_range_equal():
    variant = MockVariant("NM_000001.1:c.1_2ins[(10_10)]")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any("Please rewrite (10_10) to N[(10)]" in w for w in variant.warnings)


def test_ins_uncertain_range_without_positions():
    variant = MockVariant("NM_000001.1:c.1ins[(10_20)]")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "An insertion must be provided with the two positions"
        in w for w in variant.warnings
    )


def test_ins_uncertain_range_valid():
    variant = MockVariant("NM_000001.1:c.1_2ins[(10_20)]")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "syntactically correct"
        in w for w in variant.warnings
    )


def test_sequence_repeat_expansion():
    variant = MockVariant("NM_000001.1:c.1_2insA[3]")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble.endswith("insAAA")
    assert any("may also be written as" in w for w in variant.warnings)


def test_sequence_repeat_expansion_multiple_sections():
    variant = MockVariant("NM_000001.1:c.1_2ins[A[2];C]")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble.endswith("insAAC")
    assert any("may also be written as" in w for w in variant.warnings)


def test_sequence_repeat_without_positions():
    variant = MockVariant("NM_000001.1:c.1insA[2]")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert any(
        "An insertion must be provided with the two positions"
        in w for w in variant.warnings
    )


def test_delins_repeat_expansion():
    variant = MockVariant("NM_000001.1:c.1_2delinsA[2]")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble.endswith("delinsAA")
    assert any("may also be written as" in w for w in variant.warnings)


def test_del_repeat_expansion():
    variant = MockVariant("NM_000001.1:c.1_2delA[2]")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble.endswith("delAA")
    assert any("may also be written as" in w for w in variant.warnings)

def test_repeat_expansion_delins_multiple_groups():
    variant = MockVariant("NM_000001.1:c.1_2delins[A[2];T[3]]")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble.endswith("delinsAATTT")
    assert any("may also be written as" in w for w in variant.warnings)


def test_repeat_expansion_ins_multiple_groups():
    variant = MockVariant("NM_000001.1:c.1_2ins[G[2];C[2]]")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble.endswith("insGGCC")
    assert any("may also be written as" in w for w in variant.warnings)


def test_repeat_expansion_del_multiple_groups():
    variant = MockVariant("NM_000001.1:c.1_2del[A[3];G]")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble.endswith("delAAAG")
    assert any("may also be written as" in w for w in variant.warnings)


def test_repeat_expansion_with_N():
    variant = MockVariant("NM_000001.1:c.1_2insN[4]")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble.endswith("insNNNN")


def test_repeat_expansion_single_base():
    variant = MockVariant("NM_000001.1:c.1_2delinsG[5]")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble.endswith("delinsGGGGG")


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
