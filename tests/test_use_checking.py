import pytest
from unittest.mock import MagicMock
from VariantValidator.modules.use_checking import (
    InvalidVariantError,
    pre_parsing_global_common_mistakes,
)
from VariantValidator.modules.use_checking import refseq_common_mistakes
from VariantValidator.modules.use_checking import (
    structure_checks_g,
    structure_checks_c,
    structure_checks_n,
)
import vvhgvs
from types import SimpleNamespace


class MockRefVariant:
    def __init__(self, quibble, reftype, transcript_type="c"):
        self.quibble = quibble
        self.reftype = reftype
        self.transcript_type = transcript_type
        self.warnings = []


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

def test_refseq_nm_used_as_genomic():
    v = MockRefVariant("NM_000001.1:g.1A>G", ":g.", "c")

    assert refseq_common_mistakes(v) is True
    assert any("Did you mean NM_000001.1:c." in w for w in v.warnings)


def test_refseq_nr_used_as_coding():
    v = MockRefVariant("NR_000001.1:c.1A>G", ":c.", "n")

    assert refseq_common_mistakes(v) is True
    assert any(":n." in w for w in v.warnings)


def test_protein_reference_as_nucleotide():
    v = MockRefVariant("NP_000001.1:c.1A>G", ":c.")

    assert refseq_common_mistakes(v) is True
    assert any("Protein reference sequence" in w for w in v.warnings)


def test_coding_transcript_as_noncoding():
    v = MockRefVariant("NM_000001.1:n.1A>G", ":n.", "c")

    assert refseq_common_mistakes(v) is True
    assert any("Did you mean NM_000001.1:c." in w for w in v.warnings)


def test_nucleotide_reference_as_protein():
    v = MockRefVariant("NM_000001.1:p.Gly1Val", ":p.")

    assert refseq_common_mistakes(v) is True
    assert any("protein-level" in w.lower() for w in v.warnings)


def test_ng_c_description():
    v = MockRefVariant("NG_000001.1:c.123A>G", ":c.")

    assert refseq_common_mistakes(v) is True
    assert len(v.warnings) == 2


def test_refseq_common_mistakes_valid():
    v = MockRefVariant("NM_000001.1:c.123A>G", ":c.", "c")

    assert refseq_common_mistakes(v) is False
    assert v.warnings == []

def test_refseq_enst_used_as_genomic():
    v = MockRefVariant("ENST000001.1:g.1A>G", ":g.", "c")

    assert refseq_common_mistakes(v) is True
    assert any("Did you mean" in w for w in v.warnings)


def test_refseq_ensp_used_as_nucleotide():
    v = MockRefVariant("ENSP000001.1:c.1A>G", ":c.")

    assert refseq_common_mistakes(v) is True
    assert any("Protein reference sequence" in w for w in v.warnings)


def test_refseq_ng_as_protein():
    v = MockRefVariant("NG_000001.1:p.Arg1Gly", ":p.")

    assert refseq_common_mistakes(v) is True
    assert any("protein-level" in w.lower() for w in v.warnings)


def test_refseq_nc_as_coding():
    v = MockRefVariant("NC_000001.11:c.123A>G", ":c.")

    assert refseq_common_mistakes(v) is True
    assert len(v.warnings) == 2


def test_refseq_nr_as_genomic():
    v = MockRefVariant("NR_000001.1:g.10A>G", ":g.", "n")

    assert refseq_common_mistakes(v) is True
    assert any(":n." in w for w in v.warnings)


def test_variant_comma_after_type_autocorrect():
    variant = InvalidFormatVariant("NM_000001.1:c,123A>G")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble == "NM_000001.1:c.123A>G"
    assert any("auto" in w.lower() for w in variant.warnings)


def test_bracketed_enst_reference():
    variant = InvalidFormatVariant("(ENST00000357654.9):c.123A>G")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble == "ENST00000357654.9:c.123A>G"


def test_bracketed_lrg_reference():
    variant = InvalidFormatVariant("(LRG_1t1):c.123A>G")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble == "LRG_1t1:c.123A>G"

def test_double_colon_with_versioned_identifier():
    variant = MockVariant("NM_000001.1:abc.1:c.123A>G")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert len(variant.warnings) == 1

    assert "concatenation" in variant.warnings[0]


def test_ins_length_parentheses_without_colon():
    variant = MockVariant("ins(10)")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any("N[10]" in w for w in variant.warnings)


def test_ins_uncertain_without_colon():
    variant = MockVariant("ins[(10_20)]")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any("syntactically correct" in w for w in variant.warnings)


def test_repeat_expansion_without_colon():
    variant = MockVariant("insA[3]")

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble == "insAAA"
    assert any("may also be written as" in w for w in variant.warnings)


def test_refseq_common_mistakes_object_variant():
    class DummyHGVS:
        ac = "NM_000001.1"

        def __str__(self):
            return "NM_000001.1:g.1A>G"

    v = MockRefVariant(DummyHGVS(), ":g.", "c")
    v.quibble = DummyHGVS()

    assert refseq_common_mistakes(v) is True
    assert any("Did you mean" in w for w in v.warnings)


def test_refseq_common_mistakes_object_nc():
    class DummyHGVS:
        ac = "NC_000001.11"

        def __str__(self):
            return "NC_000001.11:c.1A>G"

    v = MockRefVariant(DummyHGVS(), ":c.")
    v.quibble = DummyHGVS()

    assert refseq_common_mistakes(v) is True

    assert len(v.warnings) == 2

def test_structure_checks_g_invalid_accession():
    variant = MagicMock()
    variant.input_parses.ac = "XM_12345"
    variant.warnings = []

    validator = MagicMock()

    assert structure_checks_g(variant, validator) is True

    assert any(
        "Invalid reference sequence identifier"
        in w
        for w in variant.warnings
    )


def test_structure_checks_g_uncertain_reference_sequence():
    variant = MagicMock()

    variant.input_parses.ac = "NC_000001.11"
    variant.input_parses.posedit.edit.ref = "NN"
    variant.warnings = []

    variant.hn.normalize.return_value = variant.input_parses

    validator = MagicMock()
    validator.vr.validate.return_value = None

    assert structure_checks_g(
        variant,
        validator,
    ) is True

    assert any(
        "UncertainSequenceError"
        in w
        for w in variant.warnings
    )


def test_structure_checks_g_normalize_error():
    variant = MagicMock()

    variant.input_parses.ac = "NC_000001.11"
    variant.input_parses.posedit.edit.ref = "A"
    variant.warnings = []

    variant.hn.normalize.side_effect = (
        vvhgvs.exceptions.HGVSError(
            "normalization failed"
        )
    )

    validator = MagicMock()
    validator.vr.validate.return_value = None

    assert structure_checks_g(
        variant,
        validator,
    ) is True

    assert any(
        "normalization failed"
        in w
        for w in variant.warnings
    )


def test_structure_checks_c_c_to_n_failure():
    variant = MagicMock()

    variant.input_parses = MagicMock()
    variant.input_parses.__str__.return_value = "NM_1:c.*1A>G"

    variant.warnings = []

    variant.evm.c_to_n.side_effect = (
        vvhgvs.exceptions.HGVSError(
            "conversion failed"
        )
    )

    validator = MagicMock()
    validator.vr.validate.return_value = None

    assert structure_checks_c(
        variant,
        validator,
    ) is True

    assert any(
        "conversion failed"
        in w
        for w in variant.warnings
    )


def test_structure_checks_n_data_not_available():
    variant = MagicMock()

    variant.input_parses = MagicMock()
    variant.input_parses.__str__.return_value = "NR_1:n.123A>G"

    variant.warnings = []

    validator = MagicMock()

    validator.vr.validate.side_effect = (
        vvhgvs.exceptions.HGVSDataNotAvailableError(
            "UTA missing"
        )
    )

    assert structure_checks_n(
        variant,
        validator,
    ) is True

    assert any(
        "UTA missing"
        in w
        for w in variant.warnings
    )


def test_structure_checks_n_invalid_interval():
    variant = MagicMock()

    variant.input_parses = MagicMock()
    variant.input_parses.__str__.return_value = "NR_1:n.20_10del"

    variant.input_parses.posedit.pos.start = 20
    variant.input_parses.posedit.pos.end = 10

    variant.warnings = []

    validator = MagicMock()

    validator.vr.validate.side_effect = (
        vvhgvs.exceptions.HGVSInvalidVariantError(
            "base start position must be <= end position"
        )
    )

    assert structure_checks_n(
        variant,
        validator,
    ) is True

    assert any(
        "Interval start position"
        in w
        for w in variant.warnings
    )


def test_structure_checks_c_invalid_interval():
    variant = MagicMock()

    variant.input_parses = MagicMock()
    variant.input_parses.__str__.return_value = "NM_1:c.20_10del"

    variant.input_parses.posedit.pos.start = 20
    variant.input_parses.posedit.pos.end = 10

    variant.warnings = []

    validator = MagicMock()

    validator.vr.validate.side_effect = (
        vvhgvs.exceptions.HGVSInvalidVariantError(
            "base start position must be <= end position"
        )
    )

    assert structure_checks_c(
        variant,
        validator,
    ) is True

    assert any(
        "Interval start position"
        in w
        for w in variant.warnings
    )


def test_structure_checks_g_validator_error():
    variant = MagicMock()

    variant.input_parses.ac = "NC_000001.11"
    variant.warnings = []

    validator = MagicMock()

    validator.vr.validate.side_effect = Exception("validation failed")

    assert structure_checks_g(
        variant,
        validator,
    ) is True

    assert any(
        "validation failed"
        in w
        for w in variant.warnings
    )

def make_input_parses(
    ac="NC_000001.11",
    hgvs_type="g",
    ref="A",
    alt="T",
    edit_type="sub",
    start=100,
    end=100,
    start_offset=0,
    end_offset=0,
):
    start_pos = SimpleNamespace(
        base=start,
        offset=start_offset,
        __str__=lambda self: str(start),
    )

    end_pos = SimpleNamespace(
        base=end,
        offset=end_offset,
        __str__=lambda self: str(end),
    )

    position = SimpleNamespace(
        start=start_pos,
        end=end_pos,
        __str__=lambda self: f"{start}_{end}",
    )

    edit = SimpleNamespace(
        ref=ref,
        alt=alt,
        type=edit_type,
    )

    posedit = SimpleNamespace(
        pos=position,
        edit=edit,
    )

    class DummyHGVS:
        def __init__(self):
            self.ac = ac
            self.type = hgvs_type
            self.posedit = posedit

        def __str__(self):
            return f"{ac}:{hgvs_type}.{start}{ref}>{alt}"

    return DummyHGVS()


def make_variant(ac="NC_000001.11", hgvs_type="g"):
    variant = MagicMock()

    variant.input_parses = make_input_parses(
        ac=ac,
        hgvs_type=hgvs_type,
    )

    variant.warnings = []
    variant.original = str(variant.input_parses)
    variant.hgvs_formatted = variant.input_parses

    variant.hn = MagicMock()
    variant.hn.normalize.side_effect = lambda x: x

    variant.evm = MagicMock()
    variant.no_replace_vm = MagicMock()
    variant.no_norm_evm = MagicMock()

    variant.primary_assembly = "GRCh38"

    return variant


def make_validator():
    validator = MagicMock()

    validator.vr = MagicMock()

    validator.hdp = MagicMock()

    validator.sf = MagicMock()
    validator.sf.fetch_seq.return_value = "A" * 5000

    validator.myevm_t_to_g = MagicMock(
        side_effect=lambda *args, **kwargs: args[0]
    )

    validator.noreplace_myevm_t_to_g = MagicMock(
        side_effect=lambda *args, **kwargs: args[0]
    )

    validator.alt_aln_method = "splign"

    return validator


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
