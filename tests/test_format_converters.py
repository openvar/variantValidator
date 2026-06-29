from unittest.mock import MagicMock, patch
from VariantValidator.modules.format_converters import final_hgvs_convert
from VariantValidator.modules.format_converters import vcf2hgvs_stage1


class MockVariant:
    def __init__(self, quibble):
        self.quibble = quibble
        self.original = quibble
        self.warnings = []
        self.write = True
        self.primary_assembly = "GRCh38"
        self.order = 1


class MockValidator:
    def __init__(self):
        # Used by final_hgvs_convert()
        self.hp = MagicMock()
        self.hp.parse_c_posedit.return_value = "C"
        self.hp.parse_g_posedit.return_value = "G"
        self.hp.parse_m_posedit.return_value = "M"
        self.hp.parse_n_posedit.return_value = "N"
        self.hp.parse_p_posedit.return_value = "P"
        self.hp.parse_r_posedit.return_value = "R"

        # Used by vcf2hgvs_stage1()
        self.batch_list = []


@patch("VariantValidator.modules.format_converters.vvhgvs.sequencevariant.SequenceVariant")
def test_final_hgvs_convert_c(mock_sv):
    validator = MockValidator()
    variant = MockVariant("NM_000001.1:c.123A>G")

    assert final_hgvs_convert(variant, validator) is False

    validator.hp.parse_c_posedit.assert_called_once_with("123A>G")
    mock_sv.assert_called_once()


@patch("VariantValidator.modules.format_converters.vvhgvs.sequencevariant.SequenceVariant")
def test_final_hgvs_convert_g(mock_sv):
    validator = MockValidator()
    variant = MockVariant("NC_000001.11:g.123A>G")

    assert final_hgvs_convert(variant, validator) is False

    validator.hp.parse_g_posedit.assert_called_once_with("123A>G")


@patch("VariantValidator.modules.format_converters.vvhgvs.sequencevariant.SequenceVariant")
def test_final_hgvs_convert_m(mock_sv):
    validator = MockValidator()
    variant = MockVariant("NC_012920.1:m.123A>G")

    assert final_hgvs_convert(variant, validator) is False

    validator.hp.parse_m_posedit.assert_called_once_with("123A>G")


@patch("VariantValidator.modules.format_converters.vvhgvs.sequencevariant.SequenceVariant")
def test_final_hgvs_convert_n(mock_sv):
    validator = MockValidator()
    variant = MockVariant("NR_000001.1:n.123A>G")

    assert final_hgvs_convert(variant, validator) is False

    validator.hp.parse_n_posedit.assert_called_once_with("123A>G")


@patch("VariantValidator.modules.format_converters.vvhgvs.sequencevariant.SequenceVariant")
def test_final_hgvs_convert_p(mock_sv):
    validator = MockValidator()
    variant = MockVariant("NP_000001.1:p.Gly12Val")

    assert final_hgvs_convert(variant, validator) is False

    validator.hp.parse_p_posedit.assert_called_once_with("Gly12Val")


@patch("VariantValidator.modules.format_converters.vvhgvs.sequencevariant.SequenceVariant")
def test_final_hgvs_convert_r(mock_sv):
    validator = MockValidator()
    variant = MockVariant("NM_000001.1:r.123a>g")

    assert final_hgvs_convert(variant, validator) is False

    validator.hp.parse_r_posedit.assert_called_once_with("123a>g")


def test_final_hgvs_convert_r_with_t():
    validator = MockValidator()
    variant = MockVariant("NM_000001.1:r.123T>G")

    assert final_hgvs_convert(variant, validator) is True

    assert any("RNA alphabet" in w for w in variant.warnings)


def test_final_hgvs_convert_invalid_type():
    validator = MockValidator()
    variant = MockVariant("NM_000001.1:x.123A>G")

    assert final_hgvs_convert(variant, validator) is True

    assert any(
        "allowed HGVS type characters" in w
        for w in variant.warnings
    )

def test_vcf_stage1_too_few_fields():
    variant = MockVariant("chr1-123")
    validator = MockValidator()

    assert vcf2hgvs_stage1(variant, validator) is False
    assert variant.quibble == "chr1-123"


def test_vcf_stage1_non_numeric_position():
    variant = MockVariant("chr1-ABC-A-T")
    validator = MockValidator()

    assert vcf2hgvs_stage1(variant, validator) is False


def test_vcf_stage1_dot_ref():
    variant = MockVariant("1-100-.-T")
    validator = MockValidator()

    assert vcf2hgvs_stage1(variant, validator) is False
    assert variant.quibble == "1:100insT"


def test_vcf_stage1_dot_alt():
    variant = MockVariant("1-100-A-.")
    validator = MockValidator()

    assert vcf2hgvs_stage1(variant, validator) is True

    assert variant.write is False
    assert len(validator.batch_list) == 2
    assert "VariantValidator has output both alternatives" in variant.warnings[-1]


def test_vcf_stage1_simple_substitution():
    variant = MockVariant("1-100-A-T")
    validator = MockValidator()

    assert vcf2hgvs_stage1(variant, validator) is False
    assert variant.quibble == "1:100A>T"


def test_vcf_stage1_invalid_bases():
    variant = MockVariant("1-100-XYZ-T")
    validator = MockValidator()

    assert vcf2hgvs_stage1(variant, validator) is False
    assert variant.quibble == "1-100-XYZ-T"


def test_vcf_stage1_assembly_prefix_removed():
    variant = MockVariant("GRCh38-1-100-A-T")
    validator = MockValidator()

    assert vcf2hgvs_stage1(variant, validator) is False
    assert variant.quibble == "1:100A>T"

def test_vcf_stage1_cnv_del():
    variant = MockVariant("1-100-200-del")
    validator = MockValidator()

    assert vcf2hgvs_stage1(variant, validator) is False
    assert variant.quibble == "1:100_200del"
    assert any("CNV identified" in warning for warning in variant.warnings)


def test_vcf_stage1_cnv_inv():
    variant = MockVariant("1-100-200-inv")
    validator = MockValidator()

    assert vcf2hgvs_stage1(variant, validator) is False
    assert variant.quibble == "1:100_200inv"
    assert any("CNV identified" in warning for warning in variant.warnings)


def test_vcf_stage1_multi_alt():
    variant = MockVariant("1-100-A-T,G")
    validator = MockValidator()

    assert vcf2hgvs_stage1(variant, validator) is False
    assert variant.quibble == "1:100A>T,G"


def test_vcf_stage1_hg19_prefix():
    variant = MockVariant("hg19-1-100-A-T")
    validator = MockValidator()

    assert vcf2hgvs_stage1(variant, validator) is False
    assert variant.quibble == "1:100A>T"


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

