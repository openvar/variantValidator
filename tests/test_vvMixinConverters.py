from unittest.mock import MagicMock, patch
import pytest
import json

from VariantValidator.modules.vvMixinConverters import Mixin
from vvhgvs.exceptions import HGVSDataNotAvailableError, HGVSError

from unittest import TestCase

from VariantValidator.validator import Validator


class TestVVMixinConvertersFunctional(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()

    def validate(self, variant, assembly="GRCh38", transcripts="all"):
        return self.vv.validate(
            variant,
            assembly,
            transcripts,
        ).format_as_dict(test=True)

    def test_allele_standard(self):
        result = self.validate("NM_000088.3:c.[600C>A];[620A>T]")

        self.assertEqual(result["flag"], "gene_variant")
        self.assertIn("NM_000088.3:c.600C>A", result)
        self.assertIn("NM_000088.3:c.620A>T", result)

    def test_allele_phased_unphased(self):
        result = self.validate("NM_000088.3:c.[600C>A];[620A>T](;)630G>A")

        self.assertEqual(result["flag"], "gene_variant")
        self.assertIn("NM_000088.3:c.600C>A", result)
        self.assertIn("NM_000088.3:c.620A>T", result)
        self.assertIn("NM_000088.3:c.630G>A", result)

    def test_allele_semicolon_parser(self):
        result = self.validate("NM_000088.3:c.600C>A(;)620A>T")

        self.assertEqual(result["flag"], "gene_variant")
        self.assertIn("NM_000088.3:c.600C>A", result)
        self.assertIn("NM_000088.3:c.620A>T", result)

    def test_allele_repeat_single_base_intronic(self):
        result = self.validate("NM_000088.3:c.[600C>A;589-1G[3]]")

        assert result["validation_warning_1"]["validation_warnings"] == ['AlleleSyntaxError: Submitted variants are out of order or their ranges overlap']

    def test_allele_repeat_single_base_intronic_nc_nm_format(self):
        result = self.validate('["NC_000017.11(NM_000088.3):r.[589G>T;600C>A]"]')

        assert result["validation_warning_1"]["validation_warnings"] == [
            "UnsupportedFormatError: RNA allele syntax variants are not currently supported. Please submit individually for additional guidance"
        ]

    def test_allele_repeat_rna(self):
        result = self.validate('["NC_000017.11(NM_000088.3):c.[600C>A;589-1G[3]]"]')

        assert result["validation_warning_1"]["validation_warnings"] == ['AlleleSyntaxError: Submitted variants are out of order or their ranges overlap']

    def test_allele_uncertain_parenthesised(self):
        result = self.validate("NM_000088.3:c.600[C>A];[(C>A)]")

        self.assertEqual(result["flag"], "gene_variant")
        self.assertIn("NM_000088.3:c.600C>A", result)
        self.assertEqual(len([k for k in result if k.startswith("NM_")]), 1)

    def test_allele_unknown(self):
        result = self.validate("NM_000088.3:c.[600C>A];[?]")

        self.assertEqual(result["flag"], "gene_variant")
        self.assertIn("NM_000088.3:c.600C>A", result)
        self.assertEqual(len([k for k in result if k.startswith("NM_")]), 1)

    def test_allele_null(self):
        result = self.validate("NM_000088.3:c.[600C>A];[0]")

        self.assertEqual(result["flag"], "gene_variant")
        self.assertIn("NM_000088.3:c.600C>A", result)
        self.assertEqual(len([k for k in result if k.startswith("NM_")]), 1)

    def test_allele_multiple_variants(self):
        result = self.validate("NM_000088.3:c.[600C>A;620A>T];[630G>A]")

        self.assertEqual(result["flag"], "gene_variant")
        self.assertIn("NM_000088.3:c.600C>A", result)
        self.assertIn("NM_000088.3:c.620A>T", result)
        self.assertIn("NM_000088.3:c.630G>A", result)

    def test_allele_repeat_only(self):
        result = self.validate("NM_000088.3:c.[589-1G[3]]")

        assert result["validation_warning_1"]["validation_warnings"] == ['RepeatSyntaxError: Ensure that the repeated sequence is included between the variant position and the number of repeat units, e.g. g.1_3ACT[20]']

    def test_allele_repeat_single_base_intronic_ordered(self):
        result = self.validate("NM_000088.3:c.[589-1G[3];600C>A]")

        self.assertEqual(result["flag"], "gene_variant")
        self.assertIn("NM_000088.3:c.589-1_590=", result)
        self.assertIn("NM_000088.3:c.600C>A", result)

    def test_noncoding_allele(self):
        variant = "NR_110010.2:n.[138del;150_150+1G[2]]"
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert "NR_110010.2:n.150_150+1=" in results.keys()
        assert "NR_110010.2:n.138del" in results.keys()


    def test_noncoding_allele_compound_accession(self):
        variant = "NC_000017.10(NR_110010.2):n.[138del;150_150+1G[2]]"
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert "NR_110010.2:n.150_150+1=" in results.keys()
        assert "NR_110010.2:n.138del" in results.keys()

    def test_mito_allele(self):
        variant = "NC_012920.1:m.[4del;7dup]"
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "mitochondrial_variant_1" in results.keys()
        assert "mitochondrial_variant_2" in results.keys()


    def test_lrg_transcript_allele(self):
        result = self.validate("LRG_1t1:c.[600C>A];[620A>T]")

        self.assertEqual(result["flag"], "gene_variant")
        self.assertIn("NM_000088.3:c.600C>A", result)
        self.assertIn("NM_000088.3:c.620A>T", result)


    def test_lrg_transcript_repeat(self):
        result = self.validate("LRG_1t1:c.[589-1G[3];600C>A]")

        self.assertEqual(result["flag"], "gene_variant")
        self.assertIn("NM_000088.3:c.589-1_590=", result)
        self.assertIn("NM_000088.3:c.600C>A", result)

    def test_lrg_gene_allele(self):
        result = self.validate("LRG1:g.[50198003C>A];[50198023A>T]")

        assert result["validation_warning_1"]["validation_warnings"] == ['VariantMappingWarning: LRG1 updated to LRG_1: LRG_1:g.[50198003C>A];[50198023A>T] automapped to equivalent RefSeq record NG_007400.1:g.[50198003C>A];[50198023A>T]', 'ReferenceMismatchError: NG_007400.1:g.50198003C>A: Variant reference (C) does not agree with reference sequence ()']


class DummyMixin(Mixin):
    pass


def make_mixin():
    mixin = DummyMixin.__new__(DummyMixin)
    mixin.vm = MagicMock()
    mixin.sf = MagicMock()
    mixin.hdp = MagicMock()
    mixin.alt_aln_method = "splign"
    mixin.hp = MagicMock()
    return mixin


def test_expand_ref_small_region():
    mixin = make_mixin()

    mixin.sf.fetch_seq.return_value = "ACGT"

    pre, post = mixin._expand_ref("NC_1", 100, 104)

    assert pre == "A"
    assert post == "T"

    mixin.sf.fetch_seq.assert_called_once_with(
        "NC_1",
        100,
        104,
    )


def test_expand_ref_large_region():
    mixin = make_mixin()

    mixin.sf.fetch_seq.side_effect = [
        "A",
        "T",
    ]

    pre, post = mixin._expand_ref("NC_1", 100, 1200)

    assert pre == "A"
    assert post == "T"

    assert mixin.sf.fetch_seq.call_count == 2


def test_coding_parses_string_variant():
    mixin = make_mixin()

    parsed = MagicMock()
    parsed.type = "c"

    mixin.hp.parse_hgvs_variant.return_value = parsed

    result = mixin.coding("NM_000001.1:c.123A>G")

    assert result is parsed

    mixin.hp.parse_hgvs_variant.assert_called_once()


def test_coding_converts_n_to_c():
    mixin = make_mixin()

    variant = MagicMock()
    variant.type = "n"

    converted = MagicMock()

    mixin.vm.n_to_c.return_value = converted

    result = mixin.coding(variant)

    assert result is converted

    mixin.vm.n_to_c.assert_called_once_with(variant)


def test_coding_returns_c_unchanged():
    mixin = make_mixin()

    variant = MagicMock()
    variant.type = "c"

    result = mixin.coding(variant)

    assert result is variant


def test_genomic_returns_existing_g_variant():
    mixin = make_mixin()

    variant = MagicMock()
    variant.type = "g"

    vv_variant = MagicMock()

    result = mixin.genomic(
        variant,
        MagicMock(),
        "GRCh38",
        vv_variant,
    )

    assert result is variant


def test_genomic_parses_g_string():
    mixin = make_mixin()

    parsed = MagicMock()

    mixin.hp.parse_hgvs_variant.return_value = parsed

    result = mixin.genomic(
        "NC_000001.11:g.123A>G",
        MagicMock(),
        "GRCh38",
        MagicMock(hn=MagicMock()),
    )

    assert result is parsed

    mixin.hp.parse_hgvs_variant.assert_called_once()


def test_myevm_g_to_t_calls_mapper():
    mixin = make_mixin()

    evm = MagicMock()
    expected = MagicMock()

    evm.g_to_t.return_value = expected

    result = mixin.myevm_g_to_t(
        evm,
        "genomic_variant",
        "NM_000001.1",
    )

    assert result is expected

    evm.g_to_t.assert_called_once_with(
        "genomic_variant",
        "NM_000001.1",
    )

def test_coding_parses_n_string():
    mixin = make_mixin()

    parsed = MagicMock()
    parsed.type = "n"

    converted = MagicMock()

    mixin.hp.parse_hgvs_variant.return_value = parsed
    mixin.vm.n_to_c.return_value = converted

    result = mixin.coding("NR_000001.1:n.123A>G")

    assert result is converted

    mixin.hp.parse_hgvs_variant.assert_called_once()
    mixin.vm.n_to_c.assert_called_once_with(parsed)


def test_coding_invalid_type_returns_none():
    mixin = make_mixin()

    variant = MagicMock()
    variant.type = "g"

    assert mixin.coding(variant) is None


def test_genomic_converts_c_variant():
    mixin = make_mixin()

    variant = MagicMock()
    variant.type = "c"

    validator = MagicMock()

    expected = MagicMock()

    mixin.myevm_t_to_g = MagicMock(return_value=expected)

    result = mixin.genomic(
        variant,
        validator,
        "GRCh38",
        MagicMock()
    )

    assert result is expected

    mixin.myevm_t_to_g.assert_called_once()


def test_genomic_converts_n_variant():
    mixin = make_mixin()

    variant = MagicMock()
    variant.type = "n"

    validator = MagicMock()

    expected = MagicMock()

    mixin.myevm_t_to_g = MagicMock(return_value=expected)

    result = mixin.genomic(
        variant,
        validator,
        "GRCh38",
        MagicMock()
    )

    assert result is expected


def test_genomic_unknown_type_returns_none():
    mixin = make_mixin()

    variant = MagicMock()
    variant.type = "p"

    assert (
        mixin.genomic(
            variant,
            MagicMock(),
            "GRCh38",
            MagicMock()
        )
        is None
    )


def test_expand_ref_single_base():
    mixin = make_mixin()

    mixin.sf.fetch_seq.return_value = "G"

    left, right = mixin._expand_ref(
        "NC_1",
        100,
        101,
    )

    assert left == "G"
    assert right == "G"


def test_myevm_g_to_t_passes_parameters():
    mixin = make_mixin()

    evm = MagicMock()

    expected = MagicMock()

    evm.g_to_t.return_value = expected

    result = mixin.myevm_g_to_t(
        evm,
        "genomic",
        "NM_000001.1"
    )

    assert result is expected

    evm.g_to_t.assert_called_once_with(
        "genomic",
        "NM_000001.1"
    )

def test_coding_parses_c_string():
    mixin = make_mixin()

    parsed = MagicMock()
    parsed.type = "c"

    mixin.hp.parse_hgvs_variant.return_value = parsed

    assert mixin.coding("NM_1:c.1A>G") is parsed


def test_coding_parses_n_string_and_converts():
    mixin = make_mixin()

    parsed = MagicMock()
    parsed.type = "n"

    converted = MagicMock()

    mixin.hp.parse_hgvs_variant.return_value = parsed
    mixin.vm.n_to_c.return_value = converted

    assert mixin.coding("NR_1:n.1A>G") is converted


def test_coding_returns_existing_c_variant():
    mixin = make_mixin()

    variant = MagicMock()
    variant.type = "c"

    assert mixin.coding(variant) is variant


def test_coding_returns_none_for_g_variant():
    mixin = make_mixin()

    variant = MagicMock()
    variant.type = "g"

    assert mixin.coding(variant) is None


def test_expand_ref_exactly_1000_bases():
    mixin = make_mixin()

    mixin.sf.fetch_seq.return_value = "A" * 1000

    left, right = mixin._expand_ref("NC_1", 1, 1001)

    assert left == "A"
    assert right == "A"


def test_expand_ref_over_1000_bases():
    mixin = make_mixin()

    mixin.sf.fetch_seq.side_effect = [
        "A",
        "T",
    ]

    left, right = mixin._expand_ref("NC_1", 1, 1500)

    assert left == "A"
    assert right == "T"

    assert mixin.sf.fetch_seq.call_count == 2


def test_genomic_returns_existing_g():
    mixin = make_mixin()

    variant = MagicMock()
    variant.type = "g"

    assert (
        mixin.genomic(
            variant,
            MagicMock(),
            "GRCh38",
            MagicMock(),
        )
        is variant
    )


def test_genomic_converts_c():
    mixin = make_mixin()

    variant = MagicMock()
    variant.type = "c"

    expected = MagicMock()

    mixin.myevm_t_to_g = MagicMock(return_value=expected)

    assert (
        mixin.genomic(
            variant,
            MagicMock(),
            "GRCh38",
            MagicMock(),
        )
        is expected
    )


def test_genomic_converts_n():
    mixin = make_mixin()

    variant = MagicMock()
    variant.type = "n"

    expected = MagicMock()

    mixin.myevm_t_to_g = MagicMock(return_value=expected)

    assert (
        mixin.genomic(
            variant,
            MagicMock(),
            "GRCh38",
            MagicMock(),
        )
        is expected
    )


def test_genomic_returns_none_for_protein():
    mixin = make_mixin()

    variant = MagicMock()
    variant.type = "p"

    assert (
        mixin.genomic(
            variant,
            MagicMock(),
            "GRCh38",
            MagicMock(),
        )
        is None
    )


def test_myevm_g_to_t():
    mixin = make_mixin()

    evm = MagicMock()

    expected = MagicMock()

    evm.g_to_t.return_value = expected

    assert (
        mixin.myevm_g_to_t(
            evm,
            "g_variant",
            "NM_1",
        )
        is expected
    )

    evm.g_to_t.assert_called_once_with(
        "g_variant",
        "NM_1",
    )

def test_noreplace_t_to_g_success():
    mixin = make_mixin()

    hgvs_c = MagicMock()

    genomic = MagicMock()

    variant = MagicMock()
    variant.evm.t_to_g.return_value = genomic

    result = mixin.noreplace_myevm_t_to_g(hgvs_c, variant)

    assert result is genomic

    variant.evm.t_to_g.assert_called_once_with(hgvs_c)


def test_noreplace_t_to_g_no_mapping_options():

    mixin = make_mixin()

    hgvs_c = MagicMock()

    variant = MagicMock()

    variant.evm.t_to_g.side_effect = HGVSError("boom")
    variant.map_dat.mapping_options.return_value = []

    with pytest.raises(HGVSDataNotAvailableError):
        mixin.noreplace_myevm_t_to_g(hgvs_c, variant)


@patch("VariantValidator.modules.vvMixinConverters.seq_data.supported_for_mapping")
def test_noreplace_search_prefers_supported(mock_supported):

    mock_supported.return_value = "1"

    mixin = make_mixin()

    hgvs_c = MagicMock()

    genomic = MagicMock()

    variant = MagicMock()

    variant.evm.t_to_g.side_effect = HGVSError("boom")

    variant.map_dat.mapping_options.return_value = [
        ("tx", "NC_000001.11", "splign", None, 1),
    ]

    mixin.vm.t_to_g.return_value = genomic

    variant.hn.normalize.side_effect = [Exception(), None]

    result = mixin.noreplace_myevm_t_to_g(hgvs_c, variant)

    assert result is genomic


def test_myevm_g_to_t_returns_mapper_result():

    mixin = make_mixin()

    mapper = MagicMock()

    expected = MagicMock()

    mapper.g_to_t.return_value = expected

    assert mixin.myevm_g_to_t(
        mapper,
        "gvar",
        "NM_1"
    ) is expected


def test_expand_ref_boundary_1001():

    mixin = make_mixin()

    mixin.sf.fetch_seq.side_effect = [
        "A",
        "C",
    ]

    left, right = mixin._expand_ref(
        "NC_1",
        1,
        1002,
    )

    assert left == "A"
    assert right == "C"

    assert mixin.sf.fetch_seq.call_count == 2


def test_expand_ref_boundary_1000():

    mixin = make_mixin()

    mixin.sf.fetch_seq.return_value = "A" * 1000

    left, right = mixin._expand_ref(
        "NC_1",
        1,
        1001,
    )

    assert left == "A"
    assert right == "A"

    mixin.sf.fetch_seq.assert_called_once()

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
