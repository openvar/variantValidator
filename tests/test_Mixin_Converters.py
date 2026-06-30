from unittest.mock import MagicMock

from VariantValidator.modules.vvMixinConverters import Mixin


class DummyMixin(Mixin):
    pass


def make_mixin():
    mixin = DummyMixin.__new__(DummyMixin)

    mixin.sf = MagicMock()
    mixin.hp = MagicMock()
    mixin.vm = MagicMock()

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
