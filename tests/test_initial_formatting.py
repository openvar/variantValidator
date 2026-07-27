from unittest.mock import MagicMock

from VariantValidator.modules.initial_formatting import (
    remove_gene_symbol_from_ref,
    initial_user_formattng,
)


def make_variant(original=None):
    variant = MagicMock()
    variant.original = original or "NM_000001.1(BRCA1):c.123A>G"
    variant.warnings = []

    variant.quibble = MagicMock()
    variant.quibble.ac = variant.original.split(":")[0]

    return variant


def make_validator():
    validator = MagicMock()
    validator.db.get_hgnc_symbol.return_value = "BRCA1"
    return validator


def test_remove_gene_symbol_from_ref_object():
    variant = make_variant()
    validator = make_validator()

    remove_gene_symbol_from_ref(variant, validator)

    assert variant.quibble.ac == "NM_000001.1"
    assert len(variant.warnings) == 1
    assert "Removing redundant gene symbol BRCA1" in variant.warnings[0]


def test_remove_gene_symbol_not_hgnc():
    variant = make_variant()
    validator = make_validator()

    validator.db.get_hgnc_symbol.return_value = "none"

    remove_gene_symbol_from_ref(variant, validator)

    assert variant.quibble.ac == "NM_000001.1(BRCA1)"
    assert variant.warnings == []


def test_remove_gene_symbol_mt_gene():
    variant = make_variant("NC_012920.1(MT-ND1):m.123A>G")
    validator = make_validator()

    validator.db.get_hgnc_symbol.return_value = "none"

    remove_gene_symbol_from_ref(variant, validator)

    assert variant.quibble.ac == "NC_012920.1"
    assert len(variant.warnings) == 1


def test_initial_user_formatting_del_warning():
    variant = make_variant("NM_000001.1:c.123_124delAG")
    validator = make_validator()

    initial_user_formattng(variant, validator)

    assert any(
        "Removing redundant reference bases" in w
        for w in variant.warnings
    )


def test_initial_user_formatting_ins_no_warning():
    variant = make_variant("NM_000001.1:c.123_124insAG")
    validator = make_validator()

    initial_user_formattng(variant, validator)

    assert not any(
        "Removing redundant reference bases" in w
        for w in variant.warnings
    )


def test_initial_user_formatting_lowercase_accession():
    variant = make_variant("nm_000001.1:c.123A>G")
    validator = make_validator()

    initial_user_formattng(variant, validator)

    assert any(
        "wrong case" in w
        for w in variant.warnings
    )


def test_initial_user_formatting_chr_accession():
    variant = make_variant("chr1:g.123A>G")
    validator = make_validator()

    initial_user_formattng(variant, validator)

    assert any(
        "no reference sequence ID" in w
        for w in variant.warnings
    )


def test_initial_user_formatting_lrg_uppercase_t():
    variant = make_variant("LRG_1T1:c.123A>G")
    validator = make_validator()

    initial_user_formattng(variant, validator)

    assert variant.quibble.ac == "LRG_1t1"
    assert any(
        "wrong case" in w
        for w in variant.warnings
    )


def test_initial_user_formatting_lowercase_lrg_prefix():
    variant = make_variant("lrg_1t1:c.123A>G")
    validator = make_validator()

    initial_user_formattng(variant, validator)

    assert any(
        "wrong case" in warning.lower()
        for warning in variant.warnings
    )


def test_initial_user_formatting_grch_reference():
    variant = make_variant("GRCh38:g.123A>G")
    validator = make_validator()

    initial_user_formattng(variant, validator)

    assert any(
        "no reference sequence ID" in warning
        for warning in variant.warnings
    )


def test_initial_user_formatting_hg_reference():
    variant = make_variant("HG19:g.123A>G")
    validator = make_validator()

    initial_user_formattng(variant, validator)

    assert any(
        "no reference sequence ID" in warning
        for warning in variant.warnings
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
