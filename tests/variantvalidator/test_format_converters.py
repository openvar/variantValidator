from unittest.mock import MagicMock, patch

import vvhgvs.exceptions

from VariantValidator.modules import format_converters
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
    batch_list = []

    assert vcf2hgvs_stage1(variant, batch_list) is False
    assert variant.quibble == "chr1-123"


def test_vcf_stage1_non_numeric_position():
    variant = MockVariant("chr1-ABC-A-T")
    batch_list = []

    assert vcf2hgvs_stage1(variant, batch_list) is False


def test_vcf_stage1_dot_ref():
    variant = MockVariant("1-100-.-T")
    batch_list = []

    assert vcf2hgvs_stage1(variant, batch_list) is False
    assert variant.quibble == "1:100insT"


def test_vcf_stage1_dot_alt():
    variant = MockVariant("1-100-A-.")
    batch_list = []

    assert vcf2hgvs_stage1(variant, batch_list) is True

    assert variant.write is False
    assert len(batch_list) == 2
    assert "VariantValidator has output both alternatives" in variant.warnings[-1]


def test_vcf_stage1_simple_substitution():
    variant = MockVariant("1-100-A-T")
    batch_list = []

    assert vcf2hgvs_stage1(variant, batch_list) is False
    assert variant.quibble == "1:100A>T"


def test_vcf_stage1_invalid_bases():
    variant = MockVariant("1-100-XYZ-T")
    batch_list = []

    assert vcf2hgvs_stage1(variant, batch_list) is False
    assert variant.quibble == "1-100-XYZ-T"


def test_vcf_stage1_assembly_prefix_removed():
    variant = MockVariant("GRCh38-1-100-A-T")
    batch_list = []

    assert vcf2hgvs_stage1(variant, batch_list) is False
    assert variant.quibble == "1:100A>T"


def test_vcf_stage1_cnv_del():
    variant = MockVariant("1-100-200-del")
    batch_list = []

    assert vcf2hgvs_stage1(variant, batch_list) is False
    assert variant.quibble == "1:100_200del"
    assert any("CNV identified" in warning for warning in variant.warnings)


def test_vcf_stage1_cnv_inv():
    variant = MockVariant("1-100-200-inv")
    batch_list = []

    assert vcf2hgvs_stage1(variant, batch_list) is False
    assert variant.quibble == "1:100_200inv"
    assert any("CNV identified" in warning for warning in variant.warnings)


def test_vcf_stage1_multi_alt():
    variant = MockVariant("1-100-A-T,G")
    batch_list = []

    assert vcf2hgvs_stage1(variant, batch_list) is False
    assert variant.quibble == "1:100A>T,G"


def test_vcf_stage1_hg19_prefix():
    variant = MockVariant("hg19-1-100-A-T")
    batch_list = []

    assert vcf2hgvs_stage1(variant, batch_list) is False
    assert variant.quibble == "1:100A>T"

@patch("VariantValidator.modules.format_converters.Variant")
def test_refseq_catch_ng_creates_transcript_variants(mock_variant):
    from VariantValidator.modules.format_converters import refseq_catch

    validator = MagicMock()

    validator.db.get_gene_symbol_from_refseq_id.return_value = "GENE1"
    validator.db.get_uta_symbol.return_value = "GENE1"

    validator.hdp.get_tx_for_gene.return_value = [
        (None, None, None, "NM_000001.1"),
        (None, None, None, "NM_000002.1"),
        (None, None, None, "ENST00000311111.2"),
    ]

    validator.select_transcripts = "mane"

    variant = MockVariant("NG_000001.1:c.123A>G")

    batch_list = []

    select_transcripts = {
        "NM_000001.1": None,
        "ENST00000311111.2": None,
    }

    assert (
        refseq_catch(
            variant,
            validator,
            select_transcripts,
            batch_list,
        )
        is True
    )

    assert variant.write is False
    assert len(batch_list) == 2

    assert mock_variant.call_count == 2

    first = mock_variant.call_args_list[0].kwargs

    assert first["primary_assembly"] == variant.primary_assembly
    assert first["order"] == variant.order
    assert first["warnings"] == variant.warnings
    assert first["quibble"] == "NG_000001.1(NM_000001.1):c.123A>G"

    second = mock_variant.call_args_list[1].kwargs

    assert second["quibble"] == "NG_000001.1(ENST00000311111.2):c.123A>G"


@patch("VariantValidator.modules.format_converters.Variant")
def test_gene_symbol_catch_creates_variants(mock_variant):
    from VariantValidator.modules.format_converters import gene_symbol_catch

    validator = MagicMock()

    validator.alt_aln_method = "splign"
    validator.select_transcripts = "select"

    validator.db.get_hgnc_symbol.return_value = "BRCA1"
    validator.db.get_uta_symbol.return_value = "BRCA1"

    validator.hdp.get_tx_for_gene.return_value = [
        (None, None, None, "NM_000001.1"),
        (None, None, None, "NM_000002.1"),
        (None, None, None, "ENST00000311111.2"),
    ]

    variant = MockVariant("BRCA1:c.123A>G")
    batch_list = []

    select_transcripts = {
        "NM_000001.1": None,
        "NM_000002.1": None,
    }

    assert gene_symbol_catch(
        variant,
        validator,
        select_transcripts,
        batch_list,
    ) is True

    # Original variant should no longer be written
    assert variant.write is False

    # One replacement per requested transcript
    assert mock_variant.call_count == 2
    assert len(batch_list) == 2

    # First replacement
    kwargs = mock_variant.call_args_list[0].kwargs
    assert kwargs["quibble"] == "NM_000001.1:c.123A>G"
    assert kwargs["primary_assembly"] == variant.primary_assembly
    assert kwargs["order"] == variant.order
    assert kwargs["warnings"] == variant.warnings

    # Second replacement
    kwargs = mock_variant.call_args_list[1].kwargs
    assert kwargs["quibble"] == "NM_000002.1:c.123A>G"
    assert kwargs["primary_assembly"] == variant.primary_assembly
    assert kwargs["order"] == variant.order

    # Warning propagated
    assert any(
        "does not allow the use of a gene symbol" in w
        for w in variant.warnings
    )


def test_gene_symbol_catch_invalid_gene():
    from VariantValidator.modules.format_converters import gene_symbol_catch

    validator = MagicMock()
    validator.db.get_hgnc_symbol.return_value = "none"

    variant = MockVariant("NOTAGENE:c.123A>G")
    batch_list = []

    assert gene_symbol_catch(
        variant,
        validator,
        {},
        batch_list,
    ) is True

    assert variant.write is True
    assert batch_list == []
    assert any("NOTAGENE" in w for w in variant.warnings)


def test_gene_symbol_catch_select_all():
    from VariantValidator.modules.format_converters import gene_symbol_catch

    validator = MagicMock()

    validator.alt_aln_method = "splign"
    validator.select_transcripts = "all"

    validator.db.get_hgnc_symbol.return_value = "BRCA1"
    validator.db.get_uta_symbol.return_value = "BRCA1"

    validator.hdp.get_tx_for_gene.return_value = [
        (None, None, None, "NM_000001.1"),
        (None, None, None, "NM_000002.1"),
    ]

    variant = MockVariant("BRCA1:c.123A>G")
    batch_list = []

    assert gene_symbol_catch(
        variant,
        validator,
        {},
        batch_list,
    ) is True

    assert variant.write is True
    assert batch_list == []
    assert any("select_transcripts=" in w for w in variant.warnings)


@patch("VariantValidator.modules.format_converters.Variant")
def test_gene_symbol_catch_genebuild_uses_enst(mock_variant):
    from VariantValidator.modules.format_converters import gene_symbol_catch

    validator = MagicMock()

    validator.alt_aln_method = "genebuild"
    validator.select_transcripts = "select"

    validator.db.get_hgnc_symbol.return_value = "BRCA1"
    validator.db.get_uta_symbol.return_value = "BRCA1"

    validator.hdp.get_tx_for_gene.return_value = [
        (None, None, None, "NM_000001.1"),
        (None, None, None, "ENST00000311111.2"),
    ]

    variant = MockVariant("BRCA1:c.123A>G")
    batch_list = []

    gene_symbol_catch(
        variant,
        validator,
        {"ENST00000311111.2": None},
        batch_list,
    )

    kwargs = mock_variant.call_args.kwargs

    assert kwargs["quibble"] == "ENST00000311111.2:c.123A>G"


def test_gene_symbol_catch_non_transcript_variant():
    from VariantValidator.modules.format_converters import gene_symbol_catch

    validator = MagicMock()

    variant = MockVariant("BRCA1:g.123A>G")
    batch_list = []

    assert gene_symbol_catch(
        variant,
        validator,
        {},
        batch_list,
    ) is False

    validator.db.get_hgnc_symbol.assert_not_called()
    assert batch_list == []

def test_gene_symbol_catch_parenthesised_gene():
    from VariantValidator.modules.format_converters import gene_symbol_catch

    validator = MagicMock()

    variant = MockVariant("BRCA1(ENST00000357654):c.123A>G")
    batch_list = []

    assert gene_symbol_catch(
        variant,
        validator,
        {},
        batch_list,
    ) is False

    validator.db.get_hgnc_symbol.assert_not_called()


def test_refseq_catch_nc_without_transcript():
    from VariantValidator.modules.format_converters import refseq_catch

    validator = MagicMock()

    variant = MockVariant("NC_000001.11:c.123A>G")
    batch_list = []

    assert refseq_catch(
        variant,
        validator,
        {},
        batch_list,
    ) is True

    assert any(
        "Unable to predict available transcripts"
        in w
        for w in variant.warnings
    )


def test_refseq_catch_ng_unknown_gene():
    from VariantValidator.modules.format_converters import refseq_catch

    validator = MagicMock()

    validator.db.get_gene_symbol_from_refseq_id.return_value = "none"

    variant = MockVariant("NG_000001.1:c.123A>G")
    batch_list = []

    assert refseq_catch(
        variant,
        validator,
        {},
        batch_list,
    ) is True

    assert any(
        "NG_(NM_)"
        in w
        for w in variant.warnings
    )


def test_refseq_catch_invalid_nested_reference():
    from VariantValidator.modules.format_converters import refseq_catch

    validator = MagicMock()

    variant = MockVariant(
        "NG_000001.1(ABC123):c.123A>G"
    )

    batch_list = []

    assert refseq_catch(
        variant,
        validator,
        {},
        batch_list,
    ) is True

    assert any(
        "NG_(NM_):c.PositionVariation"
        in w
        for w in variant.warnings
    )


def test_refseq_catch_multiple_genomic_refs():
    from VariantValidator.modules.format_converters import refseq_catch

    validator = MagicMock()

    variant = MockVariant(
        "NG_000001.1(NC_000001.11(NM_000001.1):c.123A>G"
    )

    batch_list = []

    assert refseq_catch(
        variant,
        validator,
        {},
        batch_list,
    ) is True

    assert any(
        "Multiple genomic reference sequences"
        in w
        for w in variant.warnings
    )


def test_gene_symbol_catch_empty_transcript_list():
    from VariantValidator.modules.format_converters import gene_symbol_catch

    validator = MagicMock()

    validator.alt_aln_method = "splign"
    validator.select_transcripts = "all"

    validator.db.get_hgnc_symbol.return_value = "BRCA1"
    validator.db.get_uta_symbol.return_value = "BRCA1"

    validator.hdp.get_tx_for_gene.return_value = []

    variant = MockVariant("BRCA1:c.123A>G")

    assert gene_symbol_catch(
        variant,
        validator,
        {},
        [],
    ) is True

    assert any(
        "select_transcripts="
        in w
        for w in variant.warnings
    )

def test_gene_symbol_catch_raw_select_transcripts():
    from VariantValidator.modules.format_converters import gene_symbol_catch

    validator = MagicMock()
    validator.alt_aln_method = "splign"
    validator.select_transcripts = "raw"

    validator.db.get_hgnc_symbol.return_value = "BRCA1"
    validator.db.get_uta_symbol.return_value = "BRCA1"
    validator.hdp.get_tx_for_gene.return_value = [
        (None, None, None, "NM_000001.1"),
    ]

    variant = MockVariant("BRCA1:c.123A>G")

    assert gene_symbol_catch(variant, validator, {}, []) is True
    assert variant.write is True


@patch("VariantValidator.modules.format_converters.Variant")
def test_gene_symbol_catch_duplicate_transcripts(mock_variant):
    from VariantValidator.modules.format_converters import gene_symbol_catch

    validator = MagicMock()
    validator.alt_aln_method = "splign"
    validator.select_transcripts = "select"

    validator.db.get_hgnc_symbol.return_value = "BRCA1"
    validator.db.get_uta_symbol.return_value = "BRCA1"

    validator.hdp.get_tx_for_gene.return_value = [
        (None, None, None, "NM_000001.1"),
        (None, None, None, "NM_000001.1"),
    ]

    variant = MockVariant("BRCA1:c.123A>G")
    batch = []

    gene_symbol_catch(
        variant,
        validator,
        {"NM_000001.1": None},
        batch,
    )

    assert len(batch) == 1


def test_refseq_catch_raw_mode():
    from VariantValidator.modules.format_converters import refseq_catch

    validator = MagicMock()

    validator.select_transcripts = "raw"

    validator.db.get_gene_symbol_from_refseq_id.return_value = "GENE"
    validator.db.get_uta_symbol.return_value = "GENE"

    validator.hdp.get_tx_for_gene.return_value = [
        (None, None, None, "NM_000001.1"),
    ]

    variant = MockVariant("NG_000001.1:c.123A>G")

    assert refseq_catch(
        variant,
        validator,
        {},
        [],
    ) is True

    assert variant.write is True


def test_refseq_catch_nested_valid_reference():
    from VariantValidator.modules.format_converters import refseq_catch

    validator = MagicMock()

    variant = MockVariant(
        "NG_000001.1(NM_000001.1):c.123A>G"
    )

    assert refseq_catch(
        variant,
        validator,
        {},
        [],
    ) is False


def test_gene_symbol_catch_non_coding_type():
    from VariantValidator.modules.format_converters import gene_symbol_catch

    validator = MagicMock()

    variant = MockVariant("BRCA1:p.Gly12Val")

    assert gene_symbol_catch(
        variant,
        validator,
        {},
        [],
    ) is False

    validator.db.get_hgnc_symbol.assert_not_called()


def test_refseq_catch_non_coding_type():
    from VariantValidator.modules.format_converters import refseq_catch

    validator = MagicMock()

    variant = MockVariant("NG_000001.1:g.123A>T")

    assert refseq_catch(
        variant,
        validator,
        {},
        [],
    ) is False



# ---------------------------------------------------------------------------
# initial_format_conversions
# ---------------------------------------------------------------------------

def test_initial_format_conversions_hgvs_parse_error():
    variant = MockVariant("NM_000546.6:c.bad")
    validator = MockValidator()

    with patch.object(
        format_converters,
        "vcf2hgvs_stage1",
        return_value=False,
    ), patch.object(
        format_converters,
        "vcf2hgvs_stage2",
        return_value=False,
    ), patch.object(
        format_converters,
        "gene_symbol_catch",
        return_value=False,
    ), patch.object(
        format_converters,
        "refseq_catch",
        return_value=False,
    ), patch.object(
        format_converters,
        "vcf2hgvs_stage4",
        return_value=False,
    ), patch.object(
        format_converters.use_checking,
        "pre_parsing_global_common_mistakes",
        return_value=False,
    ), patch.object(
        format_converters.methyl_syntax,
        "methyl_syntax",
    ), patch.object(
        format_converters,
        "uncertain_pos",
        return_value=False,
    ), patch.object(
        format_converters.use_checking,
        "refseq_common_mistakes",
        return_value=False,
    ), patch.object(
        format_converters,
        "intronic_converter",
    ), patch.object(
        format_converters,
        "final_hgvs_convert",
        side_effect=vvhgvs.exceptions.HGVSParseError("bad HGVS"),
    ):
        result = format_converters.initial_format_conversions(
            variant,
            validator,
            {},
            [],
        )

    assert result is True
    assert variant.warnings == [
        "HgvsSyntaxError: bad HGVS"
    ]


def test_initial_format_conversions_generic_hgvs_error():
    variant = MockVariant("NM_000546.6:c.bad")
    validator = MockValidator()

    with patch.object(
        format_converters,
        "vcf2hgvs_stage1",
        return_value=False,
    ), patch.object(
        format_converters,
        "vcf2hgvs_stage2",
        return_value=False,
    ), patch.object(
        format_converters,
        "gene_symbol_catch",
        return_value=False,
    ), patch.object(
        format_converters,
        "refseq_catch",
        return_value=False,
    ), patch.object(
        format_converters,
        "vcf2hgvs_stage4",
        return_value=False,
    ), patch.object(
        format_converters.use_checking,
        "pre_parsing_global_common_mistakes",
        return_value=False,
    ), patch.object(
        format_converters.methyl_syntax,
        "methyl_syntax",
    ), patch.object(
        format_converters,
        "uncertain_pos",
        return_value=False,
    ), patch.object(
        format_converters.use_checking,
        "refseq_common_mistakes",
        return_value=False,
    ), patch.object(
        format_converters,
        "intronic_converter",
    ), patch.object(
        format_converters,
        "final_hgvs_convert",
        side_effect=vvhgvs.exceptions.HGVSError("boom"),
    ):
        result = format_converters.initial_format_conversions(
            variant,
            validator,
            {},
            [],
        )

    assert result is True
    assert variant.warnings == [
        "HgvsParserError: Unknown error during reading of variant "
        "NM_000546.6:c.bad"
    ]

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
