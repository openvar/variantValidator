import pytest
from unittest.mock import MagicMock, patch
from types import SimpleNamespace

import vvhgvs
from vvhgvs.enums import Datum

from VariantValidator.modules.use_checking import (
    structure_checks,
    structure_checks_g,
    structure_checks_c,
    structure_checks_n,
    pre_parsing_global_common_mistakes,
    InvalidVariantError,
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
    start_datum=Datum.SEQ_START,
    end_datum=Datum.SEQ_START,
):
    start_pos = SimpleNamespace(
        base=start,
        offset=start_offset,
        datum=start_datum,
    )

    end_pos = SimpleNamespace(
        base=end,
        offset=end_offset,
        datum=end_datum,
    )

    class DummyPosition:
        def __init__(self):
            self.start = start_pos
            self.end = end_pos

        def __str__(self):
            if start == end and start_offset == end_offset:
                return str(start)
            return f"{start}_{end}"

    position = DummyPosition()

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
            return f"{self.ac}:{self.type}.{self.posedit.pos}"

    return DummyHGVS()


def make_variant(
    ac="NC_000001.11",
    hgvs_type="g",
    **kwargs,
):
    variant = MagicMock()

    variant.input_parses = make_input_parses(
        ac=ac,
        hgvs_type=hgvs_type,
        **kwargs,
    )

    variant.quibble = variant.input_parses
    variant.warnings = []
    variant.original = str(variant.input_parses)
    variant.hgvs_formatted = variant.input_parses

    variant.hn = MagicMock()
    variant.hn.normalize.side_effect = lambda x: x

    variant.evm = MagicMock()
    variant.no_norm_evm = MagicMock()
    variant.no_replace_vm = MagicMock()

    variant.primary_assembly = "GRCh38"
    variant.genomic_context_ac = None

    return variant


def make_validator():
    validator = MagicMock()

    validator.vr = MagicMock()
    validator.hdp = MagicMock()
    validator.vm = MagicMock()

    validator.sf = MagicMock()
    validator.sf.fetch_seq.return_value = "A" * 5000

    validator.myevm_t_to_g = MagicMock(
        side_effect=lambda *args, **kwargs: args[0]
    )

    validator.noreplace_myevm_t_to_g = MagicMock(
        side_effect=lambda *args, **kwargs: args[0]
    )

    validator.vm.t_to_g.side_effect = lambda *args, **kwargs: args[0]
    validator.vm.g_to_t.side_effect = lambda *args, **kwargs: args[0]

    validator.alt_aln_method = "splign"

    return validator


# ---------------------------------------------------------------------------
# structure_checks()
# ---------------------------------------------------------------------------

def test_structure_checks_dispatch_g():
    variant = make_variant()
    variant.quibble = "NC_000001.11:g.100A>G"

    validator = make_validator()

    parsed = make_input_parses(
        ac="NC_000001.11",
        hgvs_type="g",
    )

    validator.hp.parse_hgvs_variant.return_value = parsed
    validator.db.get_gene_symbol_from_transcript_id.return_value = "GENE"

    with patch(
        "VariantValidator.modules.use_checking.structure_checks_g",
        return_value=False,
    ) as mock_g:
        assert structure_checks(variant, validator) is False

    mock_g.assert_called_once_with(variant, validator)


def test_structure_checks_dispatch_c():
    variant = make_variant()
    variant.quibble = "NM_000001.1:c.1A>G"

    validator = make_validator()

    parsed = make_input_parses(
        ac="NM_000001.1",
        hgvs_type="c",
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )

    validator.hp.parse_hgvs_variant.return_value = parsed
    validator.db.get_gene_symbol_from_transcript_id.return_value = "GENE"

    with patch(
        "VariantValidator.modules.use_checking.structure_checks_c",
        return_value=False,
    ) as mock_c:
        assert structure_checks(variant, validator) is False

    mock_c.assert_called_once_with(variant, validator)


def test_structure_checks_dispatch_n():
    variant = make_variant()
    variant.quibble = "NR_000001.1:n.1A>G"

    validator = make_validator()

    parsed = make_input_parses(
        ac="NR_000001.1",
        hgvs_type="n",
    )

    validator.hp.parse_hgvs_variant.return_value = parsed
    validator.db.get_gene_symbol_from_transcript_id.return_value = "GENE"

    with patch(
        "VariantValidator.modules.use_checking.structure_checks_n",
        return_value=False,
    ) as mock_n:
        assert structure_checks(variant, validator) is False

    mock_n.assert_called_once_with(variant, validator)


def test_structure_checks_gene_symbol_none_becomes_empty():
    variant = make_variant()
    variant.quibble = "NM_000001.1:c.1A>G"

    validator = make_validator()

    parsed = make_input_parses(
        ac="NM_000001.1",
        hgvs_type="c",
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )

    validator.hp.parse_hgvs_variant.return_value = parsed
    validator.db.get_gene_symbol_from_transcript_id.return_value = "none"

    with patch(
        "VariantValidator.modules.use_checking.structure_checks_c",
        return_value=False,
    ):
        assert structure_checks(variant, validator) is False

    assert variant.gene_symbol == ""


def test_structure_checks_unknown_type():
    variant = make_variant()
    variant.quibble = "NM_000001.1:x.1"

    validator = make_validator()

    parsed = make_input_parses(
        ac="NM_000001.1",
        hgvs_type="x",
    )

    validator.hp.parse_hgvs_variant.return_value = parsed
    validator.db.get_gene_symbol_from_transcript_id.return_value = "GENE"

    assert structure_checks(variant, validator) is False


def test_structure_checks_uses_existing_hgvs_object():
    variant = make_variant(
        ac="NC_000001.11",
        hgvs_type="g",
    )
    variant.quibble = variant.input_parses

    validator = make_validator()
    validator.db.get_gene_symbol_from_transcript_id.return_value = "GENE"

    with patch(
        "VariantValidator.modules.use_checking.structure_checks_g",
        return_value=False,
    ):
        assert structure_checks(variant, validator) is False

    validator.hp.parse_hgvs_variant.assert_not_called()


# ---------------------------------------------------------------------------
# structure_checks_g()
# ---------------------------------------------------------------------------

def test_structure_checks_g_invalid_accession():
    variant = make_variant(
        ac="XM_12345",
        hgvs_type="g",
    )
    validator = make_validator()

    assert structure_checks_g(variant, validator) is True

    assert any(
        "Invalid reference sequence identifier" in warning
        for warning in variant.warnings
    )


def test_structure_checks_g_validator_error():
    variant = make_variant()
    validator = make_validator()

    validator.vr.validate.side_effect = Exception(
        "validation failed"
    )

    assert structure_checks_g(variant, validator) is True
    assert "validation failed" in variant.warnings


def test_structure_checks_g_generic_validation_error():
    variant = make_variant()
    validator = make_validator()

    validator.vr.validate.side_effect = Exception(
        "Some random validation error"
    )

    assert structure_checks_g(variant, validator) is True
    assert "Some random validation error" in variant.warnings


def test_structure_checks_g_normalizer_exception():
    variant = make_variant()
    validator = make_validator()

    validator.vr.validate.return_value = None

    variant.hn.normalize.side_effect = (
        vvhgvs.exceptions.HGVSError(
            "normalizer failed"
        )
    )

    assert structure_checks_g(variant, validator) is True
    assert "normalizer failed" in variant.warnings


def test_structure_checks_g_uncertain_sequence():
    variant = make_variant(ref="N")
    validator = make_validator()

    validator.vr.validate.return_value = None

    assert structure_checks_g(variant, validator) is True

    assert any(
        "UncertainSequenceError" in warning
        for warning in variant.warnings
    )


def test_structure_checks_g_no_reference_sequence():
    variant = make_variant(ref=None)
    validator = make_validator()

    validator.vr.validate.return_value = None

    assert structure_checks_g(variant, validator) is False
    assert variant.warnings == []


def test_structure_checks_g_insertion_length_message():
    variant = make_variant(
        start=123,
        end=123,
        alt="AA",
        edit_type="ins",
    )
    validator = make_validator()

    validator.vr.validate.side_effect = Exception(
        "insertion length must be 1"
    )

    assert structure_checks_g(variant, validator) is True

    assert any(
        "Insertion length must be 1 e.g." in warning
        for warning in variant.warnings
    )


def test_structure_checks_g_uncertain_insertion_returns_true():
    variant = make_variant(
        edit_type="ins",
    )
    validator = make_validator()

    variant.input_parses.posedit.pos = MagicMock()
    variant.input_parses.posedit.pos.__str__.return_value = "(123_124)"

    validator.vr.validate.side_effect = Exception(
        "insertion length must be 1"
    )

    assert structure_checks_g(variant, validator) is True


def test_structure_checks_g_uncertain_deletion_length():
    variant = make_variant(
        edit_type="del",
    )
    validator = make_validator()

    variant.input_parses.posedit.pos = MagicMock()
    variant.input_parses.posedit.pos.__str__.return_value = "(10_20)"

    validator.vr.validate.side_effect = Exception(
        "Length implied by coordinates must equal sequence deletion length"
    )

    assert structure_checks_g(variant, validator) is True


def test_structure_checks_g_circular_origin():
    variant = make_variant(
        ac="NC_012920.1",
    )
    validator = make_validator()

    variant.hgvs_formatted = variant.input_parses
    variant.original = str(variant.input_parses)

    validator.vr.validate.side_effect = Exception(
        "base start position must be <= end position"
    )

    assert structure_checks_g(variant, validator) is True

    assert any(
        "cannot normalize variants spanning the origin"
        in warning
        for warning in variant.warnings
    )


# ---------------------------------------------------------------------------
# structure_checks_c()
# ---------------------------------------------------------------------------

def test_structure_checks_c_unsupported_operation():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()

    validator.vr.validate.side_effect = (
        vvhgvs.exceptions.HGVSUnsupportedOperationError(
            "unsupported"
        )
    )

    assert structure_checks_c(variant, validator) is False
    assert variant.warnings == []


def test_structure_checks_c_bounds_error():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()

    validator.vr.validate.side_effect = (
        vvhgvs.exceptions.HGVSError("bounds")
    )

    assert structure_checks_c(variant, validator) is True
    assert any("bounds" in w for w in variant.warnings)


def test_structure_checks_c_invalid_interval():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start=20,
        end=10,
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()

    validator.vr.validate.side_effect = (
        vvhgvs.exceptions.HGVSInvalidVariantError(
            "base start position must be <= end position"
        )
    )

    assert structure_checks_c(variant, validator) is True

    assert any(
        "Interval start position" in warning
        for warning in variant.warnings
    )


def test_structure_checks_c_value_error_interval():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()

    validator.vr.validate.return_value = None

    validator.noreplace_myevm_t_to_g.side_effect = ValueError(
        "start > end"
    )

    assert structure_checks_c(variant, validator) is True

    assert any(
        "Interval start position" in warning
        for warning in variant.warnings
    )


def test_structure_checks_c_data_not_available():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()

    validator.vr.validate.return_value = None

    validator.noreplace_myevm_t_to_g.side_effect = (
        vvhgvs.exceptions.HGVSDataNotAvailableError(
            "UTA missing"
        )
    )

    validator.hdp.get_tx_identity_info.return_value = [
        None,
        None,
        None,
        1,
        1000,
    ]

    assert structure_checks_c(variant, validator) is True

    assert any(
        "Universal Transcript Archive" in warning
        for warning in variant.warnings
    )


def test_structure_checks_c_data_not_available_beyond_cds():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start=1200,
        end=1200,
        start_offset=1,
        end_offset=1,
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()

    validator.vr.validate.return_value = None

    validator.noreplace_myevm_t_to_g.side_effect = (
        vvhgvs.exceptions.HGVSDataNotAvailableError(
            "UTA missing"
        )
    )

    validator.hdp.get_tx_identity_info.return_value = [
        None,
        None,
        None,
        1,
        1000,
    ]

    assert structure_checks_c(variant, validator) is True

    assert any(
        "CDSError" in warning
        for warning in variant.warnings
    )


def test_structure_checks_c_g_to_t_alignment_failure():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()

    output = MagicMock()
    output.ac = "NC_000001.11"

    validator.noreplace_myevm_t_to_g.return_value = output

    variant.evm.g_to_t.side_effect = (
        vvhgvs.exceptions.HGVSError(
            "Alignment is incomplete"
        )
    )

    with patch(
        "VariantValidator.modules.use_checking."
        "hgvs_utils.incomplete_alignment_mapping_t_to_g",
        return_value=None,
    ):
        assert structure_checks_c(variant, validator) is True

    assert any(
        "Alignment is incomplete" in warning
        for warning in variant.warnings
    )


def test_structure_checks_c_invalid_variant_generic():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()

    validator.noreplace_myevm_t_to_g.side_effect = (
        vvhgvs.exceptions.HGVSInvalidVariantError(
            "Some HGVS validation error"
        )
    )

    assert structure_checks_c(variant, validator) is True
    assert "Some HGVS validation error" in variant.warnings


def test_structure_checks_c_uses_genomic_context():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()

    variant.genomic_context_ac = "NW_000001.1"

    output = MagicMock()
    output.ac = "NW_000001.1"

    check_ref_g = MagicMock()
    check_ref_g.ac = "NW_000001.1"

    no_replace_vm = MagicMock()
    no_replace_vm.t_to_g.return_value = check_ref_g
    no_replace_vm.g_to_t.return_value = variant.input_parses

    validator.vm.t_to_g.return_value = output
    variant.evm.g_to_t.return_value = variant.input_parses

    with patch(
        "vvhgvs.variantmapper.VariantMapper",
        return_value=no_replace_vm,
    ):
        result = structure_checks_c(variant, validator)

    assert result is False

    validator.vm.t_to_g.assert_called_once_with(
        variant.input_parses,
        "NW_000001.1",
    )

    no_replace_vm.t_to_g.assert_called_once()

    no_replace_vm.g_to_t.assert_called_once_with(
        check_ref_g,
        variant.input_parses.ac,
        alt_aln_method=validator.alt_aln_method,
    )


# ---------------------------------------------------------------------------
# structure_checks_n()
# ---------------------------------------------------------------------------

def test_structure_checks_n_data_not_available():
    variant = make_variant(
        ac="NR_000001.1",
        hgvs_type="n",
    )
    validator = make_validator()

    validator.vr.validate.side_effect = (
        vvhgvs.exceptions.HGVSDataNotAvailableError(
            "UTA missing"
        )
    )

    assert structure_checks_n(variant, validator) is True
    assert "UTA missing" in variant.warnings


def test_structure_checks_n_invalid_interval():
    variant = make_variant(
        ac="NR_000001.1",
        hgvs_type="n",
        start=20,
        end=10,
    )
    validator = make_validator()

    validator.vr.validate.side_effect = (
        vvhgvs.exceptions.HGVSInvalidVariantError(
            "base start position must be <= end position"
        )
    )

    assert structure_checks_n(variant, validator) is True

    assert any(
        "Interval start position" in warning
        for warning in variant.warnings
    )


def test_structure_checks_n_intronic_data_not_available():
    variant = make_variant(
        ac="NR_000001.1",
        hgvs_type="n",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
    )
    validator = make_validator()

    validator.vr.validate.return_value = None

    validator.noreplace_myevm_t_to_g.side_effect = (
        vvhgvs.exceptions.HGVSDataNotAvailableError(
            "UTA missing"
        )
    )

    assert structure_checks_n(variant, validator) is True

    assert any(
        "Universal Transcript Archive" in warning
        for warning in variant.warnings
    )


def test_structure_checks_n_intronic_value_error():
    variant = make_variant(
        ac="NR_000001.1",
        hgvs_type="n",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
    )
    validator = make_validator()

    validator.noreplace_myevm_t_to_g.side_effect = ValueError(
        "start > end"
    )

    assert structure_checks_n(variant, validator) is True

    assert any(
        "Interval start position" in warning
        for warning in variant.warnings
    )


def test_structure_checks_n_intronic_invalid_variant():
    variant = make_variant(
        ac="NR_000001.1",
        hgvs_type="n",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
    )
    validator = make_validator()

    validator.noreplace_myevm_t_to_g.side_effect = (
        vvhgvs.exceptions.HGVSInvalidVariantError(
            "Some HGVS validation error"
        )
    )

    assert structure_checks_n(variant, validator) is True
    assert "Some HGVS validation error" in variant.warnings


def test_structure_checks_n_uses_genomic_context():
    variant = make_variant(
        ac="NR_000001.1",
        hgvs_type="n",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
    )
    validator = make_validator()

    variant.genomic_context_ac = "NW_000001.1"

    output = MagicMock()
    output.ac = "NW_000001.1"

    validator.vm.t_to_g.return_value = output
    validator.vr.validate.return_value = None

    assert structure_checks_n(variant, validator) is False

    validator.vm.t_to_g.assert_called_with(
        variant.input_parses,
        "NW_000001.1",
    )

    validator.noreplace_myevm_t_to_g.assert_not_called()

def make_preparse_variant(
    original,
    quibble=None,
    format_invalid=False,
    formatted_quibble=None,
):
    variant = MagicMock()

    variant.original = original
    variant.quibble = original if quibble is None else quibble
    variant.warnings = []

    def format_quibble():
        if formatted_quibble is not None:
            variant.quibble = formatted_quibble
        return format_invalid

    variant.format_quibble.side_effect = format_quibble

    return variant


def test_preparse_numeric_only():
    variant = make_preparse_variant("  1:111111  ")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert variant.quibble == "1:111111"
    assert any(
        "only contains numeric characters" in warning
        for warning in variant.warnings
    )

    variant.format_quibble.assert_not_called()


def test_preparse_numeric_with_punctuation():
    variant = make_preparse_variant("12:30")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "only contains numeric characters" in warning
        for warning in variant.warnings
    )


def test_preparse_gene_symbol_only():
    variant = make_preparse_variant("BRCA1")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "only contains alphanumeric characters" in warning
        for warning in variant.warnings
    )


def test_preparse_alphanumeric_with_plus():
    variant = make_preparse_variant("BRCA1+TEST")

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "only contains alphanumeric characters" in warning
        for warning in variant.warnings
    )


def test_preparse_concatenated_descriptions():
    variant = make_preparse_variant(
        "NM_000001.1:c.1A>G:p.Arg1Gly"
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "concatenation" in warning
        for warning in variant.warnings
    )


def test_preparse_invalid_text_between_colons():
    variant = make_preparse_variant(
        "NM_000001.1:BAD:c.1A>G"
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert len(variant.warnings) == 2

    assert any(
        "single colon" in warning
        for warning in variant.warnings
    )

    assert any(
        "Illegal addition" in warning
        for warning in variant.warnings
    )


def test_preparse_accession_between_colons():
    variant = make_preparse_variant(
        "NM_000001.1:NM_000002.2:c.1A>G"
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "single colon" in warning
        for warning in variant.warnings
    )

    assert any(
        "Illegal addition" in warning
        for warning in variant.warnings
    )


def test_preparse_missing_substitution_arrow():
    variant = make_preparse_variant(
        "NM_000001.1:c.A123G"
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "Did you mean c.A>G?" in warning
        for warning in variant.warnings
    )


def test_preparse_format_quibble_comma_autocorrect():
    variant = make_preparse_variant(
        "NM_000001.1:c,123A>G",
        format_invalid=True,
    )

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble == "NM_000001.1:c.123A>G"

    assert any(
        "contained the , character" in warning
        for warning in variant.warnings
    )


def test_preparse_format_quibble_uppercase_type():
    variant = make_preparse_variant(
        "NM_000001.1:C.123A>G",
        format_invalid=True,
    )

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble == "NM_000001.1:c.123A>G"

    assert any(
        "wrong case" in warning
        for warning in variant.warnings
    )


def test_preparse_format_quibble_missing_dot():
    variant = make_preparse_variant(
        "NM_000001.1:c123A>G",
        format_invalid=True,
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "lacks the . character" in warning
        for warning in variant.warnings
    )


def test_preparse_strips_transcript_wrapper():
    variant = make_preparse_variant(
        "NG_000001.1(NM_000001.1):c.123A>G",
        format_invalid=True,
    )

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble == "NM_000001.1:c.123A>G"

    assert any(
        "Stripping unnecessary characters" in warning
        for warning in variant.warnings
    )


def test_preparse_invalid_format_raises():
    variant = make_preparse_variant(
        "NM_000001.1::???",
        format_invalid=True,
    )

    with pytest.raises(
        InvalidVariantError,
        match="not in an accepted format",
    ):
        pre_parsing_global_common_mistakes(variant)

    assert any(
        "not in an accepted format" in warning
        for warning in variant.warnings
    )


def test_preparse_bare_ins():
    variant = make_preparse_variant(
        "NM_000001.1:c.123_124ins"
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "inserted sequence must be provided" in warning
        for warning in variant.warnings
    )


def test_preparse_bare_ins_missing_two_positions():
    variant = make_preparse_variant(
        "NM_000001.1:c.123ins"
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "inserted sequence must be provided" in warning
        for warning in variant.warnings
    )

    assert any(
        "two positions" in warning
        for warning in variant.warnings
    )


def test_preparse_ins_parenthesised_length():
    variant = make_preparse_variant(
        "NM_000001.1:c.123_124ins(10)"
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "rewrite e.g. '(10)' to 'N[10]'" in warning
        for warning in variant.warnings
    )


def test_preparse_ins_numeric_length():
    variant = make_preparse_variant(
        "NM_000001.1:c.123_124ins10"
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "rewrite e.g. '10' to 'N[10]'" in warning
        for warning in variant.warnings
    )


def test_preparse_ins_uncertain_range_old_format():
    variant = make_preparse_variant(
        "NM_000001.1:c.123_124ins(10_20)"
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "N[(10_20)]" in warning
        for warning in variant.warnings
    )


def test_preparse_uncertain_insertion_reversed_range():
    variant = make_preparse_variant(
        "NM_000001.1:c.123_124ins[(20_10)]"
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "N[(10_20)]" in warning
        for warning in variant.warnings
    )


def test_preparse_uncertain_insertion_equal_range():
    variant = make_preparse_variant(
        "NM_000001.1:c.123_124ins[(10_10)]"
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "N[(10)]" in warning
        for warning in variant.warnings
    )


def test_preparse_valid_uncertain_insertion():
    variant = make_preparse_variant(
        "NM_000001.1:c.123_124ins[(10_20)]"
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "contains uncertainty" in warning
        for warning in variant.warnings
    )


def test_preparse_uncertain_insertion_missing_two_positions():
    variant = make_preparse_variant(
        "NM_000001.1:c.123ins[(10_20)]"
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "two positions" in warning
        for warning in variant.warnings
    )


def test_preparse_expands_simple_repeat_insertion():
    variant = make_preparse_variant(
        "NM_000001.1:c.123_124insAT[3]"
    )

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble == (
        "NM_000001.1:c.123_124insATATAT"
    )

    assert any(
        "may also be written as" in warning
        for warning in variant.warnings
    )


def test_preparse_expands_sectioned_repeat_insertion():
    variant = make_preparse_variant(
        "NM_000001.1:c.123_124ins[A[2];T[3]]"
    )

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble == (
        "NM_000001.1:c.123_124insAATTT"
    )

    assert any(
        "may also be written as" in warning
        for warning in variant.warnings
    )


def test_preparse_expands_sectioned_repeat_with_literal_section():
    variant = make_preparse_variant(
        "NM_000001.1:c.123_124ins[A[2];TG]"
    )

    assert pre_parsing_global_common_mistakes(variant) is False

    assert variant.quibble == (
        "NM_000001.1:c.123_124insAATG"
    )


def test_preparse_repeat_insertion_missing_two_positions():
    variant = make_preparse_variant(
        "NM_000001.1:c.123insA[3]"
    )

    assert pre_parsing_global_common_mistakes(variant) is False

    assert any(
        "two positions" in warning
        for warning in variant.warnings
    )


def test_preparse_substitution_with_range():
    variant = make_preparse_variant(
        "NM_000001.1:c.123_124A>G"
    )

    assert pre_parsing_global_common_mistakes(variant) is True

    assert any(
        "Base substitution (>) submitted with a reference sequence range"
        in warning
        for warning in variant.warnings
    )


def test_preparse_valid_description():
    variant = make_preparse_variant(
        "NM_000001.1:c.123A>G"
    )

    assert pre_parsing_global_common_mistakes(variant) is False
    assert variant.warnings == []

# ---------------------------------------------------------------------------
# Additional structure_checks() branch coverage
# ---------------------------------------------------------------------------

def test_structure_checks_c_check_true_intronic_bounds_remap():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start=123,
        end=123,
        start_offset=5,
        end_offset=5,
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()

    variant.warnings = [
        "Variant coordinate is beyond the bounds of the reference sequence"
    ]

    genomic = MagicMock()

    with patch(
        "VariantValidator.modules.use_checking.structure_checks_c",
        return_value=True,
    ), patch(
        "VariantValidator.modules.use_checking."
        "format_converters.remap_intronic",
    ) as remap:
        validator.myevm_t_to_g.return_value = genomic

        result = structure_checks(variant, validator)

    assert result is True
    assert variant.input_parses.posedit.pos.start.offset == 1
    assert variant.input_parses.posedit.pos.end.offset == 1

    validator.myevm_t_to_g.assert_called_once_with(
        variant.input_parses,
        variant.no_norm_evm,
        variant.primary_assembly,
        variant.hn,
        variant,
    )

    remap.assert_called_once_with(
        variant.input_parses,
        variant.input_parses,
        variant,
        validator,
    )



def test_structure_checks_c_check_true_non_intronic_no_remap():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start=123,
        end=123,
        start_offset=0,
        end_offset=0,
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()

    variant.warnings = [
        "Variant coordinate is beyond the bounds of the reference sequence"
    ]

    with patch(
        "VariantValidator.modules.use_checking.structure_checks_c",
        return_value=True,
    ), patch(
        "VariantValidator.modules.use_checking."
        "format_converters.remap_intronic",
    ) as remap:
        result = structure_checks(variant, validator)

    assert result is True
    remap.assert_not_called()


def test_structure_checks_c_check_true_without_bounds_no_remap():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start=123,
        end=123,
        start_offset=5,
        end_offset=5,
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()

    variant.warnings = ["Some other warning"]

    with patch(
        "VariantValidator.modules.use_checking.structure_checks_c",
        return_value=True,
    ), patch(
        "VariantValidator.modules.use_checking."
        "format_converters.remap_intronic",
    ) as remap:
        result = structure_checks(variant, validator)

    assert result is True
    remap.assert_not_called()


# ---------------------------------------------------------------------------
# Additional structure_checks_c() mapping coverage
# ---------------------------------------------------------------------------

def test_structure_checks_c_intronic_uses_primary_mapping_without_context():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()
    variant.genomic_context_ac = None

    genomic = MagicMock()
    genomic.ac = "NC_000001.11"

    no_replace_vm = MagicMock()
    no_replace_vm.t_to_g.return_value = genomic
    no_replace_vm.g_to_t.return_value = variant.input_parses

    validator.noreplace_myevm_t_to_g.return_value = genomic
    variant.evm.g_to_t.return_value = variant.input_parses

    with patch(
        "vvhgvs.variantmapper.VariantMapper",
        return_value=no_replace_vm,
    ):
        result = structure_checks_c(variant, validator)

    assert result is False
    validator.noreplace_myevm_t_to_g.assert_called_once_with(
        variant.input_parses,
        variant,
    )
    no_replace_vm.t_to_g.assert_called_once_with(
        variant.input_parses,
        variant.input_parses.ac,
        alt_aln_method=validator.alt_aln_method,
    )
    no_replace_vm.g_to_t.assert_called_once_with(
        genomic,
        variant.input_parses.ac,
        alt_aln_method=validator.alt_aln_method,
    )



def test_structure_checks_c_intronic_uses_genomic_context_mapping():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()
    variant.genomic_context_ac = "NW_000001.1"

    genomic = MagicMock()
    genomic.ac = "NW_000001.1"

    no_replace_vm = MagicMock()
    no_replace_vm.t_to_g.return_value = genomic
    no_replace_vm.g_to_t.return_value = variant.input_parses

    validator.vm.t_to_g.return_value = genomic
    variant.evm.g_to_t.return_value = variant.input_parses

    with patch(
        "vvhgvs.variantmapper.VariantMapper",
        return_value=no_replace_vm,
    ):
        result = structure_checks_c(variant, validator)

    assert result is False
    validator.vm.t_to_g.assert_called_once_with(
        variant.input_parses,
        variant.genomic_context_ac,
    )
    no_replace_vm.t_to_g.assert_called_once_with(
        variant.input_parses,
        variant.input_parses.ac,
        alt_aln_method=validator.alt_aln_method,
    )
    no_replace_vm.g_to_t.assert_called_once_with(
        genomic,
        variant.input_parses.ac,
        alt_aln_method=validator.alt_aln_method,
    )



def test_structure_checks_c_g_to_t_generic_error():
    variant = make_variant(
        ac="NM_000001.1",
        hgvs_type="c",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
        start_datum=Datum.CDS_START,
        end_datum=Datum.CDS_START,
    )
    validator = make_validator()

    genomic = MagicMock()
    genomic.ac = "NC_000001.11"

    validator.noreplace_myevm_t_to_g.return_value = genomic

    variant.evm.g_to_t.side_effect = (
        vvhgvs.exceptions.HGVSError(
            "generic transcript remapping failure"
        )
    )

    with patch(
        "VariantValidator.modules.use_checking."
        "hgvs_utils.incomplete_alignment_mapping_t_to_g",
        return_value=None,
    ):
        assert structure_checks_c(variant, validator) is True

    assert any(
        "generic transcript remapping failure" in warning
        for warning in variant.warnings
    )


# ---------------------------------------------------------------------------
# Additional structure_checks_n() mapping coverage
# ---------------------------------------------------------------------------

def test_structure_checks_n_intronic_uses_genomic_context():
    variant = make_variant(
        ac="NR_000001.1",
        hgvs_type="n",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
    )
    validator = make_validator()

    variant.genomic_context_ac = "NW_000001.1"

    genomic = MagicMock()
    genomic.ac = "NW_000001.1"

    validator.vm.t_to_g.return_value = genomic
    variant.evm.g_to_t.return_value = variant.input_parses

    assert structure_checks_n(variant, validator) is False

    validator.vm.t_to_g.assert_called_once_with(
        variant.input_parses,
        variant.genomic_context_ac,
    )


def test_structure_checks_n_intronic_uses_primary_mapping():
    variant = make_variant(
        ac="NR_000001.1",
        hgvs_type="n",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
    )
    validator = make_validator()

    variant.genomic_context_ac = None

    genomic = MagicMock()
    genomic.ac = "NC_000001.11"

    validator.noreplace_myevm_t_to_g.return_value = genomic

    assert structure_checks_n(variant, validator) is False

    validator.noreplace_myevm_t_to_g.assert_called_once_with(
        variant.input_parses,
        variant,
    )



def test_structure_checks_n_t_to_g_failure():
    variant = make_variant(
        ac="NR_000001.1",
        hgvs_type="n",
        start=123,
        end=123,
        start_offset=1,
        end_offset=1,
    )
    validator = make_validator()

    variant.genomic_context_ac = "NW_000001.1"

    validator.vm.t_to_g.side_effect = (
        vvhgvs.exceptions.HGVSError(
            "transcript mapping failed"
        )
    )

    with pytest.raises(
        vvhgvs.exceptions.HGVSError,
        match="transcript mapping failed",
    ):
        structure_checks_n(variant, validator)



# ---------------------------------------------------------------------------
# n. coordinate bounds
# ---------------------------------------------------------------------------

def test_structure_checks_n_negative_offset_from_base_one():
    variant = make_variant(
        ac="NR_000001.1",
        hgvs_type="n",
        start=1,
        end=1,
        start_offset=-1,
        end_offset=-1,
    )
    validator = make_validator()

    genomic = MagicMock()
    genomic.ac = "NC_000001.11"

    validator.myevm_t_to_g.return_value = genomic
    variant.hn.normalize.return_value = genomic

    with patch(
        "VariantValidator.modules.use_checking.fn.valstr",
        return_value="NR_000001.1:n.1-1",
    ):
        result = structure_checks_n(variant, validator)

    assert result is True

    assert any(
        "outside of the reference sequence" in warning
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
