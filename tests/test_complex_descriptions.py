import pytest
from unittest.mock import MagicMock, patch
import vvhgvs

from VariantValidator.modules.complex_descriptions import (
    FEInterval,
    fuzzy_ends,
    uncertain_positions,
    InvalidRangeError,
    IncompatibleTypeError,
    HgvsParseError,
    FuzzyPositionError,
)


def make_variant(quibble):
    v = MagicMock()
    v.quibble = quibble
    v.primary_assembly = "GRCh38"
    v.warnings = []
    v.hgvs_formatted = MagicMock()
    v.hgvs_formatted.posedit.pos = "?"
    return v


def make_validator():
    val = MagicMock()
    val.alt_aln_method = "splign"
    val.select_transcripts = "select"
    return val


# ---------- FEInterval ----------


def test_feinterval_format_none_start():
    iv = FEInterval(start=None, end=None)
    assert iv.format() == ""


def test_feinterval_format_same_start_end():
    pos = MagicMock()
    pos.format.return_value = "123"

    iv = FEInterval(start=pos, end=pos)

    assert iv.format() == "123"


def test_feinterval_format_range():
    start = MagicMock()
    end = MagicMock()

    start.format.return_value = "1"
    end.format.return_value = "5"

    iv = FEInterval(start=start, end=end)

    assert iv.format() == "(1)_(5)"


# ---------- fuzzy_ends ----------


def test_fuzzy_ends_unknown_start_end():
    variant = make_variant("NM_000001.1:c.(?_100)_(200_?)del")
    validator = make_validator()

    with pytest.raises(Exception):
        fuzzy_ends(variant, validator)

    assert len(variant.warnings) == 1


# ---------- uncertain_positions ----------


def test_uncertain_positions_returns_immediately():
    variant = make_variant("NM_000001.1:c.123A>G")
    validator = make_validator()

    assert uncertain_positions(variant, validator) is None


def test_uncertain_positions_invalid_missing_underscore():
    variant = make_variant("NM_000001.1:c.(1)(2)del")
    validator = make_validator()

    with pytest.raises(InvalidRangeError):
        uncertain_positions(variant, validator)


@patch("VariantValidator.modules.complex_descriptions.hgvs_obj_from_existing_edit")
def test_uncertain_positions_parse_error(mock_obj):
    variant = make_variant("NM_000001.1:c.(1_2)dup")
    validator = make_validator()

    import vvhgvs.exceptions

    mock_obj.side_effect = vvhgvs.exceptions.HGVSError("parse failed")

    with pytest.raises(HgvsParseError):
        uncertain_positions(variant, validator)

@patch("VariantValidator.modules.complex_descriptions.hgvs_obj_from_existing_edit")
def test_uncertain_positions_incompatible_type(mock_obj):
    variant = make_variant("NM_000001.1:c.(1_2)dup")
    validator = make_validator()

    parsed = MagicMock()
    parsed.posedit.pos = "1_2"
    mock_obj.return_value = parsed

    validator.vr.validate.side_effect = (
        vvhgvs.exceptions.HGVSError(
            "is not known to be compatible with variant type"
        )
    )

    with pytest.raises(IncompatibleTypeError):
        uncertain_positions(variant, validator)


@patch("VariantValidator.modules.complex_descriptions.hgvs_obj_from_existing_edit")
def test_uncertain_positions_base_start_gt_end(mock_obj):
    variant = make_variant("NM_000001.1:c.(5_2)dup")
    validator = make_validator()

    parsed = MagicMock()
    parsed.posedit.pos = "5_2"
    mock_obj.return_value = parsed

    validator.vr.validate.side_effect = (
        vvhgvs.exceptions.HGVSError(
            "base start position must be <= end position"
        )
    )

    with pytest.raises(InvalidRangeError) as exc:
        uncertain_positions(variant, validator)

    assert "base start position must be <= end position" in str(exc.value)


@patch("VariantValidator.modules.complex_descriptions.hgvs_obj_from_existing_edit")
def test_uncertain_positions_invalid_range(mock_obj):
    variant = make_variant("NM_000001.1:c.(5_8)dup")
    validator = make_validator()

    parsed = MagicMock()
    parsed.posedit.pos = "5_8"
    mock_obj.return_value = parsed

    validator.vr.validate.side_effect = (
        vvhgvs.exceptions.HGVSError(
            "some completely invalid range"
        )
    )

    with pytest.raises(InvalidRangeError) as exc:
        uncertain_positions(variant, validator)

    assert "invalid range" in str(exc.value)


@patch("VariantValidator.modules.complex_descriptions.hgvs_obj_from_existing_edit")
def test_uncertain_positions_intronic_range_passes_validation(mock_obj):
    variant = make_variant("NM_000001.1:c.(123+1_124-1)dup")
    validator = make_validator()

    parsed = MagicMock()
    parsed.posedit.pos = "123+1_124-1"
    mock_obj.return_value = parsed

    validator.vr.validate.side_effect = (
        vvhgvs.exceptions.HGVSError(
            "intronic + position"
        )
    )

    with pytest.raises(Exception):
        uncertain_positions(variant, validator)


@patch("VariantValidator.modules.complex_descriptions.hgvs_obj_from_existing_edit")
def test_uncertain_positions_success_sets_reformat_output(mock_obj):
    variant = make_variant("NM_000001.1:c.(1_2)dup")
    validator = make_validator()

    parsed = MagicMock()
    parsed.posedit.pos = "1_2"
    parsed.ac = "NM_000001.1"
    parsed.type = "c"

    mock_obj.return_value = parsed

    validator.vr.validate.return_value = None

    with patch.object(
        variant,
        "quibble",
        "NM_000001.1:c.(1_2)dup",
    ):
        try:
            uncertain_positions(variant, validator)
        except Exception:
            # Downstream mapping is expected to fail in this unit test.
            pass

    assert variant.reformat_output == "uncertain_pos"
    assert any(
        "Uncertain positions are not fully supported"
        in w
        for w in variant.warnings
    )


@patch("VariantValidator.modules.complex_descriptions.hgvs_obj_from_existing_edit")
def test_uncertain_positions_nc_sets_select(mock_obj):
    variant = make_variant("NC_000001.11:g.(100_200)del")
    validator = make_validator()

    validator.select_transcripts = "mane"

    parsed = MagicMock()
    parsed.ac = "NC_000001.11"
    parsed.type = "g"
    parsed.posedit.pos = MagicMock()

    mock_obj.return_value = parsed

    validator.vr.validate.return_value = None
    validator.relevant_transcripts.return_value = []

    uncertain_positions(variant, validator)

    assert validator.select_transcripts == "select"
    assert any(
        "Only a single transcript can be processed" in w
        for w in variant.warnings
    )


@patch("VariantValidator.modules.complex_descriptions.hgvs_obj_from_existing_edit")
def test_uncertain_positions_nc_intergenic(mock_obj):
    variant = make_variant("NC_000001.11:g.(100_200)del")
    validator = make_validator()

    parsed = MagicMock()
    parsed.ac = "NC_000001.11"
    parsed.type = "g"
    parsed.posedit.pos = MagicMock()

    mock_obj.return_value = parsed

    validator.vr.validate.return_value = None
    validator.select_transcripts = "select"
    validator.relevant_transcripts.return_value = []

    uncertain_positions(variant, validator)

    assert variant.output_type_flag == "intergenic"
    assert any(
        "Selected transcript does not span the entire range"
        in w
        for w in variant.warnings
    )


@patch("VariantValidator.modules.complex_descriptions.hgvs_obj_from_existing_edit")
def test_uncertain_positions_nc_selects_mane(mock_obj):
    variant = make_variant("NC_000001.11:g.(100_200)del")
    validator = make_validator()

    parsed = MagicMock()
    parsed.ac = "NC_000001.11"
    parsed.type = "g"
    parsed.posedit.pos = MagicMock()

    tx = MagicMock()
    tx.ac = "NM_000001.1"
    tx.type = "c"
    tx.posedit.pos = MagicMock()

    tx.__str__.return_value = "NM_000001.1:c.100_200del"

    mock_obj.side_effect = [parsed, tx]

    validator.vr.validate.return_value = None
    validator.select_transcripts = "select"
    validator.relevant_transcripts.return_value = [tx]

    validator.db.get_transcript_annotation.return_value = '{"select": "MANE"}'

    uncertain_positions(variant, validator)

    assert validator.select_transcripts == "NM_000001.1"
    assert variant.output_type_flag == "gene"


@patch("VariantValidator.modules.complex_descriptions.hgvs_obj_from_existing_edit")
def test_uncertain_positions_nc_selects_refseq(mock_obj):
    variant = make_variant("NC_000001.11:g.(100_200)del")
    validator = make_validator()

    parsed = MagicMock()
    parsed.ac = "NC_000001.11"
    parsed.type = "g"
    parsed.posedit.pos = MagicMock()

    tx = MagicMock()
    tx.ac = "NM_000001.1"
    tx.type = "c"
    tx.posedit.pos = MagicMock()
    tx.__str__.return_value = "NM_000001.1:c.100_200del"

    tx_variant = MagicMock()
    tx_variant.posedit.pos = MagicMock()

    mock_obj.side_effect = [parsed, tx_variant]

    validator.vr.validate.return_value = None
    validator.select_transcripts = "select"
    validator.relevant_transcripts.return_value = [tx]
    validator.db.get_transcript_annotation.return_value = '{"select": "RefSeq"}'

    uncertain_positions(variant, validator)

    assert validator.select_transcripts == "NM_000001.1"
    assert variant.output_type_flag == "gene"


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
