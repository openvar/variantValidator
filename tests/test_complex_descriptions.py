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
    FuzzyRangeError,
)

from unittest.mock import PropertyMock
from vvhgvs.exceptions import HGVSUnsupportedOperationError
from vvhgvs.enums import ValidationLevel

from unittest import TestCase

from VariantValidator import Validator

from VariantValidator.modules.format_converters import uncertain_pos
from VariantValidator.modules import complex_descriptions

# Functional tests
class TestComplexDescriptionsFunctional(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()

    def validate(self, variant, assembly="GRCh38", transcripts="all"):
        return self.vv.validate(
            variant,
            assembly,
            transcripts,
        ).format_as_dict(test=True)

    def test_uncertain_intronic_positions(self):
        results = self.validate(
            "NM_000546.6:c.(375+1_376-1)_(672+1_673-1)del"
        )

        assert results["NM_000546.6:c.(375+1_376-1)_(672+1_673-1)del"]["validation_warnings"] == [
            "Uncertain positions are not fully supported, however the syntax is valid"
        ]

    def test_uncertain_three_prime_utr_delins(self):
        results = self.validate(
            "NM_000546.6:c.(*1_*10)delinsA"
        )

        assert results["NM_000546.6:c.(*1_*10)delinsA"]["validation_warnings"] == [
            "Uncertain positions are not fully supported, however the syntax is valid"
        ]


    def test_uncertain_three_position_range(self):
        results = self.validate(
            "NM_000546.6:c.(100_200_300)delinsA"
        )

        assert results["validation_warning_1"]["validation_warnings"] == [
            "100_200_300 is an invalid range for accession NM_000546.6"
        ]


    def test_uncertain_out_of_order_range(self):
        results = self.validate(
            "NM_000546.6:c.(200_100)delinsA"
        )

        warnings = str(results)

        assert "base start position must be <= end position" in warnings.lower()


    def test_uncertain_single_intronic_position(self):
        results = self.validate(
            "NM_000546.6:c.(375+1)del"
        )

        assert results["NM_000546.6:c.375+1del"]["validation_warnings"] ==  [
            "Uncertain positions are not fully supported, however the syntax is valid"
        ]


    def test_uncertain_negative_intronic_position(self):
        results = self.validate(
            "NM_000546.6:c.(376-2_377-1)del"
        )

        assert results["validation_warning_1"]["validation_warnings"] ==  [
            "Uncertain positions are not fully supported, however the syntax is valid",
            "ExonBoundaryError: Position c.377-1 does not correspond with an exon boundary for transcript NM_000546.6"
        ]


    def test_lrg_transcript_uncertain_range(self):
        results = self.validate(
            "LRG_199t1:c.(100_200)del"
        )

        assert results[ "validation_warning_1"]["validation_warnings"] == [
            "Uncertain positions are not fully supported, however the syntax is valid",
            "This coding sequence variant description spans at least one intron"
        ]


    def test_lrg_genomic_uncertain_range(self):
        results = self.validate(
            "LRG_199:g.(100_200)del"
        )

        assert results["intergenic_variant_1"]["validation_warnings"] == [
            "Uncertain positions are not fully supported, however the syntax is valid",
            "NG_012232.1:g.(100_200)del automapped to genome position NC_000023.11:g.33344410_33344510del",
            "TranscriptIdentificationWarning: No individual transcripts have been identified that fully overlap the described variation in the genomic sequence. Large variants might span one or more genes and are currently only described at the genome (g.) level."
        ]

    def test_uncertain_intronic_position_not_at_exon_boundary(self):
        """Exercise the FuzzyPositionError exon-boundary branch."""
        results = self.validate(
            "NM_000546.6:c.(374+1_375+1)del"
        )

        validation = results["flag"]
        self.assertEqual(validation, "warning")

        warnings = " ".join(results["validation_warning_1"]["validation_warnings"])
        self.assertIn("ExonBoundaryError", warnings)


    def test_uncertain_non_numeric_position(self):
        """
        Exercise the second ValueError path where the coordinate cannot be
        converted to an integer after removing +/- offsets.
        """
        results = self.validate(
            "NM_000546.6:c.(?_100)del"
        )

        validation = results["flag"]
        self.assertEqual(validation, "warning")

        warnings = " ".join(results["validation_warning_1"]["validation_warnings"])
        self.assertTrue(
            "Fuzzy" in warnings or
            "unknown" in warnings or
            "Unsupported" in warnings
        )


    def test_uncertain_three_position_out_of_order(self):
        """
        Exercise the fuzzy-end ordering branch:
            if num1 <= num2 <= num3
            else:
                out of order
        """
        results = self.validate(
            "NM_000546.6:c.(200_100)delinsA"
        )

        validation = results["flag"]
        self.assertEqual(validation, "warning")

        warnings = " ".join(results["validation_warning_1"]["validation_warnings"])
        self.assertIn("base start position must be <= end position", warnings)

    def test_uncertain_5utr_single_position(self):
        result = self.validate("NM_000546.6:c.(-5)del")

        key = "NM_000546.6:c.-5del"

        self.assertEqual(result["flag"], "gene_variant")
        self.assertIn(key, result)

        data = result[key]

        self.assertEqual(
            data["hgvs_transcript_variant"],
            "NM_000546.6:c.-5del",
        )

        self.assertEqual(
            data["hgvs_predicted_protein_consequence"]["tlr"],
            "NP_000537.3:p.(=)",
        )

        self.assertIn(
            "Uncertain positions are not fully supported, however the syntax is valid",
            data["validation_warnings"],
        )

    def test_uncertain_5utr_range(self):
        result = self.validate("NM_000546.6:c.(-5_-3)del")

        key = "NM_000546.6:c.(-5_-3)del"

        self.assertEqual(result["flag"], "gene_variant")
        self.assertIn(key, result)

        data = result[key]

        warnings = data["validation_warnings"]

        self.assertIn(
            "Uncertain positions are not fully supported, however the syntax is valid",
            warnings,
        )

        self.assertTrue(
            any(
                "automapped to NM_000546.6:c.(-4_-2)del" in w
                for w in warnings
            )
        )

        self.assertEqual(
            data["hgvs_predicted_protein_consequence"]["tlr"],
            "NP_000537.3:p.(=)",
        )

    def test_uncertain_5utr_to_coding_range(self):
        result = self.validate("NM_000546.6:c.(-5_10)del")

        key = "NM_000546.6:c.(-5_10)del"

        self.assertEqual(result["flag"], "gene_variant")
        self.assertIn(key, result)

        data = result[key]

        warnings = data["validation_warnings"]

        self.assertIn(
            "Uncertain positions are not fully supported, however the syntax is valid",
            warnings,
        )

        self.assertTrue(
            any(
                "automapped to NM_000546.6:c.(-4_11)del" in w
                for w in warnings
            )
        )

        self.assertTrue(
            any(
                "ProteinTranslationError" in w
                for w in warnings
            )
        )

        self.assertEqual(
            data["hgvs_predicted_protein_consequence"]["tlr"],
            "",
        )


# Unit tests
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


# ---------- FEInterval.validate ----------

def test_feinterval_validate_valid():
    start = MagicMock()
    end = MagicMock()

    start.validate.return_value = (ValidationLevel.VALID, None)
    end.validate.return_value = (ValidationLevel.VALID, None)

    start.start = 10
    end.end = 20

    iv = FEInterval(start=start, end=end)

    assert iv.validate() == (ValidationLevel.VALID, None)


def test_feinterval_validate_invalid_order():
    start = MagicMock()
    end = MagicMock()

    start.validate.return_value = (ValidationLevel.VALID, None)
    end.validate.return_value = (ValidationLevel.VALID, None)

    start.start = 20
    end.end = 10

    iv = FEInterval(start=start, end=end)

    assert iv.validate() == (
        ValidationLevel.ERROR,
        "base start position must be <= end position",
    )


def test_feinterval_validate_warning():
    start = MagicMock()
    end = MagicMock()

    start.validate.return_value = (ValidationLevel.VALID, None)
    end.validate.return_value = (ValidationLevel.VALID, None)

    type(start).start = PropertyMock(
        side_effect=HGVSUnsupportedOperationError("unsupported")
    )
    end.end = 10

    iv = FEInterval(start=start, end=end)

    res, msg = iv.validate()

    assert res == ValidationLevel.WARNING
    assert "unsupported" in msg


# ---------- fuzzy_ends ----------

@patch("VariantValidator.modules.complex_descriptions.vv_utils.get_exon_boundary_list")
def test_fuzzy_ends_ordered_star_positions(mock_boundaries):
    variant = make_variant("NM_000001.1:c.(*1_*2)_(*3_?)del")
    validator = make_validator()

    validator.hdp.get_tx_identity_info.return_value = (
        None,
        None,
        None,
        None,
        100,
    )

    mock_boundaries.return_value = ["101", "102", "103", "NC_000001.11"]

    with pytest.raises(FuzzyRangeError):
        fuzzy_ends(variant, validator)

    assert any(
        "syntax is valid" in w
        for w in variant.warnings
    )


@patch("VariantValidator.modules.complex_descriptions.vv_utils.get_exon_boundary_list")
def test_fuzzy_ends_out_of_order_star_positions(mock_boundaries):
    variant = make_variant("NM_000001.1:c.(*5_*4)_(*3_?)del")
    validator = make_validator()

    validator.hdp.get_tx_identity_info.return_value = (
        None,
        None,
        None,
        None,
        100,
    )

    mock_boundaries.return_value = ["103", "104", "105", "NC_000001.11"]

    with pytest.raises(FuzzyRangeError):
        fuzzy_ends(variant, validator)

    assert any(
        "out of order" in w
        for w in variant.warnings
    )


def test_fuzzy_position_end_only():
    variant = make_variant("NM_000001.1:c.123A>G")

    pos = MagicMock()
    pos.__str__.return_value = "123_?"

    pos.start.__str__.return_value = "123"
    pos.end.__str__.return_value = "?"

    variant.hgvs_formatted.posedit.pos = pos

    with pytest.raises(FuzzyPositionError):
        fuzzy_ends(variant, make_validator())


def test_fuzzy_position_start_only():
    variant = make_variant("NM_000001.1:c.123A>G")

    pos = MagicMock()
    pos.__str__.return_value = "?_123"

    pos.start.__str__.return_value = "?"
    pos.end.__str__.return_value = "123"

    variant.hgvs_formatted.posedit.pos = pos

    with pytest.raises(FuzzyPositionError):
        fuzzy_ends(variant, make_validator())


def test_fuzzy_position_both():
    variant = make_variant("NM_000001.1:c.123A>G")

    pos = MagicMock()
    pos.__str__.return_value = "?_?"

    pos.start.__str__.return_value = "?"
    pos.end.__str__.return_value = "?"

    variant.hgvs_formatted.posedit.pos = pos

    with pytest.raises(FuzzyPositionError):
        fuzzy_ends(variant, make_validator())




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
