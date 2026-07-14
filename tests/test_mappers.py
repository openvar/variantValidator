from unittest import TestCase
import json

from VariantValidator.validator import Validator

class TestMappersFunctional(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()

    def validate(self, variant, assembly="GRCh38", transcripts="all"):
        return self.vv.validate(
            variant,
            assembly,
            transcripts,
        ).format_as_dict(test=True)

    def test_issue518_forward_intron_correction(self):
        result = self.validate(
            "NM_000086.2:c.791-802_1056+1445del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["NM_000086.2:c.790+532_1056+1445del"]["validation_warnings"] == [
            "ExonBoundaryError: Position c.791-802 has been updated to position to 790+532 ensuring correct HGVS numbering for transcript NM_000086.2"
        ]

    def test_issue518_reverse_intron_correction(self):
        result = self.validate(
            "NM_000088.4:c.2559+54_2560del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["NM_000088.4:c.2560-34_2561del"]["validation_warnings"] == [
            "ExonBoundaryError: Position c.2559+54 has been updated to position to 2560-35 ensuring correct HGVS numbering for transcript NM_000088.4"
        ]

    def test_issue518_double_boundary_correction(self):
        result = self.validate(
            "NM_000088.4:c.2559+54_2559+55del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["NM_000088.4:c.2560-34_2560-33del"]["validation_warnings"] == [
            "ExonBoundaryError: Position c.2559+54 has been updated to position to 2560-35 ensuring correct HGVS numbering for transcript NM_000088.4",
            "ExonBoundaryError: Position c.2559+55 has been updated to position to 2560-34 ensuring correct HGVS numbering for transcript NM_000088.4"
        ]

    def test_issue518_left_boundary_correction(self):
        result = self.validate(
            "NM_000086.2:c.790_791-802del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["NM_000086.2:c.790+1_790+533del"]["validation_warnings"] ==  [
            "ExonBoundaryError: Position c.791-802 has been updated to position to 790+532 ensuring correct HGVS numbering for transcript NM_000086.2",
            "VariantNormalizationWarning: NM_000086.2:c.790_790+532del normalized to NM_000086.2:c.790+1_790+533del"
        ]

    def test_issue518_two_reverse_boundaries(self):
        result = self.validate(
            "NM_000086.2:c.791-803_791-802del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result[ "NM_000086.2:c.790+531_790+532del"]["validation_warnings"] ==  [
            "ExonBoundaryError: Position c.791-803 has been updated to position to 790+531 ensuring correct HGVS numbering for transcript NM_000086.2",
            "ExonBoundaryError: Position c.791-802 has been updated to position to 790+532 ensuring correct HGVS numbering for transcript NM_000086.2"
        ]

    def test_issue518_invalid_start_boundary(self):
        result = self.validate(
            "NM_000088.4:c.543-1_543+1del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["validation_warning_1"]["validation_warnings"] == [
            "ExonBoundaryError: Position c.543-1 does not correspond with an exon boundary for transcript NM_000088.4"
        ]

    def test_issue518_invalid_end_boundary(self):
        result = self.validate(
            "NM_000088.4:c.2560-3_2560+3dup"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["validation_warning_1"]["validation_warnings"] == [
            "ExonBoundaryError: Position c.2560+3 does not correspond with an exon boundary for transcript NM_000088.4"
        ]

    def test_issue518_invalid_end_boundary2(self):
        result = self.validate(
            "NM_000088.4:c.2608_2608+5del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["validation_warning_1"]["validation_warnings"] == [
            "ExonBoundaryError: Position c.2608+5 does not correspond with an exon boundary for transcript NM_000088.4"
        ]

    def test_issue518_interval_error(self):
        result = self.validate(
            "NM_000088.4:c.543_453+1del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["validation_warning_1"]["validation_warnings"] == [
            "IntervalOrderError: Interval end position 453 < interval start position 543"
        ]

    def test_issue518_reverse_end_only(self):
        result = self.validate(
            "NM_000088.4:c.2559_2559+54del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["NM_000088.4:c.2559_2560-35del"]["validation_warnings"] == [
            "ExonBoundaryError: Position c.2559+54 has been updated to position to 2560-35 ensuring correct HGVS numbering for transcript NM_000088.4"
        ]