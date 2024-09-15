"""
expanded_repeats_tests
Authors: Robert Wilson (@RSWilson1) and Rebecca Locke (@rklocke)
This code runs tests on the module expanded_repeats.py to check the outputs are as expected.

It checks known edge-case HGVS compliant variant strings.
Additional functionality to add:
- Check error handling
- Check correct errors for non-HGVS compliant strings.

"""
import unittest
from unittest import TestCase
from VariantValidator.modules import expanded_repeats
from VariantValidator import Validator
vv = Validator()
vv.alt_aln_method = "splign"


class TestExpandedRepeats(unittest.TestCase):
    """Tests for the internal expanded_repeats.py module, to directly check
    that the syntax checker returns the expected results for each variant case.
    Including known edge-cases that weren't previously handled.

    Attributes
    ----------
    Variants with known strings and expected results.
    Returns
    ----------
    Number of tests completed successfully.
    """

    def test_basic_syntax_RSG(self):
        """
        Test for handling basic syntax of variant string.
        """
        variant_str = "NG_012232.1:g.4T[20]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(variant_str,  "GRCh37", "all", vv)
        my_variant.check_transcript_type()
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat(vv)
        assert formatted == "NG_012232.1:g.3_6T[20]"
        assert my_variant.variant_str == "NG_012232.1:g.4T[20]"
        assert my_variant.ref_type == "RefSeq"
        # checks correct transcript type
        assert my_variant.reference == "NG_012232.1"
        # checks correct position
        assert my_variant.variant_position == "3_6"
        # checks repeat seq
        assert my_variant.repeat_sequence == "T"
        # checks correct suffix
        assert my_variant.copy_number == "20"
        # checks number of repeats is str and correct
        assert my_variant.after_the_bracket == ""
        # checks nothing is after the bracket

    def test_basic_syntax_ENSG(self):
        """
        Test for handling basic syntax of ENSG variant string.
        """
        variant_str = "ENST00000263121.12:c.1082TCT[2]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str,  "GRCh37", "all", vv)
        my_variant.check_transcript_type()
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat(vv)
        assert formatted == "ENST00000263121.12:c.1082_1087TCT[2]"
        assert my_variant.variant_str == "ENST00000263121.12:c.1082TCT[2]"
        assert my_variant.prefix == "c"
        assert my_variant.ref_type == "Ensembl"
        # checks correct transcript type
        assert my_variant.reference == "ENST00000263121.12"
        # checks correct ref name
        assert my_variant.variant_position == "1082_1087"
        # checks correct position
        assert my_variant.repeat_sequence == "TCT"
        # checks repeat seq
        assert my_variant.copy_number == "2"
        # checks number of repeats is str and correct
        assert my_variant.after_the_bracket == ""
        # checks nothing is after the bracket

    def test_basic_syntax_NM(self):
        """
        Test for handling basic syntax with a NM_ 'c' type variant string.
        """
        variant_str = "NM_000492.4:c.1210-34TG[11]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(variant_str, "GRCh38", "all", vv)
        assert my_variant.variant_position == "1210-34"
        my_variant.check_transcript_type()
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat(vv)
        assert formatted == "NM_000492.4:c.1210-34_1210-13TG[11]"
        assert my_variant.variant_str == "NM_000492.4:c.1210-34TG[11]"
        assert my_variant.prefix == "c"
        assert my_variant.ref_type == "RefSeq"
        # checks correct transcript type
        assert my_variant.reference == "NM_000492.4"
        # checks correct ref name
        assert my_variant.variant_position == "1210-34_1210-13"
        # checks correct position
        assert my_variant.repeat_sequence == "TG"
        # checks repeat seq
        assert my_variant.copy_number == "11"
        # checks number of repeats is str and correct
        assert my_variant.after_the_bracket == ""
        # checks nothing is after the bracket

    def test_getting_full_range_from_single_pos(self):
        """
        Test to full range is calculated correctly
        """
        variant_str = "NM_003073.5:c.1085AGA[2]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str,  "GRCh37", "all", vv)
        # Includes cds offset
        self.assertEqual(expanded_repeats.TandemRepeats.get_range_from_single_pos(my_variant, vv), "1289_1297")

    def test_empty_string(self):
        """
        Test for handling empty string.
        """
        variant_str = ""
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str, "GRCh37", "all", vv)
        assert my_variant == False

    def test_simple_str_split(self):
        test_variant = "ENSG00000198947:g.1ACT[10]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(test_variant,
                                                                         "GRCh37", "all", vv)
        my_variant.simple_split_string()
        assert my_variant.begining == "ENSG00000198947"
        assert my_variant.end == ":g.1ACT[10]"


class TestCVariantsExpanded(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()
        cls.vv.testing = True

    def test_exon_boundary_single_position(self):
        variant = 'NM_004006.2:c.13-14AC[7]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert ("ExonBoundaryError: Position 13-14 does not correspond with an exon boundary for transcript "
                "NM_004006.2") in results["validation_warning_1"]["validation_warnings"]

    def test_exon_boundary_range(self):
        variant = 'NM_004006.2:c.12_13-14AC[7]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert ("ExonBoundaryError: Stated position does not correspond with an exon boundary "
                "for transcript NM_004006.2") in results["validation_warning_1"]["validation_warnings"]

    def test_intronic_single_position(self):
        variant = 'NM_000492.4:c.1210-34TG[11]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_000492.4:c.1210-34_1210-13="][
            "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000007.13:g.117188661_117188682="
        assert results["NM_000492.4:c.1210-34_1210-13="][
            "validation_warnings"] == [
            "ExpandedRepeatWarning: NM_000492.4:c.1210-34TG[11] updated to NM_000492.4:c.1210-34_1210-13TG[11]",
            "ExpandedRepeatWarning: NM_000492.4:c.1210-34_1210-13TG[11] should only be used as an annotation for the "
            "core HGVS descriptions provided",
            "NM_000492.4:c.1210-34_1210-13delinsTGTGTGTGTGTGTGTGTGTGTG automapped to NM_000492.4:c.1210-34_1210-13="
        ]

    def test_intronic_range(self):
        variant = 'NM_000492.4:c.1210-34_1210-13TG[11]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_000492.4:c.1210-34_1210-13="][
            "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000007.13:g.117188661_117188682="
        assert results["NM_000492.4:c.1210-34_1210-13="][
            "validation_warnings"] == [
            "ExpandedRepeatWarning: NM_000492.4:c.1210-34_1210-13TG[11] should only be used as an annotation for the "
            "core HGVS descriptions provided",
            "NM_000492.4:c.1210-34_1210-13delinsTGTGTGTGTGTGTGTGTGTGTG automapped to NM_000492.4:c.1210-34_1210-13="
        ]

    def test_incorrect_intronic_range(self):
        variant = 'NM_000492.4:c.1210-35_1210-14TG[11]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] == [
            "RepeatSyntaxError: The repeat sequence does not match the reference sequence at the given position "
            "1210-35_1210-14, expected TGTGTGTGTGTGTGTGTGTGTG but the reference is ATGTGTGTGTGTGTGTGTGTGT "
            "at the specified position"
        ]

    def test_exonic_single_position(self):
        variant = 'NM_003073.5:c.1085AGA[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_003073.5:c.1085_1093="][
            "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000022.10:g.24175857_24175865="
        assert results["NM_003073.5:c.1085_1093="][
            "validation_warnings"] == [
            "ExpandedRepeatWarning: NM_003073.5:c.1085AGA[3] updated to NM_003073.5:c.1085_1093AGA[3]",
            "ExpandedRepeatWarning: NM_003073.5:c.1085_1093AGA[3] should only be used as an annotation for the core "
            "HGVS descriptions provided"
        ]

    def test_exonic_range(self):
        variant = 'NM_003073.5:c.1085_1093AGA[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_003073.5:c.1085_1093="][
            "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000022.10:g.24175857_24175865="
        assert results["NM_003073.5:c.1085_1093="][
            "validation_warnings"] == [
            "ExpandedRepeatWarning: NM_003073.5:c.1085_1093AGA[3] should only be used as an annotation for the core "
            "HGVS descriptions provided"
        ]

    def test_incorrect_exonic_range(self):
        variant = 'NM_003073.5:c.1086_1094AGA[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"][
            "validation_warnings"] == [
            "RepeatSyntaxError: The provided repeat sequence AGA does not match the reference sequence GAA "
            "at the given position 1290_1292 of reference sequence NM_003073.5"
        ]

    def test_5_utr_single_pos(self):
        variant = 'NM_002024.5:c.-129CGG[10]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_002024.5:c.-129_-100="][
            "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000023.10:g.146993569_146993598="
        assert results["NM_002024.5:c.-129_-100="][
            "validation_warnings"] == [
            "ExpandedRepeatWarning: NM_002024.5:c.-129CGG[10] updated to NM_002024.5:c.-129_-100CGG[10]",
            "ExpandedRepeatWarning: NM_002024.5:c.-129_-100CGG[10] should only be used as an annotation for "
            "the core HGVS descriptions provided",
            "TranscriptVersionWarning: A more recent version of the selected reference sequence NM_002024.5 "
            "is available for genome build GRCh37 (NM_002024.6)"
        ]

    def test_5_utr_range(self):
        variant = 'NM_002024.5:c.-129_-100CGG[10]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_002024.5:c.-129_-100="][
            "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000023.10:g.146993569_146993598="
        assert results["NM_002024.5:c.-129_-100="][
            "validation_warnings"] == [
            "ExpandedRepeatWarning: NM_002024.5:c.-129_-100CGG[10] should only be used as an annotation for the core "
            "HGVS descriptions provided",
            "TranscriptVersionWarning: A more recent version of the selected reference sequence NM_002024.5 is "
            "available for genome build GRCh37 (NM_002024.6)"
        ]

    def test_5_utr_single_pos_incorrect(self):
        variant = 'NM_002024.5:c.-128CGG[10]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"][
            "validation_warnings"] == [
            "RepeatSyntaxError: The provided repeat sequence CGG does not match the reference sequence GGC at "
            "the given position 102_104 of reference sequence NM_002024.5"
        ]

    def test_gap_crossing(self):
        """Test that when the reference genome has a del in it WRT the transcript that
        the transcript version is expanded to match the repeat across this region anyway,
        and that the genome gets the correct coordinates. The same result should be
        found for transcripts described via LRG IDs, even without the gap being in the LRG
        due to the reference centric nature of our current tooling but will now tested
        later in the LRG section of the tests."""
        variant = 'NM_002111.8:c.54_110GCA[21]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_002111.8:c.54_116="]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000004.11:g.3076657_3076662dup"
        assert results["NM_002111.8:c.54_116="]["validation_warnings"] ==  [
            "ExpandedRepeatWarning: NM_002111.8:c.54_110GCA[21] updated to NM_002111.8:c.54_116GCA[21]",
            "ExpandedRepeatWarning: NM_002111.8:c.54_116GCA[21] should only be used as an "
            "annotation for the core HGVS descriptions provided"
        ]

    def test_antisense_intron_range(self):
        variant = 'NM_000088.3:c.589-1_590G[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_000088.3:c.589-1_590="][
            "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000017.10:g.48275362_48275364="
        assert results["NM_000088.3:c.589-1_590="][
            "validation_warnings"] == [
            "ExpandedRepeatWarning: NM_000088.3:c.589-1_590G[3] should only be used as an annotation for the "
            "core HGVS descriptions provided",
            "NM_000088.3:c.589-1_590delinsGGG automapped to NM_000088.3:c.589-1_590=",
            "TranscriptVersionWarning: A more recent version of the selected reference sequence NM_000088.3 "
            "is available for genome build GRCh37 (NM_000088.4)"
        ]

    def test_antisense_intron_single_pos(self):
        variant = 'NM_000088.3:c.589-1G[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_000088.3:c.589-1_590="][
                   "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000017.10:g.48275362_48275364="
        assert results["NM_000088.3:c.589-1_590="][
                   "validation_warnings"] == [
            "ExpandedRepeatWarning: NM_000088.3:c.589-1G[3] updated to NM_000088.3:c.589-1_590G[3]",
            "ExpandedRepeatWarning: NM_000088.3:c.589-1_590G[3] should only be used as an annotation for the "
            "core HGVS descriptions provided",
            "NM_000088.3:c.589-1_590delinsGGG automapped to NM_000088.3:c.589-1_590=",
            "TranscriptVersionWarning: A more recent version of the selected reference sequence NM_000088.3 "
            "is available for genome build GRCh37 (NM_000088.4)"
        ]

    def test_antisense_intron_single_pos_2(self):
        variant = 'NM_000088.3:c.589-18T[5]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_000088.3:c.589-18_589-14="][
                   "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000017.10:g.48275377_48275381="
        assert results["NM_000088.3:c.589-18_589-14="][
                   "validation_warnings"] == [
            "ExpandedRepeatWarning: NM_000088.3:c.589-18T[5] updated to NM_000088.3:c.589-18_589-14T[5]",
            "ExpandedRepeatWarning: NM_000088.3:c.589-18_589-14T[5] should only be used as an annotation for "
            "the core HGVS descriptions provided",
            "NM_000088.3:c.589-18_589-14delinsTTTTT automapped to NM_000088.3:c.589-18_589-14=",
            "TranscriptVersionWarning: A more recent version of the selected reference sequence NM_000088.3 "
            "is available for genome build GRCh37 (NM_000088.4)"
        ]

    def test_antisense_intron_single_pos_seq_inverted(self):
        variant = 'NM_000088.3:c.589-18A[5]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"][
                   "validation_warnings"] == [
            "RepeatSyntaxError: The provided repeat sequence T does not match the reference sequence A"
            " at the given position 48275381_48275381 of reference sequence NC_000017.10"
        ]


if __name__ == "__main__":
    unittest.main()


# <LICENSE>
# Copyright (C) 2016-2024 VariantValidator Contributors
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

