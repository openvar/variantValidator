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
from VariantValidator.modules import expanded_repeats


class TestExpandedRepeats(unittest.TestCase):
    """Tests for expanded_repeats.py to check the syntax checker returns
    the expected results for each variant case.
    Including known edge-cases that weren't previously handled.

    Attributes
    ----------
    Variants with known strings and expected results.
    Returns
    ----------
    Number of tests completed successfully.
    """

    def test_basic_syntax(self):
        """
        Test for handling basic syntax of variant string.
        """
        variant_str = "LRG_199:g.1ACT[20]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(variant_str, "GRCh37", "all")
        my_variant.check_transcript_type()
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat()
        assert formatted == "LRG_199:g.1ACT[20]"
        assert my_variant.variant_str == "LRG_199:g.1ACT[20]"
        assert my_variant.ref_type == "LRG"
        # checks correct transcript type
        assert my_variant.reference == "LRG_199"
        # checks correct position
        assert my_variant.variant_position == "1"
        # checks repeat seq
        assert my_variant.repeat_sequence == "ACT"
        # checks correct suffix
        assert my_variant.copy_number == "20"
        # checks number of repeats is str and correct
        assert my_variant.after_the_bracket == ""
        # checks nothing is after the bracket

    def test_basic_syntax_2(self):
        """
        Test for handling basic syntax of ENSG variant string.
        """
        variant_str = "ENSG00000198947.15:g.1ACT[10]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str, "GRCh37", "all")
        my_variant.check_transcript_type()
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat()
        assert formatted == "ENSG00000198947.15:g.1ACT[10]"
        assert my_variant.variant_str == "ENSG00000198947.15:g.1ACT[10]"
        assert my_variant.ref_type == "Ensembl"
        # checks correct transcript type
        assert my_variant.reference == "ENSG00000198947.15"
        # checks correct ref name
        assert my_variant.variant_position == "1"
        # checks correct position
        assert my_variant.repeat_sequence == "ACT"
        # checks prefix for variant
        assert my_variant.prefix == "g"
        # checks repeat seq
        assert my_variant.copy_number == "10"
        # checks number of repeats is str and correct
        assert my_variant.after_the_bracket == ""
        # checks nothing is after the bracket

    def test_basic_syntax_3(self):
        """
        Test for handling basic syntax of ENSG variant string.
        """
        variant_str = "ENST00000357033.8:c.13AC[2]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str, "GRCh37", "all")
        my_variant.check_transcript_type()
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat()
        assert formatted == "ENST00000357033.8:c.13_14insACAC"
        assert my_variant.variant_str == "ENST00000357033.8:c.13AC[2]"
        assert my_variant.prefix == "c"
        assert my_variant.ref_type == "Ensembl"
        # checks correct transcript type
        assert my_variant.reference == "ENST00000357033.8"
        # checks correct ref name
        assert my_variant.variant_position == "13"
        # checks correct position
        assert my_variant.repeat_sequence == "AC"
        # checks repeat seq
        assert my_variant.copy_number == "2"
        # checks number of repeats is str and correct
        assert my_variant.after_the_bracket == ""
        # checks nothing is after the bracket

    def test_basic_syntax_4(self):
        """
        Test for handling basic syntax of ENSG variant string.
        """
        variant_str = "ENSG00000198947.15:g.1_10A[10]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str, "GRCh37", "all")
        my_variant.check_transcript_type()
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat()
        assert formatted == "ENSG00000198947.15:g.1_10A[10]"
        assert my_variant.variant_str == "ENSG00000198947.15:g.1_10A[10]"
        assert my_variant.ref_type == "Ensembl"
        # checks correct transcript type
        assert my_variant.reference == "ENSG00000198947.15"
        # checks correct ref name
        assert my_variant.variant_position == "1_10"
        # checks correct position
        assert my_variant.repeat_sequence == "A"
        # checks prefix for variant
        assert my_variant.prefix == "g"
        # checks repeat seq
        assert my_variant.copy_number == "10"
        # checks number of repeats is str and correct
        assert my_variant.after_the_bracket == ""
        # checks nothing is after the bracket

    def test_basic_syntax_5(self):
        """
        Test for handling basic syntax of ENSG variant string.
        """
        variant_str = "LRG_199t1:c.13-25ACT[5]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str, "GRCh37", "all")
        assert my_variant.variant_position == "13_25"
        my_variant.check_transcript_type()
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat()
        assert formatted == "LRG_199t1:c.13_27ACT[5]"
        assert my_variant.variant_str == "LRG_199t1:c.13-25ACT[5]"
        assert my_variant.prefix == "c"
        assert my_variant.ref_type == "LRG"
        # checks correct transcript type
        assert my_variant.reference == "LRG_199t1"
        # checks correct ref name
        assert my_variant.variant_position == "13_27"
        # checks correct position
        assert my_variant.repeat_sequence == "ACT"
        # checks repeat seq
        assert my_variant.copy_number == "5"
        # checks number of repeats is str and correct
        assert my_variant.after_the_bracket == ""
        # checks nothing is after the bracket

    def test_getting_full_range_from_single_pos(self):
        """
        Test to full range is calculated correctly
        """
        variant_str = "LRG_199t1:c.13ACT[5]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str, "GRCh37", "all")
        self.assertEqual(expanded_repeats.TandemRepeats.get_range_from_single_pos(my_variant), "13_27")

    def test_LRG_transcript_handling(self):
        """
        Test for transcript_handling of
        LRG variant string.
        """
        # Gives LRG_199t1:c.1_60ACT[20]
        variant_str = "LRG_199t1:c.1_2ACT[20]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str, "GRCh37", "all")
        my_variant.check_transcript_type()
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat()
        assert formatted == "LRG_199t1:c.1_60ACT[20]"
        assert my_variant.variant_str == "LRG_199t1:c.1_2ACT[20]"
        assert my_variant.ref_type == "LRG"
        # checks correct transcript type
        assert my_variant.reference == "LRG_199t1"
        # checks correct ref name
        assert my_variant.variant_position == "1_60"
        # checks correct position
        assert my_variant.repeat_sequence == "ACT"
        # checks prefix for variant
        assert my_variant.prefix == "c"
        # checks repeat seq
        assert my_variant.copy_number == "20"
        # checks number of repeats is str and correct
        assert my_variant.after_the_bracket == ""
        # checks nothing is after the bracket

    def test_transcript_versions(self):
        """
        Test for handing transcript versions.
        Previous code gave an error below:
        AssertionError: Unable to identify a colon (:) in the variant
        """
        variant_str = "NM_004006.2:c.13AC[7]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str, "GRCh37", "all")
        my_variant.check_transcript_type()
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat()
        assert formatted == "NM_004006.2:c.13_14insACACACACACACAC"
        assert my_variant.variant_str == "NM_004006.2:c.13AC[7]"
        assert my_variant.ref_type == "RefSeq"
        # checks correct transcript type
        assert my_variant.reference == "NM_004006.2"
        # checks correct ref name
        assert my_variant.variant_position == "13"
        # checks correct position
        assert my_variant.repeat_sequence == "AC"
        # checks repeat seq
        assert my_variant.copy_number == "7"
        # checks number of repeats is str and correct
        assert my_variant.after_the_bracket == ""
        # checks nothing is after the bracket

    def test_throws_exception(self):
        # Test throws AssertionError if no colon included in variant
        test_variant = "NG_004006.2g.1_2act[22]"
        with self.assertRaises(AssertionError):
            expanded_repeats.TandemRepeats.parse_repeat_variant(test_variant, "GRCh37", "all")

    def test_throws_exception_2(self):
        # Test throws AssertionError if no repeat sequence is included
        test_variant = "ENST00000198947.1:c.1_2[10]"
        with self.assertRaises(AssertionError):
            expanded_repeats.TandemRepeats.parse_repeat_variant(test_variant, "GRCh37","all")

    def test_throws_exception_3(self):
        # Test throws AssertionError if allele variant
        test_variant = "LRG_199:g.[123456A>G];[345678G>C]"
        with self.assertRaises(AssertionError):
            expanded_repeats.TandemRepeats.parse_repeat_variant(test_variant, "GRCh37", "all")

    def test_throws_exception_4(self):
        # Test throws AssertionError as both sides of range not included
        test_variant = "NM_004006.2:c.1_ACT[22]"
        with self.assertRaises(AssertionError):
            expanded_repeats.TandemRepeats.parse_repeat_variant(test_variant, "GRCh37", "all")

    def test_throws_exception_5(self):
        # Test throws AssertionError if repeat number between brackets not numeric
        test_variant = "LRG_199t1:c.20A[A]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                            test_variant, "GRCh37", "all")
        my_variant.check_transcript_type()
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        with self.assertRaises(AssertionError):
            my_variant.reformat()

    def test_empty_string(self):
        """
        Test for handling empty string.
        """
        variant_str = ""
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str, "GRCh37", "all")
        assert my_variant == False

    def test_simple_str_split(self):
        test_variant = "ENSG00000198947:g.1ACT[10]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(test_variant, "GRCh37", "all")
        my_variant.simple_split_string()
        assert my_variant.begining == "ENSG00000198947"
        assert my_variant.end == ":g.1ACT[10]"

    def test_transcript_versions(self):
        """
        Test for handing - instead of _ in variant_str.
        Previous code gave an error below:
        AssertionError: Unable to identify a colon (:) in the variant
        """
        variant_str = "NM_004006.2:c.13-14AC[7]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str, "GRCh37", "all")
        my_variant.check_transcript_type()
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat()
        assert formatted == "NM_004006.2:c.13_14insACACACACACACAC"
        assert my_variant.variant_str == "NM_004006.2:c.13-14AC[7]"
        assert my_variant.variant_position == "13_14"
        # checks correct position This is the most important assert.

        # Additional checks
        assert my_variant.ref_type == "RefSeq"
        # checks correct transcript type
        assert my_variant.reference == "NM_004006.2"
        # checks correct ref name
        assert my_variant.repeat_sequence == "AC"
        # checks repeat seq
        assert my_variant.copy_number == "7"
        # checks number of repeats is str and correct
        assert my_variant.after_the_bracket == ""
        # checks nothing is after the bracket


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

