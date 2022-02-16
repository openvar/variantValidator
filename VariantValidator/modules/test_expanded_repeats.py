"""
Tandem_Repeats_tests
Authors: Robert Wilson (@RSWilson1) and Rebecca Locke (@rklocke)
This code runs tests on the module tandem_repeats.py to check the outputs are as expected.

It checks known edge-case HGVS compliant variant strings.
Additional functionality to add:
- Check error handling
- Check correct errors for non-HGVS compliant strings.

"""
import unittest
import tandem_repeats


class TestExpandedRepeats(unittest.TestCase):
    """Tests for tandem_repeats.py to check the syntax checker returns
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
        my_variant = tandem_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str, "GRCh37", "all")
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
        my_variant = tandem_repeats.TandemRepeats.parse_repeat_variant(
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
        my_variant = tandem_repeats.TandemRepeats.parse_repeat_variant(
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


    def test_transcript_versions(self):
        """
        Test for handing transcript versions.
        Previous code gave an error below:
        AssertionError: Unable to identify a colon (:) in the variant
        """
        variant_str = "NM_004006.2:c.13AC[7]"
        my_variant = tandem_repeats.TandemRepeats.parse_repeat_variant(
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
            tandem_repeats.TandemRepeats.parse_repeat_variant(test_variant, "GRCh37","all")
    
    def test_throws_exception_2(self):
        # Test throws AssertionError if no repeat sequence is included
        test_variant = "ENST00000198947.1:c.1_2[10]"
        with self.assertRaises(AssertionError):
            tandem_repeats.TandemRepeats.parse_repeat_variant(test_variant, "GRCh37","all")

    def test_throws_exception_3(self):
        # Test throws AssertionError if allele variant
        test_variant = "LRG_199:g.[123456A>G];[345678G>C]"
        with self.assertRaises(AssertionError):
            tandem_repeats.TandemRepeats.parse_repeat_variant(test_variant, "GRCh37","all")

    def test_throws_exception_4(self):
        # Test throws AssertionError as both sides of range not included
        test_variant = "NM_004006.2:c.1_ACT[22]"
        with self.assertRaises(AssertionError):
            tandem_repeats.TandemRepeats.parse_repeat_variant(test_variant, "GRCh37","all")

    def test_throws_exception_5(self):
        # Test throws AssertionError if repeat number between brackets not numeric
        test_variant = "LRG_199t1:c.20A[A]"
        my_variant = tandem_repeats.TandemRepeats.parse_repeat_variant(
                            test_variant, "GRCh37", "all")
        my_variant.check_transcript_type()
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        with self.assertRaises(AssertionError):
            my_variant.reformat()


    # def test_basic_syntax_4(self):
    #     """
    #     Test for handling basic syntax of ENSG variant string.
    #     """
    #     variant_str = "LRG_199:g.[123456A>G];[345678G>C]"
    #     my_variant = tandem_repeats.TandemRepeats.parse_repeat_variant(
    #                                 variant_str, "GRCh37", "all")
    #     my_variant.check_transcript_type()
    #     my_variant.reformat_reference()
    #     my_variant.check_genomic_or_coding()
    #     formatted = my_variant.reformat()
    #     assert my_variant.variant_str == "LRG_199:g.[123456A>G];[345678G>C]"
    #     assert raise


    def test_empty_string(self):
        """
        Test for handling empty string.
        """
        variant_str = ""
        my_variant = tandem_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str, "GRCh37", "all")
        assert my_variant == False


if __name__ == "__main__":
    unittest.main()


# <LICENSE>
# Copyright (C) 2016-2022 VariantValidator Contributors
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
