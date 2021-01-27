"""
Exon_numbering_tests

Authors: Katie Williams (@kwi11iams) and Katherine Winfield (@kjwinfield)

This code runs tests on exon_numbering.py to check the outputs are as expected

"""
import unittest # This is the automatic test runner
from exon_numbering import finds_exon_number, request_sequence, check_variant


# Define some variants to test with
test_variant_1 = "NM_007294.3:c.1067A>G"
test_variant_2 = "NM_000088.3:c.642+1G>A"
test_variant_3 = "NM_000094.3:c.6751-3_6751-2del"
test_variant_4 = "NM_000088.3:c.589G>T"  # An exon boundary, maps to 715
test_variant_5 = "NM_000088.3:c.642del"  # This should map to 768
test_variant_6 = "NM_000094.3:c.6751-3_6753del"  # Intron start, exon end


class TestExonNumbering(unittest.TestCase):
    """
    Class TextExonNumbering automates running the tests, and reports failure if
    the output is not what is expected
    """
    def test_1(self):
        test_1 = finds_exon_number(test_variant_1)
        self.assertEqual(test_1['NC_000017.10']['start_exon'],
                         '14', "Failed test variant 1")

    def test_2(self):    
        test_2 = finds_exon_number(test_variant_2)
        self.assertEqual(test_2['NC_000017.10']['start_exon'],
                         '44i', "Failed test variant 2")

    def test_3(self):
        test_3 = finds_exon_number(test_variant_3)
        self.assertEqual(test_3['NC_000003.11']['start_exon'],
                         '32i', "Failed test variant 3")

    def test_4(self):
        test_4 = finds_exon_number(test_variant_4)
        self.assertEqual(test_4['NC_000017.10']['start_exon'],
                         '44', "Failed test variant 4")

    def test_5(self):
        test_5 = finds_exon_number(test_variant_5)
        self.assertEqual(test_5['NC_000017.10']['start_exon'],
                         '44', "Failed test variant 5")

    def test_6(self):
        test_6 = finds_exon_number(test_variant_6)
        self.assertEqual(test_6['NC_000003.11']['start_exon'],
                         '32i', "Failed test variant 6")


if __name__ == "__main__":
    unittest.main()
