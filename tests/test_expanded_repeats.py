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
from VariantValidator.modules.expanded_repeats import RepeatSyntaxError
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
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat(vv)
        assert str(formatted) == "NG_012232.1:g.3_6T[20]"
        assert my_variant.variant_str == "NG_012232.1:g.4T[20]"
        # checks correct transcript ref
        assert my_variant.reference == "NG_012232.1"
        # checks correct position
        assert str(my_variant.variant_position) == "3_6"
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
        # changed from previous version of "ENST00000263121.12:c.1082TCT[2]"
        # (pre-full exon handling) after verifying that coordinates matched by
        # testing "ENST00000263121.12:c.*62_*67delinsTCTTCT" transformed into
        # "ENST00000263121.12:c.*62_*67=" with base VV
        variant_str = "ENST00000263121.12:c.*62_*67TCT[2]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str,  "GRCh37", "all", vv)
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat(vv)
        assert str(formatted) == "ENST00000263121.12:c.*62_*67TCT[2]"
        assert my_variant.variant_str == "ENST00000263121.12:c.*62_*67TCT[2]"
        assert my_variant.prefix == "c"
        assert my_variant.reference == "ENST00000263121.12"
        # checks correct ref name
        assert str(my_variant.variant_position) == "*62_*67"
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
        assert str(my_variant.variant_position) == "1210-34"
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat(vv)
        assert str(formatted) == "NM_000492.4:c.1210-34_1210-13TG[11]"
        assert my_variant.variant_str == "NM_000492.4:c.1210-34TG[11]"
        assert my_variant.prefix == "c"
        assert my_variant.reference == "NM_000492.4"
        # checks correct ref name
        assert str(my_variant.variant_position) == "1210-34_1210-13"
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
        seq_range = expanded_repeats.TandemRepeats.get_range_from_single_or_start_pos(my_variant, vv)
        self.assertEqual(str(seq_range), "1289_1297")

    def test_empty_string(self):
        """
        Test for handling empty string.
        """
        variant_str = ""
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                                    variant_str, "GRCh37", "all", vv)
        assert my_variant == False

    def test_convert_tandem_fallback(self):
        """
        Test that fallback happens when input variant is a text hgvs variant not VV variant object
        """
        variant_str = "NM_003073.5:c.1085AGA[2]"
        my_variant_data = expanded_repeats.convert_tandem(variant_str, vv, 'GRCh37', 'all')
        self.assertEqual(str(my_variant_data["position"]), "1085_1093")
        self.assertEqual(str(my_variant_data["variant"]), "NM_003073.5:c.1085_1093AGA[2]")
        self.assertEqual(my_variant_data["copy_number"],'2')
        self.assertEqual(my_variant_data["repeat_sequence"],"AGA")
        self.assertEqual(my_variant_data["reference"],"NM_003073.5")
        self.assertEqual(my_variant_data["prefix"],"c")
        self.assertEqual(my_variant_data["reference_sequence_bases"],"AGAAGAAGA")

    def test_fail_tandem_unmatched_bracket(self):
        variant_str = "NM_003073.5:c.1085AGA[2"
        with self.assertRaises(RepeatSyntaxError) as catch:
            expanded_repeats.convert_tandem(variant_str, vv, 'GRCh37', 'all')
        self.assertTrue("variant in question is missing a matched bracket pair" in \
                str(catch.exception))

    def test_fail_tandem_no_repeat(self):
        variant_str = "NM_003073.5:c.1085[2]"
        with self.assertRaises(RepeatSyntaxError) as catch:
            expanded_repeats.convert_tandem(variant_str, vv, 'GRCh37', 'all')
        self.assertTrue(
                "Ensure that the repeated sequence is included" in \
                str(catch.exception))

    def test_fail_tandem_bad_repeat_seq(self):
        variant_str = "NM_003073.5:c.1085CXXX[2]"
        with self.assertRaises(RepeatSyntaxError) as catch:
            expanded_repeats.convert_tandem(variant_str, vv, 'GRCh37', 'all')
        self.assertTrue(
            "Please ensure the repeated sequence includes only Aa, Cc, Tt, Gg, Uu" in \
            str(catch.exception))

    def test_fail_tandem_bad_repeat_length(self):
        variant_str = "NM_003073.5:c.1085_1092AGA[2]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                variant_str,  "GRCh37", "all", vv)

        with self.assertRaises(RepeatSyntaxError) as catch:
            my_variant.check_positions_given(vv)
        self.assertTrue(
            "is not a multiple of the length of the provided repeat sequence" in \
            str(catch.exception))

    def test_fail_tandem_ref_nomatch(self):
        # need a span to trigger this test
        # also needs a test for no alt at test time [0] code path
        variant_str = "NM_000492.4:c.1210-35_1210-12TG[11]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                variant_str,  "GRCh37", "all", vv)
        my_variant.get_valid_n_or_g_range_from_input_pos(vv)
        with self.assertRaises(RepeatSyntaxError) as catch:
            my_variant.check_positions_given(vv)
        self.assertTrue(
            "The repeat sequence does not match the reference sequence " in \
            str(catch.exception))
        variant_str = "NM_003073.5:c.1085_1093AGA[0]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                variant_str,  "GRCh37", "all", vv)
        my_variant.get_valid_n_or_g_range_from_input_pos(vv)
        my_variant.check_positions_given(vv)

    def test_fail_reformat_if_extra_after_the_bracket(self):
        # self.after_the_bracket needs to be true when bad variant passed
        variant_str = "NG_012232.1:g.4T[20]wrong"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                variant_str,  "GRCh37", "all", vv)
        with self.assertRaises(RepeatSyntaxError) as catch:
            formatted = my_variant.reformat(vv)
        self.assertTrue(
            "No information should be included after the number of repeat units." in \
            str(catch.exception))

    def test_fail_reformat_if_bad_repeat_characters(self):
        variant_str = "NG_012232.1:g.4CC[20]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                variant_str,  "GRCh37", "all", vv)
        with self.assertRaises(RepeatSyntaxError) as catch:
            my_variant.repeat_sequence = 'EE'
            formatted = my_variant.reformat(vv)
        self.assertTrue(
            "Please ensure the repeated sequence includes only Aa, Cc, Tt, Gg, Uu" in \
            str(catch.exception))

    def test_fail_reformat_if_not_copy_number_is_decimal(self):
        #switch to is int?
        variant_str = "NG_012232.1:g.4T[aaa]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                variant_str,  "GRCh37", "all", vv)
        with self.assertRaises(RepeatSyntaxError) as catch:
            formatted = my_variant.reformat(vv)
        self.assertTrue(
            "The number of repeat units included between square brackets must be numeric" in \
            str(catch.exception))


    #def test_pass_tandem_ref_nomatch_no_length(self):
    #    # No length for genomic match, and fail caught, should only happen for broken input
    #    # which should be caught already, or coordinates beyond the end of the seq.
    #    # This could not trigger and was removed, we may want to add a ^ / $ to the regex
    #    # however, and if so would then want to re-add this test.
    #    variant_str = "NG_012232.1:g.90000000CC[20]"
    #    my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
    #            variant_str,  "GRCh37", "all", vv)
    #    my_variant.check_positions_given(vv)

    def test_check_exon_boundaries_no_bounds_used(self):
        variant_str = "NM_003073.5:c.1085AGA[2]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                variant_str, vv, 'GRCh37', 'all')
        my_variant.check_exon_boundaries(vv)

    def test_exon_c_to_n_safe_on_n(self):
        # convert_c_to_n_coordinates should be safe to call on n/g input
        variant_str = "NG_012232.1:g.4T[5]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
                variant_str,  "GRCh37", "all", vv)
        pos = str(my_variant.variant_position)
        my_variant.convert_c_to_n_coordinates()
        assert pos == str(my_variant.variant_position)

    def test_bad_coppy_number_after_genomic_map(self):
        # test in get_range_from_single_or_start_pos
        # during remap of BaseOffsetInterval type pos with at least one .offset
        # i.e. intronic n/c
        # when number of regex matches within range of genomic span != expected
        variant_str = "NM_000492.4:c.1210-34_1210-11TG[11]"
        variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
            variant_str,'GRCh37', 'all',vv)

        with self.assertRaises(RepeatSyntaxError) as catch:
            variant.reformat(vv)
        self.assertTrue(
            "The repeat sequence does not match the expected copy number at position" in\
            str(catch.exception))

class TestCVariantsExpanded(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()
        cls.vv.testing = True

    def test_coridnates_over_cds_start(self):
        variant = 'NM_004006.2:c.-3_1A[4]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert 'NM_004006.2:c.-3_1A[4]' in results

    def test_exon_boundary_single_position(self):
        variant = 'NM_004006.2:c.13-14AC[7]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert ("ExonBoundaryError: Position 13-14 does not correspond with an exon boundary for transcript "
                "NM_004006.2") in results["validation_warning_1"]["validation_warnings"]

    def test_exon_boundary_single_position_plus(self):
        variant = 'NM_004006.2:c.13+14AC[7]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert ("ExonBoundaryError: Position 13+14 does not correspond with an exon boundary for transcript "
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
        assert results["NM_000492.4:c.1210-34_1210-13TG[11]"][
            "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000007.13:g.117188661_117188682TG[11]"
        assert results["NM_000492.4:c.1210-34_1210-13TG[11]"][
            "validation_warnings"] ==  [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NM_000492.4:c.1210-34TG[11]. The corrected description is NM_000492.4:c.1210-34_1210-13TG[11]",
            "ExpandedRepeatWarning: NM_000492.4:c.1210-34_1210-13TG[11] should only be used as an annotation for the core HGVS descriptions provided",
        ]

    def test_intronic_range(self):
        variant = 'NM_000492.4:c.1210-34_1210-13TG[11]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_000492.4:c.1210-34_1210-13TG[11]"][
            "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000007.13:g.117188661_117188682TG[11]"
        assert results["NM_000492.4:c.1210-34_1210-13TG[11]"][
            "validation_warnings"] == [
            "ExpandedRepeatWarning: NM_000492.4:c.1210-34_1210-13TG[11] should only be used as an annotation for the "
            "core HGVS descriptions provided",
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
        assert results["NM_003073.5:c.1085_1093AGA[3]"][
            "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000022.10:g.24175857_24175865AGA[3]"
        assert results["NM_003073.5:c.1085_1093AGA[3]"][
            "validation_warnings"] ==  [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NM_003073.5:c.1085AGA[3]. The corrected description is NM_003073.5:c.1085_1093AGA[3]",
            "ExpandedRepeatWarning: NM_003073.5:c.1085_1093AGA[3] should only be used as an annotation for the core HGVS descriptions provided"
        ]

    def test_exonic_range(self):
        variant = 'NM_003073.5:c.1085_1093AGA[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_003073.5:c.1085_1093AGA[3]"][
            "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000022.10:g.24175857_24175865AGA[3]"
        assert results["NM_003073.5:c.1085_1093AGA[3]"][
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
        assert results["NM_002024.5:c.-129_-100CGG[10]"][
            "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000023.10:g.146993569_146993598CGG[10]"
        assert results["NM_002024.5:c.-129_-100CGG[10]"][
            "validation_warnings"] == [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NM_002024.5:c.-129CGG[10]. The corrected description is NM_002024.5:c.-129_-100CGG[10]",
            "ExpandedRepeatWarning: NM_002024.5:c.-129_-100CGG[10] should only be used as an annotation for the core HGVS descriptions provided",
            "TranscriptVersionWarning: A more recent version of the selected reference sequence NM_002024.5 is available for genome build GRCh37 (NM_002024.6)"
        ]

    def test_5_utr_range(self):
        variant = 'NM_002024.5:c.-129_-100CGG[10]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_002024.5:c.-129_-100CGG[10]"][
            "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000023.10:g.146993569_146993598CGG[10]"
        assert results["NM_002024.5:c.-129_-100CGG[10]"][
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
        due to the reference centric nature of our current tooling but is now tested
        later in the LRG section of the tests."""
        variant = 'NM_002111.8:c.54GCA[21]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_002111.8:c.54_116GCA[21]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000004.11:g.3076606_3076662GCA[21]"
        assert results["NM_002111.8:c.54_116GCA[21]"]["validation_warnings"] ==  [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NM_002111.8:c.54GCA[21]. The corrected description is NM_002111.8:c.54_116GCA[21]",
            "ExpandedRepeatWarning: NM_002111.8:c.54_116GCA[21] should only be used as an annotation for the core HGVS descriptions provided"
        ]

    def test_antisense_intron_range(self):
        variant = 'NM_000088.3:c.589-1_590G[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_000088.3:c.589-1_590G[3]"][
            "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000017.10:g.48275362_48275364C[3]"
        assert results["NM_000088.3:c.589-1_590G[3]"][
            "validation_warnings"] == [
            "ExpandedRepeatWarning: NM_000088.3:c.589-1_590G[3] should only be used as an annotation for the "
            "core HGVS descriptions provided",
            "TranscriptVersionWarning: A more recent version of the selected reference sequence NM_000088.3 "
            "is available for genome build GRCh37 (NM_000088.4)"
        ]

    def test_antisense_intron_single_pos(self):
        variant = 'NM_000088.3:c.589-1G[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_000088.3:c.589-1_590G[3]"][
                   "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000017.10:g.48275362_48275364C[3]"
        assert results["NM_000088.3:c.589-1_590G[3]"][
                   "validation_warnings"] == [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NM_000088.3:c.589-1G[3]. The corrected description is NM_000088.3:c.589-1_590G[3]",
            "ExpandedRepeatWarning: NM_000088.3:c.589-1_590G[3] should only be used as an annotation for the core HGVS descriptions provided",
            "TranscriptVersionWarning: A more recent version of the selected reference sequence NM_000088.3 is available for genome build GRCh37 (NM_000088.4)"
        ]

    def test_antisense_intron_single_pos_2(self):
        variant = 'NM_000088.3:c.589-18T[5]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_000088.3:c.589-18_589-14T[5]"][
                   "primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000017.10:g.48275377_48275381A[5]"
        assert results["NM_000088.3:c.589-18_589-14T[5]"][
                   "validation_warnings"] == [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NM_000088.3:c.589-18T[5]. The corrected description is NM_000088.3:c.589-18_589-14T[5]",
            "ExpandedRepeatWarning: NM_000088.3:c.589-18_589-14T[5] should only be used as an annotation for the core HGVS descriptions provided",
            "TranscriptVersionWarning: A more recent version of the selected reference sequence NM_000088.3 is available for genome build GRCh37 (NM_000088.4)"
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

class TestTranscriptVariantsExpandedNvsC(TestCase):
    """
    The same tests as the C variant type specific tests above, but with genes
    that have a pair of coding and non coding transcript versions. Tests should
    pass with the same errors and/or genomic mappings for both, since target
    transcripts where selected so that each have the same first and second
    exon coordinates at least.
    # random sample pair 1  HGNC:15517 NC_000017.10 NM_022167.4 NR_110010.2
    NM_022167.4 cds 15 to 2613, +1 strand, start 48423486 end 48438546
    NR_110010.2 +1 strand start 48423486 end 48438546
    # status identical to genome ':n.' coordinates 25_31GC[3] 3_6C[3]
    random sample pair 2, -1 strand NC_000010.10 NM_001350922.2 NR_146939.2
    NM_001350922.2 cds 213 to 2313
    """
    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()
        cls.vv.testing = True

    def test_coridnates_over_cds_start_c(self):
        # It is important to test the -1 position specifically as previous versions
        # had some issues with -1 (but not -2, -3 etc.)
        variant = 'NM_022167.4:c.-1_1GA[1]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert 'NM_022167.4:c.-1_1GA[1]' in results

    def test_coridnates_over_pseudo_cds_start_n(self):
        variant = 'NR_110010.2:n.15_16GA[1]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert 'NR_110010.2:n.15_16GA[1]' in results

    def test_cordinateds_over_cds_end_c(self):
        # For this NM_001160367.2 and the following NR_027702.2 test there is a 35bp deletion
        # that breaks the translation start for the NR version of the transcript. Finding good
        # transcripts with identical sequences all through to past the CDS where one of the
        # pair was NR and the other NM and coding proved impractical, so a otherwise clean
        # start deletion was found.
        variant = 'NM_001160367.2:c.870_*1AC[1]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert 'NM_001160367.2:c.870_*1AC[1]' in results

    def test_cordinateds_over_pseudo_cds_end_n(self):
        # As mentioned n coordinates here should be -35 WRT to n equivalent coordinates for above
        variant = 'NR_027702.2:n.1032_1033AC[1]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert 'NR_027702.2:n.1032_1033AC[1]' in results

    def test_exon_boundary_single_position_c(self):
        variant = 'NM_022167.4:c.21-4GC[5]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert ("ExonBoundaryError: Position 21-4 does not correspond with an exon boundary for "
            "transcript NM_022167.4") in results["validation_warning_1"]["validation_warnings"]

    def test_exon_boundary_single_position_n(self):
        variant = 'NR_110010.2:n.36-4GC[5]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert ("ExonBoundaryError: Position 36-4 does not correspond with an exon boundary for "
            "transcript NR_110010.2") in results["validation_warning_1"]["validation_warnings"]

    def test_exon_boundary_range_c(self):
        variant = 'NM_022167.4:c.10_20-4GC[5]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert ("ExonBoundaryError: Position 10_20-4 does not correspond with an exon boundary for "
                "transcript NM_022167.4") in results["validation_warning_1"]["validation_warnings"]

    def test_exon_boundary_range_n(self):
        variant = 'NR_110010.2:n.25_35-4GC[5]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert ("ExonBoundaryError: Position 25_35-4 does not correspond with an exon boundary for "
                "transcript NR_110010.2") in results["validation_warning_1"]["validation_warnings"]

    def test_intronic_single_position_c(self):
        variant = 'NM_022167.4:c.135+1G[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_022167.4:c.135_135+1G[2]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48423636_48423637G[2]"
        assert results["NM_022167.4:c.135_135+1G[2]"][
            "validation_warnings"] ==  [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NM_022167.4:c.135+1G[2]. The corrected description is NM_022167.4:c.135_135+1G[2]",
            "ExpandedRepeatWarning: NM_022167.4:c.135_135+1G[2] should only be used as an annotation for the core HGVS descriptions provided",
        ]

    def test_intronic_single_position_n(self):
        variant = 'NR_110010.2:n.150+1G[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NR_110010.2:n.150_150+1G[2]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48423636_48423637G[2]"
        assert results["NR_110010.2:n.150_150+1G[2]"][
            "validation_warnings"] == [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NR_110010.2:n.150+1G[2]. The corrected description is NR_110010.2:n.150_150+1G[2]",
            "ExpandedRepeatWarning: NR_110010.2:n.150_150+1G[2] should only be used as an annotation for the core HGVS descriptions provided",
        ]

    def test_intronic_range_c(self):
        variant = 'NM_022167.4:c.135+6_135+7C[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_022167.4:c.135+6_135+7C[2]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48423642_48423643C[2]"
        assert results["NM_022167.4:c.135+6_135+7C[2]"][
            "validation_warnings"] == [
            "ExpandedRepeatWarning: NM_022167.4:c.135+6_135+7C[2] should only be used as an "
            "annotation for the core HGVS descriptions provided",
        ]

    def test_intronic_range_n(self):
        variant = 'NR_110010.2:n.150+6_150+7C[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NR_110010.2:n.150+6_150+7C[2]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48423642_48423643C[2]"
        assert results["NR_110010.2:n.150+6_150+7C[2]"][
            "validation_warnings"] == [
            "ExpandedRepeatWarning: NR_110010.2:n.150+6_150+7C[2] should only be used as an "
            "annotation for the core HGVS descriptions provided",
        ]

    def test_incorrect_intronic_range_c(self):
        variant = 'NM_022167.4:c.135+5_135+6C[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] == [
            "RepeatSyntaxError: The repeat sequence does not match the reference sequence at the "
            "given position 135+5_135+6, expected CC but the reference is TC at the specified position"
        ]

    def test_incorrect_intronic_range_n(self):
        variant = 'NR_110010.2:n.150+5_150+6C[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] == [
            "RepeatSyntaxError: The repeat sequence does not match the reference sequence at the "
            "given position 150+5_150+6, expected CC but the reference is TC at the specified position"
        ]

    def test_exonic_single_position_c(self):
        variant = 'NM_022167.4:c.11GC[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_022167.4:c.11_16GC[3]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48423512_48423517GC[3]"
        assert results["NM_022167.4:c.11_16GC[3]"]["validation_warnings"] ==  [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NM_022167.4:c.11GC[3]. The corrected description is NM_022167.4:c.11_16GC[3]",
            "ExpandedRepeatWarning: NM_022167.4:c.11_16GC[3] should only be used as an annotation for the core HGVS descriptions provided"
        ]

    def test_exonic_single_position_n(self):
        variant = 'NR_110010.2:n.26GC[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NR_110010.2:n.26_31GC[3]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48423512_48423517GC[3]"
        assert results["NR_110010.2:n.26_31GC[3]"]["validation_warnings"] ==  [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NR_110010.2:n.26GC[3]. The corrected description is NR_110010.2:n.26_31GC[3]",
            "ExpandedRepeatWarning: NR_110010.2:n.26_31GC[3] should only be used as an annotation for the core HGVS descriptions provided"
        ]

    def test_exonic_range_c(self):
        variant = 'NM_022167.4:c.11_16GC[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_022167.4:c.11_16GC[3]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48423512_48423517GC[3]"
        assert results["NM_022167.4:c.11_16GC[3]"]["validation_warnings"] ==  [
            "ExpandedRepeatWarning: NM_022167.4:c.11_16GC[3] should only be used as an annotation for the core HGVS descriptions provided"
        ]

    def test_exonic_range_n(self):
        variant = 'NR_110010.2:n.26_31GC[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NR_110010.2:n.26_31GC[3]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48423512_48423517GC[3]"
        assert results["NR_110010.2:n.26_31GC[3]"]["validation_warnings"] == [
            "ExpandedRepeatWarning: NR_110010.2:n.26_31GC[3] should only be used as an annotation "
            "for the core HGVS descriptions provided",
        ]

    def test_incorrect_exonic_range_c(self):
        variant = 'NM_022167.4:c.10_15GC[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] == [
            'RepeatSyntaxError: The provided repeat sequence GC does not match the reference '
            'sequence AG at the given position 25_26 of reference sequence NM_022167.4'
        ]

    def test_incorrect_exonic_range_n(self):
        variant = 'NR_110010.2:n.25_30GC[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] == [
            'RepeatSyntaxError: The provided repeat sequence GC does not match the reference '
            'sequence AG at the given position 25_26 of reference sequence NR_110010.2'
        ]

    def test_5_utr_single_pos_c(self):
        variant = 'NM_022167.4:c.-12C[3]'
        tr_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
            variant, "GRCh37", "all", vv)
        tr_variant.check_genomic_or_coding()
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_022167.4:c.-13_-11C[3]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48423489_48423491C[3]"
        assert results["NM_022167.4:c.-13_-11C[3]"]["validation_warnings"] ==  [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NM_022167.4:c.-12C[3]. The corrected description is NM_022167.4:c.-13_-11C[3]",
            "ExpandedRepeatWarning: NM_022167.4:c.-13_-11C[3] should only be used as an annotation for the core HGVS descriptions provided"
        ]

    def test_5_utr_single_pos_n(self):
        variant = 'NR_110010.2:n.3C[3]'
        tr_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
            variant, "GRCh37", "all", vv)
        tr_variant.check_genomic_or_coding()
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NR_110010.2:n.3_5C[3]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48423489_48423491C[3]"
        assert results["NR_110010.2:n.3_5C[3]"]["validation_warnings"] == [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NR_110010.2:n.3C[3]. The corrected description is NR_110010.2:n.3_5C[3]",
            "ExpandedRepeatWarning: NR_110010.2:n.3_5C[3] should only be used as an annotation for the core HGVS descriptions provided"
        ]

    def test_5_utr_range_c(self):
        variant = 'NM_022167.4:c.-13_-11C[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_022167.4:c.-13_-11C[3]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48423489_48423491C[3]"
        assert results["NM_022167.4:c.-13_-11C[3]"]["validation_warnings"] == [
            "ExpandedRepeatWarning: NM_022167.4:c.-13_-11C[3] should only be used as an annotation "
            "for the core HGVS descriptions provided",
        ]

    def test_5_utr_range_n(self):
        variant = 'NR_110010.2:n.3_5C[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NR_110010.2:n.3_5C[3]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48423489_48423491C[3]"
        assert results["NR_110010.2:n.3_5C[3]"]["validation_warnings"] == [
            "ExpandedRepeatWarning: NR_110010.2:n.3_5C[3] should only be used as an annotation for "
            "the core HGVS descriptions provided",
        ]

    def test_5_utr_single_pos_incorrect_c(self):
        variant = 'NM_022167.4:c.-12T[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] == [
            'RepeatSyntaxError: The provided repeat sequence T does not match the reference '
            'sequence C at the given position 2_2 of reference sequence NM_022167.4'
        ]

    def test_5_utr_single_pos_incorrect_c(self):
        variant = 'NM_022167.4:c.-13T[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] ==[
            'RepeatSyntaxError: The provided repeat sequence T does not match the reference '
            'sequence C at the given position 3_3 of reference sequence NM_022167.4'
        ]

    #we may want to add an additional test_gap_crossing pair here later

    def test_antisense_intron_range_c(self):
        variant = 'NM_001350922.2:c.240+14_240+21TGGG[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_001350922.2:c.240+14_240+21TGGG[2]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000010.10:g.128358789_128358796CCCA[2]"
        assert results["NM_001350922.2:c.240+14_240+21TGGG[2]"]["validation_warnings"] == [
            "ExpandedRepeatWarning: NM_001350922.2:c.240+14_240+21TGGG[2] should only be used as "
            "an annotation for the core HGVS descriptions provided",
        ]

    def test_antisense_intron_range_n(self):
        variant = 'NR_146939.2:n.453+14_453+21TGGG[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NR_146939.2:n.453+14_453+21TGGG[2]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000010.10:g.128358789_128358796CCCA[2]"
        assert results["NR_146939.2:n.453+14_453+21TGGG[2]"]["validation_warnings"] == [
            "ExpandedRepeatWarning: NR_146939.2:n.453+14_453+21TGGG[2] should only be used as an "
            "annotation for the core HGVS descriptions provided",
        ]

    def test_antisense_intron_single_pos_c(self):
        variant = 'NM_001350922.2:c.240+14TGGG[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_001350922.2:c.240+14_240+21TGGG[2]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000010.10:g.128358789_128358796CCCA[2]"
        assert results["NM_001350922.2:c.240+14_240+21TGGG[2]"]["validation_warnings"] == [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NM_001350922.2:c.240+14TGGG[2]. The corrected description is NM_001350922.2:c.240+14_240+21TGGG[2]",
            "ExpandedRepeatWarning: NM_001350922.2:c.240+14_240+21TGGG[2] should only be used as an annotation for the core HGVS descriptions provided",
        ]

    def test_antisense_intron_single_pos_n(self):
        variant = 'NR_146939.2:n.453+14TGGG[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NR_146939.2:n.453+14_453+21TGGG[2]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000010.10:g.128358789_128358796CCCA[2]"
        assert results["NR_146939.2:n.453+14_453+21TGGG[2]"]["validation_warnings"] == [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NR_146939.2:n.453+14TGGG[2]. The corrected description is NR_146939.2:n.453+14_453+21TGGG[2]",
            "ExpandedRepeatWarning: NR_146939.2:n.453+14_453+21TGGG[2] should only be used as an annotation for the core HGVS descriptions provided",
        ]

    def test_antisense_intron_single_pos_2_c(self):
        variant = 'NM_001350922.2:c.241-119TC[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_001350922.2:c.241-121_241-118TC[2]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000010.10:g.128335324_128335327GA[2]"
        assert results["NM_001350922.2:c.241-121_241-118TC[2]"]["validation_warnings"] ==  [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NM_001350922.2:c.241-119TC[2]. The corrected description is NM_001350922.2:c.241-121_241-118TC[2]",
            "ExpandedRepeatWarning: NM_001350922.2:c.241-121_241-118TC[2] should only be used as an annotation for the core HGVS descriptions provided",
        ]

    def test_antisense_intron_single_pos_2_n(self):
        variant = 'NR_146939.2:n.454-119TC[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NR_146939.2:n.454-121_454-118TC[2]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000010.10:g.128335324_128335327GA[2]"
        assert results["NR_146939.2:n.454-121_454-118TC[2]"]["validation_warnings"] == [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NR_146939.2:n.454-119TC[2]. The corrected description is NR_146939.2:n.454-121_454-118TC[2]",
            "ExpandedRepeatWarning: NR_146939.2:n.454-121_454-118TC[2] should only be used as an annotation for the core HGVS descriptions provided",
        ]

    def test_antisense_intron_single_pos_seq_inverted_c(self):
        variant = 'NM_001350922.2:c.241-119GA[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] == [
            'RepeatSyntaxError: The provided repeat sequence TC does not match the reference sequence'
            ' GA at the given position 128335324_128335325 of reference sequence NC_000010.10'
        ]

    def test_antisense_intron_single_pos_seq_inverted_n(self):
        variant = 'NR_146939.2:n.454-119GA[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] == [
            'RepeatSyntaxError: The provided repeat sequence TC does not match the reference sequence'
            ' GA at the given position 128335324_128335325 of reference sequence NC_000010.10'
        ]


class TestExpandedRepeatGenomic(TestCase):
    """
    Tests for genomic reference targeting expanded repeats, this
    test set only tests the genomic code.
    """
    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()
        cls.vv.testing = True

    def test_basic_syntax_NC(self):
        """
        Test for handling basic syntax of variant string, with NC_ genomic
        reference type variant input directly via TandemRepeats code
        Should pass if all basic syntax checks passed.
        """
        variant_str = 'NC_000023.10:g.33362721A[20]'
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
            variant_str,  "GRCh37", "all", vv)
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat(vv)
        assert str(formatted) == "NC_000023.10:g.33362721_33362724A[20]"
        assert my_variant.variant_str == "NC_000023.10:g.33362721A[20]"
        # check, in order ref ID, variant_pos, seq, copy number
        # and post-var content (should be none)
        assert my_variant.reference == "NC_000023.10"
        assert str(my_variant.variant_position) == "33362721_33362724"
        assert my_variant.repeat_sequence == "A"
        assert my_variant.copy_number == "20"
        assert my_variant.after_the_bracket == ""

    def test_NC_genomic_range(self):
        """
        Test that a working genomic input range produces an equal and valid output
        """
        variant = 'NC_000022.10:g.24175857_24175865AGA[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_003073.5:c.1085_1093AGA[3]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000022.10:g.24175857_24175865AGA[3]"
        assert results["NM_003073.5:c.1085_1093AGA[3]"]["validation_warnings"] == [
            "ExpandedRepeatWarning: NC_000022.10:g.24175857_24175865AGA[3] should only be used as "
            "an annotation for the core HGVS descriptions provided"
        ]

    def test_NC_genomic_single_position_to_span(self):
        """
        Test mapping of single location starting tandem repeat to span via full VV front end
        """
        variant = 'NC_000023.10:g.33362721A[20]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["intergenic_variant_1"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000023.10:g.33362721_33362724A[20]"
        assert "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NC_000023.10:g.33362721A[20]. The corrected description is NC_000023.10:g.33362721_33362724A[20]" in \
                results['intergenic_variant_1']["validation_warnings"]

    def test_incorrect_NC_genomic_range(self):
        """
        Test that an invalid range gives an appropriate error response
        """
        variant = 'NC_000022.10:g.24175854_24175865AGA[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] == [
            'RepeatSyntaxError: The provided repeat sequence AGA does not match the reference '
            'sequence TGG at the given position 24175854_24175856 of reference sequence NC_000022.10'
        ]

class TestExpandedRepeatGenomicToTranscript(TestCase):
    """
    Tests using a genomic reference input, the tests are the round trip
    replies of the queries in the C type transcript mapping code.
    Test positions deliberately without a proper mapping to the genome are
    excluded.
    """
    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()
        cls.vv.testing = True

    def test_intronic_single_position(self):
        """Reverse of test for 'NM_000492.4:c.1210-34TG[11]'"""
        variant = 'NC_000007.13:g.117188661TG[11]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_000492.4:c.1210-34_1210-13TG[11]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000007.13:g.117188661_117188682TG[11]"
        assert results["NM_000492.4:c.1210-34_1210-13TG[11]"]["validation_warnings"] == [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NC_000007.13:g.117188661TG[11]. The corrected description is NC_000007.13:g.117188661_117188682TG[11]",
            "ExpandedRepeatWarning: NC_000007.13:g.117188661_117188682TG[11] should only be used as an annotation for the core HGVS descriptions provided"
        ]

    def test_intronic_range(self):
        """Reverse of test for'NM_000492.4:c.1210-34_1210-13TG[11]'"""
        variant = 'NC_000007.13:g.117188661_117188682TG[11]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_000492.4:c.1210-34_1210-13TG[11]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000007.13:g.117188661_117188682TG[11]"
        assert results["NM_000492.4:c.1210-34_1210-13TG[11]"]["validation_warnings"] == [
            "ExpandedRepeatWarning: NC_000007.13:g.117188661_117188682TG[11] should only be used "
            "as an annotation for the core HGVS descriptions provided",
        ]

    def test_exonic_single_position(self):
        """Reverse of test for 'NM_003073.5:c.1085AGA[3]'"""
        variant = 'NC_000022.10:g.24175857AGA[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_003073.5:c.1085_1093AGA[3]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000022.10:g.24175857_24175865AGA[3]"
        assert results["NM_003073.5:c.1085_1093AGA[3]"]["validation_warnings"] ==  [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NC_000022.10:g.24175857AGA[3]. The corrected description is NC_000022.10:g.24175857_24175865AGA[3]",
            "ExpandedRepeatWarning: NC_000022.10:g.24175857_24175865AGA[3] should only be used as an annotation for the core HGVS descriptions provided"
        ]

    def test_exonic_range(self):
        """Reverse of test for 'NM_003073.5:c.1085_1093AGA[3]'"""
        variant = 'NC_000022.10:g.24175857_24175865AGA[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_003073.5:c.1085_1093AGA[3]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000022.10:g.24175857_24175865AGA[3]"
        assert results["NM_003073.5:c.1085_1093AGA[3]"]["validation_warnings"] == [
            "ExpandedRepeatWarning: NC_000022.10:g.24175857_24175865AGA[3] should only be used as "
            "an annotation for the core HGVS descriptions provided"
        ]

    def test_5_utr_single_pos(self):
        """Reverse of test for 'NM_002024.5:c.-129CGG[10]'"""
        variant = 'NC_000023.10:g.146993569CGG[10]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_002024.5:c.-129_-100CGG[10]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000023.10:g.146993569_146993598CGG[10]"

    def test_5_utr_range(self):
        """Reverse of test for 'NM_002024.5:c.-129_-100CGG[10]'"""
        variant = 'NC_000023.10:g.146993569_146993598CGG[10]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_002024.5:c.-129_-100CGG[10]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000023.10:g.146993569_146993598CGG[10]"
        assert results["NM_002024.5:c.-129_-100CGG[10]"]["validation_warnings"] == [
            "ExpandedRepeatWarning: NC_000023.10:g.146993569_146993598CGG[10] should only be used "
            "as an annotation for the core HGVS descriptions provided",
            "TranscriptVersionWarning: A more recent version of the selected reference sequence "
            "NM_002024.5 is available for genome build GRCh37 (NM_002024.6)"
        ]

    def test_gap_crossing(self):
        """Reverse of test for 'NM_002111.8:c.54_110GCA[21]'. Given the
        nature of the gap in the alignment between this transcript and
        the genome this represents a bad alignment case, and may change
        due to improvements later. Check that funny stuff does not happen
        without us knowing either way. """
        variant = 'NC_000004.11:g.3076606_3076662GCA[21]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        print(results.keys())
        assert results["NM_001388492.1:c.54_110GCA[21]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000004.11:g.3076606_3076662GCA[21]"
        assert (("ExpandedRepeatWarning: NC_000004.11:g.3076606_3076662GCA[21] should only be used as an annotation "
                "for the core HGVS descriptions provided") in
                results["NM_001388492.1:c.54_110GCA[21]"]["validation_warnings"])

    def test_antisense_intron_range(self):
        """Reverse of test for 'NM_000088.3:c.589-1_590G[3]'"""
        variant = 'NC_000017.10:g.48275362_48275364C[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_000088.3:c.589-1_590G[3]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48275362_48275364C[3]"
        assert results["NM_000088.3:c.589-1_590G[3]"]["validation_warnings"] == [
            "ExpandedRepeatWarning: NC_000017.10:g.48275362_48275364C[3] should only be used as an"
            " annotation for the core HGVS descriptions provided",
            "TranscriptVersionWarning: A more recent version of the selected reference sequence "
            "NM_000088.3 is available for genome build GRCh37 (NM_000088.4)",
        ]

    def test_antisense_intron_single_pos(self):
        """Reverse of test for 'NM_000088.3:c.589-1G[3]'"""
        variant = 'NC_000017.10:g.48275362C[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_000088.3:c.589-1_590G[3]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48275362_48275364C[3]"

    def test_antisense_intron_single_pos_2(self):
        """Reverse of test for 'NM_000088.3:c.589-18T[5]'"""
        variant = 'NC_000017.10:g.48275377A[5]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_000088.3:c.589-18_589-14T[5]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48275377_48275381A[5]"

class TestExpandedRepeatRefSeqGenomic(TestCase):
    """
    Tests for RefSeqGene(RSG) reference targeting expanded repeats
    Same as basic genomic mapping tests + a few known RSG transcripts
    to verify that alignments are found correctly.
    Large parts of the underlying code paths are shared, all with
    transcript mappings beyond a certain point, so more extensive RSG<->tx
    mapping tests are left out for now, until we have problematic user
    input to add.
    """
    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()
        cls.vv.testing = True

    def test_basic_syntax_RSG(self):
        """
        Test for handling basic syntax of variant string, with RSG_ genomic
        reference type variant input directly via TandemRepeats code
        Should pass if all basic syntax checks passed.
        """
        variant_str = "NG_012232.1:g.4T[20]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
            variant_str,  "GRCh37", "all", vv)
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat(vv)
        assert str(formatted) == "NG_012232.1:g.3_6T[20]"
        assert my_variant.variant_str == "NG_012232.1:g.4T[20]"
        # checks correct transcript ref
        assert my_variant.reference == "NG_012232.1"
        # checks correct position
        assert str(my_variant.variant_position) == "3_6"
        # checks repeat seq
        assert my_variant.repeat_sequence == "T"
        # checks correct suffix
        assert my_variant.copy_number == "20"
        # checks number of repeats is str and correct
        assert my_variant.after_the_bracket == ""
        # checks nothing is after the bracket

    def test_RSG_genomic_range(self):
        """
        Test that a working input range produces an equal and valid output
        """
        variant = 'NG_012232.1:g.3_6T[20]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["intergenic_variant_1"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000023.10:g.33362724_33362725insAAAAAAAAAAAAAAAA"
        assert results["intergenic_variant_1"]['hgvs_refseqgene_variant']\
            == 'NG_012232.1:g.6_7insTTTTTTTTTTTTTTTT'
        assert results["intergenic_variant_1"]["validation_warnings"] == [
            "ExpandedRepeatWarning: NG_012232.1:g.3_6T[20] should only be used as an annotation for"
            " the core HGVS descriptions provided",
            'NG_012232.1:g.6_7insTTTTTTTTTTTTTTTT automapped to genome position '
            'NC_000023.10:g.33362724_33362725insAAAAAAAAAAAAAAAA',
            'No individual transcripts have been identified that fully overlap the described variation in the genomic sequence. Large variants might span one or more genes and are currently only described at the genome (g.) level.'
        ]

    def test_RSG_genomic_single_position_to_span(self):
        """
        Test mapping of single location starting tandem repeat to span via full VV front end
        """
        variant = "NG_012232.1:g.4T[20]"
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["intergenic_variant_1"]["primary_assembly_loci"]["grch37"][
                "hgvs_genomic_description"] == "NC_000023.10:g.33362724_33362725insAAAAAAAAAAAAAAAA"
        assert results['intergenic_variant_1']["validation_warnings"] ==  [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NG_012232.1:g.4T[20]. The corrected description is NG_012232.1:g.3_6T[20]",
            "ExpandedRepeatWarning: NG_012232.1:g.3_6T[20] should only be used as an annotation for the core HGVS descriptions provided",
            "NG_012232.1:g.6_7insTTTTTTTTTTTTTTTT automapped to genome position NC_000023.10:g.33362724_33362725insAAAAAAAAAAAAAAAA",
            'No individual transcripts have been identified that fully overlap the described variation in the genomic sequence. Large variants might span one or more genes and are currently only described at the genome (g.) level.'
        ]

    def test_incorrect_RSG_genomic_range(self):
        """
        Test that an invalid range gives an appropriate error response
        """
        variant = 'NG_012232.1:g.3_6C[20]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] == [
            'RepeatSyntaxError: The provided repeat sequence C does not match the reference '
            'sequence T at the given position 3_3 of reference sequence NG_012232.1'
        ]

    def test_RSG_mapping_transcript_range(self):
        """
        Test that a working input range for a RSG mapping transcript
        produces an equal and valid RSG output in the VV results
        """
        variant = 'NM_004006.2:c.-120_-114T[7]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_004006.2:c.-120_-114T[7]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == 'NC_000023.10:g.33229543_33229549A[7]'
        assert results["NM_004006.2:c.-120_-114T[7]"]['hgvs_refseqgene_variant']\
            == 'NG_012232.1:g.133178_133184T[7]'
        assert results["NM_004006.2:c.-120_-114T[7]"]["validation_warnings"] ==  [
            'ExpandedRepeatWarning: NM_004006.2:c.-120_-114T[7] should only be used as an '
            'annotation for the core HGVS descriptions provided',
            'TranscriptVersionWarning: A more recent version of the selected reference sequence '
            'NM_004006.2 is available for genome build GRCh37 (NM_004006.3)'
            ]

    def test_RSG_mapping_transcript_single_pos(self):
        """
        Test that a working input range for a RSG mapping transcript
        produces an equal and valid RSG output in the VV results
        """
        variant = 'NM_004006.2:c.-120T[7]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_004006.2:c.-120_-114T[7]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == 'NC_000023.10:g.33229543_33229549A[7]'
        assert results["NM_004006.2:c.-120_-114T[7]"]['hgvs_refseqgene_variant']\
            == 'NG_012232.1:g.133178_133184T[7]'
        assert results["NM_004006.2:c.-120_-114T[7]"]["validation_warnings"] ==   [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NM_004006.2:c.-120T[7]. The corrected description is NM_004006.2:c.-120_-114T[7]",
            "ExpandedRepeatWarning: NM_004006.2:c.-120_-114T[7] should only be used as an annotation for the core HGVS descriptions provided",
            "TranscriptVersionWarning: A more recent version of the selected reference sequence NM_004006.2 is available for genome build GRCh37 (NM_004006.3)"
        ]

    def test_RSG_mapping_transcript_intron_single_pos(self):
        """
        Test that a working input range for an intron inside a RSG mapping transcript
        produces an equal and valid RSG output in the VV results
        """
        variant = 'NM_022167.4:c.135+1G[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_022167.4:c.135_135+1G[2]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48423636_48423637G[2]"
        assert results["NM_022167.4:c.135_135+1G[2]"]['genome_context_intronic_sequence']\
            == 'NC_000017.10(NM_022167.4):c.135_135+1='
        assert results["NM_022167.4:c.135_135+1G[2]"]['refseqgene_context_intronic_sequence']\
            == 'NG_012175.1(NM_022167.4):c.135_135+1='
        assert results["NM_022167.4:c.135_135+1G[2]"]['hgvs_refseqgene_variant']\
            == 'NG_012175.1:g.5244_5245G[2]'
        assert results["NM_022167.4:c.135_135+1G[2]"]["validation_warnings"] == [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description NM_022167.4:c.135+1G[2]. The corrected description is NM_022167.4:c.135_135+1G[2]",
            "ExpandedRepeatWarning: NM_022167.4:c.135_135+1G[2] should only be used as an annotation for the core HGVS descriptions provided",
        ]

    def test_RSG_mapping_transcript_intron_range(self):
        """
        Test that a working input range for an intron inside a RSG mapping transcript
        produces an equal and valid RSG output in the VV results
        """
        variant = 'NM_022167.4:c.135_135+1G[2]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_022167.4:c.135_135+1G[2]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000017.10:g.48423636_48423637G[2]"
        assert results["NM_022167.4:c.135_135+1G[2]"]['genome_context_intronic_sequence']\
            == 'NC_000017.10(NM_022167.4):c.135_135+1='
        assert results["NM_022167.4:c.135_135+1G[2]"]['refseqgene_context_intronic_sequence']\
            == 'NG_012175.1(NM_022167.4):c.135_135+1='
        assert results["NM_022167.4:c.135_135+1G[2]"]['hgvs_refseqgene_variant']\
            == 'NG_012175.1:g.5244_5245G[2]'
        assert results["NM_022167.4:c.135_135+1G[2]"]["validation_warnings"] == [
            "ExpandedRepeatWarning: NM_022167.4:c.135_135+1G[2] should only be used as an "
            "annotation for the core HGVS descriptions provided",
        ]


class TestExpandedRepeaLocusReferenceGenomic(TestCase):
    """
    Tests for Locus Reference Genomic targeting expanded repeats.
    This is tied to the RSG tests above the code runs LRG->RSG so if the RSG
    fails then so should the LRG.
    """
    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()
        cls.vv.testing = True

    def test_basic_syntax_LRG(self):
        """
        Test for handling basic syntax of variant string, with LRG genomic
        reference type variant input directly via TandemRepeats code
        Should pass if all basic syntax checks passed.
        """
        variant_str = "LRG_199:g.4T[20]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
            variant_str,  "GRCh37", "all", vv)
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat(vv)
        assert str(formatted) == "NG_012232.1:g.3_6T[20]"
        assert my_variant.variant_str == "LRG_199:g.4T[20]"
        # checks correct transcript ref
        assert my_variant.reference == "NG_012232.1"
        # checks correct position
        assert str(my_variant.variant_position) == "3_6"
        # checks repeat seq
        assert my_variant.repeat_sequence == "T"
        # checks correct suffix
        assert my_variant.copy_number == "20"
        # checks number of repeats is str and correct
        assert my_variant.after_the_bracket == ""
        # checks nothing is after the bracket

    def test_basic_syntax_LRG_t(self):
        """
        Test for handling basic syntax of variant string, with LRG transcript
        reference type variant input directly via TandemRepeats code
        Should pass if all basic syntax checks passed.
        """
        variant_str = "LRG_199t1:c.-120T[7]"
        my_variant = expanded_repeats.TandemRepeats.parse_repeat_variant(
            variant_str,  "GRCh37", "all", vv)
        my_variant.reformat_reference()
        my_variant.check_genomic_or_coding()
        formatted = my_variant.reformat(vv)
        assert str(formatted) == "NM_004006.2:c.-120_-114T[7]"
        assert my_variant.variant_str == "LRG_199t1:c.-120T[7]"
        # checks correct transcript ref
        assert my_variant.reference == "NM_004006.2"
        # checks correct position
        assert str(my_variant.variant_position) == "-120_-114"
        # checks repeat seq
        assert my_variant.repeat_sequence == "T"
        # checks correct suffix
        assert my_variant.copy_number == "7"
        # checks number of repeats is str and correct
        assert my_variant.after_the_bracket == ""
        # checks nothing is after the bracket

    def test_LRG_genomic_range(self):
        """
        Test that a working input range produces an equal and valid output
        """
        variant = 'LRG_199:g.3_6T[20]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] ==   [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description LRG_199:g.3_6T[20]. The corrected format would be NG_012232.1:g.3_6T[int], where int requires you to update the number of repeats"
        ]

    def test_LRG_genomic_single_position_to_span(self):
        """
        Test mapping of single location starting tandem repeat to span via full VV front end
        """
        variant = "LRG_199:g.4T[20]"
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["intergenic_variant_1"]["primary_assembly_loci"]["grch37"][
                "hgvs_genomic_description"] == "NC_000023.10:g.33362724_33362725insAAAAAAAAAAAAAAAA"
        assert results['intergenic_variant_1']["validation_warnings"] == [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description LRG_199:g.4T[20]. The corrected description is NG_012232.1:g.3_6T[20]",
            "ExpandedRepeatWarning: NG_012232.1:g.3_6T[20] should only be used as an annotation for the core HGVS descriptions provided",
            "NG_012232.1:g.6_7insTTTTTTTTTTTTTTTT automapped to genome position NC_000023.10:g.33362724_33362725insAAAAAAAAAAAAAAAA",
            "No individual transcripts have been identified that fully overlap the described variation in the genomic sequence. Large variants might span one or more genes and are currently only described at the genome (g.) level."
        ]

    def test_incorrect_LRG_genomic_range(self):
        """
        Test that an invalid range gives an appropriate error response
        """
        variant = 'LRG_199:g.3_6C[20]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] == [
            'RepeatSyntaxError: The provided repeat sequence C does not match the reference '
            'sequence T at the given position 3_3 of reference sequence NG_012232.1'
        ]

    def test_LRG_mapping_transcript_range(self):
        """
        Test that a working input range for a RSG mapping transcript
        produces an equal and valid RSG output in the VV results
        """
        variant = 'LRG_199t1:c.-120_-114T[7]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] ==  [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description LRG_199t1:c.-120_-114T[7]. The corrected format would be NM_004006.2:c.-120_-114T[int], where int requires you to update the number of repeats"
        ]

    def test_LRG_mapping_transcript_single_pos(self):
        """
        Test that a working input range for a RSG mapping transcript
        produces an equal and valid RSG output in the VV results
        """
        variant = 'LRG_199t1:c.-120T[7]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_004006.2:c.-120_-114T[7]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == 'NC_000023.10:g.33229543_33229549A[7]'
        assert results["NM_004006.2:c.-120_-114T[7]"]['hgvs_refseqgene_variant']\
            == 'NG_012232.1:g.133178_133184T[7]'
        assert results["NM_004006.2:c.-120_-114T[7]"]["validation_warnings"] == [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description LRG_199t1:c.-120T[7]. The corrected description is NM_004006.2:c.-120_-114T[7]",
            "ExpandedRepeatWarning: NM_004006.2:c.-120_-114T[7] should only be used as an annotation for the core HGVS descriptions provided",
            "TranscriptVersionWarning: A more recent version of the selected reference sequence NM_004006.2 is available for genome build GRCh37 (NM_004006.3)"
        ]

    def test_LRG_mapping_transcript_intron_single_pos(self):
        """
        Test that a working input range for an intron inside a RSG mapping transcript
        produces an equal and valid RSG output in the VV results
        """
        variant = 'LRG_199t1:c.31+9A[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_004006.2:c.31+9_31+11A[3]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000023.10:g.33229388_33229390T[3]"
        assert results["NM_004006.2:c.31+9_31+11A[3]"]['genome_context_intronic_sequence']\
            == 'NC_000023.10(NM_004006.2):c.31+9_31+11='
        assert results["NM_004006.2:c.31+9_31+11A[3]"]['refseqgene_context_intronic_sequence']\
            == 'NG_012232.1(NM_004006.2):c.31+9_31+11='
        assert results["NM_004006.2:c.31+9_31+11A[3]"]['hgvs_lrg_transcript_variant']\
            == 'LRG_199t1:c.31+9_31+11A[3]'
        assert results["NM_004006.2:c.31+9_31+11A[3]"]['hgvs_refseqgene_variant']\
            == 'NG_012232.1:g.133337_133339A[3]'
        assert results["NM_004006.2:c.31+9_31+11A[3]"]['hgvs_lrg_variant']\
            == 'LRG_199:g.133337_133339A[3]'
        assert results["NM_004006.2:c.31+9_31+11A[3]"]["validation_warnings"] ==  [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description LRG_199t1:c.31+9A[3]. The corrected description is NM_004006.2:c.31+9_31+11A[3]",
            "ExpandedRepeatWarning: NM_004006.2:c.31+9_31+11A[3] should only be used as an annotation for the core HGVS descriptions provided",
            "TranscriptVersionWarning: A more recent version of the selected reference sequence NM_004006.2 is available for genome build GRCh37 (NM_004006.3)"
        ]

    def test_LRG_mapping_transcript_intron_range(self):
        """
        Test that a working input range for an intron inside a RSG mapping transcript
        produces an equal and valid RSG output in the VV results
        """
        variant = 'LRG_199t1:c.31+9_31+11A[3]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] == [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description LRG_199t1:c.31+9_31+11A[3]. The corrected format would be NM_004006.2:c.31+9_31+11A[int], where int requires you to update the number of repeats"
        ]

    def test_gap_crossing(self):
        variant = 'LRG_763t1:c.54GCA[21]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert results["NM_002111.8:c.54_116GCA[21]"]["primary_assembly_loci"]["grch37"][
            "hgvs_genomic_description"] == "NC_000004.11:g.3076606_3076662GCA[21]"
        assert results["NM_002111.8:c.54_116GCA[21]"]["validation_warnings"] ==   [
            "ExpandedRepeatError: The coordinates for the repeat region are stated incorrectly in the submitted description LRG_763t1:c.54GCA[21]. The corrected description is NM_002111.8:c.54_116GCA[21]",
            "ExpandedRepeatWarning: NM_002111.8:c.54_116GCA[21] should only be used as an annotation for the core HGVS descriptions provided"
        ]


if __name__ == "__main__":
    unittest.main()


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

