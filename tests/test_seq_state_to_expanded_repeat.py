import unittest
from vvhgvs import exceptions
from VariantValidator import Validator
from VariantValidator.modules import seq_state_to_expanded_repeat
from VariantValidator.modules.seq_state_to_expanded_repeat import \
        convert_seq_state_to_expanded_repeat, VariantFormatError, reassemble_expanded_repeat_variant,\
        decipher_end_of_full_reference_repeated_sequence,RepeatedUnitError, quick_testfunc

class TestExpandedRepeatConversion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.validator = Validator()

    def test_expanded_repeat_polyA_variant(self):
        variant = self.validator.hp.parse("NC_000023.11:g.33344607_33344608insAAAAAAAAAAAAAAAA")
        result = convert_seq_state_to_expanded_repeat(variant, self.validator)
        self.assertEqual(str(result), "NC_000023.11:g.33344604_33344607A[20]")

    def test_expanded_repeat_polyT_variant(self):
        variant = self.validator.hp.parse("NG_012232.1:g.6_7insTTTTTTTTTTTTTTTT")
        result = convert_seq_state_to_expanded_repeat(variant, self.validator)
        self.assertEqual(str(result), "NG_012232.1:g.3_6T[20]")

    def test_invalid_variant_format_raises(self):
        # adjusted for the string-> pre-parsed hgvs object changes
        with self.assertRaises(VariantFormatError):
            variant = self.validator.hp.parse("NM_002111.8:r.54_116=")
            convert_seq_state_to_expanded_repeat(variant, self.validator)

    def test_no_seq_state_returns_none(self):
        # A valid variant format but contains a seq change
        # so can not be == to count * input repeat unit
        variant = self.validator.hp.parse("NC_000023.11:g.33344607_33344609delinsGGCT")
        with self.assertRaises(VariantFormatError):
            convert_seq_state_to_expanded_repeat(variant, self.validator)

    def test_expanded_repeat_polyT_variant_alt(self):
        variant = self.validator.hp.parse("NG_012232.1:g.6_7insTTTTTTTTTTTTTTTT")
        result = convert_seq_state_to_expanded_repeat(variant, self.validator)
        self.assertEqual(str(result), "NG_012232.1:g.3_6T[20]")

    def test_expanded_repeat_equal_variant(self):
        variant = self.validator.hp.parse("NM_002111.8:c.54_116=")
        result = convert_seq_state_to_expanded_repeat(variant, self.validator)
        self.assertEqual(str(result), "NM_002111.8:c.54_116GCA[21]")

    def test_expanded_repeat_duplication_variant(self):
        variant = self.validator.hp.parse("NM_002111.8:c.54_116dup")
        result = convert_seq_state_to_expanded_repeat(variant, self.validator)
        self.assertEqual(str(result), "NM_002111.8:c.54_116GCA[42]")

    def test_expanded_repeat_deletion_of_one_repeat(self):
        variant = self.validator.hp.parse("NM_002111.8:c.114_116del")
        result = convert_seq_state_to_expanded_repeat(variant, self.validator)
        self.assertEqual(str(result), "NM_002111.8:c.54_116GCA[21]")

    def test_expanded_repeat_deletion_of_two_repeats(self):
        variant = self.validator.hp.parse("NM_002111.8:c.111_116del")
        result = convert_seq_state_to_expanded_repeat(variant, self.validator)
        self.assertEqual(str(result), "NM_002111.8:c.54_116GCA[20]")

    def test_decipher_repeated_unit_non_repeat(self):
        """Test decipher_repeated_unit with a non-repeating sequence."""
        result = seq_state_to_expanded_repeat.decipher_repeated_unit("ACGT")
        self.assertEqual(result, "ACGT")  # Should return full sequence if no repeat found

    def test_decipher_repeated_unit_single_base_repeat(self):
        """Test decipher_repeated_unit with single-nucleotide repeat."""
        result = seq_state_to_expanded_repeat.decipher_repeated_unit("AAAAAAAA")
        self.assertEqual(result, "A")

    def test_start_of_repeat_invalid_repeated_unit(self):
        with self.assertRaises(seq_state_to_expanded_repeat.RepeatedUnitError):
            seq_state_to_expanded_repeat.decipher_start_of_full_reference_repeated_sequence("REF", "", 100, self.validator)

    def test_start_of_repeat_invalid_start_position(self):
        with self.assertRaises(seq_state_to_expanded_repeat.StartPositionError):
            seq_state_to_expanded_repeat.decipher_start_of_full_reference_repeated_sequence("REF", "A", -10, self.validator)

    def test_end_of_repeat_invalid_repeated_unit(self):
        with self.assertRaises(seq_state_to_expanded_repeat.RepeatedUnitError):
            seq_state_to_expanded_repeat.decipher_end_of_full_reference_repeated_sequence("REF", "", 100, self.validator)

    def test_end_of_repeat_invalid_start_position(self):
        with self.assertRaises(seq_state_to_expanded_repeat.StartPositionError):
            seq_state_to_expanded_repeat.decipher_end_of_full_reference_repeated_sequence("REF", "A", 0, self.validator)

    def test_expanded_repeat_coding_variant_negative_positions(self):
        variant = self.validator.hp.parse("NM_004006.2:c.-3_1=")
        expected = "NM_004006.2:c.-3_1A[4]"

        result = convert_seq_state_to_expanded_repeat(variant, validator=self.validator)
        self.assertEqual(str(result), expected)

    def test_expanded_repeat_noncoding_negative_positions(self):
        variant = self.validator.hp.parse("NR_110010.2:n.15_16=")
        expected = "NR_110010.2:n.15_16GA[1]"

        result = convert_seq_state_to_expanded_repeat(variant, validator=self.validator)
        self.assertEqual(str(result), expected)

    def test_expanded_repeat_coding_3utr_positions(self):
        variant = self.validator.hp.parse("NM_001160367.2:c.870_*1=")
        expected = "NM_001160367.2:c.870_*1AC[1]"

        result = convert_seq_state_to_expanded_repeat(variant, validator=self.validator)
        self.assertEqual(str(result), expected)

    def test_intronic_coding_sense_strand(self):
        variant = self.validator.hp.parse("NM_000492.4:c.1210-34_1210-13=")
        expected = "NM_000492.4:c.1210-34_1210-13TG[11]"

        result = convert_seq_state_to_expanded_repeat(variant, validator=self.validator, genomic_reference="NC_000007.13")
        self.assertEqual(str(result), expected)

    def test_intronic_coding_antisense_strand(self):
        variant = self.validator.hp.parse("NM_000088.3:c.589-1_590=")
        expected = "NM_000088.3:c.589-1_590G[3]"

        result = convert_seq_state_to_expanded_repeat(variant, validator=self.validator, genomic_reference="NC_000017.10")
        self.assertEqual(str(result), expected)

    # known_repeat_unit non ins, then ins, then revcomp
    def test_known_repeat_unit_nonins(self):
        variant = self.validator.hp.parse("NM_002111.8:c.111_116del")
        result = convert_seq_state_to_expanded_repeat(variant, self.validator,known_repeat_unit='GCA')
        self.assertEqual(str(result), "NM_002111.8:c.54_116GCA[20]")

    def test_known_repeat_unit_ins(self):
        # normalisation auto expands ins to dup if it can, so test as much as we can here,
        # see also test_pure_ins_variant
        variant = self.validator.hp.parse("NM_002111.8:c.116_117insGCA")
        result = convert_seq_state_to_expanded_repeat(variant, self.validator,known_repeat_unit='GCA')
        self.assertEqual(str(result), "NM_002111.8:c.54_116GCA[22]")
        variant = self.validator.hp.parse("NM_002111.8:c.53_54insGCA")
        result = convert_seq_state_to_expanded_repeat(variant, self.validator,known_repeat_unit='GCA')
        self.assertEqual(str(result), "NM_002111.8:c.54_116GCA[22]")
        # force to go through latter code as ins because ins > 22 reps, however norms to end of span
        # hence 53_54 is redundant for now, keep in for testing robustness
        variant = self.validator.hp.parse("NM_002111.8:c.53_54insGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA")
        result = convert_seq_state_to_expanded_repeat(variant, self.validator,known_repeat_unit='GCA')
        self.assertEqual(str(result), "NM_002111.8:c.54_116GCA[45]")
        variant = self.validator.hp.parse("NM_002111.8:c.116_117insGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA")
        result = convert_seq_state_to_expanded_repeat(variant, self.validator,known_repeat_unit='GCA')
        self.assertEqual(str(result), "NM_002111.8:c.54_116GCA[45]")

    def test_known_repeat_unit(self):
        variant = self.validator.hp.parse("NM_002111.8:c.111_116del")
        result = convert_seq_state_to_expanded_repeat(variant, self.validator,known_repeat_unit='TGC')
        self.assertEqual(str(result), "NM_002111.8:c.54_116GCA[20]")

    def test_bad_genomic_mapping(self):
        variant = self.validator.hp.parse("NM_000492.4:c.1210-34_1210-13=")
        with self.assertRaises(VariantFormatError):
            result = convert_seq_state_to_expanded_repeat(
                    variant, validator=self.validator,
                    genomic_reference="NC_000023.11")

    def test_no_genomic_mapping(self):
        variant = self.validator.hp.parse("NM_000492.4:c.1210-34_1210-13=")
        with self.assertRaises(VariantFormatError):
            result = convert_seq_state_to_expanded_repeat(
                    variant, validator=self.validator)

    def test_safe_null_variant(self):
        result = convert_seq_state_to_expanded_repeat(
                    '', validator=self.validator)
        assert result == ''
        result = convert_seq_state_to_expanded_repeat(
                    None, validator=self.validator)
        assert result == None

    def test_pure_ins_variant(self):
        # pure ins should not count as a HGVS Repeated Sequence to quote:
        # Repeated sequence: a sequence where, compared to a reference sequence,
        # a segment of one or more nucleotides (the repeat unit) is present
        # several times, one after the other.
        variant = self.validator.hp.parse("NM_002111.8:c.116_117insTCTCTC")
        with self.assertRaises(VariantFormatError):
            result = convert_seq_state_to_expanded_repeat(
                    variant, validator=self.validator,known_repeat_unit='TC')
        variant = self.validator.hp.parse("NM_002111.8:c.116_117insTCTCTC")
        with self.assertRaises(VariantFormatError):
            result = convert_seq_state_to_expanded_repeat(
                    variant, validator=self.validator)

    def test_reassemble_expanded_repeat_variant_input_err(self):
        # we should never encounter this, but for now test the current behaviour
        with self.assertRaises(RepeatedUnitError):
            reassemble_expanded_repeat_variant('','','','', '','','')

    def test_decipher_end_of_full_reference_repeated_sequence_end_lt_zero(self):
        # this can happen since we try from start - rep length -1 if start = 1
        # we now avoid starting end fetch from start, so as not to duplicate checks,
        # but it should still work
        end = decipher_end_of_full_reference_repeated_sequence('NM_002111.8','GCT', 1, self.validator)
        assert end == 3
        # we should never encounter this, seqfetch should already fail at an earlier point in this case,
        # but for now test the current behaviour
        end = decipher_end_of_full_reference_repeated_sequence('NM_XXXXXXX.Y','GCT', 6, self.validator)
        assert end == 6

    def test_internal_quick_test(self):
        res = quick_testfunc()
        assert str(res) == 'NM_002111.8:c.54_116GCA[21]'

    def test_bad_genomic_style_attempt(self):
        # test that bad attempts to map such as those caused by shifted or
        # changed sequence during mappings external to the expanded repeat code
        # fail as expected
        variant = self.validator.hp.parse("NM_002111.8:c.52_116=")
        with self.assertRaises(VariantFormatError):
            result = convert_seq_state_to_expanded_repeat(
                    variant, validator=self.validator, known_repeat_unit='GCA')

    def test_intronic_coding_sense_strand(self):
        # test that bad attempts to map, specifically like those caused by
        # changed intronic sequence mappings fail as expected
        variant = self.validator.hp.parse(
                "NM_000492.4:c.1210-34_1210-13CGCGCGCGCGCGCGCGCGCGCG=")

        with self.assertRaises(VariantFormatError):
            result = convert_seq_state_to_expanded_repeat(
                    variant, validator=self.validator, known_repeat_unit='TG',
                    genomic_reference="NC_000007.13")


if __name__ == "__main__":
    unittest.main()
