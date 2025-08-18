import unittest
from vvhgvs import exceptions
from VariantValidator import Validator
from VariantValidator.modules import seq_state_to_expanded_repeat
from VariantValidator.modules.seq_state_to_expanded_repeat import convert_seq_state_to_expanded_repeat, VariantFormatError

class TestExpandedRepeatConversion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.validator = Validator()

    def test_expanded_repeat_polyA_variant(self):
        variant = "NC_000023.11:g.33344607_33344608insAAAAAAAAAAAAAAAA"
        result = convert_seq_state_to_expanded_repeat(variant, self.validator)
        self.assertEqual(result, "NC_000023.11:g.33344604_33344607A[20]")

    def test_expanded_repeat_polyT_variant(self):
        variant = "NG_012232.1:g.6_7insTTTTTTTTTTTTTTTT"
        result = convert_seq_state_to_expanded_repeat(variant, self.validator)
        self.assertEqual(result, "NG_012232.1:g.3_6T[20]")

    def test_invalid_variant_format_raises(self):
        with self.assertRaises(VariantFormatError):
            convert_seq_state_to_expanded_repeat("invalid_variant_string", self.validator)

    def test_missing_colon_format_raises(self):
        # This will fail because it doesn't contain :g. or :c. etc.
        variant = "NG_012232.1g.6_7insTTTTTTTTTTTTTTTT"
        with self.assertRaises(exceptions.HGVSParseError):
            convert_seq_state_to_expanded_repeat(variant, self.validator)

    def test_no_seq_state_returns_none(self):
        # A valid variant format but no 'ins', 'del', or 'dup'
        variant = "NC_000023.11:g.33344607_33344608"
        with self.assertRaises(VariantFormatError):
            convert_seq_state_to_expanded_repeat(variant, self.validator)

    def test_expanded_repeat_polyT_variant_alt(self):
        variant = "NG_012232.1:g.6_7insTTTTTTTTTTTTTTTT"
        result = convert_seq_state_to_expanded_repeat(variant, self.validator)
        self.assertEqual(result, "NG_012232.1:g.3_6T[20]")

    def test_expanded_repeat_equal_variant(self):
        variant = "NM_002111.8:c.54_116="
        result = convert_seq_state_to_expanded_repeat(variant, self.validator)
        self.assertEqual(result, "NM_002111.8:c.54_116GCA[21]")

    def test_expanded_repeat_duplication_variant(self):
        variant = "NM_002111.8:c.54_116dup"
        result = convert_seq_state_to_expanded_repeat(variant, self.validator)
        self.assertEqual(result, "NM_002111.8:c.54_116GCA[42]")

    def test_expanded_repeat_deletion_of_one_repeat(self):
        variant = "NM_002111.8:c.114_116del"
        result = convert_seq_state_to_expanded_repeat(variant, self.validator)
        self.assertEqual(result, "NM_002111.8:c.54_116GCA[21]")

    def test_expanded_repeat_deletion_of_two_repeats(self):
        variant = "NM_002111.8:c.111_116del"
        result = convert_seq_state_to_expanded_repeat(variant, self.validator)
        self.assertEqual(result, "NM_002111.8:c.54_116GCA[20]")

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
        variant = "NM_004006.2:c.-3_1="
        expected = "NM_004006.2:c.-3_1A[4]"

        result = convert_seq_state_to_expanded_repeat(variant, validator=self.validator)
        self.assertEqual(result, expected)

    def test_expanded_repeat_noncoding_negative_positions(self):
        variant = "NR_110010.2:n.15_16="
        expected = "NR_110010.2:n.15_16GA[1]"

        result = convert_seq_state_to_expanded_repeat(variant, validator=self.validator)
        self.assertEqual(result, expected)

    def test_expanded_repeat_coding_3utr_positions(self):
        variant = "NM_001160367.2:c.870_*1="
        expected = "NM_001160367.2:c.870_*1AC[1]"

        result = convert_seq_state_to_expanded_repeat(variant, validator=self.validator)
        self.assertEqual(result, expected)

    def test_intronic_coding_sense_strand(self):
        variant = "NM_000492.4:c.1210-34_1210-13="
        expected = "NM_000492.4:c.1210-34_1210-13TG[11]"

        result = convert_seq_state_to_expanded_repeat(variant, validator=self.validator, genomic_reference="NC_000007.13")
        self.assertEqual(result, expected)

    def test_intronic_coding_antisense_strand(self):
        variant = "NM_000088.3:c.589-1_590="
        expected = "NM_000088.3:c.589-1_590G[3]"

        result = convert_seq_state_to_expanded_repeat(variant, validator=self.validator, genomic_reference="NC_000017.10")
        self.assertEqual(result, expected)

if __name__ == "__main__":
    unittest.main()
