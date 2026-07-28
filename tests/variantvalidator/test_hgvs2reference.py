import unittest

import VariantValidator


class TestHgvs2Reference(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vv = VariantValidator.Validator()

    def test_empty(self):
        result = self.vv.hgvs2ref("")

        self.assertEqual(
            list(result),
            [
                "variant",
                "start_position",
                "end_position",
                "warning",
                "sequence",
                "error",
            ],
        )
        self.assertEqual(result["error"], ": char 1: end of input")

    def test_nonsense(self):
        result = self.vv.hgvs2ref("nonsense")

        self.assertEqual(result["sequence"], "")
        self.assertEqual(result["warning"], "")
        self.assertEqual(
            result["error"],
            "nonsense: char 9: end of input",
        )

    def test_nonsense_colon(self):
        result = self.vv.hgvs2ref("non:sense")

        self.assertEqual(result["sequence"], "")
        self.assertEqual(result["warning"], "")
        self.assertEqual(
            result["error"],
            "non:sense: char 4: expected one of "
            "'c', 'g', 'm', 'n', 'p', or 'r'",
        )

    def test_nonsense_hgvs(self):
        result = self.vv.hgvs2ref("nonsense:c.34C>T")

        self.assertEqual(result["sequence"], "")
        self.assertEqual(result["warning"], "")
        self.assertIn(
            "No transcript definition for (tx_ac=nonsense)",
            result["error"],
        )

    def test_invalid_hgvs(self):
        result = self.vv.hgvs2ref("this_is_not_hgvs")

        self.assertEqual(result["sequence"], "")
        self.assertNotEqual(result["error"], "")

    def test_genomic_single_base(self):
        result = self.vv.hgvs2ref(
            "NC_000017.11:g.43071077C>T"
        )

        self.assertEqual(
            result["variant"],
            "NC_000017.11:g.43071077C>T",
        )
        self.assertEqual(result["start_position"], "43071077")
        self.assertEqual(result["end_position"], "43071077")
        self.assertEqual(result["sequence"], "T")
        self.assertEqual(result["warning"], "")
        self.assertEqual(result["error"], "")

    def test_genomic_interval(self):
        result = self.vv.hgvs2ref(
            "NC_000017.11:g.43071077_43071079del"
        )

        self.assertEqual(result["start_position"], "43071077")
        self.assertEqual(result["end_position"], "43071079")
        self.assertEqual(len(result["sequence"]), 3)
        self.assertEqual(result["warning"], "")
        self.assertEqual(result["error"], "")

    def test_coding(self):
        result = self.vv.hgvs2ref(
            "NM_007294.4:c.68_69del"
        )

        self.assertEqual(
            result["variant"],
            "NM_007294.4:c.68_69del",
        )
        self.assertEqual(result["start_position"], "68")
        self.assertEqual(result["end_position"], "69")
        self.assertEqual(len(result["sequence"]), 2)
        self.assertEqual(result["warning"], "")
        self.assertEqual(result["error"], "")

    def test_noncoding(self):
        result = self.vv.hgvs2ref(
            "NR_027676.2:n.100_102del"
        )

        self.assertEqual(
            result["variant"],
            "NR_027676.2:n.100_102del",
        )
        self.assertEqual(result["start_position"], "100")
        self.assertEqual(result["end_position"], "102")
        self.assertEqual(len(result["sequence"]), 3)
        self.assertEqual(result["warning"], "")
        self.assertEqual(result["error"], "")

    def test_coding_intronic_without_genomic_context(self):
        result = self.vv.hgvs2ref(
            "NM_007294.4:c.80+1G>A"
        )

        self.assertEqual(
            result["variant"],
            "NM_007294.4:c.80+1G>A",
        )
        self.assertEqual(result["sequence"], "")
        self.assertEqual(result["warning"], "")
        self.assertIn(
            "Unable to establish the exon structure",
            result["error"],
        )

    def test_bad_reference_accession(self):
        result = self.vv.hgvs2ref(
            "NM_999999.1:c.1A>G"
        )

        self.assertEqual(result["sequence"], "")
        self.assertNotEqual(result["error"], "")

    def test_identity_variant(self):
        result = self.vv.hgvs2ref(
            "NM_000088.4:c.589="
        )

        self.assertEqual(result["error"], "")
        self.assertNotEqual(result["sequence"], "")

    def test_multibase_interval(self):
        result = self.vv.hgvs2ref(
            "NC_000017.11:g.50198002_50198005del"
        )

        self.assertEqual(result["error"], "")
        self.assertEqual(len(result["sequence"]), 4)

    def test_protein_unsupported(self):
        result = self.vv.hgvs2ref(
            "NP_009225.1:p.Glu23Val"
        )

        self.assertEqual(result["sequence"], "")
        self.assertEqual(result["warning"], "")
        self.assertIn(
            "Unsupported HGVS reference type 'p.'",
            result["error"],
        )

    def test_rna_unsupported(self):
        result = self.vv.hgvs2ref(
            "NM_007294.4:r.68_69del"
        )

        self.assertEqual(result["sequence"], "")
        self.assertEqual(result["warning"], "")
        self.assertIn(
            "Unsupported HGVS reference type 'r.'",
            result["error"],
        )

    def test_mitochondrial_unsupported(self):
        result = self.vv.hgvs2ref(
            "NC_012920.1:m.100A>G"
        )

        self.assertEqual(result["sequence"], "")
        self.assertEqual(result["warning"], "")
        self.assertIn(
            "Unsupported HGVS reference type 'm.'",
            result["error"],
        )


if __name__ == "__main__":
    unittest.main()