import unittest

from VariantValidator.validator import Validator
from VariantValidator.modules.use_checking import (
    refseq_common_mistakes,
    refseq_type_mismatch,
)


class MockRefVariant:
    def __init__(self, quibble, reftype, transcript_type=None):
        self.quibble = quibble
        self.reftype = reftype
        self.transcript_type = transcript_type
        self.warnings = []


class MockHgvsVariant:
    def __init__(self, accession, variant_type, posedit="1A>G"):
        self.ac = accession
        self.type = variant_type
        self.posedit = posedit

    def __str__(self):
        return f"{self.ac}:{self.type}.{self.posedit}"


def make_refseq_hgvs_object(accession, variant_type, posedit="1A>G"):
    return MockHgvsVariant(accession, variant_type, posedit)


class TestUseCheckingReferenceTypeChecks(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()

    def test_refseq_common_mistakes_ignores_object(self):
        hgvs_variant = make_refseq_hgvs_object("NM_000001.1", "g")
        variant = MockRefVariant(hgvs_variant, ":g.", "c")
        self.assertFalse(refseq_common_mistakes(variant, self.vv))

    def test_refseq_common_mistakes_nm_as_genomic(self):
        variant = MockRefVariant("NM_000001.1:g.1A>G", ":g.", "c")
        self.assertTrue(refseq_common_mistakes(variant, self.vv))
        self.assertIn("Did you mean NM_000001.1:c.1A>G?", variant.warnings[0])

    def test_refseq_common_mistakes_nr_as_genomic(self):
        variant = MockRefVariant("NR_000001.1:g.10A>G", ":g.", "n")
        self.assertTrue(refseq_common_mistakes(variant, self.vv))
        self.assertIn("Did you mean NR_000001.1:n.10A>G?", variant.warnings[0])

    def test_refseq_common_mistakes_nr_as_coding(self):
        variant = MockRefVariant("NR_000001.1:c.1A>G", ":c.", "n")
        self.assertTrue(refseq_common_mistakes(variant, self.vv))
        self.assertIn("Did you mean NR_000001.1:n.1A>G?", variant.warnings[0])

    def test_refseq_common_mistakes_coding_transcript_as_noncoding(self):
        variant = MockRefVariant("NM_000001.1:n.1A>G", ":n.", "c")
        self.assertTrue(refseq_common_mistakes(variant, self.vv))
        self.assertIn("Did you mean NM_000001.1:c.1A>G?", variant.warnings[0])

    def test_refseq_common_mistakes_protein_as_nucleotide(self):
        for accession in ("NP_000001.1", "ENSP000001.1"):
            for reftype in (":c.", ":n.", ":g.", ":r."):
                with self.subTest(accession=accession, reftype=reftype):
                    variant = MockRefVariant(
                        f"{accession}{reftype}1A>G",
                        reftype,
                    )
                    self.assertTrue(
                        refseq_common_mistakes(variant, self.vv)
                    )
                    self.assertIn("ReferenceTypeError:", variant.warnings[0])

    def test_refseq_common_mistakes_nucleotide_as_protein(self):
        accessions = (
            "NM_000001.1",
            "NR_000001.1",
            "NC_000001.11",
            "NG_000001.1",
            "NT_000001.1",
            "NW_000001.1",
            "ENST000001.1",
        )
        for accession in accessions:
            with self.subTest(accession=accession):
                variant = MockRefVariant(
                    f"{accession}:p.Gly1Val",
                    ":p.",
                )
                self.assertTrue(
                    refseq_common_mistakes(variant, self.vv)
                )

    def test_refseq_common_mistakes_genomic_without_context(self):
        for accession in (
            "NC_000001.11",
            "NG_000001.1",
            "NT_000001.1",
            "NW_000001.1",
        ):
            for reftype in (":c.", ":n."):
                with self.subTest(accession=accession, reftype=reftype):
                    variant = MockRefVariant(
                        f"{accession}{reftype}123A>G",
                        reftype,
                    )
                    self.assertTrue(
                        refseq_common_mistakes(variant, self.vv)
                    )

    def test_refseq_common_mistakes_genomic_r(self):
        for accession in (
            "NC_000001.11",
            "NG_000001.1",
            "NT_000001.1",
            "NW_000001.1",
        ):
            with self.subTest(accession=accession):
                variant = MockRefVariant(
                    f"{accession}:r.123a>g",
                    ":r.",
                )
                self.assertTrue(
                    refseq_common_mistakes(variant, self.vv)
                )

    def test_refseq_common_mistakes_valid_genomic_g(self):
        for accession in (
            "NC_000001.11",
            "NG_000001.1",
            "NT_000001.1",
            "NW_000001.1",
        ):
            with self.subTest(accession=accession):
                variant = MockRefVariant(
                    f"{accession}:g.123A>G",
                    ":g.",
                )
                self.assertFalse(
                    refseq_common_mistakes(variant, self.vv)
                )

    def test_refseq_common_mistakes_valid_compound_nm_c(self):
        for genomic in (
            "NC_000001.11",
            "NG_000001.1",
            "NT_000001.1",
            "NW_000001.1",
        ):
            with self.subTest(genomic=genomic):
                variant = MockRefVariant(
                    f"{genomic}(NM_000001.1):c.123A>G",
                    ":c.",
                )
                self.assertFalse(
                    refseq_common_mistakes(variant, self.vv)
                )

    def test_refseq_common_mistakes_valid_compound_nr_n(self):
        for genomic in (
            "NC_000001.11",
            "NG_000001.1",
            "NT_000001.1",
            "NW_000001.1",
        ):
            with self.subTest(genomic=genomic):
                variant = MockRefVariant(
                    f"{genomic}(NR_000001.1):n.123A>G",
                    ":n.",
                )
                self.assertFalse(
                    refseq_common_mistakes(variant, self.vv)
                )

    def test_refseq_common_mistakes_rejects_nr_c_compound(self):
        variant = MockRefVariant(
            "NG_000001.1(NR_000001.1):c.123A>G",
            ":c.",
        )
        self.assertTrue(refseq_common_mistakes(variant, self.vv))
        self.assertIn("Did you mean", variant.warnings[0])

    def test_refseq_common_mistakes_rejects_nm_n_compound(self):
        variant = MockRefVariant(
            "NG_000001.1(NM_000001.1):n.123A>G",
            ":n.",
        )
        self.assertTrue(refseq_common_mistakes(variant, self.vv))
        self.assertIn("Did you mean", variant.warnings[0])

    def test_refseq_common_mistakes_valid_lrg_compound(self):
        for genomic in ("NC_000001.11", "NG_000001.1"):
            with self.subTest(genomic=genomic):
                variant = MockRefVariant(
                    f"{genomic}(LRG_1t1):c.123A>G",
                    ":c.",
                )
                self.assertFalse(
                    refseq_common_mistakes(variant, self.vv)
                )

    def test_refseq_common_mistakes_lrg_n_rejected(self):
        variant = MockRefVariant(
            "NG_000001.1(LRG_1t1):n.123A>G",
            ":n.",
        )
        self.assertTrue(refseq_common_mistakes(variant, self.vv))

    def test_refseq_common_mistakes_valid_transcript(self):
        variant = MockRefVariant(
            "NM_000001.1:c.123A>G",
            ":c.",
            "c",
        )
        self.assertFalse(refseq_common_mistakes(variant, self.vv))
        self.assertEqual(variant.warnings, [])

    def test_refseq_type_mismatch_ignores_string(self):
        variant = MockRefVariant(
            "NM_000001.1:g.1A>G",
            ":g.",
            "c",
        )
        self.assertFalse(refseq_type_mismatch(variant, self.vv))

    def test_refseq_type_mismatch_nm_as_genomic(self):
        hgvs_variant = make_refseq_hgvs_object("NM_000001.1", "g")
        variant = MockRefVariant(hgvs_variant, ":g.", "c")
        self.assertTrue(refseq_type_mismatch(variant, self.vv))
        self.assertIn("Did you mean NM_000001.1:c.", variant.warnings[0])

    def test_refseq_type_mismatch_nr_as_genomic(self):
        hgvs_variant = make_refseq_hgvs_object("NR_000001.1", "g")
        variant = MockRefVariant(hgvs_variant, ":g.", "n")
        self.assertTrue(refseq_type_mismatch(variant, self.vv))
        self.assertIn("Did you mean NR_000001.1:n.", variant.warnings[0])

    def test_refseq_type_mismatch_noncoding_as_coding(self):
        hgvs_variant = make_refseq_hgvs_object("NR_000001.1", "c")
        variant = MockRefVariant(hgvs_variant, ":c.", "n")
        self.assertTrue(refseq_type_mismatch(variant, self.vv))

    def test_refseq_type_mismatch_coding_as_noncoding(self):
        hgvs_variant = make_refseq_hgvs_object("NM_000001.1", "n")
        variant = MockRefVariant(hgvs_variant, ":n.", "c")
        self.assertTrue(refseq_type_mismatch(variant, self.vv))

    def test_refseq_type_mismatch_protein_as_nucleotide(self):
        for accession in ("NP_000001.1", "ENSP000001.1"):
            for variant_type in ("c", "n", "g", "r"):
                with self.subTest(
                    accession=accession,
                    variant_type=variant_type,
                ):
                    hgvs_variant = make_refseq_hgvs_object(
                        accession,
                        variant_type,
                    )
                    variant = MockRefVariant(
                        hgvs_variant,
                        f":{variant_type}.",
                    )
                    self.assertTrue(
                        refseq_type_mismatch(variant, self.vv)
                    )

    def test_refseq_type_mismatch_nucleotide_as_protein(self):
        for accession in (
            "NM_000001.1",
            "NR_000001.1",
            "NC_000001.11",
            "NG_000001.1",
            "NT_000001.1",
            "NW_000001.1",
            "ENST000001.1",
        ):
            with self.subTest(accession=accession):
                hgvs_variant = make_refseq_hgvs_object(
                    accession,
                    "p",
                    "Gly1Val",
                )
                variant = MockRefVariant(hgvs_variant, ":p.")
                self.assertTrue(
                    refseq_type_mismatch(variant, self.vv)
                )

    def test_refseq_type_mismatch_genomic_c_without_context(self):
        for accession in (
            "NC_000001.11",
            "NG_000001.1",
            "NT_000001.1",
            "NW_000001.1",
        ):
            with self.subTest(accession=accession):
                hgvs_variant = make_refseq_hgvs_object(accession, "c")
                variant = MockRefVariant(hgvs_variant, ":c.", "c")
                self.assertTrue(
                    refseq_type_mismatch(variant, self.vv)
                )

    def test_refseq_type_mismatch_valid_variants_unchanged(self):
        cases = (
            ("NM_000001.1", "c", "c"),
            ("NR_000001.1", "n", "n"),
            ("NC_000001.11", "g", "c"),
            ("NG_000001.1", "g", "c"),
            ("NT_000001.1", "g", "c"),
            ("NW_000001.1", "g", "c"),
        )
        for accession, variant_type, transcript_type in cases:
            with self.subTest(
                accession=accession,
                variant_type=variant_type,
            ):
                hgvs_variant = make_refseq_hgvs_object(
                    accession,
                    variant_type,
                )
                variant = MockRefVariant(
                    hgvs_variant,
                    f":{variant_type}.",
                    transcript_type,
                )
                self.assertFalse(
                    refseq_type_mismatch(variant, self.vv)
                )
                self.assertEqual(variant.warnings, [])


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
