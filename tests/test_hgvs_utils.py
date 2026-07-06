import unittest
from unittest.mock import Mock, patch
from VariantValidator.modules.hgvs_utils import pvcf_to_hgvs, PseudoVCF2HGVSError

from unittest import TestCase

import vvhgvs

from VariantValidator.validator import Validator
from VariantValidator.modules.hgvs_utils import (
    hgvs_dup_to_delins,
    hgvs_to_delins_hgvs,
    unset_hgvs_obj_ref,
    vcfcp_to_hgvs_obj,
    vcfcp_to_hgvsstr,
)


class TestHgvsUtilsFunctional(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()

    def validate(self, variant, assembly="GRCh38", transcripts="all"):
        return self.vv.validate(
            variant,
            assembly,
            transcripts,
        ).format_as_dict(test=True)

    def test_hybrid_genomic_deletion(self):
        result = self.validate("NC_000017.11:43045705AG>A")

        # inspect output then replace with real assertions
        print(result)

    def test_hybrid_genomic_insertion(self):
        result = self.validate("NC_000017.11:43045705A>AG")

        # inspect output then replace with real assertions
        print(result)


class TestHgvsUtils(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()

    def test_vcfcp_to_hgvsstr(self):
        start = vvhgvs.parser.Parser().parse_hgvs_variant(
            "NC_000017.11:g.43045705A>G"
        )

        hgvs = vcfcp_to_hgvsstr(
            {
                "pos": 43045705,
                "ref": "AG",
                "alt": "A",
            },
            start,
        )

        self.assertEqual(
            hgvs,
            "NC_000017.11:g.43045705_43045706delAGinsA",
        )

    def test_vcfcp_to_hgvs_obj(self):
        start = vvhgvs.parser.Parser().parse_hgvs_variant(
            "NC_000017.11:g.43045705A>G"
        )

        hgvs = vcfcp_to_hgvs_obj(
            {
                "pos": 43045705,
                "ref": "AG",
                "alt": "A",
            },
            start,
        )

        self.assertEqual(hgvs.ac, "NC_000017.11")
        self.assertEqual(hgvs.type, "g")
        self.assertEqual(hgvs.posedit.pos.start.base, 43045705)
        self.assertEqual(hgvs.posedit.pos.end.base, 43045706)
        self.assertEqual(hgvs.posedit.edit.ref, "AG")
        self.assertEqual(hgvs.posedit.edit.alt, "A")

    def test_unset_hgvs_obj_ref_delins(self):
        hp = self.vv.hp

        hgvs = hp.parse_hgvs_variant(
            "NC_000017.11:g.43045705_43045706delAGinsTT"
        )

        unset = unset_hgvs_obj_ref(hgvs)

        self.assertEqual(unset.posedit.edit.ref, "")
        self.assertEqual(unset.posedit.edit.alt, "TT")

    def test_unset_hgvs_obj_ref_duplication(self):
        hp = self.vv.hp

        hgvs = hp.parse_hgvs_variant(
            "NC_000017.11:g.43045705dupG"
        )

        unset = unset_hgvs_obj_ref(hgvs)

        self.assertEqual(unset.posedit.edit.ref, "")

    def test_unset_hgvs_obj_ref_deletion(self):
        hp = self.vv.hp

        hgvs = hp.parse_hgvs_variant(
            "NC_000017.11:g.43045705delG"
        )

        unset = unset_hgvs_obj_ref(hgvs)

        self.assertEqual(unset.posedit.edit.ref, "")

    def test_hgvs_dup_to_delins(self):
        hp = self.vv.hp

        dup = hp.parse_hgvs_variant(
            "NC_000017.11:g.43045705dupG"
        )

        delins = hgvs_dup_to_delins(dup)

        self.assertEqual(delins.posedit.edit.type, "delins")
        self.assertEqual(
            delins.posedit.edit.alt,
            delins.posedit.edit.ref * 2,
        )

    def test_hgvs_to_delins_hgvs_existing_delins(self):
        hp = self.vv.hp
        hn = Mock()
        hn.normalize.side_effect = lambda x: x

        hgvs = hp.parse_hgvs_variant(
            "NC_000017.11:g.43045705_43045706delAGinsTT"
        )

        result = hgvs_to_delins_hgvs(hgvs, hp, hn)

        self.assertEqual(result, hgvs)

    def test_hgvs_to_delins_hgvs_deletion(self):
        hp = self.vv.hp
        hn = Mock()
        hn.normalize.side_effect = lambda x: x

        hgvs = hp.parse_hgvs_variant(
            "NC_000017.11:g.43045705delG"
        )

        result = hgvs_to_delins_hgvs(hgvs, hp, hn)

        self.assertEqual(result.posedit.edit.type, "delins")
        self.assertEqual(result.posedit.edit.alt, "")

    def test_hgvs_to_delins_hgvs_duplication(self):
        hp = self.vv.hp
        hn = Mock()
        hn.normalize.side_effect = lambda x: x

        hgvs = hp.parse_hgvs_variant(
            "NC_000017.11:g.43045705dupG"
        )

        result = hgvs_to_delins_hgvs(hgvs, hp, hn)

        self.assertEqual(result.posedit.edit.type, "delins")
        self.assertEqual(
            result.posedit.edit.alt,
            result.posedit.edit.ref * 2,
        )

    def test_hgvs_to_delins_hgvs_insertion(self):
        hp = self.vv.hp
        hn = vvhgvs.normalizer.Normalizer(self.vv.hdp,
                                          cross_boundaries=False,
                                          shuffle_direction=3,
                                          alt_aln_method="splign"
                                         )

        hgvs = hp.parse_hgvs_variant(
            "NC_000017.11:g.43045705_43045706insTT"
        )

        result = hgvs_to_delins_hgvs(hgvs, hp, hn)

        self.assertEqual(result.posedit.edit.type, "delins")
        self.assertIn("TT", result.posedit.edit.alt)

    def test_pvcf_multibase_deletion(self):
        hn = vvhgvs.normalizer.Normalizer(
            self.vv.hdp,
            cross_boundaries=False,
            shuffle_direction=3,
            alt_aln_method="splign",
        )

        self.vv.hn = hn

        hgvs = pvcf_to_hgvs(
            "chr16-2099572-TC-T",
            selected_assembly="GRCh38",
            normalization_direction=3,
            reverse_normalizer=hn,
            validator=self.vv,
        )

        self.assertIsNotNone(hgvs)


    def test_pvcf_multibase_insertion(self):
        hn = vvhgvs.normalizer.Normalizer(
            self.vv.hdp,
            cross_boundaries=False,
            shuffle_direction=3,
            alt_aln_method="splign",
        )

        self.vv.hn = hn

        hgvs = pvcf_to_hgvs(
            "chr16-2099572-T-TC",
            selected_assembly="GRCh38",
            normalization_direction=3,
            reverse_normalizer=hn,
            validator=self.vv,
        )

        self.assertIsNotNone(hgvs)


class TestPVCFtoHGVS(unittest.TestCase):

    def setUp(self):
        # Mock validator object with .hn.normalize
        self.mock_validator = Mock()
        self.mock_validator.hn.normalize = Mock(side_effect=lambda x: x)
        self.mock_validator.db.get_refseq_id_from_lrg_id = Mock(side_effect=lambda x: "NM_000000.1")

        # Mock reverse_normalizer
        self.mock_reverse = Mock()
        self.mock_reverse.normalize = Mock(side_effect=lambda x: x)

    @patch('VariantValidator.modules.hgvs_utils.seq_data.to_accession')
    @patch('VariantValidator.modules.hgvs_utils.hgvs_delins_parts_to_hgvs_obj')
    def test_simple_substitution(self, mock_hgvs_obj, mock_to_accession):
        # Setup mocks
        mock_to_accession.return_value = "NM_000000.1"
        mock_hgvs_obj.return_value = "HGVS_OBJ"

        # Provide a simple pVCF string
        query = "chr1-123-A-T"

        # Call the function
        result = pvcf_to_hgvs(query, selected_assembly="GRCh38", normalization_direction=3,
                              reverse_normalizer=self.mock_reverse, validator=self.mock_validator)

        # Assertions
        self.assertEqual(result, "HGVS_OBJ")
        self.mock_validator.hn.normalize.assert_called_once()
        mock_to_accession.assert_called_once_with("1", "GRCh38")
        mock_hgvs_obj.assert_called()  # At least called once

    def test_unsupported_format(self):
        query = "chr1-123-A"  # missing ALT
        with self.assertRaises(PseudoVCF2HGVSError):
            pvcf_to_hgvs(query, selected_assembly="GRCh38", normalization_direction=3,
                          reverse_normalizer=self.mock_reverse, validator=self.mock_validator)

if __name__ == '__main__':
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
