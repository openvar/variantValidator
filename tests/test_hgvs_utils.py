import unittest
from unittest.mock import Mock, patch
from VariantValidator.modules.hgvs_utils import pvcf_to_hgvs, PseudoVCF2HGVSError
from VariantValidator.modules.hgvs_utils import pre_push_vcf_tx_g_map_fix

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
        result = self.validate(
            "NC_000017.11:43045705AG>A"
        )

        print(result)
        assert result["flag"] == "warning"
        assert result["validation_warning_1"]["validation_warnings"] == [
                 'VariantMappingWarning: NC_000017.11:g.43045705AG>A automapped to NC_000017.11:g.43045706del',
                 'ReferenceMismatchError: NC_000017.11:g.43045706delG: Variant reference (G) does not agree with reference sequence (A)']


    def test_hybrid_genomic_insertion(self):
        result = self.validate(
            "NC_000017.11:43045705A>AG"
        )

        print(result)

        assert result["flag"] == "gene_variant"
        assert "NM_001407571.1:c.5351_5352insC" in result.keys()
        assert result["NM_001407571.1:c.5351_5352insC"]["validation_warnings"] == [
            'VariantMappingWarning: NC_000017.11:g.43045705A>AG automapped to NC_000017.11:g.43045705_43045706insG']



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

    def test_pre_push_vcf_tx_g_map_fix_mapped_no_fix_required(self):
        hp = self.vv.hp

        transcript = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101del"
        )
        mapped = hp.parse_hgvs_variant(
            "NC_000017.11:g.100_101del"
        )

        var_mapper = Mock()
        normaliser = Mock()
        normaliser.normalize.return_value = mapped

        result = pre_push_vcf_tx_g_map_fix(
            transcript,
            transcript,
            "NC_000017.11",
            var_mapper,
            normaliser,
            Mock(),
            mapped,
        )

        self.assertIs(result, transcript)
        var_mapper.n_to_g.assert_not_called()
        normaliser.normalize.assert_called_once_with(mapped)

    def test_pre_push_vcf_tx_g_map_fix_maps_when_not_premapped(self):
        hp = self.vv.hp

        transcript = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101del"
        )
        mapped = hp.parse_hgvs_variant(
            "NC_000017.11:g.100_101del"
        )

        var_mapper = Mock()
        var_mapper.n_to_g.return_value = mapped

        normaliser = Mock()
        normaliser.normalize.return_value = mapped

        result = pre_push_vcf_tx_g_map_fix(
            transcript,
            transcript,
            "NC_000017.11",
            var_mapper,
            normaliser,
            Mock(),
            False,
        )

        self.assertIs(result, transcript)

        var_mapper.n_to_g.assert_called_once_with(
            transcript,
            "NC_000017.11",
        )

    def test_pre_push_vcf_tx_g_map_fix_other_invalid_variant_error(self):
        hp = self.vv.hp

        transcript = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101del"
        )
        mapped = hp.parse_hgvs_variant(
            "NC_000017.11:g.100_101del"
        )

        normaliser = Mock()
        normaliser.normalize.side_effect = (
            vvhgvs.exceptions.HGVSInvalidVariantError(
                "some other invalid variant"
            )
        )

        result = pre_push_vcf_tx_g_map_fix(
            transcript,
            transcript,
            "NC_000017.11",
            Mock(),
            normaliser,
            Mock(),
            mapped,
        )

        self.assertIs(result, transcript)

    def test_pre_push_vcf_tx_g_map_fix_gap_without_flanks(self):
        hp = self.vv.hp

        transcript = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101delAAinsTT"
        )

        genomic = hp.parse_hgvs_variant(
            "NC_000017.11:g.201_200delAAinsTT"
        )

        genomic_identity = hp.parse_hgvs_variant(
            "NC_000017.11:g.200_201="
        )

        transcript_gap = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101="
        )

        normaliser = Mock()

        normaliser.normalize.side_effect = [
            vvhgvs.exceptions.HGVSInvalidVariantError(
                "base start position must be <= end position"
            ),
            genomic_identity,
            transcript_gap,
        ]

        var_mapper = Mock()
        var_mapper.g_to_n.return_value = transcript_gap

        result = pre_push_vcf_tx_g_map_fix(
            transcript,
            transcript,
            "NC_000017.11",
            var_mapper,
            normaliser,
            Mock(),
            genomic,
        )

        self.assertEqual(
            result.posedit.pos.start.base,
            100,
        )
        self.assertEqual(
            result.posedit.pos.end.base,
            101,
        )
        self.assertEqual(result.posedit.edit.ref, "AA")
        self.assertEqual(result.posedit.edit.alt, "TT")

        var_mapper.g_to_n.assert_called_once()

    def test_pre_push_vcf_tx_g_map_fix_left_flank(self):
        hp = self.vv.hp

        transcript = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101delAAinsTT"
        )
        genomic = hp.parse_hgvs_variant(
            "NC_000017.11:g.201_200delAAinsTT"
        )
        genomic_identity = hp.parse_hgvs_variant(
            "NC_000017.11:g.200_201="
        )
        transcript_gap = hp.parse_hgvs_variant(
            "NM_000546.6:n.90_101="
        )
        left_flank = hp.parse_hgvs_variant(
            "NM_000546.6:n.90_99="
        )

        normaliser = Mock()
        normaliser.normalize.side_effect = [
            vvhgvs.exceptions.HGVSInvalidVariantError(
                "base start position must be <= end position"
            ),
            genomic_identity,
            transcript_gap,
            left_flank,
        ]

        var_mapper = Mock()
        var_mapper.g_to_n.return_value = transcript_gap

        result = pre_push_vcf_tx_g_map_fix(
            transcript,
            transcript,
            "NC_000017.11",
            var_mapper,
            normaliser,
            Mock(),
            genomic,
        )

        self.assertEqual(result.posedit.pos.start.base, 90)
        self.assertEqual(result.posedit.pos.end.base, 101)
        self.assertEqual(
            result.posedit.edit.ref,
            left_flank.posedit.edit.ref
            + transcript.posedit.edit.ref,
        )
        self.assertEqual(
            result.posedit.edit.alt,
            left_flank.posedit.edit.alt
            + transcript.posedit.edit.alt,
        )

    def test_pre_push_vcf_tx_g_map_fix_right_flank(self):
        hp = self.vv.hp

        transcript = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101delAAinsTT"
        )
        genomic = hp.parse_hgvs_variant(
            "NC_000017.11:g.201_200delAAinsTT"
        )
        genomic_identity = hp.parse_hgvs_variant(
            "NC_000017.11:g.200_201="
        )
        transcript_gap = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_110="
        )
        right_flank = hp.parse_hgvs_variant(
            "NM_000546.6:n.102_110="
        )

        normaliser = Mock()
        normaliser.normalize.side_effect = [
            vvhgvs.exceptions.HGVSInvalidVariantError(
                "base start position must be <= end position"
            ),
            genomic_identity,
            transcript_gap,
            right_flank,
        ]

        var_mapper = Mock()
        var_mapper.g_to_n.return_value = transcript_gap

        result = pre_push_vcf_tx_g_map_fix(
            transcript,
            transcript,
            "NC_000017.11",
            var_mapper,
            normaliser,
            Mock(),
            genomic,
        )

        self.assertEqual(result.posedit.pos.start.base, 100)
        self.assertEqual(result.posedit.pos.end.base, 110)
        self.assertEqual(
            result.posedit.edit.ref,
            transcript.posedit.edit.ref
            + right_flank.posedit.edit.ref,
        )
        self.assertEqual(
            result.posedit.edit.alt,
            transcript.posedit.edit.alt
            + right_flank.posedit.edit.alt,
        )

    def test_pre_push_vcf_tx_g_map_fix_both_flanks(self):
        hp = self.vv.hp

        transcript = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101delAAinsTT"
        )
        genomic = hp.parse_hgvs_variant(
            "NC_000017.11:g.201_200delAAinsTT"
        )
        genomic_identity = hp.parse_hgvs_variant(
            "NC_000017.11:g.200_201="
        )
        transcript_gap = hp.parse_hgvs_variant(
            "NM_000546.6:n.90_110="
        )
        left_flank = hp.parse_hgvs_variant(
            "NM_000546.6:n.90_99="
        )
        right_flank = hp.parse_hgvs_variant(
            "NM_000546.6:n.102_110="
        )

        normaliser = Mock()
        normaliser.normalize.side_effect = [
            vvhgvs.exceptions.HGVSInvalidVariantError(
                "base start position must be <= end position"
            ),
            genomic_identity,
            transcript_gap,
            left_flank,
            right_flank,
        ]

        var_mapper = Mock()
        var_mapper.g_to_n.return_value = transcript_gap

        result = pre_push_vcf_tx_g_map_fix(
            transcript,
            transcript,
            "NC_000017.11",
            var_mapper,
            normaliser,
            Mock(),
            genomic,
        )

        self.assertEqual(result.posedit.pos.start.base, 90)
        self.assertEqual(result.posedit.pos.end.base, 110)

        self.assertEqual(
            result.posedit.edit.ref,
            left_flank.posedit.edit.ref
            + transcript.posedit.edit.ref
            + right_flank.posedit.edit.ref,
        )
        self.assertEqual(
            result.posedit.edit.alt,
            left_flank.posedit.edit.alt
            + transcript.posedit.edit.alt
            + right_flank.posedit.edit.alt,
        )

    @patch(
        "VariantValidator.modules.hgvs_utils.hgvs_to_delins_hgvs"
    )
    def test_pre_push_vcf_tx_g_map_fix_genomic_duplication(
            self,
            mock_to_delins,
    ):
        hp = self.vv.hp

        transcript = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101delAAinsTT"
        )

        genomic = hp.parse_hgvs_variant(
            "NC_000017.11:g.201_200dup"
        )

        genomic_delins = hp.parse_hgvs_variant(
            "NC_000017.11:g.200_201delAAinsAAAA"
        )

        genomic_identity = hp.parse_hgvs_variant(
            "NC_000017.11:g.200_201="
        )

        transcript_gap = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101="
        )

        mock_to_delins.return_value = genomic_delins

        seq_fetcher = Mock()
        seq_fetcher.fetch_seq.return_value = "AA"

        normaliser = Mock()
        normaliser.normalize.side_effect = [
            vvhgvs.exceptions.HGVSInvalidVariantError(
                "base start position must be <= end position"
            ),
            genomic_identity,
            transcript_gap,
        ]

        var_mapper = Mock()
        var_mapper.g_to_n.return_value = transcript_gap

        result = pre_push_vcf_tx_g_map_fix(
            transcript,
            transcript,
            "NC_000017.11",
            var_mapper,
            normaliser,
            seq_fetcher,
            genomic,
        )

        seq_fetcher.fetch_seq.assert_called_once_with(
            "NC_000017.11",
            199,
            201,
        )

        mock_to_delins.assert_called_once()

        self.assertEqual(result.posedit.pos.start.base, 100)
        self.assertEqual(result.posedit.pos.end.base, 101)

    @patch(
        "VariantValidator.modules.hgvs_utils.hgvs_to_delins_hgvs"
    )
    def test_pre_push_vcf_tx_g_map_fix_transcript_insertion(
            self,
            mock_to_delins,
    ):
        hp = self.vv.hp

        transcript = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101insTT"
        )
        transcript_delins = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101delinsTT"
        )
        genomic = hp.parse_hgvs_variant(
            "NC_000017.11:g.201_200delinsTT"
        )
        genomic_identity = hp.parse_hgvs_variant(
            "NC_000017.11:g.200_201="
        )
        transcript_gap = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101="
        )

        mock_to_delins.return_value = transcript_delins

        normaliser = Mock()
        normaliser.normalize.side_effect = [
            vvhgvs.exceptions.HGVSInvalidVariantError(
                "base start position must be <= end position"
            ),
            genomic_identity,
            transcript_gap,
        ]

        var_mapper = Mock()
        var_mapper.g_to_n.return_value = transcript_gap

        result = pre_push_vcf_tx_g_map_fix(
            transcript,
            transcript,
            "NC_000017.11",
            var_mapper,
            normaliser,
            Mock(),
            genomic,
        )

        mock_to_delins.assert_called_once_with(
            transcript,
            None,
            normaliser,
        )

        self.assertEqual(result.posedit.pos.start.base, 100)
        self.assertEqual(result.posedit.pos.end.base, 101)


    @patch(
        "VariantValidator.modules.hgvs_utils.hgvs_to_delins_hgvs"
    )
    def test_pre_push_vcf_tx_g_map_fix_transcript_duplication(
            self,
            mock_to_delins,
    ):
        hp = self.vv.hp

        transcript = hp.parse_hgvs_variant(
            "NM_000546.6:n.100dupA"
        )
        transcript_delins = hp.parse_hgvs_variant(
            "NM_000546.6:n.100delAinsAA"
        )
        genomic = hp.parse_hgvs_variant(
            "NC_000017.11:g.201_200delinsAA"
        )
        genomic_identity = hp.parse_hgvs_variant(
            "NC_000017.11:g.200="
        )
        transcript_gap = hp.parse_hgvs_variant(
            "NM_000546.6:n.100="
        )

        mock_to_delins.return_value = transcript_delins

        normaliser = Mock()
        normaliser.normalize.side_effect = [
            vvhgvs.exceptions.HGVSInvalidVariantError(
                "base start position must be <= end position"
            ),
            genomic_identity,
            transcript_gap,
        ]

        var_mapper = Mock()
        var_mapper.g_to_n.return_value = transcript_gap

        result = pre_push_vcf_tx_g_map_fix(
            transcript,
            transcript,
            "NC_000017.11",
            var_mapper,
            normaliser,
            Mock(),
            genomic,
        )

        mock_to_delins.assert_called_once_with(
            transcript,
            None,
            normaliser,
        )

        self.assertEqual(result.posedit.pos.start.base, 100)
        self.assertEqual(result.posedit.pos.end.base, 100)


    @patch(
        "VariantValidator.modules.hgvs_utils.hgvs_to_delins_hgvs"
    )
    def test_pre_push_vcf_tx_g_map_fix_unnormalised_deletion(
            self,
            mock_to_delins,
    ):
        hp = self.vv.hp

        # The normalised representation is already a delins, but the
        # original input was a deletion. This exercises the explicit
        # unnormalised-deletion check.
        transcript = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101delAAinsTT"
        )
        unnormalised = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101delAA"
        )

        genomic = hp.parse_hgvs_variant(
            "NC_000017.11:g.201_200delAAinsTT"
        )
        genomic_identity = hp.parse_hgvs_variant(
            "NC_000017.11:g.200_201="
        )
        transcript_gap = hp.parse_hgvs_variant(
            "NM_000546.6:n.100_101="
        )

        mock_to_delins.return_value = transcript

        normaliser = Mock()
        normaliser.normalize.side_effect = [
            vvhgvs.exceptions.HGVSInvalidVariantError(
                "base start position must be <= end position"
            ),
            genomic_identity,
            transcript_gap,
        ]

        var_mapper = Mock()
        var_mapper.g_to_n.return_value = transcript_gap

        result = pre_push_vcf_tx_g_map_fix(
            transcript,
            unnormalised,
            "NC_000017.11",
            var_mapper,
            normaliser,
            Mock(),
            genomic,
        )

        mock_to_delins.assert_called_once_with(
            transcript,
            None,
            normaliser,
        )

        self.assertEqual(result.posedit.pos.start.base, 100)
        self.assertEqual(result.posedit.pos.end.base, 101)


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
