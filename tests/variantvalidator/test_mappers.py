import copy
import json
from unittest import TestCase
from unittest.mock import MagicMock, patch

import vvhgvs.exceptions

from VariantValidator.modules.mappers import (
    TranscriptMappingError,
    gene_to_transcripts,
    final_tx_to_multiple_genomic,
    transcripts_to_gene,
    MappersError
)

from VariantValidator.validator import Validator

class TestMappersFunctional(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()

    def validate(self, variant, assembly="GRCh38", transcripts="all"):
        return self.vv.validate(
            variant,
            assembly,
            transcripts,
        ).format_as_dict(test=True)

    def test_issue518_forward_intron_correction(self):
        result = self.validate(
            "NM_000086.2:c.791-802_1056+1445del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["NM_000086.2:c.790+532_1056+1445del"]["validation_warnings"] == [
            "ExonBoundaryError: Position c.791-802 has been updated to position to 790+532 ensuring correct HGVS numbering for transcript NM_000086.2"
        ]

    def test_issue518_reverse_intron_correction(self):
        result = self.validate(
            "NM_000088.4:c.2559+54_2560del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["NM_000088.4:c.2560-34_2561del"]["validation_warnings"] == [
            "ExonBoundaryError: Position c.2559+54 has been updated to position to 2560-35 ensuring correct HGVS numbering for transcript NM_000088.4"
        ]

    def test_issue518_double_boundary_correction(self):
        result = self.validate(
            "NM_000088.4:c.2559+54_2559+55del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["NM_000088.4:c.2560-34_2560-33del"]["validation_warnings"] == [
            "ExonBoundaryError: Position c.2559+54 has been updated to position to 2560-35 ensuring correct HGVS numbering for transcript NM_000088.4",
            "ExonBoundaryError: Position c.2559+55 has been updated to position to 2560-34 ensuring correct HGVS numbering for transcript NM_000088.4"
        ]

    def test_issue518_left_boundary_correction(self):
        result = self.validate(
            "NM_000086.2:c.790_791-802del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["NM_000086.2:c.790+1_790+533del"]["validation_warnings"] ==  [
            "ExonBoundaryError: Position c.791-802 has been updated to position to 790+532 ensuring correct HGVS numbering for transcript NM_000086.2",
            "VariantNormalizationWarning: NM_000086.2:c.790_790+532del normalized to NM_000086.2:c.790+1_790+533del"
        ]

    def test_issue518_two_reverse_boundaries(self):
        result = self.validate(
            "NM_000086.2:c.791-803_791-802del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result[ "NM_000086.2:c.790+531_790+532del"]["validation_warnings"] ==  [
            "ExonBoundaryError: Position c.791-803 has been updated to position to 790+531 ensuring correct HGVS numbering for transcript NM_000086.2",
            "ExonBoundaryError: Position c.791-802 has been updated to position to 790+532 ensuring correct HGVS numbering for transcript NM_000086.2"
        ]

    def test_issue518_invalid_start_boundary(self):
        result = self.validate(
            "NM_000088.4:c.543-1_543+1del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["validation_warning_1"]["validation_warnings"] == [
            "ExonBoundaryError: Position c.543-1 does not correspond with an exon boundary for transcript NM_000088.4"
        ]

    def test_issue518_invalid_end_boundary(self):
        result = self.validate(
            "NM_000088.4:c.2560-3_2560+3dup"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["validation_warning_1"]["validation_warnings"] == [
            "ExonBoundaryError: Position c.2560+3 does not correspond with an exon boundary for transcript NM_000088.4"
        ]

    def test_issue518_invalid_end_boundary2(self):
        result = self.validate(
            "NM_000088.4:c.2608_2608+5del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["validation_warning_1"]["validation_warnings"] == [
            "ExonBoundaryError: Position c.2608+5 does not correspond with an exon boundary for transcript NM_000088.4"
        ]

    def test_issue518_interval_error(self):
        result = self.validate(
            "NM_000088.4:c.543_453+1del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["validation_warning_1"]["validation_warnings"] == [
            "IntervalOrderError: Interval end position 453 < interval start position 543"
        ]

    def test_issue518_reverse_end_only(self):
        result = self.validate(
            "NM_000088.4:c.2559_2559+54del"
        )
        print(json.dumps(result, sort_keys=True, indent=4, separators=(',', ': ')))
        assert result["NM_000088.4:c.2559_2560-35del"]["validation_warnings"] == [
            "ExonBoundaryError: Position c.2559+54 has been updated to position to 2560-35 ensuring correct HGVS numbering for transcript NM_000088.4"
        ]


class TestGeneToTranscriptsUnit(TestCase):

    def setUp(self):
        self.variant = MagicMock()
        self.validator = MagicMock()
        self.batch_list = []

        self.variant.warnings = []
        self.variant.map_dat.hdp = None
        self.variant.primary_assembly = "GRCh38"
        self.variant.reverse_normalizer = MagicMock()
        self.variant.evm = MagicMock()
        self.variant.no_norm_evm = MagicMock()
        self.variant.hn = MagicMock()

        self.validator.hdp = MagicMock()
        self.validator.alt_aln_method = "splign"
        self.validator.select_transcripts = "all"

        hp = vvhgvs.parser.Parser()

        self.g_query = hp.parse_hgvs_variant(
            "NC_000001.11:g.100A>G"
        )
        self.g_test = copy.deepcopy(self.g_query)

        self.variant.hgvs_formatted = self.g_query
        self.variant.hgvs_genomic = self.g_query

        self.variant.hn.normalize.return_value = self.g_test
        self.validator.vr.validate.return_value = None

    def test_sets_hdp_when_missing(self):
        self.validator.relevant_transcripts.return_value = []

        with patch(
            "VariantValidator.modules.mappers.seq_data.supported_for_mapping",
            return_value=False,
        ):
            gene_to_transcripts(
                self.variant,
                self.validator,
                {},
                self.batch_list,
            )

        self.assertIs(
            self.variant.map_dat.hdp,
            self.validator.hdp,
        )

    def test_validation_hgvs_error(self):
        self.validator.vr.validate.side_effect = (
            vvhgvs.exceptions.HGVSInvalidVariantError(
                "bad genomic variant"
            )
        )

        result = gene_to_transcripts(
            self.variant,
            self.validator,
            {},
            self.batch_list,
        )

        self.assertTrue(result)
        self.assertIn(
            "bad genomic variant",
            self.variant.warnings,
        )

    def test_validation_key_error(self):
        self.variant.hgvs_genomic.ac = "NC_BAD"

        self.validator.vr.validate.side_effect = KeyError(
            "missing"
        )

        result = gene_to_transcripts(
            self.variant,
            self.validator,
            {},
            self.batch_list,
        )

        self.assertTrue(result)
        self.assertEqual(
            self.variant.warnings,
            [
                "Reference sequence NC_BAD is either not supported "
                "or does not exist"
            ],
        )

    def test_undefined_n_reference_warning(self):
        self.g_test.posedit.edit.ref = "ANNNT"
        self.validator.relevant_transcripts.return_value = []

        with patch(
            "VariantValidator.modules.mappers.seq_data.supported_for_mapping",
            return_value=False,
        ):
            gene_to_transcripts(
                self.variant,
                self.validator,
                {},
                self.batch_list,
            )

        self.assertIn(
            "UndefinedSequenceError: Submitted variant description cannot "
            "be fully validated because it spans a region of the reference "
            "sequence represented by base 'N' and not bases 'GATC'",
            self.variant.warnings,
        )

    def test_reference_attribute_error_is_ignored(self):
        del self.g_test.posedit.edit.ref

        self.validator.relevant_transcripts.return_value = []

        with patch(
            "VariantValidator.modules.mappers.seq_data.supported_for_mapping",
            return_value=False,
        ):
            result = gene_to_transcripts(
                self.variant,
                self.validator,
                {},
                self.batch_list,
            )

        self.assertTrue(result)

    def test_reference_type_error_is_ignored(self):
        self.g_test.posedit.edit.ref = None
        self.validator.relevant_transcripts.return_value = []

        with patch(
            "VariantValidator.modules.mappers.seq_data.supported_for_mapping",
            return_value=False,
        ):
            result = gene_to_transcripts(
                self.variant,
                self.validator,
                {},
                self.batch_list,
            )

        self.assertTrue(result)

    def test_normalized_coordinates_replace_genomic_variant(self):
        hp = vvhgvs.parser.Parser()

        self.g_test = hp.parse_hgvs_variant(
            "NC_000001.11:g.101A>G"
        )
        self.variant.hn.normalize.return_value = self.g_test

        self.validator.relevant_transcripts.return_value = []

        with patch(
            "VariantValidator.modules.mappers.seq_data.supported_for_mapping",
            return_value=False,
        ):
            gene_to_transcripts(
                self.variant,
                self.validator,
                {},
                self.batch_list,
            )

        self.assertIs(
            self.variant.hgvs_genomic,
            self.g_test,
        )
        self.assertIs(
            self.variant.hgvs_formatted,
            self.g_test,
        )

    def test_unnormalized_coordinates_keep_original_variant(self):
        self.g_test = copy.deepcopy(self.g_query)
        self.variant.hn.normalize.return_value = self.g_test

        self.validator.relevant_transcripts.return_value = []

        with patch(
            "VariantValidator.modules.mappers.seq_data.supported_for_mapping",
            return_value=False,
        ):
            gene_to_transcripts(
                self.variant,
                self.validator,
                {},
                self.batch_list,
            )

        self.assertEqual(
            str(self.variant.hgvs_genomic),
            str(self.g_query),
        )

    def test_refseqgene_mapping_hgvs_error(self):
        hp = vvhgvs.parser.Parser()

        self.g_query = hp.parse_hgvs_variant(
            "NG_000001.1:g.100A>G"
        )
        self.g_test = copy.deepcopy(self.g_query)

        self.variant.hgvs_formatted = self.g_query
        self.variant.hgvs_genomic = self.g_query
        self.variant.hn.normalize.return_value = self.g_test

        coding = hp.parse_hgvs_variant(
            "NM_000001.1:c.100A>G"
        )

        self.validator.relevant_transcripts.return_value = [
            coding
        ]

        self.validator.myevm_t_to_g.side_effect = (
            vvhgvs.exceptions.HGVSDataNotAvailableError(
                "mapping unavailable"
            )
        )

        gap_mapper = MagicMock()
        gap_mapper.gapped_g_to_c.return_value = (
            {
                "gapped_alignment_warning": "",
                "auto_info": "",
            },
            [coding],
        )

        with patch(
            "VariantValidator.modules.mappers."
            "gapped_mapping.GapMapper",
            return_value=gap_mapper,
        ):
            result = gene_to_transcripts(
                self.variant,
                self.validator,
                {},
                self.batch_list,
            )

        self.assertTrue(result)

    def test_gap_mapper_failure_raises_transcript_mapping_error(self):
        hp = vvhgvs.parser.Parser()

        coding = hp.parse_hgvs_variant(
            "NM_000001.1:c.100A>G"
        )

        self.validator.relevant_transcripts.return_value = [
            coding
        ]

        gap_mapper = MagicMock()
        gap_mapper.gapped_g_to_c.side_effect = RuntimeError(
            "gap failure"
        )

        with patch(
            "VariantValidator.modules.mappers."
            "gapped_mapping.GapMapper",
            return_value=gap_mapper,
        ):
            with self.assertRaises(TranscriptMappingError):
                gene_to_transcripts(
                    self.variant,
                    self.validator,
                    {},
                    self.batch_list,
                )

    def test_unsupported_chromosome_build(self):
        self.validator.relevant_transcripts.return_value = []

        with patch(
            "VariantValidator.modules.mappers.seq_data.supported_for_mapping",
            return_value=False,
        ):
            result = gene_to_transcripts(
                self.variant,
                self.validator,
                {},
                self.batch_list,
            )

        self.assertTrue(result)

        self.assertIn(
            "Validation will fail if the selected chromosome reference "
            "sequence does not corresponds to selected genome build.",
            self.variant.warnings[0],
        )

    def test_supported_mapping_but_validation_fails(self):
        self.validator.relevant_transcripts.return_value = []

        self.validator.vr.validate.side_effect = [
            None,
            vvhgvs.exceptions.HGVSInvalidVariantError(
                "secondary validation failed"
            ),
        ]

        with patch(
            "VariantValidator.modules.mappers.seq_data.supported_for_mapping",
            return_value=True,
        ):
            result = gene_to_transcripts(
                self.variant,
                self.validator,
                {},
                self.batch_list,
            )

        self.assertTrue(result)

        self.assertIn(
            "secondary validation failed",
            self.variant.warnings,
        )

    def test_no_transcripts_found(self):
        self.validator.relevant_transcripts.return_value = []
        self.validator.chr_to_rsg.return_value = []

        with patch(
            "VariantValidator.modules.mappers.seq_data.supported_for_mapping",
            return_value=True,
        ):
            result = gene_to_transcripts(
                self.variant,
                self.validator,
                {},
                self.batch_list,
            )

        self.assertTrue(result)

        self.assertEqual(
            self.variant.output_type_flag,
            "intergenic",
        )

        self.assertTrue(
            any(
                warning.startswith(
                    "TranscriptIdentificationWarning:"
                )
                for warning in self.variant.warnings
            )
        )

    def test_explicit_transcript_selection_finds_nothing(self):
        self.validator.select_transcripts = "NM_000001.1"
        self.validator.relevant_transcripts.return_value = []
        self.validator.chr_to_rsg.return_value = []

        with patch(
            "VariantValidator.modules.mappers.seq_data.supported_for_mapping",
            return_value=True,
        ):
            result = gene_to_transcripts(
                self.variant,
                self.validator,
                {},
                self.batch_list,
            )

        self.assertTrue(result)

        self.assertTrue(
            any(
                warning.startswith(
                    "TranscriptSelectionError:"
                )
                for warning in self.variant.warnings
            )
        )

    def test_gapped_transcripts_are_added_to_batch(self):
        hp = vvhgvs.parser.Parser()

        coding = hp.parse_hgvs_variant(
            "NM_000001.1:c.100A>G"
        )

        self.validator.relevant_transcripts.return_value = [
            coding
        ]

        gap_mapper = MagicMock()
        gap_mapper.gapped_g_to_c.return_value = (
            {
                "gapped_alignment_warning": "",
                "auto_info": "",
            },
            [coding],
        )

        with patch(
            "VariantValidator.modules.mappers."
            "gapped_mapping.GapMapper",
            return_value=gap_mapper,
        ):
            result = gene_to_transcripts(
                self.variant,
                self.validator,
                {},
                self.batch_list,
            )

        self.assertTrue(result)
        self.assertFalse(self.variant.write)
        self.assertEqual(
            len(self.batch_list),
            1,
        )

class TestTranscriptsToGeneUnit(TestCase):

    def setUp(self):
        self.variant = MagicMock()
        self.validator = MagicMock()

        self.variant.warnings = []
        self.variant.write = True
        self.variant.primary_assembly = "GRCh38"

        self.obj = MagicMock()
        self.obj.ac = "NM_000001.1"
        self.obj.rel_ac = ""
        self.variant.hgvs_formatted = self.obj
        self.variant.quibble = self.obj
        self.variant.original = "NM_000001.1:c.1A>G"

        self.validator.select_transcripts = "all"

    def test_transcript_selection_rejects_unselected_transcript(self):
        self.validator.select_transcripts = "NM_000002.1"

        with self.assertRaises(MappersError) as cm:
            transcripts_to_gene(
                self.variant,
                self.validator,
                {"NM_000002.1": {}},
            )

        self.assertIn(
            "TranscriptSelectionError:",
            str(cm.exception),
        )
        self.assertEqual(len(self.variant.warnings), 1)
        self.assertEqual(
            self.variant.warnings,
            [
                "TranscriptSelectionError: Variant "
                f"{self.variant.hgvs_formatted} is not in the list of "
                "transcripts selected for validation NM_000002.1"
            ],
        )

    def test_transcript_selection_genomic_original_suppresses_output(self):
        self.validator.select_transcripts = "NM_000002.1"
        self.variant.original = "NC_000001.11:g.100A>G"

        result = transcripts_to_gene(
            self.variant,
            self.validator,
            {"NM_000002.1": {}},
        )

        self.assertTrue(result)
        self.assertFalse(self.variant.write)
        self.assertEqual(self.variant.warnings, [])

    def test_missing_alignment_data(self):
        self.validator.myevm_t_to_g.side_effect = (
            vvhgvs.exceptions.HGVSDataNotAvailableError(
                "Alignment is incomplete ~ test"
            )
        )

        result = transcripts_to_gene(
            self.variant,
            self.validator,
            {},
        )

        self.assertTrue(result)
        self.assertTrue(
            any(
                "Full alignment data" in warning
                for warning in self.variant.warnings
            )
        )
        self.assertTrue(
            any(
                "Universal Transcript Archive" in warning
                for warning in self.variant.warnings
            )
        )
        self.assertTrue(
            any(
                "gene2transcripts" in warning
                for warning in self.variant.warnings
            )
        )

    def test_no_relevant_genomic_mapping_options(self):
        self.validator.myevm_t_to_g.side_effect = (
            vvhgvs.exceptions.HGVSDataNotAvailableError(
                "No relevant genomic mapping options"
            )
        )

        result = transcripts_to_gene(
            self.variant,
            self.validator,
            {},
        )

        self.assertTrue(result)
        self.assertTrue(
            any(
                "No relevant genomic mapping options" in warning
                for warning in self.variant.warnings
            )
        )
        self.assertTrue(
            any(
                "Universal Transcript Archive" in warning
                for warning in self.variant.warnings
            )
        )

    def test_reference_mismatch_data_error(self):
        message = (
            "NM_000001.1:c.1A>G does not agree with reference sequence"
        )

        self.validator.myevm_t_to_g.side_effect = (
            vvhgvs.exceptions.HGVSDataNotAvailableError(message)
        )

        result = transcripts_to_gene(
            self.variant,
            self.validator,
            {},
        )

        self.assertTrue(result)
        self.assertEqual(
            self.variant.warnings,
            [message],
        )

    def test_mapping_type_error(self):
        self.validator.myevm_t_to_g.side_effect = TypeError()

        result = transcripts_to_gene(
            self.variant,
            self.validator,
            {},
        )

        self.assertTrue(result)
        self.assertEqual(len(self.variant.warnings), 2)
        self.assertIn(
            "Required information for NM_000001.1",
            self.variant.warnings[0],
        )
        self.assertIn(
            "gene2transcripts",
            self.variant.warnings[1],
        )

    def test_refseqgene_sets_reset_g_origin(self):
        self.obj.rel_ac = "NG_012345.1"

        self.validator.myevm_t_to_g.side_effect = TypeError()

        result = transcripts_to_gene(
            self.variant,
            self.validator,
            {},
        )

        self.assertTrue(result)

        self.validator.myevm_t_to_g.assert_called_once()

        self.assertTrue(
            self.validator.myevm_t_to_g.call_args.kwargs[
                "reset_g_origin"
            ]
        )

    def test_non_refseqgene_does_not_reset_g_origin(self):
        self.obj.rel_ac = ""

        self.validator.myevm_t_to_g.side_effect = TypeError()

        result = transcripts_to_gene(
            self.variant,
            self.validator,
            {},
        )

        self.assertTrue(result)

        self.assertFalse(
            self.validator.myevm_t_to_g.call_args.kwargs[
                "reset_g_origin"
            ]
        )


class TestFinalTxToMultipleGenomicUnit(TestCase):

    def setUp(self):
        self.variant = MagicMock()
        self.validator = MagicMock()

        self.tx_variant = MagicMock()
        self.tx_variant.ac = "NM_000001.1"

        self.variant.hgvs_coding = self.tx_variant
        self.variant.no_norm_evm = MagicMock()
        self.variant.hn.normalize.side_effect = lambda value: value

        self.variant.genomic_g = MagicMock()
        self.variant.genomic_g.ac = "NC_000001.11"

        self.validator.alt_aln_method = "splign"

    def test_primary_liftover_uses_only_nc_splign_mappings(self):
        self.variant.map_dat.mapping_options.return_value = [
            ("NM_000001.1", "NW_000001.1", "splign"),
            ("NM_000001.1", "NC_000002.12", "genebuild"),
            ("NM_000001.1", "NC_000001.11", "splign"),
        ]

        mapped = MagicMock()
        mapped.ac = "NC_000001.11"
        mapped.posedit.edit.type = "sub"

        self.validator.myvm_t_to_g.return_value = mapped
        self.variant.map_dat.tx_exons.return_value = [
            {"alt_strand": 1}
        ]
        self.variant.map_dat.is_gapped_map.return_value = False

        result = final_tx_to_multiple_genomic(
            self.variant,
            self.validator,
            self.tx_variant,
            liftover_level="primary",
        )

        self.assertEqual(result, [mapped])

        called_alt_chromosomes = [
            call.args[1]
            for call in self.validator.myvm_t_to_g.call_args_list
        ]

        self.assertEqual(
            called_alt_chromosomes,
            ["NC_000001.11"],
        )

    def test_all_liftover_accepts_nc_nt_nw(self):
        self.variant.map_dat.mapping_options.return_value = [
            ("NM_000001.1", "NT_000001.1", "splign"),
            ("NM_000001.1", "NW_000001.1", "splign"),
        ]

        mapped = MagicMock()
        mapped.posedit.edit.type = "ins"

        self.validator.myvm_t_to_g.return_value = mapped
        self.variant.map_dat.tx_exons.return_value = [
            {"alt_strand": 1}
        ]
        self.variant.map_dat.is_gapped_map.return_value = False

        result = final_tx_to_multiple_genomic(
            self.variant,
            self.validator,
            self.tx_variant,
            liftover_level="all",
        )

        self.assertEqual(len(result), 2)

        called_alt_chromosomes = [
            call.args[1]
            for call in self.validator.myvm_t_to_g.call_args_list
        ]

        self.assertEqual(
            called_alt_chromosomes,
            [
                "NT_000001.1",
                "NW_000001.1",
            ],
        )

    def test_old_ncbi_reference_is_skipped(self):
        self.variant.map_dat.mapping_options.return_value = [
            ("NM_000001.1", "NC_001807.4", "splign"),
        ]

        result = final_tx_to_multiple_genomic(
            self.variant,
            self.validator,
            self.tx_variant,
            liftover_level="primary",
        )

        self.assertEqual(result, [])
        self.validator.myvm_t_to_g.assert_not_called()

    def test_mapping_hgvs_error_is_ignored(self):
        self.variant.map_dat.mapping_options.return_value = [
            ("NM_000001.1", "NC_000001.11", "splign"),
        ]

        self.variant.map_dat.tx_exons.return_value = [
            {"alt_strand": 1}
        ]

        self.validator.myvm_t_to_g.side_effect = (
            vvhgvs.exceptions.HGVSError("mapping failed")
        )

        result = final_tx_to_multiple_genomic(
            self.variant,
            self.validator,
            self.tx_variant,
            liftover_level="primary",
        )

        self.assertEqual(result, [])

    def test_mapping_key_error_is_ignored(self):
        self.variant.map_dat.mapping_options.return_value = [
            ("NM_000001.1", "NC_000001.11", "splign"),
        ]

        self.variant.map_dat.tx_exons.side_effect = KeyError(
            "missing exon data"
        )

        result = final_tx_to_multiple_genomic(
            self.variant,
            self.validator,
            self.tx_variant,
            liftover_level="primary",
        )

        self.assertEqual(result, [])
        self.validator.myvm_t_to_g.assert_not_called()

    def test_string_tx_variant_is_parsed(self):
        parsed = MagicMock()
        parsed.ac = "NM_000001.1"

        self.validator.hp.parse_hgvs_variant.return_value = parsed

        self.variant.map_dat.mapping_options.return_value = []

        result = final_tx_to_multiple_genomic(
            self.variant,
            self.validator,
            "NM_000001.1:c.1A>G",
            liftover_level="primary",
        )

        self.assertEqual(result, [])
        self.validator.hp.parse_hgvs_variant.assert_called_once_with(
            "NM_000001.1:c.1A>G"
        )
        self.assertIs(self.variant.hgvs_coding, parsed)

    def test_boundary_normalisation_disables_gap_compensation(self):
        self.variant.hn.normalize.side_effect = (
            vvhgvs.exceptions.HGVSUnsupportedOperationError(
                "variant spanning exon-intron boundary"
            )
        )

        self.variant.map_dat.mapping_options.return_value = []

        result = final_tx_to_multiple_genomic(
            self.variant,
            self.validator,
            self.tx_variant,
            liftover_level="primary",
        )

        self.assertEqual(result, [])

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