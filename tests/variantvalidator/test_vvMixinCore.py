import unittest
from unittest import TestCase
from unittest.mock import patch

import VariantValidator
from VariantValidator import Validator
from VariantValidator.modules.vvMixinCore import ValidatorSubmissionError
import vvhgvs


class TestTranscriptInfoFunctional(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vv = Validator()

    def _wrap_get_transcript_info(self, original):

        def wrapper(my_variant):
            print("\n==============================")
            print("_get_transcript_info() CALLED")
            print("variant:", my_variant.quibble)

            result = original(my_variant)

            print("_get_transcript_info() RETURNED", result)
            print("warnings:", my_variant.warnings)
            print("==============================\n")

            return result

        return wrapper

    @patch("VariantValidator.modules.vvDatabase.Database.in_entries")
    def test_cached_transcript_record(self, mock_in_entries):

        mock_in_entries.return_value = {
            "accession": "NM_000546.6",
            "description": "Tumour protein p53",
            "variant": "",
            "version": "NM_000546.6",
            "hgnc_symbol": "TP53",
            "uta_symbol": "TP53",
            "updated": "2026-01-01",
            "expiry": "false",
        }

        original = self.vv._get_transcript_info

        with patch.object(
            self.vv,
            "_get_transcript_info",
            wraps=self._wrap_get_transcript_info(original),
        ) as wrapped:

            results = self.vv.validate(
                "NM_000546.6:c.215C>G",
                "GRCh38",
                "all",
            ).format_as_dict(test=True)

        print("in_entries call count =", mock_in_entries.call_count)
        print(mock_in_entries.call_args_list)

        self.assertTrue(wrapped.called)
        self.assertEqual(results["flag"], "gene_variant")

    @patch("VariantValidator.modules.vvDatabase.Database.in_entries")
    def test_database_error(self, mock_in_entries):

        mock_in_entries.return_value = {
            "error": "error",
            "description": "Database exploded",
        }

        original = self.vv._get_transcript_info

        with patch.object(
            self.vv,
            "_get_transcript_info",
            wraps=self._wrap_get_transcript_info(original),
        ) as wrapped:

            results = self.vv.validate(
                "NM_000546.6:c.215C>G",
                "GRCh38",
                "all",
            ).format_as_dict(test=True)

        print("in_entries call count =", mock_in_entries.call_count)
        print(mock_in_entries.call_args_list)

        self.assertTrue(wrapped.called)

        warnings = results["validation_warning_1"]["validation_warnings"]

        self.assertIn("Database exploded", warnings)

    @patch("VariantValidator.modules.vvDatabase.Database.data_add")
    @patch("VariantValidator.modules.vvDatabase.Database.in_entries")
    def test_expired_record(self, mock_in_entries, mock_data_add):

        mock_in_entries.return_value = {
            "accession": "NM_000546.6",
            "description": "Old description",
            "variant": "",
            "version": "NM_000546.6",
            "hgnc_symbol": "TP53",
            "uta_symbol": "TP53",
            "updated": "2024-01-01",
            "expiry": "true",
        }

        mock_data_add.return_value = {
            "accession": "NM_000546.6",
            "description": "Updated description",
            "variant": "",
            "version": "NM_000546.6",
            "hgnc_symbol": "TP53",
            "uta_symbol": "TP53",
            "updated": "2026-01-01",
            "expiry": "false",
        }

        original = self.vv._get_transcript_info

        with patch.object(
            self.vv,
            "_get_transcript_info",
            wraps=self._wrap_get_transcript_info(original),
        ) as wrapped:

            self.vv.validate(
                "NM_000546.6:c.215C>G",
                "GRCh38",
                "all",
            )

        print("in_entries call count =", mock_in_entries.call_count)
        print(mock_in_entries.call_args_list)
        print("data_add call count =", mock_data_add.call_count)
        print(mock_data_add.call_args_list)

        self.assertTrue(wrapped.called)

    @patch("VariantValidator.modules.vvDatabase.Database.update_transcript_info_record")
    @patch("VariantValidator.modules.vvDatabase.Database.in_entries")
    def test_missing_record(self, mock_in_entries, mock_update):

        mock_in_entries.side_effect = [
            {
                "none": True,
            },
            {
                "accession": "NM_000546.6",
                "description": "Recovered description",
                "variant": "",
                "version": "NM_000546.6",
                "hgnc_symbol": "TP53",
                "uta_symbol": "TP53",
                "updated": "2026-01-01",
                "expiry": "false",
            },
        ]

        mock_update.return_value = None

        original = self.vv._get_transcript_info

        with patch.object(
            self.vv,
            "_get_transcript_info",
            wraps=self._wrap_get_transcript_info(original),
        ) as wrapped:

            self.vv.validate(
                "NM_000546.6:c.215C>G",
                "GRCh38",
                "all",
            )

        print("in_entries call count =", mock_in_entries.call_count)
        print(mock_in_entries.call_args_list)
        print("update_transcript_info_record call count =", mock_update.call_count)
        print(mock_update.call_args_list)

        self.assertTrue(wrapped.called)

    @patch("VariantValidator.modules.vvDatabase.Database.in_entries")
    @patch("VariantValidator.modules.vvDatabase.Database.update_transcript_info_record")
    def test_missing_record_populated_after_update(self, mock_update, mock_in_entries):

        mock_in_entries.side_effect = [
            {
                "none": "none",
            },
            {
                "accession": "NM_000546.6",
                "description": "Updated TP53",
                "variant": "",
                "version": "NM_000546.6",
                "hgnc_symbol": "TP53",
                "uta_symbol": "TP53",
                "updated": "2026-01-01",
                "expiry": "false",
            },
        ]

        results = self.vv.validate(
            "NM_000546.6:c.215C>G",
            "GRCh38",
            "all",
        ).format_as_dict(test=True)

        mock_update.assert_called_once()
        self.assertEqual(results["flag"], "gene_variant")


    @patch("VariantValidator.modules.vvDatabase.Database.in_entries")
    @patch("VariantValidator.modules.vvDatabase.Database.update_transcript_info_record")
    def test_missing_record_update_failure(self, mock_update, mock_in_entries):

        mock_in_entries.return_value = {
            "none": "none",
        }

        mock_update.side_effect = Exception("Ensembl unavailable")

        results = self.vv.validate(
            "NM_000546.6:c.215C>G",
            "GRCh38",
            "all",
        ).format_as_dict(test=True)

        warnings = results["validation_warning_1"]["validation_warnings"]

        self.assertTrue(
            any(
                "Unable to assign transcript identity records"
                in warning
                for warning in warnings
            )
        )


    @patch("VariantValidator.modules.vvDatabase.Database.in_entries")
    @patch("VariantValidator.modules.vvDatabase.Database.data_add")
    def test_expired_record_refresh_hgvs_error(self, mock_data_add, mock_in_entries):

        mock_in_entries.return_value = {
            "accession": "NM_000546.6",
            "description": "Old description",
            "variant": "",
            "version": "NM_000546.6",
            "hgnc_symbol": "TP53",
            "uta_symbol": "TP53",
            "updated": "2024-01-01",
            "expiry": "true",
        }

        mock_data_add.side_effect = vvhgvs.exceptions.HGVSError(
            "Transcript unsupported"
        )

        results = self.vv.validate(
            "NM_000546.6:c.215C>G",
            "GRCh38",
            "all",
        ).format_as_dict(test=True)

        warnings = results["validation_warning_1"]["validation_warnings"]

        self.assertTrue(
            any(
                "not currently supported"
                in warning
                for warning in warnings
            )
        )


    @patch("VariantValidator.modules.vvDatabase.Database.in_entries")
    def test_cached_record_without_hgnc_symbol(self, mock_in_entries):

        mock_in_entries.return_value = {
            "accession": "NM_000546.6",
            "description": "Tumour protein p53",
            "variant": "",
            "version": "NM_000546.6",
            "hgnc_symbol": "",
            "uta_symbol": "",
            "updated": "2026-01-01",
            "expiry": "false",
        }

        results = self.vv.validate(
            "NM_000546.6:c.215C>G",
            "GRCh38",
            "all",
        ).format_as_dict(test=True)

        self.assertEqual(results["flag"], "gene_variant")


    @patch("VariantValidator.modules.vvDatabase.Database.in_entries")
    def test_cached_record_empty_description(self, mock_in_entries):

        mock_in_entries.return_value = {
            "accession": "NM_000546.6",
            "description": "",
            "variant": "",
            "version": "NM_000546.6",
            "hgnc_symbol": "TP53",
            "uta_symbol": "TP53",
            "updated": "2026-01-01",
            "expiry": "false",
        }

        results = self.vv.validate(
            "NM_000546.6:c.215C>G",
            "GRCh38",
            "all",
        ).format_as_dict(test=True)

        self.assertEqual(results["flag"], "gene_variant")

    def test_validate_keyword_arguments(self):

        results = self.vv.validate(
            variant="NM_000546.6:c.215C>G",
            genome="GRCh38",
            select_transcripts="all",
        ).format_as_dict(test=True)

        self.assertEqual(results["flag"], "gene_variant")

    def test_validate_keyword_arguments_ensembl(self):

        results = self.vv.validate(
            variant="ENST00000269305.9:c.215C>G",
            genome="GRCh38",
            select_transcripts="all",
            transcript_set="ensembl",
        ).format_as_dict(test=True)

        self.assertEqual(results["flag"], "gene_variant")

    def test_validate_missing_variant_keyword(self):

        with self.assertRaises(ValidatorSubmissionError) as err:
            self.vv.validate(
                genome="GRCh38",
            )

        self.assertEqual(
            str(err.exception),
            "ValidatorSubmissionError: No variant descriptions submitted.",
        )

    def test_validate_missing_genome_keyword(self):

        with self.assertRaises(ValidatorSubmissionError) as err:
            self.vv.validate(
                variant="NM_000546.6:c.215C>G",
            )

        self.assertEqual(
            str(err.exception),
            "ValidatorSubmissionError: No genome build submitted.",
        )

class TestValidator(unittest.TestCase):
    """
    Going to test the Validator function with a series of different inputs/situations that aren't covered in
    test_inputs.py
    """

    @classmethod
    def setUpClass(cls):
        cls.vv = VariantValidator.Validator()

    def test_transcript_seq_nonsense(self):
        var = 'NM_015120.4:c.34C>T'
        with self.assertRaises(Exception):
            self.vv.validate(var, 'GRCh37', 'all', transcript_set='nonsense')

    # def test_transcript_seq_ensembl(self):
    #     var = 'NM_015120.4:c.34C>T'
    #     with self.assertRaises(Exception):
    #         self.vv.validate(var, 'GRCh37', 'all', transcript_set='ensembl')
    #
    #     self.assertEqual(self.vv.alt_aln_method, 'genebuild')

    def test_transcript_list(self):
        var = 'NM_015120.4:c.34C>T'

        output = self.vv.validate(var, 'GRCh37', 'Trans1').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'warning')
        assert "TranscriptSelectionError" in str(output["validation_warning_1"]["validation_warnings"])

    def test_transcript_list_realid(self):
        var = 'NM_015120.4:c.34C>T'

        output = self.vv.validate(var, 'GRCh37', 'NM_015120.4').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'gene_variant')
        self.assertEqual(list(output), ['flag', 'NM_015120.4:c.34C>T', 'metadata'])

    def test_transcript_list_real_pair(self):
        var = 'NM_015120.4:c.34C>T'

        output = self.vv.validate(var, 'GRCh37', '["NM_015120.4", "NM_015120.5"]').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'gene_variant')
        self.assertEqual(list(output), ['flag', 'NM_015120.4:c.34C>T', 'metadata'])

    def test_transcript_list_lrg(self):
        var = 'NM_015120.4:c.34C>T'

        output = self.vv.validate(var, 'GRCh37', 'LRG1').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'warning')
        assert "TranscriptSelectionError" in str(output["validation_warning_1"]["validation_warnings"])

    def test_non_ascii(self):
        var = 'NM_015120.4:c.34C>T\202'

        output = self.vv.validate(var, 'GRCh37', 'all').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'warning')
        self.assertIn('VariantSyntaxError: Submitted variant description contains an invalid character',
                      str(output['validation_warning_1']['validation_warnings']))

    def test_assembly_hg19(self):
        var = 'NM_015120.4:c.34C>T'

        out = self.vv.validate(var, 'hg19', 'all')
        for variant in out.output_list:
            self.assertEqual(variant.primary_assembly, 'GRCh37')
        output = out.format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'gene_variant')
        self.assertEqual(list(output), ['flag', 'NM_015120.4:c.34C>T', 'metadata'])

    def test_assembly_hg38(self):
        var = 'NM_015120.4:c.34C>T'

        out = self.vv.validate(var, 'hg38', 'all')
        for variant in out.output_list:
            self.assertEqual(variant.primary_assembly, 'GRCh38')
        output = out.format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'gene_variant')
        self.assertEqual(list(output), ['flag', 'NM_015120.4:c.34C>T', 'metadata'])

    def test_assembly_grch(self):
        var = 'NM_015120.4:c.34C>T'

        out = self.vv.validate(var, 'grch37', 'all')
        for variant in out.output_list:
            self.assertEqual(variant.primary_assembly, 'GRCh37')
        output = out.format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'gene_variant')
        self.assertEqual(list(output), ['flag', 'NM_015120.4:c.34C>T', 'metadata'])

    def test_assembly_invalid(self):
        var = 'NM_015120.4:c.34C>T'

        out = self.vv.validate(var, 'nonsense', 'all')
        for variant in out.output_list:
            self.assertEqual(variant.primary_assembly, 'GRCh38')
        output = out.format_as_dict()
        self.assertEqual(output['flag'], 'gene_variant')
        self.assertEqual(list(output), ['flag', 'NM_015120.4:c.34C>T', 'metadata'])
        self.assertIn('Invalid genome build has been specified',
                      str(output['NM_015120.4:c.34C>T']['validation_warnings']))

    def test_variant_invalid(self):
        var = 'NM_015120.4c.34C>T'

        output = self.vv.validate(var, 'GRCh37', 'all').format_as_dict()
        print(output)
        self.assertIn('Unable to identify a colon (:) in the variant description',
                      str(output['NM_015120.4:c.34C>T']['validation_warnings']))

    def test_variant_invalid_2(self):
        var = 'NM_015120.4:c34C>T'

        output = self.vv.validate(var, 'GRCh37', 'all').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'warning')
        self.assertIn('Unable to identify a dot (.) in the variant description NM_015120.4:c34C>T following the '
                      'reference sequence type (g,c,n,r, or p). A dot is required in HGVS variant descriptions to '
                      'separate the reference type from the variant position i.e. <accession>:<type>. e.g. :g.',
                      str(output['validation_warning_1']['validation_warnings']))

    def test_variant_con(self):
        var = 'NM_015120.4:c.34con'

        output = self.vv.validate(var, 'GRCh37', 'all').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'warning')
        self.assertIn('Conversions are no longer valid HGVS Sequence Variant Descriptions',
                      str(output['validation_warning_1']['validation_warnings']))

    def test_variant_description(self):
        var = 'NM_015120.4:c.34C>T'

        out = self.vv.validate(var, 'grch37', 'all').format_as_dict()
        self.assertNotEqual(out['NM_015120.4:c.34C>T']['transcript_description'], 'false')
        self.assertEqual(out['NM_015120.4:c.34C>T']['transcript_description'],
                         'Homo sapiens ALMS1 centrosome and basal body associated protein (ALMS1), transcript variant 1, mRNA')

    def test_variant_format(self):
        var = "NM_020812.3:c.[3190_3191delCT];[(3190_3191delCT)]"

        out = self.vv.validate(var, 'grch37', 'all').format_as_dict()
        self.assertEqual(out['flag'], 'warning')
        self.assertEqual(out['validation_warning_1']['validation_warnings'],
                         ['Unsupported format c.[3190_3191delCT];[(3190_3191delCT)]'])

    def test_variant_quotes_start(self):
        var = '"NM_015120.4:c.34C>T'

        out = self.vv.validate(var, 'GRCh37', 'all').format_as_dict()
        self.assertEqual(out['flag'], 'gene_variant')
        self.assertTrue('NM_015120.4:c.34C>T' in out.keys())

    def test_variant_quotes_end(self):
        var = 'NM_015120.4:c.34C>T"'

        out = self.vv.validate(var, 'GRCh37', 'all').format_as_dict()
        print(out)
        self.assertEqual(out['flag'], 'gene_variant')
        self.assertTrue('NM_015120.4:c.34C>T' in out.keys())

    def test_variant_quotes_both(self):
        var = '["NM_015120.4:c.34C>T"]'

        out = self.vv.validate(var, 'GRCh37', 'all').format_as_dict()
        self.assertEqual(out['flag'], 'gene_variant')
        self.assertTrue('NM_015120.4:c.34C>T' in out.keys())


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
