from unittest import TestCase
from unittest.mock import patch

from VariantValidator import Validator
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
