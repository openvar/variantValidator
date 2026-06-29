import requests
from unittest import TestCase
from unittest.mock import patch, MagicMock
from VariantValidator.modules import lovd_api

class TestLOVDApi(TestCase):

    @patch("VariantValidator.modules.lovd_api.lovd_syntax_checker.run_hgvs_checker")
    def test_run_lovd_checker_cli_success(self, mock_run_hgvs_checker):
        """Test CLI checker returning valid data."""
        mock_run_hgvs_checker.return_value = [{
            "metadata": {"library_version": "1.2.3"},
            "result": "Valid"
        }]

        variant = "c.100del"
        expected_output = {
            "data": [{
                "metadata": {"library_version": "1.2.3"},
                "result": "Valid"
            }],
            "url": f"https://api.lovd.nl/v2/checkHGVS/{variant}",
            "version": "1.2.3"
        }

        result = lovd_api.run_lovd_checker_cli(variant)
        self.assertEqual(result, expected_output)

    @patch("VariantValidator.modules.lovd_api.lovd_syntax_checker.run_hgvs_checker")
    def test_run_lovd_checker_cli_failure(self, mock_run_hgvs_checker):
        """Test CLI checker returning an error."""
        mock_run_hgvs_checker.side_effect = Exception("CLI error")

        variant = "c.100del"
        expected_output = {"lovd_api_error": "CLI check failed: CLI error"}

        result = lovd_api.run_lovd_checker_cli(variant)
        self.assertEqual(result, expected_output)

    @patch("requests.get")
    def test_run_lovd_checker_web_success(self, mock_get):
        """Test web API checker returning valid data."""
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "versions": {"library_version": "1.2.3"},
            "result": "Valid"
        }
        mock_response.status_code = 200
        mock_get.return_value = mock_response

        variant = "c.100del"
        expected_output = {
            "versions": {"library_version": "1.2.3"},
            "result": "Valid",
            "url": f"https://api.lovd.nl/v2/checkHGVS/{variant}",
            "version": "1.2.3"
        }

        result = lovd_api.run_lovd_checker_web(variant)
        self.assertEqual(result, expected_output)

    @patch("requests.get")
    def test_run_lovd_checker_web_failure(self, mock_get):
        """Test web API checker failing with a request error."""
        mock_get.side_effect = requests.RequestException("Network error")

        variant = "c.100del"
        expected_output = {"lovd_api_error": "Request failed: Network error"}

        result = lovd_api.run_lovd_checker_web(variant)
        self.assertEqual(result, expected_output)

    @patch("requests.get")
    def test_run_lovd_checker_web_invalid_json(self, mock_get):
        """Test web API checker failing with invalid JSON."""
        mock_response = MagicMock()
        mock_response.json.side_effect = ValueError("Invalid JSON")
        mock_get.return_value = mock_response

        variant = "c.100del"
        expected_output = {"lovd_api_error": "Unexpected error: Invalid JSON"}

        result = lovd_api.run_lovd_checker_web(variant)
        self.assertEqual(result, expected_output)

    @patch("VariantValidator.modules.lovd_api.run_lovd_checker_web")
    @patch("VariantValidator.modules.lovd_api.run_lovd_checker_cli")
    def test_lovd_syntax_check_cli_success_does_not_call_web(
        self,
        mock_run_cli,
        mock_run_web,
    ):
        """A successful CLI call should not fall back to the web API."""
        mock_run_cli.return_value = {
            "result": "Valid"
        }

        result = lovd_api.lovd_syntax_check("c.100del")

        mock_run_cli.assert_called_once_with(
            "c.100del",
            is_a_gene=False,
        )
        mock_run_web.assert_not_called()

        self.assertEqual(
            result,
            {
                "result": "Valid"
            },
        )

    @patch("VariantValidator.modules.lovd_api.lovd_syntax_checker.run_hgvs_checker")
    def test_run_lovd_checker_cli_gene_success(self, mock_run_hgvs_checker):
        """Test CLI checker for gene symbols."""
        mock_run_hgvs_checker.return_value = [{
            "metadata": {"library_version": "1.2.3"},
            "result": "Valid"
        }]

        result = lovd_api.run_lovd_checker_cli(
            "BRCA1",
            is_a_gene=True,
        )

        self.assertEqual(
            result["url"],
            "https://api.lovd.nl/v2/checkGene/BRCA1",
        )
        self.assertEqual(result["version"], "1.2.3")


    @patch("requests.get")
    def test_run_lovd_checker_web_gene_not_supported(self, mock_get):
        """Gene symbols are not currently supported by the web API."""
        result = lovd_api.run_lovd_checker_web(
            "BRCA1",
            is_a_gene=True,
        )

        self.assertEqual(
            result,
            {
                "lovd_api_error":
                    "Unsupported value: Web API is currently not configured to support gene symbols"
            },
        )
        mock_get.assert_not_called()


    @patch("VariantValidator.modules.lovd_api.run_lovd_checker_web")
    @patch("VariantValidator.modules.lovd_api.run_lovd_checker_cli")
    def test_lovd_syntax_check_gene_fallback(
        self,
        mock_run_cli,
        mock_run_web,
    ):
        """CLI failure for a gene should fall back to the web handler."""
        mock_run_cli.return_value = {
            "lovd_api_error": "CLI failed"
        }

        mock_run_web.return_value = {
            "lovd_api_error":
                "Unsupported value: Web API is currently not configured to support gene symbols"
        }

        result = lovd_api.lovd_syntax_check(
            "BRCA1",
            is_a_gene=True,
        )

        mock_run_cli.assert_called_once_with(
            "BRCA1",
            is_a_gene=True,
        )
        mock_run_web.assert_called_once_with(
            "BRCA1",
            is_a_gene=True,
        )

        self.assertEqual(
            result,
            {
                "lovd_api_error":
                    "Unsupported value: Web API is currently not configured to support gene symbols"
            },
        )


    def test_remove_double_quotes_scalar(self):
        """Objects that are not containers should be returned unchanged."""
        self.assertEqual(lovd_api.remove_double_quotes(123), 123)
        self.assertEqual(lovd_api.remove_double_quotes(None), None)
        self.assertEqual(lovd_api.remove_double_quotes(True), True)

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