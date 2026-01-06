import requests
from unittest import TestCase
from unittest.mock import patch, MagicMock
from VariantValidator.modules import lovd_api
from VariantValidator.bin import lovd_syntax_checker


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

    @patch("VariantValidator.modules.lovd_api.run_lovd_checker_cli", side_effect=Exception("CLI failure"))
    @patch("VariantValidator.modules.lovd_api.run_lovd_checker_web")
    def test_lovd_syntax_check_fallback(self, mock_run_web, mock_run_cli):
        """Test that the syntax check falls back to the web API if CLI fails."""
        mock_run_web.return_value = {
            "versions": {"library_version": "1.2.3"},
            "result": "Valid"
        }

        variant = "c.100del"
        expected_output = {
            "versions": {"library_version": "1.2.3"},
            "result": "Valid"
        }

        result = lovd_api.lovd_syntax_check(variant)
        self.assertEqual(result, expected_output)

    def test_lovd_syntax_check_disabled(self):
        """Test when the LOVD check is explicitly disabled."""
        variant = "c.100del"
        expected_output = {"lovd_api_error": "Do LOVD syntax check set to False"}

        result = lovd_api.lovd_syntax_check(variant, do_lovd_check=False)
        self.assertEqual(result, expected_output)

    def test_remove_double_quotes(self):
        """Test the removal of double quotes from different data types."""
        input_data = {
            "string": '"Hello"',
            "list": ['"Hello"', '"World"'],
            "tuple": ('"Hello"', '"World"'),
            "set": {'"Hello"', '"World"'},
            "dict": {"key": '"Value"'}
        }

        expected_output = {
            "string": "Hello",
            "list": ["Hello", "World"],
            "tuple": ("Hello", "World"),
            "set": {"Hello", "World"},
            "dict": {"key": "Value"}
        }

        result = lovd_api.remove_double_quotes(input_data)
        self.assertEqual(result, expected_output)

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