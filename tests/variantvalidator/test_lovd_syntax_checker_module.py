from unittest import TestCase
from unittest.mock import MagicMock
from VariantValidator.bin import lovd_syntax_checker
import subprocess
import importlib
from unittest.mock import patch
import VariantValidator.bin.lovd_syntax_checker as checker

class TestLovdSyntaxChecker(TestCase):

    @patch("VariantValidator.bin.lovd_syntax_checker.os.path.exists", return_value=False)
    def test_missing_php_script(self, mock_exists):
        with self.assertRaises(FileNotFoundError):
            lovd_syntax_checker.run_hgvs_checker("c.100del")

    @patch("VariantValidator.bin.lovd_syntax_checker.os.path.exists", return_value=True)
    @patch("VariantValidator.bin.lovd_syntax_checker.subprocess.run")
    def test_gene_mode(self, mock_run, mock_exists):

        result = MagicMock()
        result.stdout = '[{"result":"Valid"}]'

        meta = MagicMock()
        meta.stdout = '[{"library_version":"1.2.3"}]'

        mock_run.side_effect = [result, meta]

        output = lovd_syntax_checker.run_hgvs_checker(
            "BRCA1",
            is_a_gene=True
        )

        self.assertEqual(output[0]["result"], "Valid")

        command = mock_run.call_args_list[0][0][0]
        self.assertIn("gene:BRCA1", command)

    @patch("VariantValidator.bin.lovd_syntax_checker.time.sleep")
    @patch("VariantValidator.bin.lovd_syntax_checker.os.path.exists", return_value=True)
    @patch("VariantValidator.bin.lovd_syntax_checker.subprocess.run")
    def test_subprocess_retry_failure(
        self,
        mock_run,
        mock_exists,
        mock_sleep,
    ):

        mock_run.side_effect = subprocess.CalledProcessError(
            1,
            "cmd",
            stderr="boom"
        )

        with self.assertRaises(RuntimeError):
            lovd_syntax_checker.run_hgvs_checker("c.100del")

        self.assertEqual(mock_sleep.call_count, 2)

    @patch("VariantValidator.bin.lovd_syntax_checker.time.sleep")
    @patch("VariantValidator.bin.lovd_syntax_checker.os.path.exists", return_value=True)
    @patch("VariantValidator.bin.lovd_syntax_checker.subprocess.run")
    def test_invalid_json_retry_failure(
        self,
        mock_run,
        mock_exists,
        mock_sleep,
    ):

        result = MagicMock()
        result.stdout = "not json"

        meta = MagicMock()
        meta.stdout = '[{"library_version":"1.2.3"}]'

        mock_run.side_effect = [
            result, meta,
            result, meta,
            result, meta,
        ]

        with self.assertRaises(RuntimeError):
            lovd_syntax_checker.run_hgvs_checker("c.100del")

    @patch("VariantValidator.bin.lovd_syntax_checker.time.sleep")
    @patch("VariantValidator.bin.lovd_syntax_checker.os.path.exists", return_value=True)
    @patch("VariantValidator.bin.lovd_syntax_checker.subprocess.run")
    def test_unexpected_exception_retry_failure(
        self,
        mock_run,
        mock_exists,
        mock_sleep,
    ):

        mock_run.side_effect = OSError("unexpected")

        with self.assertRaises(RuntimeError):
            lovd_syntax_checker.run_hgvs_checker("c.100del")

    @patch("VariantValidator.bin.lovd_syntax_checker.pkg_resources.files")
    def test_get_installation_path_module_not_found(self, mock_files):
        mock_files.side_effect = ModuleNotFoundError()

        with self.assertRaises(FileNotFoundError):
            lovd_syntax_checker.get_installation_path()

    @patch("configparser.ConfigParser")
    def test_php_path_from_config(self, mock_parser):
        config = {
            "php": {
                "php_path": "/custom/php"
            }
        }

        mock_instance = mock_parser.return_value
        mock_instance.read.return_value = None
        mock_instance.__getitem__.return_value.get.return_value = "/custom/php"
        mock_instance.__getitem__.return_value.__getitem__.return_value = "/custom/php"

        reloaded = importlib.reload(checker)

        self.assertEqual(reloaded.PHP_PATH, "/custom/php")

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