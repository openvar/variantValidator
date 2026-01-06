import os
import subprocess
from unittest.mock import patch, MagicMock
import importlib
import pytest

MODULE = "VariantValidator.bin.setup_lovd_syntax_checker"


@pytest.fixture
def reload_setup_module():
    """Reload the module fresh for each test."""
    if MODULE in importlib.sys.modules:
        del importlib.sys.modules[MODULE]
    return importlib.import_module(MODULE)


def test_setup_lovd_creates_directories_and_downloads(reload_setup_module):
    """Test that directories are created and files downloaded, chmod applied, subprocess run."""

    # Patch side-effect functions
    with patch("VariantValidator.bin.setup_lovd_syntax_checker.os.makedirs") as mock_makedirs, \
         patch("VariantValidator.bin.setup_lovd_syntax_checker.urllib.request.urlretrieve") as mock_urlretrieve, \
         patch("VariantValidator.bin.setup_lovd_syntax_checker.os.chmod") as mock_chmod, \
         patch("VariantValidator.bin.setup_lovd_syntax_checker.subprocess.run") as mock_run:

        # Make subprocess.run succeed
        mock_run.return_value = MagicMock(returncode=0)

        # Reload module to apply mocks
        importlib.reload(reload_setup_module)
        reload_setup_module.setup_lovd()

        # Check directories created
        LOVD_DIR = reload_setup_module.LOVD_DIR
        WEB_DIR = reload_setup_module.WEB_DIR
        CACHE_DIR = reload_setup_module.CACHE_DIR

        mock_makedirs.assert_any_call(LOVD_DIR, exist_ok=True)
        mock_makedirs.assert_any_call(WEB_DIR, exist_ok=True)
        mock_makedirs.assert_any_call(CACHE_DIR, exist_ok=True)

        # Check files downloaded
        mock_urlretrieve.assert_any_call(reload_setup_module.HGVS_CHECKER_URL, reload_setup_module.PHP_SCRIPT)
        mock_urlretrieve.assert_any_call(reload_setup_module.AJAX_URL, reload_setup_module.AJAX_SCRIPT)
        mock_urlretrieve.assert_any_call(reload_setup_module.INDEX_URL, reload_setup_module.INDEX_SCRIPT)
        mock_urlretrieve.assert_any_call(reload_setup_module.UPDATE_URL, reload_setup_module.UPDATE_SCRIPT)

        # Check chmod applied to files and cache dir
        mock_chmod.assert_any_call(reload_setup_module.PHP_SCRIPT, 0o755)
        mock_chmod.assert_any_call(reload_setup_module.CACHE_DIR, 0o777)

        # Check PHP subprocess called
        mock_run.assert_called_with(["php", "-f", reload_setup_module.UPDATE_SCRIPT], check=True)


def test_download_file_failure(reload_setup_module):
    """Test download_file prints error on exception."""

    with patch("VariantValidator.bin.setup_lovd_syntax_checker.urllib.request.urlretrieve", side_effect=Exception("fail")), \
         patch("builtins.print") as mock_print:
        reload_setup_module.download_file("http://example.com/file.php", "/tmp/file.php")
        mock_print.assert_any_call("Failed to download http://example.com/file.php: fail")


def test_subprocess_failure_warning(reload_setup_module):
    """Test subprocess.CalledProcessError is caught and warning printed."""

    with patch("VariantValidator.bin.setup_lovd_syntax_checker.subprocess.run", side_effect=subprocess.CalledProcessError(1, "php")), \
         patch("builtins.print") as mock_print, \
         patch("VariantValidator.bin.setup_lovd_syntax_checker.os.makedirs"), \
         patch("VariantValidator.bin.setup_lovd_syntax_checker.urllib.request.urlretrieve"), \
         patch("VariantValidator.bin.setup_lovd_syntax_checker.os.chmod"):
        reload_setup_module.setup_lovd()
        mock_print.assert_any_call("Failed to run PHP cache update: Command 'php' returned non-zero exit status 1.")

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
