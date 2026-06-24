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
    """Test subprocess.CalledProcessError is caught, retried, and warning printed."""

    with patch(
        "VariantValidator.bin.setup_lovd_syntax_checker.subprocess.run",
        side_effect=subprocess.CalledProcessError(1, "php"),
    ), patch(
        "builtins.print"
    ) as mock_print, patch(
        "VariantValidator.bin.setup_lovd_syntax_checker.os.makedirs"
    ), patch(
        "VariantValidator.bin.setup_lovd_syntax_checker.urllib.request.urlretrieve"
    ), patch(
        "VariantValidator.bin.setup_lovd_syntax_checker.os.chmod"
    ):
        reload_setup_module.setup_lovd()

        # Assert retry message is printed
        mock_print.assert_any_call(
            "PHP cache update failed (likely due to PHP memory limits). "
            "Retrying with increased PHP memory..."
        )

        # Assert final failure warning is printed
        mock_print.assert_any_call(
            "Failed to update the LOVD HGVS cache even after increasing "
            "PHP memory limit.\n"
            "Please check your PHP installation."
        )

def test_get_installation_path_module_not_found(reload_setup_module):
    with patch(
        "VariantValidator.bin.setup_lovd_syntax_checker.pkg_resources.files",
        side_effect=ModuleNotFoundError(),
    ), patch(
        "VariantValidator.bin.setup_lovd_syntax_checker.sys.exit",
        side_effect=SystemExit,
    ) as mock_exit, patch(
        "builtins.print"
    ) as mock_print:

        with pytest.raises(SystemExit):
            reload_setup_module.get_installation_path()

        mock_print.assert_called_with(
            "Error: VariantValidator package is not installed."
        )
        mock_exit.assert_called_once_with(1)


def test_cache_permission_warning(reload_setup_module):
    with patch(
        "VariantValidator.bin.setup_lovd_syntax_checker.os.makedirs"
    ), patch(
        "VariantValidator.bin.setup_lovd_syntax_checker.download_file"
    ), patch(
        "VariantValidator.bin.setup_lovd_syntax_checker.run_php_cache_update"
    ), patch(
        "VariantValidator.bin.setup_lovd_syntax_checker.os.chmod",
        side_effect=Exception("permission denied"),
    ), patch(
        "builtins.print"
    ) as mock_print:

        reload_setup_module.setup_lovd()

        mock_print.assert_any_call(
            f"Warning: could not change permissions on "
            f"{reload_setup_module.CACHE_DIR}: permission denied"
        )


def test_run_php_cache_update_prints_details(reload_setup_module):
    error = subprocess.CalledProcessError(1, "php")

    with patch(
        "VariantValidator.bin.setup_lovd_syntax_checker.subprocess.run",
        side_effect=[error, error],
    ), patch(
        "builtins.print"
    ) as mock_print:

        reload_setup_module.run_php_cache_update()

        mock_print.assert_any_call(
            "Failed to update the LOVD HGVS cache even after increasing "
            "PHP memory limit.\n"
            "Please check your PHP installation."
        )

        mock_print.assert_any_call(f"Details: {error}")

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
