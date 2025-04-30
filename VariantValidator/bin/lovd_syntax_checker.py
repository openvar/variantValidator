import os
import subprocess
import json
import importlib.resources as pkg_resources

def get_installation_path():
    """Determine the correct installation path inside the installed VariantValidator package."""
    try:
        package_path = str(pkg_resources.files("VariantValidator"))  # Get installed package path
        return os.path.join(package_path, "php", "lovd_hgvs_checker")
    except ModuleNotFoundError:
        raise FileNotFoundError("Error: VariantValidator package is not installed.")

LOVD_DIR = get_installation_path()
PHP_SCRIPT = os.path.join(LOVD_DIR, "HGVS.php")

def run_hgvs_checker(variant):
    """Run the LOVD HGVS Syntax Checker with the given variant using PHP."""
    if not os.path.exists(PHP_SCRIPT):
        raise FileNotFoundError(
            f"HGVS Syntax Checker is not installed. Expected at: {PHP_SCRIPT}"
        )

    try:
        result = subprocess.run(
            ["php", "-f", PHP_SCRIPT, variant],
            capture_output=True,
            text=True,
            check=True
        )

        result_meta = subprocess.run(
            ["php", "-f", PHP_SCRIPT, "getVersions"],
            capture_output=True,
            text=True,
            check=True
        )
        json_result = json.loads(result.stdout.strip())
        json_meta = json.loads(result_meta.stdout.strip())
        json_result[0]["metadata"] = json_meta[0]
        return json_result # Parse JSON output
    except subprocess.CalledProcessError as e:
        raise f"Error running HGVS Checker: {e.stderr.strip()}"
    except json.JSONDecodeError:
        raise f"Invalid JSON output: {result.stdout.strip()}"
    except Exception:
        raise

if __name__ == "__main__":
    test_variant = "c.100del"
    output = run_hgvs_checker(test_variant)
    print(output)
    test_variant = "Dmd\\"
    output = run_hgvs_checker(test_variant)
    print(output)


# <LICENSE>
# Copyright (C) 2016-2025 VariantValidator Contributors
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