import subprocess
import json
import time
import threading
import importlib.resources as pkg_resources
import os
from configparser import ConfigParser

CONFIG_DIR = os.path.join(os.path.expanduser('~'), '.variantvalidator')

# Change settings based on config
config = ConfigParser()
config.read(CONFIG_DIR)

# Set PHP path if specified in config
PHP_PATH = None
try:
    if config['php']['php_path']:
        PHP_PATH = config['php']['php_path']
except Exception:
    pass

# Lock to ensure thread-safe subprocess access
lovd_cli_lock = threading.Lock()

def get_installation_path():
    """Determine the correct installation path inside the installed VariantValidator package."""
    try:
        package_path = str(pkg_resources.files("VariantValidator"))  # Get installed package path
        return os.path.join(package_path, "php", "lovd_hgvs_checker")
    except ModuleNotFoundError:
        raise FileNotFoundError("Error: VariantValidator package is not installed.")

LOVD_DIR = get_installation_path()
PHP_SCRIPT = os.path.join(LOVD_DIR, "HGVS.php")

MAX_RETRIES = 3
RETRY_DELAY = 0.2  # 200 milliseconds

def run_hgvs_checker(variant, is_a_gene=False):
    """Run the LOVD HGVS Syntax Checker with quick retries on failure, using thread lock."""
    if not os.path.exists(PHP_SCRIPT):
        raise FileNotFoundError(f"HGVS Syntax Checker is not installed. Expected at: {PHP_SCRIPT}")

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with lovd_cli_lock:
                if is_a_gene is True:
                    if PHP_PATH is not None:
                        command = [PHP_PATH, "-f", PHP_SCRIPT, f"gene:{variant}"]
                    else:
                        command = ["php", "-f", PHP_SCRIPT, f"gene:{variant}"]
                else:
                    if PHP_PATH is not None:
                        command = [PHP_PATH, "-f", PHP_SCRIPT, variant]
                    else:
                        command = ["php", "-f", PHP_SCRIPT, variant]

                result = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    check=True
                )

                if PHP_PATH is not None:
                    result_meta = subprocess.run(
                        [PHP_PATH, "-f", PHP_SCRIPT, "getVersions"],
                        capture_output=True,
                        text=True,
                        check=True
                    )
                else:
                    result_meta = subprocess.run(
                        ["php", "-f", PHP_SCRIPT, "getVersions"],
                        capture_output=True,
                        text=True,
                        check=True
                    )

            json_result = json.loads(result.stdout.strip())
            json_meta = json.loads(result_meta.stdout.strip())
            json_result[0]["metadata"] = json_meta[0]
            return json_result

        except subprocess.CalledProcessError as e:
            if attempt == MAX_RETRIES:
                raise RuntimeError(f"HGVS Checker error after {MAX_RETRIES} retries: {e.stderr.strip()}")
            time.sleep(RETRY_DELAY)

        except json.JSONDecodeError:
            if attempt == MAX_RETRIES:
                raise RuntimeError(f"Invalid JSON from HGVS Checker after {MAX_RETRIES} retries: {result.stdout.strip()}")
            time.sleep(RETRY_DELAY)

        except Exception as e:
            if attempt == MAX_RETRIES:
                raise RuntimeError(f"Unexpected error in HGVS Checker after {MAX_RETRIES} retries: {e}")
            time.sleep(RETRY_DELAY)

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
