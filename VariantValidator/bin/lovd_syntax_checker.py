import subprocess
import json
import time
import threading
import importlib.resources as pkg_resources
import os
from configparser import ConfigParser

CONFIG_DIR = os.path.join(os.path.expanduser("~"), ".variantvalidator")

# Load config
config = ConfigParser()
config.read(CONFIG_DIR)

# PHP binary path
PHP_PATH = "php"
try:
    if config["php"].get("php_path"):
        PHP_PATH = config["php"]["php_path"]
except Exception:
    pass

# PHP memory limit (env > config > default)
PHP_MEMORY_LIMIT = os.environ.get(
    "VV_PHP_MEMORY_LIMIT",
    config.get("php", "memory_limit", fallback="1G"),
)

# Base PHP command (no paths hard-coded)
PHP_BASE_CMD = [PHP_PATH, "-d", f"memory_limit={PHP_MEMORY_LIMIT}"]

# Lock to ensure thread-safe subprocess access
lovd_cli_lock = threading.Lock()


def get_installation_path():
    """Determine the correct installation path inside the installed VariantValidator package."""
    try:
        package_path = str(pkg_resources.files("VariantValidator"))
        return os.path.join(package_path, "php", "lovd_hgvs_checker")
    except ModuleNotFoundError:
        raise FileNotFoundError("Error: VariantValidator package is not installed.")


LOVD_DIR = get_installation_path()
PHP_SCRIPT = os.path.join(LOVD_DIR, "HGVS.php")

MAX_RETRIES = 3
RETRY_DELAY = 0.2  # 200 ms


def run_hgvs_checker(variant, is_a_gene=False):
    """
    Run the LOVD HGVS Syntax Checker with retries and increased PHP memory.
    Thread-safe and path-agnostic.
    """
    if not os.path.exists(PHP_SCRIPT):
        raise FileNotFoundError(
            f"HGVS Syntax Checker is not installed. Expected at: {PHP_SCRIPT}"
        )

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with lovd_cli_lock:
                if is_a_gene:
                    command = (
                        PHP_BASE_CMD
                        + ["-f", PHP_SCRIPT, f"gene:{variant}"]
                    )
                else:
                    command = (
                        PHP_BASE_CMD
                        + ["-f", PHP_SCRIPT, variant]
                    )

                result = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    check=True,
                )

                result_meta = subprocess.run(
                    PHP_BASE_CMD + ["-f", PHP_SCRIPT, "getVersions"],
                    capture_output=True,
                    text=True,
                    check=True,
                )

            json_result = json.loads(result.stdout.strip())
            json_meta = json.loads(result_meta.stdout.strip())
            json_result[0]["metadata"] = json_meta[0]
            return json_result

        except subprocess.CalledProcessError as e:
            if attempt == MAX_RETRIES:
                raise RuntimeError(
                    f"HGVS Checker error after {MAX_RETRIES} retries:\n{e.stderr.strip()}"
                )
            time.sleep(RETRY_DELAY)

        except json.JSONDecodeError:
            if attempt == MAX_RETRIES:
                raise RuntimeError(
                    f"Invalid JSON from HGVS Checker after {MAX_RETRIES} retries:\n"
                    f"{result.stdout.strip()}"
                )
            time.sleep(RETRY_DELAY)

        except Exception as e:
            if attempt == MAX_RETRIES:
                raise RuntimeError(
                    f"Unexpected error in HGVS Checker after {MAX_RETRIES} retries: {e}"
                )
            time.sleep(RETRY_DELAY)


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
