import os
import urllib.request
import subprocess
import importlib.resources as pkg_resources
import sys


# Define installation path dynamically using importlib
def get_installation_path():
    try:
        package_path = str(pkg_resources.files("VariantValidator"))
        return os.path.join(package_path, "php", "lovd_hgvs_checker")
    except ModuleNotFoundError:
        print("Error: VariantValidator package is not installed.")
        sys.exit(1)


# Define directories
LOVD_DIR = get_installation_path()
WEB_DIR = os.path.join(LOVD_DIR, "web")
CACHE_DIR = os.path.join(LOVD_DIR, "cache")

# Define file paths
PHP_SCRIPT = os.path.join(LOVD_DIR, "HGVS.php")
AJAX_SCRIPT = os.path.join(WEB_DIR, "ajax.php")
INDEX_SCRIPT = os.path.join(WEB_DIR, "index.php")
UPDATE_SCRIPT = os.path.join(CACHE_DIR, "update.php")

# Define URLs
HGVS_CHECKER_URL = (
    "https://raw.githubusercontent.com/LOVDnl/HGVS-syntax-checker/main/HGVS.php"
)
AJAX_URL = (
    "https://raw.githubusercontent.com/LOVDnl/HGVS-syntax-checker/main/web/ajax.php"
)
INDEX_URL = (
    "https://raw.githubusercontent.com/LOVDnl/HGVS-syntax-checker/main/web/index.php"
)
UPDATE_URL = (
    "https://raw.githubusercontent.com/LOVDnl/HGVS-syntax-checker/main/cache/update.php"
)


def download_file(url, dest):
    try:
        print(f"Downloading {url} -> {dest}...")
        urllib.request.urlretrieve(url, dest)
        os.chmod(dest, 0o755)
        print(f"Downloaded: {dest}")
    except Exception as e:
        print(f"Failed to download {url}: {e}")


def run_php_cache_update():
    """
    Run the LOVD HGVS cache updater.
    Retry with an increased PHP memory limit if the first attempt fails.
    """
    print("Running PHP cache updater...")

    # First attempt: default PHP settings
    try:
        subprocess.run(["php", "-f", UPDATE_SCRIPT], check=True)
        print("Cache updated successfully.")
        return
    except subprocess.CalledProcessError:
        print(
            "PHP cache update failed (likely due to PHP memory limits). "
            "Retrying with increased PHP memory..."
        )

    # Second attempt: increase PHP memory limit
    try:
        subprocess.run(
            ["php", "-d", "memory_limit=1G", "-f", UPDATE_SCRIPT],
            check=True,
        )
        print("Cache updated successfully with increased PHP memory.")
    except subprocess.CalledProcessError as e:
        print(
            "Failed to update the LOVD HGVS cache even after increasing "
            "PHP memory limit.\n"
            "Please check your PHP installation."
        )
        print(f"Details: {e}")


def setup_lovd():
    os.makedirs(LOVD_DIR, exist_ok=True)
    os.makedirs(WEB_DIR, exist_ok=True)
    os.makedirs(CACHE_DIR, exist_ok=True)

    download_file(HGVS_CHECKER_URL, PHP_SCRIPT)
    download_file(AJAX_URL, AJAX_SCRIPT)
    download_file(INDEX_URL, INDEX_SCRIPT)
    download_file(UPDATE_URL, UPDATE_SCRIPT)

    try:
        os.chmod(CACHE_DIR, 0o777)
    except Exception as e:
        print(f"Warning: could not change permissions on {CACHE_DIR}: {e}")

    run_php_cache_update()

    print("Setup complete! HGVS Syntax Checker is ready to use.")


if __name__ == "__main__":
    setup_lovd()


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
