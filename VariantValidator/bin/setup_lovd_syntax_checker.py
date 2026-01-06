import os
import urllib.request
import subprocess
import importlib.resources as pkg_resources

# Define installation path dynamically using importlib
def get_installation_path():
    try:
        package_path = str(pkg_resources.files("VariantValidator"))
        return os.path.join(package_path, "php", "lovd_hgvs_checker")
    except ModuleNotFoundError:
        print("Error: VariantValidator package is not installed.")
        exit(1)

# Define directories
LOVD_DIR = get_installation_path()
WEB_DIR = os.path.join(LOVD_DIR, "web")
CACHE_DIR = os.path.join(LOVD_DIR, "cache")  # <-- NEW

# Define file paths
PHP_SCRIPT = os.path.join(LOVD_DIR, "HGVS.php")
AJAX_SCRIPT = os.path.join(WEB_DIR, "ajax.php")
INDEX_SCRIPT = os.path.join(WEB_DIR, "index.php")
UPDATE_SCRIPT = os.path.join(CACHE_DIR, "update.php")  # <-- NEW

# Define URLs
HGVS_CHECKER_URL = "https://raw.githubusercontent.com/LOVDnl/HGVS-syntax-checker/main/HGVS.php"
AJAX_URL = "https://raw.githubusercontent.com/LOVDnl/HGVS-syntax-checker/main/web/ajax.php"
INDEX_URL = "https://raw.githubusercontent.com/LOVDnl/HGVS-syntax-checker/main/web/index.php"
UPDATE_URL = "https://raw.githubusercontent.com/LOVDnl/HGVS-syntax-checker/main/cache/update.php"  # <-- NEW

def download_file(url, dest):
    try:
        print(f"Downloading {url} -> {dest}...")
        urllib.request.urlretrieve(url, dest)
        os.chmod(dest, 0o755)
        print(f"Downloaded: {dest}")
    except Exception as e:
        print(f"Failed to download {url}: {e}")

def setup_lovd():
    os.makedirs(LOVD_DIR, exist_ok=True)
    os.makedirs(WEB_DIR, exist_ok=True)
    os.makedirs(CACHE_DIR, exist_ok=True)

    download_file(HGVS_CHECKER_URL, PHP_SCRIPT)
    download_file(AJAX_URL, AJAX_SCRIPT)
    download_file(INDEX_URL, INDEX_SCRIPT)
    download_file(UPDATE_URL, UPDATE_SCRIPT)

    try:
        os.chmod(CACHE_DIR, 0o777)  # Make directory writable
    except Exception as e:
        print(f"Warning: could not change permissions on {CACHE_DIR}: {e}")

    print("Running PHP cache updater...")
    try:
        subprocess.run(["php", "-f", UPDATE_SCRIPT], check=True)
        print("Cache updated successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Failed to run PHP cache update: {e}")

    print("Setup complete! HGVS Syntax Checker is ready to use.")

if __name__ == "__main__":
    setup_lovd()

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
