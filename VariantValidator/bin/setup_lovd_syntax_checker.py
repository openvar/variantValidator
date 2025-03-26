import os
import urllib.request
import importlib.resources as pkg_resources

# Define installation path dynamically using importlib
def get_installation_path():
    """Determine the correct installation path for the script inside VariantValidator."""
    try:
        # Get the installed package directory
        package_path = str(pkg_resources.files("VariantValidator"))  # Works with installed packages
        return os.path.join(package_path, "php", "lovd_hgvs_checker")
    except ModuleNotFoundError:
        print("Error: VariantValidator package is not installed.")
        exit(1)

# Define installation directories
LOVD_DIR = get_installation_path()
WEB_DIR = os.path.join(LOVD_DIR, "web")

# Define file paths
PHP_SCRIPT = os.path.join(LOVD_DIR, "HGVS.php")
AJAX_SCRIPT = os.path.join(WEB_DIR, "ajax.php")
INDEX_SCRIPT = os.path.join(WEB_DIR, "index.php")

# Define URLs
HGVS_CHECKER_URL = "https://raw.githubusercontent.com/LOVDnl/HGVS-syntax-checker/main/HGVS.php"
AJAX_URL = "https://raw.githubusercontent.com/LOVDnl/HGVS-syntax-checker/main/web/ajax.php"
INDEX_URL = "https://raw.githubusercontent.com/LOVDnl/HGVS-syntax-checker/main/web/index.php"

def download_file(url, dest):
    """Download a file from a URL and save it to the specified destination."""
    try:
        print(f"Downloading {url} -> {dest}...")
        urllib.request.urlretrieve(url, dest)
        os.chmod(dest, 0o755)  # Make the file executable if needed
        print(f"Downloaded: {dest}")
    except Exception as e:
        print(f"Failed to download {url}: {e}")

def setup_lovd():
    """Setup LOVD HGVS Syntax Checker by downloading required files."""
    os.makedirs(LOVD_DIR, exist_ok=True)
    os.makedirs(WEB_DIR, exist_ok=True)

    download_file(HGVS_CHECKER_URL, PHP_SCRIPT)
    download_file(AJAX_URL, AJAX_SCRIPT)
    download_file(INDEX_URL, INDEX_SCRIPT)

    print("Setup complete! HGVS Syntax Checker is ready to use.")

if __name__ == "__main__":
    setup_lovd()
