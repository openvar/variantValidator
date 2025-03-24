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
        return json.loads(result.stdout.strip())  # Parse JSON output
    except subprocess.CalledProcessError as e:
        return f"Error running HGVS Checker: {e.stderr.strip()}"
    except json.JSONDecodeError:
        return f"Invalid JSON output: {result.stdout.strip()}"

if __name__ == "__main__":
    test_variant = "NC_000023.11:g.32893361G>A"
    output = run_hgvs_checker(test_variant)
    print(output)
