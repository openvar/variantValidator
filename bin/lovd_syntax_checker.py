import os
import subprocess
import json

LOVD_DIR = os.path.expanduser("~/.lovd_hgvs_checker")
PHP_SCRIPT = os.path.join(LOVD_DIR, "HGVS.php")

def run_hgvs_checker(variant):
    """Run the LOVD HGVS Syntax Checker with the given variant using the -f option."""
    if not os.path.exists(PHP_SCRIPT):
        raise FileNotFoundError("HGVS Syntax Checker is not installed. "
                                "Please ensure the HGVS.php script is located in the ~/.lovd_hgvs_checker/ directory.")

    try:
        result = subprocess.run(
            ["php", "-f", PHP_SCRIPT, variant],
            capture_output=True,
            text=True,
            check=True
        )
        # Attempt to parse the output as JSON
        return json.loads(result.stdout.strip())
    except subprocess.CalledProcessError as e:
        return f"Error running HGVS Checker: {e.stderr.strip()}"
    except json.JSONDecodeError:
        return f"Invalid JSON output: {result.stdout.strip()}"

if __name__ == "__main__":
    test_variant = "NC_000023.11:g.32893361G>A"
    output = run_hgvs_checker(test_variant)
    print(output)
