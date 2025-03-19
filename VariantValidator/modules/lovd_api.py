import requests
import subprocess
import json

def run_lovd_checker_cli(variant):
    try:
        result = subprocess.run(
            ["lovd_syntax_checker", variant],  # Calls the script
            capture_output=True,
            text=True,
            check=True
        )
        return json.loads(result.stdout.strip())  # Parse the JSON output
    except subprocess.CalledProcessError:
        raise
    except json.JSONDecodeError:
        raise

def run_lovd_checker_web(variant_description):
    base_url = "https://api.lovd.nl/v2/checkHGVS"
    url = f"{base_url}/{variant_description}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        json_data = response.json()
        json_data["url"] = url
        json_data = remove_double_quotes(json_data)
        return json_data
    except requests.RequestException:
        raise


def lovd_syntax_check(variant_description, do_lovd_check=True):
    if do_lovd_check is False:
        return {"lovd_api_error": f"Do LOVD syntax check set to {do_lovd_check}"}

    # Try the cli first
    try:
        json_data = run_lovd_checker_cli(variant_description)
        return json_data
    except subprocess.CalledProcessError:





def remove_double_quotes(obj):
    if isinstance(obj, str):
        return obj.replace('"', '')  # Remove double quotes from strings
    elif isinstance(obj, dict):
        return {k: remove_double_quotes(v) for k, v in obj.items()}  # Process dict recursively
    elif isinstance(obj, list):
        return [remove_double_quotes(i) for i in obj]  # Process lists recursively
    elif isinstance(obj, tuple):
        return tuple(remove_double_quotes(i) for i in obj)  # Process tuples recursively
    elif isinstance(obj, set):
        return {remove_double_quotes(i) for i in obj}  # Process sets recursively
    return obj  # Return unchanged if not a str, dict, list, tuple, or set

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