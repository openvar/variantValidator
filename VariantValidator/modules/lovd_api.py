import requests
import logging
from VariantValidator.bin import lovd_syntax_checker

logger = logging.getLogger(__name__)


def run_lovd_checker_cli(variant):
    """Runs the LOVD syntax checker via CLI."""
    base_url = "https://api.lovd.nl/v2/checkHGVS"
    url = f"{base_url}/{variant}"

    try:
        result = lovd_syntax_checker.run_hgvs_checker(variant)[0]
        result = {"data": [result]}
        result["url"] = url
        result["version"] = result["data"][0]["metadata"]["library_version"]
        return result  # Ensure it returns a dictionary
    except Exception as e:
        logger.error(f"Error running LOVD checker CLI: {e}")
        return {"lovd_api_error": f"CLI check failed: {e}"}


def run_lovd_checker_web(variant_description):
    """Runs the LOVD syntax checker via the web API."""
    base_url = "https://api.lovd.nl/v2/checkHGVS"
    url = f"{base_url}/{variant_description}"

    try:
        response = requests.get(url)
        response.raise_for_status()
        json_data = response.json()
        json_data["url"] = url
        json_data["version"] = json_data["versions"]["library_version"]
        return remove_double_quotes(json_data)
    except requests.RequestException as e:
        return {"lovd_api_error": f"Request failed: {e}"}
    except Exception as e:
        return {"lovd_api_error": f"Unexpected error: {e}"}


def lovd_syntax_check(variant_description, do_lovd_check=True):
    """Performs LOVD syntax check using CLI first, then falls back to web API if necessary."""
    if not do_lovd_check:
        return {"lovd_api_error": f"Do LOVD syntax check set to {do_lovd_check}"}

    json_data = None  # Ensure json_data is always defined

    try:
        json_data = run_lovd_checker_cli(variant_description)
        if "lovd_api_error" in json_data:  # Fallback if CLI fails
            raise ValueError(json_data["lovd_api_error"])
    except Exception as e:
        logger.error(f"Error running LOVD checker CLI: {e}")
        json_data = run_lovd_checker_web(variant_description)

    json_data = remove_double_quotes(json_data)

    # Ensure the final return value is always a dictionary
    if not isinstance(json_data, dict):
        return {"lovd_api_error": "Unexpected output format"}

    return json_data


def remove_double_quotes(obj):
    """Recursively removes double quotes from all strings in a structure."""
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
