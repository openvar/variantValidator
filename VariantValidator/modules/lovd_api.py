import logging

import requests

from VariantValidator.bin import lovd_syntax_checker


logger = logging.getLogger(__name__)


class LovdApiFlowException(Exception):
    pass


def run_lovd_checker_cli(variant, is_a_gene=False):
    """Runs the LOVD syntax checker via CLI."""
    if not is_a_gene:
        base_url = "https://api.lovd.nl/v2/checkHGVS"
    else:
        base_url = "https://api.lovd.nl/v2/checkGene"

    url = f"{base_url}/{variant}"

    logger.info(
        "Calling LOVD CLI with: %s",
        variant,
    )

    try:
        result = lovd_syntax_checker.run_hgvs_checker(
            variant,
            is_a_gene,
        )[0]
        result = {"data": [result]}
        result["url"] = url
        result["version"] = result["data"][0]["metadata"]["library_version"]

        logger.info(
            "Called LOVD CLI successfully with: %s",
            variant,
        )

        return result

    except Exception as exc:
        logger.error(
            "Error running LOVD checker CLI: %s",
            exc,
        )
        return {
            "lovd_api_error": f"CLI check failed: {exc}"
        }


def run_lovd_checker_web(variant_description, is_a_gene=False):
    """Runs the LOVD syntax checker via the web API."""
    if not is_a_gene:
        base_url = "https://api.lovd.nl/v2/checkHGVS"
    else:
        base_url = "https://api.lovd.nl/v2/checkGene"

    url = f"{base_url}/{variant_description}"

    logger.info(
        "Calling LOVD API with: %s",
        variant_description,
    )

    try:
        if is_a_gene:
            raise LovdApiFlowException(
                "Web API is currently not configured to support gene symbols"
            )

        response = requests.get(url)
        response.raise_for_status()

        json_data = response.json()
        json_data["url"] = url
        json_data["version"] = json_data["versions"]["library_version"]

        logger.info(
            "Called LOVD Web API successfully with: %s",
            variant_description,
        )

        return remove_double_quotes(json_data)

    except requests.RequestException as exc:
        return {
            "lovd_api_error": f"Request failed: {exc}"
        }

    except LovdApiFlowException as exc:
        return {
            "lovd_api_error": f"Unsupported value: {exc}"
        }

    except Exception as exc:
        return {
            "lovd_api_error": f"Unexpected error: {exc}"
        }


def lovd_syntax_check(
        variant_description,
        do_lovd_check=True,
        is_a_gene=False,
):
    """
    Perform LOVD syntax check using CLI first, then fall back to the web API.
    """
    if not do_lovd_check:
        return {
            "lovd_api_error":
                f"Do LOVD syntax check set to {do_lovd_check}"
        }

    if ":p." in variant_description:
        raise LovdApiFlowException(
            "Protein-level variant descriptions are not supported: "
            f"{variant_description}"
        )

    if ":r." in variant_description:
        raise LovdApiFlowException(
            "RNA-level variant descriptions are not supported: "
            f"{variant_description}"
        )

    try:
        json_data = run_lovd_checker_cli(
            variant_description,
            is_a_gene=is_a_gene,
        )

        if "lovd_api_error" in json_data:
            raise ValueError(
                json_data["lovd_api_error"]
            )

    except Exception as exc:
        logger.error(
            "Error running LOVD checker CLI: %s",
            exc,
        )
        json_data = run_lovd_checker_web(
            variant_description,
            is_a_gene=is_a_gene,
        )

    json_data = remove_double_quotes(json_data)

    if not isinstance(json_data, dict):
        return {
            "lovd_api_error": "Unexpected output format"
        }

    return json_data


def remove_double_quotes(obj):
    """Recursively remove double quotes from strings in a structure."""
    if isinstance(obj, str):
        return obj.replace('"', '')

    if isinstance(obj, dict):
        return {
            key: remove_double_quotes(value)
            for key, value in obj.items()
        }

    if isinstance(obj, list):
        return [
            remove_double_quotes(item)
            for item in obj
        ]

    if isinstance(obj, tuple):
        return tuple(
            remove_double_quotes(item)
            for item in obj
        )

    if isinstance(obj, set):
        return {
            remove_double_quotes(item)
            for item in obj
        }

    return obj


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
