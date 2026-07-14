#!/usr/bin/env python3

"""
VariantValidator LOVD HGVS Syntax Checker Command Line Interface.

The LOVD HGVS Syntax Checker CLI provides convenient access to the
LOVD HGVS Syntax Checker through the VariantValidator Python API.

The CLI performs HGVS syntax checking of sequence variant descriptions
or gene symbols and returns the results as formatted JSON.

The CLI first attempts to use the locally installed LOVD HGVS Syntax
Checker. If the local installation is unavailable, it automatically
falls back to the LOVD web API.

Features
--------
* Check HGVS sequence variant descriptions
* Check gene symbols
* Automatic fallback to the LOVD web API
* Pretty-printed JSON output
* Configurable logging

Unlike the VariantValidator validation engine, this CLI does not require
the VariantValidator databases or configuration to be installed.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys

import VariantValidator

from VariantValidator import logger as vvlogger
from VariantValidator.modules import lovd_api


LOGGER = logging.getLogger(__name__)


EXIT_SUCCESS = 0
EXIT_INPUT_ERROR = 2
EXIT_UNEXPECTED_ERROR = 5


EXAMPLES = """
Examples
--------

Check an HGVS variant

    lovd-hgvs-syntax-checker -q NM_000059.4:c.7790G>A

Check a genomic variant

    lovd-hgvs-syntax-checker -q NC_000013.11:g.32355250G>A

Check a gene symbol

    lovd-hgvs-syntax-checker -q BRCA2 --gene

Write JSON output to a file

    lovd-hgvs-syntax-checker \\
        -q NM_000059.4:c.7790G>A \\
        -o result.json
"""


# ----------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------


def configure_logging(level: str | None) -> None:
    """
    Configure logging for the LOVD HGVS Syntax Checker CLI.

    Parameters
    ----------
    level
        Console logging level supplied on the command line.
        If None, the logging configuration defined in
        VariantValidator.settings is used.
    """

    vvlogger.configure_logging(
        console_level=level,
    )

    LOGGER.debug("Logging initialised")

# ----------------------------------------------------------------------
# Argument parser
# ----------------------------------------------------------------------


class CustomFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawTextHelpFormatter,
):
    """
    Preserve help formatting while displaying default values.
    """
    pass


def build_parser() -> argparse.ArgumentParser:
    """
    Build the command-line argument parser.
    """

    parser = argparse.ArgumentParser(
        prog="lovd-hgvs-syntax-checker",
        description=(
            "Check HGVS sequence variant descriptions or gene symbols "
            "using the LOVD HGVS Syntax Checker."
        ),
        formatter_class=CustomFormatter,
        epilog=EXAMPLES,
    )

    #
    # Query
    #

    parser.add_argument(
        "-q",
        "--query",
        required=True,
        metavar="QUERY",
        help="""
HGVS sequence variant description or gene symbol.

Examples

  NM_000059.4:c.7790G>A
  NC_000013.11:g.32355250G>A
  BRCA2
""",
    )

    #
    # Options
    #

    options = parser.add_argument_group(
        "Options"
    )

    options.add_argument(
        "-g",
        "--gene",
        action="store_true",
        help=(
            "Interpret QUERY as a gene symbol rather "
            "than an HGVS variant description."
        ),
    )

    #
    # Output
    #

    output = parser.add_argument_group(
        "Output"
    )

    output.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        default=sys.stdout,
        metavar="FILE",
        help=(
            "Write formatted JSON to FILE "
            "(default: stdout)."
        ),
    )

    #
    # Logging
    #

    logging_group = parser.add_argument_group(
        "Logging"
    )

    logging_group.add_argument(
        "--log-level",
        default="WARNING",
        choices=[
            "DEBUG",
            "INFO",
            "WARNING",
            "ERROR",
            "CRITICAL",
        ],
        help="Console logging level.",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"VariantValidator {VariantValidator.__version__}",
    )

    return parser

# ----------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------


class Writer:
    """
    Write LOVD HGVS Syntax Checker results as formatted JSON.
    """

    def __init__(
        self,
        handle,
    ) -> None:

        self.handle = handle

    def write(
        self,
        result: dict,
    ) -> None:
        """
        Write formatted JSON.
        """

        json.dump(
            result,
            self.handle,
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
        )

        self.handle.write("\n")

    def close(
        self,
    ) -> None:
        """
        Flush the output stream.
        """

        self.handle.flush()


# ----------------------------------------------------------------------
# Query execution
# ----------------------------------------------------------------------


def run_query(
    query: str,
    is_a_gene: bool,
) -> dict:
    """
    Execute a LOVD HGVS Syntax Checker query.
    """

    LOGGER.info(
        "Running LOVD HGVS Syntax Checker"
    )

    LOGGER.info(
        "Query    : %s",
        query,
    )

    LOGGER.info(
        "Gene mode: %s",
        is_a_gene,
    )

    result = lovd_api.lovd_syntax_check(
        variant_description=query,
        is_a_gene=is_a_gene,
    )

    LOGGER.info(
        "Query complete"
    )

    return result

# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------


def main() -> int:
    """
    LOVD HGVS Syntax Checker command-line entry point.
    """

    parser = build_parser()

    args = parser.parse_args()

    configure_logging(
        args.log_level,
    )

    try:

        writer = Writer(
            args.output,
        )

        result = run_query(
            query=args.query,
            is_a_gene=args.gene,
        )

        writer.write(
            result,
        )

        writer.close()

        return EXIT_SUCCESS

    except KeyboardInterrupt:

        LOGGER.error(
            "Interrupted by user."
        )

        return EXIT_UNEXPECTED_ERROR

    except (
        FileNotFoundError,
        ValueError,
    ) as exc:

        LOGGER.error(
            "%s",
            exc,
        )

        print(
            f"Error: {exc}",
            file=sys.stderr,
        )

        return EXIT_INPUT_ERROR

    except Exception as exc:

        LOGGER.exception(
            "LOVD HGVS Syntax Checker failed."
        )

        print(
            f"Unexpected error: {exc}",
            file=sys.stderr,
        )

        return EXIT_UNEXPECTED_ERROR


if __name__ == "__main__":
    raise SystemExit(main())


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