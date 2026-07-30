#!/usr/bin/env python3

"""
VariantValidator HGVS2Reference Command Line Interface.

The HGVS2Reference CLI provides convenient access to the
VariantValidator `hgvs2ref()` method for retrieving the reference
sequence corresponding to an HGVS variant description.

Features
--------
* Query genomic, coding and non-coding HGVS variants
* Support transcript variants with explicit genomic context
* Pretty-printed JSON output
* Configurable logging
"""

from __future__ import annotations

import argparse
import json
import logging
import sys

import VariantValidator

from VariantValidator import logger as vvlogger


LOGGER = logging.getLogger(__name__)


EXIT_SUCCESS = 0
EXIT_INPUT_ERROR = 2
EXIT_UNEXPECTED_ERROR = 5


EXAMPLES = """
Examples
--------

Common usage
~~~~~~~~~~~~

Query a genomic variant

    hgvs2reference -q NC_000017.11:g.7676594C>T

Query a coding variant

    hgvs2reference -q NM_000546.6:c.215C>G

Query a transcript variant with genomic context

    hgvs2reference -q NC_000017.11(NM_000546.6):c.375+1G>A


Output
~~~~~~

Write JSON to a file

    hgvs2reference -q NM_000546.6:c.215C>G -o results.json
"""

# ----------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------


def configure_logging(level: str | None) -> None:
    """
    Configure logging for the HGVS2Reference CLI.

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
        prog="hgvs2reference",
        description=(
            "Retrieve the reference sequence corresponding "
            "to an HGVS variant description."
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
        metavar="HGVS",
        help="""
HGVS variant description.

Supported formats

  • Genomic (g.)
  • Coding (c.)
  • Non-coding transcript (n.)
  • Transcript with explicit genomic context
    (NC_(NM_) / NC_(NR_))

See the examples below.
""",
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
            "Write JSON output to FILE "
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

    #
    # Version
    #

    parser.add_argument(
        "--version",
        action="version",
        version=(
            f"hgvs2reference "
            f"(VariantValidator {VariantValidator.__version__})"
        ),
    )

    return parser

# ----------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------


class Writer:
    """
    Write HGVS2Reference results as JSON.
    """

    def __init__(
        self,
        handle,
    ) -> None:

        self.handle = handle

    def write(
        self,
        result,
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

    def close(self) -> None:
        """
        Flush the output stream.
        """

        self.handle.flush()


# ----------------------------------------------------------------------
# Query execution
# ----------------------------------------------------------------------


def run_query(
    validator: VariantValidator.Validator,
    query: str,
) -> dict:
    """
    Execute a HGVS2Reference query.
    """

    LOGGER.info(
        "Initialising HGVS2Reference"
    )

    LOGGER.info(
        "Processing %s",
        query,
    )

    result = validator.hgvs2ref(
        query,
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
    HGVS2Reference command-line entry point.
    """

    parser = build_parser()

    args = parser.parse_args()

    configure_logging(
        args.log_level,
    )

    try:

        validator = VariantValidator.Validator()

        writer = Writer(
            args.output,
        )

        result = run_query(
            validator=validator,
            query=args.query,
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

    except ValueError as exc:

        LOGGER.error(
            "%s",
            exc,
        )

        print(
            f"Error: {exc}",
            file=sys.stderr,
        )

        return EXIT_INPUT_ERROR

    except Exception:

        LOGGER.exception(
            "HGVS2Reference failed."
        )

        print(
            "Unexpected error. "
            "Run with --log-level DEBUG for additional information.",
            file=sys.stderr,
        )

        return EXIT_UNEXPECTED_ERROR


if __name__ == "__main__":
    raise SystemExit(main())


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
