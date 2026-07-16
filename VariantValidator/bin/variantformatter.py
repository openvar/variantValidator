#!/usr/bin/env python3

"""
VariantFormatter Command Line Interface.

The VariantFormatter CLI provides a convenient interface to the
VariantFormatter Python API for converting HGVS sequence variants
between genomic, transcript and protein representations.

Features
--------
* Single or multiple variant formatting
* JSON or pipe-delimited input
* Read variants from file
* Configurable transcript selection
* Optional genomic syntax checking only
* Optional genomic liftover
* JSON output
* Logging suitable for debugging
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
import time

from VariantFormatter.simpleVariantFormatter import (
    SimpleVariantFormatter,
)

from VariantValidator import logger as vvlogger
import VariantFormatter


LOGGER = logging.getLogger(__name__)


EXIT_SUCCESS = 0
EXIT_INPUT_ERROR = 2
EXIT_UNEXPECTED_ERROR = 5

EXAMPLES = """
Examples
--------

Common usage
~~~~~~~~~~~~

Format a genomic HGVS variant

    variantformatter -v "NC_000017.10:g.48275363C>A"

Format a pseudo-VCF variant

    variantformatter -v "17-50198002-C-A"

    variantformatter -v "17:50198002:C:A"

Format using Ensembl transcripts

    variantformatter -v "NC_000017.10:g.48275363C>A" --transcript-model ensembl

Format using all transcript models

    variantformatter -v "NC_000017.10:g.48275363C>A" --transcript-model all


Transcript selection
~~~~~~~~~~~~~~~~~~~~

Use one of the built-in transcript selection options

    variantformatter -v "NC_000017.10:g.48275363C>A" --select-transcripts mane_select

    variantformatter -v "NC_000017.10:g.48275363C>A" --select-transcripts all

Return specific transcripts

    variantformatter -v "NC_000017.10:g.48275363C>A" --select-transcripts NM_000093.4

    variantformatter -v "NC_000017.10:g.48275363C>A" --select-transcripts "NM_000093.4|NM_001278074.1|NM_000093.3"


Formatting options
~~~~~~~~~~~~~~~~~~

Check genomic syntax only

    variantformatter -v "NC_000017.10:g.48275363C>A" --check-only

Generate lifted-over genomic representations

    variantformatter -v "NC_000017.10:g.48275363C>A" --liftover


Input formats
~~~~~~~~~~~~~

Multiple command-line arguments

    variantformatter -v "NC_000017.10:g.48275363C>A" "17-50198002-C-A"

JSON array (recommended)

    variantformatter -v '[
        "NC_000017.10:g.48275363C>A",
        "17-50198002-C-A"
    ]'

Pipe-delimited list

    variantformatter -v "NC_000017.10:g.48275363C>A|17-50198002-C-A"

Text file

    variantformatter -v @/path/to/variants.txt

JSON file

    variantformatter -v @/path/to/variants.json


Output
~~~~~~

Write formatted JSON to a file

    variantformatter -v "NC_000017.10:g.48275363C>A" -o results.json
"""


TRANSCRIPT_KEYWORDS = {
    "mane",
    "mane_select",
    "all",
    "raw",
    "select",
}

# ----------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------


def configure_logging(level: str | None) -> None:
    """
    Configure logging for the VariantFormatter CLI.

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
# Input parsing
# ----------------------------------------------------------------------


def _load_json_file(path: Path) -> list[str]:
    """
    Load a JSON array from disk.
    """

    with path.open() as handle:
        data = json.load(handle)

    if not isinstance(data, list):
        raise ValueError(
            f"{path} does not contain a JSON array."
        )

    return [str(item) for item in data]


def _parse_value(value: str) -> list[str]:
    """
    Parse a string into a list.

    Accepted formats

    * single value
    * pipe-delimited string
    * JSON array
    """

    value = value.strip()

    if not value:
        return []

    if value.startswith("["):
        try:
            parsed = json.loads(value)
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"Invalid JSON input:\n{exc}"
            ) from exc

        if not isinstance(parsed, list):
            raise ValueError(
                "JSON input must be an array."
            )

        return [str(item) for item in parsed]

    if "|" in value:
        return [
            item.strip()
            for item in value.split("|")
            if item.strip()
        ]

    return [value]


def parse_variants(values: list[str]) -> list[str]:
    """
    Parse variant input supplied on the command line.

    Supports

    * repeated arguments
    * JSON array
    * pipe-delimited string
    * @filename
    """

    variants: list[str] = []

    for value in values:

        if value.startswith("@"):

            path = Path(value[1:])

            LOGGER.info(
                "Loading variants from %s",
                path,
            )

            if not path.exists():
                raise FileNotFoundError(path)

            if path.suffix.lower() in {".json", ".jsn"}:
                variants.extend(
                    _load_json_file(path)
                )
                continue

            with path.open() as handle:

                for line in handle:

                    line = line.strip()

                    if not line:
                        continue

                    if line.startswith("#"):
                        continue

                    variants.append(line)

            continue

        variants.extend(
            _parse_value(value)
        )

    if not variants:
        raise ValueError(
            "No variants were supplied."
        )

    LOGGER.info(
        "Parsed %d variants",
        len(variants),
    )

    return variants


def parse_transcripts(value: str) -> str | list[str]:
    """
    Parse transcript selection.

    Returns either

    * keyword
    * list of transcript accessions
    """

    value = value.strip()

    if value in TRANSCRIPT_KEYWORDS:
        return value

    parsed = _parse_value(value)

    if not parsed:
        raise ValueError(
            "Transcript selection is empty."
        )

    return parsed

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
        prog="variantformatter",
        description=(
            "Format HGVS sequence variants using VariantFormatter."
        ),
        epilog=EXAMPLES,
        formatter_class=CustomFormatter,
    )

    #
    # Required arguments
    #

    parser.add_argument(
        "-v",
        "--variant",
        nargs="+",
        required=True,
        metavar="VARIANT",
        help="""
Variant(s) to format.

Supported inputs include:

  • Genomic HGVS
  • Pseudo-VCF
  • VCF records
  • Multiple variants
  • JSON arrays
  • @text file
  • @JSON file

See the examples below for details.
""",
    )

    #
    # Reference genome
    #

    parser.add_argument(
        "-g",
        "--genome",
        "--assembly",
        default="GRCh38",
        choices=[
            "GRCh37",
            "GRCh38",
            "hg19",
            "hg38",
        ],
        help=(
            "Reference genome assembly "
            "(default: %(default)s)."
        ),
    )

    #
    # Transcript selection
    #

    parser.add_argument(
        "-t",
        "--select-transcripts",
        default="mane_select",
        metavar="TRANSCRIPTS",
        help="""
Transcript selection.

Keywords

  mane
  mane_select
  all
  raw
  select

Or supply a JSON array or pipe-delimited list.

See the examples below.
""",
    )

    #
    # Transcript model
    #

    parser.add_argument(
        "--transcript-model",
        choices=[
            "refseq",
            "ensembl",
            "all",
        ],
        default=None,
        metavar="MODEL",
        help="""
Transcript model to use.

Choices

    refseq
    ensembl
    all
""",
    )

    #
    # Formatting options
    #

    formatting = parser.add_argument_group(
        "Formatting options"
    )

    formatting.add_argument(
        "--check-only",
        action="store_true",
        help=(
            "Validate genomic HGVS syntax only. "
            "Skip transcript and protein mapping."
        ),
    )

    formatting.add_argument(
        "--liftover",
        action="store_true",
        help=(
            "Generate equivalent representations "
            "on compatible genome assemblies."
        ),
    )

    #
    # Output
    #

    output = parser.add_argument_group(
        "Output options"
    )

    output.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        default=sys.stdout,
        metavar="FILE",
        help=(
            "Output file "
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
        help=(
            "Logging verbosity "
            "(default: %(default)s)."
        ),
    )

    #
    # Miscellaneous
    #

    parser.add_argument(
        "--version",
        action="version",
        version=f"VariantFormatter {VariantFormatter.__version__}",
    )

    return parser

# ----------------------------------------------------------------------
# Output formatting
# ----------------------------------------------------------------------


class Writer:
    """
    Write VariantFormatter results as JSON.
    """

    def __init__(
            self,
            handle,
    ) -> None:

        self.handle = handle

    def write(self, result: dict) -> None:
        """
        Write a formatted JSON document.
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
# Formatting
# ----------------------------------------------------------------------


def run_formatting(
    formatter: SimpleVariantFormatter,
    variants: list[str],
    args: argparse.Namespace,
) -> dict:
    """
    Execute VariantFormatter.

    Returns
    -------
    dict
        Formatted variant data.
    """

    LOGGER.info(
        "Initialising formatting"
    )

    LOGGER.info(
        "Assembly             : %s",
        args.genome,
    )

    LOGGER.info(
        "Transcript model     : %s",
        args.transcript_model,
    )

    LOGGER.info(
        "Transcript selection : %s",
        args.select_transcripts,
    )

    LOGGER.info(
        "Check only           : %s",
        args.check_only,
    )

    LOGGER.info(
        "Liftover             : %s",
        args.liftover,
    )

    LOGGER.info(
        "Number of variants   : %d",
        len(variants),
    )

    #
    # Submit all variants in a single request.
    #

    batch_input = json.dumps(variants)

    LOGGER.debug(
        "Submitting %d variant(s)",
        len(variants),
    )

    result = formatter.format(
        batch_input=batch_input,
        genome_build=args.genome,
        transcript_model=args.transcript_model,
        specify_transcripts=(
            json.dumps(args.select_transcripts)
            if isinstance(args.select_transcripts, list)
            else args.select_transcripts
        ),
        checkOnly=args.check_only,
        liftover=args.liftover,
    )

    LOGGER.debug(
        "Formatting returned %d variant(s)",
        len(result),
    )

    LOGGER.info(
        "Formatting complete"
    )

    return result

# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------


def main() -> int:
    """
    VariantFormatter command-line entry point.
    """

    parser = build_parser()

    try:

        args = parser.parse_args()

        configure_logging(
            args.log_level,
        )

        start = time.perf_counter()

        LOGGER.info(
            "VariantFormatter %s",
            VariantFormatter.__version__,
        )

        LOGGER.debug(
            "Python %s",
            sys.version.split()[0],
        )

        formatter = SimpleVariantFormatter()

        LOGGER.debug(
            "Formatter initialised"
        )

        variants = parse_variants(
            args.variant,
        )

        args.select_transcripts = (
            parse_transcripts(
                args.select_transcripts,
            )
        )

        writer = Writer(
            handle=args.output,
        )

        result = run_formatting(
            formatter=formatter,
            variants=variants,
            args=args,
        )

        writer.write(result)

        writer.close()

        LOGGER.info(
            "Finished successfully"
        )

        LOGGER.info(
            "Completed in %.2f seconds",
            time.perf_counter() - start,
        )

        return EXIT_SUCCESS

    except KeyboardInterrupt:

        LOGGER.error(
            "Interrupted by user."
        )

        return EXIT_UNEXPECTED_ERROR

    except FileNotFoundError as exc:

        LOGGER.error(
            "Unable to locate %s",
            exc,
        )

        return EXIT_INPUT_ERROR

    except ValueError as exc:

        LOGGER.error(
            "%s",
            exc,
        )

        return EXIT_INPUT_ERROR

    except Exception:

        LOGGER.exception(
            "Unexpected error"
        )

        return EXIT_UNEXPECTED_ERROR

    finally:

        try:

            if (
                "args" in locals()
                and args.output is not sys.stdout
            ):
                args.output.close()

        except Exception:
            pass

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
