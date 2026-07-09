#!/usr/bin/env python3

"""
VariantValidator Gene2Transcripts Command Line Interface.

The Gene2Transcripts CLI provides convenient access to the
VariantValidator gene2transcripts API for retrieving transcript
information associated with genes and transcript accessions.

Features
--------
* Query gene symbols, transcript accessions or HGNC identifiers
* Single or multiple queries
* JSON or pipe-delimited input
* Read queries from file
* Select transcript subsets
* RefSeq or Ensembl transcript models
* Optional LOVD syntax checking
* Optional HGNC web lookup bypass
* Optional genomic span suppression
* Pretty-printed JSON output
* Configurable logging
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import VariantValidator

from VariantValidator import logger as vvlogger


LOGGER = logging.getLogger(__name__)


EXIT_SUCCESS = 0
EXIT_INPUT_ERROR = 2
EXIT_UNEXPECTED_ERROR = 5


TRANSCRIPT_KEYWORDS = {
    "mane",
    "mane_select",
    "all",
    "raw",
    "select",
}

EXAMPLES = """
Examples
--------

Common usage
~~~~~~~~~~~~

Query a gene symbol

    gene2transcripts -q COL1A1

Query a transcript accession

    gene2transcripts -q NM_000088.4

Query an HGNC identifier

    gene2transcripts -q HGNC:2197


Transcript selection
~~~~~~~~~~~~~~~~~~~~

Use one of the built-in transcript selection options

    gene2transcripts -q BRAF --select-transcripts mane_select

    gene2transcripts -q BRCA2 --select-transcripts all

Return a specific transcript

    gene2transcripts -q BRCA2 --select-transcripts NM_000059.4

Return multiple transcripts

    gene2transcripts -q BRCA2 --select-transcripts "NM_000059.3|NM_000059.4"


Transcript models
~~~~~~~~~~~~~~~~~

Use RefSeq transcripts

    gene2transcripts -q COL1A1 --transcript-model refseq

Use Ensembl transcripts

    gene2transcripts -q COL1A1 --transcript-model ensembl


Options
~~~~~~~

Disable HGNC web searches

    gene2transcripts -q COL1A1 --no-web-searches

Omit genomic span information

    gene2transcripts -q COL1A1 --no-genomic-spans

Enable LOVD syntax checking

    gene2transcripts -q C3orf52 --lovd-syntax-check


Input formats
~~~~~~~~~~~~~

Multiple command-line arguments

    gene2transcripts -q COL1A1 BRCA2 P3H1

JSON array (recommended)

    gene2transcripts -q '[
        "COL1A1",
        "BRCA2",
        "P3H1"
    ]'

Pipe-delimited list

    gene2transcripts -q "COL1A1|BRCA2|P3H1"

Text file

    gene2transcripts -q @/path/to/queries.txt

JSON file

    gene2transcripts -q @/path/to/queries.json


Output
~~~~~~

Write JSON to a file

    gene2transcripts -q COL1A1 -o results.json
"""

# ----------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------


def configure_logging(level: str | None) -> None:
    """
    Configure logging for the Gene2Transcripts CLI.

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


def parse_queries(values: list[str]) -> list[str]:
    """
    Parse query input supplied on the command line.

    Supports

    * repeated arguments
    * JSON array
    * pipe-delimited string
    * @filename
    """

    queries: list[str] = []

    for value in values:

        if value.startswith("@"):

            path = Path(value[1:])

            LOGGER.info(
                "Loading queries from %s",
                path,
            )

            if not path.exists():
                raise FileNotFoundError(path)

            if path.suffix.lower() in {
                ".json",
                ".jsn",
            }:
                queries.extend(
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

                    queries.append(line)

            continue

        queries.extend(
            _parse_value(value)
        )

    if not queries:
        raise ValueError(
            "No queries were supplied."
        )

    LOGGER.info(
        "Parsed %d queries",
        len(queries),
    )

    return queries


def parse_transcripts(
    value: str,
) -> str | list[str]:
    """
    Parse transcript selection.
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
        prog="gene2transcripts",
        description=(
            "Retrieve transcript information for genes "
            "and transcript accessions."
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
        nargs="+",
        required=True,
        metavar="QUERY",
        help="""
Gene or transcript query.

Supported formats

  • Gene symbol
  • Transcript accession
  • HGNC identifier
  • Multiple arguments
  • JSON array (recommended)
  • Pipe-delimited string
  • @text file
  • @JSON file

See the examples below.
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
            "Reference genome assembly."
        ),
    )

    #
    # Transcript selection
    #

    parser.add_argument(
        "-t",
        "--select-transcripts",
        default=None,
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
        ],
        default="refseq",
        metavar="MODEL",
        help="""
Transcript model.

Choices

    refseq
    ensembl
""",
    )

    #
    # Options
    #

    options = parser.add_argument_group(
        "Options"
    )

    options.add_argument(
        "--no-web-searches",
        action="store_true",
        help=(
            "Disable HGNC web lookups for "
            "gene symbol correction."
        ),
    )

    options.add_argument(
        "--no-genomic-spans",
        action="store_true",
        help=(
            "Omit exon boundary and genomic "
            "alignment information."
        ),
    )

    options.add_argument(
        "--lovd-syntax-check",
        action="store_true",
        help=(
            "Enable LOVD syntax checking."
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
    Write Gene2Transcripts results as JSON.
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
    queries: list[str],
    args: argparse.Namespace,
) -> dict:
    """
    Execute one or more Gene2Transcripts queries.
    """

    LOGGER.info(
        "Initialising Gene2Transcripts"
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
        "Number of queries    : %d",
        len(queries),
    )

    batch = json.dumps(queries)

    result = validator.gene2transcripts(
        batch,
        select_transcripts=args.select_transcripts,
        transcript_set=args.transcript_model,
        genome_build=args.genome,
        bypass_web_searches=args.no_web_searches,
        bypass_genomic_spans=args.no_genomic_spans,
        lovd_syntax_check=args.lovd_syntax_check,
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
    Gene2Transcripts command-line entry point.
    """

    parser = build_parser()

    args = parser.parse_args()

    configure_logging(
        args.log_level,
    )

    try:

        queries = parse_queries(
            args.query,
        )

        if args.select_transcripts is not None:
            args.select_transcripts = (
                parse_transcripts(
                    args.select_transcripts,
                )
            )

        validator = VariantValidator.Validator()

        writer = Writer(
            args.output,
        )

        result = run_query(
            validator=validator,
            queries=queries,
            args=args,
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
            "Gene2Transcripts failed."
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
