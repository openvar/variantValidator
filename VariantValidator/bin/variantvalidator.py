#!/usr/bin/env python3

"""
VariantValidator Command Line Interface.

The VariantValidator CLI provides a convenient interface to the
VariantValidator Python API for validating HGVS sequence variants.

Features
--------
* Single or multiple variant validation
* JSON or pipe-delimited input
* Read variants from file
* Configurable transcript selection
* Automatic LOVD HGVS syntax checking
* Automatic genomic liftover
* Multiple output formats
* Logging suitable for debugging
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
import time

from VariantValidator import Validator
from VariantValidator import logger as vvlogger
from VariantValidator.version import __version__


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

Validate a genomic variant (default transcript selection: mane_select)

    variantvalidator -v "NC_000017.10:g.48275363C>A"

Validate a genomic variant against all compatible transcripts

    variantvalidator -v "NC_000017.10:g.48275363C>A" --select-transcripts all

Validate a transcript variant

    variantvalidator -v "NM_000088.4:c.589G>T" --select-transcripts all

Validate variants from a text file and write a TSV report

    variantvalidator -v @/path/to/variants.txt --select-transcripts all -o /path/to/results.tsv

Validate variants from a JSON file and write JSON results

    variantvalidator -v @/path/to/variants.json --select-transcripts all -f json -o /path/to/results.json


Transcript selection
~~~~~~~~~~~~~~~~~~~~

Use one of the built-in transcript selection options

    variantvalidator -v "NC_000017.10:g.48275363C>A" --select-transcripts mane_select

    variantvalidator -v "NC_000017.10:g.48275363C>A" --select-transcripts mane

    variantvalidator -v "NC_000017.10:g.48275363C>A" --select-transcripts all

    variantvalidator -v "NC_000017.10:g.48275363C>A" --select-transcripts raw

Select specific transcripts

    variantvalidator -v "NC_000017.10:g.48275363C>A" --select-transcripts '["NM_000088.4","NM_000088.3"]'

Use Ensembl transcripts

    variantvalidator -v "ENST00000269305.9:c.215C>G" --select-transcripts all --transcript-set ensembl


Input formats
~~~~~~~~~~~~~

Multiple command-line arguments

    variantvalidator -v "NM_000546.6:c.215C>G" "NM_004006.3:c.274G>T" --select-transcripts all

JSON array (recommended)

    variantvalidator -v '[
        "NM_000546.6:c.215C>G",
        "NM_004006.3:c.274G>T"
    ]' --select-transcripts all

Pipe-delimited list

    variantvalidator -v "NM_000546.6:c.215C>G|NM_004006.3:c.274G>T" --select-transcripts all

Text file

    variantvalidator -v @/path/to/variants.txt

    One variant per line.
    Blank lines and lines beginning with '#' are ignored.

JSON file

    variantvalidator -v @/path/to/variants.json

    JSON array of variant strings.
"""

# ----------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------


def configure_logging(level: str | None) -> None:
    """
    Configure logging for the VariantValidator CLI.

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
        prog="variantvalidator",
        description=(
            "Validate HGVS sequence variants using VariantValidator."
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
Variant(s) to validate.

Supported formats

  • Single variant
  • Multiple arguments
  • JSON array (recommended)
  • Pipe-delimited string
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
    # Transcript database
    #

    parser.add_argument(
        "--transcript-set",
        choices=[
            "refseq",
            "ensembl"
        ],
        default="refseq",
        metavar="SET",
        help="""
Restrict validation to a transcript source.

Choices

    refseq
    ensembl
""",
    )

    #
    # Validation behaviour
    #

    validation = parser.add_argument_group(
        "Validation options"
    )

    validation.add_argument(
        "--no-liftover",
        dest="liftover_level",
        action="store_false",
        default=True,
        help=(
            "Disable liftover between genome builds."
        ),
    )

    validation.add_argument(
        "--no-lovd-syntax-check",
        dest="lovd_syntax_check",
        action="store_false",
        default=True,
        help=(
            "Disable the LOVD HGVS syntax checker."
        ),
    )

    validation.add_argument(
        "--shorthand-vcf",
        action="store_true",
        help=(
            "Convert large deletions "
            "and inversions to shorthand "
            "VCF-style."
        ),
    )

    #
    # Output
    #

    output = parser.add_argument_group(
        "Output options"
    )

    output.add_argument(
        "-f",
        "--output-format",
        choices=[
            "json",
            "table",
        ],
        default=None,
        help="""
Output format.

If omitted

    stdout -> json

    output file -> table

table
    Human-readable tab-separated table.

json
    JSON document.
""",
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

    output.add_argument(
        "-m",
        "--meta",
        action="store_true",
        help="Include metadata.",
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
        version=f"VariantValidator {__version__}",
    )

    return parser

# ----------------------------------------------------------------------
# Output formatting
# ----------------------------------------------------------------------

class Writer:
    """
    Write VariantValidator validation results.

    Formats
    ------
    json
        One JSON document per validation.

    table
        Human-readable table.

    Notes
    -----
    Table output can either repeat the header for every validation
    (interactive console use) or write it only once (TSV file output).
    """

    def __init__(
            self,
            handle,
            output_format: str = "table",
            with_meta: bool = False,
            repeat_header: bool = True,
    ) -> None:

        self.handle = handle
        self.output_format = output_format.lower()
        self.with_meta = with_meta
        self.repeat_header = repeat_header

        self._header_written = False

    def write(self, validation) -> None:

        if self.output_format == "json":
            self._write_json(validation)
        elif self.output_format == "table":
            self._write_table(validation)
        else:
            raise ValueError(
                f"Unknown output format '{self.output_format}'."
            )

    def close(self) -> None:
        self.handle.flush()

    @staticmethod
    def _clean(value) -> str:

        if value is None or value == "":
            return "None"

        if value is True:
            return "True"

        if value is False:
            return "False"

        return str(value)

    def _write_json(self, validation) -> None:

        json.dump(
            json.loads(
                validation.format_as_json(
                    with_meta=self.with_meta
                )
            ),
            self.handle,
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
        )

        self.handle.write("\n")

        self.handle.write("\n")

    def _write_table(self, validation) -> None:

        table = validation.format_as_table(
            with_meta=self.with_meta
        )

        for row_number, row in enumerate(table):

            #
            # First row is the validation heading.
            # Keep it for interactive output only.
            #
            if row_number == 0:

                if not self.repeat_header:
                    continue

                self.handle.write(str(row))
                self.handle.write("\n")
                continue

            #
            # Second row contains column names.
            #
            if row_number == 1:

                if (
                    self._header_written
                    and not self.repeat_header
                ):
                    continue

                self._header_written = True

            if isinstance(row, list):

                line = "\t".join(
                    self._clean(v)
                    for v in row
                )

            else:

                line = self._clean(row)

            self.handle.write(line)
            self.handle.write("\n")

        if self.repeat_header:
            self.handle.write("\n")


# ----------------------------------------------------------------------
# Validation
# ----------------------------------------------------------------------


def run_validation(
    validator: Validator,
    variants: list[str],
    args: argparse.Namespace,
):
    """
    Execute VariantValidator.

    Yields
    ------
    ValidationOutput
        Validation results.
    """

    LOGGER.info(
        "Initialising validation"
    )

    LOGGER.info(
        "Assembly             : %s",
        args.genome,
    )

    LOGGER.info(
        "Transcript selection : %s",
        args.select_transcripts,
    )

    LOGGER.info(
        "Transcript set       : %s",
        args.transcript_set,
    )

    LOGGER.info(
        "Number of variants    : %d",
        len(variants),
    )

    #
    # Submit all variants in a single request.
    #
    request = json.dumps(variants)

    LOGGER.debug(
        "Submitting %d variant(s)",
        len(variants),
    )

    result = validator.validate(
        batch_variant=request,
        selected_assembly=args.genome,
        select_transcripts=(
            json.dumps(args.select_transcripts)
            if isinstance(args.select_transcripts, list)
            else args.select_transcripts
        ),
        transcript_set=args.transcript_set,
        liftover_level=args.liftover_level,
        lovd_syntax_check=args.lovd_syntax_check,
        shorthand_vcf=args.shorthand_vcf,
    )

    yield result

    LOGGER.info(
        "Validation complete"
    )


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------


def main() -> int:
    """
    VariantValidator command-line entry point.
    """

    parser = build_parser()

    try:

        args = parser.parse_args()

        if args.output_format is None:

            if args.output is sys.stdout:
                args.output_format = "json"
            else:
                args.output_format = "table"

        configure_logging(
            args.log_level,
        )

        start = time.perf_counter()

        LOGGER.info(
            "VariantValidator %s",
            __version__,
        )

        LOGGER.debug(
            "Python %s",
            sys.version.split()[0],
        )

        validator = Validator()

        LOGGER.debug(
            "Validator initialised"
        )

        variants = parse_variants(
            args.variant,
        )

        args.select_transcripts = (
            parse_transcripts(
                args.select_transcripts,
            )
        )

        #
        # Choose output style
        #

        writer = Writer(
            handle=args.output,
            output_format=args.output_format,
            with_meta=args.meta if args.output is sys.stdout else True,
            repeat_header=args.output is sys.stdout,
        )

        #
        # Stream results directly to the writer
        #

        for result in run_validation(
                validator=validator,
                variants=variants,
                args=args,
        ):
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