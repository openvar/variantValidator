"""
Tests for the VariantValidator command-line interface.
"""

import json
import sys
import io

import pytest

from VariantValidator.bin import variantvalidator


# ----------------------------------------------------------------------
# Test variants
# ----------------------------------------------------------------------

GRCH38_GENOMIC = "NC_000017.11:g.50198002C>A"
GRCH37_GENOMIC = "NC_000017.10:g.48275363C>A"

NM4 = "NM_000088.4:c.589G>T"
NM3 = "NM_000088.3:c.589G>T"
NM2 = "NM_000088.2:c.589G>T"

ENST10 = "ENST00000225964.10:c.589G>T"
ENST5 = "ENST00000225964.5:c.589G>T"

ENST_NC = "ENST00000495677.1:n.316G>T"


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def run_cli(
    monkeypatch,
    capsys,
    *args,
    parse_output=True,
):
    """
    Execute the VariantValidator CLI.
    """

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantvalidator",
            "--log-level",
            "CRITICAL",
            *args,
        ],
    )

    exit_code = variantvalidator.main()

    captured = capsys.readouterr()

    assert exit_code == 0

    if not parse_output:
        return captured

    return json.loads(captured.out)


def get_primary_entry(results):
    """
    Return the first validation entry, ignoring metadata keys.
    """

    for key, value in results.items():

        if key in {
            "flag",
            "metadata",
        }:
            continue

        if isinstance(value, dict):
            return key, value

    pytest.fail("No validation entry found")


# ----------------------------------------------------------------------
# Basic validation
# ----------------------------------------------------------------------

def test_grch38_genomic(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        GRCH38_GENOMIC,
    )

    assert results["flag"] == "gene_variant"

    _, entry = get_primary_entry(results)

    assert entry["submitted_variant"] == GRCH38_GENOMIC
    assert entry["gene_symbol"] == "COL1A1"
    assert entry["hgvs_transcript_variant"] == NM4


def test_grch37_genomic(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        GRCH37_GENOMIC,
        "--genome",
        "GRCh37",
    )

    assert results["flag"] == "gene_variant"

    _, entry = get_primary_entry(results)

    assert entry["submitted_variant"] == GRCH37_GENOMIC
    assert entry["gene_symbol"] == "COL1A1"
    assert entry["hgvs_transcript_variant"] == NM4


def test_refseq_transcript_latest(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        NM4,
        "--select-transcripts",
        "all",
    )

    assert results["flag"] == "gene_variant"

    _, entry = get_primary_entry(results)

    assert entry["submitted_variant"] == NM4
    assert entry["gene_symbol"] == "COL1A1"
    assert entry["hgvs_transcript_variant"] == NM4


def test_refseq_transcript_previous(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        NM3,
        "--select-transcripts",
        "all",
    )

    assert results["flag"] == "gene_variant"

    _, entry = get_primary_entry(results)

    assert entry["submitted_variant"] == NM3
    assert entry["gene_symbol"] == "COL1A1"
    assert entry["hgvs_transcript_variant"] == NM3

# ----------------------------------------------------------------------
# Input formats
# ----------------------------------------------------------------------

def test_multiple_arguments(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        NM4,
        NM3,
        "--select-transcripts",
        "all",
    )

    assert results["flag"] == "gene_variant"
    assert NM4 in results
    assert NM3 in results


def test_pipe_delimited(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        f"{NM4}|{NM3}",
        "--select-transcripts",
        "all",
    )

    assert results["flag"] == "gene_variant"
    assert NM4 in results
    assert NM3 in results


def test_json_array(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        json.dumps([NM4, NM3]),
        "--select-transcripts",
        "all",
    )

    assert results["flag"] == "gene_variant"
    assert NM4 in results
    assert NM3 in results


def test_text_file(tmp_path, monkeypatch, capsys):

    variants = tmp_path / "variants.txt"

    variants.write_text(
        f"{NM4}\n"
        f"{NM3}\n"
    )

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        f"@{variants}",
        "--select-transcripts",
        "all",
    )

    assert results["flag"] == "gene_variant"
    assert NM4 in results
    assert NM3 in results


def test_json_file(tmp_path, monkeypatch, capsys):

    variants = tmp_path / "variants.json"

    variants.write_text(
        json.dumps([NM4, NM3])
    )

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        f"@{variants}",
        "--select-transcripts",
        "all",
    )

    assert results["flag"] == "gene_variant"
    assert NM4 in results
    assert NM3 in results

# ----------------------------------------------------------------------
# Transcript selection
# ----------------------------------------------------------------------

def test_select_transcripts_mane(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        GRCH38_GENOMIC,
        "--select-transcripts",
        "mane",
    )

    assert results["flag"] == "gene_variant"

    _, entry = get_primary_entry(results)

    assert entry["gene_symbol"] == "COL1A1"
    assert entry["hgvs_transcript_variant"] == NM4


def test_select_transcripts_mane_select(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        GRCH38_GENOMIC,
    )

    assert results["flag"] == "gene_variant"

    _, entry = get_primary_entry(results)

    assert entry["gene_symbol"] == "COL1A1"
    assert entry["hgvs_transcript_variant"] == NM4


def test_select_transcripts_all(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        GRCH38_GENOMIC,
        "--select-transcripts",
        "all",
    )

    assert results["flag"] == "gene_variant"

    assert NM4 in results
    assert NM3 not in results


def test_select_transcripts_raw(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        GRCH38_GENOMIC,
        "--genome",
        "GRCh37",
        "--select-transcripts",
        "raw",
    )

    assert results["flag"] == "gene_variant"

    assert NM4 in results
    assert NM3 in results
    assert NM2 in results


def test_select_specific_transcripts(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        GRCH38_GENOMIC,
        "--select-transcripts",
        '["NM_000088.4","NM_000088.3"]',
    )

    assert results["flag"] == "gene_variant"

    assert NM4 in results
    assert NM3 in results


def test_transcript_set_ensembl(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        ENST10,
        "--transcript-set",
        "ensembl",
        "--select-transcripts",
        "all",
    )

    assert results["flag"] == "gene_variant"

    _, entry = get_primary_entry(results)

    assert entry["gene_symbol"] == "COL1A1"
    assert entry["hgvs_transcript_variant"] == ENST10


def test_transcript_set_refseq(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        NM4,
        "--transcript-set",
        "refseq",
        "--select-transcripts",
        "all",
    )

    assert results["flag"] == "gene_variant"

    _, entry = get_primary_entry(results)

    assert entry["gene_symbol"] == "COL1A1"
    assert entry["hgvs_transcript_variant"] == NM4

# ----------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------

def test_output_json_file(tmp_path, monkeypatch, capsys):

    output = tmp_path / "results.json"

    run_cli(
        monkeypatch,
        capsys,
        "-v",
        NM4,
        "-f",
        "json",
        "-o",
        str(output),
        parse_output=False,
    )

    assert output.exists()

    results = json.loads(output.read_text())

    assert results["flag"] == "gene_variant"

    _, entry = get_primary_entry(results)

    assert entry["submitted_variant"] == NM4
    assert entry["gene_symbol"] == "COL1A1"


def test_output_table_file(tmp_path, monkeypatch, capsys):

    output = tmp_path / "results.tsv"

    run_cli(
        monkeypatch,
        capsys,
        "-v",
        NM4,
        "-f",
        "table",
        "-o",
        str(output),
        parse_output=False,
    )

    assert output.exists()

    text = output.read_text()

    assert "HGVS_transcript" in text
    assert "COL1A1" in text
    assert NM4 in text


def test_stdout_json(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        NM4,
        "-f",
        "json",
    )

    assert results["flag"] == "gene_variant"

    _, entry = get_primary_entry(results)

    assert entry["submitted_variant"] == NM4


def test_stdout_table(monkeypatch, capsys):

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantvalidator",
            "--log-level",
            "CRITICAL",
            "-v",
            NM4,
            "-f",
            "table",
        ],
    )

    exit_code = variantvalidator.main()

    captured = capsys.readouterr()

    assert exit_code == 0
    assert "HGVS_transcript" in captured.out
    assert "COL1A1" in captured.out
    assert NM4 in captured.out


def test_output_with_metadata(tmp_path, monkeypatch, capsys):

    output = tmp_path / "results.json"

    run_cli(
        monkeypatch,
        capsys,
        "-v",
        NM4,
        "-m",
        "-f",
        "json",
        "-o",
        str(output),
        parse_output=False,
    )

    results = json.loads(output.read_text())

    assert "metadata" in results
    assert "variantvalidator_version" in results["metadata"]

# ----------------------------------------------------------------------
# Validation options
# ----------------------------------------------------------------------

def test_no_liftover(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        GRCH38_GENOMIC,
        "--liftover-level",
        "false",
    )

    assert results["flag"] == "gene_variant"

    _, entry = get_primary_entry(results)

    assert entry["submitted_variant"] == GRCH38_GENOMIC
    assert entry["gene_symbol"] == "COL1A1"


def test_no_lovd_syntax_check(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        GRCH38_GENOMIC,
        "--no-lovd-syntax-check",
    )

    assert results["flag"] == "gene_variant"

    _, entry = get_primary_entry(results)

    assert entry["submitted_variant"] == GRCH38_GENOMIC
    assert entry["gene_symbol"] == "COL1A1"


def test_shorthand_vcf(monkeypatch, capsys):

    results = run_cli(
        monkeypatch,
        capsys,
        "-v",
        GRCH38_GENOMIC,
        "--shorthand-vcf",
    )

    assert results["flag"] == "gene_variant"

    _, entry = get_primary_entry(results)

    assert entry["submitted_variant"] == GRCH38_GENOMIC
    assert entry["gene_symbol"] == "COL1A1"


# ----------------------------------------------------------------------
# Miscellaneous
# ----------------------------------------------------------------------

def test_version(monkeypatch, capsys):

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantvalidator",
            "--version",
        ],
    )

    with pytest.raises(SystemExit) as exc:
        variantvalidator.main()

    assert exc.value.code == 0

    captured = capsys.readouterr()

    assert "VariantValidator" in captured.out


def test_help(monkeypatch, capsys):

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantvalidator",
            "--help",
        ],
    )

    with pytest.raises(SystemExit) as exc:
        variantvalidator.main()

    assert exc.value.code == 0

    captured = capsys.readouterr()

    assert "usage:" in captured.out
    assert "--variant" in captured.out
    assert "--select-transcripts" in captured.out

# ----------------------------------------------------------------------
# Error handling
# ----------------------------------------------------------------------

def test_missing_input_file(monkeypatch):

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantvalidator",
            "--log-level",
            "CRITICAL",
            "-v",
            "@/this/file/does/not/exist.txt",
        ],
    )

    assert (
        variantvalidator.main()
        == variantvalidator.EXIT_INPUT_ERROR
    )


def test_invalid_json_argument(monkeypatch):

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantvalidator",
            "--log-level",
            "CRITICAL",
            "-v",
            '["NM_000088.4:c.589G>T"',
        ],
    )

    assert (
        variantvalidator.main()
        == variantvalidator.EXIT_INPUT_ERROR
    )


def test_invalid_json_file(tmp_path, monkeypatch):

    infile = tmp_path / "variants.json"

    infile.write_text(
        '{"variant":"NM_000088.4:c.589G>T"}'
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantvalidator",
            "--log-level",
            "CRITICAL",
            "-v",
            f"@{infile}",
        ],
    )

    assert (
        variantvalidator.main()
        == variantvalidator.EXIT_INPUT_ERROR
    )


def test_empty_text_file(tmp_path, monkeypatch):

    infile = tmp_path / "variants.txt"

    infile.write_text("")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantvalidator",
            "--log-level",
            "CRITICAL",
            "-v",
            f"@{infile}",
        ],
    )

    assert (
        variantvalidator.main()
        == variantvalidator.EXIT_INPUT_ERROR
    )


def test_empty_json_array(monkeypatch):

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantvalidator",
            "--log-level",
            "CRITICAL",
            "-v",
            "[]",
        ],
    )

    assert (
        variantvalidator.main()
        == variantvalidator.EXIT_INPUT_ERROR
    )


def test_invalid_transcript_json(monkeypatch):

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantvalidator",
            "--log-level",
            "CRITICAL",
            "-v",
            NM4,
            "--select-transcripts",
            '["NM_000088.4"',
        ],
    )

    assert (
        variantvalidator.main()
        == variantvalidator.EXIT_INPUT_ERROR
    )


def test_keyboard_interrupt(monkeypatch):

    monkeypatch.setattr(
        variantvalidator,
        "run_validation",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            KeyboardInterrupt
        ),
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantvalidator",
            "--log-level",
            "CRITICAL",
            "-v",
            NM4,
        ],
    )

    assert (
        variantvalidator.main()
        == variantvalidator.EXIT_UNEXPECTED_ERROR
    )


def test_unexpected_exception(monkeypatch):

    monkeypatch.setattr(
        variantvalidator,
        "run_validation",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            RuntimeError("Boom")
        ),
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantvalidator",
            "--log-level",
            "CRITICAL",
            "-v",
            NM4,
        ],
    )

    assert (
        variantvalidator.main()
        == variantvalidator.EXIT_UNEXPECTED_ERROR
    )

def test_unknown_writer_format():

    writer = variantvalidator.Writer(
        handle=io.StringIO(),
        output_format="unknown",
    )

    with pytest.raises(ValueError):
        writer.write(None)

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
