"""
Tests for VariantValidator.bin.gene2transcripts.
"""

import argparse
import json
import sys

import pytest

import VariantValidator
from VariantValidator.bin import gene2transcripts

def test_parse_single_query():

    parsed = gene2transcripts._parse_value(
        "COL1A1"
    )

    assert parsed == [
        "COL1A1"
    ]


def test_parse_pipe_delimited():

    parsed = gene2transcripts._parse_value(
        "COL1A1|BRCA2|P3H1"
    )

    assert parsed == [
        "COL1A1",
        "BRCA2",
        "P3H1",
    ]


def test_parse_json_array():

    parsed = gene2transcripts._parse_value(
        '["COL1A1","BRCA2","P3H1"]'
    )

    assert parsed == [
        "COL1A1",
        "BRCA2",
        "P3H1",
    ]


def test_parse_empty_string():

    assert gene2transcripts._parse_value("") == []


def test_parse_invalid_json():

    with pytest.raises(ValueError):
        gene2transcripts._parse_value(
            '["COL1A1"'
        )

def test_parse_queries_multiple():

    parsed = gene2transcripts.parse_queries(
        [
            "COL1A1",
            "BRCA2",
        ]
    )

    assert parsed == [
        "COL1A1",
        "BRCA2",
    ]


def test_parse_queries_pipe():

    parsed = gene2transcripts.parse_queries(
        [
            "COL1A1|BRCA2|P3H1",
        ]
    )

    assert parsed == [
        "COL1A1",
        "BRCA2",
        "P3H1",
    ]


def test_parse_queries_empty():

    with pytest.raises(ValueError):
        gene2transcripts.parse_queries([])


def test_parse_queries_missing_file():

    with pytest.raises(FileNotFoundError):
        gene2transcripts.parse_queries(
            ["@this_file_does_not_exist.txt"]
        )

def test_parse_queries_text_file(tmp_path):

    infile = tmp_path / "queries.txt"

    infile.write_text(
        "\n".join(
            [
                "COL1A1",
                "",
                "# comment",
                "BRCA2",
            ]
        )
    )

    parsed = gene2transcripts.parse_queries(
        [f"@{infile}"]
    )

    assert parsed == [
        "COL1A1",
        "BRCA2",
    ]


def test_parse_queries_json_file(tmp_path):

    infile = tmp_path / "queries.json"

    infile.write_text(
        json.dumps(
            [
                "COL1A1",
                "BRCA2",
            ]
        )
    )

    parsed = gene2transcripts.parse_queries(
        [f"@{infile}"]
    )

    assert parsed == [
        "COL1A1",
        "BRCA2",
    ]


def test_load_json_file_invalid(tmp_path):

    infile = tmp_path / "queries.json"

    infile.write_text(
        json.dumps(
            {
                "query": "COL1A1"
            }
        )
    )

    with pytest.raises(ValueError):
        gene2transcripts._load_json_file(
            infile
        )

@pytest.mark.parametrize(
    "keyword",
    [
        "mane",
        "mane_select",
        "all",
        "raw",
        "select",
    ],
)
def test_parse_transcript_keyword(keyword):

    assert (
        gene2transcripts.parse_transcripts(keyword)
        == keyword
    )


def test_parse_transcript_pipe():

    parsed = gene2transcripts.parse_transcripts(
        "NM_000059.3|NM_000059.4"
    )

    assert parsed == [
        "NM_000059.3",
        "NM_000059.4",
    ]


def test_parse_transcript_json():

    parsed = gene2transcripts.parse_transcripts(
        '["NM_000059.3","NM_000059.4"]'
    )

    assert parsed == [
        "NM_000059.3",
        "NM_000059.4",
    ]


def test_parse_transcript_empty():

    with pytest.raises(ValueError):
        gene2transcripts.parse_transcripts("")

def test_writer_stdout(capsys):

    writer = gene2transcripts.Writer(
        handle=sys.stdout,
    )

    writer.write(
        {
            "COL1A1": {
                "current_symbol": "COL1A1",
            }
        }
    )

    captured = capsys.readouterr()

    assert '"COL1A1"' in captured.out
    assert '"current_symbol"' in captured.out


def test_writer_file(tmp_path):

    outfile = tmp_path / "output.json"

    with outfile.open("w") as handle:

        writer = gene2transcripts.Writer(
            handle=handle,
        )

        writer.write(
            {
                "COL1A1": {
                    "current_symbol": "COL1A1",
                }
            }
        )

        writer.close()

    data = json.loads(
        outfile.read_text()
    )

    assert (
        data["COL1A1"]["current_symbol"]
        == "COL1A1"
    )


def test_writer_close():

    from io import StringIO

    handle = StringIO()

    writer = gene2transcripts.Writer(
        handle=handle,
    )

    writer.close()

    assert handle.closed is False

class FakeValidator:

    def __init__(self):
        self.args = None
        self.kwargs = None

    def gene2transcripts(
        self,
        *args,
        **kwargs,
    ):

        self.args = args
        self.kwargs = kwargs

        return {
            "COL1A1": {
                "current_symbol": "COL1A1",
            }
        }


def test_run_query():

    validator = FakeValidator()

    args = argparse.Namespace(
        genome="GRCh38",
        transcript_model="refseq",
        select_transcripts="mane_select",
        no_web_searches=False,
        no_genomic_spans=False,
        lovd_syntax_check=False,
    )

    result = gene2transcripts.run_query(
        validator=validator,
        queries=[
            "COL1A1",
        ],
        args=args,
    )

    assert "COL1A1" in result

    assert validator.args[0] == '["COL1A1"]'

    assert validator.kwargs["genome_build"] == "GRCh38"

    assert validator.kwargs["transcript_set"] == "refseq"

    assert (
        validator.kwargs["select_transcripts"]
        == "mane_select"
    )

    assert (
        validator.kwargs["bypass_web_searches"]
        is False
    )

    assert (
        validator.kwargs["bypass_genomic_spans"]
        is False
    )

    assert (
        validator.kwargs["lovd_syntax_check"]
        is False
    )


def test_run_query_all_options():

    validator = FakeValidator()

    args = argparse.Namespace(
        genome="GRCh37",
        transcript_model="ensembl",
        select_transcripts=[
            "NM_000059.3",
            "NM_000059.4",
        ],
        no_web_searches=True,
        no_genomic_spans=True,
        lovd_syntax_check=True,
    )

    gene2transcripts.run_query(
        validator=validator,
        queries=[
            "COL1A1",
            "BRCA2",
        ],
        args=args,
    )

    assert validator.args[0] == (
        '["COL1A1", "BRCA2"]'
    )

    assert (
        validator.kwargs["transcript_set"]
        == "ensembl"
    )

    assert (
        validator.kwargs["select_transcripts"]
        == [
            "NM_000059.3",
            "NM_000059.4",
        ]
    )

    assert (
        validator.kwargs["bypass_web_searches"]
        is True
    )

    assert (
        validator.kwargs["bypass_genomic_spans"]
        is True
    )

    assert (
        validator.kwargs["lovd_syntax_check"]
        is True
    )

def test_main_success(tmp_path, monkeypatch):

    outfile = tmp_path / "results.json"

    monkeypatch.setattr(
        VariantValidator,
        "Validator",
        FakeValidator,
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gene2transcripts",
            "-q",
            "COL1A1",
            "-o",
            str(outfile),
        ],
    )

    rc = gene2transcripts.main()

    assert rc == gene2transcripts.EXIT_SUCCESS

    data = json.loads(
        outfile.read_text()
    )

    assert "COL1A1" in data


def test_main_bad_input(monkeypatch):

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gene2transcripts",
        ],
    )

    with pytest.raises(SystemExit):
        gene2transcripts.main()


def test_main_input_error(tmp_path, monkeypatch):

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gene2transcripts",
            "-q",
            "@missing.txt",
        ],
    )

    rc = gene2transcripts.main()

    assert rc == (
        gene2transcripts.EXIT_INPUT_ERROR
    )


def test_main_keyboard_interrupt(monkeypatch):

    class InterruptValidator:

        def __init__(self):
            raise KeyboardInterrupt

    monkeypatch.setattr(
        VariantValidator,
        "Validator",
        InterruptValidator,
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gene2transcripts",
            "-q",
            "COL1A1",
        ],
    )

    rc = gene2transcripts.main()

    assert rc == (
        gene2transcripts.EXIT_UNEXPECTED_ERROR
    )


def test_main_unexpected_exception(monkeypatch):

    class BrokenValidator:

        def __init__(self):
            raise RuntimeError(
                "Boom"
            )

    monkeypatch.setattr(
        VariantValidator,
        "Validator",
        BrokenValidator,
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gene2transcripts",
            "-q",
            "COL1A1",
        ],
    )

    rc = gene2transcripts.main()

    assert rc == (
        gene2transcripts.EXIT_UNEXPECTED_ERROR
    )

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
