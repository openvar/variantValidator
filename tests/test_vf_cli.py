"""
Tests for VariantFormatter.bin.variantformatter.
"""

import json
import sys
import argparse

import pytest

from VariantValidator.bin import variantformatter


class FakeFormatter:

    def __init__(self):
        self.kwargs = None

    def format(self, **kwargs):

        self.kwargs = kwargs

        return {
            "17-50198002-C-A": {
                "status": "ok",
            }
        }

def test_parse_single_variant():

    parsed = variantformatter._parse_value(
        "NC_000017.10:g.48275363C>A"
    )

    assert parsed == [
        "NC_000017.10:g.48275363C>A"
    ]


def test_parse_pipe_delimited():

    parsed = variantformatter._parse_value(
        "17-50198002-C-A|17-50198003-G-T"
    )

    assert parsed == [
        "17-50198002-C-A",
        "17-50198003-G-T",
    ]


def test_parse_json_array():

    parsed = variantformatter._parse_value(
        '["17-50198002-C-A","17-50198003-G-T"]'
    )

    assert parsed == [
        "17-50198002-C-A",
        "17-50198003-G-T",
    ]


def test_parse_empty_string():

    assert variantformatter._parse_value("") == []


def test_parse_invalid_json():

    with pytest.raises(ValueError):
        variantformatter._parse_value(
            '["17-50198002-C-A"'
        )

def test_parse_variants_multiple():

    parsed = variantformatter.parse_variants(
        [
            "17-50198002-C-A",
            "17-50198003-G-T",
        ]
    )

    assert parsed == [
        "17-50198002-C-A",
        "17-50198003-G-T",
    ]


def test_parse_variants_pipe():

    parsed = variantformatter.parse_variants(
        [
            "17-50198002-C-A|17-50198003-G-T",
        ]
    )

    assert parsed == [
        "17-50198002-C-A",
        "17-50198003-G-T",
    ]


def test_parse_variants_empty():

    with pytest.raises(ValueError):
        variantformatter.parse_variants([])


def test_parse_variants_missing_file():

    with pytest.raises(FileNotFoundError):
        variantformatter.parse_variants(
            ["@this_file_does_not_exist.txt"]
        )

def test_parse_variants_text_file(tmp_path):

    infile = tmp_path / "variants.txt"

    infile.write_text(
        "\n".join(
            [
                "17-50198002-C-A",
                "",
                "# comment",
                "17-50198003-G-T",
            ]
        )
    )

    parsed = variantformatter.parse_variants(
        [f"@{infile}"]
    )

    assert parsed == [
        "17-50198002-C-A",
        "17-50198003-G-T",
    ]


def test_parse_variants_json_file(tmp_path):

    infile = tmp_path / "variants.json"

    infile.write_text(
        json.dumps(
            [
                "17-50198002-C-A",
                "17-50198003-G-T",
            ]
        )
    )

    parsed = variantformatter.parse_variants(
        [f"@{infile}"]
    )

    assert parsed == [
        "17-50198002-C-A",
        "17-50198003-G-T",
    ]


def test_load_json_file_invalid(tmp_path):

    infile = tmp_path / "variants.json"

    infile.write_text(
        json.dumps(
            {
                "variant": "17-50198002-C-A"
            }
        )
    )

    with pytest.raises(ValueError):
        variantformatter._load_json_file(infile)

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
        variantformatter.parse_transcripts(keyword)
        == keyword
    )


def test_parse_transcript_pipe():

    parsed = variantformatter.parse_transcripts(
        "NM_000093.4|NM_001278074.1"
    )

    assert parsed == [
        "NM_000093.4",
        "NM_001278074.1",
    ]


def test_parse_transcript_json():

    parsed = variantformatter.parse_transcripts(
        '["NM_000093.4","NM_001278074.1"]'
    )

    assert parsed == [
        "NM_000093.4",
        "NM_001278074.1",
    ]


def test_parse_transcript_empty():

    with pytest.raises(ValueError):
        variantformatter.parse_transcripts("")


def test_writer_stdout(capsys):

    writer = variantformatter.Writer(
        handle=sys.stdout,
    )

    writer.write(
        {
            "variant": {
                "gene": "COL1A1",
            }
        }
    )

    captured = capsys.readouterr()

    assert '"variant"' in captured.out
    assert '"gene"' in captured.out


def test_writer_file(tmp_path):

    outfile = tmp_path / "output.json"

    with outfile.open("w") as handle:

        writer = variantformatter.Writer(
            handle=handle,
        )

        writer.write(
            {
                "variant": {
                    "gene": "COL1A1",
                }
            }
        )

        writer.close()

    data = json.loads(
        outfile.read_text()
    )

    assert data["variant"]["gene"] == "COL1A1"


def test_writer_close():

    from io import StringIO

    handle = StringIO()

    writer = variantformatter.Writer(
        handle=handle,
    )

    writer.close()

    assert handle.closed is False


def test_run_formatting():

    formatter = FakeFormatter()

    args = argparse.Namespace(
        genome="GRCh38",
        transcript_model="refseq",
        select_transcripts="mane_select",
        check_only=False,
        liftover=False,
    )

    result = variantformatter.run_formatting(
        formatter=formatter,
        variants=[
            "17-50198002-C-A",
        ],
        args=args,
    )

    assert "17-50198002-C-A" in result

    assert formatter.kwargs["batch_input"] == (
        '["17-50198002-C-A"]'
    )

    assert formatter.kwargs["genome_build"] == "GRCh38"

    assert formatter.kwargs["transcript_model"] == "refseq"

    assert formatter.kwargs["specify_transcripts"] == "mane_select"

    assert formatter.kwargs["checkOnly"] is False

    assert formatter.kwargs["liftover"] is False


def test_run_formatting_transcript_list():

    formatter = FakeFormatter()

    args = argparse.Namespace(
        genome="GRCh37",
        transcript_model="ensembl",
        select_transcripts=[
            "NM_000093.4",
            "NM_001278074.1",
        ],
        check_only=True,
        liftover=True,
    )

    variantformatter.run_formatting(
        formatter=formatter,
        variants=[
            "17-50198002-C-A",
        ],
        args=args,
    )

    assert (
        formatter.kwargs["specify_transcripts"]
        ==
        '["NM_000093.4", "NM_001278074.1"]'
    )

    assert formatter.kwargs["checkOnly"] is True

    assert formatter.kwargs["liftover"] is True


def test_main_success(tmp_path, monkeypatch):

    outfile = tmp_path / "results.json"

    monkeypatch.setattr(
        variantformatter,
        "SimpleVariantFormatter",
        FakeFormatter,
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantformatter",
            "-v",
            "17-50198002-C-A",
            "-o",
            str(outfile),
        ],
    )

    rc = variantformatter.main()

    assert rc == variantformatter.EXIT_SUCCESS

    data = json.loads(outfile.read_text())

    assert "17-50198002-C-A" in data


def test_main_bad_input(monkeypatch):

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantformatter",
        ],
    )

    with pytest.raises(SystemExit):
        variantformatter.main()


def test_main_keyboard_interrupt(monkeypatch):

    class InterruptFormatter:

        def __init__(self):
            raise KeyboardInterrupt

    monkeypatch.setattr(
        variantformatter,
        "SimpleVariantFormatter",
        InterruptFormatter,
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantformatter",
            "-v",
            "17-50198002-C-A",
        ],
    )

    rc = variantformatter.main()

    assert rc == variantformatter.EXIT_UNEXPECTED_ERROR


def test_main_unexpected_exception(monkeypatch):

    class BrokenFormatter:

        def __init__(self):
            raise RuntimeError("Boom")

    monkeypatch.setattr(
        variantformatter,
        "SimpleVariantFormatter",
        BrokenFormatter,
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "variantformatter",
            "-v",
            "17-50198002-C-A",
        ],
    )

    rc = variantformatter.main()

    assert rc == variantformatter.EXIT_UNEXPECTED_ERROR

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
