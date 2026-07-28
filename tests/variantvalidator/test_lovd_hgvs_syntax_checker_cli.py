"""
Tests for VariantValidator.bin.lovd_hgvs_syntax_checker.
"""

import json
import sys
from io import StringIO
from unittest.mock import patch

import pytest

from VariantValidator.bin import lovd_hgvs_syntax_checker


# ----------------------------------------------------------------------
# Writer
# ----------------------------------------------------------------------


def test_writer_stdout(capsys):

    writer = lovd_hgvs_syntax_checker.Writer(
        handle=sys.stdout,
    )

    writer.write(
        {
            "result": "Valid",
            "version": "1.2.3",
        }
    )

    captured = capsys.readouterr()

    assert '"result"' in captured.out
    assert '"Valid"' in captured.out
    assert '"version"' in captured.out


def test_writer_file(tmp_path):

    outfile = tmp_path / "result.json"

    with outfile.open("w") as handle:

        writer = lovd_hgvs_syntax_checker.Writer(
            handle=handle,
        )

        writer.write(
            {
                "result": "Valid",
                "version": "1.2.3",
            }
        )

        writer.close()

    data = json.loads(
        outfile.read_text()
    )

    assert data["result"] == "Valid"
    assert data["version"] == "1.2.3"


def test_writer_close():

    handle = StringIO()

    writer = lovd_hgvs_syntax_checker.Writer(
        handle=handle,
    )

    writer.close()

    assert handle.closed is False


# ----------------------------------------------------------------------
# run_query
# ----------------------------------------------------------------------


@patch(
    "VariantValidator.modules.lovd_api.lovd_syntax_check"
)
def test_run_query_variant(
    mock_check,
):

    mock_check.return_value = {
        "result": "Valid",
    }

    result = lovd_hgvs_syntax_checker.run_query(
        query="NM_000059.4:c.7790G>A",
        is_a_gene=False,
    )

    assert result == {
        "result": "Valid",
    }

    mock_check.assert_called_once_with(
        variant_description="NM_000059.4:c.7790G>A",
        is_a_gene=False,
    )


@patch(
    "VariantValidator.modules.lovd_api.lovd_syntax_check"
)
def test_run_query_gene(
    mock_check,
):

    mock_check.return_value = {
        "result": "Valid",
    }

    result = lovd_hgvs_syntax_checker.run_query(
        query="BRCA2",
        is_a_gene=True,
    )

    assert result == {
        "result": "Valid",
    }

    mock_check.assert_called_once_with(
        variant_description="BRCA2",
        is_a_gene=True,
    )

# ----------------------------------------------------------------------
# Parser
# ----------------------------------------------------------------------


def test_build_parser_defaults():

    parser = lovd_hgvs_syntax_checker.build_parser()

    args = parser.parse_args(
        [
            "-q",
            "NM_000059.4:c.7790G>A",
        ]
    )

    assert (
        args.query
        == "NM_000059.4:c.7790G>A"
    )

    assert args.gene is False

    assert (
        args.log_level
        == "WARNING"
    )

    assert args.output is sys.stdout


def test_build_parser_gene():

    parser = lovd_hgvs_syntax_checker.build_parser()

    args = parser.parse_args(
        [
            "-q",
            "BRCA2",
            "--gene",
        ]
    )

    assert args.query == "BRCA2"

    assert args.gene is True


def test_build_parser_output(tmp_path):

    outfile = tmp_path / "result.json"

    parser = lovd_hgvs_syntax_checker.build_parser()

    args = parser.parse_args(
        [
            "-q",
            "BRCA2",
            "-o",
            str(outfile),
        ]
    )

    assert args.query == "BRCA2"

    assert args.output.name == str(outfile)

    args.output.close()


def test_build_parser_version():

    parser = lovd_hgvs_syntax_checker.build_parser()

    with pytest.raises(SystemExit):

        parser.parse_args(
            [
                "--version",
            ]
        )

# ----------------------------------------------------------------------
# main()
# ----------------------------------------------------------------------


@patch(
    "VariantValidator.modules.lovd_api.lovd_syntax_check"
)
def test_main_success_stdout(
    mock_check,
    monkeypatch,
    capsys,
):

    mock_check.return_value = {
        "result": "Valid",
    }

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "lovd-hgvs-syntax-checker",
            "-q",
            "NM_000059.4:c.7790G>A",
        ],
    )

    rc = lovd_hgvs_syntax_checker.main()

    assert rc == (
        lovd_hgvs_syntax_checker.EXIT_SUCCESS
    )

    captured = capsys.readouterr()

    data = json.loads(
        captured.out
    )

    assert data["result"] == "Valid"


@patch(
    "VariantValidator.modules.lovd_api.lovd_syntax_check"
)
def test_main_success_file(
    mock_check,
    tmp_path,
    monkeypatch,
):

    mock_check.return_value = {
        "result": "Valid",
    }

    outfile = tmp_path / "result.json"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "lovd-hgvs-syntax-checker",
            "-q",
            "NM_000059.4:c.7790G>A",
            "-o",
            str(outfile),
        ],
    )

    rc = lovd_hgvs_syntax_checker.main()

    assert rc == (
        lovd_hgvs_syntax_checker.EXIT_SUCCESS
    )

    data = json.loads(
        outfile.read_text()
    )

    assert data["result"] == "Valid"


def test_main_bad_input(
    monkeypatch,
):

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "lovd-hgvs-syntax-checker",
        ],
    )

    with pytest.raises(
        SystemExit
    ):
        lovd_hgvs_syntax_checker.main()


@patch(
    "VariantValidator.modules.lovd_api.lovd_syntax_check"
)
def test_main_input_error(
    mock_check,
    monkeypatch,
):

    mock_check.side_effect = ValueError(
        "Invalid HGVS"
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "lovd-hgvs-syntax-checker",
            "-q",
            "bad_variant",
        ],
    )

    rc = lovd_hgvs_syntax_checker.main()

    assert rc == (
        lovd_hgvs_syntax_checker.EXIT_INPUT_ERROR
    )


@patch(
    "VariantValidator.modules.lovd_api.lovd_syntax_check"
)
def test_main_keyboard_interrupt(
    mock_check,
    monkeypatch,
):

    mock_check.side_effect = KeyboardInterrupt

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "lovd-hgvs-syntax-checker",
            "-q",
            "NM_000059.4:c.7790G>A",
        ],
    )

    rc = lovd_hgvs_syntax_checker.main()

    assert rc == (
        lovd_hgvs_syntax_checker.EXIT_UNEXPECTED_ERROR
    )


@patch(
    "VariantValidator.modules.lovd_api.lovd_syntax_check"
)
def test_main_unexpected_exception(
    mock_check,
    monkeypatch,
):

    mock_check.side_effect = RuntimeError(
        "Boom"
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "lovd-hgvs-syntax-checker",
            "-q",
            "NM_000059.4:c.7790G>A",
        ],
    )

    rc = lovd_hgvs_syntax_checker.main()

    assert rc == (
        lovd_hgvs_syntax_checker.EXIT_UNEXPECTED_ERROR
    )

@patch("VariantValidator.logger.configure_logging")
def test_configure_logging(mock_configure):

    lovd_hgvs_syntax_checker.configure_logging(
        "INFO"
    )

    mock_configure.assert_called_once_with(
        console_level="INFO",
    )


@patch("VariantValidator.logger.configure_logging")
def test_configure_logging_default(mock_configure):

    lovd_hgvs_syntax_checker.configure_logging(
        None
    )

    mock_configure.assert_called_once_with(
        console_level=None,
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
