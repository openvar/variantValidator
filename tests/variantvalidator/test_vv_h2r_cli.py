"""
Tests for VariantValidator.bin.hgvs2reference.
"""

import json
import sys
from io import StringIO

import pytest

import VariantValidator
from VariantValidator.bin import hgvs2reference


def test_writer_stdout(capsys):

    writer = hgvs2reference.Writer(
        handle=sys.stdout,
    )

    writer.write(
        {
            "variant": "NM_000546.6:c.215C>G",
            "sequence": "ACTG",
            "error": "",
        }
    )

    captured = capsys.readouterr()

    assert '"variant"' in captured.out
    assert '"sequence"' in captured.out
    assert '"ACTG"' in captured.out


def test_writer_file(tmp_path):

    outfile = tmp_path / "output.json"

    with outfile.open("w") as handle:

        writer = hgvs2reference.Writer(
            handle=handle,
        )

        writer.write(
            {
                "variant": "NM_000546.6:c.215C>G",
                "sequence": "ACTG",
                "error": "",
            }
        )

        writer.close()

    data = json.loads(
        outfile.read_text()
    )

    assert (
        data["variant"]
        == "NM_000546.6:c.215C>G"
    )

    assert (
        data["sequence"]
        == "ACTG"
    )


def test_writer_close():

    handle = StringIO()

    writer = hgvs2reference.Writer(
        handle=handle,
    )

    writer.close()

    assert handle.closed is False


class FakeValidator:

    def __init__(self):
        self.query = None

    def hgvs2ref(
        self,
        query,
    ):

        self.query = query

        return {
            "variant": query,
            "start_position": "215",
            "end_position": "215",
            "sequence": "ACTG",
            "warning": "",
            "error": "",
        }


def test_run_query():

    validator = FakeValidator()

    result = hgvs2reference.run_query(
        validator=validator,
        query="NM_000546.6:c.215C>G",
    )

    assert (
        validator.query
        == "NM_000546.6:c.215C>G"
    )

    assert (
        result["variant"]
        == "NM_000546.6:c.215C>G"
    )

    assert (
        result["sequence"]
        == "ACTG"
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
            "hgvs2reference",
            "-q",
            "NM_000546.6:c.215C>G",
            "-o",
            str(outfile),
        ],
    )

    rc = hgvs2reference.main()

    assert (
        rc
        == hgvs2reference.EXIT_SUCCESS
    )

    data = json.loads(
        outfile.read_text()
    )

    assert (
        data["variant"]
        == "NM_000546.6:c.215C>G"
    )


def test_main_bad_input(monkeypatch):

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "hgvs2reference",
        ],
    )

    with pytest.raises(SystemExit):
        hgvs2reference.main()


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
            "hgvs2reference",
            "-q",
            "NM_000546.6:c.215C>G",
        ],
    )

    rc = hgvs2reference.main()

    assert (
        rc
        == hgvs2reference.EXIT_UNEXPECTED_ERROR
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
            "hgvs2reference",
            "-q",
            "NM_000546.6:c.215C>G",
        ],
    )

    rc = hgvs2reference.main()

    assert (
        rc
        == hgvs2reference.EXIT_UNEXPECTED_ERROR
    )


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later