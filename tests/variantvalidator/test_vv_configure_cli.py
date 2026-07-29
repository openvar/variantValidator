"""
Tests for VariantValidator.bin.vv_configure.
"""

import configparser
from pathlib import Path

import pytest

from VariantValidator.bin import vv_configure


@pytest.fixture
def temp_config(tmp_path, monkeypatch):

    config_file = tmp_path / ".variantvalidator"

    monkeypatch.setenv(
        "VARIANTVALIDATOR_TEST_CONFIG",
        str(config_file),
    )

    yield config_file

"""
Tests for VariantValidator.bin.vv_configure.
"""

import configparser
from pathlib import Path

import pytest

from VariantValidator.bin import vv_configure


@pytest.fixture
def temp_config(tmp_path, monkeypatch):

    config_file = tmp_path / ".variantvalidator"

    monkeypatch.setenv(
        "VARIANTVALIDATOR_TEST_CONFIG",
        str(config_file),
    )

    yield config_file

def test_write_configuration(temp_config):

    config = configparser.ConfigParser()

    config["logging"] = {
        "console": "INFO",
    }

    vv_configure.write_configuration(
        config,
        temp_config,
    )

    assert temp_config.exists()

    check = configparser.ConfigParser()
    check.read(temp_config)

    assert check["logging"]["console"] == "INFO"

def test_load_existing_configuration(temp_config):

    config = configparser.ConfigParser()

    config["logging"] = {
        "console": "ERROR",
    }

    with open(temp_config, "w") as fh:
        config.write(fh)

    loaded, filename, new_file = (
        vv_configure.load_configuration()
    )

    assert new_file is False
    assert filename == temp_config
    assert loaded["logging"]["console"] == "ERROR"

def test_configure_changes_value(monkeypatch):

    config = configparser.ConfigParser()

    config["logging"] = {
        "console": "INFO",
    }

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: "WARNING",
    )

    changed = vv_configure.configure(config)

    assert changed is True
    assert config["logging"]["console"] == "WARNING"

def test_configure_blank_input(monkeypatch):

    config = configparser.ConfigParser()

    config["logging"] = {
        "console": "INFO",
    }

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: "",
    )

    changed = vv_configure.configure(config)

    assert changed is False
    assert config["logging"]["console"] == "INFO"

def test_configure_single_section(monkeypatch):

    config = configparser.ConfigParser()

    config["logging"] = {
        "console": "INFO",
    }

    config["mysql"] = {
        "host": "localhost",
    }

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: "WARNING",
    )

    vv_configure.configure(
        config,
        section="logging",
    )

    assert config["logging"]["console"] == "WARNING"
    assert config["mysql"]["host"] == "localhost"

def test_empty_configuration():

    config = configparser.ConfigParser()

    changed = vv_configure.configure(config)

    assert changed is False

def test_default_config_file_exists():

    path = vv_configure.default_config_file()

    assert path.name == "default.ini"

import sys


def test_main_help(capsys, monkeypatch):

    monkeypatch.setattr(
        sys,
        "argv",
        ["vv_configure", "--help"],
    )

    with pytest.raises(SystemExit) as exc:
        vv_configure.main()

    assert exc.value.code == 0

    captured = capsys.readouterr()

    assert "Configure VariantValidator" in captured.out


def test_main_new_configuration(
        temp_config,
        monkeypatch,
        capsys,
):

    monkeypatch.setattr(
        sys,
        "argv",
        ["vv_configure"],
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: "",
    )

    rc = vv_configure.main()

    captured = capsys.readouterr()

    assert rc == 0
    assert temp_config.exists()
    assert "Configuration written to" in captured.out


def test_main_existing_configuration_no_changes(
        temp_config,
        monkeypatch,
        capsys,
):

    config = configparser.ConfigParser()

    config["logging"] = {
        "console": "INFO",
    }

    with open(temp_config, "w") as handle:
        config.write(handle)

    monkeypatch.setattr(
        sys,
        "argv",
        ["vv_configure"],
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: "",
    )

    rc = vv_configure.main()

    captured = capsys.readouterr()

    assert rc == 0
    assert "No changes made." in captured.out


def test_main_existing_configuration_changed(
        temp_config,
        monkeypatch,
        capsys,
):

    config = configparser.ConfigParser()

    config["logging"] = {
        "console": "INFO",
    }

    with open(temp_config, "w") as handle:
        config.write(handle)

    monkeypatch.setattr(
        sys,
        "argv",
        ["vv_configure"],
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: "DEBUG",
    )

    rc = vv_configure.main()

    captured = capsys.readouterr()

    assert rc == 0
    assert "Configuration written to" in captured.out


def test_main_single_section(
        temp_config,
        monkeypatch,
        capsys,
):

    config = configparser.ConfigParser()

    config["logging"] = {
        "console": "INFO",
    }

    config["mysql"] = {
        "host": "localhost",
    }

    with open(temp_config, "w") as handle:
        config.write(handle)

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "vv_configure",
            "--section",
            "logging",
        ],
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: "",
    )

    rc = vv_configure.main()

    captured = capsys.readouterr()

    assert rc == 0
    assert "Section: logging" in captured.out
    assert "Section: mysql" not in captured.out

# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
