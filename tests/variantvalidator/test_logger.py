"""
Tests for VariantValidator.logger.
"""

import copy
import logging

import pytest

from VariantValidator import settings
from VariantValidator.logger import configure_logging


@pytest.fixture(autouse=True)
def reset_logging_config():
    """
    Run each test with the default logging levels and restore
    the incoming logging configuration afterwards.
    """
    original = copy.deepcopy(settings.LOGGING_CONFIG)

    settings.LOGGING_CONFIG["handlers"]["console"]["level"] = (
        settings.CONSOLE_LEVEL
    )
    settings.LOGGING_CONFIG["handlers"]["file"]["level"] = (
        settings.FILE_LEVEL
    )

    yield

    settings.LOGGING_CONFIG = original


def test_configure_logging_defaults():
    """
    Default logging configuration is applied from settings.
    """
    configure_logging()

    logger = logging.getLogger("VariantValidator")

    assert logger.hasHandlers()

    assert (
        settings.LOGGING_CONFIG["handlers"]["console"]["level"]
        == settings.CONSOLE_LEVEL
    )

    assert (
        settings.LOGGING_CONFIG["handlers"]["file"]["level"]
        == settings.FILE_LEVEL
    )


def test_console_level_override():
    """
    Console logging level can be overridden.
    """
    configure_logging(console_level="WARNING")

    assert (
        settings.LOGGING_CONFIG["handlers"]["console"]["level"]
        == "WARNING"
    )

    assert (
        settings.LOGGING_CONFIG["handlers"]["file"]["level"]
        == settings.FILE_LEVEL
    )


def test_file_level_override():
    """
    File logging level can be overridden.
    """
    configure_logging(file_level="INFO")

    assert (
        settings.LOGGING_CONFIG["handlers"]["file"]["level"]
        == "INFO"
    )

    assert (
        settings.LOGGING_CONFIG["handlers"]["console"]["level"]
        == settings.CONSOLE_LEVEL
    )


def test_both_levels_override():
    """
    Console and file logging levels can both be overridden.
    """
    configure_logging(
        console_level="ERROR",
        file_level="CRITICAL",
    )

    assert (
        settings.LOGGING_CONFIG["handlers"]["console"]["level"]
        == "ERROR"
    )

    assert (
        settings.LOGGING_CONFIG["handlers"]["file"]["level"]
        == "CRITICAL"
    )


def test_rotating_file_handler_configuration():
    """
    File handler should be converted to a RotatingFileHandler.
    """
    configure_logging()

    handler = settings.LOGGING_CONFIG["handlers"]["file"]

    assert (
        handler["class"]
        == "logging.handlers.RotatingFileHandler"
    )

    assert handler["level"] == settings.FILE_LEVEL
    assert handler["filename"] == settings.LOG_FILE
    assert handler["mode"] == "a"
    assert handler["maxBytes"] == 500000
    assert handler["backupCount"] == 2
    assert handler["formatter"] == "detailed"


def test_console_handler_configuration():
    """
    Console handler configuration should be preserved.
    """
    configure_logging()

    handler = settings.LOGGING_CONFIG["handlers"]["console"]

    assert handler["class"] == "logging.StreamHandler"
    assert handler["formatter"] == "simple"
    assert handler["level"] == settings.CONSOLE_LEVEL


def test_logger_configuration():
    """
    VariantValidator logger should exist and use configured handlers.
    """
    configure_logging()

    logger = logging.getLogger("VariantValidator")

    assert logger.level == logging.DEBUG
    assert logger.propagate is False
    assert logger.hasHandlers()
    assert len(logger.handlers) == 2


def test_file_handler_filename_preserved():
    """
    Log file path should remain unchanged.
    """
    configure_logging()

    assert (
        settings.LOGGING_CONFIG["handlers"]["file"]["filename"]
        == settings.LOG_FILE
    )


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later