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
    Restore the original logging configuration after every test.
    """
    original = copy.deepcopy(settings.LOGGING_CONFIG)

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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>
