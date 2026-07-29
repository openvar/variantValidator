import importlib
import os
from unittest import TestCase
from unittest.mock import patch

from VariantValidator import settings
from VariantValidator.settings import (
    LOG_FILE,
    LOGGING_CONFIG,
    vvDB_GET_CACHE,
    vvDB_GET_CACHE_SIZE,
)


CONFIG_DIR = settings.get_config_dir()


class TestSettings(TestCase):
    def test_config_dir_exists(self):
        assert os.path.exists(CONFIG_DIR)

    def test_log_file_exists(self):
        assert os.path.exists(LOG_FILE)

    def test_logging_config_structure(self):
        assert isinstance(LOGGING_CONFIG, dict)
        assert "version" in LOGGING_CONFIG
        assert "formatters" in LOGGING_CONFIG
        assert "handlers" in LOGGING_CONFIG
        assert "loggers" in LOGGING_CONFIG

    def test_logging_config_handlers(self):
        handlers = LOGGING_CONFIG.get("handlers", {})
        assert "console" in handlers
        assert "file" in handlers

    def test_logging_config_console_handler(self):
        console_handler = LOGGING_CONFIG["handlers"].get(
            "console",
            {},
        )
        assert (
            console_handler.get("class")
            == "logging.StreamHandler"
        )
        assert console_handler.get("formatter") == "simple"

    def test_logging_config_file_handler(self):
        file_handler = LOGGING_CONFIG["handlers"].get(
            "file",
            {},
        )
        assert (
            file_handler.get("class")
            == "logging.handlers.RotatingFileHandler"
        )
        assert file_handler.get("filename") == LOG_FILE
        assert file_handler.get("mode") == "a"
        assert file_handler.get("formatter") == "detailed"

    def test_logging_config_loggers(self):
        loggers = LOGGING_CONFIG.get("loggers", {})
        assert "VariantValidator" in loggers

    def test_logging_config_variantvalidator_logger(self):
        vv_logger = LOGGING_CONFIG["loggers"].get(
            "VariantValidator",
            {},
        )
        assert vv_logger.get("handlers") == [
            "console",
            "file",
        ]
        assert vv_logger.get("propagate") is False

    def test_vvdb_get_cache_is_boolean(self):
        assert isinstance(vvDB_GET_CACHE, bool)

    def test_vvdb_get_cache_size_is_positive_integer(self):
        assert isinstance(vvDB_GET_CACHE_SIZE, int)
        assert vvDB_GET_CACHE_SIZE > 0

    def test_vvdb_get_cache_environment_true(self):
        values = (
            "true",
            "TRUE",
            "1",
            "yes",
            "YES",
        )

        for value in values:
            with self.subTest(value=value):
                with patch.dict(
                    os.environ,
                    {"VV_DB_GET_CACHE": value},
                    clear=False,
                ):
                    importlib.reload(settings)
                    assert settings.vvDB_GET_CACHE is True

                importlib.reload(settings)

    def test_vvdb_get_cache_environment_false(self):
        values = (
            "false",
            "FALSE",
            "0",
            "no",
            "NO",
        )

        for value in values:
            with self.subTest(value=value):
                with patch.dict(
                    os.environ,
                    {"VV_DB_GET_CACHE": value},
                    clear=False,
                ):
                    importlib.reload(settings)
                    assert settings.vvDB_GET_CACHE is False

                importlib.reload(settings)

    def test_vvdb_get_cache_size_environment(self):
        with patch.dict(
            os.environ,
            {"VV_DB_GET_CACHE_SIZE": "12345"},
            clear=False,
        ):
            importlib.reload(settings)
            assert settings.vvDB_GET_CACHE_SIZE == 12345

        importlib.reload(settings)


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
