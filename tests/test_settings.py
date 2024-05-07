import os
from unittest import TestCase
from VariantValidator.settings import CONFIG_DIR, LOG_FILE, LOGGING_CONFIG


class TestSettings(TestCase):
    def test_config_dir_exists(self):
        assert os.path.exists(CONFIG_DIR)

    def test_log_file_exists(self):
        assert os.path.exists(LOG_FILE)

    def test_logging_config_structure(self):
        assert isinstance(LOGGING_CONFIG, dict)
        assert 'version' in LOGGING_CONFIG
        assert 'formatters' in LOGGING_CONFIG
        assert 'handlers' in LOGGING_CONFIG
        assert 'loggers' in LOGGING_CONFIG

    def test_logging_config_handlers(self):
        handlers = LOGGING_CONFIG.get('handlers', {})
        assert 'console' in handlers
        assert 'file' in handlers

    def test_logging_config_console_handler(self):
        console_handler = LOGGING_CONFIG['handlers'].get('console', {})
        assert console_handler.get('class') == 'logging.StreamHandler'
        assert console_handler.get('formatter') == 'simple'

    def test_logging_config_file_handler(self):
        file_handler = LOGGING_CONFIG['handlers'].get('file', {})
        assert file_handler.get('class') == 'logging.FileHandler'
        assert file_handler.get('filename') == LOG_FILE
        assert file_handler.get('mode') == 'a'
        assert file_handler.get('formatter') == 'detailed'

    def test_logging_config_loggers(self):
        loggers = LOGGING_CONFIG.get('loggers', {})
        assert 'VariantValidator' in loggers

    def test_logging_config_variantvalidator_logger(self):
        vv_logger = LOGGING_CONFIG['loggers'].get('VariantValidator', {})
        assert vv_logger.get('handlers') == ['console', 'file']
        assert vv_logger.get('propagate') == 'no'


# <LICENSE>
# Copyright (C) 2016-2024 VariantValidator Contributors
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
