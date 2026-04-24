import os
from configparser import ConfigParser

CONFIG_DIR = os.path.join(os.path.expanduser('~'), '.variantvalidator')

# Change settings based on config
config = ConfigParser()
config.read(CONFIG_DIR)

try:
    if config['logging']['file_name']:
        LOG_FILE = config['logging']['file_name']
    else:
        LOG_FILE = os.path.join(os.path.expanduser('~'), '.vv_errorlog')
except KeyError:
    LOG_FILE = os.path.join(os.path.expanduser('~'), '.vv_errorlog')

LOGGING_CONFIG = {
    'version': 1,
    'formatters': {
        'simple': {
            'class': 'logging.Formatter',
            # Console: clean, readable
            # Example:
            # 2026-04-24 17:11:02 INFO     VariantValidator.validate Started validation
            'format': '%(asctime)s %(levelname)-8s %(name)s %(message)s',
            'datefmt': '%Y-%m-%d %H:%M:%S',
        },
        'detailed': {
            'class': 'logging.Formatter',
            # File: richer but still human
            # Example:
            # 2026-04-24 17:11:02 VariantValidator.validate validate (line 412) INFO Started validation
            'format': (
                '%(asctime)s %(name)s %(funcName)s '
                '(line %(lineno)d) %(levelname)-8s %(message)s'
            ),
            'datefmt': '%Y-%m-%d %H:%M:%S',
        }
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'level': 'DEBUG',
            'formatter': 'simple'
        },
        'file': {
            'class': 'logging.FileHandler',
            'level': 'ERROR',
            'filename': LOG_FILE,
            'mode': 'a',
            'formatter': 'detailed',
        },
    },
    'loggers': {
        'VariantValidator': {
            'level': 'DEBUG',
            'handlers': ['console', 'file'],
            'propagate': 'no',
        }
    }
}

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
