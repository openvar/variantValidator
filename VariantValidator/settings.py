import os
from configparser import ConfigParser

CONFIG_DIR = os.path.join(os.path.expanduser('~'), '.variantvalidator')

config = ConfigParser()
config.read(CONFIG_DIR)

# ----------------------------------------
# DEFAULT LOG FILE
# ----------------------------------------
DEFAULT_LOG = os.path.join(
    os.path.expanduser('~'),
    '.vv_errorlog'
)

# ----------------------------------------
# FILE LOCATION (SAFE)
# ----------------------------------------
if config.has_section('logging'):
    file_name = config.get('logging', 'file_name', fallback=None)
else:
    file_name = None

LOG_FILE = file_name if file_name else DEFAULT_LOG

# Normalise path (handles ~ and relative paths)
LOG_FILE = os.path.abspath(os.path.expanduser(LOG_FILE))


# ----------------------------------------
# ENSURE LOG DIRECTORY EXISTS
# ----------------------------------------
log_dir = os.path.dirname(LOG_FILE)
if log_dir:
    os.makedirs(log_dir, exist_ok=True)


# ----------------------------------------
# LOG LEVELS (SAFE)
# ----------------------------------------
CONSOLE_LEVEL = config.get('logging', 'console', fallback='DEBUG').upper()
FILE_LEVEL = config.get('logging', 'file', fallback='ERROR').upper()


# ----------------------------------------
# LOGGING CONFIG
# ----------------------------------------
LOGGING_CONFIG = {
    'version': 1,

    'formatters': {
        'simple': {
            'format': '%(asctime)s | %(levelname)-8s | %(name)s | %(message)s',
            'datefmt': '%Y-%m-%d %H:%M:%S',
        },
        'detailed': {
            'format': (
                '%(asctime)s | %(levelname)-8s | %(name)s | '
                '%(filename)s:%(lineno)d | %(message)s'
            ),
            'datefmt': '%Y-%m-%d %H:%M:%S',
        }
    },

    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'level': CONSOLE_LEVEL,
            'formatter': 'simple'
        },
        'file': {
            'class': 'logging.FileHandler',
            'level': FILE_LEVEL,
            'filename': LOG_FILE,
            'mode': 'a',
            'formatter': 'detailed',
        },
    },

    'loggers': {
        'VariantValidator': {
            'level': 'DEBUG',
            'handlers': ['console', 'file'],
            'propagate': False,
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
