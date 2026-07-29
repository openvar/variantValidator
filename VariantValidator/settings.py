import os
from configparser import ConfigParser

config = ConfigParser()

vvDB_GET_CACHE = True
vvDB_GET_CACHE_SIZE = 20000

if "VV_DB_GET_CACHE" in os.environ:
    vvDB_GET_CACHE = os.environ["VV_DB_GET_CACHE"].lower() in (
        "true",
        "1",
        "yes",
    )

if "VV_DB_GET_CACHE_SIZE" in os.environ:
    vvDB_GET_CACHE_SIZE = int(
        os.environ["VV_DB_GET_CACHE_SIZE"]
    )

def get_config_dir():
    if 'VARIANTVALIDATOR_TEST_CONFIG' in os.environ:
        return os.environ['VARIANTVALIDATOR_TEST_CONFIG']

    return os.path.join(
        os.path.expanduser('~'),
        '.variantvalidator'
    )

# Read config
config.read(get_config_dir())


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
LOG_FILE = os.path.abspath(os.path.expanduser(LOG_FILE))

# ----------------------------------------
# ENSURE LOG DIRECTORY EXISTS
# ----------------------------------------
log_dir = os.path.dirname(LOG_FILE)
if log_dir:
    os.makedirs(log_dir, exist_ok=True)

# ----------------------------------------
# GLOBAL LOGGING SWITCH
# ----------------------------------------
logging_enabled = config.get('logging', 'log', fallback='true').lower() not in ('false', '0', 'no', 'off')

# ----------------------------------------
# LOG LEVELS
# ----------------------------------------
CONSOLE_LEVEL = config.get('logging', 'console', fallback='DEBUG').upper()
FILE_LEVEL = config.get('logging', 'file', fallback='ERROR').upper()

# If logging disabled → silence everything
if not logging_enabled:
    CONSOLE_LEVEL = 'CRITICAL'
    FILE_LEVEL = 'CRITICAL'

# ----------------------------------------
# LOGGING CONFIG
# ----------------------------------------
LOGGING_CONFIG = {
    'version': 1,

    'formatters': {
        'simple': {
            'format': (
                '%(asctime)s | %(levelname)-8s | %(name)s | '
                '%(filename)s:%(lineno)d | %(message)s'
            ),
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


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
