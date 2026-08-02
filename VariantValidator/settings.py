"""
===============================================================================
VariantValidator global settings
===============================================================================

This module defines the global configuration used throughout
VariantValidator.

Configuration values may originate from one of four sources, listed below
in order of precedence.

Configuration precedence
------------------------

1. Runtime overrides
   Applications may override the configured logging levels for the current
   process using:

       VariantValidator.logger.configure_logging(
           console_level=...,
           file_level=...
       )

   These overrides affect only the current process and do not modify the
   configuration file or the default settings.

2. Environment variable overrides
   Selected settings may be overridden using environment variables,
   allowing deployment-specific configuration without modifying the source
   code or user configuration file.

3. VariantValidator configuration file
   User configuration is read from:

       ~/.variantvalidator

   unless an alternative configuration file is specified using:

       VARIANTVALIDATOR_TEST_CONFIG

4. Built-in defaults
   Default values defined in this module are used whenever no higher
   precedence configuration is supplied.

===============================================================================
Cache configuration
===============================================================================

VariantValidator uses several optional in-memory caches to improve
performance by reducing repeated database lookups and sequence retrieval
operations. Each cache may be enabled or disabled independently and its
maximum size configured.

DBGet lookup cache
------------------
Caches deterministic database lookup results returned by DBGet helper
methods.

Settings:

    vvDB_GET_CACHE
        Enable or disable the DBGet lookup cache.

    vvDB_GET_CACHE_SIZE
        Maximum number of cached lookup results.

Environment overrides:

    VV_DB_GET_CACHE
    VV_DB_GET_CACHE_SIZE


SeqFetcher cache
----------------
Caches recently retrieved sequence fragments to avoid repeated sequence
fetch operations.

Settings:

    SEQFETCHER_CACHE
        Enable or disable the SeqFetcher cache.

    SEQFETCHER_CACHE_SIZE
        Maximum number of cached sequence fragments.

Environment overrides:

    SEQFETCHER_CACHE
    SEQFETCHER_CACHE_SIZE


HGVS data-provider cache
------------------------
Caches HGVS data-provider (HDP) lookups performed through the local HGVS
interface.

Settings:

    vvHGVS_HDP_CACHE
        Enable or disable the local HGVS HDP cache.

    vvHGVS_HDP_CACHE_SIZE
        Maximum number of cached HGVS data-provider lookups.

Environment overrides:

    VV_HGVS_HDP_CACHE
    VV_HGVS_HDP_CACHE_SIZE

===============================================================================
Configuration file
===============================================================================

VariantValidator reads user configuration from a configuration file
located at:

    ~/.variantvalidator

This file typically contains:

    • Database connection settings
    • Logging configuration
    • Other user-specific configuration values

The configuration file location may be overridden by setting:

    VARIANTVALIDATOR_TEST_CONFIG

===============================================================================
Logging configuration
===============================================================================

Logging behaviour is configured using the [logging] section of the
VariantValidator configuration file.

Supported configuration options:

    log
        Enable or disable logging globally.

    console
        Console logging level.

    file
        Log file logging level.

    file_name
        Path to the log file.

If no log file is specified, VariantValidator writes to:

    ~/.vv_errorlog

Applications may override the configured console and/or file logging
levels for the current process by calling:

    VariantValidator.logger.configure_logging(
        console_level=...,
        file_level=...
    )

without modifying the configuration file or the default settings defined
in this module.

===============================================================================
"""

import os
from configparser import ConfigParser

# =============================================================================
# Cache configuration
#
# VariantValidator uses several optional in-memory caches to reduce repeated
# database lookups and sequence retrieval operations. These caches improve
# validation performance for repeated queries while allowing memory usage to be
# tuned for different deployment environments. Each cache can be enabled or
# disabled independently, and cache sizes may be overridden using environment
# variables.
# =============================================================================

# DBGet lookup cache settings.
vvDB_GET_CACHE = True
vvDB_GET_CACHE_SIZE = 20000

# SeqFetcher sequence cache settings.
SEQFETCHER_CACHE = True
SEQFETCHER_CACHE_SIZE = 32768

# vvHGVS HDP cache settings.
vvHGVS_HDP_CACHE = True
vvHGVS_HDP_CACHE_SIZE = 1000


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

if "SEQFETCHER_CACHE" in os.environ:
    SEQFETCHER_CACHE = os.environ["SEQFETCHER_CACHE"].lower() in (
        "true",
        "1",
        "yes",
    )

if "SEQFETCHER_CACHE_SIZE" in os.environ:
    SEQFETCHER_CACHE_SIZE = int(
        os.environ["SEQFETCHER_CACHE_SIZE"]
    )

if "vvHGVS_HDP_CACHE" in os.environ:
    vvHGVS_HDP_CACHE = os.environ["vvHGVS_HDP_CACHE"].lower() in (
        "true",
        "1",
        "yes",
    )

if "vvHGVS_HDP_CACHE_SIZE" in os.environ:
    vvHGVS_HDP_CACHE_SIZE = int(
        os.environ["vvHGVS_HDP_CACHE_SIZE"]
    )


# ----------------------------------------
# READ VV CONFIG FILE
# ----------------------------------------
config = ConfigParser()
def get_config_dir():
    if 'VARIANTVALIDATOR_TEST_CONFIG' in os.environ:
        return os.environ['VARIANTVALIDATOR_TEST_CONFIG']

    return os.path.join(
        os.path.expanduser('~'),
        '.variantvalidator'
    )

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
