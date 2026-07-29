import logging.config

from VariantValidator import settings


def configure_logging(
    console_level=None,
    file_level=None,
):
    """
    Configure VariantValidator logging.

    Parameters
    ----------
    console_level : str, optional
        Override the console logging level.

    file_level : str, optional
        Override the file logging level.

    Notes
    -----
    If no overrides are supplied, the logging configuration from
    VariantValidator.settings is used unchanged.

    Applications such as the CLI may override the console or file
    logging level without modifying the default configuration.
    """

    #
    # STEP 1: Configure RotatingFileHandler
    #
    if file_level is None:
        f_level = settings.LOGGING_CONFIG["handlers"]["file"]["level"]
    else:
        f_level = file_level.upper()

    settings.LOGGING_CONFIG["handlers"]["file"] = {
        "class": "logging.handlers.RotatingFileHandler",
        "level": f_level,
        "filename": settings.LOGGING_CONFIG["handlers"]["file"]["filename"],
        "mode": "a",
        "maxBytes": 500000,
        "backupCount": 2,
        "formatter": "detailed",
    }

    #
    # STEP 2: Configure console handler
    #
    if console_level is None:
        c_level = settings.LOGGING_CONFIG["handlers"]["console"]["level"]
    else:
        c_level = console_level.upper()

    settings.LOGGING_CONFIG["handlers"]["console"]["level"] = c_level

    #
    # STEP 3: Apply logging configuration
    #
    logging.config.dictConfig(settings.LOGGING_CONFIG)


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later