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