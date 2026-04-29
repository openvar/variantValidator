import logging.config
from VariantValidator import settings


# --------------------------------------------------
# STEP 1: Use RotatingFileHandler
# --------------------------------------------------
settings.LOGGING_CONFIG['handlers']['file'] = {
    'class': 'logging.handlers.RotatingFileHandler',
    'level': settings.LOGGING_CONFIG['handlers']['file']['level'],
    'filename': settings.LOGGING_CONFIG['handlers']['file']['filename'],
    'mode': 'a',
    'maxBytes': 500000,
    'backupCount': 2,
    'formatter': 'detailed',
}


# --------------------------------------------------
# STEP 2: Apply logging configuration
# --------------------------------------------------
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>
