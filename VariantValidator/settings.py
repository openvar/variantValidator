import os

CONFIG_DIR = os.path.join(os.path.expanduser('~'), '.variantvalidator')

LOG_FILE = os.path.join(os.path.expanduser('~'), '.vv_errorlog')

LOGGING_CONFIG = {
    'version': 1,
    'formatters': {
        'simple': {
            'class': 'logging.Formatter',
            'format': '%(levelname)s: %(message)s'
        },
        'detailed': {
            'class': 'logging.Formatter',
            'format': '%(asctime)s %(name)-5s %(funcName)-10s (line %(lineno)d) %(levelname)-8s %(message)s'
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
            'level': 'ERROR',
            'handlers': ['console', 'file'],
            'propagate': 'no',
        }
    }
}

# <LICENSE>
# Copyright (C) 2016-2021 VariantValidator Contributors
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
