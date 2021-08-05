import logging.config
from configparser import ConfigParser
from . import settings

# Change settings based on config
config = ConfigParser()
config.read(settings.CONFIG_DIR)

if config['logging'].getboolean('log') is True:
    settings.LOGGING_CONFIG['handlers']['console']['level'] = config['logging']['console'].upper()
    settings.LOGGING_CONFIG['handlers']['file']['level'] = config['logging']['file'].upper()

    logging.config.dictConfig(settings.LOGGING_CONFIG)
else:
    logging.getLogger('VariantValidator').addHandler(logging.NullHandler())

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
