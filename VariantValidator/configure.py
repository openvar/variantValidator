import os
import shutil
import configparser
from .settings import CONFIG_DIR


def read_configuration():
    config = configparser.ConfigParser()
    config.read(CONFIG_DIR)

    if config['mysql']['user'] == 'USERNAME' or config['mysql']['password'] == 'PASSWORD':
        print("MySQL username and password have not been updated from default.")
        exit_with_message()

    if config['postgres']['user'] == 'USERNAME' or config['postgres']['password'] == 'PASSWORD':
        print("PostgreSQL username and password have not been updated from default.")
        exit_with_message()

    if config['seqrepo']['location'] == '/PATH/TO/SEQREPO':
        print("Seqrepo directory location has not been updated from default.")
        exit_with_message()


def exit_with_message():
    print("Please edit your configuration file %s" % CONFIG_DIR)
    print()
    raise SystemExit


if os.path.exists(CONFIG_DIR):
    read_configuration()
else:
    print("*-----------------------------*")
    print("| Welcome to VariantValidator |")
    print("*-----------------------------*")
    shutil.copyfile(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'configuration',
                                 'default.ini'), CONFIG_DIR)
    print("A configuration file has been copied into your home directory (%s)." % CONFIG_DIR)
    print("Please edit this file with your database connection settings prior to continuing.")
    print("Items that must be changed are highlighted in capitals.")
    print()
    raise SystemExit

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
