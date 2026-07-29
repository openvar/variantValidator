import os
import shutil
import configparser
from VariantValidator import settings


def read_configuration():
    config = configparser.ConfigParser()
    config.read(settings.get_config_dir())

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
    print("Please edit your configuration file %s" % settings.get_config_dir())
    print()
    raise SystemExit


if os.path.exists(settings.get_config_dir()):
    read_configuration()
else:
    print("*-----------------------------*")
    print("| Welcome to VariantValidator |")
    print("*-----------------------------*")
    shutil.copyfile(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'configuration',
                                 'default.ini'), settings.get_config_dir())
    print("A configuration file has been copied into your home directory (%s)." % settings.get_config_dir())
    print("Please edit this file with your database connection settings prior to continuing.")
    print("Items that must be changed are highlighted in capitals.")
    print()
    raise SystemExit

# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
