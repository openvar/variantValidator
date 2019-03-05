import os
import shutil
import configparser

CONFIG_DIR = os.path.join(os.path.expanduser('~'), '.variantvalidator')


def read_configuration():
    print("Going to read configuration here")
    config = configparser.ConfigParser()
    config.read(CONFIG_DIR)

    if config['mysql']['user'] == 'USERNAME' or config['mysql']['password'] == 'PASSWORD':
        print("MySQL username and password have not been updated from default.")
        exit_with_message()

    if config['postgres']['user'] == 'USERNAME' or config['postgres']['password'] == 'PASSWORD':
        print("PostgreSQL username and password have not been updated from default.")
        exit_with_message()

    if config['seqrepo']['location'] == 'PATH/TO/SEQREPO':
        print("Seqrepo directory location has not been updated from default.")
        exit_with_message()


def exit_with_message():
    print("Please edit your configuration file %s" % CONFIG_DIR)
    print()
    raise SystemExit


if os.path.exists(CONFIG_DIR):
    print("Configuration already set up for this user")
    read_configuration()
else:
    print("*-----------------------------*")
    print("| Welcome to VariantValidator |")
    print("*-----------------------------*")
    shutil.copyfile(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'configuration', 'default.ini'), CONFIG_DIR)
    print("A configuration file has been copied into your home directory (%s)." % CONFIG_DIR)
    print("Please edit this file with your database connection settings prior to continuing.")
    print("Items that must be changed are highlighted in capitals.")
    print()
    raise SystemExit

