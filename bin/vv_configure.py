#! /usr/bin/env python
from __future__ import print_function
import argparse
import os
import configparser
import pkgutil


def find_root():
    package = pkgutil.get_loader('VariantValidator')
    path = os.path.dirname(os.path.dirname(package.get_filename()))
    return path


def read_settings():
    root = find_root()
    settings_file = os.path.join(root, 'VariantValidator', 'settings.py')
    with open(settings_file) as f:
        values = {}
        exec(f.read(), {}, values)
        return values


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--section', choices=['mysql', 'seqrepo', 'postgres', 'logging', 'EntrezID', 'liftover'],
                        nargs='?', help='Optional choice of section to configure')

    args = parser.parse_args()

    settings = read_settings()
    newfile = False

    if os.path.exists(settings['CONFIG_DIR']):
        readfile = settings['CONFIG_DIR']
    else:
        root = find_root()
        readfile = os.path.join(root, 'configuration', 'default.ini')
        newfile = True

    config = configparser.ConfigParser()
    config.read(readfile)

    values_changed = False

    for section in config.sections():
        if not newfile and args.section and args.section != section:
            continue
        print('Section:', section)
        for name, value in config.items(section):
            print("{} [{}]: ".format(name, value), end="")
            newval = input()
            if newval != '':
                config.set(section, name, newval.strip())
                values_changed = True

        print()

    if newfile or values_changed:
        with open(settings['CONFIG_DIR'], 'w') as fh:
            config.write(fh)
