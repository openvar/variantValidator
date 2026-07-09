#!/usr/bin/env python3

"""
Interactive configuration utility for VariantValidator.
"""

from __future__ import annotations

import argparse
import configparser
from importlib import resources
from pathlib import Path

from VariantValidator import settings


def default_config_file():
    """
    Return the packaged default configuration file.
    """
    return (
        resources.files("VariantValidator")
        / "configuration"
        / "default.ini"
    )


def load_configuration() -> tuple[configparser.ConfigParser, str, bool]:
    """
    Load the VariantValidator configuration.

    Returns
    -------
    config
        Loaded ConfigParser object.

    config_path
        User configuration file.

    new_file
        True if the user configuration does not yet exist.
    """

    config_path = Path(settings.get_config_dir())

    config = configparser.ConfigParser()

    if config_path.exists():
        config.read(config_path)
        new_file = False
    else:
        with resources.as_file(default_config_file()) as default_file:
            config.read(default_file)
        new_file = True

    return config, config_path, new_file


def configure(
    config: configparser.ConfigParser,
    section: str | None = None,
) -> bool:
    """
    Interactively edit configuration values.

    Returns
    -------
    bool
        True if any values were modified.
    """

    values_changed = False

    for current_section in config.sections():

        if section is not None and current_section != section:
            continue

        print(f"Section: {current_section}")

        for name, value in config.items(current_section):

            response = input(
                f"{name} [{value}]: "
            ).strip()

            if response:
                config.set(
                    current_section,
                    name,
                    response,
                )
                values_changed = True

        print()

    return values_changed


def write_configuration(
    config: configparser.ConfigParser,
    filename,
) -> None:
    """
    Write the configuration to disk.
    """

    with open(filename, "w") as handle:
        config.write(handle)


def main() -> int:
    """
    VariantValidator configuration utility.
    """

    parser = argparse.ArgumentParser(
        description="Configure VariantValidator."
    )

    parser.add_argument(
        "-s",
        "--section",
        choices=[
            "mysql",
            "seqrepo",
            "postgres",
            "logging",
            "EntrezID",
        ],
        help="Configure a single section.",
    )

    args = parser.parse_args()

    config, config_file, new_file = load_configuration()

    changed = configure(
        config,
        section=args.section,
    )

    if new_file or changed:
        write_configuration(
            config,
            config_file,
        )
        print(f"Configuration written to {config_file}")

    else:
        print("No changes made.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


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