import unittest
from VariantValidator.version import _is_released_version
from unittest.mock import patch
import importlib.metadata
import warnings


class TestVersionFetching(unittest.TestCase):

    @patch('importlib.metadata.version')
    def test_version_fetching_package_not_found(self, mock_version):
        # Set up the mock to raise PackageNotFoundError
        mock_version.side_effect = importlib.metadata.PackageNotFoundError

        # Capture the warning using the warnings module
        with warnings.catch_warnings(record=True) as w:
            # Run the code from VariantValidator.version.py
            exec(open("VariantValidator/version.py").read(), globals())

            # Check that the warning was issued
            self.assertTrue(any(issubclass(warn.category, Warning) and "can't get __version__" in str(warn.message)
                                for warn in w))

        # Check that _is_released_version is False
        self.assertFalse(_is_released_version)


if __name__ == '__main__':
    unittest.main()

# <LICENSE>
# Copyright (C) 2016-2024 VariantValidator Contributors
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