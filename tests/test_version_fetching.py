import unittest
import importlib
from unittest.mock import patch
import warnings


class TestVersionFetching(unittest.TestCase):

    @patch("importlib.metadata.version")
    def test_version_fetching_package_not_found(self, mock_version):
        mock_version.side_effect = importlib.metadata.PackageNotFoundError

        with warnings.catch_warnings(record=True) as w:
            import VariantValidator.version as version
            importlib.reload(version)

            self.assertTrue(
                any(
                    "can't get __version__" in str(warn.message)
                    for warn in w
                )
            )

            self.assertIsNone(version.__version__)
            self.assertFalse(version._is_released_version)

    @patch("importlib.metadata.version", return_value="3.2.1")
    def test_version_fetching_release_version(self, mock_version):
        import VariantValidator.version as version
        importlib.reload(version)

        self.assertEqual(version.__version__, "3.2.1")
        self.assertTrue(version._is_released_version)

    @patch("importlib.metadata.version", return_value="3.2.1.dev1")
    def test_version_fetching_dev_version(self, mock_version):
        import VariantValidator.version as version
        importlib.reload(version)

        self.assertEqual(version.__version__, "3.2.1.dev1")
        self.assertFalse(version._is_released_version)

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
