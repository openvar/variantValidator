import importlib
import unittest
from unittest.mock import patch


class TestVariantFormatterVersion(unittest.TestCase):

    @patch("importlib.metadata.version", return_value="3.2.1")
    def test_version_fetching_release_version(self, mock_version):
        import VariantFormatter
        importlib.reload(VariantFormatter)

        self.assertEqual(VariantFormatter.__version__, "3.2.1")
        self.assertTrue(VariantFormatter._is_released_version)
        mock_version.assert_called_with("VariantFormatter")

    @patch("importlib.metadata.version", return_value="3.2.1.dev1")
    def test_version_fetching_dev_version(self, mock_version):
        import VariantFormatter
        importlib.reload(VariantFormatter)

        self.assertEqual(VariantFormatter.__version__, "3.2.1.dev1")
        self.assertFalse(VariantFormatter._is_released_version)
        mock_version.assert_called_with("VariantFormatter")


if __name__ == "__main__":
    unittest.main()


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
