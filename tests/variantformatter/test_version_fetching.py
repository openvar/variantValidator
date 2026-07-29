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


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
