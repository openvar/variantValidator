import importlib
import unittest
import warnings
from unittest.mock import patch


class TestVersionFetching(unittest.TestCase):

    @patch("importlib.metadata.version")
    def test_version_fetching_package_not_found(self, mock_version):
        mock_version.side_effect = importlib.metadata.PackageNotFoundError

        with warnings.catch_warnings(record=True) as caught_warnings:
            import VariantValidator.version as version
            importlib.reload(version)

            self.assertTrue(
                any(
                    "can't get __version__" in str(warning.message)
                    for warning in caught_warnings
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


if __name__ == "__main__":
    unittest.main()


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
