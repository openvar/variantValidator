import unittest
from VariantValidator.version import __version__ as variant_version, _is_released_version
from unittest.mock import patch, MagicMock
import importlib.metadata
import re
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
