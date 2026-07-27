from unittest import TestCase
from unittest.mock import MagicMock

from VariantValidator.modules.methyl_syntax import methyl_syntax


class TestMethylSyntax(TestCase):

    def test_no_methylation_syntax(self):
        variant = MagicMock()
        variant.quibble = "NM_000001.1:c.123A>G"

        result = methyl_syntax(variant)

        self.assertIsNone(result)
        self.assertEqual(
            variant.quibble,
            "NM_000001.1:c.123A>G",
        )

    def test_gom(self):
        variant = MagicMock()
        variant.quibble = "NM_000001.1:c.123|gom"

        result = methyl_syntax(variant)

        self.assertIs(result, variant)
        self.assertEqual(
            variant.reformat_output,
            "|gom",
        )
        self.assertEqual(
            variant.quibble,
            "NM_000001.1:c.123=",
        )

    def test_lom(self):
        variant = MagicMock()
        variant.quibble = "NM_000001.1:c.123|lom"

        result = methyl_syntax(variant)

        self.assertIs(result, variant)
        self.assertEqual(
            variant.reformat_output,
            "|lom",
        )
        self.assertEqual(
            variant.quibble,
            "NM_000001.1:c.123=",
        )

    def test_met(self):
        variant = MagicMock()
        variant.quibble = "NM_000001.1:c.123|met=0.5"

        result = methyl_syntax(variant)

        self.assertIs(result, variant)
        self.assertEqual(
            variant.reformat_output,
            "|met=",
        )
        self.assertEqual(
            variant.quibble,
            "NM_000001.1:c.123=",
        )

    def test_unrecognised_pipe_syntax(self):
        variant = MagicMock()
        variant.quibble = "NM_000001.1:c.123|foo"

        result = methyl_syntax(variant)

        self.assertIsNone(result)
        self.assertEqual(
            variant.quibble,
            "NM_000001.1:c.123|foo",
        )


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
