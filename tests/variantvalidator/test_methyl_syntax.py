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


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
