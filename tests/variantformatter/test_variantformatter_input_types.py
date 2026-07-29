from VariantFormatter import simpleVariantFormatter
import VariantValidator
vfo = VariantValidator.Validator()


class TestVFvariantsInputs(object):
    @classmethod
    def setup_class(cls):
        vfo.testing = True

    def test_hybrid_syntax_1(self):
        results = simpleVariantFormatter.format("chr17:50198002C>A", 'GRCh38',
                                                                 'all', "all", True, False, testing=False)
        print(results)
        assert 'chr17:50198002C>A' in results.keys()
        assert 'NC_000017.11:g.50198002C>A' in results["chr17:50198002C>A"].keys()

    def test_hybrid_syntax_2(self):
        results = simpleVariantFormatter.format("17:50198002C>A", 'GRCh38',
                                                                 'all', "all", True, False, testing=False)
        print(results)
        assert '17:50198002C>A' in results.keys()
        assert 'NC_000017.11:g.50198002C>A' in results["17:50198002C>A"].keys()

# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
