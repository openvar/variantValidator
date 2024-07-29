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
