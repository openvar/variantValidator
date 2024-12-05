from VariantValidator import Validator
from unittest import TestCase


class TestAltAlignments(TestCase):

    @classmethod
    def setup_class(cls):
        cls.vv = Validator()
        cls.vv.testing = True

    def test_sub(self):
        variant = 'NW_011332691.1(NM_012234.6):c.335+1G>C'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_012234.6:c.335+1G>C" in list(results.keys())
        assert results["NM_012234.6:c.335+1G>C"][
                   "genome_context_intronic_sequence"] == "NW_011332691.1(NM_012234.6):c.335+1G>C"
        assert results["NM_012234.6:c.335+1G>C"]["alt_genomic_loci"][0]["grch38"][
            "hgvs_genomic_description"] == "NW_011332691.1:g.53828C>G"

    def test_del(self):
        variant = 'NW_011332691.1(NM_012234.6):c.335+1del'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_012234.6:c.335+1del" in list(results.keys())
        assert "NW_011332691.1(NM_012234.6):c.335+1del" in results[
            "NM_012234.6:c.335+1del"]["genome_context_intronic_sequence"]
        assert "NW_011332691.1:g.53828del" in results["NM_012234.6:c.335+1del"]["alt_genomic_loci"][0]["grch38"][
            "hgvs_genomic_description"]

    def test_delins(self):
        variant = "NW_011332691.1(NM_012234.6):c.335+1delinsAT"
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_012234.6:c.335+1delinsAT" in list(results.keys())
        assert "NW_011332691.1(NM_012234.6):c.335+1delinsAT" in results[
            "NM_012234.6:c.335+1delinsAT"]["genome_context_intronic_sequence"]
        assert "NW_011332691.1:g.53828delinsAT" in results["NM_012234.6:c.335+1delinsAT"]["alt_genomic_loci"][
            0]["grch38"]["hgvs_genomic_description"]

    def test_inv(self):
        variant = "NW_011332691.1(NM_012234.6):c.335+1_335+6inv"
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_012234.6:c.335+1_335+6inv" in list(results.keys())
        assert "NW_011332691.1(NM_012234.6):c.335+1_335+6inv" in results[
            "NM_012234.6:c.335+1_335+6inv"]["genome_context_intronic_sequence"]
        assert "NW_011332691.1:g.53823_53828inv" in results["NM_012234.6:c.335+1_335+6inv"]["alt_genomic_loci"][0][
            "grch38"]["hgvs_genomic_description"]

    def test_dup(self):
        variant = "NW_011332691.1(NM_012234.6):c.335+1dup"
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_012234.6:c.335+1dup" in list(results.keys())
        assert "NW_011332691.1(NM_012234.6):c.335+1dup" in results[
            "NM_012234.6:c.335+1dup"]["genome_context_intronic_sequence"]
        assert "NW_011332691.1:g.53828dup" in results["NM_012234.6:c.335+1dup"]["alt_genomic_loci"][0]["grch38"][
            "hgvs_genomic_description"]

    def test_identity(self):
        variant = "NW_011332691.1(NM_012234.6):c.335+1_335+6="
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_012234.6:c.335+1_335+6=" in list(results.keys())
        assert "NW_011332691.1(NM_012234.6):c.335+1_335+6=" in results[
            "NM_012234.6:c.335+1_335+6="]["genome_context_intronic_sequence"]
        assert "NW_011332691.1:g.53823_53828=" in results["NM_012234.6:c.335+1_335+6="]["alt_genomic_loci"][0]["grch38"][
            "hgvs_genomic_description"]

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



