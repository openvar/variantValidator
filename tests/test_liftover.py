from VariantFormatter import simpleVariantFormatter
from unittest import TestCase


class TestLiftover(TestCase):

    def test_gene_variant(self):
        result = simpleVariantFormatter.format('NC_000023.10:g.591732del', 'GRCh37', 'refseq', "raw", False, True)
        assert result["NC_000023.10:g.591732del"]["NC_000023.10:g.591732del"]["g_hgvs"] == "NC_000023.10:g.591732del"
        assert result["NC_000023.10:g.591732del"]["NC_000023.10:g.591732del"]["p_vcf"] == "X:591731:TA:T"
        assert result["NC_000023.10:g.591732del"]["NC_000023.10:g.591732del"]["hgvs_t_and_p"]["NM_000451.3"]["alt_genomic_loci"][0]["grch37"]["hgvs_genomic_description"] == "NC_000024.9:g.541732del"

    def test_mito_variant(self):
        result = simpleVariantFormatter.format('NC_012920.1:g.100del', 'GRCh37', 'refseq', None, False, True)
        assert result["NC_012920.1:g.100del"]["NC_012920.1:g.100del"]["g_hgvs"] == "NC_012920.1:m.101del"
        assert result["NC_012920.1:g.100del"]["NC_012920.1:g.100del"]["p_vcf"] == "M:99:TG:T"


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
