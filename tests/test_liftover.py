from VariantFormatter import simpleVariantFormatter
from unittest import TestCase


class TestLiftover(TestCase):

    def test_gene_variant(self):
        result = simpleVariantFormatter.format('NC_000023.10:g.591732del', 'GRCh37', 'refseq', None, False, True)
        assert result["NC_000023.10:g.591732del"]["NC_000023.10:g.591732del"]["g_hgvs"] == "NC_000023.10:g.591732del"
        assert result["NC_000023.10:g.591732del"]["NC_000023.10:g.591732del"]["p_vcf"] == "X:591731:TA:T"
        assert result["NC_000023.10:g.591732del"]["NC_000023.10:g.591732del"]["hgvs_t_and_p"]["NM_000451.3"]["alt_genomic_loci"][0]["grch37"]["hgvs_genomic_description"] == "NC_000024.9:g.541732del"

    def test_mito_variant(self):
        result = simpleVariantFormatter.format('NC_012920.1:g.100del', 'GRCh37', 'refseq', None, False, True)
        assert result["NC_012920.1:g.100del"]["NC_012920.1:g.100del"]["g_hgvs"] == "NC_012920.1:m.101del"
        assert result["NC_012920.1:g.100del"]["NC_012920.1:g.100del"]["p_vcf"] == "M:99:TG:T"


