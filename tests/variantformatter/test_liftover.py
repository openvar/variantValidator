from unittest import TestCase

from VariantFormatter import simpleVariantFormatter


class TestLiftover(TestCase):

    def test_gene_variant_par_vf(self):
        result = simpleVariantFormatter.format(
            'NC_000023.10:g.591732del',
            'GRCh37',
            'refseq',
            "raw",
            False,
            True
        )
        assert result["NC_000023.10:g.591732del"]["NC_000023.10:g.591732del"]["g_hgvs"] == "NC_000023.10:g.591732del"
        assert result["NC_000023.10:g.591732del"]["NC_000023.10:g.591732del"]["p_vcf"] == "X:591731:TA:T"
        assert result["NC_000023.10:g.591732del"]["NC_000023.10:g.591732del"]["hgvs_t_and_p"]["NM_000451.3"]["alt_genomic_loci"][0]["grch37"]["hgvs_genomic_description"] == "NC_000024.9:g.541732del"
        assert result["NC_000023.10:g.591732del"]["NC_000023.10:g.591732del"]["genomic_variant_warnings"] == "ParRegionWarning: Variant is located in a pseudoautosomal region (PAR) of the X and Y chromosomes, so the Y context description has been moved to alt_genomic_loci"

        result = simpleVariantFormatter.format(
            'NC_000024.9:g.541732del',
            'GRCh37',
            'refseq',
            "raw",
            False,
            True
        )
        assert result["NC_000024.9:g.541732del"]["NC_000024.9:g.541732del"]["g_hgvs"] == "NC_000024.9:g.541732del"
        assert result["NC_000024.9:g.541732del"]["NC_000024.9:g.541732del"]["hgvs_t_and_p"]["NM_000451.3"]["alt_genomic_loci"][0]["grch37"]["hgvs_genomic_description"] == "NC_000024.9:g.541732del"
        assert result["NC_000024.9:g.541732del"]["NC_000024.9:g.541732del"]["genomic_variant_warnings"] == "ParRegionWarning: Variant is located in a pseudoautosomal region (PAR) of the X and Y chromosomes, so the Y context description has been moved to alt_genomic_loci"

    def test_gene_variant_par_vf_non_par(self):
        result = simpleVariantFormatter.format(
            'NC_000024.9:g.9197998C>T',
            'GRCh37',
            'refseq',
            'mane_select',
            False,
            'False'
        )
        assert result["NC_000024.9:g.9197998C>T"]["NC_000024.9:g.9197998C>T"]["g_hgvs"] == "NC_000024.9:g.9197998C>T"

    def test_mito_variant_vf(self):
        result = simpleVariantFormatter.format(
            'NC_012920.1:g.100del',
            'GRCh37',
            'refseq',
            None,
            False,
            True
        )
        assert result["NC_012920.1:g.100del"]["NC_012920.1:g.100del"]["g_hgvs"] == "NC_012920.1:m.101del"
        assert result["NC_012920.1:g.100del"]["NC_012920.1:g.100del"]["p_vcf"] == "M:99:TG:T"


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
