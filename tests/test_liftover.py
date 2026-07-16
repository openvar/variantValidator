from VariantFormatter import simpleVariantFormatter
import VariantValidator
vval = VariantValidator.Validator()

from unittest import TestCase


class TestLiftover(TestCase):

    def test_gene_variant_par_vf(self):
        result = simpleVariantFormatter.format('NC_000023.10:g.591732del', 'GRCh37', 'refseq', "raw", False, True)
        assert result["NC_000023.10:g.591732del"]["NC_000023.10:g.591732del"]["g_hgvs"] == "NC_000023.10:g.591732del"
        assert result["NC_000023.10:g.591732del"]["NC_000023.10:g.591732del"]["p_vcf"] == "X:591731:TA:T"
        assert result["NC_000023.10:g.591732del"]["NC_000023.10:g.591732del"]["hgvs_t_and_p"]["NM_000451.3"]["alt_genomic_loci"][0]["grch37"]["hgvs_genomic_description"] == "NC_000024.9:g.541732del"
        assert result["NC_000023.10:g.591732del"]["NC_000023.10:g.591732del"]["genomic_variant_warnings"] == "ParRegionWarning: Variant is located in a pseudoautosomal region (PAR) of the X and Y chromosomes, so the Y context description has been moved to alt_genomic_loci"

        result = simpleVariantFormatter.format('NC_000024.9:g.541732del', 'GRCh37', 'refseq', "raw", False, True)
        assert result["NC_000024.9:g.541732del"]["NC_000024.9:g.541732del"]["g_hgvs"] == "NC_000024.9:g.541732del"
        assert result["NC_000024.9:g.541732del"]["NC_000024.9:g.541732del"]["hgvs_t_and_p"]["NM_000451.3"]["alt_genomic_loci"][0]["grch37"]["hgvs_genomic_description"] == "NC_000024.9:g.541732del"
        assert result["NC_000024.9:g.541732del"]["NC_000024.9:g.541732del"]["genomic_variant_warnings"] == "ParRegionWarning: Variant is located in a pseudoautosomal region (PAR) of the X and Y chromosomes, so the Y context description has been moved to alt_genomic_loci"


    def test_gene_variant_par_vf_non_par(self):
        result = simpleVariantFormatter.format('NC_000024.9:g.9197998C>T', 'GRCh37', 'refseq', 'mane_select', False, 'False')
        assert result["NC_000024.9:g.9197998C>T"]["NC_000024.9:g.9197998C>T"]["g_hgvs"] == "NC_000024.9:g.9197998C>T"

    def test_mito_variant_vf(self):
        result = simpleVariantFormatter.format('NC_012920.1:g.100del', 'GRCh37', 'refseq', None, False, True)
        assert result["NC_012920.1:g.100del"]["NC_012920.1:g.100del"]["g_hgvs"] == "NC_012920.1:m.101del"
        assert result["NC_012920.1:g.100del"]["NC_012920.1:g.100del"]["p_vcf"] == "M:99:TG:T"

    def test_gene_variant_par_vv(self):
        variant = 'NC_000023.10:g.591732del'
        genome_build = 'GRCh38'
        select_transcripts = 'all'
        result = vval.validate(variant, genome_build, select_transcripts, transcript_set="refseq").format_as_dict(test=True)
        assert result["NM_000451.4:c.100del"]["primary_assembly_loci"]["grch38"]["hgvs_genomic_description"] == "NC_000023.11:g.630997del"
        assert result["NM_000451.4:c.100del"]["primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000023.10:g.591732del"
        assert {
                "grch38": {
                    "hgvs_genomic_description": "NC_000024.10:g.630997del",
                    "vcf": {
                        "alt": "T",
                        "chr": "Y",
                        "pos": "630996",
                        "ref": "TA"
                    }
                }
            } in result["NM_000451.4:c.100del"]["alt_genomic_loci"]
        assert {
                "grch37": {
                    "hgvs_genomic_description": "NC_000024.9:g.541732del",
                    "vcf": {
                        "alt": "T",
                        "chr": "Y",
                        "pos": "541731",
                        "ref": "TA"
                    }
                }
            } in result["NM_000451.4:c.100del"]["alt_genomic_loci"]
        assert ("ParRegionWarning: Variant is located in a pseudoautosomal region (PAR) of the X and Y chromosomes, "
                "so the Y context description has been moved to alt_genomic_loci") in result[
            "NM_000451.4:c.100del"]["validation_warnings"]

        variant = 'NC_000024.10:g.630997del'
        genome_build = 'GRCh38'
        select_transcripts = 'all'
        result = vval.validate(variant, genome_build, select_transcripts, transcript_set="refseq").format_as_dict(test=True)
        assert result["NM_000451.4:c.100del"]["primary_assembly_loci"]["grch38"]["hgvs_genomic_description"] == "NC_000023.11:g.630997del"
        assert result["NM_000451.4:c.100del"]["primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000023.10:g.591732del"
        assert {
                "grch38": {
                    "hgvs_genomic_description": "NC_000024.10:g.630997del",
                    "vcf": {
                        "alt": "T",
                        "chr": "Y",
                        "pos": "630996",
                        "ref": "TA"
                    }
                }
            } in result["NM_000451.4:c.100del"]["alt_genomic_loci"]
        assert {
                "grch37": {
                    "hgvs_genomic_description": "NC_000024.9:g.541732del",
                    "vcf": {
                        "alt": "T",
                        "chr": "Y",
                        "pos": "541731",
                        "ref": "TA"
                    }
                }
            } in result["NM_000451.4:c.100del"]["alt_genomic_loci"]
        assert ("ParRegionWarning: Variant is located in a pseudoautosomal region (PAR) of the X and Y chromosomes, "
                "so the Y context description has been moved to alt_genomic_loci") in result[
            "NM_000451.4:c.100del"]["validation_warnings"]

    def test_gene_variant_non_par_vv(self):
        variant = 'NC_000024.9:g.9197998C>T'
        genome_build = 'GRCh37'
        select_transcripts = 'all'
        result = vval.validate(variant, genome_build, select_transcripts, transcript_set="refseq").format_as_dict(test=True)
        assert result["NM_001243721.2:c.911C>T"]["primary_assembly_loci"]["grch37"]["hgvs_genomic_description"] == "NC_000024.9:g.9197998C>T"


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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>
