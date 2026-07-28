from unittest import TestCase

import VariantValidator


vval = VariantValidator.Validator()


class TestLiftover(TestCase):

    def test_gene_variant_par_vv(self):
        variant = 'NC_000023.10:g.591732del'
        genome_build = 'GRCh38'
        select_transcripts = 'all'
        result = vval.validate(
            variant,
            genome_build,
            select_transcripts,
            transcript_set="refseq"
        ).format_as_dict(test=True)

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
        assert (
            "ParRegionWarning: Variant is located in a pseudoautosomal region (PAR) of the X and Y chromosomes, "
            "so the Y context description has been moved to alt_genomic_loci"
        ) in result["NM_000451.4:c.100del"]["validation_warnings"]

        variant = 'NC_000024.10:g.630997del'
        genome_build = 'GRCh38'
        select_transcripts = 'all'
        result = vval.validate(
            variant,
            genome_build,
            select_transcripts,
            transcript_set="refseq"
        ).format_as_dict(test=True)

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
        assert (
            "ParRegionWarning: Variant is located in a pseudoautosomal region (PAR) of the X and Y chromosomes, "
            "so the Y context description has been moved to alt_genomic_loci"
        ) in result["NM_000451.4:c.100del"]["validation_warnings"]

    def test_gene_variant_non_par_vv(self):
        variant = 'NC_000024.9:g.9197998C>T'
        genome_build = 'GRCh37'
        select_transcripts = 'all'
        result = vval.validate(
            variant,
            genome_build,
            select_transcripts,
            transcript_set="refseq"
        ).format_as_dict(test=True)

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

