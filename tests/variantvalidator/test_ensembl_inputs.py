from unittest import TestCase

from VariantValidator import Validator


class TestVariantsEnsembl(TestCase):
    @classmethod
    def setup_class(cls):
        cls.vv = Validator()

    def _validate(self, variant, build, transcript="all"):
        return self.vv.validate(
            variant, build, transcript, transcript_set="ensembl"
        ).format_as_dict(test=True)

    @staticmethod
    def _check(results, key, submitted, gene, protein, *, gene_ids=None,
               genome_context="", loci=None, refs=None):
        assert results["flag"] == "gene_variant"
        assert key in results
        result = results[key]

        assert result["submitted_variant"] == submitted
        assert result["gene_symbol"] == gene
        assert result["hgvs_transcript_variant"] == key
        assert result["genome_context_intronic_sequence"] == genome_context

        if gene_ids is not None:
            assert result["gene_ids"] == gene_ids
        if protein is not None:
            assert result["hgvs_predicted_protein_consequence"] == protein
        if loci is not None:
            assert result["alt_genomic_loci"] == []

            for assembly, expected in loci.items():
                assert result["primary_assembly_loci"][assembly] == expected

        if refs is not None:
            assert result["reference_sequence_records"] == refs

    def test_variant1(self):
        variant = "ENST00000225964.10:c.589-1GG>G"
        key = "ENST00000225964.10:c.590del"
        results = self._validate(variant, "GRCh38")

        self._check(
            results, key, variant, "COL1A1",
            {
                "tlr": "ENSP00000225964.6:p.(Gly197ValfsTer68)",
                "slr": "ENSP00000225964.6:p.(G197Vfs*68)",
            },
            gene_ids={
                "hgnc_id": "HGNC:2197", "entrez_gene_id": "1277",
                "ucsc_id": "uc002iqm.4", "omim_id": ["120150"],
            },
            loci={
                "hg38": {
                    "hgvs_genomic_description": "NC_000017.11:g.50198003del",
                    "vcf": {"chr": "chr17", "pos": "50198000", "ref": "AC", "alt": "A"},
                },
                "grch38": {
                    "hgvs_genomic_description": "NC_000017.11:g.50198003del",
                    "vcf": {"chr": "17", "pos": "50198000", "ref": "AC", "alt": "A"},
                },
            },
            refs={
                "transcript": "https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000225964.10",
                "protein": "https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000225964.6",
            },
        )

    def test_variant2(self):
        variant = "ENST00000371817.8:c.5071A>T"
        results = self._validate(variant, "GRCh38")

        self._check(
            results, variant, variant, "COL5A1",
            {
                "tlr": "ENSP00000360882.3:p.(Arg1691Ter)",
                "slr": "ENSP00000360882.3:p.(R1691*)",
            },
            gene_ids={
                "hgnc_id": "HGNC:2209", "entrez_gene_id": "1289",
                "ucsc_id": "uc004cfe.5", "omim_id": ["120215"],
            },
            loci={
                "hg38": {
                    "hgvs_genomic_description": "NC_000009.12:g.134829979A>T",
                    "vcf": {"chr": "chr9", "pos": "134829979", "ref": "A", "alt": "T"},
                },
                "grch38": {
                    "hgvs_genomic_description": "NC_000009.12:g.134829979A>T",
                    "vcf": {"chr": "9", "pos": "134829979", "ref": "A", "alt": "T"},
                },
            },
            refs={
                "transcript": "https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000371817.8",
                "protein": "https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000360882.3",
            },
        )

    def test_variant3(self):
        variant = "ENST00000269305.9:c.652_654del"
        results = self._validate(variant, "GRCh38")

        self._check(
            results, variant, variant, "TP53",
            {
                "tlr": "ENSP00000269305.4:p.(Val218del)",
                "slr": "ENSP00000269305.4:p.(V218del)",
            },
            gene_ids={
                "hgnc_id": "HGNC:11998", "entrez_gene_id": "7157",
                "ucsc_id": "uc060aur.1", "omim_id": ["191170"],
            },
            loci={
                "hg38": {
                    "hgvs_genomic_description": "NC_000017.11:g.7674883_7674885del",
                    "vcf": {"chr": "chr17", "pos": "7674876", "ref": "GCAC", "alt": "G"},
                },
                "grch38": {
                    "hgvs_genomic_description": "NC_000017.11:g.7674883_7674885del",
                    "vcf": {"chr": "17", "pos": "7674876", "ref": "GCAC", "alt": "G"},
                },
            },
            refs={
                "transcript": "https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000269305.9",
                "protein": "https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000269305.4",
            },
        )

    def test_variant4(self):
        variant = "ENST00000296388.10:c.2055+18G>A"
        results = self._validate(variant, "GRCh38")

        self._check(
            results, variant, variant, "P3H1",
            {"tlr": "ENSP00000296388.5:p.?", "slr": "ENSP00000296388.5:p.?"},
            gene_ids={
                "hgnc_id": "HGNC:19316", "entrez_gene_id": "64175",
                "ucsc_id": "", "omim_id": ["610339"],
            },
            genome_context="NC_000001.11(ENST00000296388.10):c.2055+18G>A",
            loci={
                "hg38": {
                    "hgvs_genomic_description": "NC_000001.11:g.42747254C>T",
                    "vcf": {"chr": "chr1", "pos": "42747254", "ref": "C", "alt": "T"},
                },
                "grch38": {
                    "hgvs_genomic_description": "NC_000001.11:g.42747254C>T",
                    "vcf": {"chr": "1", "pos": "42747254", "ref": "C", "alt": "T"},
                },
            },
            refs={
                "transcript": "https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000296388.10",
                "protein": "https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000296388.5",
            },
        )

    def test_variant5(self):
        variant = "ENST00000357654.9:c.301+1G>C"
        results = self._validate(variant, "GRCh38")

        self._check(
            results, variant, variant, "BRCA1",
            {"tlr": "ENSP00000350283.3:p.?", "slr": "ENSP00000350283.3:p.?"},
            gene_ids={
                "hgnc_id": "HGNC:1100", "entrez_gene_id": "672",
                "ucsc_id": "uc002ict.4", "omim_id": ["113705"],
            },
            genome_context="NC_000017.11(ENST00000357654.9):c.301+1G>C",
            loci={
                "hg38": {
                    "hgvs_genomic_description": "NC_000017.11:g.43104867C>G",
                    "vcf": {"chr": "chr17", "pos": "43104867", "ref": "C", "alt": "G"},
                },
                "grch38": {
                    "hgvs_genomic_description": "NC_000017.11:g.43104867C>G",
                    "vcf": {"chr": "17", "pos": "43104867", "ref": "C", "alt": "G"},
                },
            },
            refs={
                "transcript": "https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000357654.9",
                "protein": "https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000350283.3",
            },
        )

    def test_variant6(self):
        variant = "NC_000013.10:g.32929387T>C"
        results = self._validate(variant, "GRCh37")

        expected = {
            "ENST00000380152.3:c.7397T>C": "ENSP00000369497.3",
            "ENST00000544455.1:c.7397T>C": "ENSP00000439902.1",
        }

        for key, protein_id in expected.items():
            self._check(
                results, key, variant, "BRCA2",
                {
                    "tlr": f"{protein_id}:p.(Val2466Ala)",
                    "slr": f"{protein_id}:p.(V2466A)",
                },
                gene_ids={
                    "hgnc_id": "HGNC:1101", "entrez_gene_id": "675",
                    "ucsc_id": "uc001uub.2", "omim_id": ["600185"],
                },
                loci={
                    "hg19": {
                        "hgvs_genomic_description": "NC_000013.10:g.32929387T>C",
                        "vcf": {"chr": "chr13", "pos": "32929387", "ref": "T", "alt": "C"},
                    },
                    "grch37": {
                        "hgvs_genomic_description": "NC_000013.10:g.32929387T>C",
                        "vcf": {"chr": "13", "pos": "32929387", "ref": "T", "alt": "C"},
                    },
                },
            )

    def test_variant7(self):
        variant = "11-5248232-T-A"
        key = "ENST00000330597.3:c.*127A>T"
        results = self._validate(variant, "GRCh38")

        self._check(
            results, key, variant, "HBG1",
            {"tlr": "ENSP00000327431.3:p.(=)", "slr": "ENSP00000327431.3:p.(=)"},
            gene_ids={
                "hgnc_id": "HGNC:4831", "entrez_gene_id": "3047",
                "ucsc_id": "uc001mah.2", "omim_id": ["142200"],
            },
            loci={
                "hg38": {
                    "hgvs_genomic_description": "NC_000011.10:g.5248232T>A",
                    "vcf": {"chr": "chr11", "pos": "5248232", "ref": "T", "alt": "A"},
                },
                "grch38": {
                    "hgvs_genomic_description": "NC_000011.10:g.5248232T>A",
                    "vcf": {"chr": "11", "pos": "5248232", "ref": "T", "alt": "A"},
                },
            },
            refs={
                "transcript": "https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000330597.3",
                "protein": "https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000327431.3",
            },
        )

    def test_variant9(self):
        variant = "6-32012992-CG-C"
        results = self._validate(variant, "GRCh37")

        expected = {
            "ENST00000375247.2:c.10711del": (
                "ENSP00000364396.2",
                "ENSP00000364396.2:p.(Arg3571AlafsTer91)",
                "ENSP00000364396.2:p.(R3571Afs*91)",
                "https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000375247.2",
            ),
            "ENST00000451343.1:c.4del": (
                "ENSP00000407685.1",
                "ENSP00000407685.1:p.(Arg2AlafsTer91)",
                "ENSP00000407685.1:p.(R2Afs*91)",
                "https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000451343.1",
            ),
        }

        for key, (_, tlr, slr, transcript_url) in expected.items():
            self._check(
                results, key, variant, "TNXB",
                {"tlr": tlr, "slr": slr},
                gene_ids={
                    "hgnc_id": "HGNC:11976", "entrez_gene_id": "7148",
                    "ucsc_id": "uc063nnw.1", "omim_id": ["600985"],
                },
                loci={
                    "hg19": {
                        "hgvs_genomic_description": "NC_000006.11:g.32012993del",
                        "vcf": {"chr": "chr6", "pos": "32012992", "ref": "CG", "alt": "C"},
                    },
                    "grch37": {
                        "hgvs_genomic_description": "NC_000006.11:g.32012993del",
                        "vcf": {"chr": "6", "pos": "32012992", "ref": "CG", "alt": "C"},
                    },
                },
                refs={
                    "transcript": transcript_url,
                    "protein": f"https://grch37.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p={tlr.split(':')[0]}",
                },
            )

    def test_variant10(self):
        variant = "ENST00000298552.9:c.363+1dupG"
        key = "ENST00000298552.9:c.363+1dup"
        results = self._validate(variant, "GRCh38")

        self._check(
            results, key, variant, "TSC1",
            {"tlr": "ENSP00000298552.3:p.?", "slr": "ENSP00000298552.3:p.?"},
            gene_ids={
                "hgnc_id": "HGNC:12362", "entrez_gene_id": "7248",
                "ucsc_id": "uc004cca.3", "omim_id": ["605284"],
            },
            genome_context="NC_000009.12(ENST00000298552.9):c.363+1dup",
            loci={
                "hg38": {
                    "hgvs_genomic_description": "NC_000009.12:g.132925587dup",
                    "vcf": {"chr": "chr9", "pos": "132925585", "ref": "A", "alt": "AC"},
                },
                "grch38": {
                    "hgvs_genomic_description": "NC_000009.12:g.132925587dup",
                    "vcf": {"chr": "9", "pos": "132925585", "ref": "A", "alt": "AC"},
                },
            },
            refs={
                "transcript": "https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000298552.9",
                "protein": "https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000298552.3",
            },
        )

    def test_variant11(self):
        variant = "NC_000016.10:g.2099572TG>T"
        key = "ENST00000262304.9:c.10050+71del"
        results = self._validate(variant, "GRCh38")

        self._check(
            results, key, variant, "PKD1",
            {"tlr": "ENSP00000262304.4:p.?", "slr": "ENSP00000262304.4:p.?"},
            genome_context="NC_000016.10(ENST00000262304.9):c.10050+71del",
            loci={
                "hg38": {
                    "hgvs_genomic_description": "NC_000016.10:g.2099574del",
                    "vcf": {"chr": "chr16", "pos": "2099572", "ref": "TG", "alt": "T"},
                },
            },
        )

    def _check_gapped(self, results, key, variant, gene, protein):
        self._check(
            results, key, variant, gene, protein,
            gene_ids={
                "hgnc_id": "HGNC:6717", "entrez_gene_id": "8425",
                "ucsc_id": "uc032hxp.2", "omim_id": ["604710"],
            },
            loci={
                "hg19": {
                    "hgvs_genomic_description": "NC_000019.9:g.41123095dup",
                    "vcf": {"chr": "chr19", "pos": "41123093", "ref": "A", "alt": "AG"},
                },
                "grch37": {
                    "hgvs_genomic_description": "NC_000019.9:g.41123095dup",
                    "vcf": {"chr": "19", "pos": "41123093", "ref": "A", "alt": "AG"},
                },
            },
        )
        result = results[key]
        assert result["refseqgene_context_intronic_sequence"] == ""
        assert result["hgvs_lrg_transcript_variant"] == ""
        assert result["hgvs_lrg_variant"] == ""

    def test_variant12(self):
        variant = "19-41123094-G-GG"
        key = "ENST00000396819.3:c.3033_3034insGGT"
        results = self._validate(variant, "GRCh37")
        self._check_gapped(
            results, key, variant, "LTBP4",
            {
                "slr": "ENSP00000380031.3:p.(Q1011_Y1012insG)",
                "tlr": "ENSP00000380031.3:p.(Gln1011_Tyr1012insGly)",
            },
        )

    def test_variant12b(self):
        variant = "ENST00000396819.3:c.3033_3034insGGT"
        results = self._validate(variant, "GRCh37")
        self._check_gapped(
            results, variant, variant, "LTBP4",
            {
                "slr": "ENSP00000380031.3:p.(Q1011_Y1012insG)",
                "tlr": "ENSP00000380031.3:p.(Gln1011_Tyr1012insGly)",
            },
        )

    def _check_nr2e3(self, results, key, submitted):
        self._check(
            results, key, submitted, "NR2E3",
            {"tlr": "", "slr": ""},
            loci={
                "hg19": {
                    "hgvs_genomic_description": "NC_000015.9:g.72105933del",
                    "vcf": {"alt": "A", "chr": "chr15", "pos": "72105928", "ref": "AC"},
                },
                "grch37": {
                    "hgvs_genomic_description": "NC_000015.9:g.72105933del",
                    "vcf": {"alt": "A", "chr": "15", "pos": "72105928", "ref": "AC"},
                },
            },
        )

    def test_variant13(self):
        variant = "15-72105928-AC-A"
        key = "ENST00000398840.2:n.1133_1141="
        results = self._validate(variant, "GRCh37", "ENST00000398840.2")
        self._check_nr2e3(results, key, variant)

    def test_variant13b(self):
        variant = "ENST00000398840.2:n.1133_1141="
        results = self._validate(variant, "GRCh37", "ENST00000398840.2")
        self._check_nr2e3(results, variant, variant)

    def test_variant14(self):
        variant = "NC_000002.11:g.95847041_95847043GCG="
        key = "ENST00000340539.5:c.468_470="
        results = self._validate(variant, "GRCh37", "ENST00000340539.5")

        self._check(
            results, key, variant, "ZNF2",
            {
                "slr": "ENSP00000345392.5:p.(L156_R157=)",
                "tlr": "ENSP00000345392.5:p.(Leu156_Arg157=)",
            },
        )

        for assembly in ("hg19", "grch37"):
            assert (
                results[key]["primary_assembly_loci"][assembly]
                ["hgvs_genomic_description"]
                == "NC_000002.11:g.95847041_95847043="
            )


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 (or at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
