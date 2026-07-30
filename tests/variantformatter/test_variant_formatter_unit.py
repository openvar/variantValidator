"""
High-value branch coverage tests for VariantFormatter.variantformatter.
"""

import VariantValidator

from VariantFormatter.variantformatter import FormatVariant


class TestVariantFormatterBranches:

    @classmethod
    def setup_class(cls):
        cls.validator = VariantValidator.Validator()

    def _transcripts(self, result, variant):
        return result.stucture_data()[variant]["hgvs_t_and_p"]

    def test_select_transcripts_mane_select(self):
        variant = "NC_000017.11:g.50198003del"

        result = FormatVariant(
            variant,
            "GRCh38",
            self.validator,
            specify_transcripts="mane_select",
        )

        transcripts = self._transcripts(result, variant)

        assert transcripts

        assert "NM_000088.4" in transcripts
        assert "ENST00000225964.10" in transcripts

        assert transcripts["NM_000088.4"]["select_status"]["mane_select"]
        assert transcripts["ENST00000225964.10"]["select_status"]["mane_select"]

    def test_select_transcripts_select(self):
        variant = "NC_000017.11:g.50198003del"

        result = FormatVariant(
            variant,
            "GRCh38",
            self.validator,
            specify_transcripts="select",
        )

        transcripts = self._transcripts(result, variant)

        assert transcripts
        assert "NM_000088.4" in transcripts

    def test_select_transcripts_mane(self):
        variant = "NC_000017.11:g.50198003del"

        result = FormatVariant(
            variant,
            "GRCh38",
            self.validator,
            specify_transcripts="mane",
        )

        transcripts = self._transcripts(result, variant)

        assert transcripts
        assert "NM_000088.4" in transcripts

    def test_json_transcript_selection(self):
        variant = "NC_000017.11:g.50198003del"

        result = FormatVariant(
            variant,
            "GRCh38",
            self.validator,
            specify_transcripts='["NM_000088.4"]',
        )

        transcripts = self._transcripts(result, variant)

        assert list(transcripts) == ["NM_000088.4"]

    def test_variantvalidator_primary_loci_structure(self):
        variant = "NC_000017.11:g.50198003del"

        result = FormatVariant(
            variant,
            "GRCh38",
            self.validator,
            specify_transcripts="mane_select",
            legacy_genomic_structure=False,
        )

        transcript = next(
            iter(self._transcripts(result, variant).values())
        )

        primary = transcript["primary_assembly_loci"]

        assert isinstance(primary, dict)
        assert "grch38" in primary

        grch38 = primary["grch38"]

        assert "hgvs_genomic_description" in grch38
        assert "vcf" in grch38

    def test_liftover_variantvalidator_structure(self):
        variant = "NC_000017.11:g.50198003del"

        result = FormatVariant(
            variant,
            "GRCh38",
            self.validator,
            specify_transcripts="mane_select",
            legacy_genomic_structure=False,
            liftover=True,
        )

        transcript = next(
            iter(self._transcripts(result, variant).values())
        )

        assert "primary_assembly_loci" in transcript
        assert "alt_genomic_loci" in transcript

    def test_protein_descriptions_generated(self):
        variant = "NC_000017.11:g.50198003del"

        result = FormatVariant(
            variant,
            "GRCh38",
            self.validator,
            specify_transcripts='["NM_000088.4"]',
        )

        transcript = self._transcripts(result, variant)["NM_000088.4"]

        assert transcript["t_hgvs"] is not None
        assert transcript["p_hgvs_tlc"] is not None
        assert transcript["p_hgvs_slc"] is not None

    def test_check_only_returns_genomic_description(self):
        variant = "NC_000017.11:g.50198003del"

        result = FormatVariant(
            variant,
            "GRCh38",
            self.validator,
            checkOnly=True,
        )

        assert result.genomic_descriptions.g_hgvs is not None

    def test_legacy_and_new_structure_are_both_supported(self):
        variant = "NC_000017.11:g.50198003del"

        legacy = FormatVariant(
            variant,
            "GRCh38",
            self.validator,
            specify_transcripts='["NM_000088.4"]',
            legacy_genomic_structure=True,
        )

        modern = FormatVariant(
            variant,
            "GRCh38",
            self.validator,
            specify_transcripts='["NM_000088.4"]',
            legacy_genomic_structure=False,
        )

        legacy_tx = self._transcripts(legacy, variant)["NM_000088.4"]
        modern_tx = self._transcripts(modern, variant)["NM_000088.4"]

        assert "primary_assembly_loci" in legacy_tx
        assert "primary_assembly_loci" in modern_tx

    def test_refseq_only_json_filter(self):
        variant = "NC_000017.11:g.50198003del"

        result = FormatVariant(
            variant,
            "GRCh38",
            self.validator,
            specify_transcripts='["NM_000088.4"]',
        )

        transcripts = self._transcripts(result, variant)

        assert len(transcripts) == 1
        assert "NM_000088.4" in transcripts

    def test_ensembl_only_json_filter(self):
        variant = "NC_000017.11:g.50198003del"

        result = FormatVariant(
            variant,
            "GRCh38",
            self.validator,
            specify_transcripts='["ENST00000225964.10"]',
        )

        transcripts = self._transcripts(result, variant)

        assert len(transcripts) == 1
        assert "ENST00000225964.10" in transcripts

    def test_primary_loci_contains_expected_fields(self):
        variant = "NC_000017.11:g.50198003del"

        result = FormatVariant(
            variant,
            "GRCh38",
            self.validator,
            specify_transcripts='["NM_000088.4"]',
            legacy_genomic_structure=False,
        )

        transcript = self._transcripts(result, variant)["NM_000088.4"]

        grch38 = transcript["primary_assembly_loci"]["grch38"]

        assert "hgvs_genomic_description" in grch38
        assert "vcf" in grch38

        assert grch38["hgvs_genomic_description"] == variant

        vcf = grch38["vcf"]

        assert vcf["chr"] == "17"
        assert vcf["ref"] == "AC"
        assert vcf["alt"] == "A"

    def test_liftover_adds_alt_loci(self):
        variant = "NC_000017.11:g.50198003del"

        result = FormatVariant(
            variant,
            "GRCh38",
            self.validator,
            specify_transcripts='["NM_000088.4"]',
            legacy_genomic_structure=False,
            liftover=True,
        )

        transcript = self._transcripts(result, variant)["NM_000088.4"]

        assert "alt_genomic_loci" in transcript
        assert transcript["alt_genomic_loci"] is not None
