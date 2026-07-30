"""
Tests for VariantFormatter.gapGenes.
"""

from VariantValidator import Validator
from VariantValidator.modules.transcript_map_data import TranscriptMapData
from VariantFormatter.gapGenes import (
    _get_alt_aln_method,
    fully_normalize,
)


class TestGapGenes:

    @classmethod
    def setup_class(cls):
        cls.validator = Validator()

        # Configure the validator for RefSeq mappings.
        cls.validator.alt_aln_method = "splign"
        cls.validator.create_additional_normalizers_and_mappers()

        cls.map_data = TranscriptMapData(cls.validator.hdp)

    def test_get_alt_aln_method_refseq(self):
        assert _get_alt_aln_method("NM_000088.4") == "splign"

    def test_get_alt_aln_method_ensembl(self):
        assert _get_alt_aln_method("ENST00000225964.10") == "genebuild"

    def test_fully_normalize_non_gap_gene(self):
        genomic = self.validator.hp.parse_hgvs_variant(
            "NC_000017.11:g.50198003del"
        )

        transcript = self.validator.vm.g_to_t(
            genomic,
            "NM_000088.4",
        )

        result = fully_normalize(
            hgvs_tx=transcript,
            hgvs_genomic=genomic,
            hn=self.validator.hn,
            reverse_normalizer=self.validator.reverse_hn,
            vm=self.validator.vm,
            vfo=self.validator,
            map_dat=self.map_data,
        )

        assert result is not None
        assert result.ac == transcript.ac
        assert result.type == "c"
        assert result.posedit is not None

    def test_fully_normalize_is_idempotent(self):
        genomic = self.validator.hp.parse_hgvs_variant(
            "NC_000017.11:g.50198003del"
        )

        transcript = self.validator.vm.g_to_t(
            genomic,
            "NM_000088.4",
        )

        first = fully_normalize(
            hgvs_tx=transcript,
            hgvs_genomic=genomic,
            hn=self.validator.hn,
            reverse_normalizer=self.validator.reverse_hn,
            vm=self.validator.vm,
            vfo=self.validator,
            map_dat=self.map_data,
        )

        second = fully_normalize(
            hgvs_tx=first,
            hgvs_genomic=genomic,
            hn=self.validator.hn,
            reverse_normalizer=self.validator.reverse_hn,
            vm=self.validator.vm,
            vfo=self.validator,
            map_dat=self.map_data,
        )

        assert str(first) == str(second)

    def test_fully_normalize_ensembl_transcript(self):
        """
        Exercise the Ensembl/genebuild mapping path.
        """

        genomic = self.validator.hp.parse_hgvs_variant(
            "NC_000017.11:g.50198003del"
        )

        transcript = self.validator.vm.g_to_t(
            genomic,
            "ENST00000225964.10",
            alt_aln_method="genebuild",
        )

        result = fully_normalize(
            hgvs_tx=transcript,
            hgvs_genomic=genomic,
            hn=self.validator.hn,
            reverse_normalizer=self.validator.reverse_hn,
            vm=self.validator.vm,
            vfo=self.validator,
            map_dat=self.map_data,
        )

        assert result is not None
        assert result.ac == transcript.ac
        assert result.type == "c"
        assert result.posedit is not None

    def test_fully_normalize_gap_gene(self):
        """
        Exercise a genuine gap-compensated mapping.
        """

        genomic = self.validator.hp.parse_hgvs_variant(
            "NC_000015.9:g.72105933del"
        )

        transcript = self.validator.vm.g_to_t(
            genomic,
            "NM_014249.3",
        )

        result = fully_normalize(
            hgvs_tx=transcript,
            hgvs_genomic=genomic,
            hn=self.validator.hn,
            reverse_normalizer=self.validator.reverse_hn,
            vm=self.validator.vm,
            vfo=self.validator,
            map_dat=self.map_data,
        )

        assert result is not None
        assert result.ac == transcript.ac
        assert result.type == "c"

    def test_fully_normalize_patch_alignment(self):
        """
        Exercise a transcript with alternate patch/NW mappings.
        """

        genomic = self.validator.hp.parse_hgvs_variant(
            "NC_000012.11:g.122064777C>A"
        )

        transcript = self.validator.vm.g_to_t(
            genomic,
            "NM_032790.3",
        )

        result = fully_normalize(
            hgvs_tx=transcript,
            hgvs_genomic=genomic,
            hn=self.validator.hn,
            reverse_normalizer=self.validator.reverse_hn,
            vm=self.validator.vm,
            vfo=self.validator,
            map_dat=self.map_data,
        )

        assert result is not None
        assert result.ac == transcript.ac
        assert result.type == "c"

# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later