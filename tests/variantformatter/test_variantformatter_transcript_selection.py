from VariantFormatter import simpleVariantFormatter
import VariantValidator
import pytest


class TestVFvariantsTranscriptSelection(object):
    @classmethod
    def setup_class(cls):
        cls.vfo = VariantValidator.Validator()
        cls.vfo.testing = False

    def test_transcript_selection_raw(self):
        results = simpleVariantFormatter.format(
                'NC_000005.10:g.140114829del',
                'GRCh38', 'refseq', "raw",
                False, True, testing=False,
                validator=self.vfo)
        print(results)
        assert 'NC_000005.10:g.140114829del' in results.keys()
        assert 'NM_005859.3' in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()
        assert 'NM_005859.4' in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()
        assert 'NM_005859.5' in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()

    def test_transcript_selection_all(self):
        results = simpleVariantFormatter.format(
                'NC_000005.10:g.140114829del',
                'GRCh38', 'refseq', "all",
                False, True, testing=False,
                validator=self.vfo)
        print(results)
        assert 'NC_000005.10:g.140114829del' in results.keys()
        assert 'NM_005859.3' not in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()
        assert 'NM_005859.4' not in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()
        assert 'NM_005859.5' in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()

    def test_transcript_selection_mane_select(self):
        results = simpleVariantFormatter.format(
                'NC_000005.10:g.140114829del',
                'GRCh38', 'refseq', "mane_select",
                False, True, testing=False,
                validator=self.vfo)
        print(results)
        assert 'NC_000005.10:g.140114829del' in results.keys()
        assert 'NM_005859.3' not in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()
        assert 'NM_005859.4' not in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()
        assert 'NM_005859.5' in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()

    def test_transcript_selection_select(self):
        results = simpleVariantFormatter.format(
                'NC_000005.10:g.140114829del',
                'GRCh38', 'refseq', "select",
                False, True, testing=False,
                validator=self.vfo)
        print(results)
        assert 'NC_000005.10:g.140114829del' in results.keys()
        assert 'NM_005859.3' not in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()
        assert 'NM_005859.4' in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()
        assert 'NM_005859.5' in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()

    def test_transcript_selection_nm(self):
        results = simpleVariantFormatter.format(
                'NC_000005.10:g.140114829del',
                'GRCh38', 'refseq', "NM_005859.4",
                False, True, testing=False,
                validator=self.vfo)
        print(results)
        assert 'NC_000005.10:g.140114829del' in results.keys()
        assert 'NM_005859.3' not in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()
        assert 'NM_005859.4' in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()
        assert 'NM_005859.5' not in results['NC_000005.10:g.140114829del'][
            'NC_000005.10:g.140114829del']['hgvs_t_and_p'].keys()

    def test_transcript_selection_mane(self):
        results = simpleVariantFormatter.format(
                'NC_000007.14:g.140924703T>C',
                'GRCh38', 'refseq', "mane",
                False, True, testing=False,
                validator=self.vfo)
        print(results)
        assert 'NC_000007.14:g.140924703T>C' in results.keys()
        assert 'NM_004333.6' in results['NC_000007.14:g.140924703T>C'][
            'NC_000007.14:g.140924703T>C']['hgvs_t_and_p'].keys()
        assert 'NM_001374258.1' in results['NC_000007.14:g.140924703T>C'][
            'NC_000007.14:g.140924703T>C']['hgvs_t_and_p'].keys()
        assert 'NM_001354609.1' not in results['NC_000007.14:g.140924703T>C'][
            'NC_000007.14:g.140924703T>C']['hgvs_t_and_p'].keys()
        assert 'NM_001354609.2' not in results['NC_000007.14:g.140924703T>C'][
            'NC_000007.14:g.140924703T>C']['hgvs_t_and_p'].keys()
        assert 'NM_001374244.1' not in results['NC_000007.14:g.140924703T>C'][
            'NC_000007.14:g.140924703T>C']['hgvs_t_and_p'].keys()
        assert 'NM_001378467.1' not in results['NC_000007.14:g.140924703T>C'][
            'NC_000007.14:g.140924703T>C']['hgvs_t_and_p'].keys()

    def test_format_keyword_arguments(self):
        formatter = simpleVariantFormatter.SimpleVariantFormatter()
        formatter.validator = self.vfo

        results = formatter.format(
            variant="NC_000005.10:g.140114829del",
            genome="GRCh38",
            transcript_model="refseq",
            select_transcripts="all",
            liftover_level=True,
        )

        assert "NC_000005.10:g.140114829del" in results

    def test_format_missing_variant_keyword(self):
        formatter = simpleVariantFormatter.SimpleVariantFormatter()
        formatter.validator = self.vfo

        with pytest.raises(
                simpleVariantFormatter.FormatterSubmissionError,
                match="Variant is required"):
            formatter.format(
                genome="GRCh38",
            )

    def test_format_missing_genome_keyword(self):
        formatter = simpleVariantFormatter.SimpleVariantFormatter()
        formatter.validator = self.vfo

        with pytest.raises(
                simpleVariantFormatter.FormatterSubmissionError,
                match="Genome build is required"):
            formatter.format(
                variant="NC_000005.10:g.140114829del",
            )

# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later