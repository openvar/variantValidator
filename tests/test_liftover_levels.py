from VariantValidator import Validator
from VariantValidator.modules.vvMixinCore import ValidatorSubmissionError
from VariantFormatter import simpleVariantFormatter
from VariantFormatter.simpleVariantFormatter import FormatterSubmissionError
from pprint import pprint

from unittest import TestCase


class TestLiftoverLevel(TestCase):

    @classmethod
    def setup_class(cls):
        cls.vv = Validator()
        cls.vf = simpleVariantFormatter.SimpleVariantFormatter(testing=True)

    VARIANT = "NM_001040114.1:c.3055_3056inv"
    KEY = "NM_001040114.1:c.3055_3056inv"

    GENOMIC_VARIANT = "NC_000016.10:g.15738651_15738652inv"
    GENOMIC_KEY = "NC_000016.10:g.15738651_15738652inv"
    TRANSCRIPT = "NM_002474.3"

    def test_liftover_true(self):
        results = self.vv.validate(
            self.VARIANT,
            "GRCh38",
            "all",
            liftover_level=True
        ).format_as_dict(test=True)

        entry = results[self.KEY]

        assert len(entry["alt_genomic_loci"]) > 0

        assert "grch37" in entry["primary_assembly_loci"]
        assert "hg19" in entry["primary_assembly_loci"]
        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

    def test_liftover_true_string(self):
        results = self.vv.validate(
            self.VARIANT,
            "GRCh38",
            "all",
            liftover_level="True"
        ).format_as_dict(test=True)

        entry = results[self.KEY]

        assert len(entry["alt_genomic_loci"]) > 0

        assert "grch37" in entry["primary_assembly_loci"]
        assert "hg19" in entry["primary_assembly_loci"]
        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

    def test_liftover_primary(self):
        results = self.vv.validate(
            self.VARIANT,
            "GRCh38",
            "all",
            liftover_level="primary"
        ).format_as_dict(test=True)

        entry = results[self.KEY]

        assert entry["alt_genomic_loci"] == []

        assert "grch37" in entry["primary_assembly_loci"]
        assert "hg19" in entry["primary_assembly_loci"]
        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

    def test_liftover_false(self):
        results = self.vv.validate(
            self.VARIANT,
            "GRCh38",
            "all",
            liftover_level=False
        ).format_as_dict(test=True)

        entry = results[self.KEY]

        assert entry["alt_genomic_loci"] == []

        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

        assert "grch37" not in entry["primary_assembly_loci"]
        assert "hg19" not in entry["primary_assembly_loci"]

    def test_liftover_false_string(self):
        results = self.vv.validate(
            self.VARIANT,
            "GRCh38",
            "all",
            liftover_level="False"
        ).format_as_dict(test=True)

        entry = results[self.KEY]

        assert entry["alt_genomic_loci"] == []

        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

        assert "grch37" not in entry["primary_assembly_loci"]
        assert "hg19" not in entry["primary_assembly_loci"]

    def test_liftover_none(self):
        results = self.vv.validate(
            self.VARIANT,
            "GRCh38",
            "all",
            liftover_level=None
        ).format_as_dict(test=True)

        entry = results[self.KEY]

        assert entry["alt_genomic_loci"] == []

        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

        assert "grch37" not in entry["primary_assembly_loci"]
        assert "hg19" not in entry["primary_assembly_loci"]

    def test_liftover_invalid_string(self):
        with self.assertRaises(ValidatorSubmissionError):
            self.vv.validate(
                self.VARIANT,
                "GRCh38",
                "all",
                liftover_level="banana"
            )


    def test_liftover_integer_1(self):
        results = self.vv.validate(
            self.VARIANT,
            "GRCh38",
            "all",
            liftover_level=1
        ).format_as_dict(test=True)

        entry = results[self.KEY]

        assert len(entry["alt_genomic_loci"]) > 0

        assert "grch37" in entry["primary_assembly_loci"]
        assert "hg19" in entry["primary_assembly_loci"]
        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]


    def test_liftover_integer_0(self):
        results = self.vv.validate(
            self.VARIANT,
            "GRCh38",
            "all",
            liftover_level=0
        ).format_as_dict(test=True)

        entry = results[self.KEY]

        assert entry["alt_genomic_loci"] == []

        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

        assert "grch37" not in entry["primary_assembly_loci"]
        assert "hg19" not in entry["primary_assembly_loci"]


    # ------------------------------------------------------------------
    # VariantFormatter object API
    # ------------------------------------------------------------------

    from pprint import pprint

    def _vf_entry(self, liftover_level):
        result = self.vf.format(
            variant=self.GENOMIC_VARIANT,
            genome="GRCh38",
            transcript_model="refseq",
            select_transcripts="mane_select",
            liftover_level=liftover_level
        )

        tx = result[self.GENOMIC_KEY][self.GENOMIC_KEY]["hgvs_t_and_p"][self.TRANSCRIPT]

        print("\n==============================")
        print("liftover_level =", liftover_level)
        print("==============================")

        pprint(tx)

        print("\nKeys:")
        pprint(tx.keys())

        print("\nalt_genomic_loci:")
        pprint(tx.get("alt_genomic_loci"))

        print("\nprimary_assembly_loci:")
        pprint(tx.get("primary_assembly_loci"))

        return tx

    def test_formatter_missing_variant(self):
        with self.assertRaises(FormatterSubmissionError):
            self.vf.format(
                genome="GRCh37"
            )

    def test_formatter_missing_genome(self):
        with self.assertRaises(FormatterSubmissionError):
            self.vf.format(
                variant=self.GENOMIC_VARIANT
            )

    def test_formatter_invalid_liftover(self):
        with self.assertRaises(FormatterSubmissionError):
            self.vf.format(
                variant=self.GENOMIC_VARIANT,
                genome="GRCh37",
                transcript_model="refseq",
                select_transcripts="mane_select",
                liftover_level="banana"
            )

    def test_formatter_liftover_true(self):
        entry = self._vf_entry(True)

        assert len(entry["alt_genomic_loci"]) > 0

        assert "grch37" in entry["primary_assembly_loci"]
        assert "hg19" in entry["primary_assembly_loci"]
        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

    def test_formatter_liftover_true_string(self):
        entry = self._vf_entry("True")

        assert len(entry["alt_genomic_loci"]) > 0

        assert "grch37" in entry["primary_assembly_loci"]
        assert "hg19" in entry["primary_assembly_loci"]
        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

    def test_formatter_liftover_integer_1(self):
        entry = self._vf_entry(1)

        assert len(entry["alt_genomic_loci"]) > 0

        assert "grch37" in entry["primary_assembly_loci"]
        assert "hg19" in entry["primary_assembly_loci"]
        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

    def test_formatter_liftover_primary(self):
        entry = self._vf_entry("primary")

        assert entry["alt_genomic_loci"] == []

        assert "grch37" in entry["primary_assembly_loci"]
        assert "hg19" in entry["primary_assembly_loci"]
        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

    def test_formatter_liftover_false(self):
        entry = self._vf_entry(False)

        assert entry["alt_genomic_loci"] == []

        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

        assert "grch37" not in entry["primary_assembly_loci"]
        assert "hg19" not in entry["primary_assembly_loci"]

    def test_formatter_liftover_false_string(self):
        entry = self._vf_entry("False")

        assert entry["alt_genomic_loci"] == []

        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

        assert "grch37" not in entry["primary_assembly_loci"]
        assert "hg19" not in entry["primary_assembly_loci"]

    def test_formatter_liftover_none(self):
        entry = self._vf_entry(None)

        assert entry["alt_genomic_loci"] == []

        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

        assert "grch37" not in entry["primary_assembly_loci"]
        assert "hg19" not in entry["primary_assembly_loci"]

    def test_formatter_liftover_integer_0(self):
        entry = self._vf_entry(0)

        assert entry["alt_genomic_loci"] == []

        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

        assert "grch37" not in entry["primary_assembly_loci"]
        assert "hg19" not in entry["primary_assembly_loci"]


    # ------------------------------------------------------------------
    # VariantFormatter legacy API
    # ------------------------------------------------------------------

    def test_formatter_legacy_missing_variant(self):
        with self.assertRaises(FormatterSubmissionError):
            simpleVariantFormatter.format(
                genome="GRCh38"
            )


    def test_formatter_legacy_missing_genome(self):
        with self.assertRaises(FormatterSubmissionError):
            simpleVariantFormatter.format(
                variant=self.GENOMIC_VARIANT
            )


    def test_formatter_legacy_invalid_liftover(self):
        with self.assertRaises(FormatterSubmissionError):
            simpleVariantFormatter.format(
                variant=self.GENOMIC_VARIANT,
                genome="GRCh38",
                transcript_model="refseq",
                select_transcripts="mane_select",
                liftover_level="banana"
            )


    def test_formatter_legacy_liftover_true(self):
        result = simpleVariantFormatter.format(
            variant=self.GENOMIC_VARIANT,
            genome="GRCh38",
            transcript_model="refseq",
            select_transcripts="mane_select",
            liftover_level=True
        )

        entry = result[self.GENOMIC_KEY][self.GENOMIC_KEY]["hgvs_t_and_p"][self.TRANSCRIPT]

        assert len(entry["alt_genomic_loci"]) > 0

        assert "grch37" in entry["primary_assembly_loci"]
        assert "hg19" in entry["primary_assembly_loci"]
        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]


    def test_formatter_legacy_liftover_true_string(self):
        result = simpleVariantFormatter.format(
            variant=self.GENOMIC_VARIANT,
            genome="GRCh38",
            transcript_model="refseq",
            select_transcripts="mane_select",
            liftover_level="True"
        )

        entry = result[self.GENOMIC_KEY][self.GENOMIC_KEY]["hgvs_t_and_p"][self.TRANSCRIPT]

        assert len(entry["alt_genomic_loci"]) > 0

        assert "grch37" in entry["primary_assembly_loci"]
        assert "hg19" in entry["primary_assembly_loci"]
        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]


    def test_formatter_legacy_liftover_integer_1(self):
        result = simpleVariantFormatter.format(
            variant=self.GENOMIC_VARIANT,
            genome="GRCh38",
            transcript_model="refseq",
            select_transcripts="mane_select",
            liftover_level=1
        )

        entry = result[self.GENOMIC_KEY][self.GENOMIC_KEY]["hgvs_t_and_p"][self.TRANSCRIPT]

        assert len(entry["alt_genomic_loci"]) > 0

        assert "grch37" in entry["primary_assembly_loci"]
        assert "hg19" in entry["primary_assembly_loci"]
        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]


    def test_formatter_legacy_liftover_primary(self):
        result = simpleVariantFormatter.format(
            variant=self.GENOMIC_VARIANT,
            genome="GRCh38",
            transcript_model="refseq",
            select_transcripts="mane_select",
            liftover_level="primary"
        )

        entry = result[self.GENOMIC_KEY][self.GENOMIC_KEY]["hgvs_t_and_p"][self.TRANSCRIPT]

        assert entry["alt_genomic_loci"] == []

        assert "grch37" in entry["primary_assembly_loci"]
        assert "hg19" in entry["primary_assembly_loci"]
        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]


    def test_formatter_legacy_liftover_false(self):
        result = simpleVariantFormatter.format(
            variant=self.GENOMIC_VARIANT,
            genome="GRCh38",
            transcript_model="refseq",
            select_transcripts="mane_select",
            liftover_level=False
        )

        entry = result[self.GENOMIC_KEY][self.GENOMIC_KEY]["hgvs_t_and_p"][self.TRANSCRIPT]

        assert entry["alt_genomic_loci"] == []

        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

        assert "grch37" not in entry["primary_assembly_loci"]
        assert "hg19" not in entry["primary_assembly_loci"]


    def test_formatter_legacy_liftover_false_string(self):
        result = simpleVariantFormatter.format(
            variant=self.GENOMIC_VARIANT,
            genome="GRCh38",
            transcript_model="refseq",
            select_transcripts="mane_select",
            liftover_level="False"
        )

        entry = result[self.GENOMIC_KEY][self.GENOMIC_KEY]["hgvs_t_and_p"][self.TRANSCRIPT]

        assert entry["alt_genomic_loci"] == []

        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

        assert "grch37" not in entry["primary_assembly_loci"]
        assert "hg19" not in entry["primary_assembly_loci"]


    def test_formatter_legacy_liftover_integer_0(self):
        result = simpleVariantFormatter.format(
            variant=self.GENOMIC_VARIANT,
            genome="GRCh38",
            transcript_model="refseq",
            select_transcripts="mane_select",
            liftover_level=0
        )

        entry = result[self.GENOMIC_KEY][self.GENOMIC_KEY]["hgvs_t_and_p"][self.TRANSCRIPT]

        assert entry["alt_genomic_loci"] == []

        assert "grch38" in entry["primary_assembly_loci"]
        assert "hg38" in entry["primary_assembly_loci"]

        assert "grch37" not in entry["primary_assembly_loci"]
        assert "hg19" not in entry["primary_assembly_loci"]

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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>
