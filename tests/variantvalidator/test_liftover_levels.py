from unittest import TestCase

from VariantValidator import Validator
from VariantValidator.modules.vvMixinCore import ValidatorSubmissionError


class TestLiftoverLevel(TestCase):

    @classmethod
    def setup_class(cls):
        cls.vv = Validator()

    VARIANT = "NM_001040114.1:c.3055_3056inv"
    KEY = "NM_001040114.1:c.3055_3056inv"

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
