from VariantValidator import Validator
from unittest import TestCase


class TestMethylCases(TestCase):

    @classmethod
    def setup_class(cls):
        cls.vv = Validator()
        cls.vv.testing = True

    def test_gom(self):
        variant = 'NC_000011.10:g.1999904_1999946|gom'
        data_dict = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)

        # Check if "|gom" is present in the hgvs_genomic_description for each variant
        assert "|gom" in data_dict["NM_001400176.1:c.498-11637_498-11595|gom"]["alt_genomic_loci"][0]["grch38"][
            "hgvs_genomic_description"], "Assertion failed: '|gom' not found in grch38 genomic description of the first variant."
        assert "|gom" in data_dict["NM_001400176.1:c.498-11637_498-11595|gom"]["alt_genomic_loci"][1]["hg38"][
            "hgvs_genomic_description"], "Assertion failed: '|gom' not found in hg38 genomic description of the first variant."

        assert "|gom" in data_dict["NR_131224.1:n.249+1272_249+1314|gom"]["alt_genomic_loci"][0]["grch38"][
            "hgvs_genomic_description"], "Assertion failed: '|gom' not found in grch38 genomic description of the second variant."
        assert "|gom" in data_dict["NR_131224.1:n.249+1272_249+1314|gom"]["alt_genomic_loci"][1]["hg38"][
            "hgvs_genomic_description"], "Assertion failed: '|gom' not found in hg38 genomic description of the second variant."

    def test_lom(self):
        variant = 'NC_000011.10:g.1999904_1999946|lom'
        data_dict = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)

        # Check if "|lom" is present in the hgvs_genomic_description for each variant
        assert "|lom" in data_dict["NM_001400176.1:c.498-11637_498-11595|lom"]["alt_genomic_loci"][0]["grch38"][
            "hgvs_genomic_description"], "Assertion failed: '|lom' not found in grch38 genomic description of the first variant."
        assert "|lom" in data_dict["NM_001400176.1:c.498-11637_498-11595|lom"]["alt_genomic_loci"][1]["hg38"][
            "hgvs_genomic_description"], "Assertion failed: '|lom' not found in hg38 genomic description of the first variant."

        assert "|lom" in data_dict["NR_131224.1:n.249+1272_249+1314|lom"]["alt_genomic_loci"][0]["grch38"][
            "hgvs_genomic_description"], "Assertion failed: '|lom' not found in grch38 genomic description of the second variant."
        assert "|lom" in data_dict["NR_131224.1:n.249+1272_249+1314|lom"]["alt_genomic_loci"][1]["hg38"][
            "hgvs_genomic_description"], "Assertion failed: '|lom' not found in hg38 genomic description of the second variant."

    def test_meteq(self):
        variant = 'NC_000011.10:g.1999904_1999946|met='
        data_dict = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        # Assuming your dictionary is stored in a variable named 'data_dict'

        # Check if "|met=" is present in the hgvs_genomic_description for each variant
        assert "|met=" in data_dict["NM_001400176.1:c.498-11637_498-11595|met="]["alt_genomic_loci"][0]["grch38"][
            "hgvs_genomic_description"], "Assertion failed: '|met=' not found in grch38 genomic description of the first variant."
        assert "|met=" in data_dict["NM_001400176.1:c.498-11637_498-11595|met="]["alt_genomic_loci"][1]["hg38"][
            "hgvs_genomic_description"], "Assertion failed: '|met=' not found in hg38 genomic description of the first variant."

        assert "|met=" in data_dict["NR_131224.1:n.249+1272_249+1314|met="]["alt_genomic_loci"][0]["grch38"][
            "hgvs_genomic_description"], "Assertion failed: '|met=' not found in grch38 genomic description of the second variant."
        assert "|met=" in data_dict["NR_131224.1:n.249+1272_249+1314|met="]["alt_genomic_loci"][1]["hg38"][
            "hgvs_genomic_description"], "Assertion failed: '|met=' not found in hg38 genomic description of the second variant."

