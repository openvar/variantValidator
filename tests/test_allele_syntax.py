from VariantValidator import Validator
from unittest import TestCase


class TestAlleleSyntax(TestCase):

    @classmethod
    def setup_class(cls):
        cls.vv = Validator()
        cls.vv.testing = True

    def test_variant1(self):
        variant = 'LRG_199t1:c.[2376G>C];[3103del]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_004006.2:c.2376G>C" in results.keys()
        assert "NM_004006.2:c.3103del" in results.keys()

    def test_variant2(self):
        variant = 'LRG_199t1:c.[4358_4359del;4361_4372del]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_004006.2:c.4358_4359del" in results.keys()
        assert "NM_004006.2:c.4362_4373del" in results.keys()

    def test_variant3(self):
        variant = 'LRG_199t1:c.2376G>C(;)3103del'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_004006.2:c.2376G>C" in results.keys()
        assert "NM_004006.2:c.3103del" in results.keys()

    def test_variant4(self):
        variant = 'LRG_199t1:c.2376[G>C];[(G>C)]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_004006.2:c.2376G>C" in results.keys()

    def test_variant5(self):
        variant = 'LRG_199t1:c.[2376G>C];[?]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_004006.2:c.2376G>C" in results.keys()
        # Note: Needs an update with the incorporation of uncertain positions

    def test_variant6(self):
        variant = 'LRG_199t1:c.[296T>G;476T=];[476T=](;)1083A>C'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_004006.2:c.1083A>C" in results.keys()
        assert "NM_004006.2:c.296T>G" in results.keys()
        assert "NM_004006.2:c.476=" in results.keys()

    def test_variant7(self):
        variant = 'LRG_199t1:c.[296T>G];[476T>C](;)1083A>C(;)1406del'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_004006.2:c.1083A>C" in results.keys()
        assert "NM_004006.2:c.1408del" in results.keys()
        assert "NM_004006.2:c.296T>G" in results.keys()
        assert "NM_004006.2:c.476T>C" in results.keys()

    def test_variant8(self):
        variant = 'LRG_199t1:c.[976-20T>A;976-17_976-1dup]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert "AlleleSyntaxError: Intronic variants can only be validated if a genomic/gene reference sequence " \
               "is also provided e.g. NC_000017.11(NM_000088.3):c.589-1G>T" in \
               results["validation_warning_1"]["validation_warnings"]

    def test_variant9(self):
        variant = 'NC_000023.10(LRG_199t1):c.[976-20T>A;976-17_976-1dup]'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)
        assert "Reference sequence LRG_199t1 updated to NM_004006.2" in results[
            "validation_warning_1"]["validation_warnings"]
        assert "ExonBoundaryError: Position c.976-17 does not correspond with an exon boundary " \
               "for transcript NM_004006.2" in results[
            "validation_warning_1"]["validation_warnings"]

    def test_variant10(self):
        variant = 'NM_004006.2:c.[145C>T;147C>G]'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert ("AlleleSyntaxError: Variants [145C>T;147C>G] should be merged into NM_004006.2:c.145_147delinsTGG" in
                results["validation_warning_1"]["validation_warnings"])

    def test_variant11(self):
        variant = 'NM_000059.4:c.[1916dup;1929del]'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert ("AlleleSyntaxError: Merging variants [1916dup;1929del] restores the original reading frame, "
                "so should be described as NM_000059.4:c.1917_1929delinsTGCATTCTTCTGT" in
                results["validation_warning_1"]["validation_warnings"])

