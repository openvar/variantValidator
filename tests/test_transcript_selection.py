from VariantValidator import Validator
from unittest import TestCase


class TestTranscriptSelection(TestCase):

    @classmethod
    def setup_class(cls):
        cls.vv = Validator()
        cls.vv.testing = True

    def test_all_transcripts(self):
        variant = 'NC_000017.11:g.7676594T>A'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert len(results.keys()) >= 20

    def test_selected_tx(self):
        variant = 'NC_000017.11:g.7676594T>A'
        results = self.vv.validate(variant, 'GRCh38', 'NM_000546.6|NM_000546.5').format_as_dict(test=True)
        print(results)
        assert len(results.keys()) == 4

    def test_all_select(self):
        variant = 'NC_000017.11:g.7676594T>A'
        results = self.vv.validate(variant, 'GRCh38', 'select').format_as_dict(test=True)
        print(results)
        assert len(results.keys()) >= 4

    def test_mane_select(self):
        variant = 'NC_000017.11:g.7676594T>A'
        results = self.vv.validate(variant, 'GRCh38', 'mane_select').format_as_dict(test=True)
        print(results)
        assert len(results.keys()) == 3

    def test_mane(self):
        variant = 'NC_000008.11:g.22165406A>T'
        results = self.vv.validate(variant, 'GRCh38', 'mane').format_as_dict(test=True)
        print(results)
        assert len(results.keys()) >= 4

    def test_latest_transcripts(self):
        variant = 'NC_000010.10:g.11330462A>G'
        self.vv.testing = True
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        self.vv.testing = False
        results2 = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        assert len(results.keys()) > len(results2.keys())
        self.vv.testing = True
        assert self.vv.testing is True

