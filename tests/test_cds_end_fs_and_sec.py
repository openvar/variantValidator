from VariantValidator import Validator
from unittest import TestCase


class TestCDSEndFS(TestCase):

    @classmethod
    def setup_class(cls):
        cls.vv = Validator()
        cls.vv.testing = True

    def test_deletion_1(self):
        variant = 'NM_001005242.3:c.2510_*1del'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_001005242.3:c.2510_*1del" in results.keys()
        assert results["NM_001005242.3:c.2510_*1del"]["hgvs_predicted_protein_consequence"
               ]["tlr"] == "NP_001005242.2:p.(Asp837GlyfsTer3)"

    def test_deletion_2(self):
        variant = 'NM_000074.2:c.782_*2del'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_000074.2:c.782_*2del" in results.keys()
        assert results["NM_000074.2:c.782_*2del"]["hgvs_predicted_protein_consequence"
               ]["tlr"] == "NP_000065.1:p.(Leu261GlnfsTer50)"

    def test_deletion_3(self):
        variant = 'NM_015474.3:c.1273_*2del'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_015474.3:c.1273_*2del" in results.keys()
        assert results["NM_015474.3:c.1273_*2del"]["hgvs_predicted_protein_consequence"
               ]["tlr"] == "NP_056289.2:p.(Asn425ValfsTer48)"

    def test_deletion_4(self):
        variant = 'NM_000330.3:c.658_*7del'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_000330.3:c.658_*7del" in results.keys()
        assert results["NM_000330.3:c.658_*7del"]["hgvs_predicted_protein_consequence"
               ]["tlr"] == "NP_000321.1:p.(Val220LeufsTer9)"

    def test_deletion_5(self):
        variant = 'NM_144573.3:c.2026_*1del'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_144573.3:c.2026_*1del" in results.keys()
        assert results["NM_144573.3:c.2026_*1del"]["hgvs_predicted_protein_consequence"
               ]["tlr"] == "NP_653174.3:p.(Ter676HisextTer8)"

    def test_deletion_6(self):
        variant = 'NM_000277.1:c.1357_*2del'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_000277.1:c.1357_*2del" in results.keys()
        assert results["NM_000277.1:c.1357_*2del"]["hgvs_predicted_protein_consequence"
               ]["tlr"] == "NP_000268.1:p.(Ter453ProextTer33)"

    def test_inversion_1(self):
        variant = 'NM_001005242.3:c.2510_*1inv'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_001005242.3:c.2510_*1inv" in results.keys()
        assert results["NM_001005242.3:c.2510_*1inv"]["hgvs_predicted_protein_consequence"
               ]["tlr"] == "NP_001005242.2:p.(Asp837AlafsTer3)"

    def test_issue_815a(self):
        results = self.vv.validate('NM_020451.3:c.447dup', 'GRCh38', 'all', liftover_level=True).format_as_dict(test=True)
        assert "NM_020451.3:c.447dup" in results.keys()
        assert results["NM_020451.3:c.447dup"]["hgvs_predicted_protein_consequence"]['slr'] == "NP_065184.2:p.(D150Ufs*2)"
        assert results["NM_020451.3:c.447dup"]["hgvs_predicted_protein_consequence"]['tlr'] == "NP_065184.2:p.(Asp150SecfsTer2)"

    def test_issue_815b(self):
        results = self.vv.validate('NM_020451.3:c.407del', 'GRCh38', 'all', liftover_level=True).format_as_dict(test=True)
        assert "NM_020451.3:c.407del" in results.keys()
        assert results["NM_020451.3:c.407del"]["hgvs_predicted_protein_consequence"]['slr'] == "NP_065184.2:p.(S136*)"
        assert results["NM_020451.3:c.407del"]["hgvs_predicted_protein_consequence"]['tlr'] == "NP_065184.2:p.(Ser136Ter)"

    def test_issue_815c(self):
        results = self.vv.validate('NM_020451.3:c.532del', 'GRCh38', 'all', liftover_level=True).format_as_dict(test=True)
        assert "NM_020451.3:c.532del" in results.keys()
        assert results["NM_020451.3:c.532del"]["hgvs_predicted_protein_consequence"]['slr'] == "NP_065184.2:p.(L178*)"
        assert results["NM_020451.3:c.532del"]["hgvs_predicted_protein_consequence"]['tlr'] == "NP_065184.2:p.(Leu178Ter)"

    def test_issue_815d(self):
        results = self.vv.validate('NM_020451.3:c.406_407insG', 'GRCh38', 'all', liftover_level=True).format_as_dict(test=True)
        assert "NM_020451.3:c.406_407insG" in results.keys()
        assert results["NM_020451.3:c.406_407insG"]["hgvs_predicted_protein_consequence"]['slr'] == "NP_065184.2:p.(S136Cfs*16)"
        assert results["NM_020451.3:c.406_407insG"]["hgvs_predicted_protein_consequence"]['tlr'] == "NP_065184.2:p.(Ser136CysfsTer16)"

    def test_issue_815e(self):
        # test that variants affected by selenocysteine edits inform users
        results = self.vv.validate('NM_020451.3:c.406_407insG', 'GRCh38', 'all', liftover_level=True).format_as_dict(test=True)
        print(results)
        assert "NM_020451.3:c.406_407insG" in results.keys()
        assert results["NM_020451.3:c.406_407insG"]["hgvs_predicted_protein_consequence"]['slr'] == "NP_065184.2:p.(S136Cfs*16)"
        assert results["NM_020451.3:c.406_407insG"]["hgvs_predicted_protein_consequence"]['tlr'] == "NP_065184.2:p.(Ser136CysfsTer16)"

        sec_warn_found = False
        for warn in results["NM_020451.3:c.406_407insG"]["validation_warnings"]:
            if warn.startswith('ProteinTranslationInfo: Sel'):
                sec_warn_found = True
        assert sec_warn_found

    def test_delins_1(self):
        variant = 'NM_001005242.3:c.2510_*1delinsTTT'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_001005242.3:c.2510_*1delinsTTT" in results.keys()
        assert results["NM_001005242.3:c.2510_*1delinsTTT"]["hgvs_predicted_protein_consequence"
               ]["tlr"] == "NP_001005242.2:p.(Asp837ValfsTer2)"

    def test_delins_2(self):
        variant = 'NM_001005242.3:c.2510_*1delinsTT'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_001005242.3:c.2510_*1delinsTT" in results.keys()
        assert results["NM_001005242.3:c.2510_*1delinsTT"]["hgvs_predicted_protein_consequence"
               ]["tlr"] == "NP_001005242.2:p.(Asp837ValfsTer49)"

    def test_delins_3(self):
        variant = 'NM_001005242.3:c.2510_*1delinsT'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_001005242.3:c.2510_*1delinsT" in results.keys()
        assert results["NM_001005242.3:c.2510_*1delinsT"]["hgvs_predicted_protein_consequence"
               ]["tlr"] == "NP_001005242.2:p.(Asp837ValfsTer43)"