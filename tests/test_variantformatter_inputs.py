import VariantFormatter
import VariantFormatter.variantformatter as vf
import VariantFormatter.simpleVariantFormatter
import VariantValidator
vfo = VariantValidator.Validator()


class TestVFvariantsAuto(object):
    @classmethod
    def setup_class(cls):
        VariantFormatter.__version__
        vfo.testing = True

    def test_variant1(self):
        variant = 'NC_000019.10:g.50378563_50378564insTAC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000019.10:g.50378563_50378564insTAC' in results.keys()
        assert results['NC_000019.10:g.50378563_50378564insTAC']['p_vcf'] is None
        assert results['NC_000019.10:g.50378563_50378564insTAC']['g_hgvs'] is None
        assert results['NC_000019.10:g.50378563_50378564insTAC']['genomic_variant_error'] == 'chromosome ID NC_000019.10 is not associated with genome build GRCh37'
        assert results['NC_000019.10:g.50378563_50378564insTAC']['hgvs_t_and_p'] is None

    def test_variant2(self):
        variant = '11-5248232-A-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '11-5248232-A-T' in results.keys()
        assert results['11-5248232-A-T']['p_vcf'] is None
        assert results['11-5248232-A-T']['g_hgvs'] is None
        assert results['11-5248232-A-T']['genomic_variant_error'] == 'NC_000011.9:g.5248232A>T: Variant reference (A) does not agree with reference sequence (T)'
        assert results['11-5248232-A-T']['hgvs_t_and_p'] is None

    def test_variant3(self):
        variant = 'NC_000012.11:g.122064777A>C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000012.11:g.122064777A>C' in results.keys()
        assert results['NC_000012.11:g.122064777A>C']['p_vcf'] is None
        assert results['NC_000012.11:g.122064777A>C']['g_hgvs'] is None
        assert results['NC_000012.11:g.122064777A>C']['genomic_variant_error'] == 'NC_000012.11:g.122064777A>C: Variant reference (A) does not agree with reference sequence (C)'
        assert results['NC_000012.11:g.122064777A>C']['hgvs_t_and_p'] is None

    def test_variant4(self):
        variant = 'NC_000002.11:g.73613030C>T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000002.11:g.73613030C>T' in results.keys()
        assert results['NC_000002.11:g.73613030C>T']['p_vcf'] == '2:73613030:C:T'
        assert results['NC_000002.11:g.73613030C>T']['g_hgvs'] == 'NC_000002.11:g.73613030C>T'
        assert results['NC_000002.11:g.73613030C>T']['genomic_variant_error'] is None
        assert 'NM_015120.4' in results['NC_000002.11:g.73613030C>T']['hgvs_t_and_p'].keys()
        assert results['NC_000002.11:g.73613030C>T']['hgvs_t_and_p']['NM_015120.4']['t_hgvs'] == 'NM_015120.4:c.34C>T'
        assert results['NC_000002.11:g.73613030C>T']['hgvs_t_and_p']['NM_015120.4']['p_hgvs_tlc'] == 'NP_055935.4:p.(Leu12=)'
        assert results['NC_000002.11:g.73613030C>T']['hgvs_t_and_p']['NM_015120.4']['transcript_variant_error'] is None

    def test_variant5(self):
        variant = 'NC_000023.10:g.33229673A>T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.10:g.33229673A>T' in results.keys()
        assert results['NC_000023.10:g.33229673A>T']['p_vcf'] == 'X:33229673:A:T'
        assert results['NC_000023.10:g.33229673A>T']['g_hgvs'] == 'NC_000023.10:g.33229673A>T'
        assert results['NC_000023.10:g.33229673A>T']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.10:g.33229673A>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.10:g.33229673A>T']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.7+127703T>A'
        assert results['NC_000023.10:g.33229673A>T']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.10:g.33229673A>T']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.10:g.33229673A>T']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.10:g.33229673A>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.10:g.33229673A>T']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.-244T>A'
        assert results['NC_000023.10:g.33229673A>T']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.10:g.33229673A>T']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.10:g.33229673A>T']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None

    def test_variant6(self):
        variant = 'NC_000017.10:g.48279242G>T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000017.10:g.48279242G>T' in results.keys()
        assert results['NC_000017.10:g.48279242G>T']['p_vcf'] == '17:48279242:G:T'
        assert results['NC_000017.10:g.48279242G>T']['g_hgvs'] == 'NC_000017.10:g.48279242G>T'
        assert results['NC_000017.10:g.48279242G>T']['genomic_variant_error'] is None
        assert results['NC_000017.10:g.48279242G>T']['hgvs_t_and_p'] == {'intergenic': {'alt_genomic_loci': None}}

    def test_variant7(self):
        variant = 'NC_000017.10:g.48261457_48261463TTATGTT='
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000017.10:g.48261457_48261463TTATGTT=' in results.keys()
        assert results['NC_000017.10:g.48261457_48261463TTATGTT=']['p_vcf'] == '17:48261457:TTATGTT:TTATGTT'
        assert results['NC_000017.10:g.48261457_48261463TTATGTT=']['g_hgvs'] == 'NC_000017.10:g.48261457_48261463='
        assert results['NC_000017.10:g.48261457_48261463TTATGTT=']['genomic_variant_error'] is None
        assert 'NM_000088.3' in results['NC_000017.10:g.48261457_48261463TTATGTT=']['hgvs_t_and_p'].keys()
        assert results['NC_000017.10:g.48261457_48261463TTATGTT=']['hgvs_t_and_p']['NM_000088.3']['t_hgvs'] == 'NM_000088.3:c.*1400_*1406='
        assert results['NC_000017.10:g.48261457_48261463TTATGTT=']['hgvs_t_and_p']['NM_000088.3']['p_hgvs_tlc'] == 'NP_000079.2:p.?'
        assert results['NC_000017.10:g.48261457_48261463TTATGTT=']['hgvs_t_and_p']['NM_000088.3']['p_hgvs_slc'] == 'NP_000079.2:p.?'
        assert results['NC_000017.10:g.48261457_48261463TTATGTT=']['hgvs_t_and_p']['NM_000088.3']['transcript_variant_error'] is None

    def test_variant8(self):
        variant = 'NC_000017.10:g.48275363C>A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000017.10:g.48275363C>A' in results.keys()
        assert results['NC_000017.10:g.48275363C>A']['p_vcf'] == '17:48275363:C:A'
        assert results['NC_000017.10:g.48275363C>A']['g_hgvs'] == 'NC_000017.10:g.48275363C>A'
        assert results['NC_000017.10:g.48275363C>A']['genomic_variant_error'] is None
        assert 'NM_000088.3' in results['NC_000017.10:g.48275363C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000017.10:g.48275363C>A']['hgvs_t_and_p']['NM_000088.3']['t_hgvs'] == 'NM_000088.3:c.589G>T'
        assert results['NC_000017.10:g.48275363C>A']['hgvs_t_and_p']['NM_000088.3']['p_hgvs_tlc'] == 'NP_000079.2:p.(Gly197Cys)'
        assert results['NC_000017.10:g.48275363C>A']['hgvs_t_and_p']['NM_000088.3']['p_hgvs_slc'] == 'NP_000079.2:p.(G197C)'
        assert results['NC_000017.10:g.48275363C>A']['hgvs_t_and_p']['NM_000088.3']['transcript_variant_error'] is None

    def test_variant9(self):
        variant = '11-5248232-T-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '11-5248232-T-A' in results.keys()
        assert results['11-5248232-T-A']['p_vcf'] == '11-5248232-T-A'
        assert results['11-5248232-T-A']['g_hgvs'] == 'NC_000011.9:g.5248232T>A'
        assert results['11-5248232-T-A']['genomic_variant_error'] is None
        assert 'NM_000518.5' in results['11-5248232-T-A']['hgvs_t_and_p'].keys()
        assert results['11-5248232-T-A']['hgvs_t_and_p']['NM_000518.5']['t_hgvs'] == 'NM_000518.5:c.20A>T'
        assert results['11-5248232-T-A']['hgvs_t_and_p']['NM_000518.5']['p_hgvs_tlc'] == 'NP_000509.1:p.(Glu7Val)'
        assert results['11-5248232-T-A']['hgvs_t_and_p']['NM_000518.5']['p_hgvs_slc'] == 'NP_000509.1:p.(E7V)'
        assert results['11-5248232-T-A']['hgvs_t_and_p']['NM_000518.5']['transcript_variant_error'] is None
        assert 'NM_000518.4' in results['11-5248232-T-A']['hgvs_t_and_p'].keys()
        assert results['11-5248232-T-A']['hgvs_t_and_p']['NM_000518.4']['t_hgvs'] == 'NM_000518.4:c.20A>T'
        assert results['11-5248232-T-A']['hgvs_t_and_p']['NM_000518.4']['p_hgvs_tlc'] == 'NP_000509.1:p.(Glu7Val)'
        assert results['11-5248232-T-A']['hgvs_t_and_p']['NM_000518.4']['p_hgvs_slc'] == 'NP_000509.1:p.(E7V)'
        assert results['11-5248232-T-A']['hgvs_t_and_p']['NM_000518.4']['transcript_variant_error'] is None

    def test_variant10(self):
        variant = '1-150550916-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-150550916-G-A' in results.keys()
        assert results['1-150550916-G-A']['p_vcf'] == '1-150550916-G-A'
        assert results['1-150550916-G-A']['g_hgvs'] == 'NC_000001.10:g.150550916G>A'
        assert results['1-150550916-G-A']['genomic_variant_error'] is None
        assert 'NM_001197320.1' in results['1-150550916-G-A']['hgvs_t_and_p'].keys()
        assert results['1-150550916-G-A']['hgvs_t_and_p']['NM_001197320.1']['t_hgvs'] == 'NM_001197320.1:c.281C>T'
        assert results['1-150550916-G-A']['hgvs_t_and_p']['NM_001197320.1']['p_hgvs_tlc'] == 'NP_001184249.1:p.(Ser94Phe)'
        assert results['1-150550916-G-A']['hgvs_t_and_p']['NM_001197320.1']['p_hgvs_slc'] == 'NP_001184249.1:p.(S94F)'
        assert results['1-150550916-G-A']['hgvs_t_and_p']['NM_001197320.1']['transcript_variant_error'] is None
        assert 'NM_021960.4' in results['1-150550916-G-A']['hgvs_t_and_p'].keys()
        assert results['1-150550916-G-A']['hgvs_t_and_p']['NM_021960.4']['t_hgvs'] == 'NM_021960.4:c.740C>T'
        assert results['1-150550916-G-A']['hgvs_t_and_p']['NM_021960.4']['p_hgvs_tlc'] == 'NP_068779.1:p.(Ser247Phe)'
        assert results['1-150550916-G-A']['hgvs_t_and_p']['NM_021960.4']['p_hgvs_slc'] == 'NP_068779.1:p.(S247F)'
        assert results['1-150550916-G-A']['hgvs_t_and_p']['NM_021960.4']['transcript_variant_error'] is None
        assert 'NM_182763.2' in results['1-150550916-G-A']['hgvs_t_and_p'].keys()
        assert results['1-150550916-G-A']['hgvs_t_and_p']['NM_182763.2']['t_hgvs'] == "NM_182763.2:c.688+403C>T"
        assert results['1-150550916-G-A']['hgvs_t_and_p']['NM_182763.2']['p_hgvs_tlc'] == "NP_877495.1:p.?"
        assert results['1-150550916-G-A']['hgvs_t_and_p']['NM_182763.2']['p_hgvs_slc'] == "NP_877495.1:p.?"
        assert results['1-150550916-G-A']['hgvs_t_and_p']['NM_182763.2']['transcript_variant_error'] is None

    def test_variant11(self):
        variant = 'HSCHR6_MHC_SSTO_CTG1-3852542-C-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'HSCHR6_MHC_SSTO_CTG1-3852542-C-G' in results.keys()
        assert results['HSCHR6_MHC_SSTO_CTG1-3852542-C-G']['p_vcf'] == 'HSCHR6_MHC_SSTO_CTG1-3852542-C-G'
        assert results['HSCHR6_MHC_SSTO_CTG1-3852542-C-G']['g_hgvs'] == 'NT_167249.1:g.3852542C>G'
        assert results['HSCHR6_MHC_SSTO_CTG1-3852542-C-G']['genomic_variant_error'] is None
        assert 'NM_021983.4' in results['HSCHR6_MHC_SSTO_CTG1-3852542-C-G']['hgvs_t_and_p'].keys()
        assert results['HSCHR6_MHC_SSTO_CTG1-3852542-C-G']['hgvs_t_and_p']['NM_021983.4']['t_hgvs'] == 'NM_021983.4:c.490G>C'
        assert results['HSCHR6_MHC_SSTO_CTG1-3852542-C-G']['hgvs_t_and_p']['NM_021983.4']['p_hgvs_tlc'] == 'NP_068818.4:p.(Gly164Arg)'
        assert results['HSCHR6_MHC_SSTO_CTG1-3852542-C-G']['hgvs_t_and_p']['NM_021983.4']['p_hgvs_slc'] == 'NP_068818.4:p.(G164R)'
        assert results['HSCHR6_MHC_SSTO_CTG1-3852542-C-G']['hgvs_t_and_p']['NM_021983.4']['transcript_variant_error'] is None

    def test_variant12(self):
        variant = 'NC_000013.10:g.32929387T>C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000013.10:g.32929387T>C' in results.keys()
        assert results['NC_000013.10:g.32929387T>C']['p_vcf'] == '13:32929387:T:C'
        assert results['NC_000013.10:g.32929387T>C']['g_hgvs'] == 'NC_000013.10:g.32929387T>C'
        assert results['NC_000013.10:g.32929387T>C']['genomic_variant_error'] is None
        assert 'NM_000059.3' in results['NC_000013.10:g.32929387T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000013.10:g.32929387T>C']['hgvs_t_and_p']['NM_000059.3']['t_hgvs'] == 'NM_000059.3:c.7397='
        assert results['NC_000013.10:g.32929387T>C']['hgvs_t_and_p']['NM_000059.3']['p_hgvs_tlc'] == 'NP_000050.2:p.(Ala2466=)'
        assert results['NC_000013.10:g.32929387T>C']['hgvs_t_and_p']['NM_000059.3']['p_hgvs_slc'] == 'NP_000050.2:p.(A2466=)'
        assert results['NC_000013.10:g.32929387T>C']['hgvs_t_and_p']['NM_000059.3']['transcript_variant_error'] is None

    def test_variant13(self):
        variant = '19-41123094-G-GG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '19-41123094-G-GG' in results.keys()
        assert results['19-41123094-G-GG']['p_vcf'] == '19-41123093-A-AG'
        assert results['19-41123094-G-GG']['g_hgvs'] == 'NC_000019.9:g.41123095dup'
        assert results['19-41123094-G-GG']['genomic_variant_error'] is None
        assert 'NM_001042544.1' in results['19-41123094-G-GG']['hgvs_t_and_p'].keys()
        assert results['19-41123094-G-GG']['hgvs_t_and_p']['NM_001042544.1']['t_hgvs'] == 'NM_001042544.1:c.3233_3235='
        assert results['19-41123094-G-GG']['hgvs_t_and_p']['NM_001042544.1']['p_hgvs_tlc'] == 'NP_001036009.1:p.(Gln1078_Gly1079=)'
        assert results['19-41123094-G-GG']['hgvs_t_and_p']['NM_001042544.1']['transcript_variant_error'] is None
        assert 'NM_003573.2' in results['19-41123094-G-GG']['hgvs_t_and_p'].keys()
        assert results['19-41123094-G-GG']['hgvs_t_and_p']['NM_003573.2']['t_hgvs'] == 'NM_003573.2:c.3122_3124='
        assert results['19-41123094-G-GG']['hgvs_t_and_p']['NM_003573.2']['p_hgvs_tlc'] == 'NP_003564.2:p.(Gln1041_Gly1042=)'
        assert results['19-41123094-G-GG']['hgvs_t_and_p']['NM_003573.2']['transcript_variant_error'] is None
        assert 'NM_001042545.1' in results['19-41123094-G-GG']['hgvs_t_and_p'].keys()
        assert results['19-41123094-G-GG']['hgvs_t_and_p']['NM_001042545.1']['t_hgvs'] == 'NM_001042545.1:c.3032_3034='
        assert results['19-41123094-G-GG']['hgvs_t_and_p']['NM_001042545.1']['p_hgvs_tlc'] == 'NP_001036010.1:p.(Gln1011_Gly1012=)'
        assert results['19-41123094-G-GG']['hgvs_t_and_p']['NM_001042545.1']['transcript_variant_error'] is None

    def test_variant14(self):
        variant = '15-72105928-AC-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '15-72105928-AC-A' in results.keys()
        assert results['15-72105928-AC-A']['p_vcf'] == '15-72105928-AC-A'
        assert results['15-72105928-AC-A']['g_hgvs'] == 'NC_000015.9:g.72105933del'
        assert results['15-72105928-AC-A']['genomic_variant_error'] is None
        assert 'NM_014249.3' in results['15-72105928-AC-A']['hgvs_t_and_p'].keys()
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_014249.3']['t_hgvs'] == 'NM_014249.3:c.947_948='
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_014249.3']['p_hgvs_tlc'] == 'NP_055064.1:p.(Asp316=)'
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_014249.3']['p_hgvs_slc'] == 'NP_055064.1:p.(D316=)'
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_014249.3']['transcript_variant_error'] is None
        assert 'NM_014249.2' in results['15-72105928-AC-A']['hgvs_t_and_p'].keys()
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_014249.2']['t_hgvs'] == 'NM_014249.2:c.947_948='
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_014249.2']['p_hgvs_tlc'] == 'NP_055064.1:p.(Asp316=)'
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_014249.2']['p_hgvs_slc'] == 'NP_055064.1:p.(D316=)'
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_014249.2']['transcript_variant_error'] is None
        assert 'NM_016346.3' in results['15-72105928-AC-A']['hgvs_t_and_p'].keys()
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_016346.3']['t_hgvs'] == 'NM_016346.3:c.947_948='
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_016346.3']['p_hgvs_tlc'] == 'NP_057430.1:p.(Asp316=)'
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_016346.3']['p_hgvs_slc'] == 'NP_057430.1:p.(D316=)'
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_016346.3']['transcript_variant_error'] is None
        assert 'NM_016346.2' in results['15-72105928-AC-A']['hgvs_t_and_p'].keys()
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_016346.2']['t_hgvs'] == 'NM_016346.2:c.947_948='
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_016346.2']['p_hgvs_tlc'] == 'NP_057430.1:p.(Asp316=)'
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_016346.2']['p_hgvs_slc'] == 'NP_057430.1:p.(D316=)'
        assert results['15-72105928-AC-A']['hgvs_t_and_p']['NM_016346.2']['transcript_variant_error'] is None

    def test_variant15(self):
        variant = '12-122064773-CCCGCCA-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '12-122064773-CCCGCCA-C' in results.keys()
        assert results['12-122064773-CCCGCCA-C']['p_vcf'] == '12-122064773-CCCGCCA-C'
        assert results['12-122064773-CCCGCCA-C']['g_hgvs'] == 'NC_000012.11:g.122064785_122064790del'
        assert results['12-122064773-CCCGCCA-C']['genomic_variant_error'] is None
        assert 'NM_032790.3' in results['12-122064773-CCCGCCA-C']['hgvs_t_and_p'].keys()
        assert results['12-122064773-CCCGCCA-C']['hgvs_t_and_p']['NM_032790.3']['t_hgvs'] == 'NM_032790.3:c.126_127='
        assert results['12-122064773-CCCGCCA-C']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_tlc'] == 'NP_116179.2:p.(Ala42_Pro43=)'
        assert results['12-122064773-CCCGCCA-C']['hgvs_t_and_p']['NM_032790.3']['transcript_variant_error'] is None

    def test_variant16(self):
        variant = '12-122064774-CCGCCA-CCGCCA'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '12-122064774-CCGCCA-CCGCCA' in results.keys()
        assert results['12-122064774-CCGCCA-CCGCCA']['p_vcf'] == '12-122064774-CCGCCA-CCGCCA'
        assert results['12-122064774-CCGCCA-CCGCCA']['g_hgvs'] == 'NC_000012.11:g.122064774_122064779='
        assert results['12-122064774-CCGCCA-CCGCCA']['genomic_variant_error'] is None
        assert 'NM_032790.3' in results['12-122064774-CCGCCA-CCGCCA']['hgvs_t_and_p'].keys()
        assert results['12-122064774-CCGCCA-CCGCCA']['hgvs_t_and_p']['NM_032790.3']['t_hgvs'] == 'NM_032790.3:c.132_137dup'
        assert results['12-122064774-CCGCCA-CCGCCA']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_tlc'] == 'NP_116179.2:p.(Pro46_Pro47dup)'
        assert results['12-122064774-CCGCCA-CCGCCA']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_slc'] == 'NP_116179.2:p.(P46_P47dup)'
        assert results['12-122064774-CCGCCA-CCGCCA']['hgvs_t_and_p']['NM_032790.3']['transcript_variant_error'] is None

    # def test_variant17(self):
    #     variant = '12-122064773-CCCGCCACCGCCACCGC-CCCGCCACCGCCGCCGTC'
    #     results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
    #     results = results.stucture_data()
    #     print(results)
    #     assert '12-122064773-CCCGCCACCGCCACCGC-CCCGCCACCGCCGCCGTC' in results.keys()
    #     assert results['12-122064773-CCCGCCACCGCCACCGC-CCCGCCACCGCCGCCGTC']['p_vcf'] == '12-122064773-CCCGCCACCGCCACCGC-CCCGCCACCGCCGCCGTC'
    #     assert results['12-122064773-CCCGCCACCGCCACCGC-CCCGCCACCGCCGCCGTC']['g_hgvs'] == 'NC_000012.11:g.122064785_122064788delinsGCCGT'
    #     assert results['12-122064773-CCCGCCACCGCCACCGC-CCCGCCACCGCCGCCGTC']['genomic_variant_error'] is None
    #     assert 'NM_032790.3' in results['12-122064773-CCCGCCACCGCCACCGC-CCCGCCACCGCCGCCGTC']['hgvs_t_and_p'].keys()
    #     assert results['12-122064773-CCCGCCACCGCCACCGC-CCCGCCACCGCCGCCGTC']['hgvs_t_and_p']['NM_032790.3']['t_hgvs'] == 'NM_032790.3:c.132_135delinsGCCGT'
    #     assert results['12-122064773-CCCGCCACCGCCACCGC-CCCGCCACCGCCGCCGTC']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_tlc'] == 'NP_116179.2:p.(Pro46SerfsTer42)'
    #     assert results['12-122064773-CCCGCCACCGCCACCGC-CCCGCCACCGCCGCCGTC']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_slc'] == 'NP_116179.2:p.(P46Sfs*42)'
    #     assert results['12-122064773-CCCGCCACCGCCACCGC-CCCGCCACCGCCGCCGTC']['hgvs_t_and_p']['NM_032790.3']['transcript_variant_error'] is None

    def test_variant18(self):
        variant = 'NC_000012.11:g.122064777C>A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000012.11:g.122064777C>A' in results.keys()
        assert results['NC_000012.11:g.122064777C>A']['p_vcf'] == '12:122064777:C:A'
        assert results['NC_000012.11:g.122064777C>A']['g_hgvs'] == 'NC_000012.11:g.122064777C>A'
        assert results['NC_000012.11:g.122064777C>A']['genomic_variant_error'] is None
        assert 'NM_032790.3' in results['NC_000012.11:g.122064777C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000012.11:g.122064777C>A']['hgvs_t_and_p']['NM_032790.3']['t_hgvs'] == "NM_032790.3:c.129_130insACACCG"
        assert results['NC_000012.11:g.122064777C>A']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_tlc'] == "NP_116179.2:p.(Pro43_Pro44insThrPro)"
        assert results['NC_000012.11:g.122064777C>A']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_slc'] == "NP_116179.2:p.(P43_P44insTP)"
        assert results['NC_000012.11:g.122064777C>A']['hgvs_t_and_p']['NM_032790.3']['transcript_variant_error'] is None

    def test_variant19(self):
        variant = 'NC_000012.11:g.122064776delG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000012.11:g.122064776delG' in results.keys()
        assert results['NC_000012.11:g.122064776delG']['p_vcf'] == '12:122064775:CG:C'
        assert results['NC_000012.11:g.122064776delG']['g_hgvs'] == 'NC_000012.11:g.122064776del'
        assert results['NC_000012.11:g.122064776delG']['genomic_variant_error'] is None
        assert 'NM_032790.3' in results['NC_000012.11:g.122064776delG']['hgvs_t_and_p'].keys()
        assert results['NC_000012.11:g.122064776delG']['hgvs_t_and_p']['NM_032790.3']['t_hgvs'] == "NM_032790.3:c.128_129insCCACC"
        assert results['NC_000012.11:g.122064776delG']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_tlc'] == "NP_116179.2:p.(Pro44HisfsTer22)"
        assert results['NC_000012.11:g.122064776delG']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_slc'] == "NP_116179.2:p.(P44Hfs*22)"
        assert results['NC_000012.11:g.122064776delG']['hgvs_t_and_p']['NM_032790.3']['transcript_variant_error'] is None

    def test_variant20(self):
        variant = 'NC_000012.11:g.122064776dupG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000012.11:g.122064776dupG' in results.keys()
        assert results['NC_000012.11:g.122064776dupG']['p_vcf'] == '12:122064775:C:CG'
        assert results['NC_000012.11:g.122064776dupG']['g_hgvs'] == 'NC_000012.11:g.122064776dup'
        assert results['NC_000012.11:g.122064776dupG']['genomic_variant_error'] is None
        assert 'NM_032790.3' in results['NC_000012.11:g.122064776dupG']['hgvs_t_and_p'].keys()
        assert results['NC_000012.11:g.122064776dupG']['hgvs_t_and_p']['NM_032790.3']['t_hgvs'] == "NM_032790.3:c.129_130insGCCACCG"
        assert results['NC_000012.11:g.122064776dupG']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_tlc'] == "NP_116179.2:p.(Pro44AlafsTer46)"
        assert results['NC_000012.11:g.122064776dupG']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_slc'] == "NP_116179.2:p.(P44Afs*46)"
        assert results['NC_000012.11:g.122064776dupG']['hgvs_t_and_p']['NM_032790.3']['transcript_variant_error'] is None

    def test_variant21(self):
        variant = 'NC_000012.11:g.122064776_122064777insTTT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000012.11:g.122064776_122064777insTTT' in results.keys()
        assert results['NC_000012.11:g.122064776_122064777insTTT']['p_vcf'] == '12:122064776:G:GTTT'
        assert results['NC_000012.11:g.122064776_122064777insTTT']['g_hgvs'] == 'NC_000012.11:g.122064776_122064777insTTT'
        assert results['NC_000012.11:g.122064776_122064777insTTT']['genomic_variant_error'] is None
        assert 'NM_032790.3' in results['NC_000012.11:g.122064776_122064777insTTT']['hgvs_t_and_p'].keys()
        assert results['NC_000012.11:g.122064776_122064777insTTT']['hgvs_t_and_p']['NM_032790.3']['t_hgvs'] == "NM_032790.3:c.129_130insTTTCCACCG"
        assert results['NC_000012.11:g.122064776_122064777insTTT']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_tlc'] == "NP_116179.2:p.(Pro43_Pro44insPheProPro)"
        assert results['NC_000012.11:g.122064776_122064777insTTT']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_slc'] == "NP_116179.2:p.(P43_P44insFPP)"
        assert results['NC_000012.11:g.122064776_122064777insTTT']['hgvs_t_and_p']['NM_032790.3']['transcript_variant_error'] is None

    def test_variant22(self):
        variant = 'NC_000012.11:g.122064772_122064775del'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000012.11:g.122064772_122064775del' in results.keys()
        assert results['NC_000012.11:g.122064772_122064775del']['p_vcf'] == '12:122064771:GCCCC:G'
        assert results['NC_000012.11:g.122064772_122064775del']['g_hgvs'] == 'NC_000012.11:g.122064772_122064775del'
        assert results['NC_000012.11:g.122064772_122064775del']['genomic_variant_error'] is None
        assert 'NM_032790.3' in results['NC_000012.11:g.122064772_122064775del']['hgvs_t_and_p'].keys()
        assert results['NC_000012.11:g.122064772_122064775del']['hgvs_t_and_p']['NM_032790.3']['t_hgvs'] == "NM_032790.3:c.125_126delinsGCCA"
        assert results['NC_000012.11:g.122064772_122064775del']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_tlc'] == "NP_116179.2:p.(Ala42GlyfsTer23)"
        assert results['NC_000012.11:g.122064772_122064775del']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_slc'] == "NP_116179.2:p.(A42Gfs*23)"
        assert results['NC_000012.11:g.122064772_122064775del']['hgvs_t_and_p']['NM_032790.3']['transcript_variant_error'] is None

    def test_variant23(self):
        variant = 'NC_000012.11:g.122064772_122064775dup'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000012.11:g.122064772_122064775dup' in results.keys()
        assert results['NC_000012.11:g.122064772_122064775dup']['p_vcf'] == '12:122064771:G:GCCCC'
        assert results['NC_000012.11:g.122064772_122064775dup']['g_hgvs'] == 'NC_000012.11:g.122064772_122064775dup'
        assert results['NC_000012.11:g.122064772_122064775dup']['genomic_variant_error'] is None
        assert 'NM_032790.3' in results['NC_000012.11:g.122064772_122064775dup']['hgvs_t_and_p'].keys()
        assert results['NC_000012.11:g.122064772_122064775dup']['hgvs_t_and_p']['NM_032790.3']['t_hgvs'] == "NM_032790.3:c.128_129insCCCCGCCACC"
        assert results['NC_000012.11:g.122064772_122064775dup']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_tlc'] == "NP_116179.2:p.(Pro45AlafsTer46)"
        assert results['NC_000012.11:g.122064772_122064775dup']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_slc'] == "NP_116179.2:p.(P45Afs*46)"
        assert results['NC_000012.11:g.122064772_122064775dup']['hgvs_t_and_p']['NM_032790.3']['transcript_variant_error'] is None

    def test_variant24(self):
        variant = 'NC_000012.11:g.122064773_122064774insTTTT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000012.11:g.122064773_122064774insTTTT' in results.keys()
        assert results['NC_000012.11:g.122064773_122064774insTTTT']['p_vcf'] == '12:122064773:C:CTTTT'
        assert results['NC_000012.11:g.122064773_122064774insTTTT']['g_hgvs'] == 'NC_000012.11:g.122064773_122064774insTTTT'
        assert results['NC_000012.11:g.122064773_122064774insTTTT']['genomic_variant_error'] is None
        assert 'NM_032790.3' in results['NC_000012.11:g.122064773_122064774insTTTT']['hgvs_t_and_p'].keys()
        assert results['NC_000012.11:g.122064773_122064774insTTTT']['hgvs_t_and_p']['NM_032790.3']['t_hgvs'] == "NM_032790.3:c.126_127insTTTTCCGCCA"
        assert results['NC_000012.11:g.122064773_122064774insTTTT']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_tlc'] == "NP_116179.2:p.(Pro43PhefsTer48)"
        assert results['NC_000012.11:g.122064773_122064774insTTTT']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_slc'] == "NP_116179.2:p.(P43Ffs*48)"
        assert results['NC_000012.11:g.122064773_122064774insTTTT']['hgvs_t_and_p']['NM_032790.3']['transcript_variant_error'] is None

    def test_variant25(self):
        variant = 'NC_000012.11:g.122064772_122064777del'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000012.11:g.122064772_122064777del' in results.keys()
        assert results['NC_000012.11:g.122064772_122064777del']['p_vcf'] == '12:122064771:GCCCCGC:G'
        assert results['NC_000012.11:g.122064772_122064777del']['g_hgvs'] == 'NC_000012.11:g.122064773_122064778del'
        assert results['NC_000012.11:g.122064772_122064777del']['genomic_variant_error'] == "NC_000012.11:g.122064772_122064777del updated to NC_000012.11:g.122064773_122064778del"
        assert 'NM_032790.3' in results['NC_000012.11:g.122064772_122064777del']['hgvs_t_and_p'].keys()
        assert results['NC_000012.11:g.122064772_122064777del']['hgvs_t_and_p']['NM_032790.3']['t_hgvs'] == "NM_032790.3:c.126C>A"
        assert results['NC_000012.11:g.122064772_122064777del']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_tlc'] == "NP_116179.2:p.(Ala42=)"
        assert results['NC_000012.11:g.122064772_122064777del']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_slc'] == "NP_116179.2:p.(A42=)"
        assert results['NC_000012.11:g.122064772_122064777del']['hgvs_t_and_p']['NM_032790.3']['transcript_variant_error'] is None

    def test_variant26(self):
        variant = 'NC_000012.11:g.122064772_122064777dup'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000012.11:g.122064772_122064777dup' in results.keys()
        assert results['NC_000012.11:g.122064772_122064777dup']['p_vcf'] == '12:122064771:G:GCCCCGC'
        assert results['NC_000012.11:g.122064772_122064777dup']['g_hgvs'] == 'NC_000012.11:g.122064773_122064778dup'
        assert results['NC_000012.11:g.122064772_122064777dup']['genomic_variant_error'] == "NC_000012.11:g.122064772_122064777dup updated to NC_000012.11:g.122064773_122064778dup"
        assert 'NM_032790.3' in results['NC_000012.11:g.122064772_122064777dup']['hgvs_t_and_p'].keys()
        assert results['NC_000012.11:g.122064772_122064777dup']['hgvs_t_and_p']['NM_032790.3']['t_hgvs'] == "NM_032790.3:c.131_132insCCCGCCACCGCC"
        assert results['NC_000012.11:g.122064772_122064777dup']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_tlc'] == "NP_116179.2:p.(Pro44_Pro47dup)"
        assert results['NC_000012.11:g.122064772_122064777dup']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_slc'] == "NP_116179.2:p.(P44_P47dup)"
        assert results['NC_000012.11:g.122064772_122064777dup']['hgvs_t_and_p']['NM_032790.3']['transcript_variant_error'] is None

    def test_variant27(self):
        variant = 'NC_000012.11:g.122064779_122064782dup'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000012.11:g.122064779_122064782dup' in results.keys()
        assert results['NC_000012.11:g.122064779_122064782dup']['p_vcf'] == '12:122064778:C:CACCG'
        assert results['NC_000012.11:g.122064779_122064782dup']['g_hgvs'] == 'NC_000012.11:g.122064779_122064782dup'
        assert results['NC_000012.11:g.122064779_122064782dup']['genomic_variant_error'] is None
        assert 'NM_032790.3' in results['NC_000012.11:g.122064779_122064782dup']['hgvs_t_and_p'].keys()
        assert results['NC_000012.11:g.122064779_122064782dup']['hgvs_t_and_p']['NM_032790.3']['t_hgvs'] == "NM_032790.3:c.135_136insACCGCCACCG"
        assert results['NC_000012.11:g.122064779_122064782dup']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_tlc'] == "NP_116179.2:p.(Pro46ThrfsTer45)"
        assert results['NC_000012.11:g.122064779_122064782dup']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_slc'] == "NP_116179.2:p.(P46Tfs*45)"
        assert results['NC_000012.11:g.122064779_122064782dup']['hgvs_t_and_p']['NM_032790.3']['transcript_variant_error'] is None

    def test_variant28(self):
        variant = 'NC_000012.11:g.122064772_122064782del'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000012.11:g.122064772_122064782del' in results.keys()
        assert results['NC_000012.11:g.122064772_122064782del']['p_vcf'] == '12:122064770:GGCCCCGCCACC:G'
        assert results['NC_000012.11:g.122064772_122064782del']['g_hgvs'] == 'NC_000012.11:g.122064774_122064784del'
        assert results['NC_000012.11:g.122064772_122064782del']['genomic_variant_error'] == "NC_000012.11:g.122064772_122064782del updated to NC_000012.11:g.122064774_122064784del"
        assert 'NM_032790.3' in results['NC_000012.11:g.122064772_122064782del']['hgvs_t_and_p'].keys()
        assert results['NC_000012.11:g.122064772_122064782del']['hgvs_t_and_p']['NM_032790.3']['t_hgvs'] == "NM_032790.3:c.126_127insA"
        assert results['NC_000012.11:g.122064772_122064782del']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_tlc'] == "NP_116179.2:p.(Pro43ThrfsTer45)"
        assert results['NC_000012.11:g.122064772_122064782del']['hgvs_t_and_p']['NM_032790.3']['p_hgvs_slc'] == "NP_116179.2:p.(P43Tfs*45)"
        assert results['NC_000012.11:g.122064772_122064782del']['hgvs_t_and_p']['NM_032790.3']['transcript_variant_error'] is None

    def test_variant29(self):
        variant = 'NC_000002.11:g.95847041_95847043GCG='
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000002.11:g.95847041_95847043GCG=' in results.keys()
        assert results['NC_000002.11:g.95847041_95847043GCG=']['p_vcf'] == '2:95847041:GCG:GCG'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['g_hgvs'] == 'NC_000002.11:g.95847041_95847043='
        assert results['NC_000002.11:g.95847041_95847043GCG=']['genomic_variant_error'] is None
        assert 'NM_021088.2' in results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'].keys()
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_021088.2']['t_hgvs'] == 'NM_021088.2:c.471_473dup'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_021088.2']['p_hgvs_tlc'] == 'NP_066574.2:p.(Arg159dup)'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_021088.2']['p_hgvs_slc'] == 'NP_066574.2:p.(R159dup)'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_021088.2']['transcript_variant_error'] is None
        assert 'NM_021088.3' in results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'].keys()
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_021088.3']['t_hgvs'] == 'NM_021088.3:c.471_473dup'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_021088.3']['p_hgvs_tlc'] == 'NP_066574.2:p.(Arg159dup)'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_021088.3']['p_hgvs_slc'] == 'NP_066574.2:p.(R159dup)'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_021088.3']['transcript_variant_error'] is None
        assert 'NM_001291605.1' in results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'].keys()
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001291605.1']['t_hgvs'] == 'NM_001291605.1:c.510_512dup'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001291605.1']['p_hgvs_tlc'] == 'NP_001278534.1:p.(Arg172dup)'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001291605.1']['p_hgvs_slc'] == 'NP_001278534.1:p.(R172dup)'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001291605.1']['transcript_variant_error'] is None
        assert 'NM_001291604.1' in results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'].keys()
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001291604.1']['t_hgvs'] == 'NM_001291604.1:c.231_233dup'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001291604.1']['p_hgvs_tlc'] == 'NP_001278533.1:p.(Arg79dup)'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001291604.1']['p_hgvs_slc'] == 'NP_001278533.1:p.(R79dup)'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001291604.1']['transcript_variant_error'] is None
        assert 'NM_001017396.1' in results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'].keys()
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001017396.1']['t_hgvs'] == 'NM_001017396.1:c.345_347dup'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001017396.1']['p_hgvs_tlc'] == 'NP_001017396.1:p.(Arg117dup)'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001017396.1']['p_hgvs_slc'] == 'NP_001017396.1:p.(R117dup)'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001017396.1']['transcript_variant_error'] is None
        assert 'NM_001017396.2' in results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'].keys()
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001017396.2']['t_hgvs'] == 'NM_001017396.2:c.345_347dup'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001017396.2']['p_hgvs_tlc'] == 'NP_001017396.1:p.(Arg117dup)'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001017396.2']['p_hgvs_slc'] == 'NP_001017396.1:p.(R117dup)'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001017396.2']['transcript_variant_error'] is None
        assert 'NM_001282398.1' in results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'].keys()
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001282398.1']['t_hgvs'] == 'NM_001282398.1:c.357_359dup'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001282398.1']['p_hgvs_tlc'] == 'NP_001269327.1:p.(Arg121dup)'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001282398.1']['p_hgvs_slc'] == 'NP_001269327.1:p.(R121dup)'
        assert results['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p']['NM_001282398.1']['transcript_variant_error'] is None

    def test_variant30(self):
        variant = 'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG='
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=' in results.keys()
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['p_vcf'] == '17:5286863:AGTGTTTGGAATTTTCTGTTCATATAG:AGTGTTTGGAATTTTCTGTTCATATAG'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['g_hgvs'] == 'NC_000017.10:g.5286863_5286889='
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['genomic_variant_error'] is None
        assert 'NM_001291581.1' in results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'].keys()
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_001291581.1']['t_hgvs'] == 'NM_001291581.1:c.*344_*368dup'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_001291581.1']['p_hgvs_tlc'] == 'NP_001278510.1:p.?'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_001291581.1']['p_hgvs_slc'] == 'NP_001278510.1:p.?'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_001291581.1']['transcript_variant_error'] is None
        assert 'NM_004703.4' in results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'].keys()
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_004703.4']['t_hgvs'] == 'NM_004703.4:c.*344_*368dup'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_004703.4']['p_hgvs_tlc'] == 'NP_004694.2:p.?'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_004703.4']['p_hgvs_slc'] == 'NP_004694.2:p.?'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_004703.4']['transcript_variant_error'] is None
        assert 'NM_001083585.2' in results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'].keys()
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_001083585.2']['t_hgvs'] == 'NM_001083585.2:c.*344_*368dup'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_001083585.2']['p_hgvs_tlc'] == 'NP_001077054.1:p.?'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_001083585.2']['p_hgvs_slc'] == 'NP_001077054.1:p.?'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_001083585.2']['transcript_variant_error'] is None
        assert 'NM_004703.5' in results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'].keys()
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_004703.5']['t_hgvs'] == 'NM_004703.5:c.*344_*368dup'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_004703.5']['p_hgvs_tlc'] == 'NP_004694.2:p.?'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_004703.5']['p_hgvs_slc'] == 'NP_004694.2:p.?'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_004703.5']['transcript_variant_error'] is None
        assert 'NM_001083585.1' in results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'].keys()
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_001083585.1']['t_hgvs'] == 'NM_001083585.1:c.*344_*368dup'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_001083585.1']['p_hgvs_tlc'] == 'NP_001077054.1:p.?'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_001083585.1']['p_hgvs_slc'] == 'NP_001077054.1:p.?'
        assert results['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p']['NM_001083585.1']['transcript_variant_error'] is None

    def test_variant31(self):
        variant = 'NC_000003.11:g.14561629_14561630GC='
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000003.11:g.14561629_14561630GC=' in results.keys()
        assert results['NC_000003.11:g.14561629_14561630GC=']['p_vcf'] == '3:14561629:GC:GC'
        assert results['NC_000003.11:g.14561629_14561630GC=']['g_hgvs'] == 'NC_000003.11:g.14561629_14561630='
        assert results['NC_000003.11:g.14561629_14561630GC=']['genomic_variant_error'] is None
        assert 'NM_001080423.2' in results['NC_000003.11:g.14561629_14561630GC=']['hgvs_t_and_p'].keys()
        assert results['NC_000003.11:g.14561629_14561630GC=']['hgvs_t_and_p']['NM_001080423.2']['t_hgvs'] == 'NM_001080423.2:c.1311del'
        assert results['NC_000003.11:g.14561629_14561630GC=']['hgvs_t_and_p']['NM_001080423.2']['p_hgvs_tlc'] == 'NP_001073892.2:p.(Ser438GlnfsTer4)'
        assert results['NC_000003.11:g.14561629_14561630GC=']['hgvs_t_and_p']['NM_001080423.2']['p_hgvs_slc'] == 'NP_001073892.2:p.(S438Qfs*4)'
        assert results['NC_000003.11:g.14561629_14561630GC=']['hgvs_t_and_p']['NM_001080423.2']['transcript_variant_error'] is None
        assert 'NM_001080423.3' in results['NC_000003.11:g.14561629_14561630GC=']['hgvs_t_and_p'].keys()
        assert results['NC_000003.11:g.14561629_14561630GC=']['hgvs_t_and_p']['NM_001080423.3']['t_hgvs'] == 'NM_001080423.3:c.1020del'
        assert results['NC_000003.11:g.14561629_14561630GC=']['hgvs_t_and_p']['NM_001080423.3']['p_hgvs_tlc'] == 'NP_001073892.3:p.(Ser341GlnfsTer4)'
        assert results['NC_000003.11:g.14561629_14561630GC=']['hgvs_t_and_p']['NM_001080423.3']['p_hgvs_slc'] == 'NP_001073892.3:p.(S341Qfs*4)'
        assert results['NC_000003.11:g.14561629_14561630GC=']['hgvs_t_and_p']['NM_001080423.3']['transcript_variant_error'] is None

    def test_variant32(self):
        variant = 'NC_000003.11:g.14561629_14561630insG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000003.11:g.14561629_14561630insG' in results.keys()
        assert results['NC_000003.11:g.14561629_14561630insG']['p_vcf'] == '3:14561627:A:AG'
        assert results['NC_000003.11:g.14561629_14561630insG']['g_hgvs'] == 'NC_000003.11:g.14561629dup'
        assert results['NC_000003.11:g.14561629_14561630insG']['genomic_variant_error'] is None
        assert 'NM_001080423.2' in results['NC_000003.11:g.14561629_14561630insG']['hgvs_t_and_p'].keys()
        assert results['NC_000003.11:g.14561629_14561630insG']['hgvs_t_and_p']['NM_001080423.2']['t_hgvs'] == 'NM_001080423.2:c.1308_1311='
        assert results['NC_000003.11:g.14561629_14561630insG']['hgvs_t_and_p']['NM_001080423.2']['p_hgvs_tlc'] == 'NP_001073892.2:p.(Arg436_Pro437=)'
        assert results['NC_000003.11:g.14561629_14561630insG']['hgvs_t_and_p']['NM_001080423.2']['transcript_variant_error'] is None
        assert 'NM_001080423.3' in results['NC_000003.11:g.14561629_14561630insG']['hgvs_t_and_p'].keys()
        assert results['NC_000003.11:g.14561629_14561630insG']['hgvs_t_and_p']['NM_001080423.3']['t_hgvs'] == 'NM_001080423.3:c.1017_1020='
        assert results['NC_000003.11:g.14561629_14561630insG']['hgvs_t_and_p']['NM_001080423.3']['p_hgvs_tlc'] == 'NP_001073892.3:p.(Arg339_Pro340=)'
        assert results['NC_000003.11:g.14561629_14561630insG']['hgvs_t_and_p']['NM_001080423.3']['transcript_variant_error'] is None

    def test_variant33(self):
        variant = 'NC_000004.11:g.140811111_140811122del'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000004.11:g.140811111_140811122del' in results.keys()
        assert results['NC_000004.11:g.140811111_140811122del']['p_vcf'] == '4:140811063:TTGCTGCTGCTGC:T'
        assert results['NC_000004.11:g.140811111_140811122del']['g_hgvs'] == 'NC_000004.11:g.140811111_140811122del'
        assert results['NC_000004.11:g.140811111_140811122del']['genomic_variant_error'] is None
        assert 'NM_018717.5' in results['NC_000004.11:g.140811111_140811122del']['hgvs_t_and_p'].keys()
        assert results['NC_000004.11:g.140811111_140811122del']['hgvs_t_and_p']['NM_018717.5']['t_hgvs'] == 'NM_018717.5:c.1515_1526del'
        assert results['NC_000004.11:g.140811111_140811122del']['hgvs_t_and_p']['NM_018717.5']['p_hgvs_tlc'] == 'NP_061187.3:p.(Gln507_Gln510del)'
        assert results['NC_000004.11:g.140811111_140811122del']['hgvs_t_and_p']['NM_018717.5']['p_hgvs_slc'] == 'NP_061187.3:p.(Q507_Q510del)'
        assert results['NC_000004.11:g.140811111_140811122del']['hgvs_t_and_p']['NM_018717.5']['transcript_variant_error'] is None
        assert 'NM_018717.4' in results['NC_000004.11:g.140811111_140811122del']['hgvs_t_and_p'].keys()
        assert results['NC_000004.11:g.140811111_140811122del']['hgvs_t_and_p']['NM_018717.4']['t_hgvs'] == 'NM_018717.4:c.1466_1468='
        assert results['NC_000004.11:g.140811111_140811122del']['hgvs_t_and_p']['NM_018717.4']['p_hgvs_tlc'] == 'NP_061187.2:p.(Gln489_Gln490=)'
        assert results['NC_000004.11:g.140811111_140811122del']['hgvs_t_and_p']['NM_018717.4']['transcript_variant_error'] is None

    def test_variant34(self):
        variant = 'NC_000004.11:g.140811111_140811122CTGCTGCTGCTG='
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=' in results.keys()
        assert results['NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=']['p_vcf'] == '4:140811111:CTGCTGCTGCTG:CTGCTGCTGCTG'
        assert results['NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=']['g_hgvs'] == 'NC_000004.11:g.140811111_140811122='
        assert results['NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=']['genomic_variant_error'] is None
        assert 'NM_018717.5' in results['NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=']['hgvs_t_and_p'].keys()
        assert results['NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=']['hgvs_t_and_p']['NM_018717.5']['t_hgvs'] == 'NM_018717.5:c.1468_1479='
        assert results['NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=']['hgvs_t_and_p']['NM_018717.5']['p_hgvs_tlc'] == 'NP_061187.3:p.(Gln490_Gln493=)'
        assert results['NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=']['hgvs_t_and_p']['NM_018717.5']['transcript_variant_error'] is None
        assert 'NM_018717.4' in results['NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=']['hgvs_t_and_p'].keys()
        assert results['NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=']['hgvs_t_and_p']['NM_018717.4']['t_hgvs'] == 'NM_018717.4:c.1503_1514dup'
        assert results['NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=']['hgvs_t_and_p']['NM_018717.4']['p_hgvs_tlc'] == 'NP_061187.2:p.(Gln503_Gln506dup)'
        assert results['NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=']['hgvs_t_and_p']['NM_018717.4']['p_hgvs_slc'] == 'NP_061187.2:p.(Q503_Q506dup)'
        assert results['NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=']['hgvs_t_and_p']['NM_018717.4']['transcript_variant_error'] is None

    def test_variant35(self):
        variant = 'NC_000004.11:g.140811117_140811122del'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000004.11:g.140811117_140811122del' in results.keys()
        assert results['NC_000004.11:g.140811117_140811122del']['p_vcf'] == '4:140811063:TTGCTGC:T'
        assert results['NC_000004.11:g.140811117_140811122del']['g_hgvs'] == 'NC_000004.11:g.140811117_140811122del'
        assert results['NC_000004.11:g.140811117_140811122del']['genomic_variant_error'] is None
        assert 'NM_018717.5' in results['NC_000004.11:g.140811117_140811122del']['hgvs_t_and_p'].keys()
        assert results['NC_000004.11:g.140811117_140811122del']['hgvs_t_and_p']['NM_018717.5']['t_hgvs'] == 'NM_018717.5:c.1521_1526del'
        assert results['NC_000004.11:g.140811117_140811122del']['hgvs_t_and_p']['NM_018717.5']['p_hgvs_tlc'] == 'NP_061187.3:p.(Gln509_Gln510del)'
        assert results['NC_000004.11:g.140811117_140811122del']['hgvs_t_and_p']['NM_018717.5']['p_hgvs_slc'] == 'NP_061187.3:p.(Q509_Q510del)'
        assert results['NC_000004.11:g.140811117_140811122del']['hgvs_t_and_p']['NM_018717.5']['transcript_variant_error'] is None
        assert 'NM_018717.4' in results['NC_000004.11:g.140811117_140811122del']['hgvs_t_and_p'].keys()
        assert results['NC_000004.11:g.140811117_140811122del']['hgvs_t_and_p']['NM_018717.4']['t_hgvs'] == 'NM_018717.4:c.1509_1514dup'
        assert results['NC_000004.11:g.140811117_140811122del']['hgvs_t_and_p']['NM_018717.4']['p_hgvs_tlc'] == 'NP_061187.2:p.(Gln505_Gln506dup)'
        assert results['NC_000004.11:g.140811117_140811122del']['hgvs_t_and_p']['NM_018717.4']['p_hgvs_slc'] == 'NP_061187.2:p.(Q505_Q506dup)'
        assert results['NC_000004.11:g.140811117_140811122del']['hgvs_t_and_p']['NM_018717.4']['transcript_variant_error'] is None

    def test_variant36(self):
        variant = 'NC_000004.11:g.140811111_140811117del'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000004.11:g.140811111_140811117del' in results.keys()
        assert results['NC_000004.11:g.140811111_140811117del']['p_vcf'] == '4:140811110:GCTGCTGC:G'
        assert results['NC_000004.11:g.140811111_140811117del']['g_hgvs'] == 'NC_000004.11:g.140811111_140811117del'
        assert results['NC_000004.11:g.140811111_140811117del']['genomic_variant_error'] is None
        assert 'NM_018717.5' in results['NC_000004.11:g.140811111_140811117del']['hgvs_t_and_p'].keys()
        assert results['NC_000004.11:g.140811111_140811117del']['hgvs_t_and_p']['NM_018717.5']['t_hgvs'] == 'NM_018717.5:c.1473_1479del'
        assert results['NC_000004.11:g.140811111_140811117del']['hgvs_t_and_p']['NM_018717.5']['p_hgvs_tlc'] == 'NP_061187.3:p.(Gln491HisfsTer29)'
        assert results['NC_000004.11:g.140811111_140811117del']['hgvs_t_and_p']['NM_018717.5']['p_hgvs_slc'] == 'NP_061187.3:p.(Q491Hfs*29)'
        assert results['NC_000004.11:g.140811111_140811117del']['hgvs_t_and_p']['NM_018717.5']['transcript_variant_error'] is None
        assert 'NM_018717.4' in results['NC_000004.11:g.140811111_140811117del']['hgvs_t_and_p'].keys()
        assert results['NC_000004.11:g.140811111_140811117del']['hgvs_t_and_p']['NM_018717.4']['t_hgvs'] == "NM_018717.4:c.1468_1472dup"
        assert results['NC_000004.11:g.140811111_140811117del']['hgvs_t_and_p']['NM_018717.4']['p_hgvs_tlc'] == "NP_061187.2:p.(Gln491HisfsTer29)"
        assert results['NC_000004.11:g.140811111_140811117del']['hgvs_t_and_p']['NM_018717.4']['p_hgvs_slc'] == "NP_061187.2:p.(Q491Hfs*29)"
        assert results['NC_000004.11:g.140811111_140811117del']['hgvs_t_and_p']['NM_018717.4']['transcript_variant_error'] is None

    def test_variant37(self):
        variant = 'NC_000004.11:g.140811117C>A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000004.11:g.140811117C>A' in results.keys()
        assert results['NC_000004.11:g.140811117C>A']['p_vcf'] == '4:140811117:C:A'
        assert results['NC_000004.11:g.140811117C>A']['g_hgvs'] == 'NC_000004.11:g.140811117C>A'
        assert results['NC_000004.11:g.140811117C>A']['genomic_variant_error'] is None
        assert 'NM_018717.5' in results['NC_000004.11:g.140811117C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000004.11:g.140811117C>A']['hgvs_t_and_p']['NM_018717.5']['t_hgvs'] == 'NM_018717.5:c.1473G>T'
        assert results['NC_000004.11:g.140811117C>A']['hgvs_t_and_p']['NM_018717.5']['p_hgvs_tlc'] == 'NP_061187.3:p.(Gln491His)'
        assert results['NC_000004.11:g.140811117C>A']['hgvs_t_and_p']['NM_018717.5']['p_hgvs_slc'] == 'NP_061187.3:p.(Q491H)'
        assert results['NC_000004.11:g.140811117C>A']['hgvs_t_and_p']['NM_018717.5']['transcript_variant_error'] is None
        assert 'NM_018717.4' in results['NC_000004.11:g.140811117C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000004.11:g.140811117C>A']['hgvs_t_and_p']['NM_018717.4']['t_hgvs'] == "NM_018717.4:c.1472_1473insTCAGCAGCAGCA"
        assert results['NC_000004.11:g.140811117C>A']['hgvs_t_and_p']['NM_018717.4']['p_hgvs_tlc'] == "NP_061187.2:p.(Gln490_Gln491insHisGlnGlnGln)"
        assert results['NC_000004.11:g.140811117C>A']['hgvs_t_and_p']['NM_018717.4']['p_hgvs_slc'] == "NP_061187.2:p.(Q490_Q491insHQQQ)"
        assert results['NC_000004.11:g.140811117C>A']['hgvs_t_and_p']['NM_018717.4']['transcript_variant_error'] is None

    def test_variant38(self):
        variant = 'NC_000002.11:g.73675227_73675228insCTC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000002.11:g.73675227_73675228insCTC' in results.keys()
        assert results['NC_000002.11:g.73675227_73675228insCTC']['p_vcf'] == '2:73675227:T:TCTC'
        assert results['NC_000002.11:g.73675227_73675228insCTC']['g_hgvs'] == 'NC_000002.11:g.73675228_73675230dup'
        assert results['NC_000002.11:g.73675227_73675228insCTC']['genomic_variant_error'] == "NC_000002.11:g.73675227_73675228insCTC updated to NC_000002.11:g.73675228_73675230dup"
        assert 'NM_015120.4' in results['NC_000002.11:g.73675227_73675228insCTC']['hgvs_t_and_p'].keys()
        assert results['NC_000002.11:g.73675227_73675228insCTC']['hgvs_t_and_p']['NM_015120.4']['t_hgvs'] == 'NM_015120.4:c.1573_1579='
        assert results['NC_000002.11:g.73675227_73675228insCTC']['hgvs_t_and_p']['NM_015120.4']['p_hgvs_tlc'] == 'NP_055935.4:p.(Ser525_Leu527=)'
        assert results['NC_000002.11:g.73675227_73675228insCTC']['hgvs_t_and_p']['NM_015120.4']['transcript_variant_error'] is None

    def test_variant39(self):
        variant = '9-136132908-T-TC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '9-136132908-T-TC' in results.keys()
        assert results['9-136132908-T-TC']['p_vcf'] == '9-136132908-T-TC'
        assert results['9-136132908-T-TC']['g_hgvs'] == 'NC_000009.11:g.136132908_136132909insC'
        assert results['9-136132908-T-TC']['genomic_variant_error'] is None
        assert 'NM_020469.2' in results['9-136132908-T-TC']['hgvs_t_and_p'].keys()
        assert results['9-136132908-T-TC']['hgvs_t_and_p']['NM_020469.2']['t_hgvs'] == 'NM_020469.2:c.260_262='
        assert results['9-136132908-T-TC']['hgvs_t_and_p']['NM_020469.2']['p_hgvs_tlc'] == 'NP_065202.2:p.(Val87_Thr88=)'
        assert results['9-136132908-T-TC']['hgvs_t_and_p']['NM_020469.2']['transcript_variant_error'] is None

    def test_variant40(self):
        variant = '9-136132908-TAC-TCA'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '9-136132908-TAC-TCA' in results.keys()
        assert results['9-136132908-TAC-TCA']['p_vcf'] == '9-136132909-AC-CA'
        assert results['9-136132908-TAC-TCA']['g_hgvs'] == 'NC_000009.11:g.136132909_136132910delinsCA'
        assert results['9-136132908-TAC-TCA']['genomic_variant_error'] is None
        assert 'NM_020469.2' in results['9-136132908-TAC-TCA']['hgvs_t_and_p'].keys()
        assert results['9-136132908-TAC-TCA']['hgvs_t_and_p']['NM_020469.2']['t_hgvs'] == 'NM_020469.2:c.259del'
        assert results['9-136132908-TAC-TCA']['hgvs_t_and_p']['NM_020469.2']['p_hgvs_tlc'] == 'NP_065202.2:p.(Val87Ter)'
        assert results['9-136132908-TAC-TCA']['hgvs_t_and_p']['NM_020469.2']['p_hgvs_slc'] == 'NP_065202.2:p.(V87*)'
        assert results['9-136132908-TAC-TCA']['hgvs_t_and_p']['NM_020469.2']['transcript_variant_error'] is None

    def test_variant41(self):
        variant = '9-136132908-TA-TA'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '9-136132908-TA-TA' in results.keys()
        assert results['9-136132908-TA-TA']['p_vcf'] == '9-136132908-TA-TA'
        assert results['9-136132908-TA-TA']['g_hgvs'] == 'NC_000009.11:g.136132908_136132909='
        assert results['9-136132908-TA-TA']['genomic_variant_error'] is None
        assert 'NM_020469.2' in results['9-136132908-TA-TA']['hgvs_t_and_p'].keys()
        assert results['9-136132908-TA-TA']['hgvs_t_and_p']['NM_020469.2']['t_hgvs'] == 'NM_020469.2:c.261del'
        assert results['9-136132908-TA-TA']['hgvs_t_and_p']['NM_020469.2']['p_hgvs_tlc'] == 'NP_065202.2:p.(Thr88ProfsTer31)'
        assert results['9-136132908-TA-TA']['hgvs_t_and_p']['NM_020469.2']['p_hgvs_slc'] == 'NP_065202.2:p.(T88Pfs*31)'
        assert results['9-136132908-TA-TA']['hgvs_t_and_p']['NM_020469.2']['transcript_variant_error'] is None

    def test_variant42(self):
        variant = 'NC_012920.1:m.1011C>T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_012920.1:m.1011C>T' in results.keys()
        assert results['NC_012920.1:m.1011C>T']['p_vcf'] == 'M:1011:C:T'
        assert results['NC_012920.1:m.1011C>T']['g_hgvs'] == 'NC_012920.1:m.1011C>T'
        assert results['NC_012920.1:m.1011C>T']['hgvs_t_and_p'] == {'intergenic': {'alt_genomic_loci': None}}

    def test_variant43(self):
        variant = 'NC_000006.11:g.90403795G='
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000006.11:g.90403795G=' in results.keys()
        assert results['NC_000006.11:g.90403795G=']['p_vcf'] == '6:90403795:G:G'
        assert results['NC_000006.11:g.90403795G=']['g_hgvs'] == 'NC_000006.11:g.90403795='
        assert results['NC_000006.11:g.90403795G=']['genomic_variant_error'] is None
        assert 'NM_014611.2' in results['NC_000006.11:g.90403795G=']['hgvs_t_and_p'].keys()
        assert results['NC_000006.11:g.90403795G=']['hgvs_t_and_p']['NM_014611.2']['t_hgvs'] == 'NM_014611.2:c.9879='
        assert results['NC_000006.11:g.90403795G=']['hgvs_t_and_p']['NM_014611.2']['p_hgvs_tlc'] == 'NP_055426.1:p.(Val3293=)'
        assert results['NC_000006.11:g.90403795G=']['hgvs_t_and_p']['NM_014611.2']['p_hgvs_slc'] == 'NP_055426.1:p.(V3293=)'
        assert results['NC_000006.11:g.90403795G=']['hgvs_t_and_p']['NM_014611.2']['transcript_variant_error'] is None
        assert 'NM_014611.1' in results['NC_000006.11:g.90403795G=']['hgvs_t_and_p'].keys()
        assert results['NC_000006.11:g.90403795G=']['hgvs_t_and_p']['NM_014611.1']['t_hgvs'] == 'NM_014611.1:c.9879T>C'
        assert results['NC_000006.11:g.90403795G=']['hgvs_t_and_p']['NM_014611.1']['p_hgvs_tlc'] == 'NP_055426.1:p.(Val3293=)'
        assert results['NC_000006.11:g.90403795G=']['hgvs_t_and_p']['NM_014611.1']['p_hgvs_slc'] == 'NP_055426.1:p.(V3293=)'
        assert results['NC_000006.11:g.90403795G=']['hgvs_t_and_p']['NM_014611.1']['transcript_variant_error'] is None

    def test_variant44(self):
        variant = 'NC_000005.9:g.35058667_35058668AG='
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000005.9:g.35058667_35058668AG=' in results.keys()
        assert results['NC_000005.9:g.35058667_35058668AG=']['p_vcf'] == '5:35058667:AG:AG'
        assert results['NC_000005.9:g.35058667_35058668AG=']['g_hgvs'] == 'NC_000005.9:g.35058667_35058668='
        assert results['NC_000005.9:g.35058667_35058668AG=']['genomic_variant_error'] is None
        assert 'NM_001204314.1' in results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p'].keys()
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204314.1']['t_hgvs'] == 'NM_001204314.1:c.*6523_*6524='
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204314.1']['p_hgvs_tlc'] == 'NP_001191243.1:p.?'
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204314.1']['p_hgvs_slc'] == 'NP_001191243.1:p.?'
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204314.1']['transcript_variant_error'] is None
        assert 'NM_001204314.2' in results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p'].keys()
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204314.2']['t_hgvs'] == 'NM_001204314.2:c.*6528del'
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204314.2']['p_hgvs_tlc'] == 'NP_001191243.1:p.?'
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204314.2']['p_hgvs_slc'] == 'NP_001191243.1:p.?'
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204314.2']['transcript_variant_error'] is None
        assert 'NM_000949.5' in results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p'].keys()
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_000949.5']['t_hgvs'] == 'NM_000949.5:c.*6523_*6524='
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_000949.5']['p_hgvs_tlc'] == 'NP_000940.1:p.?'
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_000949.5']['p_hgvs_slc'] == 'NP_000940.1:p.?'
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_000949.5']['transcript_variant_error'] is None
        assert 'NM_001204317.1' in results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p'].keys()
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204317.1']['t_hgvs'] == 'NM_001204317.1:c.856-9155_856-9154='
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204317.1']['p_hgvs_tlc'] == 'NP_001191246.1:p.?'
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204317.1']['p_hgvs_slc'] == 'NP_001191246.1:p.?'
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204317.1']['transcript_variant_error'] is None
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204316.1']['t_hgvs'] == 'NM_001204316.1:c.1009+7383_1009+7384='
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204316.1']['p_hgvs_tlc'] == 'NP_001191245.1:p.?'
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204316.1']['p_hgvs_slc'] == 'NP_001191245.1:p.?'
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204316.1']['transcript_variant_error'] is None
        assert 'NR_037910.1' in results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p'].keys()
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NR_037910.1']['t_hgvs'] == 'NR_037910.1:n.828-9155_828-9154='
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NR_037910.1']['p_hgvs_tlc'] is None
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NR_037910.1']['p_hgvs_slc'] is None
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NR_037910.1']['transcript_variant_error'] is None
        assert 'NM_001204318.1' in results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p'].keys()
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204318.1']['t_hgvs'] == 'NM_001204318.1:c.686-9155_686-9154='
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204318.1']['p_hgvs_tlc'] == 'NP_001191247.1:p.?'
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204318.1']['p_hgvs_slc'] == 'NP_001191247.1:p.?'
        assert results['NC_000005.9:g.35058667_35058668AG=']['hgvs_t_and_p']['NM_001204318.1']['transcript_variant_error'] is None

    def test_variant45(self):
        variant = '2-73675227-TCTC-TCTCCTC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '2-73675227-TCTC-TCTCCTC' in results.keys()
        assert results['2-73675227-TCTC-TCTCCTC']['p_vcf'] == '2-73675227-T-TCTC'
        assert results['2-73675227-TCTC-TCTCCTC']['g_hgvs'] == 'NC_000002.11:g.73675228_73675230dup'
        assert results['2-73675227-TCTC-TCTCCTC']['genomic_variant_error'] is None
        assert 'NM_015120.4' in results['2-73675227-TCTC-TCTCCTC']['hgvs_t_and_p'].keys()
        assert results['2-73675227-TCTC-TCTCCTC']['hgvs_t_and_p']['NM_015120.4']['t_hgvs'] == 'NM_015120.4:c.1573_1579='
        assert results['2-73675227-TCTC-TCTCCTC']['hgvs_t_and_p']['NM_015120.4']['p_hgvs_tlc'] == 'NP_055935.4:p.(Ser525_Leu527=)'
        assert results['2-73675227-TCTC-TCTCCTC']['hgvs_t_and_p']['NM_015120.4']['transcript_variant_error'] is None

    def test_variant46(self):
        variant = '2-73675227-TC-TC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '2-73675227-TC-TC' in results.keys()
        assert results['2-73675227-TC-TC']['p_vcf'] == '2-73675227-TC-TC'
        assert results['2-73675227-TC-TC']['g_hgvs'] == 'NC_000002.11:g.73675227_73675228='
        assert results['2-73675227-TC-TC']['genomic_variant_error'] is None
        assert 'NM_015120.4' in results['2-73675227-TC-TC']['hgvs_t_and_p'].keys()
        assert results['2-73675227-TC-TC']['hgvs_t_and_p']['NM_015120.4']['t_hgvs'] == 'NM_015120.4:c.1577_1579del'
        assert results['2-73675227-TC-TC']['hgvs_t_and_p']['NM_015120.4']['p_hgvs_tlc'] == 'NP_055935.4:p.(Pro526del)'
        assert results['2-73675227-TC-TC']['hgvs_t_and_p']['NM_015120.4']['p_hgvs_slc'] == 'NP_055935.4:p.(P526del)'
        assert results['2-73675227-TC-TC']['hgvs_t_and_p']['NM_015120.4']['transcript_variant_error'] is None

    def test_variant47(self):
        variant = '3-14561627-AG-AGG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '3-14561627-AG-AGG' in results.keys()
        assert results['3-14561627-AG-AGG']['p_vcf'] == '3-14561627-A-AG'
        assert results['3-14561627-AG-AGG']['g_hgvs'] == 'NC_000003.11:g.14561629dup'
        assert results['3-14561627-AG-AGG']['genomic_variant_error'] is None
        assert 'NM_001080423.2' in results['3-14561627-AG-AGG']['hgvs_t_and_p'].keys()
        assert results['3-14561627-AG-AGG']['hgvs_t_and_p']['NM_001080423.2']['t_hgvs'] == 'NM_001080423.2:c.1308_1311='
        assert results['3-14561627-AG-AGG']['hgvs_t_and_p']['NM_001080423.2']['p_hgvs_tlc'] == 'NP_001073892.2:p.(Arg436_Pro437=)'
        assert results['3-14561627-AG-AGG']['hgvs_t_and_p']['NM_001080423.2']['p_hgvs_slc'] == 'NP_001073892.2:p.(R436_P437=)'
        assert results['3-14561627-AG-AGG']['hgvs_t_and_p']['NM_001080423.2']['transcript_variant_error'] is None
        assert 'NM_001080423.3' in results['3-14561627-AG-AGG']['hgvs_t_and_p'].keys()
        assert results['3-14561627-AG-AGG']['hgvs_t_and_p']['NM_001080423.3']['t_hgvs'] == 'NM_001080423.3:c.1017_1020='
        assert results['3-14561627-AG-AGG']['hgvs_t_and_p']['NM_001080423.3']['p_hgvs_tlc'] == 'NP_001073892.3:p.(Arg339_Pro340=)'
        assert results['3-14561627-AG-AGG']['hgvs_t_and_p']['NM_001080423.3']['p_hgvs_slc'] == 'NP_001073892.3:p.(R339_P340=)'
        assert results['3-14561627-AG-AGG']['hgvs_t_and_p']['NM_001080423.3']['transcript_variant_error'] is None

    def test_variant48(self):
        variant = '3-14561630-CC-CC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '3-14561630-CC-CC' in results.keys()
        assert results['3-14561630-CC-CC']['p_vcf'] == '3-14561630-CC-CC'
        assert results['3-14561630-CC-CC']['g_hgvs'] == 'NC_000003.11:g.14561630_14561631='
        assert results['3-14561630-CC-CC']['genomic_variant_error'] is None
        assert 'NM_001080423.2' in results['3-14561630-CC-CC']['hgvs_t_and_p'].keys()
        assert results['3-14561630-CC-CC']['hgvs_t_and_p']['NM_001080423.2']['t_hgvs'] == 'NM_001080423.2:c.1311del'
        assert results['3-14561630-CC-CC']['hgvs_t_and_p']['NM_001080423.2']['p_hgvs_tlc'] == 'NP_001073892.2:p.(Ser438GlnfsTer4)'
        assert results['3-14561630-CC-CC']['hgvs_t_and_p']['NM_001080423.2']['p_hgvs_slc'] == 'NP_001073892.2:p.(S438Qfs*4)'
        assert results['3-14561630-CC-CC']['hgvs_t_and_p']['NM_001080423.2']['transcript_variant_error'] is None
        assert 'NM_001080423.3' in results['3-14561630-CC-CC']['hgvs_t_and_p'].keys()
        assert results['3-14561630-CC-CC']['hgvs_t_and_p']['NM_001080423.3']['t_hgvs'] == 'NM_001080423.3:c.1020del'
        assert results['3-14561630-CC-CC']['hgvs_t_and_p']['NM_001080423.3']['p_hgvs_tlc'] == 'NP_001073892.3:p.(Ser341GlnfsTer4)'
        assert results['3-14561630-CC-CC']['hgvs_t_and_p']['NM_001080423.3']['p_hgvs_slc'] == 'NP_001073892.3:p.(S341Qfs*4)'
        assert results['3-14561630-CC-CC']['hgvs_t_and_p']['NM_001080423.3']['transcript_variant_error'] is None

    def test_variant49(self):
        variant = '6-90403795-G-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '6-90403795-G-G' in results.keys()
        assert results['6-90403795-G-G']['p_vcf'] == '6-90403795-G-G'
        assert results['6-90403795-G-G']['g_hgvs'] == 'NC_000006.11:g.90403795='
        assert results['6-90403795-G-G']['genomic_variant_error'] is None
        assert 'NM_014611.2' in results['6-90403795-G-G']['hgvs_t_and_p'].keys()
        assert results['6-90403795-G-G']['hgvs_t_and_p']['NM_014611.2']['t_hgvs'] == 'NM_014611.2:c.9879='
        assert results['6-90403795-G-G']['hgvs_t_and_p']['NM_014611.2']['p_hgvs_tlc'] == 'NP_055426.1:p.(Val3293=)'
        assert results['6-90403795-G-G']['hgvs_t_and_p']['NM_014611.2']['p_hgvs_slc'] == 'NP_055426.1:p.(V3293=)'
        assert results['6-90403795-G-G']['hgvs_t_and_p']['NM_014611.2']['transcript_variant_error'] is None
        assert 'NM_014611.1' in results['6-90403795-G-G']['hgvs_t_and_p'].keys()
        assert results['6-90403795-G-G']['hgvs_t_and_p']['NM_014611.1']['t_hgvs'] == 'NM_014611.1:c.9879T>C'
        assert results['6-90403795-G-G']['hgvs_t_and_p']['NM_014611.1']['p_hgvs_tlc'] == 'NP_055426.1:p.(Val3293=)'
        assert results['6-90403795-G-G']['hgvs_t_and_p']['NM_014611.1']['p_hgvs_slc'] == 'NP_055426.1:p.(V3293=)'
        assert results['6-90403795-G-G']['hgvs_t_and_p']['NM_014611.1']['transcript_variant_error'] is None

    def test_variant50(self):
        variant = '6-90403795-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '6-90403795-G-A' in results.keys()
        assert results['6-90403795-G-A']['p_vcf'] == '6-90403795-G-A'
        assert results['6-90403795-G-A']['g_hgvs'] == 'NC_000006.11:g.90403795G>A'
        assert results['6-90403795-G-A']['genomic_variant_error'] is None
        assert 'NM_014611.2' in results['6-90403795-G-A']['hgvs_t_and_p'].keys()
        assert results['6-90403795-G-A']['hgvs_t_and_p']['NM_014611.2']['t_hgvs'] == 'NM_014611.2:c.9879C>T'
        assert results['6-90403795-G-A']['hgvs_t_and_p']['NM_014611.2']['p_hgvs_tlc'] == 'NP_055426.1:p.(Val3293=)'
        assert results['6-90403795-G-A']['hgvs_t_and_p']['NM_014611.2']['p_hgvs_slc'] == 'NP_055426.1:p.(V3293=)'
        assert results['6-90403795-G-A']['hgvs_t_and_p']['NM_014611.2']['transcript_variant_error'] is None
        assert 'NM_014611.1' in results['6-90403795-G-A']['hgvs_t_and_p'].keys()
        assert results['6-90403795-G-A']['hgvs_t_and_p']['NM_014611.1']['t_hgvs'] == 'NM_014611.1:c.9879='
        assert results['6-90403795-G-A']['hgvs_t_and_p']['NM_014611.1']['p_hgvs_tlc'] == 'NP_055426.1:p.(Val3293=)'
        assert results['6-90403795-G-A']['hgvs_t_and_p']['NM_014611.1']['p_hgvs_slc'] == 'NP_055426.1:p.(V3293=)'
        assert results['6-90403795-G-A']['hgvs_t_and_p']['NM_014611.1']['transcript_variant_error'] is None

    def test_variant51(self):
        variant = '6-32012992-CG-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '6-32012992-CG-C' in results.keys()
        assert results['6-32012992-CG-C']['p_vcf'] == '6-32012992-CG-C'
        assert results['6-32012992-CG-C']['g_hgvs'] == 'NC_000006.11:g.32012993del'
        assert results['6-32012992-CG-C']['genomic_variant_error'] is None
        assert 'NM_019105.7' in results['6-32012992-CG-C']['hgvs_t_and_p'].keys()
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_019105.7']['t_hgvs'] == 'NM_019105.7:c.10711del'
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_019105.7']['p_hgvs_tlc'] == 'NP_061978.6:p.(Arg3571AlafsTer91)'
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_019105.7']['p_hgvs_slc'] == 'NP_061978.6:p.(R3571Afs*91)'
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_019105.7']['transcript_variant_error'] is None
        assert 'NM_019105.6' in results['6-32012992-CG-C']['hgvs_t_and_p'].keys()
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_019105.6']['t_hgvs'] == 'NM_019105.6:c.10711del'
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_019105.6']['p_hgvs_tlc'] == 'NP_061978.6:p.(Arg3571AlafsTer91)'
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_019105.6']['p_hgvs_slc'] == 'NP_061978.6:p.(R3571Afs*91)'
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_019105.6']['transcript_variant_error'] is None
        assert 'NM_032470.3' in results['6-32012992-CG-C']['hgvs_t_and_p'].keys()
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_032470.3']['t_hgvs'] == 'NM_032470.3:c.4del'
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_032470.3']['p_hgvs_tlc'] == 'NP_115859.2:p.(Arg2AlafsTer91)'
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_032470.3']['p_hgvs_slc'] == 'NP_115859.2:p.(R2Afs*91)'
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_032470.3']['transcript_variant_error'] is None
        assert 'NM_001365276.1' in results['6-32012992-CG-C']['hgvs_t_and_p'].keys()
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_001365276.1']['t_hgvs'] == 'NM_001365276.1:c.10717del'
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_001365276.1']['p_hgvs_tlc'] == 'NP_001352205.1:p.(Arg3573AlafsTer91)'
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_001365276.1']['p_hgvs_slc'] == 'NP_001352205.1:p.(R3573Afs*91)'
        assert results['6-32012992-CG-C']['hgvs_t_and_p']['NM_001365276.1']['transcript_variant_error'] is None

    def test_variant52(self):
        variant = '17-48275363-C-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-48275363-C-A' in results.keys()
        assert results['17-48275363-C-A']['p_vcf'] == '17-48275363-C-A'
        assert results['17-48275363-C-A']['g_hgvs'] == 'NC_000017.10:g.48275363C>A'
        assert results['17-48275363-C-A']['genomic_variant_error'] is None
        assert 'NM_000088.3' in results['17-48275363-C-A']['hgvs_t_and_p'].keys()
        assert results['17-48275363-C-A']['hgvs_t_and_p']['NM_000088.3']['t_hgvs'] == 'NM_000088.3:c.589G>T'
        assert results['17-48275363-C-A']['hgvs_t_and_p']['NM_000088.3']['p_hgvs_tlc'] == 'NP_000079.2:p.(Gly197Cys)'
        assert results['17-48275363-C-A']['hgvs_t_and_p']['NM_000088.3']['p_hgvs_slc'] == 'NP_000079.2:p.(G197C)'
        assert results['17-48275363-C-A']['hgvs_t_and_p']['NM_000088.3']['transcript_variant_error'] is None

    def test_variant53(self):
        variant = '17-48275364-C-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-48275364-C-A' in results.keys()
        assert results['17-48275364-C-A']['p_vcf'] == '17-48275364-C-A'
        assert results['17-48275364-C-A']['g_hgvs'] == 'NC_000017.10:g.48275364C>A'
        assert results['17-48275364-C-A']['genomic_variant_error'] is None
        assert 'NM_000088.3' in results['17-48275364-C-A']['hgvs_t_and_p'].keys()
        assert results['17-48275364-C-A']['hgvs_t_and_p']['NM_000088.3']['t_hgvs'] == 'NM_000088.3:c.589-1G>T'
        assert results['17-48275364-C-A']['hgvs_t_and_p']['NM_000088.3']['p_hgvs_tlc'] == 'NP_000079.2:p.?'
        assert results['17-48275364-C-A']['hgvs_t_and_p']['NM_000088.3']['p_hgvs_slc'] == 'NP_000079.2:p.?'
        assert results['17-48275364-C-A']['hgvs_t_and_p']['NM_000088.3']['transcript_variant_error'] is None

    def test_variant54(self):
        variant = '17-48275359-GGA-TCC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-48275359-GGA-TCC' in results.keys()
        assert results['17-48275359-GGA-TCC']['p_vcf'] == '17-48275359-GGA-TCC'
        assert results['17-48275359-GGA-TCC']['g_hgvs'] == 'NC_000017.10:g.48275359_48275361inv'
        assert results['17-48275359-GGA-TCC']['genomic_variant_error'] is None
        assert 'NM_000088.3' in results['17-48275359-GGA-TCC']['hgvs_t_and_p'].keys()
        assert results['17-48275359-GGA-TCC']['hgvs_t_and_p']['NM_000088.3']['t_hgvs'] == 'NM_000088.3:c.591_593inv'
        assert results['17-48275359-GGA-TCC']['hgvs_t_and_p']['NM_000088.3']['p_hgvs_tlc'] == 'NP_000079.2:p.(Pro198Asp)'
        assert results['17-48275359-GGA-TCC']['hgvs_t_and_p']['NM_000088.3']['p_hgvs_slc'] == 'NP_000079.2:p.(P198D)'
        assert results['17-48275359-GGA-TCC']['hgvs_t_and_p']['NM_000088.3']['transcript_variant_error'] is None

    def test_variant55(self):
        variant = '7-94039128-CTTG-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '7-94039128-CTTG-C' in results.keys()
        assert results['7-94039128-CTTG-C']['p_vcf'] == '7-94039128-CTTG-C'
        assert results['7-94039128-CTTG-C']['g_hgvs'] == 'NC_000007.13:g.94039133_94039135del'
        assert results['7-94039128-CTTG-C']['genomic_variant_error'] is None
        assert 'NM_000089.3' in results['7-94039128-CTTG-C']['hgvs_t_and_p'].keys()
        assert results['7-94039128-CTTG-C']['hgvs_t_and_p']['NM_000089.3']['t_hgvs'] == 'NM_000089.3:c.1035_1035+2del'
        assert results['7-94039128-CTTG-C']['hgvs_t_and_p']['NM_000089.3']['p_hgvs_tlc'] == 'NP_000080.2:p.?'
        assert results['7-94039128-CTTG-C']['hgvs_t_and_p']['NM_000089.3']['p_hgvs_slc'] == 'NP_000080.2:p.?'
        assert results['7-94039128-CTTG-C']['hgvs_t_and_p']['NM_000089.3']['transcript_variant_error'] is None

    def test_variant56(self):
        variant = '9-135800972-AC-ACC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '9-135800972-AC-ACC' in results.keys()
        assert results['9-135800972-AC-ACC']['p_vcf'] == '9-135800972-A-AC'
        assert results['9-135800972-AC-ACC']['g_hgvs'] == 'NC_000009.11:g.135800974dup'
        assert results['9-135800972-AC-ACC']['genomic_variant_error'] is None
        assert 'NM_001162427.1' in results['9-135800972-AC-ACC']['hgvs_t_and_p'].keys()
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_001162427.1']['t_hgvs'] == 'NM_001162427.1:c.210+1615dup'
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_001162427.1']['p_hgvs_tlc'] == 'NP_001155899.1:p.?'
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_001162427.1']['p_hgvs_slc'] == 'NP_001155899.1:p.?'
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_001162427.1']['transcript_variant_error'] is None
        assert 'NM_001362177.1' in results['9-135800972-AC-ACC']['hgvs_t_and_p'].keys()
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_001362177.1']['t_hgvs'] == 'NM_001362177.1:c.-1+1dup'
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_001362177.1']['p_hgvs_tlc'] == 'NP_001349106.1:p.?'
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_001362177.1']['p_hgvs_slc'] == 'NP_001349106.1:p.?'
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_001362177.1']['transcript_variant_error'] is None
        assert 'NM_001162426.1' in results['9-135800972-AC-ACC']['hgvs_t_and_p'].keys()
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_001162426.1']['t_hgvs'] == 'NM_001162426.1:c.363+1dup'
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_001162426.1']['p_hgvs_tlc'] == 'NP_001155898.1:p.?'
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_001162426.1']['p_hgvs_slc'] == 'NP_001155898.1:p.?'
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_001162426.1']['transcript_variant_error'] is None
        assert 'NM_000368.4' in results['9-135800972-AC-ACC']['hgvs_t_and_p'].keys()
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_000368.4']['t_hgvs'] == 'NM_000368.4:c.363+1dup'
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_000368.4']['p_hgvs_tlc'] == 'NP_000359.1:p.?'
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_000368.4']['p_hgvs_slc'] == 'NP_000359.1:p.?'
        assert results['9-135800972-AC-ACC']['hgvs_t_and_p']['NM_000368.4']['transcript_variant_error'] is None

    def test_variant57(self):
        variant = '1-43212925-C-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-43212925-C-T' in results.keys()
        assert results['1-43212925-C-T']['p_vcf'] == '1-43212925-C-T'
        assert results['1-43212925-C-T']['g_hgvs'] == 'NC_000001.10:g.43212925C>T'
        assert results['1-43212925-C-T']['genomic_variant_error'] is None
        assert 'NM_001146289.1' in results['1-43212925-C-T']['hgvs_t_and_p'].keys()
        assert results['1-43212925-C-T']['hgvs_t_and_p']['NM_001146289.1']['t_hgvs'] == 'NM_001146289.1:c.2073G>A'
        assert results['1-43212925-C-T']['hgvs_t_and_p']['NM_001146289.1']['p_hgvs_tlc'] == 'NP_001139761.1:p.(Ala691=)'
        assert results['1-43212925-C-T']['hgvs_t_and_p']['NM_001146289.1']['p_hgvs_slc'] == 'NP_001139761.1:p.(A691=)'
        assert results['1-43212925-C-T']['hgvs_t_and_p']['NM_001146289.1']['transcript_variant_error'] is None
        assert 'NM_022356.3' in results['1-43212925-C-T']['hgvs_t_and_p'].keys()
        assert results['1-43212925-C-T']['hgvs_t_and_p']['NM_022356.3']['t_hgvs'] == 'NM_022356.3:c.2055+18G>A'
        assert results['1-43212925-C-T']['hgvs_t_and_p']['NM_022356.3']['p_hgvs_tlc'] == 'NP_071751.3:p.?'
        assert results['1-43212925-C-T']['hgvs_t_and_p']['NM_022356.3']['p_hgvs_slc'] == 'NP_071751.3:p.?'
        assert results['1-43212925-C-T']['hgvs_t_and_p']['NM_022356.3']['transcript_variant_error'] is None
        assert 'NM_001243246.1' in results['1-43212925-C-T']['hgvs_t_and_p'].keys()
        assert results['1-43212925-C-T']['hgvs_t_and_p']['NM_001243246.1']['t_hgvs'] == 'NM_001243246.1:c.2073G>A'
        assert results['1-43212925-C-T']['hgvs_t_and_p']['NM_001243246.1']['p_hgvs_tlc'] == 'NP_001230175.1:p.(Ala691=)'
        assert results['1-43212925-C-T']['hgvs_t_and_p']['NM_001243246.1']['p_hgvs_slc'] == 'NP_001230175.1:p.(A691=)'
        assert results['1-43212925-C-T']['hgvs_t_and_p']['NM_001243246.1']['transcript_variant_error'] is None

    def test_variant58(self):
        variant = 'HG987_PATCH-355171-C-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'HG987_PATCH-355171-C-A' in results.keys()
        assert results['HG987_PATCH-355171-C-A']['p_vcf'] == 'HG987_PATCH-355171-C-A'
        assert results['HG987_PATCH-355171-C-A']['g_hgvs'] == 'NW_003315950.2:g.355171C>A'
        assert results['HG987_PATCH-355171-C-A']['genomic_variant_error'] is None
        assert 'NM_001194958.2' in results['HG987_PATCH-355171-C-A']['hgvs_t_and_p'].keys()
        assert results['HG987_PATCH-355171-C-A']['hgvs_t_and_p']['NM_001194958.2']['t_hgvs'] == 'NM_001194958.2:c.20C>A'
        assert results['HG987_PATCH-355171-C-A']['hgvs_t_and_p']['NM_001194958.2']['p_hgvs_tlc'] == 'NP_001181887.2:p.(Ala7Asp)'
        assert results['HG987_PATCH-355171-C-A']['hgvs_t_and_p']['NM_001194958.2']['p_hgvs_slc'] == 'NP_001181887.2:p.(A7D)'
        assert results['HG987_PATCH-355171-C-A']['hgvs_t_and_p']['NM_001194958.2']['transcript_variant_error'] is None

    def test_variant59(self):
        variant = '20-43252915-T-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '20-43252915-T-C' in results.keys()
        assert results['20-43252915-T-C']['p_vcf'] == '20-43252915-T-C'
        assert results['20-43252915-T-C']['g_hgvs'] == 'NC_000020.10:g.43252915T>C'
        assert results['20-43252915-T-C']['genomic_variant_error'] is None
        assert 'NM_001322050.1' in results['20-43252915-T-C']['hgvs_t_and_p'].keys()
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_001322050.1']['t_hgvs'] == 'NM_001322050.1:c.129A>G'
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_001322050.1']['p_hgvs_tlc'] == 'NP_001308979.1:p.(Val43=)'
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_001322050.1']['p_hgvs_slc'] == 'NP_001308979.1:p.(V43=)'
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_001322050.1']['transcript_variant_error'] is None
        assert 'NR_136160.1' in results['20-43252915-T-C']['hgvs_t_and_p'].keys()
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NR_136160.1']['t_hgvs'] == 'NR_136160.1:n.685A>G'
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NR_136160.1']['p_hgvs_tlc'] is None
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NR_136160.1']['p_hgvs_slc'] is None
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NR_136160.1']['transcript_variant_error'] is None
        assert 'NM_000022.3' in results['20-43252915-T-C']['hgvs_t_and_p'].keys()
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_000022.3']['t_hgvs'] == 'NM_000022.3:c.534A>G'
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_000022.3']['p_hgvs_tlc'] == 'NP_000013.2:p.(Val178=)'
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_000022.3']['p_hgvs_slc'] == 'NP_000013.2:p.(V178=)'
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_000022.3']['transcript_variant_error'] is None
        assert 'NM_000022.2' in results['20-43252915-T-C']['hgvs_t_and_p'].keys()
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_000022.2']['t_hgvs'] == 'NM_000022.2:c.534A>G'
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_000022.2']['p_hgvs_tlc'] == 'NP_000013.2:p.(Val178=)'
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_000022.2']['p_hgvs_slc'] == 'NP_000013.2:p.(V178=)'
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_000022.2']['transcript_variant_error'] is None
        assert 'NM_001322051.1' in results['20-43252915-T-C']['hgvs_t_and_p'].keys()
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_001322051.1']['t_hgvs'] == 'NM_001322051.1:c.534A>G'
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_001322051.1']['p_hgvs_tlc'] == 'NP_001308980.1:p.(Val178=)'
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_001322051.1']['p_hgvs_slc'] == 'NP_001308980.1:p.(V178=)'
        assert results['20-43252915-T-C']['hgvs_t_and_p']['NM_001322051.1']['transcript_variant_error'] is None

    def test_variant60(self):
        variant = '1-216219781-A-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-216219781-A-C' in results.keys()
        assert results['1-216219781-A-C']['p_vcf'] == '1-216219781-A-C'
        assert results['1-216219781-A-C']['g_hgvs'] == 'NC_000001.10:g.216219781A>C'
        assert results['1-216219781-A-C']['genomic_variant_error'] is None
        assert 'NM_206933.2' in results['1-216219781-A-C']['hgvs_t_and_p'].keys()
        assert results['1-216219781-A-C']['hgvs_t_and_p']['NM_206933.2']['t_hgvs'] == 'NM_206933.2:c.6317C>G'
        assert results['1-216219781-A-C']['hgvs_t_and_p']['NM_206933.2']['p_hgvs_tlc'] == 'NP_996816.2:p.(Thr2106Arg)'
        assert results['1-216219781-A-C']['hgvs_t_and_p']['NM_206933.2']['p_hgvs_slc'] == 'NP_996816.2:p.(T2106R)'
        assert results['1-216219781-A-C']['hgvs_t_and_p']['NM_206933.2']['transcript_variant_error'] is None

    def test_variant61(self):
        variant = 'NC_000005.9:g.35058665_35058666CA='
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000005.9:g.35058665_35058666CA=' in results.keys()
        assert results['NC_000005.9:g.35058665_35058666CA=']['p_vcf'] == '5:35058665:CA:CA'
        assert results['NC_000005.9:g.35058665_35058666CA=']['g_hgvs'] == 'NC_000005.9:g.35058665_35058666='
        assert results['NC_000005.9:g.35058665_35058666CA=']['genomic_variant_error'] is None
        assert 'NM_001204314.1' in results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p'].keys()
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204314.1']['t_hgvs'] == 'NM_001204314.1:c.*6525_*6526='
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204314.1']['p_hgvs_tlc'] == 'NP_001191243.1:p.?'
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204314.1']['p_hgvs_slc'] == 'NP_001191243.1:p.?'
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204314.1']['transcript_variant_error'] is None
        assert 'NM_001204314.2' in results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p'].keys()
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204314.2']['t_hgvs'] == 'NM_001204314.2:c.*6528_*6529='
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204314.2']['p_hgvs_tlc'] == 'NP_001191243.1:p.?'
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204314.2']['p_hgvs_slc'] == 'NP_001191243.1:p.?'
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204314.2']['transcript_variant_error'] is None
        assert 'NM_000949.5' in results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p'].keys()
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_000949.5']['t_hgvs'] == 'NM_000949.5:c.*6525_*6526='
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_000949.5']['p_hgvs_tlc'] == 'NP_000940.1:p.?'
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_000949.5']['p_hgvs_slc'] == 'NP_000940.1:p.?'
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_000949.5']['transcript_variant_error'] is None
        assert 'NM_001204317.1' in results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p'].keys()
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204317.1']['t_hgvs'] == 'NM_001204317.1:c.856-9153_856-9152='
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204317.1']['p_hgvs_tlc'] == 'NP_001191246.1:p.?'
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204317.1']['p_hgvs_slc'] == 'NP_001191246.1:p.?'
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204317.1']['transcript_variant_error'] is None
        assert 'NM_000949.6' in results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p'].keys()
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_000949.6']['t_hgvs'] == 'NM_000949.6:c.*6528_*6529='
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_000949.6']['p_hgvs_tlc'] == 'NP_000940.1:p.?'
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_000949.6']['p_hgvs_slc'] == 'NP_000940.1:p.?'
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_000949.6']['transcript_variant_error'] is None
        assert 'NM_001204316.1' in results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p'].keys()
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204316.1']['t_hgvs'] == 'NM_001204316.1:c.1009+7385_1009+7386='
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204316.1']['p_hgvs_tlc'] == 'NP_001191245.1:p.?'
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204316.1']['p_hgvs_slc'] == 'NP_001191245.1:p.?'
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204316.1']['transcript_variant_error'] is None
        assert 'NR_037910.1' in results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p'].keys()
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NR_037910.1']['t_hgvs'] == 'NR_037910.1:n.828-9153_828-9152='
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NR_037910.1']['p_hgvs_tlc'] is None
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NR_037910.1']['p_hgvs_slc'] is None
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NR_037910.1']['transcript_variant_error'] is None
        assert 'NM_001204318.1' in results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p'].keys()
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204318.1']['t_hgvs'] == 'NM_001204318.1:c.686-9153_686-9152='
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204318.1']['p_hgvs_tlc'] == 'NP_001191247.1:p.?'
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204318.1']['p_hgvs_slc'] == 'NP_001191247.1:p.?'
        assert results['NC_000005.9:g.35058665_35058666CA=']['hgvs_t_and_p']['NM_001204318.1']['transcript_variant_error'] is None

    def test_variant62(self):
        variant = 'NC_000002.11:g.73675227_73675229delTCTinsTCTCTC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000002.11:g.73675227_73675229delTCTinsTCTCTC' in results.keys()
        assert results['NC_000002.11:g.73675227_73675229delTCTinsTCTCTC']['p_vcf'] == '2:73675229:T:TCTC'
        assert results['NC_000002.11:g.73675227_73675229delTCTinsTCTCTC']['g_hgvs'] == 'NC_000002.11:g.73675231_73675232insCCT'
        assert results['NC_000002.11:g.73675227_73675229delTCTinsTCTCTC']['genomic_variant_error'] == "NC_000002.11:g.73675227_73675229delTCTinsTCTCTC updated to NC_000002.11:g.73675231_73675232insCCT"
        assert 'NM_015120.4' in results['NC_000002.11:g.73675227_73675229delTCTinsTCTCTC']['hgvs_t_and_p'].keys()
        assert results['NC_000002.11:g.73675227_73675229delTCTinsTCTCTC']['hgvs_t_and_p']['NM_015120.4']['t_hgvs'] == 'NM_015120.4:c.1580_1581insCCT'
        assert results['NC_000002.11:g.73675227_73675229delTCTinsTCTCTC']['hgvs_t_and_p']['NM_015120.4']['p_hgvs_tlc'] == 'NP_055935.4:p.(Leu527dup)'
        assert results['NC_000002.11:g.73675227_73675229delTCTinsTCTCTC']['hgvs_t_and_p']['NM_015120.4']['p_hgvs_slc'] == 'NP_055935.4:p.(L527dup)'
        assert results['NC_000002.11:g.73675227_73675229delTCTinsTCTCTC']['hgvs_t_and_p']['NM_015120.4']['transcript_variant_error'] is None

    def test_variant63(self):
        variant = 'X-122318386-A-AGG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'X-122318386-A-AGG' in results.keys()
        assert results['X-122318386-A-AGG']['p_vcf'] == 'X-122318386-A-AGG'
        assert results['X-122318386-A-AGG']['g_hgvs'] == 'NC_000023.10:g.122318386_122318387insGG'
        assert results['X-122318386-A-AGG']['genomic_variant_error'] is None
        assert 'NM_000828.4' in results['X-122318386-A-AGG']['hgvs_t_and_p'].keys()
        assert results['X-122318386-A-AGG']['hgvs_t_and_p']['NM_000828.4']['t_hgvs'] == 'NM_000828.4:c.-2dup'
        assert results['X-122318386-A-AGG']['hgvs_t_and_p']['NM_000828.4']['p_hgvs_tlc'] == 'NP_000819.3:p.?'
        assert results['X-122318386-A-AGG']['hgvs_t_and_p']['NM_000828.4']['p_hgvs_slc'] == 'NP_000819.3:p.?'
        assert results['X-122318386-A-AGG']['hgvs_t_and_p']['NM_000828.4']['transcript_variant_error'] is None
        assert 'NM_007325.4' in results['X-122318386-A-AGG']['hgvs_t_and_p'].keys()
        assert results['X-122318386-A-AGG']['hgvs_t_and_p']['NM_007325.4']['t_hgvs'] == 'NM_007325.4:c.-2dup'
        assert results['X-122318386-A-AGG']['hgvs_t_and_p']['NM_007325.4']['p_hgvs_tlc'] == 'NP_015564.4:p.?'
        assert results['X-122318386-A-AGG']['hgvs_t_and_p']['NM_007325.4']['p_hgvs_slc'] == 'NP_015564.4:p.?'
        assert results['X-122318386-A-AGG']['hgvs_t_and_p']['NM_007325.4']['transcript_variant_error'] is None
        assert 'NM_001256743.1' in results['X-122318386-A-AGG']['hgvs_t_and_p'].keys()
        assert results['X-122318386-A-AGG']['hgvs_t_and_p']['NM_001256743.1']['t_hgvs'] == 'NM_001256743.1:c.-2dup'
        assert results['X-122318386-A-AGG']['hgvs_t_and_p']['NM_001256743.1']['p_hgvs_tlc'] == 'NP_001243672.1:p.?'
        assert results['X-122318386-A-AGG']['hgvs_t_and_p']['NM_001256743.1']['p_hgvs_slc'] == 'NP_001243672.1:p.?'
        assert results['X-122318386-A-AGG']['hgvs_t_and_p']['NM_001256743.1']['transcript_variant_error'] is None

    def test_variant64(self):
        variant = 'X-122318386-A-AT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'X-122318386-A-AT' in results.keys()
        assert results['X-122318386-A-AT']['p_vcf'] == 'X-122318386-A-AT'
        assert results['X-122318386-A-AT']['g_hgvs'] == 'NC_000023.10:g.122318386_122318387insT'
        assert results['X-122318386-A-AT']['genomic_variant_error'] is None
        assert 'NM_000828.4' in results['X-122318386-A-AT']['hgvs_t_and_p'].keys()
        assert results['X-122318386-A-AT']['hgvs_t_and_p']['NM_000828.4']['t_hgvs'] == 'NM_000828.4:c.-2G>T'
        assert results['X-122318386-A-AT']['hgvs_t_and_p']['NM_000828.4']['p_hgvs_tlc'] == 'NP_000819.3:p.?'
        assert results['X-122318386-A-AT']['hgvs_t_and_p']['NM_000828.4']['p_hgvs_slc'] == 'NP_000819.3:p.?'
        assert results['X-122318386-A-AT']['hgvs_t_and_p']['NM_000828.4']['transcript_variant_error'] is None
        assert 'NM_007325.4' in results['X-122318386-A-AT']['hgvs_t_and_p'].keys()
        assert results['X-122318386-A-AT']['hgvs_t_and_p']['NM_007325.4']['t_hgvs'] == 'NM_007325.4:c.-2G>T'
        assert results['X-122318386-A-AT']['hgvs_t_and_p']['NM_007325.4']['p_hgvs_tlc'] == 'NP_015564.4:p.?'
        assert results['X-122318386-A-AT']['hgvs_t_and_p']['NM_007325.4']['p_hgvs_slc'] == 'NP_015564.4:p.?'
        assert results['X-122318386-A-AT']['hgvs_t_and_p']['NM_007325.4']['transcript_variant_error'] is None
        assert 'NM_001256743.1' in results['X-122318386-A-AT']['hgvs_t_and_p'].keys()
        assert results['X-122318386-A-AT']['hgvs_t_and_p']['NM_001256743.1']['t_hgvs'] == 'NM_001256743.1:c.-2G>T'
        assert results['X-122318386-A-AT']['hgvs_t_and_p']['NM_001256743.1']['p_hgvs_tlc'] == 'NP_001243672.1:p.?'
        assert results['X-122318386-A-AT']['hgvs_t_and_p']['NM_001256743.1']['p_hgvs_slc'] == 'NP_001243672.1:p.?'
        assert results['X-122318386-A-AT']['hgvs_t_and_p']['NM_001256743.1']['transcript_variant_error'] is None

    def test_variant65(self):
        variant = '15-72105929-C-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '15-72105929-C-C' in results.keys()
        assert results['15-72105929-C-C']['p_vcf'] == '15-72105929-C-C'
        assert results['15-72105929-C-C']['g_hgvs'] == 'NC_000015.9:g.72105929='
        assert results['15-72105929-C-C']['genomic_variant_error'] is None
        assert 'NM_014249.3' in results['15-72105929-C-C']['hgvs_t_and_p'].keys()
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_014249.3']['t_hgvs'] == 'NM_014249.3:c.951dup'
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_014249.3']['p_hgvs_tlc'] == 'NP_055064.1:p.(Thr318HisfsTer23)'
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_014249.3']['p_hgvs_slc'] == 'NP_055064.1:p.(T318Hfs*23)'
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_014249.3']['transcript_variant_error'] is None
        assert 'NM_014249.2' in results['15-72105929-C-C']['hgvs_t_and_p'].keys()
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_014249.2']['t_hgvs'] == 'NM_014249.2:c.951dup'
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_014249.2']['p_hgvs_tlc'] == 'NP_055064.1:p.(Thr318HisfsTer23)'
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_014249.2']['p_hgvs_slc'] == 'NP_055064.1:p.(T318Hfs*23)'
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_014249.2']['transcript_variant_error'] is None
        assert 'NM_016346.3' in results['15-72105929-C-C']['hgvs_t_and_p'].keys()
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_016346.3']['t_hgvs'] == 'NM_016346.3:c.951dup'
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_016346.3']['p_hgvs_tlc'] == 'NP_057430.1:p.(Thr318HisfsTer23)'
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_016346.3']['p_hgvs_slc'] == 'NP_057430.1:p.(T318Hfs*23)'
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_016346.3']['transcript_variant_error'] is None
        assert 'NM_016346.2' in results['15-72105929-C-C']['hgvs_t_and_p'].keys()
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_016346.2']['t_hgvs'] == 'NM_016346.2:c.951dup'
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_016346.2']['p_hgvs_tlc'] == 'NP_057430.1:p.(Thr318HisfsTer23)'
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_016346.2']['p_hgvs_slc'] == 'NP_057430.1:p.(T318Hfs*23)'
        assert results['15-72105929-C-C']['hgvs_t_and_p']['NM_016346.2']['transcript_variant_error'] is None

    def test_variant66(self):
        variant = '15-72105928-AC-ATT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '15-72105928-AC-ATT' in results.keys()
        assert results['15-72105928-AC-ATT']['p_vcf'] == '15-72105929-C-TT'
        assert results['15-72105928-AC-ATT']['g_hgvs'] == 'NC_000015.9:g.72105929delinsTT'
        assert results['15-72105928-AC-ATT']['genomic_variant_error'] is None
        assert 'NM_014249.3' in results['15-72105928-AC-ATT']['hgvs_t_and_p'].keys()
        assert results['15-72105928-AC-ATT']['hgvs_t_and_p']['NM_014249.3']['t_hgvs'] == "NM_014249.3:c.947_948insTT"
        assert results['15-72105928-AC-ATT']['hgvs_t_and_p']['NM_014249.3']['p_hgvs_tlc'] == "NP_055064.1:p.(Pro317SerfsTer8)"
        assert results['15-72105928-AC-ATT']['hgvs_t_and_p']['NM_014249.3']['p_hgvs_slc'] == "NP_055064.1:p.(P317Sfs*8)"
        assert results['15-72105928-AC-ATT']['hgvs_t_and_p']['NM_014249.3']['transcript_variant_error'] is None
        assert 'NM_014249.2' in results['15-72105928-AC-ATT']['hgvs_t_and_p'].keys()
        assert results['15-72105928-AC-ATT']['hgvs_t_and_p']['NM_014249.2']['t_hgvs'] == "NM_014249.2:c.947_948insTT"
        assert results['15-72105928-AC-ATT']['hgvs_t_and_p']['NM_014249.2']['p_hgvs_tlc'] == "NP_055064.1:p.(Pro317SerfsTer8)"
        assert results['15-72105928-AC-ATT']['hgvs_t_and_p']['NM_014249.2']['p_hgvs_slc'] == "NP_055064.1:p.(P317Sfs*8)"
        assert results['15-72105928-AC-ATT']['hgvs_t_and_p']['NM_014249.2']['transcript_variant_error'] is None

    def test_variant67(self):
        variant = '15-72105928-ACC-ATT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '15-72105928-ACC-ATT' in results.keys()
        assert results['15-72105928-ACC-ATT']['p_vcf'] == '15-72105929-CC-TT'
        assert results['15-72105928-ACC-ATT']['g_hgvs'] == 'NC_000015.9:g.72105929_72105930delinsTT'
        assert results['15-72105928-ACC-ATT']['genomic_variant_error'] is None
        assert 'NM_014249.3' in results['15-72105928-ACC-ATT']['hgvs_t_and_p'].keys()
        assert results['15-72105928-ACC-ATT']['hgvs_t_and_p']['NM_014249.3']['t_hgvs'] == "NM_014249.3:c.947_948insTT"
        assert results['15-72105928-ACC-ATT']['hgvs_t_and_p']['NM_014249.3']['p_hgvs_tlc'] == "NP_055064.1:p.(Pro317SerfsTer8)"
        assert results['15-72105928-ACC-ATT']['hgvs_t_and_p']['NM_014249.3']['p_hgvs_slc'] == "NP_055064.1:p.(P317Sfs*8)"
        assert results['15-72105928-ACC-ATT']['hgvs_t_and_p']['NM_014249.3']['transcript_variant_error'] is None
        assert 'NM_014249.2' in results['15-72105928-ACC-ATT']['hgvs_t_and_p'].keys()
        assert results['15-72105928-ACC-ATT']['hgvs_t_and_p']['NM_014249.2']['t_hgvs'] == "NM_014249.2:c.947_948insTT"
        assert results['15-72105928-ACC-ATT']['hgvs_t_and_p']['NM_014249.2']['p_hgvs_tlc'] == "NP_055064.1:p.(Pro317SerfsTer8)"
        assert results['15-72105928-ACC-ATT']['hgvs_t_and_p']['NM_014249.2']['p_hgvs_slc'] == "NP_055064.1:p.(P317Sfs*8)"
        assert results['15-72105928-ACC-ATT']['hgvs_t_and_p']['NM_014249.2']['transcript_variant_error'] is None

    def test_variant68(self):
        variant = '15-72105927-GACC-GTT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '15-72105927-GACC-GTT' in results.keys()
        assert results['15-72105927-GACC-GTT']['p_vcf'] == '15-72105928-ACC-TT'
        assert results['15-72105927-GACC-GTT']['g_hgvs'] == 'NC_000015.9:g.72105928_72105930delinsTT'
        assert results['15-72105927-GACC-GTT']['genomic_variant_error'] is None
        assert 'NM_014249.3' in results['15-72105927-GACC-GTT']['hgvs_t_and_p'].keys()
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_014249.3']['t_hgvs'] == 'NM_014249.3:c.947delinsTT'
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_014249.3']['p_hgvs_tlc'] == 'NP_055064.1:p.(Asp316ValfsTer25)'
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_014249.3']['p_hgvs_slc'] == 'NP_055064.1:p.(D316Vfs*25)'
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_014249.3']['transcript_variant_error'] is None
        assert 'NM_014249.2' in results['15-72105927-GACC-GTT']['hgvs_t_and_p'].keys()
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_014249.2']['t_hgvs'] == 'NM_014249.2:c.947delinsTT'
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_014249.2']['p_hgvs_tlc'] == 'NP_055064.1:p.(Asp316ValfsTer25)'
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_014249.2']['p_hgvs_slc'] == 'NP_055064.1:p.(D316Vfs*25)'
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_014249.2']['transcript_variant_error'] is None
        assert 'NM_016346.3' in results['15-72105927-GACC-GTT']['hgvs_t_and_p'].keys()
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_016346.3']['t_hgvs'] == 'NM_016346.3:c.947delinsTT'
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_016346.3']['p_hgvs_tlc'] == 'NP_057430.1:p.(Asp316ValfsTer25)'
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_016346.3']['p_hgvs_slc'] == 'NP_057430.1:p.(D316Vfs*25)'
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_016346.3']['transcript_variant_error'] is None
        assert 'NM_016346.2' in results['15-72105927-GACC-GTT']['hgvs_t_and_p'].keys()
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_016346.2']['t_hgvs'] == 'NM_016346.2:c.947delinsTT'
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_016346.2']['p_hgvs_tlc'] == 'NP_057430.1:p.(Asp316ValfsTer25)'
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_016346.2']['p_hgvs_slc'] == 'NP_057430.1:p.(D316Vfs*25)'
        assert results['15-72105927-GACC-GTT']['hgvs_t_and_p']['NM_016346.2']['transcript_variant_error'] is None

    def test_variant69(self):
        variant = '19-41123093-A-AG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '19-41123093-A-AG' in results.keys()
        assert results['19-41123093-A-AG']['p_vcf'] == '19-41123093-A-AG'
        assert results['19-41123093-A-AG']['g_hgvs'] == 'NC_000019.9:g.41123095dup'
        assert results['19-41123093-A-AG']['genomic_variant_error'] is None
        assert 'NM_001042544.1' in results['19-41123093-A-AG']['hgvs_t_and_p'].keys()
        assert results['19-41123093-A-AG']['hgvs_t_and_p']['NM_001042544.1']['t_hgvs'] == 'NM_001042544.1:c.3233_3235='
        assert results['19-41123093-A-AG']['hgvs_t_and_p']['NM_001042544.1']['p_hgvs_tlc'] == 'NP_001036009.1:p.(Gln1078_Gly1079=)'
        assert results['19-41123093-A-AG']['hgvs_t_and_p']['NM_001042544.1']['transcript_variant_error'] is None
        assert 'NM_003573.2' in results['19-41123093-A-AG']['hgvs_t_and_p'].keys()
        assert results['19-41123093-A-AG']['hgvs_t_and_p']['NM_003573.2']['t_hgvs'] == 'NM_003573.2:c.3122_3124='
        assert results['19-41123093-A-AG']['hgvs_t_and_p']['NM_003573.2']['p_hgvs_tlc'] == 'NP_003564.2:p.(Gln1041_Gly1042=)'
        assert results['19-41123093-A-AG']['hgvs_t_and_p']['NM_003573.2']['transcript_variant_error'] is None
        assert 'NM_001042545.1' in results['19-41123093-A-AG']['hgvs_t_and_p'].keys()
        assert results['19-41123093-A-AG']['hgvs_t_and_p']['NM_001042545.1']['t_hgvs'] == 'NM_001042545.1:c.3032_3034='
        assert results['19-41123093-A-AG']['hgvs_t_and_p']['NM_001042545.1']['p_hgvs_tlc'] == 'NP_001036010.1:p.(Gln1011_Gly1012=)'
        assert results['19-41123093-A-AG']['hgvs_t_and_p']['NM_001042545.1']['transcript_variant_error'] is None

    def test_variant70(self):
        variant = '19-41123093-A-AT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '19-41123093-A-AT' in results.keys()
        assert results['19-41123093-A-AT']['p_vcf'] == '19-41123093-A-AT'
        assert results['19-41123093-A-AT']['g_hgvs'] == 'NC_000019.9:g.41123093_41123094insT'
        assert results['19-41123093-A-AT']['genomic_variant_error'] is None
        assert 'NM_001042544.1' in results['19-41123093-A-AT']['hgvs_t_and_p'].keys()
        assert results['19-41123093-A-AT']['hgvs_t_and_p']['NM_001042544.1']['t_hgvs'] == 'NM_001042544.1:c.3234G>T'
        assert results['19-41123093-A-AT']['hgvs_t_and_p']['NM_001042544.1']['p_hgvs_tlc'] == 'NP_001036009.1:p.(Gln1078His)'
        assert results['19-41123093-A-AT']['hgvs_t_and_p']['NM_001042544.1']['p_hgvs_slc'] == 'NP_001036009.1:p.(Q1078H)'
        assert results['19-41123093-A-AT']['hgvs_t_and_p']['NM_001042544.1']['transcript_variant_error'] is None
        assert 'NM_003573.2' in results['19-41123093-A-AT']['hgvs_t_and_p'].keys()
        assert results['19-41123093-A-AT']['hgvs_t_and_p']['NM_003573.2']['t_hgvs'] == 'NM_003573.2:c.3123G>T'
        assert results['19-41123093-A-AT']['hgvs_t_and_p']['NM_003573.2']['p_hgvs_tlc'] == 'NP_003564.2:p.(Gln1041His)'
        assert results['19-41123093-A-AT']['hgvs_t_and_p']['NM_003573.2']['p_hgvs_slc'] == 'NP_003564.2:p.(Q1041H)'
        assert results['19-41123093-A-AT']['hgvs_t_and_p']['NM_003573.2']['transcript_variant_error'] is None
        assert 'NM_001042545.1' in results['19-41123093-A-AT']['hgvs_t_and_p'].keys()
        assert results['19-41123093-A-AT']['hgvs_t_and_p']['NM_001042545.1']['t_hgvs'] == 'NM_001042545.1:c.3033G>T'
        assert results['19-41123093-A-AT']['hgvs_t_and_p']['NM_001042545.1']['p_hgvs_tlc'] == 'NP_001036010.1:p.(Gln1011His)'
        assert results['19-41123093-A-AT']['hgvs_t_and_p']['NM_001042545.1']['p_hgvs_slc'] == 'NP_001036010.1:p.(Q1011H)'
        assert results['19-41123093-A-AT']['hgvs_t_and_p']['NM_001042545.1']['transcript_variant_error'] is None

    def test_variant71(self):
        variant = '19-41123093-AG-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '19-41123093-AG-A' in results.keys()
        assert results['19-41123093-AG-A']['p_vcf'] == '19-41123093-AG-A'
        assert results['19-41123093-AG-A']['g_hgvs'] == 'NC_000019.9:g.41123095del'
        assert results['19-41123093-AG-A']['genomic_variant_error'] is None
        assert 'NM_001042544.1' in results['19-41123093-AG-A']['hgvs_t_and_p'].keys()
        assert results['19-41123093-AG-A']['hgvs_t_and_p']['NM_001042544.1']['t_hgvs'] == 'NM_001042544.1:c.3235_3236del'
        assert results['19-41123093-AG-A']['hgvs_t_and_p']['NM_001042544.1']['p_hgvs_tlc'] == 'NP_001036009.1:p.(Gly1079LeufsTer17)'
        assert results['19-41123093-AG-A']['hgvs_t_and_p']['NM_001042544.1']['p_hgvs_slc'] == 'NP_001036009.1:p.(G1079Lfs*17)'
        assert results['19-41123093-AG-A']['hgvs_t_and_p']['NM_001042544.1']['transcript_variant_error'] is None
        assert 'NM_003573.2' in results['19-41123093-AG-A']['hgvs_t_and_p'].keys()
        assert results['19-41123093-AG-A']['hgvs_t_and_p']['NM_003573.2']['t_hgvs'] == 'NM_003573.2:c.3124_3125del'
        assert results['19-41123093-AG-A']['hgvs_t_and_p']['NM_003573.2']['p_hgvs_tlc'] == 'NP_003564.2:p.(Gly1042LeufsTer17)'
        assert results['19-41123093-AG-A']['hgvs_t_and_p']['NM_003573.2']['p_hgvs_slc'] == 'NP_003564.2:p.(G1042Lfs*17)'
        assert results['19-41123093-AG-A']['hgvs_t_and_p']['NM_003573.2']['transcript_variant_error'] is None
        assert 'NM_001042545.1' in results['19-41123093-AG-A']['hgvs_t_and_p'].keys()
        assert results['19-41123093-AG-A']['hgvs_t_and_p']['NM_001042545.1']['t_hgvs'] == 'NM_001042545.1:c.3034_3035del'
        assert results['19-41123093-AG-A']['hgvs_t_and_p']['NM_001042545.1']['p_hgvs_tlc'] == 'NP_001036010.1:p.(Gly1012LeufsTer17)'
        assert results['19-41123093-AG-A']['hgvs_t_and_p']['NM_001042545.1']['p_hgvs_slc'] == 'NP_001036010.1:p.(G1012Lfs*17)'
        assert results['19-41123093-AG-A']['hgvs_t_and_p']['NM_001042545.1']['transcript_variant_error'] is None

    def test_variant72(self):
        variant = '19-41123093-AG-AG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '19-41123093-AG-AG' in results.keys()
        assert results['19-41123093-AG-AG']['p_vcf'] == '19-41123093-AG-AG'
        assert results['19-41123093-AG-AG']['g_hgvs'] == 'NC_000019.9:g.41123093_41123094='
        assert results['19-41123093-AG-AG']['genomic_variant_error'] is None
        assert 'NM_001042544.1' in results['19-41123093-AG-AG']['hgvs_t_and_p'].keys()
        assert results['19-41123093-AG-AG']['hgvs_t_and_p']['NM_001042544.1']['t_hgvs'] == 'NM_001042544.1:c.3236del'
        assert results['19-41123093-AG-AG']['hgvs_t_and_p']['NM_001042544.1']['p_hgvs_tlc'] == 'NP_001036009.1:p.(Gly1079ValfsTer14)'
        assert results['19-41123093-AG-AG']['hgvs_t_and_p']['NM_001042544.1']['p_hgvs_slc'] == 'NP_001036009.1:p.(G1079Vfs*14)'
        assert results['19-41123093-AG-AG']['hgvs_t_and_p']['NM_001042544.1']['transcript_variant_error'] is None
        assert 'NM_003573.2' in results['19-41123093-AG-AG']['hgvs_t_and_p'].keys()
        assert results['19-41123093-AG-AG']['hgvs_t_and_p']['NM_003573.2']['t_hgvs'] == 'NM_003573.2:c.3125del'
        assert results['19-41123093-AG-AG']['hgvs_t_and_p']['NM_003573.2']['p_hgvs_tlc'] == 'NP_003564.2:p.(Gly1042ValfsTer14)'
        assert results['19-41123093-AG-AG']['hgvs_t_and_p']['NM_003573.2']['p_hgvs_slc'] == 'NP_003564.2:p.(G1042Vfs*14)'
        assert results['19-41123093-AG-AG']['hgvs_t_and_p']['NM_003573.2']['transcript_variant_error'] is None
        assert 'NM_001042545.1' in results['19-41123093-AG-AG']['hgvs_t_and_p'].keys()
        assert results['19-41123093-AG-AG']['hgvs_t_and_p']['NM_001042545.1']['t_hgvs'] == 'NM_001042545.1:c.3035del'
        assert results['19-41123093-AG-AG']['hgvs_t_and_p']['NM_001042545.1']['p_hgvs_tlc'] == 'NP_001036010.1:p.(Gly1012ValfsTer14)'
        assert results['19-41123093-AG-AG']['hgvs_t_and_p']['NM_001042545.1']['p_hgvs_slc'] == 'NP_001036010.1:p.(G1012Vfs*14)'
        assert results['19-41123093-AG-AG']['hgvs_t_and_p']['NM_001042545.1']['transcript_variant_error'] is None

    def test_variant73(self):
        variant = '1-5935162-A-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-5935162-A-T' in results.keys()
        assert results['1-5935162-A-T']['p_vcf'] == '1-5935162-A-T'
        assert results['1-5935162-A-T']['g_hgvs'] == 'NC_000001.10:g.5935162A>T'
        assert results['1-5935162-A-T']['genomic_variant_error'] is None
        assert 'NM_001291594.1' in results['1-5935162-A-T']['hgvs_t_and_p'].keys()
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_001291594.1']['t_hgvs'] == 'NM_001291594.1:c.1282-2T>A'
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_001291594.1']['p_hgvs_tlc'] == 'NP_001278523.1:p.?'
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_001291594.1']['p_hgvs_slc'] == 'NP_001278523.1:p.?'
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_001291594.1']['transcript_variant_error'] is None
        assert 'NM_015102.3' in results['1-5935162-A-T']['hgvs_t_and_p'].keys()
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_015102.3']['t_hgvs'] == 'NM_015102.3:c.2818-2T>A'
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_015102.3']['p_hgvs_tlc'] == 'NP_055917.1:p.?'
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_015102.3']['p_hgvs_slc'] == 'NP_055917.1:p.?'
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_015102.3']['transcript_variant_error'] is None
        assert 'NM_001291593.1' in results['1-5935162-A-T']['hgvs_t_and_p'].keys()
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_001291593.1']['t_hgvs'] == 'NM_001291593.1:c.1279-2T>A'
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_001291593.1']['p_hgvs_tlc'] == 'NP_001278522.1:p.?'
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_001291593.1']['p_hgvs_slc'] == 'NP_001278522.1:p.?'
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_001291593.1']['transcript_variant_error'] is None
        assert 'NR_111987.1' in results['1-5935162-A-T']['hgvs_t_and_p'].keys()
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NR_111987.1']['t_hgvs'] == 'NR_111987.1:n.3633-2T>A'
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NR_111987.1']['p_hgvs_tlc'] is None
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NR_111987.1']['p_hgvs_slc'] is None
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NR_111987.1']['transcript_variant_error'] is None
        assert 'NM_015102.4' in results['1-5935162-A-T']['hgvs_t_and_p'].keys()
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_015102.4']['t_hgvs'] == 'NM_015102.4:c.2818-2T>A'
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_015102.4']['p_hgvs_tlc'] == 'NP_055917.1:p.?'
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_015102.4']['p_hgvs_slc'] == 'NP_055917.1:p.?'
        assert results['1-5935162-A-T']['hgvs_t_and_p']['NM_015102.4']['transcript_variant_error'] is None

    def test_variant74(self):
        variant = '1-12065948-C-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-12065948-C-T' in results.keys()
        assert results['1-12065948-C-T']['p_vcf'] == '1-12065948-C-T'
        assert results['1-12065948-C-T']['g_hgvs'] == 'NC_000001.10:g.12065948C>T'
        assert results['1-12065948-C-T']['genomic_variant_error'] is None
        assert 'NM_014874.3' in results['1-12065948-C-T']['hgvs_t_and_p'].keys()
        assert results['1-12065948-C-T']['hgvs_t_and_p']['NM_014874.3']['t_hgvs'] == 'NM_014874.3:c.1676C>T'
        assert results['1-12065948-C-T']['hgvs_t_and_p']['NM_014874.3']['p_hgvs_tlc'] == 'NP_055689.1:p.(Pro559Leu)'
        assert results['1-12065948-C-T']['hgvs_t_and_p']['NM_014874.3']['p_hgvs_slc'] == 'NP_055689.1:p.(P559L)'
        assert results['1-12065948-C-T']['hgvs_t_and_p']['NM_014874.3']['transcript_variant_error'] is None
        assert 'NM_001127660.1' in results['1-12065948-C-T']['hgvs_t_and_p'].keys()
        assert results['1-12065948-C-T']['hgvs_t_and_p']['NM_001127660.1']['t_hgvs'] == 'NM_001127660.1:c.1676C>T'
        assert results['1-12065948-C-T']['hgvs_t_and_p']['NM_001127660.1']['p_hgvs_tlc'] == 'NP_001121132.1:p.(Pro559Leu)'
        assert results['1-12065948-C-T']['hgvs_t_and_p']['NM_001127660.1']['p_hgvs_slc'] == 'NP_001121132.1:p.(P559L)'
        assert results['1-12065948-C-T']['hgvs_t_and_p']['NM_001127660.1']['transcript_variant_error'] is None

    def test_variant75(self):
        variant = '1-46655125-CTCAC-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-46655125-CTCAC-C' in results.keys()
        assert results['1-46655125-CTCAC-C']['p_vcf'] == '1-46655121-GTCAC-G'
        assert results['1-46655125-CTCAC-C']['g_hgvs'] == 'NC_000001.10:g.46655126_46655129del'
        assert results['1-46655125-CTCAC-C']['genomic_variant_error'] is None
        assert 'NM_001290130.1' in results['1-46655125-CTCAC-C']['hgvs_t_and_p'].keys()
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_001290130.1']['t_hgvs'] == 'NM_001290130.1:c.1466+5_1466+8del'
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_001290130.1']['p_hgvs_tlc'] == 'NP_001277059.1:p.?'
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_001290130.1']['p_hgvs_slc'] == 'NP_001277059.1:p.?'
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_001290130.1']['transcript_variant_error'] is None
        assert 'NM_001290129.1' in results['1-46655125-CTCAC-C']['hgvs_t_and_p'].keys()
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_001290129.1']['t_hgvs'] == 'NM_001290129.1:c.1829+5_1829+8del'
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_001290129.1']['p_hgvs_tlc'] == 'NP_001277058.1:p.?'
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_001290129.1']['p_hgvs_slc'] == 'NP_001277058.1:p.?'
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_001290129.1']['transcript_variant_error'] is None
        assert 'NM_017739.3' in results['1-46655125-CTCAC-C']['hgvs_t_and_p'].keys()
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_017739.3']['t_hgvs'] == 'NM_017739.3:c.1895+5_1895+8del'
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_017739.3']['p_hgvs_tlc'] == 'NP_060209.3:p.?'
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_017739.3']['p_hgvs_slc'] == 'NP_060209.3:p.?'
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_017739.3']['transcript_variant_error'] is None
        assert 'NM_001243766.1' in results['1-46655125-CTCAC-C']['hgvs_t_and_p'].keys()
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_001243766.1']['t_hgvs'] == 'NM_001243766.1:c.1869+31_1869+34del'
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_001243766.1']['p_hgvs_tlc'] == 'NP_001230695.1:p.?'
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_001243766.1']['p_hgvs_slc'] == 'NP_001230695.1:p.?'
        assert results['1-46655125-CTCAC-C']['hgvs_t_and_p']['NM_001243766.1']['transcript_variant_error'] is None

    def test_variant76(self):
        variant = '1-68912523-TGAGCCAGAG-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-68912523-TGAGCCAGAG-T' in results.keys()
        assert results['1-68912523-TGAGCCAGAG-T']['p_vcf'] == '1-68912523-TGAGCCAGAG-T'
        assert results['1-68912523-TGAGCCAGAG-T']['g_hgvs'] == 'NC_000001.10:g.68912525_68912533del'
        assert results['1-68912523-TGAGCCAGAG-T']['genomic_variant_error'] is None
        assert 'NM_000329.2' in results['1-68912523-TGAGCCAGAG-T']['hgvs_t_and_p'].keys()
        assert results['1-68912523-TGAGCCAGAG-T']['hgvs_t_and_p']['NM_000329.2']['t_hgvs'] == 'NM_000329.2:c.106_114del'
        assert results['1-68912523-TGAGCCAGAG-T']['hgvs_t_and_p']['NM_000329.2']['p_hgvs_tlc'] == 'NP_000320.1:p.(Leu36_Leu38del)'
        assert results['1-68912523-TGAGCCAGAG-T']['hgvs_t_and_p']['NM_000329.2']['p_hgvs_slc'] == 'NP_000320.1:p.(L36_L38del)'
        assert results['1-68912523-TGAGCCAGAG-T']['hgvs_t_and_p']['NM_000329.2']['transcript_variant_error'] is None

    def test_variant77(self):
        variant = '1-68912526-GCCAGAG-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-68912526-GCCAGAG-G' in results.keys()
        assert results['1-68912526-GCCAGAG-G']['p_vcf'] == '1-68912523-TGAGCCA-T'
        assert results['1-68912526-GCCAGAG-G']['g_hgvs'] == 'NC_000001.10:g.68912527_68912532del'
        assert results['1-68912526-GCCAGAG-G']['genomic_variant_error'] is None
        assert 'NM_000329.2' in results['1-68912526-GCCAGAG-G']['hgvs_t_and_p'].keys()
        assert results['1-68912526-GCCAGAG-G']['hgvs_t_and_p']['NM_000329.2']['t_hgvs'] == 'NM_000329.2:c.109_114del'
        assert results['1-68912526-GCCAGAG-G']['hgvs_t_and_p']['NM_000329.2']['p_hgvs_tlc'] == 'NP_000320.1:p.(Trp37_Leu38del)'
        assert results['1-68912526-GCCAGAG-G']['hgvs_t_and_p']['NM_000329.2']['p_hgvs_slc'] == 'NP_000320.1:p.(W37_L38del)'
        assert results['1-68912526-GCCAGAG-G']['hgvs_t_and_p']['NM_000329.2']['transcript_variant_error'] is None

    def test_variant78(self):
        variant = '1-109817590-G-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-109817590-G-T' in results.keys()
        assert results['1-109817590-G-T']['p_vcf'] == '1-109817590-G-T'
        assert results['1-109817590-G-T']['g_hgvs'] == 'NC_000001.10:g.109817590G>T'
        assert results['1-109817590-G-T']['genomic_variant_error'] is None
        assert 'NM_001408.2' in results['1-109817590-G-T']['hgvs_t_and_p'].keys()
        assert results['1-109817590-G-T']['hgvs_t_and_p']['NM_001408.2']['t_hgvs'] == 'NM_001408.2:c.*919G>T'
        assert results['1-109817590-G-T']['hgvs_t_and_p']['NM_001408.2']['p_hgvs_tlc'] == 'NP_001399.1:p.?'
        assert results['1-109817590-G-T']['hgvs_t_and_p']['NM_001408.2']['p_hgvs_slc'] == 'NP_001399.1:p.?'
        assert results['1-109817590-G-T']['hgvs_t_and_p']['NM_001408.2']['transcript_variant_error'] is None

    def test_variant79(self):
        variant = '1-145597475-GAAGT-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-145597475-GAAGT-G' in results.keys()
        assert results['1-145597475-GAAGT-G']['p_vcf'] == '1-145597475-GAAGT-G'
        assert results['1-145597475-GAAGT-G']['g_hgvs'] == 'NC_000001.10:g.145597477_145597480del'
        assert results['1-145597475-GAAGT-G']['genomic_variant_error'] is None
        # Deprecated in VVTA
        # assert 'NM_006468.7' in results['1-145597475-GAAGT-G']['hgvs_t_and_p'].keys()
        # assert results['1-145597475-GAAGT-G']['hgvs_t_and_p']['NM_006468.7']['t_hgvs'] == "NM_006468.7:c.1070+35_1070+38del"
        # assert results['1-145597475-GAAGT-G']['hgvs_t_and_p']['NM_006468.7']['p_hgvs_tlc'] == "NP_006459.3:p.?"
        # assert results['1-145597475-GAAGT-G']['hgvs_t_and_p']['NM_006468.7']['p_hgvs_slc'] == "NP_006459.3:p.?"
        # assert results['1-145597475-GAAGT-G']['hgvs_t_and_p']['NM_006468.7']['transcript_variant_error'] is None
        assert 'NM_006468.6' in results['1-145597475-GAAGT-G']['hgvs_t_and_p'].keys()
        assert results['1-145597475-GAAGT-G']['hgvs_t_and_p']['NM_006468.6']['t_hgvs'] == "NM_006468.6:c.1070+35_1070+38del"
        assert results['1-145597475-GAAGT-G']['hgvs_t_and_p']['NM_006468.6']['p_hgvs_tlc'] == "NP_006459.3:p.?"
        assert results['1-145597475-GAAGT-G']['hgvs_t_and_p']['NM_006468.6']['p_hgvs_slc'] == "NP_006459.3:p.?"
        assert results['1-145597475-GAAGT-G']['hgvs_t_and_p']['NM_006468.6']['transcript_variant_error'] is None
        assert 'NM_001303456.1' in results['1-145597475-GAAGT-G']['hgvs_t_and_p'].keys()
        assert results['1-145597475-GAAGT-G']['hgvs_t_and_p']['NM_001303456.1']['t_hgvs'] == "NM_001303456.1:c.1109+35_1109+38del"
        assert results['1-145597475-GAAGT-G']['hgvs_t_and_p']['NM_001303456.1']['p_hgvs_tlc'] == "NP_001290385.1:p.?"
        assert results['1-145597475-GAAGT-G']['hgvs_t_and_p']['NM_001303456.1']['p_hgvs_slc'] == "NP_001290385.1:p.?"
        assert results['1-145597475-GAAGT-G']['hgvs_t_and_p']['NM_001303456.1']['transcript_variant_error'] is None

    def test_variant80(self):
        variant = '1-153791300-CTG-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-153791300-CTG-C' in results.keys()
        assert results['1-153791300-CTG-C']['p_vcf'] == '1-153791300-CTG-C'
        assert results['1-153791300-CTG-C']['g_hgvs'] == 'NC_000001.10:g.153791302_153791303del'
        assert results['1-153791300-CTG-C']['genomic_variant_error'] is None
        assert 'NM_020699.3' in results['1-153791300-CTG-C']['hgvs_t_and_p'].keys()
        assert results['1-153791300-CTG-C']['hgvs_t_and_p']['NM_020699.3']['t_hgvs'] == 'NM_020699.3:c.562_563del'
        assert results['1-153791300-CTG-C']['hgvs_t_and_p']['NM_020699.3']['p_hgvs_tlc'] == 'NP_065750.1:p.(Gln188GlufsTer36)'
        assert results['1-153791300-CTG-C']['hgvs_t_and_p']['NM_020699.3']['p_hgvs_slc'] == 'NP_065750.1:p.(Q188Efs*36)'
        assert results['1-153791300-CTG-C']['hgvs_t_and_p']['NM_020699.3']['transcript_variant_error'] is None
        assert 'NM_020699.2' in results['1-153791300-CTG-C']['hgvs_t_and_p'].keys()
        assert results['1-153791300-CTG-C']['hgvs_t_and_p']['NM_020699.2']['t_hgvs'] == 'NM_020699.2:c.562_563del'
        assert results['1-153791300-CTG-C']['hgvs_t_and_p']['NM_020699.2']['p_hgvs_tlc'] == 'NP_065750.1:p.(Gln188GlufsTer36)'
        assert results['1-153791300-CTG-C']['hgvs_t_and_p']['NM_020699.2']['p_hgvs_slc'] == 'NP_065750.1:p.(Q188Efs*36)'
        assert results['1-153791300-CTG-C']['hgvs_t_and_p']['NM_020699.2']['transcript_variant_error'] is None

    def test_variant81(self):
        variant = '1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC' in results.keys()
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['g_hgvs'] == 'NC_000001.10:g.156104667_156104690delinsCCCC'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['genomic_variant_error'] is None
        assert 'NM_001257374.2' in results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p'].keys()
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001257374.2']['t_hgvs'] == 'NM_001257374.2:c.375_398delinsCCCC'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001257374.2']['p_hgvs_tlc'] == 'NP_001244303.1:p.(Glu126ProfsTer9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001257374.2']['p_hgvs_slc'] == 'NP_001244303.1:p.(E126Pfs*9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001257374.2']['transcript_variant_error'] is None
        assert 'NM_001257374.1' in results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p'].keys()
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001257374.1']['t_hgvs'] == 'NM_001257374.1:c.375_398delinsCCCC'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001257374.1']['p_hgvs_tlc'] == 'NP_001244303.1:p.(Glu126ProfsTer9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001257374.1']['p_hgvs_slc'] == 'NP_001244303.1:p.(E126Pfs*9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001257374.1']['transcript_variant_error'] is None
        assert 'NM_001282626.1' in results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p'].keys()
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001282626.1']['t_hgvs'] == 'NM_001282626.1:c.711_734delinsCCCC'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001282626.1']['p_hgvs_tlc'] == 'NP_001269555.1:p.(Glu238ProfsTer9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001282626.1']['p_hgvs_slc'] == 'NP_001269555.1:p.(E238Pfs*9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001282626.1']['transcript_variant_error'] is None
        assert 'NM_005572.3' in results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p'].keys()
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_005572.3']['t_hgvs'] == 'NM_005572.3:c.711_734delinsCCCC'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_005572.3']['p_hgvs_tlc'] == 'NP_005563.1:p.(Glu238ProfsTer9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_005572.3']['p_hgvs_slc'] == 'NP_005563.1:p.(E238Pfs*9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_005572.3']['transcript_variant_error'] is None
        assert 'NM_170707.3' in results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p'].keys()
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_170707.3']['t_hgvs'] == 'NM_170707.3:c.711_734delinsCCCC'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_170707.3']['p_hgvs_tlc'] == 'NP_733821.1:p.(Glu238ProfsTer9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_170707.3']['p_hgvs_slc'] == 'NP_733821.1:p.(E238Pfs*9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_170707.3']['transcript_variant_error'] is None
        assert 'NM_170708.3' in results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p'].keys()
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_170708.3']['t_hgvs'] == 'NM_170708.3:c.711_734delinsCCCC'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_170708.3']['p_hgvs_tlc'] == 'NP_733822.1:p.(Glu238ProfsTer9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_170708.3']['p_hgvs_slc'] == 'NP_733822.1:p.(E238Pfs*9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_170708.3']['transcript_variant_error'] is None
        assert 'NM_001282625.1' in results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p'].keys()
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001282625.1']['t_hgvs'] == 'NM_001282625.1:c.711_734delinsCCCC'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001282625.1']['p_hgvs_tlc'] == 'NP_001269554.1:p.(Glu238ProfsTer9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001282625.1']['p_hgvs_slc'] == 'NP_001269554.1:p.(E238Pfs*9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001282625.1']['transcript_variant_error'] is None
        assert 'NM_001282624.1' in results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p'].keys()
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001282624.1']['t_hgvs'] == 'NM_001282624.1:c.468_491delinsCCCC'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001282624.1']['p_hgvs_tlc'] == 'NP_001269553.1:p.(Glu157ProfsTer9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001282624.1']['p_hgvs_slc'] == 'NP_001269553.1:p.(E157Pfs*9)'
        assert results['1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC']['hgvs_t_and_p']['NM_001282624.1']['transcript_variant_error'] is None

    def test_variant82(self):
        variant = '1-156108541-G-GG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-156108541-G-GG' in results.keys()
        assert results['1-156108541-G-GG']['p_vcf'] == '1-156108540-C-CG'
        assert results['1-156108541-G-GG']['g_hgvs'] == 'NC_000001.10:g.156108541dup'
        assert results['1-156108541-G-GG']['genomic_variant_error'] is None
        assert 'NM_170707.3' in results['1-156108541-G-GG']['hgvs_t_and_p'].keys()
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_170707.3']['t_hgvs'] == 'NM_170707.3:c.1961dup'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_170707.3']['p_hgvs_tlc'] == 'NP_733821.1:p.(Thr655AsnfsTer49)'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_170707.3']['p_hgvs_slc'] == 'NP_733821.1:p.(T655Nfs*49)'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_170707.3']['transcript_variant_error'] is None
        assert 'NM_001257374.2' in results['1-156108541-G-GG']['hgvs_t_and_p'].keys()
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_001257374.2']['t_hgvs'] == 'NM_001257374.2:c.1625dup'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_001257374.2']['p_hgvs_tlc'] == 'NP_001244303.1:p.(Thr543AsnfsTer90)'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_001257374.2']['p_hgvs_slc'] == 'NP_001244303.1:p.(T543Nfs*90)'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_001257374.2']['transcript_variant_error'] is None
        assert 'NM_001257374.1' in results['1-156108541-G-GG']['hgvs_t_and_p'].keys()
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_001257374.1']['t_hgvs'] == 'NM_001257374.1:c.1625dup'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_001257374.1']['p_hgvs_tlc'] == 'NP_001244303.1:p.(Thr543AsnfsTer90)'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_001257374.1']['p_hgvs_slc'] == 'NP_001244303.1:p.(T543Nfs*90)'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_001257374.1']['transcript_variant_error'] is None
        assert 'NM_170708.3' in results['1-156108541-G-GG']['hgvs_t_and_p'].keys()
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_170708.3']['t_hgvs'] == 'NM_170708.3:c.1871dup'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_170708.3']['p_hgvs_tlc'] == 'NP_733822.1:p.(Thr625AsnfsTer49)'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_170708.3']['p_hgvs_slc'] == 'NP_733822.1:p.(T625Nfs*49)'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_170708.3']['transcript_variant_error'] is None
        assert 'NM_001282626.1' in results['1-156108541-G-GG']['hgvs_t_and_p'].keys()
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_001282626.1']['t_hgvs'] == 'NM_001282626.1:c.1818+143dup'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_001282626.1']['p_hgvs_tlc'] == 'NP_001269555.1:p.?'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_001282626.1']['p_hgvs_slc'] == 'NP_001269555.1:p.?'
        assert results['1-156108541-G-GG']['hgvs_t_and_p']['NM_001282626.1']['transcript_variant_error'] is None

    def test_variant83(self):
        variant = '1-161279695-T-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-161279695-T-A' in results.keys()
        assert results['1-161279695-T-A']['p_vcf'] == '1-161279695-T-A'
        assert results['1-161279695-T-A']['g_hgvs'] == 'NC_000001.10:g.161279695T>A'
        assert results['1-161279695-T-A']['genomic_variant_error'] is None
        assert 'NM_000530.7' in results['1-161279695-T-A']['hgvs_t_and_p'].keys()
        assert results['1-161279695-T-A']['hgvs_t_and_p']['NM_000530.7']['t_hgvs'] == 'NM_000530.7:c.1A>T'
        assert results['1-161279695-T-A']['hgvs_t_and_p']['NM_000530.7']['p_hgvs_tlc'] == 'NP_000521.2:p.(Met1?)'
        assert results['1-161279695-T-A']['hgvs_t_and_p']['NM_000530.7']['p_hgvs_slc'] == 'NP_000521.2:p.(M1?)'
        assert results['1-161279695-T-A']['hgvs_t_and_p']['NM_000530.7']['transcript_variant_error'] is None
        assert 'NM_000530.6' in results['1-161279695-T-A']['hgvs_t_and_p'].keys()
        assert results['1-161279695-T-A']['hgvs_t_and_p']['NM_000530.6']['t_hgvs'] == 'NM_000530.6:c.1A>T'
        assert results['1-161279695-T-A']['hgvs_t_and_p']['NM_000530.6']['p_hgvs_tlc'] == 'NP_000521.2:p.(Met1?)'
        assert results['1-161279695-T-A']['hgvs_t_and_p']['NM_000530.6']['p_hgvs_slc'] == 'NP_000521.2:p.(M1?)'
        assert results['1-161279695-T-A']['hgvs_t_and_p']['NM_000530.6']['transcript_variant_error'] is None
        assert 'NM_001315491.1' in results['1-161279695-T-A']['hgvs_t_and_p'].keys()
        assert results['1-161279695-T-A']['hgvs_t_and_p']['NM_001315491.1']['t_hgvs'] == 'NM_001315491.1:c.1A>T'
        assert results['1-161279695-T-A']['hgvs_t_and_p']['NM_001315491.1']['p_hgvs_tlc'] == 'NP_001302420.1:p.(Met1?)'
        assert results['1-161279695-T-A']['hgvs_t_and_p']['NM_001315491.1']['p_hgvs_slc'] == 'NP_001302420.1:p.(M1?)'
        assert results['1-161279695-T-A']['hgvs_t_and_p']['NM_001315491.1']['transcript_variant_error'] is None

    def test_variant84(self):
        variant = '1-169519049-T-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-169519049-T-T' in results.keys()
        assert results['1-169519049-T-T']['p_vcf'] == '1-169519049-T-T'
        assert results['1-169519049-T-T']['g_hgvs'] == 'NC_000001.10:g.169519049='
        assert results['1-169519049-T-T']['genomic_variant_error'] is None
        assert 'NM_000130.4' in results['1-169519049-T-T']['hgvs_t_and_p'].keys()
        assert results['1-169519049-T-T']['hgvs_t_and_p']['NM_000130.4']['t_hgvs'] == 'NM_000130.4:c.1601G>A'
        assert results['1-169519049-T-T']['hgvs_t_and_p']['NM_000130.4']['p_hgvs_tlc'] == 'NP_000121.2:p.(Arg534Gln)'
        assert results['1-169519049-T-T']['hgvs_t_and_p']['NM_000130.4']['p_hgvs_slc'] == 'NP_000121.2:p.(R534Q)'
        assert results['1-169519049-T-T']['hgvs_t_and_p']['NM_000130.4']['transcript_variant_error'] is None

    def test_variant85(self):
        variant = '1-226125468-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '1-226125468-G-A' in results.keys()
        assert results['1-226125468-G-A']['p_vcf'] == '1-226125468-G-A'
        assert results['1-226125468-G-A']['g_hgvs'] == 'NC_000001.10:g.226125468G>A'
        assert results['1-226125468-G-A']['genomic_variant_error'] is None
        assert 'NM_003240.4' in results['1-226125468-G-A']['hgvs_t_and_p'].keys()
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_003240.4']['t_hgvs'] == 'NM_003240.4:c.774C>T'
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_003240.4']['p_hgvs_tlc'] == 'NP_003231.2:p.(Thr258=)'
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_003240.4']['p_hgvs_slc'] == 'NP_003231.2:p.(T258=)'
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_003240.4']['transcript_variant_error'] is None
        assert 'NM_003240.3' in results['1-226125468-G-A']['hgvs_t_and_p'].keys()
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_003240.3']['t_hgvs'] == 'NM_003240.3:c.774C>T'
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_003240.3']['p_hgvs_tlc'] == 'NP_003231.2:p.(Thr258=)'
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_003240.3']['p_hgvs_slc'] == 'NP_003231.2:p.(T258=)'
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_003240.3']['transcript_variant_error'] is None
        assert 'NM_001172425.2' in results['1-226125468-G-A']['hgvs_t_and_p'].keys()
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_001172425.2']['t_hgvs'] == 'NM_001172425.2:c.672C>T'
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_001172425.2']['p_hgvs_tlc'] == 'NP_001165896.1:p.(Thr224=)'
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_001172425.2']['p_hgvs_slc'] == 'NP_001165896.1:p.(T224=)'
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_001172425.2']['transcript_variant_error'] is None
        assert 'NM_001172425.1' in results['1-226125468-G-A']['hgvs_t_and_p'].keys()
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_001172425.1']['t_hgvs'] == 'NM_001172425.1:c.672C>T'
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_001172425.1']['p_hgvs_tlc'] == 'NP_001165896.1:p.(Thr224=)'
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_001172425.1']['p_hgvs_slc'] == 'NP_001165896.1:p.(T224=)'
        assert results['1-226125468-G-A']['hgvs_t_and_p']['NM_001172425.1']['transcript_variant_error'] is None

    def test_variant86(self):
        variant = '10-89623035-CGCA-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '10-89623035-CGCA-C' in results.keys()
        assert results['10-89623035-CGCA-C']['p_vcf'] == '10-89623035-CGCA-C'
        assert results['10-89623035-CGCA-C']['g_hgvs'] == 'NC_000010.10:g.89623039_89623041del'
        assert results['10-89623035-CGCA-C']['genomic_variant_error'] is None
        assert 'NM_001126049.1' in results['10-89623035-CGCA-C']['hgvs_t_and_p'].keys()
        assert results['10-89623035-CGCA-C']['hgvs_t_and_p']['NM_001126049.1']['t_hgvs'] == 'NM_001126049.1:c.-794_-792del'
        assert results['10-89623035-CGCA-C']['hgvs_t_and_p']['NM_001126049.1']['p_hgvs_tlc'] == 'NP_001119521.1:p.?'
        assert results['10-89623035-CGCA-C']['hgvs_t_and_p']['NM_001126049.1']['p_hgvs_slc'] == 'NP_001119521.1:p.?'
        assert results['10-89623035-CGCA-C']['hgvs_t_and_p']['NM_001126049.1']['transcript_variant_error'] is None

    def test_variant87(self):
        variant = '11-62457852-C-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '11-62457852-C-A' in results.keys()
        assert results['11-62457852-C-A']['p_vcf'] == '11-62457852-C-A'
        assert results['11-62457852-C-A']['g_hgvs'] == 'NC_000011.9:g.62457852C>A'
        assert results['11-62457852-C-A']['genomic_variant_error'] is None
        assert 'NM_001122955.3' in results['11-62457852-C-A']['hgvs_t_and_p'].keys()
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NM_001122955.3']['t_hgvs'] == 'NM_001122955.3:c.1376G>T'
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NM_001122955.3']['p_hgvs_tlc'] == 'NP_001116427.1:p.(Cys459Phe)'
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NM_001122955.3']['p_hgvs_slc'] == 'NP_001116427.1:p.(C459F)'
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NM_001122955.3']['transcript_variant_error'] is None
        assert 'NR_037946.1' in results['11-62457852-C-A']['hgvs_t_and_p'].keys()
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NR_037946.1']['t_hgvs'] == 'NR_037946.1:n.3896G>T'
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NR_037946.1']['p_hgvs_tlc'] is None
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NR_037946.1']['p_hgvs_slc'] is None
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NR_037946.1']['transcript_variant_error'] is None
        assert 'NM_032667.6' in results['11-62457852-C-A']['hgvs_t_and_p'].keys()
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NM_032667.6']['t_hgvs'] == 'NM_032667.6:c.1184G>T'
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NM_032667.6']['p_hgvs_tlc'] == 'NP_116056.3:p.(Cys395Phe)'
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NM_032667.6']['p_hgvs_slc'] == 'NP_116056.3:p.(C395F)'
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NM_032667.6']['transcript_variant_error'] is None
        assert 'NM_001130702.2' in results['11-62457852-C-A']['hgvs_t_and_p'].keys()
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NM_001130702.2']['t_hgvs'] == 'NM_001130702.2:c.*178G>T'
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NM_001130702.2']['p_hgvs_tlc'] == 'NP_001124174.2:p.?'
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NM_001130702.2']['p_hgvs_slc'] == 'NP_001124174.2:p.?'
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NM_001130702.2']['transcript_variant_error'] is None
        assert 'NR_037949.1' in results['11-62457852-C-A']['hgvs_t_and_p'].keys()
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NR_037949.1']['t_hgvs'] == 'NR_037949.1:n.1984G>T'
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NR_037949.1']['p_hgvs_tlc'] is None
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NR_037949.1']['p_hgvs_slc'] is None
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NR_037949.1']['transcript_variant_error'] is None
        assert 'NR_037948.1' in results['11-62457852-C-A']['hgvs_t_and_p'].keys()
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NR_037948.1']['t_hgvs'] == 'NR_037948.1:n.1978G>T'
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NR_037948.1']['p_hgvs_tlc'] is None
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NR_037948.1']['p_hgvs_slc'] is None
        assert results['11-62457852-C-A']['hgvs_t_and_p']['NR_037948.1']['transcript_variant_error'] is None

    def test_variant88(self):
        variant = '11-108178710-A-AT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '11-108178710-A-AT' in results.keys()
        assert results['11-108178710-A-AT']['p_vcf'] == '11-108178710-A-AT'
        assert results['11-108178710-A-AT']['g_hgvs'] == 'NC_000011.9:g.108178710_108178711insT'
        assert results['11-108178710-A-AT']['genomic_variant_error'] is None
        assert 'NM_001351834.1' in results['11-108178710-A-AT']['hgvs_t_and_p'].keys()
        assert results['11-108178710-A-AT']['hgvs_t_and_p']['NM_001351834.1']['t_hgvs'] == 'NM_001351834.1:c.5761_5762insT'
        assert results['11-108178710-A-AT']['hgvs_t_and_p']['NM_001351834.1']['p_hgvs_tlc'] == 'NP_001338763.1:p.(Arg1921MetfsTer9)'
        assert results['11-108178710-A-AT']['hgvs_t_and_p']['NM_001351834.1']['p_hgvs_slc'] == 'NP_001338763.1:p.(R1921Mfs*9)'
        assert results['11-108178710-A-AT']['hgvs_t_and_p']['NM_001351834.1']['transcript_variant_error'] is None
        assert 'NM_000051.3' in results['11-108178710-A-AT']['hgvs_t_and_p'].keys()
        assert results['11-108178710-A-AT']['hgvs_t_and_p']['NM_000051.3']['t_hgvs'] == 'NM_000051.3:c.5761_5762insT'
        assert results['11-108178710-A-AT']['hgvs_t_and_p']['NM_000051.3']['p_hgvs_tlc'] == 'NP_000042.3:p.(Arg1921MetfsTer9)'
        assert results['11-108178710-A-AT']['hgvs_t_and_p']['NM_000051.3']['p_hgvs_slc'] == 'NP_000042.3:p.(R1921Mfs*9)'
        assert results['11-108178710-A-AT']['hgvs_t_and_p']['NM_000051.3']['transcript_variant_error'] is None

    def test_variant89(self):
        variant = '11-111735981-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '11-111735981-G-A' in results.keys()
        assert results['11-111735981-G-A']['p_vcf'] == '11-111735981-G-A'
        assert results['11-111735981-G-A']['g_hgvs'] == 'NC_000011.9:g.111735981G>A'
        assert results['11-111735981-G-A']['genomic_variant_error'] is None
        assert 'NM_001352422.1' in results['11-111735981-G-A']['hgvs_t_and_p'].keys()
        assert results['11-111735981-G-A']['hgvs_t_and_p']['NM_001352422.1']['t_hgvs'] == "NM_001352422.1:c.-326-7C>T"
        assert results['11-111735981-G-A']['hgvs_t_and_p']['NM_001352422.1']['p_hgvs_tlc'] == "NP_001339351.1:p.?"
        assert results['11-111735981-G-A']['hgvs_t_and_p']['NM_001352422.1']['p_hgvs_slc'] == "NP_001339351.1:p.?"
        assert results['11-111735981-G-A']['hgvs_t_and_p']['NM_001352422.1']['transcript_variant_error'] is None
        assert 'NM_001352415.1' in results['11-111735981-G-A']['hgvs_t_and_p'].keys()
        assert results['11-111735981-G-A']['hgvs_t_and_p']['NM_001352415.1']['t_hgvs'] == "NM_001352415.1:c.-108-7C>T"
        assert results['11-111735981-G-A']['hgvs_t_and_p']['NM_001352415.1']['p_hgvs_tlc'] == "NP_001339344.1:p.?"
        assert results['11-111735981-G-A']['hgvs_t_and_p']['NM_001352415.1']['p_hgvs_slc'] == "NP_001339344.1:p.?"
        assert results['11-111735981-G-A']['hgvs_t_and_p']['NM_001352415.1']['transcript_variant_error'] is None
        assert 'NM_001352410.1' in results['11-111735981-G-A']['hgvs_t_and_p'].keys()
        assert results['11-111735981-G-A']['hgvs_t_and_p']['NM_001352410.1']['t_hgvs'] == "NM_001352410.1:c.-108-7C>T"
        assert results['11-111735981-G-A']['hgvs_t_and_p']['NM_001352410.1']['p_hgvs_tlc'] == "NP_001339339.1:p.?"
        assert results['11-111735981-G-A']['hgvs_t_and_p']['NM_001352410.1']['p_hgvs_slc'] == "NP_001339339.1:p.?"
        assert results['11-111735981-G-A']['hgvs_t_and_p']['NM_001352410.1']['transcript_variant_error'] is None
    def test_variant90(self):
        variant = '12-11023080-C-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '12-11023080-C-A' in results.keys()
        assert results['12-11023080-C-A']['p_vcf'] == '12-11023080-C-A'
        assert results['12-11023080-C-A']['g_hgvs'] == 'NC_000012.11:g.11023080C>A'
        assert results['12-11023080-C-A']['genomic_variant_error'] is None
        # Deprecated in VVTA
        # assert 'NR_037918.2' in results['12-11023080-C-A']['hgvs_t_and_p'].keys()
        # assert results['12-11023080-C-A']['hgvs_t_and_p']['NR_037918.2']['t_hgvs'] == 'NR_037918.2:n.1184+11736G>T'
        # assert results['12-11023080-C-A']['hgvs_t_and_p']['NR_037918.2']['p_hgvs_tlc'] is None
        # assert results['12-11023080-C-A']['hgvs_t_and_p']['NR_037918.2']['p_hgvs_slc'] is None
        # assert results['12-11023080-C-A']['hgvs_t_and_p']['NR_037918.2']['transcript_variant_error'] is None

    def test_variant91(self):
        variant = '12-22018712-TC-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '12-22018712-TC-T' in results.keys()
        assert results['12-22018712-TC-T']['p_vcf'] == '12-22018712-TC-T'
        assert results['12-22018712-TC-T']['g_hgvs'] == 'NC_000012.11:g.22018713del'
        assert results['12-22018712-TC-T']['genomic_variant_error'] is None
        assert 'NM_020297.3' in results['12-22018712-TC-T']['hgvs_t_and_p'].keys()
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_020297.3']['t_hgvs'] == 'NM_020297.3:c.2199-1302del'
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_020297.3']['p_hgvs_tlc'] == 'NP_064693.2:p.?'
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_020297.3']['p_hgvs_slc'] == 'NP_064693.2:p.?'
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_020297.3']['transcript_variant_error'] is None
        assert 'NM_020297.2' in results['12-22018712-TC-T']['hgvs_t_and_p'].keys()
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_020297.2']['t_hgvs'] == 'NM_020297.2:c.2199-1302del'
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_020297.2']['p_hgvs_tlc'] == 'NP_064693.2:p.?'
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_020297.2']['p_hgvs_slc'] == 'NP_064693.2:p.?'
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_020297.2']['transcript_variant_error'] is None
        assert 'NM_005691.2' in results['12-22018712-TC-T']['hgvs_t_and_p'].keys()
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_005691.2']['t_hgvs'] == 'NM_005691.2:c.2199-1302del'
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_005691.2']['p_hgvs_tlc'] == 'NP_005682.2:p.?'
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_005691.2']['p_hgvs_slc'] == 'NP_005682.2:p.?'
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_005691.2']['transcript_variant_error'] is None
        assert 'NM_005691.3' in results['12-22018712-TC-T']['hgvs_t_and_p'].keys()
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_005691.3']['t_hgvs'] == 'NM_005691.3:c.2199-1302del'
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_005691.3']['p_hgvs_tlc'] == 'NP_005682.2:p.?'
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_005691.3']['p_hgvs_slc'] == 'NP_005682.2:p.?'
        assert results['12-22018712-TC-T']['hgvs_t_and_p']['NM_005691.3']['transcript_variant_error'] is None

    def test_variant92(self):
        variant = '12-52912946-T-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '12-52912946-T-C' in results.keys()
        assert results['12-52912946-T-C']['p_vcf'] == '12-52912946-T-C'
        assert results['12-52912946-T-C']['g_hgvs'] == 'NC_000012.11:g.52912946T>C'
        assert results['12-52912946-T-C']['genomic_variant_error'] is None
        assert 'NM_000424.3' in results['12-52912946-T-C']['hgvs_t_and_p'].keys()
        assert results['12-52912946-T-C']['hgvs_t_and_p']['NM_000424.3']['t_hgvs'] == 'NM_000424.3:c.556-2A>G'
        assert results['12-52912946-T-C']['hgvs_t_and_p']['NM_000424.3']['p_hgvs_tlc'] == 'NP_000415.2:p.?'
        assert results['12-52912946-T-C']['hgvs_t_and_p']['NM_000424.3']['p_hgvs_slc'] == 'NP_000415.2:p.?'
        assert results['12-52912946-T-C']['hgvs_t_and_p']['NM_000424.3']['transcript_variant_error'] is None

    def test_variant93(self):
        variant = '12-103234292-TC-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '12-103234292-TC-T' in results.keys()
        assert results['12-103234292-TC-T']['p_vcf'] == '12-103234292-TC-T'
        assert results['12-103234292-TC-T']['g_hgvs'] == 'NC_000012.11:g.103234294del'
        assert results['12-103234292-TC-T']['genomic_variant_error'] is None
        assert 'NM_000277.2' in results['12-103234292-TC-T']['hgvs_t_and_p'].keys()
        assert results['12-103234292-TC-T']['hgvs_t_and_p']['NM_000277.2']['t_hgvs'] == 'NM_000277.2:c.1200del'
        assert results['12-103234292-TC-T']['hgvs_t_and_p']['NM_000277.2']['p_hgvs_tlc'] == 'NP_000268.1:p.(Asn401ThrfsTer51)'
        assert results['12-103234292-TC-T']['hgvs_t_and_p']['NM_000277.2']['p_hgvs_slc'] == 'NP_000268.1:p.(N401Tfs*51)'
        assert results['12-103234292-TC-T']['hgvs_t_and_p']['NM_000277.2']['transcript_variant_error'] is None
        assert 'NM_000277.1' in results['12-103234292-TC-T']['hgvs_t_and_p'].keys()
        assert results['12-103234292-TC-T']['hgvs_t_and_p']['NM_000277.1']['t_hgvs'] == 'NM_000277.1:c.1200del'
        assert results['12-103234292-TC-T']['hgvs_t_and_p']['NM_000277.1']['p_hgvs_tlc'] == 'NP_000268.1:p.(Asn401ThrfsTer51)'
        assert results['12-103234292-TC-T']['hgvs_t_and_p']['NM_000277.1']['p_hgvs_slc'] == 'NP_000268.1:p.(N401Tfs*51)'
        assert results['12-103234292-TC-T']['hgvs_t_and_p']['NM_000277.1']['transcript_variant_error'] is None
        assert 'NM_001354304.1' in results['12-103234292-TC-T']['hgvs_t_and_p'].keys()
        assert results['12-103234292-TC-T']['hgvs_t_and_p']['NM_001354304.1']['t_hgvs'] == 'NM_001354304.1:c.1200del'
        assert results['12-103234292-TC-T']['hgvs_t_and_p']['NM_001354304.1']['p_hgvs_tlc'] == 'NP_001341233.1:p.(Asn401ThrfsTer51)'
        assert results['12-103234292-TC-T']['hgvs_t_and_p']['NM_001354304.1']['p_hgvs_slc'] == 'NP_001341233.1:p.(N401Tfs*51)'
        assert results['12-103234292-TC-T']['hgvs_t_and_p']['NM_001354304.1']['transcript_variant_error'] is None

    def test_variant94(self):
        variant = '12-103311124-T-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '12-103311124-T-C' in results.keys()
        assert results['12-103311124-T-C']['p_vcf'] == '12-103311124-T-C'
        assert results['12-103311124-T-C']['g_hgvs'] == 'NC_000012.11:g.103311124T>C'
        assert results['12-103311124-T-C']['genomic_variant_error'] is None
        assert 'NM_000277.2' in results['12-103311124-T-C']['hgvs_t_and_p'].keys()
        assert results['12-103311124-T-C']['hgvs_t_and_p']['NM_000277.2']['t_hgvs'] == 'NM_000277.2:c.-216A>G'
        assert results['12-103311124-T-C']['hgvs_t_and_p']['NM_000277.2']['p_hgvs_tlc'] == 'NP_000268.1:p.?'
        assert results['12-103311124-T-C']['hgvs_t_and_p']['NM_000277.2']['p_hgvs_slc'] == 'NP_000268.1:p.?'
        assert results['12-103311124-T-C']['hgvs_t_and_p']['NM_000277.2']['transcript_variant_error'] is None
        assert 'NM_000277.1' in results['12-103311124-T-C']['hgvs_t_and_p'].keys()
        assert results['12-103311124-T-C']['hgvs_t_and_p']['NM_000277.1']['t_hgvs'] == 'NM_000277.1:c.-215A>G'
        assert results['12-103311124-T-C']['hgvs_t_and_p']['NM_000277.1']['p_hgvs_tlc'] == 'NP_000268.1:p.?'
        assert results['12-103311124-T-C']['hgvs_t_and_p']['NM_000277.1']['p_hgvs_slc'] == 'NP_000268.1:p.?'
        assert results['12-103311124-T-C']['hgvs_t_and_p']['NM_000277.1']['transcript_variant_error'] is None
        assert 'NM_001354304.1' in results['12-103311124-T-C']['hgvs_t_and_p'].keys()
        assert results['12-103311124-T-C']['hgvs_t_and_p']['NM_001354304.1']['t_hgvs'] == "NM_001354304.1:c.-95-121A>G"
        assert results['12-103311124-T-C']['hgvs_t_and_p']['NM_001354304.1']['p_hgvs_tlc'] == "NP_001341233.1:p.?"
        assert results['12-103311124-T-C']['hgvs_t_and_p']['NM_001354304.1']['p_hgvs_slc'] == "NP_001341233.1:p.?"
        assert results['12-103311124-T-C']['hgvs_t_and_p']['NM_001354304.1']['transcript_variant_error'] is None

    def test_variant95(self):
        variant = '12-111064166-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '12-111064166-G-A' in results.keys()
        assert results['12-111064166-G-A']['p_vcf'] == '12-111064166-G-A'
        assert results['12-111064166-G-A']['g_hgvs'] == 'NC_000012.11:g.111064166G>A'
        assert results['12-111064166-G-A']['genomic_variant_error'] is None
        assert 'NM_001082538.2' in results['12-111064166-G-A']['hgvs_t_and_p'].keys()
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001082538.2']['t_hgvs'] == 'NM_001082538.2:c.342-1G>A'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001082538.2']['p_hgvs_tlc'] == 'NP_001076007.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001082538.2']['p_hgvs_slc'] == 'NP_001076007.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001082538.2']['transcript_variant_error'] is None
        assert 'NM_001173976.1' in results['12-111064166-G-A']['hgvs_t_and_p'].keys()
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001173976.1']['t_hgvs'] == 'NM_001173976.1:c.162-1G>A'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001173976.1']['p_hgvs_tlc'] == 'NP_001167447.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001173976.1']['p_hgvs_slc'] == 'NP_001167447.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001173976.1']['transcript_variant_error'] is None
        assert 'NM_001082537.2' in results['12-111064166-G-A']['hgvs_t_and_p'].keys()
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001082537.2']['t_hgvs'] == 'NM_001082537.2:c.342-1G>A'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001082537.2']['p_hgvs_tlc'] == 'NP_001076006.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001082537.2']['p_hgvs_slc'] == 'NP_001076006.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001082537.2']['transcript_variant_error'] is None
        assert 'NM_001319680.1' in results['12-111064166-G-A']['hgvs_t_and_p'].keys()
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001319680.1']['t_hgvs'] == 'NM_001319680.1:c.342-1G>A'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001319680.1']['p_hgvs_tlc'] == 'NP_001306609.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001319680.1']['p_hgvs_slc'] == 'NP_001306609.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001319680.1']['transcript_variant_error'] is None
        assert 'NM_001319681.1' in results['12-111064166-G-A']['hgvs_t_and_p'].keys()
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001319681.1']['t_hgvs'] == 'NM_001319681.1:c.-366-1G>A'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001319681.1']['p_hgvs_tlc'] == 'NP_001306610.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001319681.1']['p_hgvs_slc'] == 'NP_001306610.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001319681.1']['transcript_variant_error'] is None
        assert 'NM_024549.5' in results['12-111064166-G-A']['hgvs_t_and_p'].keys()
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_024549.5']['t_hgvs'] == 'NM_024549.5:c.342-1G>A'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_024549.5']['p_hgvs_tlc'] == 'NP_078825.2:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_024549.5']['p_hgvs_slc'] == 'NP_078825.2:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_024549.5']['transcript_variant_error'] is None
        assert 'NM_001173975.2' in results['12-111064166-G-A']['hgvs_t_and_p'].keys()
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001173975.2']['t_hgvs'] == 'NM_001173975.2:c.174-1G>A'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001173975.2']['p_hgvs_tlc'] == 'NP_001167446.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001173975.2']['p_hgvs_slc'] == 'NP_001167446.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001173975.2']['transcript_variant_error'] is None
        assert 'NR_135088.1' in results['12-111064166-G-A']['hgvs_t_and_p'].keys()
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NR_135088.1']['t_hgvs'] == 'NR_135088.1:n.559-1G>A'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NR_135088.1']['p_hgvs_tlc'] is None
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NR_135088.1']['p_hgvs_slc'] is None
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NR_135088.1']['transcript_variant_error'] is None
        assert 'NM_001319682.1' in results['12-111064166-G-A']['hgvs_t_and_p'].keys()
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001319682.1']['t_hgvs'] == 'NM_001319682.1:c.174-1G>A'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001319682.1']['p_hgvs_tlc'] == 'NP_001306611.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001319682.1']['p_hgvs_slc'] == 'NP_001306611.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001319682.1']['transcript_variant_error'] is None
        assert 'NM_001173975.1' in results['12-111064166-G-A']['hgvs_t_and_p'].keys()
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001173975.1']['t_hgvs'] == 'NM_001173975.1:c.174-1G>A'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001173975.1']['p_hgvs_tlc'] == 'NP_001167446.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001173975.1']['p_hgvs_slc'] == 'NP_001167446.1:p.?'
        assert results['12-111064166-G-A']['hgvs_t_and_p']['NM_001173975.1']['transcript_variant_error'] is None

    def test_variant96(self):
        variant = '12-123738430-CA-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '12-123738430-CA-C' in results.keys()
        assert results['12-123738430-CA-C']['p_vcf'] == '12-123738430-CA-C'
        assert results['12-123738430-CA-C']['g_hgvs'] == 'NC_000012.11:g.123738431del'
        assert results['12-123738430-CA-C']['genomic_variant_error'] is None
        assert 'NM_001194995.1' in results['12-123738430-CA-C']['hgvs_t_and_p'].keys()
        assert results['12-123738430-CA-C']['hgvs_t_and_p']['NM_001194995.1']['t_hgvs'] == 'NM_001194995.1:c.210del'
        assert results['12-123738430-CA-C']['hgvs_t_and_p']['NM_001194995.1']['p_hgvs_tlc'] == 'NP_001181924.1:p.(Gly72AlafsTer13)'
        assert results['12-123738430-CA-C']['hgvs_t_and_p']['NM_001194995.1']['p_hgvs_slc'] == 'NP_001181924.1:p.(G72Afs*13)'
        assert results['12-123738430-CA-C']['hgvs_t_and_p']['NM_001194995.1']['transcript_variant_error'] is None
        assert 'NM_152269.4' in results['12-123738430-CA-C']['hgvs_t_and_p'].keys()
        assert results['12-123738430-CA-C']['hgvs_t_and_p']['NM_152269.4']['t_hgvs'] == 'NM_152269.4:c.210del'
        assert results['12-123738430-CA-C']['hgvs_t_and_p']['NM_152269.4']['p_hgvs_tlc'] == 'NP_689482.1:p.(Gly72AlafsTer13)'
        assert results['12-123738430-CA-C']['hgvs_t_and_p']['NM_152269.4']['p_hgvs_slc'] == 'NP_689482.1:p.(G72Afs*13)'
        assert results['12-123738430-CA-C']['hgvs_t_and_p']['NM_152269.4']['transcript_variant_error'] is None
        assert 'NM_001143905.2' in results['12-123738430-CA-C']['hgvs_t_and_p'].keys()
        assert results['12-123738430-CA-C']['hgvs_t_and_p']['NM_001143905.2']['t_hgvs'] == 'NM_001143905.2:c.210del'
        assert results['12-123738430-CA-C']['hgvs_t_and_p']['NM_001143905.2']['p_hgvs_tlc'] == 'NP_001137377.1:p.(Gly72AlafsTer13)'
        assert results['12-123738430-CA-C']['hgvs_t_and_p']['NM_001143905.2']['p_hgvs_slc'] == 'NP_001137377.1:p.(G72Afs*13)'
        assert results['12-123738430-CA-C']['hgvs_t_and_p']['NM_001143905.2']['transcript_variant_error'] is None

    def test_variant97(self):
        variant = '13-31789169-CT-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '13-31789169-CT-C' in results.keys()
        assert results['13-31789169-CT-C']['p_vcf'] == '13-31789169-CT-C'
        assert results['13-31789169-CT-C']['g_hgvs'] == 'NC_000013.10:g.31789183del'
        assert results['13-31789169-CT-C']['genomic_variant_error'] is None
        assert 'NM_194318.3' in results['13-31789169-CT-C']['hgvs_t_and_p'].keys()
        assert results['13-31789169-CT-C']['hgvs_t_and_p']['NM_194318.3']['t_hgvs'] == 'NM_194318.3:c.71-5del'
        assert results['13-31789169-CT-C']['hgvs_t_and_p']['NM_194318.3']['p_hgvs_tlc'] == 'NP_919299.3:p.?'
        assert results['13-31789169-CT-C']['hgvs_t_and_p']['NM_194318.3']['p_hgvs_slc'] == 'NP_919299.3:p.?'
        assert results['13-31789169-CT-C']['hgvs_t_and_p']['NM_194318.3']['transcript_variant_error'] is None

    def test_variant98(self):
        variant = '14-62187287-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '14-62187287-G-A' in results.keys()
        assert results['14-62187287-G-A']['p_vcf'] == '14-62187287-G-A'
        assert results['14-62187287-G-A']['g_hgvs'] == 'NC_000014.8:g.62187287G>A'
        assert results['14-62187287-G-A']['genomic_variant_error'] is None
        assert 'NM_001243084.1' in results['14-62187287-G-A']['hgvs_t_and_p'].keys()
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NM_001243084.1']['t_hgvs'] == 'NM_001243084.1:c.295G>A'
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NM_001243084.1']['p_hgvs_tlc'] == 'NP_001230013.1:p.(Ala99Thr)'
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NM_001243084.1']['p_hgvs_slc'] == 'NP_001230013.1:p.(A99T)'
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NM_001243084.1']['transcript_variant_error'] is None
        assert 'NM_181054.2' in results['14-62187287-G-A']['hgvs_t_and_p'].keys()
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NM_181054.2']['t_hgvs'] == 'NM_181054.2:c.223G>A'
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NM_181054.2']['p_hgvs_tlc'] == 'NP_851397.1:p.(Ala75Thr)'
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NM_181054.2']['p_hgvs_slc'] == 'NP_851397.1:p.(A75T)'
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NM_181054.2']['transcript_variant_error'] is None
        assert 'NR_144368.1' in results['14-62187287-G-A']['hgvs_t_and_p'].keys()
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NR_144368.1']['t_hgvs'] == 'NR_144368.1:n.214-3552C>T'
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NR_144368.1']['p_hgvs_tlc'] is None
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NR_144368.1']['p_hgvs_slc'] is None
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NR_144368.1']['transcript_variant_error'] is None
        assert 'NM_001530.3' in results['14-62187287-G-A']['hgvs_t_and_p'].keys()
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NM_001530.3']['t_hgvs'] == 'NM_001530.3:c.223G>A'
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NM_001530.3']['p_hgvs_tlc'] == 'NP_001521.1:p.(Ala75Thr)'
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NM_001530.3']['p_hgvs_slc'] == 'NP_001521.1:p.(A75T)'
        assert results['14-62187287-G-A']['hgvs_t_and_p']['NM_001530.3']['transcript_variant_error'] is None

    def test_variant99(self):
        variant = '14-62188231-TT-GA'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '14-62188231-TT-GA' in results.keys()
        assert results['14-62188231-TT-GA']['p_vcf'] == '14-62188231-TT-GA'
        assert results['14-62188231-TT-GA']['g_hgvs'] == 'NC_000014.8:g.62188231_62188232delinsGA'
        assert results['14-62188231-TT-GA']['genomic_variant_error'] is None
        assert 'NM_001243084.1' in results['14-62188231-TT-GA']['hgvs_t_and_p'].keys()
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NM_001243084.1']['t_hgvs'] == 'NM_001243084.1:c.303_304delinsGA'
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NM_001243084.1']['p_hgvs_tlc'] == 'NP_001230013.1:p.(Asp101_Leu102delinsGluMet)'
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NM_001243084.1']['p_hgvs_slc'] == 'NP_001230013.1:p.(D101_L102delinsEM)'
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NM_001243084.1']['transcript_variant_error'] is None
        assert 'NM_181054.2' in results['14-62188231-TT-GA']['hgvs_t_and_p'].keys()
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NM_181054.2']['t_hgvs'] == 'NM_181054.2:c.231_232delinsGA'
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NM_181054.2']['p_hgvs_tlc'] == 'NP_851397.1:p.(Asp77_Leu78delinsGluMet)'
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NM_181054.2']['p_hgvs_slc'] == 'NP_851397.1:p.(D77_L78delinsEM)'
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NM_181054.2']['transcript_variant_error'] is None
        assert 'NR_144368.1' in results['14-62188231-TT-GA']['hgvs_t_and_p'].keys()
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NR_144368.1']['t_hgvs'] == 'NR_144368.1:n.214-4497_214-4496delinsTC'
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NR_144368.1']['p_hgvs_tlc'] is None
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NR_144368.1']['p_hgvs_slc'] is None
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NR_144368.1']['transcript_variant_error'] is None
        assert 'NM_001530.3' in results['14-62188231-TT-GA']['hgvs_t_and_p'].keys()
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NM_001530.3']['t_hgvs'] == 'NM_001530.3:c.231_232delinsGA'
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NM_001530.3']['p_hgvs_tlc'] == 'NP_001521.1:p.(Asp77_Leu78delinsGluMet)'
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NM_001530.3']['p_hgvs_slc'] == 'NP_001521.1:p.(D77_L78delinsEM)'
        assert results['14-62188231-TT-GA']['hgvs_t_and_p']['NM_001530.3']['transcript_variant_error'] is None

    def test_variant100(self):
        variant = '14-63174827-C-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '14-63174827-C-A' in results.keys()
        assert results['14-63174827-C-A']['p_vcf'] == '14-63174827-C-A'
        assert results['14-63174827-C-A']['g_hgvs'] == 'NC_000014.8:g.63174827C>A'
        assert results['14-63174827-C-A']['genomic_variant_error'] is None
        assert 'NM_172375.2' in results['14-63174827-C-A']['hgvs_t_and_p'].keys()
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_172375.2']['t_hgvs'] == 'NM_172375.2:c.*333G>T'
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_172375.2']['p_hgvs_tlc'] == 'NP_758963.1:p.?'
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_172375.2']['p_hgvs_slc'] == 'NP_758963.1:p.?'
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_172375.2']['transcript_variant_error'] is None
        assert 'NM_139318.3' in results['14-63174827-C-A']['hgvs_t_and_p'].keys()
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_139318.3']['t_hgvs'] == 'NM_139318.3:c.2366G>T'
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_139318.3']['p_hgvs_tlc'] == 'NP_647479.2:p.(Gly789Val)'
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_139318.3']['p_hgvs_slc'] == 'NP_647479.2:p.(G789V)'
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_139318.3']['transcript_variant_error'] is None
        assert 'NM_139318.4' in results['14-63174827-C-A']['hgvs_t_and_p'].keys()
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_139318.4']['t_hgvs'] == 'NM_139318.4:c.2366G>T'
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_139318.4']['p_hgvs_tlc'] == 'NP_647479.2:p.(Gly789Val)'
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_139318.4']['p_hgvs_slc'] == 'NP_647479.2:p.(G789V)'
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_139318.4']['transcript_variant_error'] is None
        assert 'NM_172375.1' in results['14-63174827-C-A']['hgvs_t_and_p'].keys()
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_172375.1']['t_hgvs'] == 'NM_172375.1:c.*333G>T'
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_172375.1']['p_hgvs_tlc'] == 'NP_758963.1:p.?'
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_172375.1']['p_hgvs_slc'] == 'NP_758963.1:p.?'
        assert results['14-63174827-C-A']['hgvs_t_and_p']['NM_172375.1']['transcript_variant_error'] is None

    def test_variant101(self):
        variant = '15-42680000-CA-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '15-42680000-CA-C' in results.keys()
        assert results['15-42680000-CA-C']['p_vcf'] == '15-42680000-CA-C'
        assert results['15-42680000-CA-C']['g_hgvs'] == 'NC_000015.9:g.42680002del'
        assert results['15-42680000-CA-C']['genomic_variant_error'] is None
        assert 'NM_173087.1' in results['15-42680000-CA-C']['hgvs_t_and_p'].keys()
        assert results['15-42680000-CA-C']['hgvs_t_and_p']['NM_173087.1']['t_hgvs'] == 'NM_173087.1:c.550del'
        assert results['15-42680000-CA-C']['hgvs_t_and_p']['NM_173087.1']['p_hgvs_tlc'] == 'NP_775110.1:p.(Thr184ArgfsTer36)'
        assert results['15-42680000-CA-C']['hgvs_t_and_p']['NM_173087.1']['p_hgvs_slc'] == 'NP_775110.1:p.(T184Rfs*36)'
        assert results['15-42680000-CA-C']['hgvs_t_and_p']['NM_173087.1']['transcript_variant_error'] is None
        assert 'NM_024344.1' in results['15-42680000-CA-C']['hgvs_t_and_p'].keys()
        assert results['15-42680000-CA-C']['hgvs_t_and_p']['NM_024344.1']['t_hgvs'] == 'NM_024344.1:c.550del'
        assert results['15-42680000-CA-C']['hgvs_t_and_p']['NM_024344.1']['p_hgvs_tlc'] == 'NP_077320.1:p.(Thr184ArgfsTer36)'
        assert results['15-42680000-CA-C']['hgvs_t_and_p']['NM_024344.1']['p_hgvs_slc'] == 'NP_077320.1:p.(T184Rfs*36)'
        assert results['15-42680000-CA-C']['hgvs_t_and_p']['NM_024344.1']['transcript_variant_error'] is None
        assert 'NM_000070.2' in results['15-42680000-CA-C']['hgvs_t_and_p'].keys()
        assert results['15-42680000-CA-C']['hgvs_t_and_p']['NM_000070.2']['t_hgvs'] == 'NM_000070.2:c.550del'
        assert results['15-42680000-CA-C']['hgvs_t_and_p']['NM_000070.2']['p_hgvs_tlc'] == 'NP_000061.1:p.(Thr184ArgfsTer36)'
        assert results['15-42680000-CA-C']['hgvs_t_and_p']['NM_000070.2']['p_hgvs_slc'] == 'NP_000061.1:p.(T184Rfs*36)'
        assert results['15-42680000-CA-C']['hgvs_t_and_p']['NM_000070.2']['transcript_variant_error'] is None

    def test_variant102(self):
        variant = '15-42680000-CA-CAA'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '15-42680000-CA-CAA' in results.keys()
        assert results['15-42680000-CA-CAA']['p_vcf'] == '15-42680000-C-CA'
        assert results['15-42680000-CA-CAA']['g_hgvs'] == 'NC_000015.9:g.42680002dup'
        assert results['15-42680000-CA-CAA']['genomic_variant_error'] is None
        assert 'NM_173087.1' in results['15-42680000-CA-CAA']['hgvs_t_and_p'].keys()
        assert results['15-42680000-CA-CAA']['hgvs_t_and_p']['NM_173087.1']['t_hgvs'] == 'NM_173087.1:c.550dup'
        assert results['15-42680000-CA-CAA']['hgvs_t_and_p']['NM_173087.1']['p_hgvs_tlc'] == 'NP_775110.1:p.(Thr184AsnfsTer16)'
        assert results['15-42680000-CA-CAA']['hgvs_t_and_p']['NM_173087.1']['p_hgvs_slc'] == 'NP_775110.1:p.(T184Nfs*16)'
        assert results['15-42680000-CA-CAA']['hgvs_t_and_p']['NM_173087.1']['transcript_variant_error'] is None
        assert 'NM_024344.1' in results['15-42680000-CA-CAA']['hgvs_t_and_p'].keys()
        assert results['15-42680000-CA-CAA']['hgvs_t_and_p']['NM_024344.1']['t_hgvs'] == 'NM_024344.1:c.550dup'
        assert results['15-42680000-CA-CAA']['hgvs_t_and_p']['NM_024344.1']['p_hgvs_tlc'] == 'NP_077320.1:p.(Thr184AsnfsTer16)'
        assert results['15-42680000-CA-CAA']['hgvs_t_and_p']['NM_024344.1']['p_hgvs_slc'] == 'NP_077320.1:p.(T184Nfs*16)'
        assert results['15-42680000-CA-CAA']['hgvs_t_and_p']['NM_024344.1']['transcript_variant_error'] is None
        assert 'NM_000070.2' in results['15-42680000-CA-CAA']['hgvs_t_and_p'].keys()
        assert results['15-42680000-CA-CAA']['hgvs_t_and_p']['NM_000070.2']['t_hgvs'] == 'NM_000070.2:c.550dup'
        assert results['15-42680000-CA-CAA']['hgvs_t_and_p']['NM_000070.2']['p_hgvs_tlc'] == 'NP_000061.1:p.(Thr184AsnfsTer16)'
        assert results['15-42680000-CA-CAA']['hgvs_t_and_p']['NM_000070.2']['p_hgvs_slc'] == 'NP_000061.1:p.(T184Nfs*16)'
        assert results['15-42680000-CA-CAA']['hgvs_t_and_p']['NM_000070.2']['transcript_variant_error'] is None

    def test_variant103(self):
        variant = '15-42703179-T-TTCA'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '15-42703179-T-TTCA' in results.keys()
        assert results['15-42703179-T-TTCA']['p_vcf'] == '15-42703179-T-TTCA'
        assert results['15-42703179-T-TTCA']['g_hgvs'] == 'NC_000015.9:g.42703179_42703180insTCA'
        assert results['15-42703179-T-TTCA']['genomic_variant_error'] is None
        assert 'NM_173087.1' in results['15-42703179-T-TTCA']['hgvs_t_and_p'].keys()
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173087.1']['t_hgvs'] == 'NM_173087.1:c.2085_2086insTCA'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173087.1']['p_hgvs_tlc'] == 'NP_775110.1:p.(Val695_Arg696insSer)'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173087.1']['p_hgvs_slc'] == 'NP_775110.1:p.(V695_R696insS)'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173087.1']['transcript_variant_error'] is None
        assert 'NM_173088.1' in results['15-42703179-T-TTCA']['hgvs_t_and_p'].keys()
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173088.1']['t_hgvs'] == 'NM_173088.1:c.825_826insTCA'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173088.1']['p_hgvs_tlc'] == 'NP_775111.1:p.(Val275_Arg276insSer)'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173088.1']['p_hgvs_slc'] == 'NP_775111.1:p.(V275_R276insS)'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173088.1']['transcript_variant_error'] is None
        assert 'NM_173089.1' in results['15-42703179-T-TTCA']['hgvs_t_and_p'].keys()
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173089.1']['t_hgvs'] == 'NM_173089.1:c.366_367insTCA'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173089.1']['p_hgvs_tlc'] == 'NP_775112.1:p.(Val122_Arg123insSer)'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173089.1']['p_hgvs_slc'] == 'NP_775112.1:p.(V122_R123insS)'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173089.1']['transcript_variant_error'] is None
        assert 'NM_173090.1' in results['15-42703179-T-TTCA']['hgvs_t_and_p'].keys()
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173090.1']['t_hgvs'] == 'NM_173090.1:c.366_367insTCA'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173090.1']['p_hgvs_tlc'] == 'NP_775113.1:p.(Val122_Arg123insSer)'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173090.1']['p_hgvs_slc'] == 'NP_775113.1:p.(V122_R123insS)'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_173090.1']['transcript_variant_error'] is None
        assert 'NM_000070.2' in results['15-42703179-T-TTCA']['hgvs_t_and_p'].keys()
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_000070.2']['t_hgvs'] == 'NM_000070.2:c.2361_2362insTCA'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_000070.2']['p_hgvs_tlc'] == 'NP_000061.1:p.(Val787_Arg788insSer)'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_000070.2']['p_hgvs_slc'] == 'NP_000061.1:p.(V787_R788insS)'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_000070.2']['transcript_variant_error'] is None
        assert 'NM_024344.1' in results['15-42703179-T-TTCA']['hgvs_t_and_p'].keys()
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_024344.1']['t_hgvs'] == 'NM_024344.1:c.2343_2344insTCA'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_024344.1']['p_hgvs_tlc'] == 'NP_077320.1:p.(Val781_Arg782insSer)'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_024344.1']['p_hgvs_slc'] == 'NP_077320.1:p.(V781_R782insS)'
        assert results['15-42703179-T-TTCA']['hgvs_t_and_p']['NM_024344.1']['transcript_variant_error'] is None

    def test_variant104(self):
        variant = '15-42703179-TAG-TTCATCT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '15-42703179-TAG-TTCATCT' in results.keys()
        assert results['15-42703179-TAG-TTCATCT']['p_vcf'] == '15-42703180-AG-TCATCT'
        assert results['15-42703179-TAG-TTCATCT']['g_hgvs'] == 'NC_000015.9:g.42703180_42703181delinsTCATCT'
        assert results['15-42703179-TAG-TTCATCT']['genomic_variant_error'] is None
        assert 'NM_173087.1' in results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p'].keys()
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173087.1']['t_hgvs'] == 'NM_173087.1:c.2086_2087delinsTCATCT'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173087.1']['p_hgvs_tlc'] == 'NP_775110.1:p.(Arg696SerfsTer14)'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173087.1']['p_hgvs_slc'] == 'NP_775110.1:p.(R696Sfs*14)'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173087.1']['transcript_variant_error'] is None
        assert 'NM_173088.1' in results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p'].keys()
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173088.1']['t_hgvs'] == 'NM_173088.1:c.826_827delinsTCATCT'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173088.1']['p_hgvs_tlc'] == 'NP_775111.1:p.(Arg276SerfsTer14)'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173088.1']['p_hgvs_slc'] == 'NP_775111.1:p.(R276Sfs*14)'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173088.1']['transcript_variant_error'] is None
        assert 'NM_173089.1' in results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p'].keys()
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173089.1']['t_hgvs'] == 'NM_173089.1:c.367_368delinsTCATCT'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173089.1']['p_hgvs_tlc'] == 'NP_775112.1:p.(Arg123SerfsTer14)'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173089.1']['p_hgvs_slc'] == 'NP_775112.1:p.(R123Sfs*14)'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173089.1']['transcript_variant_error'] is None
        assert 'NM_173090.1' in results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p'].keys()
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173090.1']['t_hgvs'] == 'NM_173090.1:c.367_368delinsTCATCT'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173090.1']['p_hgvs_tlc'] == 'NP_775113.1:p.(Arg123SerfsTer14)'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173090.1']['p_hgvs_slc'] == 'NP_775113.1:p.(R123Sfs*14)'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_173090.1']['transcript_variant_error'] is None
        assert 'NM_000070.2' in results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p'].keys()
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_000070.2']['t_hgvs'] == 'NM_000070.2:c.2362_2363delinsTCATCT'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_000070.2']['p_hgvs_tlc'] == 'NP_000061.1:p.(Arg788SerfsTer14)'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_000070.2']['p_hgvs_slc'] == 'NP_000061.1:p.(R788Sfs*14)'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_000070.2']['transcript_variant_error'] is None
        assert 'NM_024344.1' in results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p'].keys()
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_024344.1']['t_hgvs'] == 'NM_024344.1:c.2344_2345delinsTCATCT'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_024344.1']['p_hgvs_tlc'] == 'NP_077320.1:p.(Arg782SerfsTer14)'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_024344.1']['p_hgvs_slc'] == 'NP_077320.1:p.(R782Sfs*14)'
        assert results['15-42703179-TAG-TTCATCT']['hgvs_t_and_p']['NM_024344.1']['transcript_variant_error'] is None

    def test_variant105(self):
        variant = '15-48782203-C-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '15-48782203-C-T' in results.keys()
        assert results['15-48782203-C-T']['p_vcf'] == '15-48782203-C-T'
        assert results['15-48782203-C-T']['g_hgvs'] == 'NC_000015.9:g.48782203C>T'
        assert results['15-48782203-C-T']['genomic_variant_error'] is None
        assert 'NM_000138.4' in results['15-48782203-C-T']['hgvs_t_and_p'].keys()
        assert results['15-48782203-C-T']['hgvs_t_and_p']['NM_000138.4']['t_hgvs'] == 'NM_000138.4:c.2927G>A'
        assert results['15-48782203-C-T']['hgvs_t_and_p']['NM_000138.4']['p_hgvs_tlc'] == 'NP_000129.3:p.(Arg976His)'
        assert results['15-48782203-C-T']['hgvs_t_and_p']['NM_000138.4']['p_hgvs_slc'] == 'NP_000129.3:p.(R976H)'
        assert results['15-48782203-C-T']['hgvs_t_and_p']['NM_000138.4']['transcript_variant_error'] is None

    def test_variant106(self):
        variant = '15-72105929-CC-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '15-72105929-CC-C' in results.keys()
        assert results['15-72105929-CC-C']['p_vcf'] == '15-72105928-AC-A'
        assert results['15-72105929-CC-C']['g_hgvs'] == 'NC_000015.9:g.72105933del'
        assert results['15-72105929-CC-C']['genomic_variant_error'] is None
        assert 'NM_014249.3' in results['15-72105929-CC-C']['hgvs_t_and_p'].keys()
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_014249.3']['t_hgvs'] == 'NM_014249.3:c.947_948='
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_014249.3']['p_hgvs_tlc'] == 'NP_055064.1:p.(Asp316=)'
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_014249.3']['p_hgvs_slc'] == 'NP_055064.1:p.(D316=)'
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_014249.3']['transcript_variant_error'] is None
        assert 'NM_014249.2' in results['15-72105929-CC-C']['hgvs_t_and_p'].keys()
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_014249.2']['t_hgvs'] == 'NM_014249.2:c.947_948='
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_014249.2']['p_hgvs_tlc'] == 'NP_055064.1:p.(Asp316=)'
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_014249.2']['p_hgvs_slc'] == 'NP_055064.1:p.(D316=)'
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_014249.2']['transcript_variant_error'] is None
        assert 'NM_016346.3' in results['15-72105929-CC-C']['hgvs_t_and_p'].keys()
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_016346.3']['t_hgvs'] == 'NM_016346.3:c.947_948='
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_016346.3']['p_hgvs_tlc'] == 'NP_057430.1:p.(Asp316=)'
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_016346.3']['p_hgvs_slc'] == 'NP_057430.1:p.(D316=)'
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_016346.3']['transcript_variant_error'] is None
        assert 'NM_016346.2' in results['15-72105929-CC-C']['hgvs_t_and_p'].keys()
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_016346.2']['t_hgvs'] == 'NM_016346.2:c.947_948='
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_016346.2']['p_hgvs_tlc'] == 'NP_057430.1:p.(Asp316=)'
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_016346.2']['p_hgvs_slc'] == 'NP_057430.1:p.(D316=)'
        assert results['15-72105929-CC-C']['hgvs_t_and_p']['NM_016346.2']['transcript_variant_error'] is None

    def test_variant107(self):
        variant = '15-89873415-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '15-89873415-G-A' in results.keys()
        assert results['15-89873415-G-A']['p_vcf'] == '15-89873415-G-A'
        assert results['15-89873415-G-A']['g_hgvs'] == 'NC_000015.9:g.89873415G>A'
        assert results['15-89873415-G-A']['genomic_variant_error'] is None
        assert 'NM_001126131.1' in results['15-89873415-G-A']['hgvs_t_and_p'].keys()
        assert results['15-89873415-G-A']['hgvs_t_and_p']['NM_001126131.1']['t_hgvs'] == 'NM_001126131.1:c.752C>T'
        assert results['15-89873415-G-A']['hgvs_t_and_p']['NM_001126131.1']['p_hgvs_tlc'] == 'NP_001119603.1:p.(Thr251Ile)'
        assert results['15-89873415-G-A']['hgvs_t_and_p']['NM_001126131.1']['p_hgvs_slc'] == 'NP_001119603.1:p.(T251I)'
        assert results['15-89873415-G-A']['hgvs_t_and_p']['NM_001126131.1']['transcript_variant_error'] is None
        assert 'NM_002693.2' in results['15-89873415-G-A']['hgvs_t_and_p'].keys()
        assert results['15-89873415-G-A']['hgvs_t_and_p']['NM_002693.2']['t_hgvs'] == 'NM_002693.2:c.752C>T'
        assert results['15-89873415-G-A']['hgvs_t_and_p']['NM_002693.2']['p_hgvs_tlc'] == 'NP_002684.1:p.(Thr251Ile)'
        assert results['15-89873415-G-A']['hgvs_t_and_p']['NM_002693.2']['p_hgvs_slc'] == 'NP_002684.1:p.(T251I)'
        assert results['15-89873415-G-A']['hgvs_t_and_p']['NM_002693.2']['transcript_variant_error'] is None

    def test_variant108(self):
        variant = '16-2103394-C-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-2103394-C-T' in results.keys()
        assert results['16-2103394-C-T']['p_vcf'] == '16-2103394-C-T'
        assert results['16-2103394-C-T']['g_hgvs'] == 'NC_000016.9:g.2103394C>T'
        assert results['16-2103394-C-T']['genomic_variant_error'] is None
        assert 'NM_001077183.1' in results['16-2103394-C-T']['hgvs_t_and_p'].keys()
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001077183.1']['t_hgvs'] == 'NM_001077183.1:c.277C>T'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001077183.1']['p_hgvs_tlc'] == 'NP_001070651.1:p.(Arg93Trp)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001077183.1']['p_hgvs_slc'] == 'NP_001070651.1:p.(R93W)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001077183.1']['transcript_variant_error'] is None
        assert 'NM_001318831.1' in results['16-2103394-C-T']['hgvs_t_and_p'].keys()
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318831.1']['t_hgvs'] == 'NM_001318831.1:c.-1-2803C>T'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318831.1']['p_hgvs_tlc'] == 'NP_001305760.1:p.?'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318831.1']['p_hgvs_slc'] == 'NP_001305760.1:p.?'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318831.1']['transcript_variant_error'] is None
        assert 'NM_001318827.1' in results['16-2103394-C-T']['hgvs_t_and_p'].keys()
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318827.1']['t_hgvs'] == 'NM_001318827.1:c.226-903C>T'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318827.1']['p_hgvs_tlc'] == 'NP_001305756.1:p.?'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318827.1']['p_hgvs_slc'] == 'NP_001305756.1:p.?'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318827.1']['transcript_variant_error'] is None
        assert 'NM_001114382.1' in results['16-2103394-C-T']['hgvs_t_and_p'].keys()
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001114382.1']['t_hgvs'] == 'NM_001114382.1:c.277C>T'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001114382.1']['p_hgvs_tlc'] == 'NP_001107854.1:p.(Arg93Trp)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001114382.1']['p_hgvs_slc'] == 'NP_001107854.1:p.(R93W)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001114382.1']['transcript_variant_error'] is None
        assert 'NM_001114382.2' in results['16-2103394-C-T']['hgvs_t_and_p'].keys()
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001114382.2']['t_hgvs'] == 'NM_001114382.2:c.277C>T'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001114382.2']['p_hgvs_tlc'] == 'NP_001107854.1:p.(Arg93Trp)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001114382.2']['p_hgvs_slc'] == 'NP_001107854.1:p.(R93W)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001114382.2']['transcript_variant_error'] is None
        assert 'NM_001318829.1' in results['16-2103394-C-T']['hgvs_t_and_p'].keys()
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318829.1']['t_hgvs'] == 'NM_001318829.1:c.130C>T'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318829.1']['p_hgvs_tlc'] == 'NP_001305758.1:p.(Arg44Trp)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318829.1']['p_hgvs_slc'] == 'NP_001305758.1:p.(R44W)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318829.1']['transcript_variant_error'] is None
        assert 'NM_001318832.1' in results['16-2103394-C-T']['hgvs_t_and_p'].keys()
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318832.1']['t_hgvs'] == 'NM_001318832.1:c.310C>T'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318832.1']['p_hgvs_tlc'] == 'NP_001305761.1:p.(Arg104Trp)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318832.1']['p_hgvs_slc'] == 'NP_001305761.1:p.(R104W)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001318832.1']['transcript_variant_error'] is None
        assert 'NM_001363528.1' in results['16-2103394-C-T']['hgvs_t_and_p'].keys()
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001363528.1']['t_hgvs'] == 'NM_001363528.1:c.277C>T'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001363528.1']['p_hgvs_tlc'] == 'NP_001350457.1:p.(Arg93Trp)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001363528.1']['p_hgvs_slc'] == 'NP_001350457.1:p.(R93W)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001363528.1']['transcript_variant_error'] is None
        assert 'NM_001077183.2' in results['16-2103394-C-T']['hgvs_t_and_p'].keys()
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001077183.2']['t_hgvs'] == 'NM_001077183.2:c.277C>T'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001077183.2']['p_hgvs_tlc'] == 'NP_001070651.1:p.(Arg93Trp)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001077183.2']['p_hgvs_slc'] == 'NP_001070651.1:p.(R93W)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_001077183.2']['transcript_variant_error'] is None
        assert 'NM_021055.2' in results['16-2103394-C-T']['hgvs_t_and_p'].keys()
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_021055.2']['t_hgvs'] == 'NM_021055.2:c.277C>T'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_021055.2']['p_hgvs_tlc'] == 'NP_066399.2:p.(Arg93Trp)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_021055.2']['p_hgvs_slc'] == 'NP_066399.2:p.(R93W)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_021055.2']['transcript_variant_error'] is None
        assert 'NM_000548.4' in results['16-2103394-C-T']['hgvs_t_and_p'].keys()
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_000548.4']['t_hgvs'] == 'NM_000548.4:c.277C>T'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_000548.4']['p_hgvs_tlc'] == 'NP_000539.2:p.(Arg93Trp)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_000548.4']['p_hgvs_slc'] == 'NP_000539.2:p.(R93W)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_000548.4']['transcript_variant_error'] is None
        assert 'NM_000548.3' in results['16-2103394-C-T']['hgvs_t_and_p'].keys()
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_000548.3']['t_hgvs'] == 'NM_000548.3:c.277C>T'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_000548.3']['p_hgvs_tlc'] == 'NP_000539.2:p.(Arg93Trp)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_000548.3']['p_hgvs_slc'] == 'NP_000539.2:p.(R93W)'
        assert results['16-2103394-C-T']['hgvs_t_and_p']['NM_000548.3']['transcript_variant_error'] is None

    def test_variant109(self):
        variant = '16-3779300-C-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-3779300-C-G' in results.keys()
        assert results['16-3779300-C-G']['p_vcf'] == '16-3779300-C-G'
        assert results['16-3779300-C-G']['g_hgvs'] == 'NC_000016.9:g.3779300C>G'
        assert results['16-3779300-C-G']['genomic_variant_error'] is None
        assert 'NM_001079846.1' in results['16-3779300-C-G']['hgvs_t_and_p'].keys()
        assert results['16-3779300-C-G']['hgvs_t_and_p']['NM_001079846.1']['t_hgvs'] == 'NM_001079846.1:c.5634G>C'
        assert results['16-3779300-C-G']['hgvs_t_and_p']['NM_001079846.1']['p_hgvs_tlc'] == 'NP_001073315.1:p.(Met1878Ile)'
        assert results['16-3779300-C-G']['hgvs_t_and_p']['NM_001079846.1']['p_hgvs_slc'] == 'NP_001073315.1:p.(M1878I)'
        assert results['16-3779300-C-G']['hgvs_t_and_p']['NM_001079846.1']['transcript_variant_error'] is None
        assert 'NM_004380.2' in results['16-3779300-C-G']['hgvs_t_and_p'].keys()
        assert results['16-3779300-C-G']['hgvs_t_and_p']['NM_004380.2']['t_hgvs'] == 'NM_004380.2:c.5748G>C'
        assert results['16-3779300-C-G']['hgvs_t_and_p']['NM_004380.2']['p_hgvs_tlc'] == 'NP_004371.2:p.(Met1916Ile)'
        assert results['16-3779300-C-G']['hgvs_t_and_p']['NM_004380.2']['p_hgvs_slc'] == 'NP_004371.2:p.(M1916I)'
        assert results['16-3779300-C-G']['hgvs_t_and_p']['NM_004380.2']['transcript_variant_error'] is None

    def test_variant110(self):
        variant = '16-5128843-C-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-5128843-C-G' in results.keys()
        assert results['16-5128843-C-G']['p_vcf'] == '16-5128843-C-G'
        assert results['16-5128843-C-G']['g_hgvs'] == 'NC_000016.9:g.5128843C>G'
        assert results['16-5128843-C-G']['genomic_variant_error'] is None
        assert 'NM_019109.4' in results['16-5128843-C-G']['hgvs_t_and_p'].keys()
        assert results['16-5128843-C-G']['hgvs_t_and_p']['NM_019109.4']['t_hgvs'] == 'NM_019109.4:c.826C>G'
        assert results['16-5128843-C-G']['hgvs_t_and_p']['NM_019109.4']['p_hgvs_tlc'] == 'NP_061982.3:p.(Arg276Gly)'
        assert results['16-5128843-C-G']['hgvs_t_and_p']['NM_019109.4']['p_hgvs_slc'] == 'NP_061982.3:p.(R276G)'
        assert results['16-5128843-C-G']['hgvs_t_and_p']['NM_019109.4']['transcript_variant_error'] is None
        assert 'NM_001330504.1' in results['16-5128843-C-G']['hgvs_t_and_p'].keys()
        assert results['16-5128843-C-G']['hgvs_t_and_p']['NM_001330504.1']['t_hgvs'] == 'NM_001330504.1:c.493C>G'
        assert results['16-5128843-C-G']['hgvs_t_and_p']['NM_001330504.1']['p_hgvs_tlc'] == 'NP_001317433.1:p.(Arg165Gly)'
        assert results['16-5128843-C-G']['hgvs_t_and_p']['NM_001330504.1']['p_hgvs_slc'] == 'NP_001317433.1:p.(R165G)'
        assert results['16-5128843-C-G']['hgvs_t_and_p']['NM_001330504.1']['transcript_variant_error'] is None

    def test_variant111(self):
        variant = '16-74808559-C-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-74808559-C-T' in results.keys()
        assert results['16-74808559-C-T']['p_vcf'] == '16-74808559-C-T'
        assert results['16-74808559-C-T']['g_hgvs'] == 'NC_000016.9:g.74808559C>T'
        assert results['16-74808559-C-T']['genomic_variant_error'] is None
        assert 'NM_024306.4' in results['16-74808559-C-T']['hgvs_t_and_p'].keys()
        assert results['16-74808559-C-T']['hgvs_t_and_p']['NM_024306.4']['t_hgvs'] == 'NM_024306.4:c.95G>A'
        assert results['16-74808559-C-T']['hgvs_t_and_p']['NM_024306.4']['p_hgvs_tlc'] == 'NP_077282.3:p.(Arg32His)'
        assert results['16-74808559-C-T']['hgvs_t_and_p']['NM_024306.4']['p_hgvs_slc'] == 'NP_077282.3:p.(R32H)'
        assert results['16-74808559-C-T']['hgvs_t_and_p']['NM_024306.4']['transcript_variant_error'] is None

    def test_variant112(self):
        variant = '16-89574804-C-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-89574804-C-A' in results.keys()
        assert results['16-89574804-C-A']['p_vcf'] == '16-89574804-C-A'
        assert results['16-89574804-C-A']['g_hgvs'] == 'NC_000016.9:g.89574804C>A'
        assert results['16-89574804-C-A']['genomic_variant_error'] is None
        assert 'NM_003119.3' in results['16-89574804-C-A']['hgvs_t_and_p'].keys()
        assert results['16-89574804-C-A']['hgvs_t_and_p']['NM_003119.3']['t_hgvs'] == 'NM_003119.3:c.-22C>A'
        assert results['16-89574804-C-A']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_tlc'] == 'NP_003110.1:p.?'
        assert results['16-89574804-C-A']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_slc'] == 'NP_003110.1:p.?'
        assert results['16-89574804-C-A']['hgvs_t_and_p']['NM_003119.3']['transcript_variant_error'] is None
        assert 'NM_001363850.1' in results['16-89574804-C-A']['hgvs_t_and_p'].keys()
        assert results['16-89574804-C-A']['hgvs_t_and_p']['NM_001363850.1']['t_hgvs'] == 'NM_001363850.1:c.-22C>A'
        assert results['16-89574804-C-A']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_tlc'] == 'NP_001350779.1:p.?'
        assert results['16-89574804-C-A']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_slc'] == 'NP_001350779.1:p.?'
        assert results['16-89574804-C-A']['hgvs_t_and_p']['NM_001363850.1']['transcript_variant_error'] is None
        assert 'NM_199367.2' in results['16-89574804-C-A']['hgvs_t_and_p'].keys()
        assert results['16-89574804-C-A']['hgvs_t_and_p']['NM_199367.2']['t_hgvs'] == 'NM_199367.2:c.-22C>A'
        assert results['16-89574804-C-A']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_tlc'] == 'NP_955399.1:p.?'
        assert results['16-89574804-C-A']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_slc'] == 'NP_955399.1:p.?'
        assert results['16-89574804-C-A']['hgvs_t_and_p']['NM_199367.2']['transcript_variant_error'] is None

    def test_variant113(self):
        variant = '16-89574826-A-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-89574826-A-C' in results.keys()
        assert results['16-89574826-A-C']['p_vcf'] == '16-89574826-A-C'
        assert results['16-89574826-A-C']['g_hgvs'] == 'NC_000016.9:g.89574826A>C'
        assert results['16-89574826-A-C']['genomic_variant_error'] is None
        assert 'NM_003119.3' in results['16-89574826-A-C']['hgvs_t_and_p'].keys()
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_003119.3']['t_hgvs'] == 'NM_003119.3:c.1A>C'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_tlc'] == 'NP_003110.1:p.(Met1?)'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_slc'] == 'NP_003110.1:p.(M1?)'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_003119.3']['transcript_variant_error'] is None
        assert 'NM_001363850.1' in results['16-89574826-A-C']['hgvs_t_and_p'].keys()
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_001363850.1']['t_hgvs'] == 'NM_001363850.1:c.1A>C'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_tlc'] == 'NP_001350779.1:p.(Met1?)'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_slc'] == 'NP_001350779.1:p.(M1?)'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_001363850.1']['transcript_variant_error'] is None
        assert 'NM_003119.2' in results['16-89574826-A-C']['hgvs_t_and_p'].keys()
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_003119.2']['t_hgvs'] == 'NM_003119.2:c.1A>C'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_tlc'] == 'NP_003110.1:p.(Met1?)'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_slc'] == 'NP_003110.1:p.(M1?)'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_003119.2']['transcript_variant_error'] is None
        assert 'NM_199367.1' in results['16-89574826-A-C']['hgvs_t_and_p'].keys()
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_199367.1']['t_hgvs'] == 'NM_199367.1:c.1A>C'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_199367.1']['p_hgvs_tlc'] == 'NP_955399.1:p.(Met1?)'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_199367.1']['p_hgvs_slc'] == 'NP_955399.1:p.(M1?)'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_199367.1']['transcript_variant_error'] is None
        assert 'NM_199367.2' in results['16-89574826-A-C']['hgvs_t_and_p'].keys()
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_199367.2']['t_hgvs'] == 'NM_199367.2:c.1A>C'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_tlc'] == 'NP_955399.1:p.(Met1?)'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_slc'] == 'NP_955399.1:p.(M1?)'
        assert results['16-89574826-A-C']['hgvs_t_and_p']['NM_199367.2']['transcript_variant_error'] is None

    def test_variant114(self):
        variant = '16-89574914-G-GT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-89574914-G-GT' in results.keys()
        assert results['16-89574914-G-GT']['p_vcf'] == '16-89574914-G-GT'
        assert results['16-89574914-G-GT']['g_hgvs'] == 'NC_000016.9:g.89574915dup'
        assert results['16-89574914-G-GT']['genomic_variant_error'] is None
        assert 'NM_003119.3' in results['16-89574914-G-GT']['hgvs_t_and_p'].keys()
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_003119.3']['t_hgvs'] == 'NM_003119.3:c.90dup'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_tlc'] == 'NP_003110.1:p.(Pro31SerfsTer43)'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_slc'] == 'NP_003110.1:p.(P31Sfs*43)'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_003119.3']['transcript_variant_error'] is None
        assert 'NM_001363850.1' in results['16-89574914-G-GT']['hgvs_t_and_p'].keys()
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_001363850.1']['t_hgvs'] == 'NM_001363850.1:c.90dup'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_tlc'] == 'NP_001350779.1:p.(Pro31SerfsTer43)'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_slc'] == 'NP_001350779.1:p.(P31Sfs*43)'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_001363850.1']['transcript_variant_error'] is None
        assert 'NM_003119.2' in results['16-89574914-G-GT']['hgvs_t_and_p'].keys()
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_003119.2']['t_hgvs'] == 'NM_003119.2:c.90dup'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_tlc'] == 'NP_003110.1:p.(Pro31SerfsTer43)'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_slc'] == 'NP_003110.1:p.(P31Sfs*43)'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_003119.2']['transcript_variant_error'] is None
        assert 'NM_199367.1' in results['16-89574914-G-GT']['hgvs_t_and_p'].keys()
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_199367.1']['t_hgvs'] == 'NM_199367.1:c.90dup'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_199367.1']['p_hgvs_tlc'] == 'NP_955399.1:p.(Pro31SerfsTer43)'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_199367.1']['p_hgvs_slc'] == 'NP_955399.1:p.(P31Sfs*43)'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_199367.1']['transcript_variant_error'] is None
        assert 'NM_199367.2' in results['16-89574914-G-GT']['hgvs_t_and_p'].keys()
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_199367.2']['t_hgvs'] == 'NM_199367.2:c.90dup'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_tlc'] == 'NP_955399.1:p.(Pro31SerfsTer43)'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_slc'] == 'NP_955399.1:p.(P31Sfs*43)'
        assert results['16-89574914-G-GT']['hgvs_t_and_p']['NM_199367.2']['transcript_variant_error'] is None

    def test_variant115(self):
        variant = '16-89574916-C-CGTC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-89574916-C-CGTC' in results.keys()
        assert results['16-89574916-C-CGTC']['p_vcf'] == '16-89574913-A-AGTC'
        assert results['16-89574916-C-CGTC']['g_hgvs'] == 'NC_000016.9:g.89574914_89574916dup'
        assert results['16-89574916-C-CGTC']['genomic_variant_error'] is None
        assert 'NM_003119.3' in results['16-89574916-C-CGTC']['hgvs_t_and_p'].keys()
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_003119.3']['t_hgvs'] == 'NM_003119.3:c.89_91dup'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_tlc'] == 'NP_003110.1:p.(Ser30_Pro31insArg)'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_slc'] == 'NP_003110.1:p.(S30_P31insR)'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_003119.3']['transcript_variant_error'] is None
        assert 'NM_001363850.1' in results['16-89574916-C-CGTC']['hgvs_t_and_p'].keys()
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_001363850.1']['t_hgvs'] == 'NM_001363850.1:c.89_91dup'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_tlc'] == 'NP_001350779.1:p.(Ser30_Pro31insArg)'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_slc'] == 'NP_001350779.1:p.(S30_P31insR)'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_001363850.1']['transcript_variant_error'] is None
        assert 'NM_003119.2' in results['16-89574916-C-CGTC']['hgvs_t_and_p'].keys()
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_003119.2']['t_hgvs'] == 'NM_003119.2:c.89_91dup'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_tlc'] == 'NP_003110.1:p.(Ser30_Pro31insArg)'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_slc'] == 'NP_003110.1:p.(S30_P31insR)'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_003119.2']['transcript_variant_error'] is None
        assert 'NM_199367.1' in results['16-89574916-C-CGTC']['hgvs_t_and_p'].keys()
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_199367.1']['t_hgvs'] == 'NM_199367.1:c.89_91dup'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_199367.1']['p_hgvs_tlc'] == 'NP_955399.1:p.(Ser30_Pro31insArg)'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_199367.1']['p_hgvs_slc'] == 'NP_955399.1:p.(S30_P31insR)'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_199367.1']['transcript_variant_error'] is None
        assert 'NM_199367.2' in results['16-89574916-C-CGTC']['hgvs_t_and_p'].keys()
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_199367.2']['t_hgvs'] == 'NM_199367.2:c.89_91dup'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_tlc'] == 'NP_955399.1:p.(Ser30_Pro31insArg)'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_slc'] == 'NP_955399.1:p.(S30_P31insR)'
        assert results['16-89574916-C-CGTC']['hgvs_t_and_p']['NM_199367.2']['transcript_variant_error'] is None

    def test_variant116(self):
        variant = '16-89575009-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-89575009-G-A' in results.keys()
        assert results['16-89575009-G-A']['p_vcf'] == '16-89575009-G-A'
        assert results['16-89575009-G-A']['g_hgvs'] == 'NC_000016.9:g.89575009G>A'
        assert results['16-89575009-G-A']['genomic_variant_error'] is None
        assert 'NM_003119.3' in results['16-89575009-G-A']['hgvs_t_and_p'].keys()
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_003119.3']['t_hgvs'] == 'NM_003119.3:c.183+1G>A'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_tlc'] == 'NP_003110.1:p.?'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_slc'] == 'NP_003110.1:p.?'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_003119.3']['transcript_variant_error'] is None
        assert 'NM_001363850.1' in results['16-89575009-G-A']['hgvs_t_and_p'].keys()
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_001363850.1']['t_hgvs'] == 'NM_001363850.1:c.183+1G>A'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_tlc'] == 'NP_001350779.1:p.?'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_slc'] == 'NP_001350779.1:p.?'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_001363850.1']['transcript_variant_error'] is None
        assert 'NM_003119.2' in results['16-89575009-G-A']['hgvs_t_and_p'].keys()
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_003119.2']['t_hgvs'] == 'NM_003119.2:c.183+1G>A'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_tlc'] == 'NP_003110.1:p.?'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_slc'] == 'NP_003110.1:p.?'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_003119.2']['transcript_variant_error'] is None
        assert 'NM_199367.1' in results['16-89575009-G-A']['hgvs_t_and_p'].keys()
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_199367.1']['t_hgvs'] == 'NM_199367.1:c.183+1G>A'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_199367.1']['p_hgvs_tlc'] == 'NP_955399.1:p.?'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_199367.1']['p_hgvs_slc'] == 'NP_955399.1:p.?'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_199367.1']['transcript_variant_error'] is None
        assert 'NM_199367.2' in results['16-89575009-G-A']['hgvs_t_and_p'].keys()
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_199367.2']['t_hgvs'] == 'NM_199367.2:c.183+1G>A'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_tlc'] == 'NP_955399.1:p.?'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_slc'] == 'NP_955399.1:p.?'
        assert results['16-89575009-G-A']['hgvs_t_and_p']['NM_199367.2']['transcript_variant_error'] is None

    def test_variant117(self):
        variant = '16-89576896-A-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-89576896-A-C' in results.keys()
        assert results['16-89576896-A-C']['p_vcf'] == '16-89576896-A-C'
        assert results['16-89576896-A-C']['g_hgvs'] == 'NC_000016.9:g.89576896A>C'
        assert results['16-89576896-A-C']['genomic_variant_error'] is None
        assert 'NM_003119.3' in results['16-89576896-A-C']['hgvs_t_and_p'].keys()
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_003119.3']['t_hgvs'] == 'NM_003119.3:c.184-2A>C'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_tlc'] == 'NP_003110.1:p.?'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_slc'] == 'NP_003110.1:p.?'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_003119.3']['transcript_variant_error'] is None
        assert 'NM_001363850.1' in results['16-89576896-A-C']['hgvs_t_and_p'].keys()
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_001363850.1']['t_hgvs'] == 'NM_001363850.1:c.184-2A>C'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_tlc'] == 'NP_001350779.1:p.?'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_slc'] == 'NP_001350779.1:p.?'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_001363850.1']['transcript_variant_error'] is None
        assert 'NM_003119.2' in results['16-89576896-A-C']['hgvs_t_and_p'].keys()
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_003119.2']['t_hgvs'] == 'NM_003119.2:c.184-2A>C'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_tlc'] == 'NP_003110.1:p.?'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_slc'] == 'NP_003110.1:p.?'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_003119.2']['transcript_variant_error'] is None
        assert 'NM_199367.1' in results['16-89576896-A-C']['hgvs_t_and_p'].keys()
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_199367.1']['t_hgvs'] == 'NM_199367.1:c.184-2A>C'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_199367.1']['p_hgvs_tlc'] == 'NP_955399.1:p.?'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_199367.1']['p_hgvs_slc'] == 'NP_955399.1:p.?'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_199367.1']['transcript_variant_error'] is None
        assert 'NM_199367.2' in results['16-89576896-A-C']['hgvs_t_and_p'].keys()
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_199367.2']['t_hgvs'] == 'NM_199367.2:c.184-2A>C'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_tlc'] == 'NP_955399.1:p.?'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_slc'] == 'NP_955399.1:p.?'
        assert results['16-89576896-A-C']['hgvs_t_and_p']['NM_199367.2']['transcript_variant_error'] is None

    def test_variant118(self):
        variant = '16-89576931-G-GTG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-89576931-G-GTG' in results.keys()
        assert results['16-89576931-G-GTG']['p_vcf'] == '16-89576929-T-TTG'
        assert results['16-89576931-G-GTG']['g_hgvs'] == 'NC_000016.9:g.89576930_89576931dup'
        assert results['16-89576931-G-GTG']['genomic_variant_error'] is None
        assert 'NM_003119.3' in results['16-89576931-G-GTG']['hgvs_t_and_p'].keys()
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_003119.3']['t_hgvs'] == 'NM_003119.3:c.216_217dup'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_tlc'] == 'NP_003110.1:p.(Glu73ValfsTer9)'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_slc'] == 'NP_003110.1:p.(E73Vfs*9)'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_003119.3']['transcript_variant_error'] is None
        assert 'NM_001363850.1' in results['16-89576931-G-GTG']['hgvs_t_and_p'].keys()
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_001363850.1']['t_hgvs'] == 'NM_001363850.1:c.216_217dup'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_tlc'] == 'NP_001350779.1:p.(Glu73ValfsTer9)'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_slc'] == 'NP_001350779.1:p.(E73Vfs*9)'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_001363850.1']['transcript_variant_error'] is None
        assert 'NM_003119.2' in results['16-89576931-G-GTG']['hgvs_t_and_p'].keys()
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_003119.2']['t_hgvs'] == 'NM_003119.2:c.216_217dup'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_tlc'] == 'NP_003110.1:p.(Glu73ValfsTer9)'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_slc'] == 'NP_003110.1:p.(E73Vfs*9)'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_003119.2']['transcript_variant_error'] is None
        assert 'NM_199367.1' in results['16-89576931-G-GTG']['hgvs_t_and_p'].keys()
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_199367.1']['t_hgvs'] == 'NM_199367.1:c.216_217dup'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_199367.1']['p_hgvs_tlc'] == 'NP_955399.1:p.(Glu73ValfsTer9)'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_199367.1']['p_hgvs_slc'] == 'NP_955399.1:p.(E73Vfs*9)'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_199367.1']['transcript_variant_error'] is None
        assert 'NM_199367.2' in results['16-89576931-G-GTG']['hgvs_t_and_p'].keys()
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_199367.2']['t_hgvs'] == 'NM_199367.2:c.216_217dup'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_tlc'] == 'NP_955399.1:p.(Glu73ValfsTer9)'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_slc'] == 'NP_955399.1:p.(E73Vfs*9)'
        assert results['16-89576931-G-GTG']['hgvs_t_and_p']['NM_199367.2']['transcript_variant_error'] is None

    def test_variant119(self):
        variant = '16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C' in results.keys()
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['p_vcf'] == '16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['g_hgvs'] == 'NC_000016.9:g.89598370_89598395del'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['genomic_variant_error'] is None
        assert 'NM_003119.3' in results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p'].keys()
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_003119.3']['t_hgvs'] == 'NM_003119.3:c.1046_1071del'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_tlc'] == 'NP_003110.1:p.(Gly349AlafsTer38)'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_slc'] == 'NP_003110.1:p.(G349Afs*38)'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_003119.3']['transcript_variant_error'] is None
        assert 'NM_001363850.1' in results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p'].keys()
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_001363850.1']['t_hgvs'] == 'NM_001363850.1:c.1046_1071del'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_tlc'] == 'NP_001350779.1:p.(Gly349AlafsTer38)'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_slc'] == 'NP_001350779.1:p.(G349Afs*38)'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_001363850.1']['transcript_variant_error'] is None
        assert 'NM_003119.2' in results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p'].keys()
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_003119.2']['t_hgvs'] == 'NM_003119.2:c.1046_1071del'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_tlc'] == 'NP_003110.1:p.(Gly349AlafsTer38)'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_slc'] == 'NP_003110.1:p.(G349Afs*38)'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_003119.2']['transcript_variant_error'] is None
        assert 'NM_199367.1' in results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p'].keys()
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_199367.1']['t_hgvs'] == 'NM_199367.1:c.1046_1071del'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_199367.1']['p_hgvs_tlc'] == 'NP_955399.1:p.(Gly349AlafsTer38)'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_199367.1']['p_hgvs_slc'] == 'NP_955399.1:p.(G349Afs*38)'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_199367.1']['transcript_variant_error'] is None
        assert 'NM_199367.2' in results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p'].keys()
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_199367.2']['t_hgvs'] == 'NM_199367.2:c.1046_1071del'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_tlc'] == 'NP_955399.1:p.(Gly349AlafsTer38)'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_199367.2']['p_hgvs_slc'] == 'NP_955399.1:p.(G349Afs*38)'
        assert results['16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C']['hgvs_t_and_p']['NM_199367.2']['transcript_variant_error'] is None

    def test_variant120(self):
        variant = '16-89613064-AGGAGAGGCG-AT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-89613064-AGGAGAGGCG-AT' in results.keys()
        assert results['16-89613064-AGGAGAGGCG-AT']['p_vcf'] == '16-89613065-GGAGAGGCG-T'
        assert results['16-89613064-AGGAGAGGCG-AT']['g_hgvs'] == 'NC_000016.9:g.89613065_89613073delinsT'
        assert results['16-89613064-AGGAGAGGCG-AT']['genomic_variant_error'] is None
        assert 'NM_003119.3' in results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p'].keys()
        assert results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p']['NM_003119.3']['t_hgvs'] == 'NM_003119.3:c.1450-1_1457delinsT'
        assert results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_tlc'] == 'NP_003110.1:p.?'
        assert results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_slc'] == 'NP_003110.1:p.?'
        assert results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p']['NM_003119.3']['transcript_variant_error'] is None
        assert 'NM_001363850.1' in results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p'].keys()
        assert results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p']['NM_001363850.1']['t_hgvs'] == 'NM_001363850.1:c.1450-1_1457delinsT'
        assert results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_tlc'] == 'NP_001350779.1:p.?'
        assert results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_slc'] == 'NP_001350779.1:p.?'
        assert results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p']['NM_001363850.1']['transcript_variant_error'] is None
        assert 'NM_003119.2' in results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p'].keys()
        assert results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p']['NM_003119.2']['t_hgvs'] == 'NM_003119.2:c.1450-1_1457delinsT'
        assert results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_tlc'] == 'NP_003110.1:p.?'
        assert results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_slc'] == 'NP_003110.1:p.?'
        assert results['16-89613064-AGGAGAGGCG-AT']['hgvs_t_and_p']['NM_003119.2']['transcript_variant_error'] is None

    def test_variant121(self):
        variant = '16-89613069-AGGCGGGAGA-AT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-89613069-AGGCGGGAGA-AT' in results.keys()
        assert results['16-89613069-AGGCGGGAGA-AT']['p_vcf'] == '16-89613070-GGCGGGAGA-T'
        assert results['16-89613069-AGGCGGGAGA-AT']['g_hgvs'] == 'NC_000016.9:g.89613070_89613078delinsT'
        assert results['16-89613069-AGGCGGGAGA-AT']['genomic_variant_error'] is None
        assert 'NM_003119.3' in results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p'].keys()
        assert results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p']['NM_003119.3']['t_hgvs'] == 'NM_003119.3:c.1454_1462delinsT'
        assert results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_tlc'] == 'NP_003110.1:p.(Arg485IlefsTer3)'
        assert results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_slc'] == 'NP_003110.1:p.(R485Ifs*3)'
        assert results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p']['NM_003119.3']['transcript_variant_error'] is None
        assert 'NM_001363850.1' in results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p'].keys()
        assert results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p']['NM_001363850.1']['t_hgvs'] == 'NM_001363850.1:c.1454_1462delinsT'
        assert results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_tlc'] == 'NP_001350779.1:p.(Arg485IlefsTer3)'
        assert results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_slc'] == 'NP_001350779.1:p.(R485Ifs*3)'
        assert results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p']['NM_001363850.1']['transcript_variant_error'] is None
        assert 'NM_003119.2' in results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p'].keys()
        assert results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p']['NM_003119.2']['t_hgvs'] == 'NM_003119.2:c.1454_1462delinsT'
        assert results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_tlc'] == 'NP_003110.1:p.(Arg485IlefsTer3)'
        assert results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_slc'] == 'NP_003110.1:p.(R485Ifs*3)'
        assert results['16-89613069-AGGCGGGAGA-AT']['hgvs_t_and_p']['NM_003119.2']['transcript_variant_error'] is None

    def test_variant122(self):
        variant = '16-89613145-C-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '16-89613145-C-T' in results.keys()
        assert results['16-89613145-C-T']['p_vcf'] == '16-89613145-C-T'
        assert results['16-89613145-C-T']['g_hgvs'] == 'NC_000016.9:g.89613145C>T'
        assert results['16-89613145-C-T']['genomic_variant_error'] is None
        assert 'NM_003119.3' in results['16-89613145-C-T']['hgvs_t_and_p'].keys()
        assert results['16-89613145-C-T']['hgvs_t_and_p']['NM_003119.3']['t_hgvs'] == 'NM_003119.3:c.1529C>T'
        assert results['16-89613145-C-T']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_tlc'] == 'NP_003110.1:p.(Ala510Val)'
        assert results['16-89613145-C-T']['hgvs_t_and_p']['NM_003119.3']['p_hgvs_slc'] == 'NP_003110.1:p.(A510V)'
        assert results['16-89613145-C-T']['hgvs_t_and_p']['NM_003119.3']['transcript_variant_error'] is None
        assert 'NM_001363850.1' in results['16-89613145-C-T']['hgvs_t_and_p'].keys()
        assert results['16-89613145-C-T']['hgvs_t_and_p']['NM_001363850.1']['t_hgvs'] == 'NM_001363850.1:c.1529C>T'
        assert results['16-89613145-C-T']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_tlc'] == 'NP_001350779.1:p.(Ala510Val)'
        assert results['16-89613145-C-T']['hgvs_t_and_p']['NM_001363850.1']['p_hgvs_slc'] == 'NP_001350779.1:p.(A510V)'
        assert results['16-89613145-C-T']['hgvs_t_and_p']['NM_001363850.1']['transcript_variant_error'] is None
        assert 'NM_003119.2' in results['16-89613145-C-T']['hgvs_t_and_p'].keys()
        assert results['16-89613145-C-T']['hgvs_t_and_p']['NM_003119.2']['t_hgvs'] == 'NM_003119.2:c.1529C>T'
        assert results['16-89613145-C-T']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_tlc'] == 'NP_003110.1:p.(Ala510Val)'
        assert results['16-89613145-C-T']['hgvs_t_and_p']['NM_003119.2']['p_hgvs_slc'] == 'NP_003110.1:p.(A510V)'
        assert results['16-89613145-C-T']['hgvs_t_and_p']['NM_003119.2']['transcript_variant_error'] is None

    def test_variant123(self):
        variant = '17-7578194-GCAC-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-7578194-GCAC-G' in results.keys()
        assert results['17-7578194-GCAC-G']['p_vcf'] == '17-7578194-GCAC-G'
        assert results['17-7578194-GCAC-G']['g_hgvs'] == 'NC_000017.10:g.7578201_7578203del'
        assert results['17-7578194-GCAC-G']['genomic_variant_error'] is None
        assert 'NM_001276697.1' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276697.1']['t_hgvs'] == 'NM_001276697.1:c.175_177del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276697.1']['p_hgvs_tlc'] == 'NP_001263626.1:p.(Val59del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276697.1']['p_hgvs_slc'] == 'NP_001263626.1:p.(V59del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276697.1']['transcript_variant_error'] is None
        assert 'NM_001276699.1' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276699.1']['t_hgvs'] == 'NM_001276699.1:c.175_177del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276699.1']['p_hgvs_tlc'] == 'NP_001263628.1:p.(Val59del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276699.1']['p_hgvs_slc'] == 'NP_001263628.1:p.(V59del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276699.1']['transcript_variant_error'] is None
        assert 'NM_001276696.1' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276696.1']['t_hgvs'] == 'NM_001276696.1:c.535_537del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276696.1']['p_hgvs_tlc'] == 'NP_001263625.1:p.(Val179del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276696.1']['p_hgvs_slc'] == 'NP_001263625.1:p.(V179del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276696.1']['transcript_variant_error'] is None
        assert 'NM_001276760.1' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276760.1']['t_hgvs'] == 'NM_001276760.1:c.535_537del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276760.1']['p_hgvs_tlc'] == 'NP_001263689.1:p.(Val179del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276760.1']['p_hgvs_slc'] == 'NP_001263689.1:p.(V179del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276760.1']['transcript_variant_error'] is None
        assert 'NM_001126116.1' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126116.1']['t_hgvs'] == 'NM_001126116.1:c.256_258del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126116.1']['p_hgvs_tlc'] == 'NP_001119588.1:p.(Val86del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126116.1']['p_hgvs_slc'] == 'NP_001119588.1:p.(V86del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126116.1']['transcript_variant_error'] is None
        assert 'NM_001276761.1' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276761.1']['t_hgvs'] == 'NM_001276761.1:c.535_537del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276761.1']['p_hgvs_tlc'] == 'NP_001263690.1:p.(Val179del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276761.1']['p_hgvs_slc'] == 'NP_001263690.1:p.(V179del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276761.1']['transcript_variant_error'] is None
        assert 'NM_001126113.2' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126113.2']['t_hgvs'] == 'NM_001126113.2:c.652_654del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126113.2']['p_hgvs_tlc'] == 'NP_001119585.1:p.(Val218del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126113.2']['p_hgvs_slc'] == 'NP_001119585.1:p.(V218del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126113.2']['transcript_variant_error'] is None
        assert 'NM_001126117.1' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126117.1']['t_hgvs'] == 'NM_001126117.1:c.256_258del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126117.1']['p_hgvs_tlc'] == 'NP_001119589.1:p.(Val86del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126117.1']['p_hgvs_slc'] == 'NP_001119589.1:p.(V86del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126117.1']['transcript_variant_error'] is None
        assert 'NM_001276695.1' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276695.1']['t_hgvs'] == 'NM_001276695.1:c.535_537del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276695.1']['p_hgvs_tlc'] == 'NP_001263624.1:p.(Val179del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276695.1']['p_hgvs_slc'] == 'NP_001263624.1:p.(V179del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276695.1']['transcript_variant_error'] is None
        assert 'NM_001276698.1' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276698.1']['t_hgvs'] == 'NM_001276698.1:c.175_177del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276698.1']['p_hgvs_tlc'] == 'NP_001263627.1:p.(Val59del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276698.1']['p_hgvs_slc'] == 'NP_001263627.1:p.(V59del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001276698.1']['transcript_variant_error'] is None
        assert 'NM_001126112.2' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126112.2']['t_hgvs'] == 'NM_001126112.2:c.652_654del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126112.2']['p_hgvs_tlc'] == 'NP_001119584.1:p.(Val218del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126112.2']['p_hgvs_slc'] == 'NP_001119584.1:p.(V218del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126112.2']['transcript_variant_error'] is None
        assert 'NM_001126114.2' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126114.2']['t_hgvs'] == 'NM_001126114.2:c.652_654del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126114.2']['p_hgvs_tlc'] == 'NP_001119586.1:p.(Val218del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126114.2']['p_hgvs_slc'] == 'NP_001119586.1:p.(V218del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126114.2']['transcript_variant_error'] is None
        assert 'NM_000546.5' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_000546.5']['t_hgvs'] == 'NM_000546.5:c.652_654del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_000546.5']['p_hgvs_tlc'] == 'NP_000537.3:p.(Val218del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_000546.5']['p_hgvs_slc'] == 'NP_000537.3:p.(V218del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_000546.5']['transcript_variant_error'] is None
        assert 'NM_001126115.1' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126115.1']['t_hgvs'] == 'NM_001126115.1:c.256_258del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126115.1']['p_hgvs_tlc'] == 'NP_001119587.1:p.(Val86del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126115.1']['p_hgvs_slc'] == 'NP_001119587.1:p.(V86del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126115.1']['transcript_variant_error'] is None
        assert 'NM_001126118.1' in results['17-7578194-GCAC-G']['hgvs_t_and_p'].keys()
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126118.1']['t_hgvs'] == 'NM_001126118.1:c.535_537del'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126118.1']['p_hgvs_tlc'] == 'NP_001119590.1:p.(Val179del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126118.1']['p_hgvs_slc'] == 'NP_001119590.1:p.(V179del)'
        assert results['17-7578194-GCAC-G']['hgvs_t_and_p']['NM_001126118.1']['transcript_variant_error'] is None

    def test_variant124(self):
        variant = '17-7578523-T-TG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-7578523-T-TG' in results.keys()
        assert results['17-7578523-T-TG']['p_vcf'] == '17-7578523-T-TG'
        assert results['17-7578523-T-TG']['g_hgvs'] == 'NC_000017.10:g.7578525dup'
        assert results['17-7578523-T-TG']['genomic_variant_error'] is None
        assert 'NM_001276697.1' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276697.1']['t_hgvs'] == 'NM_001276697.1:c.-72dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276697.1']['p_hgvs_tlc'] == 'NP_001263626.1:p.?'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276697.1']['p_hgvs_slc'] == 'NP_001263626.1:p.?'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276697.1']['transcript_variant_error'] is None
        assert 'NM_001276699.1' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276699.1']['t_hgvs'] == 'NM_001276699.1:c.-72dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276699.1']['p_hgvs_tlc'] == 'NP_001263628.1:p.?'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276699.1']['p_hgvs_slc'] == 'NP_001263628.1:p.?'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276699.1']['transcript_variant_error'] is None
        assert 'NM_001276696.1' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276696.1']['t_hgvs'] == 'NM_001276696.1:c.289dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276696.1']['p_hgvs_tlc'] == 'NP_001263625.1:p.(Gln97ProfsTer13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276696.1']['p_hgvs_slc'] == 'NP_001263625.1:p.(Q97Pfs*13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276696.1']['transcript_variant_error'] is None
        assert 'NM_001276760.1' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276760.1']['t_hgvs'] == 'NM_001276760.1:c.289dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276760.1']['p_hgvs_tlc'] == 'NP_001263689.1:p.(Gln97ProfsTer13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276760.1']['p_hgvs_slc'] == 'NP_001263689.1:p.(Q97Pfs*13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276760.1']['transcript_variant_error'] is None
        assert 'NM_001126116.1' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126116.1']['t_hgvs'] == 'NM_001126116.1:c.10dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126116.1']['p_hgvs_tlc'] == 'NP_001119588.1:p.(Gln4ProfsTer13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126116.1']['p_hgvs_slc'] == 'NP_001119588.1:p.(Q4Pfs*13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126116.1']['transcript_variant_error'] is None
        assert 'NM_001276761.1' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276761.1']['t_hgvs'] == 'NM_001276761.1:c.289dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276761.1']['p_hgvs_tlc'] == 'NP_001263690.1:p.(Gln97ProfsTer13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276761.1']['p_hgvs_slc'] == 'NP_001263690.1:p.(Q97Pfs*13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276761.1']['transcript_variant_error'] is None
        assert 'NM_001126113.2' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126113.2']['t_hgvs'] == 'NM_001126113.2:c.406dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126113.2']['p_hgvs_tlc'] == 'NP_001119585.1:p.(Gln136ProfsTer13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126113.2']['p_hgvs_slc'] == 'NP_001119585.1:p.(Q136Pfs*13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126113.2']['transcript_variant_error'] is None
        assert 'NM_001126117.1' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126117.1']['t_hgvs'] == 'NM_001126117.1:c.10dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126117.1']['p_hgvs_tlc'] == 'NP_001119589.1:p.(Gln4ProfsTer13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126117.1']['p_hgvs_slc'] == 'NP_001119589.1:p.(Q4Pfs*13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126117.1']['transcript_variant_error'] is None
        assert 'NM_001276695.1' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276695.1']['t_hgvs'] == 'NM_001276695.1:c.289dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276695.1']['p_hgvs_tlc'] == 'NP_001263624.1:p.(Gln97ProfsTer13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276695.1']['p_hgvs_slc'] == 'NP_001263624.1:p.(Q97Pfs*13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276695.1']['transcript_variant_error'] is None
        assert 'NM_001276698.1' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276698.1']['t_hgvs'] == 'NM_001276698.1:c.-72dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276698.1']['p_hgvs_tlc'] == 'NP_001263627.1:p.?'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276698.1']['p_hgvs_slc'] == 'NP_001263627.1:p.?'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001276698.1']['transcript_variant_error'] is None
        assert 'NM_001126112.2' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126112.2']['t_hgvs'] == 'NM_001126112.2:c.406dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126112.2']['p_hgvs_tlc'] == 'NP_001119584.1:p.(Gln136ProfsTer13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126112.2']['p_hgvs_slc'] == 'NP_001119584.1:p.(Q136Pfs*13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126112.2']['transcript_variant_error'] is None
        assert 'NM_001126114.2' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126114.2']['t_hgvs'] == 'NM_001126114.2:c.406dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126114.2']['p_hgvs_tlc'] == 'NP_001119586.1:p.(Gln136ProfsTer13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126114.2']['p_hgvs_slc'] == 'NP_001119586.1:p.(Q136Pfs*13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126114.2']['transcript_variant_error'] is None
        assert 'NM_000546.5' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_000546.5']['t_hgvs'] == 'NM_000546.5:c.406dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_000546.5']['p_hgvs_tlc'] == 'NP_000537.3:p.(Gln136ProfsTer13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_000546.5']['p_hgvs_slc'] == 'NP_000537.3:p.(Q136Pfs*13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_000546.5']['transcript_variant_error'] is None
        assert 'NM_001126115.1' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126115.1']['t_hgvs'] == 'NM_001126115.1:c.10dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126115.1']['p_hgvs_tlc'] == 'NP_001119587.1:p.(Gln4ProfsTer13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126115.1']['p_hgvs_slc'] == 'NP_001119587.1:p.(Q4Pfs*13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126115.1']['transcript_variant_error'] is None
        assert 'NM_001126118.1' in results['17-7578523-T-TG']['hgvs_t_and_p'].keys()
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126118.1']['t_hgvs'] == 'NM_001126118.1:c.289dup'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126118.1']['p_hgvs_tlc'] == 'NP_001119590.1:p.(Gln97ProfsTer13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126118.1']['p_hgvs_slc'] == 'NP_001119590.1:p.(Q97Pfs*13)'
        assert results['17-7578523-T-TG']['hgvs_t_and_p']['NM_001126118.1']['transcript_variant_error'] is None

    def test_variant125(self):
        variant = '17-17119692-A-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-17119692-A-C' in results.keys()
        assert results['17-17119692-A-C']['p_vcf'] == '17-17119692-A-C'
        assert results['17-17119692-A-C']['g_hgvs'] == 'NC_000017.10:g.17119692A>C'
        assert results['17-17119692-A-C']['genomic_variant_error'] is None
        assert 'NM_001353230.1' in results['17-17119692-A-C']['hgvs_t_and_p'].keys()
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_001353230.1']['t_hgvs'] == 'NM_001353230.1:c.1300+2T>G'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_001353230.1']['p_hgvs_tlc'] == 'NP_001340159.1:p.?'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_001353230.1']['p_hgvs_slc'] == 'NP_001340159.1:p.?'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_001353230.1']['transcript_variant_error'] is None
        assert 'NM_144997.5' in results['17-17119692-A-C']['hgvs_t_and_p'].keys()
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_144997.5']['t_hgvs'] == 'NM_144997.5:c.1300+2T>G'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_144997.5']['p_hgvs_tlc'] == 'NP_659434.2:p.?'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_144997.5']['p_hgvs_slc'] == 'NP_659434.2:p.?'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_144997.5']['transcript_variant_error'] is None
        assert 'NM_001353231.1' in results['17-17119692-A-C']['hgvs_t_and_p'].keys()
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_001353231.1']['t_hgvs'] == 'NM_001353231.1:c.1300+2T>G'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_001353231.1']['p_hgvs_tlc'] == 'NP_001340160.1:p.?'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_001353231.1']['p_hgvs_slc'] == 'NP_001340160.1:p.?'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_001353231.1']['transcript_variant_error'] is None
        assert 'NM_001353229.1' in results['17-17119692-A-C']['hgvs_t_and_p'].keys()
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_001353229.1']['t_hgvs'] == 'NM_001353229.1:c.1354+2T>G'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_001353229.1']['p_hgvs_tlc'] == 'NP_001340158.1:p.?'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_001353229.1']['p_hgvs_slc'] == 'NP_001340158.1:p.?'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_001353229.1']['transcript_variant_error'] is None
        assert 'NM_144997.6' in results['17-17119692-A-C']['hgvs_t_and_p'].keys()
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_144997.6']['t_hgvs'] == 'NM_144997.6:c.1300+2T>G'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_144997.6']['p_hgvs_tlc'] == 'NP_659434.2:p.?'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_144997.6']['p_hgvs_slc'] == 'NP_659434.2:p.?'
        assert results['17-17119692-A-C']['hgvs_t_and_p']['NM_144997.6']['transcript_variant_error'] is None

    def test_variant126(self):
        variant = '17-41197588-GGACA-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-41197588-GGACA-G' in results.keys()
        assert results['17-41197588-GGACA-G']['p_vcf'] == '17-41197588-GGACA-G'
        assert results['17-41197588-GGACA-G']['g_hgvs'] == 'NC_000017.10:g.41197590_41197593del'
        assert results['17-41197588-GGACA-G']['genomic_variant_error'] is None
        assert 'NM_007294.3' in results['17-41197588-GGACA-G']['hgvs_t_and_p'].keys()
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007294.3']['t_hgvs'] == 'NM_007294.3:c.*103_*106del'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007294.3']['p_hgvs_tlc'] == 'NP_009225.1:p.?'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007294.3']['p_hgvs_slc'] == 'NP_009225.1:p.?'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007294.3']['transcript_variant_error'] is None
        assert 'NM_007300.3' in results['17-41197588-GGACA-G']['hgvs_t_and_p'].keys()
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007300.3']['t_hgvs'] == 'NM_007300.3:c.*103_*106del'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007300.3']['p_hgvs_tlc'] == 'NP_009231.2:p.?'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007300.3']['p_hgvs_slc'] == 'NP_009231.2:p.?'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007300.3']['transcript_variant_error'] is None
        assert 'NM_007299.3' in results['17-41197588-GGACA-G']['hgvs_t_and_p'].keys()
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007299.3']['t_hgvs'] == 'NM_007299.3:c.*209_*212del'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007299.3']['p_hgvs_tlc'] == 'NP_009230.2:p.?'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007299.3']['p_hgvs_slc'] == 'NP_009230.2:p.?'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007299.3']['transcript_variant_error'] is None
        assert 'NM_007298.3' in results['17-41197588-GGACA-G']['hgvs_t_and_p'].keys()
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007298.3']['t_hgvs'] == 'NM_007298.3:c.*103_*106del'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007298.3']['p_hgvs_tlc'] == 'NP_009229.2:p.?'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007298.3']['p_hgvs_slc'] == 'NP_009229.2:p.?'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007298.3']['transcript_variant_error'] is None
        assert 'NR_027676.1' in results['17-41197588-GGACA-G']['hgvs_t_and_p'].keys()
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NR_027676.1']['t_hgvs'] == 'NR_027676.1:n.5831_5834del'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NR_027676.1']['p_hgvs_tlc'] is None
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NR_027676.1']['p_hgvs_slc'] is None
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NR_027676.1']['transcript_variant_error'] is None
        assert 'NM_007297.3' in results['17-41197588-GGACA-G']['hgvs_t_and_p'].keys()
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007297.3']['t_hgvs'] == 'NM_007297.3:c.*103_*106del'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007297.3']['p_hgvs_tlc'] == 'NP_009228.2:p.?'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007297.3']['p_hgvs_slc'] == 'NP_009228.2:p.?'
        assert results['17-41197588-GGACA-G']['hgvs_t_and_p']['NM_007297.3']['transcript_variant_error'] is None

    def test_variant127(self):
        variant = '17-41256884-C-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-41256884-C-G' in results.keys()
        assert results['17-41256884-C-G']['p_vcf'] == '17-41256884-C-G'
        assert results['17-41256884-C-G']['g_hgvs'] == 'NC_000017.10:g.41256884C>G'
        assert results['17-41256884-C-G']['genomic_variant_error'] is None
        assert 'NM_007294.3' in results['17-41256884-C-G']['hgvs_t_and_p'].keys()
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007294.3']['t_hgvs'] == 'NM_007294.3:c.301+1G>C'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007294.3']['p_hgvs_tlc'] == 'NP_009225.1:p.?'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007294.3']['p_hgvs_slc'] == 'NP_009225.1:p.?'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007294.3']['transcript_variant_error'] is None
        assert 'NM_007300.3' in results['17-41256884-C-G']['hgvs_t_and_p'].keys()
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007300.3']['t_hgvs'] == 'NM_007300.3:c.301+1G>C'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007300.3']['p_hgvs_tlc'] == 'NP_009231.2:p.?'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007300.3']['p_hgvs_slc'] == 'NP_009231.2:p.?'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007300.3']['transcript_variant_error'] is None
        assert 'NM_007299.3' in results['17-41256884-C-G']['hgvs_t_and_p'].keys()
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007299.3']['t_hgvs'] == 'NM_007299.3:c.301+1G>C'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007299.3']['p_hgvs_tlc'] == 'NP_009230.2:p.?'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007299.3']['p_hgvs_slc'] == 'NP_009230.2:p.?'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007299.3']['transcript_variant_error'] is None
        assert 'NM_007298.3' in results['17-41256884-C-G']['hgvs_t_and_p'].keys()
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007298.3']['t_hgvs'] == 'NM_007298.3:c.301+1G>C'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007298.3']['p_hgvs_tlc'] == 'NP_009229.2:p.?'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007298.3']['p_hgvs_slc'] == 'NP_009229.2:p.?'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007298.3']['transcript_variant_error'] is None
        assert 'NR_027676.1' in results['17-41256884-C-G']['hgvs_t_and_p'].keys()
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NR_027676.1']['t_hgvs'] == 'NR_027676.1:n.440+1G>C'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NR_027676.1']['p_hgvs_tlc'] is None
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NR_027676.1']['p_hgvs_slc'] is None
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NR_027676.1']['transcript_variant_error'] is None
        assert 'NM_007297.3' in results['17-41256884-C-G']['hgvs_t_and_p'].keys()
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007297.3']['t_hgvs'] == 'NM_007297.3:c.160+1G>C'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007297.3']['p_hgvs_tlc'] == 'NP_009228.2:p.?'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007297.3']['p_hgvs_slc'] == 'NP_009228.2:p.?'
        assert results['17-41256884-C-G']['hgvs_t_and_p']['NM_007297.3']['transcript_variant_error'] is None

    def test_variant128(self):
        variant = '17-42991428-C-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-42991428-C-A' in results.keys()
        assert results['17-42991428-C-A']['p_vcf'] == '17-42991428-C-A'
        assert results['17-42991428-C-A']['g_hgvs'] == 'NC_000017.10:g.42991428C>A'
        assert results['17-42991428-C-A']['genomic_variant_error'] is None
        assert 'NM_001131019.2' in results['17-42991428-C-A']['hgvs_t_and_p'].keys()
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_001131019.2']['t_hgvs'] == 'NM_001131019.2:c.490G>T'
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_001131019.2']['p_hgvs_tlc'] == 'NP_001124491.1:p.(Glu164Ter)'
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_001131019.2']['p_hgvs_slc'] == 'NP_001124491.1:p.(E164*)'
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_001131019.2']['transcript_variant_error'] is None
        assert 'NM_001242376.1' in results['17-42991428-C-A']['hgvs_t_and_p'].keys()
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_001242376.1']['t_hgvs'] == 'NM_001242376.1:c.490G>T'
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_001242376.1']['p_hgvs_tlc'] == 'NP_001229305.1:p.(Glu164Ter)'
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_001242376.1']['p_hgvs_slc'] == 'NP_001229305.1:p.(E164*)'
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_001242376.1']['transcript_variant_error'] is None
        assert 'NM_002055.4' in results['17-42991428-C-A']['hgvs_t_and_p'].keys()
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_002055.4']['t_hgvs'] == 'NM_002055.4:c.490G>T'
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_002055.4']['p_hgvs_tlc'] == 'NP_002046.1:p.(Glu164Ter)'
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_002055.4']['p_hgvs_slc'] == 'NP_002046.1:p.(E164*)'
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_002055.4']['transcript_variant_error'] is None
        assert 'NM_001363846.1' in results['17-42991428-C-A']['hgvs_t_and_p'].keys()
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_001363846.1']['t_hgvs'] == 'NM_001363846.1:c.490G>T'
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_001363846.1']['p_hgvs_tlc'] == 'NP_001350775.1:p.(Glu164Ter)'
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_001363846.1']['p_hgvs_slc'] == 'NP_001350775.1:p.(E164*)'
        assert results['17-42991428-C-A']['hgvs_t_and_p']['NM_001363846.1']['transcript_variant_error'] is None

    def test_variant129(self):
        variant = '17-48252809-A-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-48252809-A-T' in results.keys()
        assert results['17-48252809-A-T']['p_vcf'] == '17-48252809-A-T'
        assert results['17-48252809-A-T']['g_hgvs'] == 'NC_000017.10:g.48252809A>T'
        assert results['17-48252809-A-T']['genomic_variant_error'] is None
        assert 'NM_001135697.2' in results['17-48252809-A-T']['hgvs_t_and_p'].keys()
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_001135697.2']['t_hgvs'] == 'NM_001135697.2:c.*11A>T'
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_001135697.2']['p_hgvs_tlc'] == 'NP_001129169.1:p.?'
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_001135697.2']['p_hgvs_slc'] == 'NP_001129169.1:p.?'
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_001135697.2']['transcript_variant_error'] is None
        assert 'NM_000023.2' in results['17-48252809-A-T']['hgvs_t_and_p'].keys()
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_000023.2']['t_hgvs'] == 'NM_000023.2:c.*11A>T'
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_000023.2']['p_hgvs_tlc'] == 'NP_000014.1:p.?'
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_000023.2']['p_hgvs_slc'] == 'NP_000014.1:p.?'
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_000023.2']['transcript_variant_error'] is None
        assert 'NM_001135697.1' in results['17-48252809-A-T']['hgvs_t_and_p'].keys()
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_001135697.1']['t_hgvs'] == 'NM_001135697.1:c.*11A>T'
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_001135697.1']['p_hgvs_tlc'] == 'NP_001129169.1:p.?'
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_001135697.1']['p_hgvs_slc'] == 'NP_001129169.1:p.?'
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_001135697.1']['transcript_variant_error'] is None
        assert 'NR_135553.1' in results['17-48252809-A-T']['hgvs_t_and_p'].keys()
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NR_135553.1']['t_hgvs'] == 'NR_135553.1:n.1022A>T'
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NR_135553.1']['p_hgvs_tlc'] is None
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NR_135553.1']['p_hgvs_slc'] is None
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NR_135553.1']['transcript_variant_error'] is None
        assert 'NM_000023.3' in results['17-48252809-A-T']['hgvs_t_and_p'].keys()
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_000023.3']['t_hgvs'] == 'NM_000023.3:c.*11A>T'
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_000023.3']['p_hgvs_tlc'] == 'NP_000014.1:p.?'
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_000023.3']['p_hgvs_slc'] == 'NP_000014.1:p.?'
        assert results['17-48252809-A-T']['hgvs_t_and_p']['NM_000023.3']['transcript_variant_error'] is None

    def test_variant130(self):
        variant = '17-62022709-G-GTC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-62022709-G-GTC' in results.keys()
        assert results['17-62022709-G-GTC']['p_vcf'] == '17-62022709-G-GTC'
        assert results['17-62022709-G-GTC']['g_hgvs'] == 'NC_000017.10:g.62022710_62022711dup'
        assert results['17-62022709-G-GTC']['genomic_variant_error'] is None
        assert 'NM_000334.4' in results['17-62022709-G-GTC']['hgvs_t_and_p'].keys()
        assert results['17-62022709-G-GTC']['hgvs_t_and_p']['NM_000334.4']['t_hgvs'] == 'NM_000334.4:c.3720+9_3720+10dup'
        assert results['17-62022709-G-GTC']['hgvs_t_and_p']['NM_000334.4']['p_hgvs_tlc'] == 'NP_000325.4:p.?'
        assert results['17-62022709-G-GTC']['hgvs_t_and_p']['NM_000334.4']['p_hgvs_slc'] == 'NP_000325.4:p.?'
        assert results['17-62022709-G-GTC']['hgvs_t_and_p']['NM_000334.4']['transcript_variant_error'] is None

    def test_variant131(self):
        variant = '17-62022711-C-CT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-62022711-C-CT' in results.keys()
        assert results['17-62022711-C-CT']['p_vcf'] == '17-62022711-C-CT'
        assert results['17-62022711-C-CT']['g_hgvs'] == 'NC_000017.10:g.62022711_62022712insT'
        assert results['17-62022711-C-CT']['genomic_variant_error'] is None
        assert 'NM_000334.4' in results['17-62022711-C-CT']['hgvs_t_and_p'].keys()
        assert results['17-62022711-C-CT']['hgvs_t_and_p']['NM_000334.4']['t_hgvs'] == 'NM_000334.4:c.3720+8_3720+9insA'
        assert results['17-62022711-C-CT']['hgvs_t_and_p']['NM_000334.4']['p_hgvs_tlc'] == 'NP_000325.4:p.?'
        assert results['17-62022711-C-CT']['hgvs_t_and_p']['NM_000334.4']['p_hgvs_slc'] == 'NP_000325.4:p.?'
        assert results['17-62022711-C-CT']['hgvs_t_and_p']['NM_000334.4']['transcript_variant_error'] is None

    def test_variant132(self):
        variant = '17-62023005-G-GGC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-62023005-G-GGC' in results.keys()
        assert results['17-62023005-G-GGC']['p_vcf'] == '17-62023005-G-GGC'
        assert results['17-62023005-G-GGC']['g_hgvs'] == 'NC_000017.10:g.62023005_62023006insGC'
        assert results['17-62023005-G-GGC']['genomic_variant_error'] is None
        assert 'NM_000334.4' in results['17-62023005-G-GGC']['hgvs_t_and_p'].keys()
        assert results['17-62023005-G-GGC']['hgvs_t_and_p']['NM_000334.4']['t_hgvs'] == 'NM_000334.4:c.3442-8_3442-7insGC'
        assert results['17-62023005-G-GGC']['hgvs_t_and_p']['NM_000334.4']['p_hgvs_tlc'] == 'NP_000325.4:p.?'
        assert results['17-62023005-G-GGC']['hgvs_t_and_p']['NM_000334.4']['p_hgvs_slc'] == 'NP_000325.4:p.?'
        assert results['17-62023005-G-GGC']['hgvs_t_and_p']['NM_000334.4']['transcript_variant_error'] is None

    def test_variant133(self):
        variant = '17-62023006-C-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-62023006-C-A' in results.keys()
        assert results['17-62023006-C-A']['p_vcf'] == '17-62023006-C-A'
        assert results['17-62023006-C-A']['g_hgvs'] == 'NC_000017.10:g.62023006C>A'
        assert results['17-62023006-C-A']['genomic_variant_error'] is None
        assert 'NM_000334.4' in results['17-62023006-C-A']['hgvs_t_and_p'].keys()
        assert results['17-62023006-C-A']['hgvs_t_and_p']['NM_000334.4']['t_hgvs'] == 'NM_000334.4:c.3442-8G>T'
        assert results['17-62023006-C-A']['hgvs_t_and_p']['NM_000334.4']['p_hgvs_tlc'] == 'NP_000325.4:p.?'
        assert results['17-62023006-C-A']['hgvs_t_and_p']['NM_000334.4']['p_hgvs_slc'] == 'NP_000325.4:p.?'
        assert results['17-62023006-C-A']['hgvs_t_and_p']['NM_000334.4']['transcript_variant_error'] is None

    def test_variant134(self):
        variant = '17-62034787-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '17-62034787-G-A' in results.keys()
        assert results['17-62034787-G-A']['p_vcf'] == '17-62034787-G-A'
        assert results['17-62034787-G-A']['g_hgvs'] == 'NC_000017.10:g.62034787G>A'
        assert results['17-62034787-G-A']['genomic_variant_error'] is None
        assert 'NM_000334.4' in results['17-62034787-G-A']['hgvs_t_and_p'].keys()
        assert results['17-62034787-G-A']['hgvs_t_and_p']['NM_000334.4']['t_hgvs'] == 'NM_000334.4:c.2111C>T'
        assert results['17-62034787-G-A']['hgvs_t_and_p']['NM_000334.4']['p_hgvs_tlc'] == 'NP_000325.4:p.(Thr704Met)'
        assert results['17-62034787-G-A']['hgvs_t_and_p']['NM_000334.4']['p_hgvs_slc'] == 'NP_000325.4:p.(T704M)'
        assert results['17-62034787-G-A']['hgvs_t_and_p']['NM_000334.4']['transcript_variant_error'] is None

    def test_variant135(self):
        variant = '18-24128261-GTCCTCC-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '18-24128261-GTCCTCC-G' in results.keys()
        assert results['18-24128261-GTCCTCC-G']['p_vcf'] == '18-24128261-GTCCTCC-G'
        assert results['18-24128261-GTCCTCC-G']['g_hgvs'] == 'NC_000018.9:g.24128273_24128278del'
        assert results['18-24128261-GTCCTCC-G']['genomic_variant_error'] is None
        assert 'NM_001258221.1' in results['18-24128261-GTCCTCC-G']['hgvs_t_and_p'].keys()
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001258221.1']['t_hgvs'] == 'NM_001258221.1:c.-16+1426_-16+1431del'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001258221.1']['p_hgvs_tlc'] == 'NP_001245150.1:p.?'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001258221.1']['p_hgvs_slc'] == 'NP_001245150.1:p.?'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001258221.1']['transcript_variant_error'] is None
        assert 'NM_001142730.2' in results['18-24128261-GTCCTCC-G']['hgvs_t_and_p'].keys()
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001142730.2']['t_hgvs'] == 'NM_001142730.2:c.234_239del'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001142730.2']['p_hgvs_tlc'] == 'NP_001136202.1:p.(Glu78_Glu79del)'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001142730.2']['p_hgvs_slc'] == 'NP_001136202.1:p.(E78_E79del)'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001142730.2']['transcript_variant_error'] is None
        assert 'NM_001351443.1' in results['18-24128261-GTCCTCC-G']['hgvs_t_and_p'].keys()
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001351443.1']['t_hgvs'] == 'NM_001351443.1:c.-16+941_-16+946del'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001351443.1']['p_hgvs_tlc'] == 'NP_001338372.1:p.?'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001351443.1']['p_hgvs_slc'] == 'NP_001338372.1:p.?'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001351443.1']['transcript_variant_error'] is None
        assert 'NM_001136205.2' in results['18-24128261-GTCCTCC-G']['hgvs_t_and_p'].keys()
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001136205.2']['t_hgvs'] == 'NM_001136205.2:c.-16+588_-16+593del'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001136205.2']['p_hgvs_tlc'] == 'NP_001129677.1:p.?'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001136205.2']['p_hgvs_slc'] == 'NP_001129677.1:p.?'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001136205.2']['transcript_variant_error'] is None
        assert 'NM_198991.3' in results['18-24128261-GTCCTCC-G']['hgvs_t_and_p'].keys()
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_198991.3']['t_hgvs'] == 'NM_198991.3:c.-15-47053_-15-47048del'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_198991.3']['p_hgvs_tlc'] == 'NP_945342.1:p.?'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_198991.3']['p_hgvs_slc'] == 'NP_945342.1:p.?'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_198991.3']['transcript_variant_error'] is None
        assert 'NM_001258222.2' in results['18-24128261-GTCCTCC-G']['hgvs_t_and_p'].keys()
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001258222.2']['t_hgvs'] == 'NM_001258222.2:c.10-47053_10-47048del'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001258222.2']['p_hgvs_tlc'] == 'NP_001245151.1:p.?'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001258222.2']['p_hgvs_slc'] == 'NP_001245151.1:p.?'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001258222.2']['transcript_variant_error'] is None
        assert 'NM_001258222.1' in results['18-24128261-GTCCTCC-G']['hgvs_t_and_p'].keys()
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001258222.1']['t_hgvs'] == 'NM_001258222.1:c.10-47053_10-47048del'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001258222.1']['p_hgvs_tlc'] == 'NP_001245151.1:p.?'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001258222.1']['p_hgvs_slc'] == 'NP_001245151.1:p.?'
        assert results['18-24128261-GTCCTCC-G']['hgvs_t_and_p']['NM_001258222.1']['transcript_variant_error'] is None

    def test_variant136(self):
        variant = '19-15291774-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '19-15291774-G-A' in results.keys()
        assert results['19-15291774-G-A']['p_vcf'] == '19-15291774-G-A'
        assert results['19-15291774-G-A']['g_hgvs'] == 'NC_000019.9:g.15291774G>A'
        assert results['19-15291774-G-A']['genomic_variant_error'] is None
        assert 'NM_000435.2' in results['19-15291774-G-A']['hgvs_t_and_p'].keys()
        assert results['19-15291774-G-A']['hgvs_t_and_p']['NM_000435.2']['t_hgvs'] == 'NM_000435.2:c.2992C>T'
        assert results['19-15291774-G-A']['hgvs_t_and_p']['NM_000435.2']['p_hgvs_tlc'] == 'NP_000426.2:p.(Gln998Ter)'
        assert results['19-15291774-G-A']['hgvs_t_and_p']['NM_000435.2']['p_hgvs_slc'] == 'NP_000426.2:p.(Q998*)'
        assert results['19-15291774-G-A']['hgvs_t_and_p']['NM_000435.2']['transcript_variant_error'] is None

    # Test removed. No longer intergenic in VVTA
    # def test_variant137(self):
    # 	variant = '19-15311794-A-G'
    # 	results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
    # 	results = results.stucture_data()
    # 	print(results)
    #
    # 	assert '19-15311794-A-G' in results.keys()
    # 	assert results['19-15311794-A-G']['p_vcf'] == '19-15311794-A-G'
    # 	assert results['19-15311794-A-G']['g_hgvs'] == 'NC_000019.9:g.15311794A>G'
    # 	assert results['19-15311794-A-G']['genomic_variant_error'] is None
    # 	assert results['19-15311794-A-G']['hgvs_t_and_p'] == {'intergenic': {'alt_genomic_loci': None}}

    def test_variant138(self):
        variant = '19-39076592-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '19-39076592-G-A' in results.keys()
        assert results['19-39076592-G-A']['p_vcf'] == '19-39076592-G-A'
        assert results['19-39076592-G-A']['g_hgvs'] == 'NC_000019.9:g.39076592G>A'
        assert results['19-39076592-G-A']['genomic_variant_error'] is None
        assert 'NM_001042723.1' in results['19-39076592-G-A']['hgvs_t_and_p'].keys()
        assert results['19-39076592-G-A']['hgvs_t_and_p']['NM_001042723.1']['t_hgvs'] == 'NM_001042723.1:c.14803G>A'
        assert results['19-39076592-G-A']['hgvs_t_and_p']['NM_001042723.1']['p_hgvs_tlc'] == 'NP_001036188.1:p.(Ala4935Thr)'
        assert results['19-39076592-G-A']['hgvs_t_and_p']['NM_001042723.1']['p_hgvs_slc'] == 'NP_001036188.1:p.(A4935T)'
        assert results['19-39076592-G-A']['hgvs_t_and_p']['NM_001042723.1']['transcript_variant_error'] is None
        assert 'NM_000540.2' in results['19-39076592-G-A']['hgvs_t_and_p'].keys()
        assert results['19-39076592-G-A']['hgvs_t_and_p']['NM_000540.2']['t_hgvs'] == 'NM_000540.2:c.14818G>A'
        assert results['19-39076592-G-A']['hgvs_t_and_p']['NM_000540.2']['p_hgvs_tlc'] == 'NP_000531.2:p.(Ala4940Thr)'
        assert results['19-39076592-G-A']['hgvs_t_and_p']['NM_000540.2']['p_hgvs_slc'] == 'NP_000531.2:p.(A4940T)'
        assert results['19-39076592-G-A']['hgvs_t_and_p']['NM_000540.2']['transcript_variant_error'] is None

    def test_variant139(self):
        variant = '2-50149352-T-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '2-50149352-T-C' in results.keys()
        assert results['2-50149352-T-C']['p_vcf'] == '2-50149352-T-C'
        assert results['2-50149352-T-C']['g_hgvs'] == 'NC_000002.11:g.50149352T>C'
        assert results['2-50149352-T-C']['genomic_variant_error'] is None
        assert 'NM_001330083.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330083.1']['t_hgvs'] == 'NM_001330083.1:c.4089A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330083.1']['p_hgvs_tlc'] == 'NP_001317012.1:p.(Pro1363=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330083.1']['p_hgvs_slc'] == 'NP_001317012.1:p.(P1363=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330083.1']['transcript_variant_error'] is None
        assert 'NM_001330082.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330082.1']['t_hgvs'] == 'NM_001330082.1:c.4221A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330082.1']['p_hgvs_tlc'] == 'NP_001317011.1:p.(Pro1407=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330082.1']['p_hgvs_slc'] == 'NP_001317011.1:p.(P1407=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330082.1']['transcript_variant_error'] is None
        assert 'NM_001330095.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330095.1']['t_hgvs'] == 'NM_001330095.1:c.4113A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330095.1']['p_hgvs_tlc'] == 'NP_001317024.1:p.(Pro1371=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330095.1']['p_hgvs_slc'] == 'NP_001317024.1:p.(P1371=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330095.1']['transcript_variant_error'] is None
        assert 'NM_001330097.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330097.1']['t_hgvs'] == 'NM_001330097.1:c.1050A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330097.1']['p_hgvs_tlc'] == 'NP_001317026.1:p.(Pro350=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330097.1']['p_hgvs_slc'] == 'NP_001317026.1:p.(P350=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330097.1']['transcript_variant_error'] is None
        assert 'NM_001330085.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330085.1']['t_hgvs'] == 'NM_001330085.1:c.4227A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330085.1']['p_hgvs_tlc'] == 'NP_001317014.1:p.(Pro1409=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330085.1']['p_hgvs_slc'] == 'NP_001317014.1:p.(P1409=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330085.1']['transcript_variant_error'] is None
        assert 'NM_001330084.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330084.1']['t_hgvs'] == 'NM_001330084.1:c.4188A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330084.1']['p_hgvs_tlc'] == 'NP_001317013.1:p.(Pro1396=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330084.1']['p_hgvs_slc'] == 'NP_001317013.1:p.(P1396=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330084.1']['transcript_variant_error'] is None
        assert 'NM_001330094.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330094.1']['t_hgvs'] == 'NM_001330094.1:c.4233A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330094.1']['p_hgvs_tlc'] == 'NP_001317023.1:p.(Pro1411=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330094.1']['p_hgvs_slc'] == 'NP_001317023.1:p.(P1411=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330094.1']['transcript_variant_error'] is None
        assert 'NM_138735.2' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_138735.2']['t_hgvs'] == 'NM_138735.2:c.1059A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_138735.2']['p_hgvs_tlc'] == 'NP_620072.1:p.(Pro353=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_138735.2']['p_hgvs_slc'] == 'NP_620072.1:p.(P353=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_138735.2']['transcript_variant_error'] is None
        assert 'NM_138735.4' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_138735.4']['t_hgvs'] == 'NM_138735.4:c.1059A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_138735.4']['p_hgvs_tlc'] == 'NP_620072.1:p.(Pro353=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_138735.4']['p_hgvs_slc'] == 'NP_620072.1:p.(P353=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_138735.4']['transcript_variant_error'] is None
        assert 'NM_004801.4' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_004801.4']['t_hgvs'] == 'NM_004801.4:c.4164A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_004801.4']['p_hgvs_tlc'] == 'NP_004792.1:p.(Pro1388=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_004801.4']['p_hgvs_slc'] == 'NP_004792.1:p.(P1388=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_004801.4']['transcript_variant_error'] is None
        assert 'NM_004801.5' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_004801.5']['t_hgvs'] == 'NM_004801.5:c.4164A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_004801.5']['p_hgvs_tlc'] == 'NP_004792.1:p.(Pro1388=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_004801.5']['p_hgvs_slc'] == 'NP_004792.1:p.(P1388=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_004801.5']['transcript_variant_error'] is None
        assert 'NM_001330077.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330077.1']['t_hgvs'] == 'NM_001330077.1:c.4230A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330077.1']['p_hgvs_tlc'] == 'NP_001317006.1:p.(Pro1410=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330077.1']['p_hgvs_slc'] == 'NP_001317006.1:p.(P1410=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330077.1']['transcript_variant_error'] is None
        assert 'NM_001135659.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001135659.1']['t_hgvs'] == 'NM_001135659.1:c.4374A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001135659.1']['p_hgvs_tlc'] == 'NP_001129131.1:p.(Pro1458=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001135659.1']['p_hgvs_slc'] == 'NP_001129131.1:p.(P1458=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001135659.1']['transcript_variant_error'] is None
        assert 'NM_001135659.2' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001135659.2']['t_hgvs'] == 'NM_001135659.2:c.4374A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001135659.2']['p_hgvs_tlc'] == 'NP_001129131.1:p.(Pro1458=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001135659.2']['p_hgvs_slc'] == 'NP_001129131.1:p.(P1458=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001135659.2']['transcript_variant_error'] is None
        assert 'NM_001330087.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330087.1']['t_hgvs'] == 'NM_001330087.1:c.4053A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330087.1']['p_hgvs_tlc'] == 'NP_001317016.1:p.(Pro1351=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330087.1']['p_hgvs_slc'] == 'NP_001317016.1:p.(P1351=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330087.1']['transcript_variant_error'] is None
        assert 'NM_001330086.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330086.1']['t_hgvs'] == 'NM_001330086.1:c.4245A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330086.1']['p_hgvs_tlc'] == 'NP_001317015.1:p.(Pro1415=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330086.1']['p_hgvs_slc'] == 'NP_001317015.1:p.(P1415=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330086.1']['transcript_variant_error'] is None
        assert 'NM_001330096.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330096.1']['t_hgvs'] == 'NM_001330096.1:c.4044A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330096.1']['p_hgvs_tlc'] == 'NP_001317025.1:p.(Pro1348=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330096.1']['p_hgvs_slc'] == 'NP_001317025.1:p.(P1348=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330096.1']['transcript_variant_error'] is None
        assert 'NM_001330092.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330092.1']['t_hgvs'] == 'NM_001330092.1:c.1149A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330092.1']['p_hgvs_tlc'] == 'NP_001317021.1:p.(Pro383=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330092.1']['p_hgvs_slc'] == 'NP_001317021.1:p.(P383=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330092.1']['transcript_variant_error'] is None
        assert 'NM_001330093.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330093.1']['t_hgvs'] == 'NM_001330093.1:c.4251A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330093.1']['p_hgvs_tlc'] == 'NP_001317022.1:p.(Pro1417=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330093.1']['p_hgvs_slc'] == 'NP_001317022.1:p.(P1417=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330093.1']['transcript_variant_error'] is None
        assert 'NM_001330088.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330088.1']['t_hgvs'] == 'NM_001330088.1:c.4074A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330088.1']['p_hgvs_tlc'] == 'NP_001317017.1:p.(Pro1358=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330088.1']['p_hgvs_slc'] == 'NP_001317017.1:p.(P1358=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330088.1']['transcript_variant_error'] is None
        assert 'NM_001320156.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320156.1']['t_hgvs'] == 'NM_001320156.1:c.159A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320156.1']['p_hgvs_tlc'] == 'NP_001307085.1:p.(Pro53=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320156.1']['p_hgvs_slc'] == 'NP_001307085.1:p.(P53=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320156.1']['transcript_variant_error'] is None
        assert 'NM_001320156.3' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320156.3']['t_hgvs'] == 'NM_001320156.3:c.159A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320156.3']['p_hgvs_tlc'] == 'NP_001307085.1:p.(Pro53=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320156.3']['p_hgvs_slc'] == 'NP_001307085.1:p.(P53=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320156.3']['transcript_variant_error'] is None
        assert 'NM_001320157.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320157.1']['t_hgvs'] == 'NM_001320157.1:c.150A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320157.1']['p_hgvs_tlc'] == 'NP_001307086.1:p.(Pro50=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320157.1']['p_hgvs_slc'] == 'NP_001307086.1:p.(P50=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320157.1']['transcript_variant_error'] is None
        assert 'NM_001320157.3' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320157.3']['t_hgvs'] == 'NM_001320157.3:c.150A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320157.3']['p_hgvs_tlc'] == 'NP_001307086.1:p.(Pro50=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320157.3']['p_hgvs_slc'] == 'NP_001307086.1:p.(P50=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001320157.3']['transcript_variant_error'] is None
        assert 'NM_001330091.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330091.1']['t_hgvs'] == 'NM_001330091.1:c.1140A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330091.1']['p_hgvs_tlc'] == 'NP_001317020.1:p.(Pro380=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330091.1']['p_hgvs_slc'] == 'NP_001317020.1:p.(P380=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330091.1']['transcript_variant_error'] is None
        assert 'NM_001330078.1' in results['2-50149352-T-C']['hgvs_t_and_p'].keys()
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330078.1']['t_hgvs'] == 'NM_001330078.1:c.4254A>G'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330078.1']['p_hgvs_tlc'] == 'NP_001317007.1:p.(Pro1418=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330078.1']['p_hgvs_slc'] == 'NP_001317007.1:p.(P1418=)'
        assert results['2-50149352-T-C']['hgvs_t_and_p']['NM_001330078.1']['transcript_variant_error'] is None

    def test_variant140(self):
        variant = '2-50847195-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '2-50847195-G-A' in results.keys()
        assert results['2-50847195-G-A']['p_vcf'] == '2-50847195-G-A'
        assert results['2-50847195-G-A']['g_hgvs'] == 'NC_000002.11:g.50847195G>A'
        assert results['2-50847195-G-A']['genomic_variant_error'] is None
        assert 'NM_001330088.1' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330088.1']['t_hgvs'] == 'NM_001330088.1:c.1231C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330088.1']['p_hgvs_tlc'] == 'NP_001317017.1:p.(Pro411Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330088.1']['p_hgvs_slc'] == 'NP_001317017.1:p.(P411S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330088.1']['transcript_variant_error'] is None
        assert 'NM_001135659.1' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001135659.1']['t_hgvs'] == 'NM_001135659.1:c.1405C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001135659.1']['p_hgvs_tlc'] == 'NP_001129131.1:p.(Pro469Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001135659.1']['p_hgvs_slc'] == 'NP_001129131.1:p.(P469S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001135659.1']['transcript_variant_error'] is None
        assert 'NM_001135659.2' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001135659.2']['t_hgvs'] == 'NM_001135659.2:c.1405C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001135659.2']['p_hgvs_tlc'] == 'NP_001129131.1:p.(Pro469Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001135659.2']['p_hgvs_slc'] == 'NP_001129131.1:p.(P469S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001135659.2']['transcript_variant_error'] is None
        assert 'NM_001330096.1' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330096.1']['t_hgvs'] == 'NM_001330096.1:c.1201C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330096.1']['p_hgvs_tlc'] == 'NP_001317025.1:p.(Pro401Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330096.1']['p_hgvs_slc'] == 'NP_001317025.1:p.(P401S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330096.1']['transcript_variant_error'] is None
        assert 'NM_001330083.1' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330083.1']['t_hgvs'] == 'NM_001330083.1:c.1246C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330083.1']['p_hgvs_tlc'] == 'NP_001317012.1:p.(Pro416Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330083.1']['p_hgvs_slc'] == 'NP_001317012.1:p.(P416S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330083.1']['transcript_variant_error'] is None
        assert 'NM_001330082.1' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330082.1']['t_hgvs'] == 'NM_001330082.1:c.1261C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330082.1']['p_hgvs_tlc'] == 'NP_001317011.1:p.(Pro421Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330082.1']['p_hgvs_slc'] == 'NP_001317011.1:p.(P421S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330082.1']['transcript_variant_error'] is None
        assert 'NM_001330087.1' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330087.1']['t_hgvs'] == 'NM_001330087.1:c.1201C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330087.1']['p_hgvs_tlc'] == 'NP_001317016.1:p.(Pro401Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330087.1']['p_hgvs_slc'] == 'NP_001317016.1:p.(P401S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330087.1']['transcript_variant_error'] is None
        assert 'NM_001330095.1' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330095.1']['t_hgvs'] == 'NM_001330095.1:c.1261C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330095.1']['p_hgvs_tlc'] == 'NP_001317024.1:p.(Pro421Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330095.1']['p_hgvs_slc'] == 'NP_001317024.1:p.(P421S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330095.1']['transcript_variant_error'] is None
        assert 'NM_004801.4' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_004801.4']['t_hgvs'] == 'NM_004801.4:c.1285C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_004801.4']['p_hgvs_tlc'] == 'NP_004792.1:p.(Pro429Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_004801.4']['p_hgvs_slc'] == 'NP_004792.1:p.(P429S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_004801.4']['transcript_variant_error'] is None
        assert 'NM_004801.5' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_004801.5']['t_hgvs'] == 'NM_004801.5:c.1285C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_004801.5']['p_hgvs_tlc'] == 'NP_004792.1:p.(Pro429Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_004801.5']['p_hgvs_slc'] == 'NP_004792.1:p.(P429S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_004801.5']['transcript_variant_error'] is None
        assert 'NM_001330077.1' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330077.1']['t_hgvs'] == 'NM_001330077.1:c.1261C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330077.1']['p_hgvs_tlc'] == 'NP_001317006.1:p.(Pro421Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330077.1']['p_hgvs_slc'] == 'NP_001317006.1:p.(P421S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330077.1']['transcript_variant_error'] is None
        assert 'NM_001330093.1' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330093.1']['t_hgvs'] == 'NM_001330093.1:c.1282C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330093.1']['p_hgvs_tlc'] == 'NP_001317022.1:p.(Pro428Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330093.1']['p_hgvs_slc'] == 'NP_001317022.1:p.(P428S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330093.1']['transcript_variant_error'] is None
        assert 'NM_001330086.1' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330086.1']['t_hgvs'] == 'NM_001330086.1:c.1285C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330086.1']['p_hgvs_tlc'] == 'NP_001317015.1:p.(Pro429Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330086.1']['p_hgvs_slc'] == 'NP_001317015.1:p.(P429S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330086.1']['transcript_variant_error'] is None
        assert 'NM_001330078.1' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330078.1']['t_hgvs'] == 'NM_001330078.1:c.1285C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330078.1']['p_hgvs_tlc'] == 'NP_001317007.1:p.(Pro429Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330078.1']['p_hgvs_slc'] == 'NP_001317007.1:p.(P429S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330078.1']['transcript_variant_error'] is None
        assert 'NM_001330085.1' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330085.1']['t_hgvs'] == 'NM_001330085.1:c.1285C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330085.1']['p_hgvs_tlc'] == 'NP_001317014.1:p.(Pro429Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330085.1']['p_hgvs_slc'] == 'NP_001317014.1:p.(P429S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330085.1']['transcript_variant_error'] is None
        assert 'NM_001330084.1' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330084.1']['t_hgvs'] == 'NM_001330084.1:c.1246C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330084.1']['p_hgvs_tlc'] == 'NP_001317013.1:p.(Pro416Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330084.1']['p_hgvs_slc'] == 'NP_001317013.1:p.(P416S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330084.1']['transcript_variant_error'] is None
        assert 'NM_001330094.1' in results['2-50847195-G-A']['hgvs_t_and_p'].keys()
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330094.1']['t_hgvs'] == 'NM_001330094.1:c.1273C>T'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330094.1']['p_hgvs_tlc'] == 'NP_001317023.1:p.(Pro425Ser)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330094.1']['p_hgvs_slc'] == 'NP_001317023.1:p.(P425S)'
        assert results['2-50847195-G-A']['hgvs_t_and_p']['NM_001330094.1']['transcript_variant_error'] is None

    def test_variant141(self):
        variant = '2-71825797-C-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '2-71825797-C-G' in results.keys()
        assert results['2-71825797-C-G']['p_vcf'] == '2-71825797-C-G'
        assert results['2-71825797-C-G']['g_hgvs'] == 'NC_000002.11:g.71825797C>G'
        assert results['2-71825797-C-G']['genomic_variant_error'] is None
        assert 'NM_001130455.1' in results['2-71825797-C-G']['hgvs_t_and_p'].keys()
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130455.1']['t_hgvs'] == 'NM_001130455.1:c.3627C>G'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130455.1']['p_hgvs_tlc'] == 'NP_001123927.1:p.(Ile1209Met)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130455.1']['p_hgvs_slc'] == 'NP_001123927.1:p.(I1209M)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130455.1']['transcript_variant_error'] is None
        assert 'NM_001130982.1' in results['2-71825797-C-G']['hgvs_t_and_p'].keys()
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130982.1']['t_hgvs'] == 'NM_001130982.1:c.3720C>G'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130982.1']['p_hgvs_tlc'] == 'NP_001124454.1:p.(Ile1240Met)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130982.1']['p_hgvs_slc'] == 'NP_001124454.1:p.(I1240M)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130982.1']['transcript_variant_error'] is None
        assert 'NM_001130983.1' in results['2-71825797-C-G']['hgvs_t_and_p'].keys()
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130983.1']['t_hgvs'] == 'NM_001130983.1:c.3627C>G'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130983.1']['p_hgvs_tlc'] == 'NP_001124455.1:p.(Ile1209Met)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130983.1']['p_hgvs_slc'] == 'NP_001124455.1:p.(I1209M)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130983.1']['transcript_variant_error'] is None
        assert 'NM_001130977.1' in results['2-71825797-C-G']['hgvs_t_and_p'].keys()
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130977.1']['t_hgvs'] == 'NM_001130977.1:c.3582C>G'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130977.1']['p_hgvs_tlc'] == 'NP_001124449.1:p.(Ile1194Met)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130977.1']['p_hgvs_slc'] == 'NP_001124449.1:p.(I1194M)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130977.1']['transcript_variant_error'] is None
        assert 'NM_001130984.1' in results['2-71825797-C-G']['hgvs_t_and_p'].keys()
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130984.1']['t_hgvs'] == 'NM_001130984.1:c.3585C>G'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130984.1']['p_hgvs_tlc'] == 'NP_001124456.1:p.(Ile1195Met)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130984.1']['p_hgvs_slc'] == 'NP_001124456.1:p.(I1195M)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130984.1']['transcript_variant_error'] is None
        assert 'NM_001130976.1' in results['2-71825797-C-G']['hgvs_t_and_p'].keys()
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130976.1']['t_hgvs'] == 'NM_001130976.1:c.3582C>G'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130976.1']['p_hgvs_tlc'] == 'NP_001124448.1:p.(Ile1194Met)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130976.1']['p_hgvs_slc'] == 'NP_001124448.1:p.(I1194M)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130976.1']['transcript_variant_error'] is None
        assert 'NM_001130985.1' in results['2-71825797-C-G']['hgvs_t_and_p'].keys()
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130985.1']['t_hgvs'] == 'NM_001130985.1:c.3678C>G'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130985.1']['p_hgvs_tlc'] == 'NP_001124457.1:p.(Ile1226Met)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130985.1']['p_hgvs_slc'] == 'NP_001124457.1:p.(I1226M)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130985.1']['transcript_variant_error'] is None
        assert 'NM_001130980.1' in results['2-71825797-C-G']['hgvs_t_and_p'].keys()
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130980.1']['t_hgvs'] == 'NM_001130980.1:c.3675C>G'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130980.1']['p_hgvs_tlc'] == 'NP_001124452.1:p.(Ile1225Met)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130980.1']['p_hgvs_slc'] == 'NP_001124452.1:p.(I1225M)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130980.1']['transcript_variant_error'] is None
        assert 'NM_001130987.1' in results['2-71825797-C-G']['hgvs_t_and_p'].keys()
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130987.1']['t_hgvs'] == 'NM_001130987.1:c.3678C>G'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130987.1']['p_hgvs_tlc'] == 'NP_001124459.1:p.(Ile1226Met)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130987.1']['p_hgvs_slc'] == 'NP_001124459.1:p.(I1226M)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130987.1']['transcript_variant_error'] is None
        assert 'NM_001130978.1' in results['2-71825797-C-G']['hgvs_t_and_p'].keys()
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130978.1']['t_hgvs'] == 'NM_001130978.1:c.3624C>G'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130978.1']['p_hgvs_tlc'] == 'NP_001124450.1:p.(Ile1208Met)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130978.1']['p_hgvs_slc'] == 'NP_001124450.1:p.(I1208M)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130978.1']['transcript_variant_error'] is None
        assert 'NM_001130981.1' in results['2-71825797-C-G']['hgvs_t_and_p'].keys()
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130981.1']['t_hgvs'] == 'NM_001130981.1:c.3675C>G'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130981.1']['p_hgvs_tlc'] == 'NP_001124453.1:p.(Ile1225Met)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130981.1']['p_hgvs_slc'] == 'NP_001124453.1:p.(I1225M)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130981.1']['transcript_variant_error'] is None
        assert 'NM_001130979.1' in results['2-71825797-C-G']['hgvs_t_and_p'].keys()
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130979.1']['t_hgvs'] == 'NM_001130979.1:c.3717C>G'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130979.1']['p_hgvs_tlc'] == 'NP_001124451.1:p.(Ile1239Met)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130979.1']['p_hgvs_slc'] == 'NP_001124451.1:p.(I1239M)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130979.1']['transcript_variant_error'] is None
        assert 'NM_001130986.1' in results['2-71825797-C-G']['hgvs_t_and_p'].keys()
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130986.1']['t_hgvs'] == 'NM_001130986.1:c.3585C>G'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130986.1']['p_hgvs_tlc'] == 'NP_001124458.1:p.(Ile1195Met)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130986.1']['p_hgvs_slc'] == 'NP_001124458.1:p.(I1195M)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_001130986.1']['transcript_variant_error'] is None
        assert 'NM_003494.3' in results['2-71825797-C-G']['hgvs_t_and_p'].keys()
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_003494.3']['t_hgvs'] == 'NM_003494.3:c.3624C>G'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_003494.3']['p_hgvs_tlc'] == 'NP_003485.1:p.(Ile1208Met)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_003494.3']['p_hgvs_slc'] == 'NP_003485.1:p.(I1208M)'
        assert results['2-71825797-C-G']['hgvs_t_and_p']['NM_003494.3']['transcript_variant_error'] is None

    def test_variant142(self):
        variant = '2-166179712-G-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '2-166179712-G-C' in results.keys()
        assert results['2-166179712-G-C']['p_vcf'] == '2-166179712-G-C'
        assert results['2-166179712-G-C']['g_hgvs'] == 'NC_000002.11:g.166179712G>C'
        assert results['2-166179712-G-C']['genomic_variant_error'] is None
        assert 'NM_001040142.1' in results['2-166179712-G-C']['hgvs_t_and_p'].keys()
        assert results['2-166179712-G-C']['hgvs_t_and_p']['NM_001040142.1']['t_hgvs'] == 'NM_001040142.1:c.1718G>C'
        assert results['2-166179712-G-C']['hgvs_t_and_p']['NM_001040142.1']['p_hgvs_tlc'] == 'NP_001035232.1:p.(Ser573Thr)'
        assert results['2-166179712-G-C']['hgvs_t_and_p']['NM_001040142.1']['p_hgvs_slc'] == 'NP_001035232.1:p.(S573T)'
        assert results['2-166179712-G-C']['hgvs_t_and_p']['NM_001040142.1']['transcript_variant_error'] is None
        assert 'NM_001040143.1' in results['2-166179712-G-C']['hgvs_t_and_p'].keys()
        assert results['2-166179712-G-C']['hgvs_t_and_p']['NM_001040143.1']['t_hgvs'] == 'NM_001040143.1:c.1718G>C'
        assert results['2-166179712-G-C']['hgvs_t_and_p']['NM_001040143.1']['p_hgvs_tlc'] == 'NP_001035233.1:p.(Ser573Thr)'
        assert results['2-166179712-G-C']['hgvs_t_and_p']['NM_001040143.1']['p_hgvs_slc'] == 'NP_001035233.1:p.(S573T)'
        assert results['2-166179712-G-C']['hgvs_t_and_p']['NM_001040143.1']['transcript_variant_error'] is None
        assert 'NM_021007.2' in results['2-166179712-G-C']['hgvs_t_and_p'].keys()
        assert results['2-166179712-G-C']['hgvs_t_and_p']['NM_021007.2']['t_hgvs'] == 'NM_021007.2:c.1718G>C'
        assert results['2-166179712-G-C']['hgvs_t_and_p']['NM_021007.2']['p_hgvs_tlc'] == 'NP_066287.2:p.(Ser573Thr)'
        assert results['2-166179712-G-C']['hgvs_t_and_p']['NM_021007.2']['p_hgvs_slc'] == 'NP_066287.2:p.(S573T)'
        assert results['2-166179712-G-C']['hgvs_t_and_p']['NM_021007.2']['transcript_variant_error'] is None

    def test_variant143(self):
        variant = '2-166183371-A-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '2-166183371-A-G' in results.keys()
        assert results['2-166183371-A-G']['p_vcf'] == '2-166183371-A-G'
        assert results['2-166183371-A-G']['g_hgvs'] == 'NC_000002.11:g.166183371A>G'
        assert results['2-166183371-A-G']['genomic_variant_error'] is None
        assert 'NM_001040142.1' in results['2-166183371-A-G']['hgvs_t_and_p'].keys()
        assert results['2-166183371-A-G']['hgvs_t_and_p']['NM_001040142.1']['t_hgvs'] == 'NM_001040142.1:c.2026A>G'
        assert results['2-166183371-A-G']['hgvs_t_and_p']['NM_001040142.1']['p_hgvs_tlc'] == 'NP_001035232.1:p.(Thr676Ala)'
        assert results['2-166183371-A-G']['hgvs_t_and_p']['NM_001040142.1']['p_hgvs_slc'] == 'NP_001035232.1:p.(T676A)'
        assert results['2-166183371-A-G']['hgvs_t_and_p']['NM_001040142.1']['transcript_variant_error'] is None
        assert 'NM_001040143.1' in results['2-166183371-A-G']['hgvs_t_and_p'].keys()
        assert results['2-166183371-A-G']['hgvs_t_and_p']['NM_001040143.1']['t_hgvs'] == 'NM_001040143.1:c.2026A>G'
        assert results['2-166183371-A-G']['hgvs_t_and_p']['NM_001040143.1']['p_hgvs_tlc'] == 'NP_001035233.1:p.(Thr676Ala)'
        assert results['2-166183371-A-G']['hgvs_t_and_p']['NM_001040143.1']['p_hgvs_slc'] == 'NP_001035233.1:p.(T676A)'
        assert results['2-166183371-A-G']['hgvs_t_and_p']['NM_001040143.1']['transcript_variant_error'] is None
        assert 'NM_021007.2' in results['2-166183371-A-G']['hgvs_t_and_p'].keys()
        assert results['2-166183371-A-G']['hgvs_t_and_p']['NM_021007.2']['t_hgvs'] == 'NM_021007.2:c.2026A>G'
        assert results['2-166183371-A-G']['hgvs_t_and_p']['NM_021007.2']['p_hgvs_tlc'] == 'NP_066287.2:p.(Thr676Ala)'
        assert results['2-166183371-A-G']['hgvs_t_and_p']['NM_021007.2']['p_hgvs_slc'] == 'NP_066287.2:p.(T676A)'
        assert results['2-166183371-A-G']['hgvs_t_and_p']['NM_021007.2']['transcript_variant_error'] is None

    def test_variant144(self):
        variant = '2-166929889-GTCCAGGTCCT-GAC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '2-166929889-GTCCAGGTCCT-GAC' in results.keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['p_vcf'] == '2-166929890-TCCAGGTCCT-AC'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['g_hgvs'] == 'NC_000002.11:g.166929890_166929899delinsAC'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['genomic_variant_error'] is None
        assert 'NM_001353952.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353952.1']['t_hgvs'] == 'NM_001353952.1:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353952.1']['p_hgvs_tlc'] == 'NP_001340881.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353952.1']['p_hgvs_slc'] == 'NP_001340881.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353952.1']['transcript_variant_error'] is None
        assert 'NM_001353954.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353954.1']['t_hgvs'] == 'NM_001353954.1:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353954.1']['p_hgvs_tlc'] == 'NP_001340883.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353954.1']['p_hgvs_slc'] == 'NP_001340883.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353954.1']['transcript_variant_error'] is None
        assert 'NM_001353948.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353948.1']['t_hgvs'] == 'NM_001353948.1:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353948.1']['p_hgvs_tlc'] == 'NP_001340877.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353948.1']['p_hgvs_slc'] == 'NP_001340877.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353948.1']['transcript_variant_error'] is None
        assert 'NM_001353958.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353958.1']['t_hgvs'] == 'NM_001353958.1:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353958.1']['p_hgvs_tlc'] == 'NP_001340887.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353958.1']['p_hgvs_slc'] == 'NP_001340887.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353958.1']['transcript_variant_error'] is None
        assert 'NM_001353955.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353955.1']['t_hgvs'] == 'NM_001353955.1:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353955.1']['p_hgvs_tlc'] == 'NP_001340884.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353955.1']['p_hgvs_slc'] == 'NP_001340884.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353955.1']['transcript_variant_error'] is None
        assert 'NM_001353950.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353950.1']['t_hgvs'] == 'NM_001353950.1:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353950.1']['p_hgvs_tlc'] == 'NP_001340879.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353950.1']['p_hgvs_slc'] == 'NP_001340879.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353950.1']['transcript_variant_error'] is None
        assert 'NM_001165964.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165964.1']['t_hgvs'] == 'NM_001165964.1:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165964.1']['p_hgvs_tlc'] == 'NP_001159436.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165964.1']['p_hgvs_slc'] == 'NP_001159436.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165964.1']['transcript_variant_error'] is None
        assert 'NM_001353951.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353951.1']['t_hgvs'] == 'NM_001353951.1:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353951.1']['p_hgvs_tlc'] == 'NP_001340880.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353951.1']['p_hgvs_slc'] == 'NP_001340880.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353951.1']['transcript_variant_error'] is None
        assert 'NM_006920.4' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_006920.4']['t_hgvs'] == 'NM_006920.4:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_006920.4']['p_hgvs_tlc'] == 'NP_008851.3:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_006920.4']['p_hgvs_slc'] == 'NP_008851.3:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_006920.4']['transcript_variant_error'] is None
        assert 'NM_001353961.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353961.1']['t_hgvs'] == 'NM_001353961.1:c.-2193_-2184delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353961.1']['p_hgvs_tlc'] == 'NP_001340890.1:p.?'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353961.1']['p_hgvs_slc'] == 'NP_001340890.1:p.?'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353961.1']['transcript_variant_error'] is None
        assert 'NM_001353949.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353949.1']['t_hgvs'] == 'NM_001353949.1:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353949.1']['p_hgvs_tlc'] == 'NP_001340878.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353949.1']['p_hgvs_slc'] == 'NP_001340878.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353949.1']['transcript_variant_error'] is None
        assert 'NM_001353960.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353960.1']['t_hgvs'] == 'NM_001353960.1:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353960.1']['p_hgvs_tlc'] == 'NP_001340889.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353960.1']['p_hgvs_slc'] == 'NP_001340889.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353960.1']['transcript_variant_error'] is None
        assert 'NM_006920.5' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_006920.5']['t_hgvs'] == 'NM_006920.5:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_006920.5']['p_hgvs_tlc'] == 'NP_008851.3:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_006920.5']['p_hgvs_slc'] == 'NP_008851.3:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_006920.5']['transcript_variant_error'] is None
        assert 'NR_148667.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NR_148667.1']['t_hgvs'] == 'NR_148667.1:n.638_647delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NR_148667.1']['p_hgvs_tlc'] is None
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NR_148667.1']['p_hgvs_slc'] is None
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NR_148667.1']['transcript_variant_error'] is None
        assert 'NM_001165963.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165963.1']['t_hgvs'] == 'NM_001165963.1:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165963.1']['p_hgvs_tlc'] == 'NP_001159435.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165963.1']['p_hgvs_slc'] == 'NP_001159435.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165963.1']['transcript_variant_error'] is None
        assert 'NM_001202435.2' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001202435.2']['t_hgvs'] == 'NM_001202435.2:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001202435.2']['p_hgvs_tlc'] == 'NP_001189364.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001202435.2']['p_hgvs_slc'] == 'NP_001189364.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001202435.2']['transcript_variant_error'] is None
        assert 'NM_001353957.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353957.1']['t_hgvs'] == 'NM_001353957.1:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353957.1']['p_hgvs_tlc'] == 'NP_001340886.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353957.1']['p_hgvs_slc'] == 'NP_001340886.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001353957.1']['transcript_variant_error'] is None
        assert 'NM_001165963.2' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165963.2']['t_hgvs'] == 'NM_001165963.2:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165963.2']['p_hgvs_tlc'] == 'NP_001159435.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165963.2']['p_hgvs_slc'] == 'NP_001159435.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165963.2']['transcript_variant_error'] is None
        assert 'NM_001202435.1' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001202435.1']['t_hgvs'] == 'NM_001202435.1:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001202435.1']['p_hgvs_tlc'] == 'NP_001189364.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001202435.1']['p_hgvs_slc'] == 'NP_001189364.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001202435.1']['transcript_variant_error'] is None
        assert 'NM_001165964.2' in results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p'].keys()
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165964.2']['t_hgvs'] == 'NM_001165964.2:c.233_242delinsGT'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165964.2']['p_hgvs_tlc'] == 'NP_001159436.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165964.2']['p_hgvs_slc'] == 'NP_001159436.1:p.(E78Gfs*7)'
        assert results['2-166929889-GTCCAGGTCCT-GAC']['hgvs_t_and_p']['NM_001165964.2']['transcript_variant_error'] is None

    def test_variant145(self):
        variant = '2-166929891-CCAGGTCCT-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '2-166929891-CCAGGTCCT-C' in results.keys()
        assert results['2-166929891-CCAGGTCCT-C']['p_vcf'] == '2-166929891-CCAGGTCCT-C'
        assert results['2-166929891-CCAGGTCCT-C']['g_hgvs'] == 'NC_000002.11:g.166929893_166929900del'
        assert results['2-166929891-CCAGGTCCT-C']['genomic_variant_error'] is None
        assert 'NM_001353952.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353952.1']['t_hgvs'] == 'NM_001353952.1:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353952.1']['p_hgvs_tlc'] == 'NP_001340881.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353952.1']['p_hgvs_slc'] == 'NP_001340881.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353952.1']['transcript_variant_error'] is None
        assert 'NM_001353954.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353954.1']['t_hgvs'] == 'NM_001353954.1:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353954.1']['p_hgvs_tlc'] == 'NP_001340883.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353954.1']['p_hgvs_slc'] == 'NP_001340883.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353954.1']['transcript_variant_error'] is None
        assert 'NM_001353948.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353948.1']['t_hgvs'] == 'NM_001353948.1:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353948.1']['p_hgvs_tlc'] == 'NP_001340877.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353948.1']['p_hgvs_slc'] == 'NP_001340877.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353948.1']['transcript_variant_error'] is None
        assert 'NM_001353958.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353958.1']['t_hgvs'] == 'NM_001353958.1:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353958.1']['p_hgvs_tlc'] == 'NP_001340887.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353958.1']['p_hgvs_slc'] == 'NP_001340887.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353958.1']['transcript_variant_error'] is None
        assert 'NM_001353955.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353955.1']['t_hgvs'] == 'NM_001353955.1:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353955.1']['p_hgvs_tlc'] == 'NP_001340884.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353955.1']['p_hgvs_slc'] == 'NP_001340884.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353955.1']['transcript_variant_error'] is None
        assert 'NM_001353950.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353950.1']['t_hgvs'] == 'NM_001353950.1:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353950.1']['p_hgvs_tlc'] == 'NP_001340879.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353950.1']['p_hgvs_slc'] == 'NP_001340879.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353950.1']['transcript_variant_error'] is None
        assert 'NM_001165964.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165964.1']['t_hgvs'] == 'NM_001165964.1:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165964.1']['p_hgvs_tlc'] == 'NP_001159436.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165964.1']['p_hgvs_slc'] == 'NP_001159436.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165964.1']['transcript_variant_error'] is None
        assert 'NM_001353951.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353951.1']['t_hgvs'] == 'NM_001353951.1:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353951.1']['p_hgvs_tlc'] == 'NP_001340880.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353951.1']['p_hgvs_slc'] == 'NP_001340880.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353951.1']['transcript_variant_error'] is None
        assert 'NM_006920.4' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_006920.4']['t_hgvs'] == 'NM_006920.4:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_006920.4']['p_hgvs_tlc'] == 'NP_008851.3:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_006920.4']['p_hgvs_slc'] == 'NP_008851.3:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_006920.4']['transcript_variant_error'] is None
        assert 'NM_001353961.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353961.1']['t_hgvs'] == 'NM_001353961.1:c.-2193_-2186del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353961.1']['p_hgvs_tlc'] == 'NP_001340890.1:p.?'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353961.1']['p_hgvs_slc'] == 'NP_001340890.1:p.?'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353961.1']['transcript_variant_error'] is None
        assert 'NM_001353949.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353949.1']['t_hgvs'] == 'NM_001353949.1:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353949.1']['p_hgvs_tlc'] == 'NP_001340878.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353949.1']['p_hgvs_slc'] == 'NP_001340878.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353949.1']['transcript_variant_error'] is None
        assert 'NM_001353960.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353960.1']['t_hgvs'] == 'NM_001353960.1:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353960.1']['p_hgvs_tlc'] == 'NP_001340889.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353960.1']['p_hgvs_slc'] == 'NP_001340889.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353960.1']['transcript_variant_error'] is None
        assert 'NM_006920.5' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_006920.5']['t_hgvs'] == 'NM_006920.5:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_006920.5']['p_hgvs_tlc'] == 'NP_008851.3:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_006920.5']['p_hgvs_slc'] == 'NP_008851.3:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_006920.5']['transcript_variant_error'] is None
        assert 'NR_148667.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NR_148667.1']['t_hgvs'] == 'NR_148667.1:n.638_645del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NR_148667.1']['p_hgvs_tlc'] is None
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NR_148667.1']['p_hgvs_slc'] is None
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NR_148667.1']['transcript_variant_error'] is None
        assert 'NM_001165963.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165963.1']['t_hgvs'] == 'NM_001165963.1:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165963.1']['p_hgvs_tlc'] == 'NP_001159435.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165963.1']['p_hgvs_slc'] == 'NP_001159435.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165963.1']['transcript_variant_error'] is None
        assert 'NM_001202435.2' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001202435.2']['t_hgvs'] == 'NM_001202435.2:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001202435.2']['p_hgvs_tlc'] == 'NP_001189364.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001202435.2']['p_hgvs_slc'] == 'NP_001189364.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001202435.2']['transcript_variant_error'] is None
        assert 'NM_001353957.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353957.1']['t_hgvs'] == 'NM_001353957.1:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353957.1']['p_hgvs_tlc'] == 'NP_001340886.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353957.1']['p_hgvs_slc'] == 'NP_001340886.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001353957.1']['transcript_variant_error'] is None
        assert 'NM_001165963.2' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165963.2']['t_hgvs'] == 'NM_001165963.2:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165963.2']['p_hgvs_tlc'] == 'NP_001159435.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165963.2']['p_hgvs_slc'] == 'NP_001159435.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165963.2']['transcript_variant_error'] is None
        assert 'NM_001202435.1' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001202435.1']['t_hgvs'] == 'NM_001202435.1:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001202435.1']['p_hgvs_tlc'] == 'NP_001189364.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001202435.1']['p_hgvs_slc'] == 'NP_001189364.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001202435.1']['transcript_variant_error'] is None
        assert 'NM_001165964.2' in results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p'].keys()
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165964.2']['t_hgvs'] == 'NM_001165964.2:c.233_240del'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165964.2']['p_hgvs_tlc'] == 'NP_001159436.1:p.(Glu78GlyfsTer7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165964.2']['p_hgvs_slc'] == 'NP_001159436.1:p.(E78Gfs*7)'
        assert results['2-166929891-CCAGGTCCT-C']['hgvs_t_and_p']['NM_001165964.2']['transcript_variant_error'] is None

    def test_variant146(self):
        variant = '2-179393504-G-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '2-179393504-G-T' in results.keys()
        assert results['2-179393504-G-T']['p_vcf'] == '2-179393504-G-T'
        assert results['2-179393504-G-T']['g_hgvs'] == 'NC_000002.11:g.179393504G>T'
        assert results['2-179393504-G-T']['genomic_variant_error'] is None
        assert 'NM_001256850.1' in results['2-179393504-G-T']['hgvs_t_and_p'].keys()
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_001256850.1']['t_hgvs'] == 'NM_001256850.1:c.102051C>A'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_001256850.1']['p_hgvs_tlc'] == 'NP_001243779.1:p.(Ser34017Arg)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_001256850.1']['p_hgvs_slc'] == 'NP_001243779.1:p.(S34017R)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_001256850.1']['transcript_variant_error'] is None
        assert 'NM_133437.3' in results['2-179393504-G-T']['hgvs_t_and_p'].keys()
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133437.3']['t_hgvs'] == 'NM_133437.3:c.80355C>A'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133437.3']['p_hgvs_tlc'] == 'NP_597681.3:p.(Ser26785Arg)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133437.3']['p_hgvs_slc'] == 'NP_597681.3:p.(S26785R)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133437.3']['transcript_variant_error'] is None
        assert 'NM_133437.4' in results['2-179393504-G-T']['hgvs_t_and_p'].keys()
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133437.4']['t_hgvs'] == 'NM_133437.4:c.80355C>A'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133437.4']['p_hgvs_tlc'] == 'NP_597681.4:p.(Ser26785Arg)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133437.4']['p_hgvs_slc'] == 'NP_597681.4:p.(S26785R)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133437.4']['transcript_variant_error'] is None
        assert 'NM_001267550.1' in results['2-179393504-G-T']['hgvs_t_and_p'].keys()
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_001267550.1']['t_hgvs'] == 'NM_001267550.1:c.106974C>A'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_001267550.1']['p_hgvs_tlc'] == 'NP_001254479.1:p.(Ser35658Arg)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_001267550.1']['p_hgvs_slc'] == 'NP_001254479.1:p.(S35658R)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_001267550.1']['transcript_variant_error'] is None
        assert 'NM_001267550.2' in results['2-179393504-G-T']['hgvs_t_and_p'].keys()
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_001267550.2']['t_hgvs'] == 'NM_001267550.2:c.106974C>A'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_001267550.2']['p_hgvs_tlc'] == 'NP_001254479.2:p.(Ser35658Arg)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_001267550.2']['p_hgvs_slc'] == 'NP_001254479.2:p.(S35658R)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_001267550.2']['transcript_variant_error'] is None
        assert 'NM_133378.4' in results['2-179393504-G-T']['hgvs_t_and_p'].keys()
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133378.4']['t_hgvs'] == 'NM_133378.4:c.99270C>A'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133378.4']['p_hgvs_tlc'] == 'NP_596869.4:p.(Ser33090Arg)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133378.4']['p_hgvs_slc'] == 'NP_596869.4:p.(S33090R)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133378.4']['transcript_variant_error'] is None
        assert 'NR_038272.1' in results['2-179393504-G-T']['hgvs_t_and_p'].keys()
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NR_038272.1']['t_hgvs'] == 'NR_038272.1:n.219+5141G>T'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NR_038272.1']['p_hgvs_tlc'] is None
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NR_038272.1']['p_hgvs_slc'] is None
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NR_038272.1']['transcript_variant_error'] is None
        assert 'NM_003319.4' in results['2-179393504-G-T']['hgvs_t_and_p'].keys()
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_003319.4']['t_hgvs'] == 'NM_003319.4:c.79779C>A'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_003319.4']['p_hgvs_tlc'] == 'NP_003310.4:p.(Ser26593Arg)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_003319.4']['p_hgvs_slc'] == 'NP_003310.4:p.(S26593R)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_003319.4']['transcript_variant_error'] is None
        assert 'NR_038271.1' in results['2-179393504-G-T']['hgvs_t_and_p'].keys()
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NR_038271.1']['t_hgvs'] == 'NR_038271.1:n.446+5141G>T'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NR_038271.1']['p_hgvs_tlc'] is None
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NR_038271.1']['p_hgvs_slc'] is None
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NR_038271.1']['transcript_variant_error'] is None
        assert 'NM_133432.3' in results['2-179393504-G-T']['hgvs_t_and_p'].keys()
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133432.3']['t_hgvs'] == 'NM_133432.3:c.80154C>A'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133432.3']['p_hgvs_tlc'] == 'NP_597676.3:p.(Ser26718Arg)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133432.3']['p_hgvs_slc'] == 'NP_597676.3:p.(S26718R)'
        assert results['2-179393504-G-T']['hgvs_t_and_p']['NM_133432.3']['transcript_variant_error'] is None

    def test_variant147(self):
        variant = '2-185803444-TGCAGCTGCTGCAGCTGCAGCTGCA-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '2-185803444-TGCAGCTGCTGCAGCTGCAGCTGCA-T' in results.keys()
        assert results['2-185803444-TGCAGCTGCTGCAGCTGCAGCTGCA-T']['p_vcf'] == '2-185803444-TGCAGCTGCTGCAGCTGCAGCTGCA-T'
        assert results['2-185803444-TGCAGCTGCTGCAGCTGCAGCTGCA-T']['g_hgvs'] == 'NC_000002.11:g.185803447_185803470del'
        assert results['2-185803444-TGCAGCTGCTGCAGCTGCAGCTGCA-T']['genomic_variant_error'] is None
        assert 'NM_194250.1' in results['2-185803444-TGCAGCTGCTGCAGCTGCAGCTGCA-T']['hgvs_t_and_p'].keys()
        assert results['2-185803444-TGCAGCTGCTGCAGCTGCAGCTGCA-T']['hgvs_t_and_p']['NM_194250.1']['t_hgvs'] == 'NM_194250.1:c.3324_3347del'
        assert results['2-185803444-TGCAGCTGCTGCAGCTGCAGCTGCA-T']['hgvs_t_and_p']['NM_194250.1']['p_hgvs_tlc'] == 'NP_919226.1:p.(Ala1112_Ala1119del)'
        assert results['2-185803444-TGCAGCTGCTGCAGCTGCAGCTGCA-T']['hgvs_t_and_p']['NM_194250.1']['p_hgvs_slc'] == 'NP_919226.1:p.(A1112_A1119del)'
        assert results['2-185803444-TGCAGCTGCTGCAGCTGCAGCTGCA-T']['hgvs_t_and_p']['NM_194250.1']['transcript_variant_error'] is None

    def test_variant148(self):
        variant = '2-201950249-G-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '2-201950249-G-T' in results.keys()
        assert results['2-201950249-G-T']['p_vcf'] == '2-201950249-G-T'
        assert results['2-201950249-G-T']['g_hgvs'] == 'NC_000002.11:g.201950249G>T'
        assert results['2-201950249-G-T']['genomic_variant_error'] is None
        assert 'NM_002491.2' in results['2-201950249-G-T']['hgvs_t_and_p'].keys()
        assert results['2-201950249-G-T']['hgvs_t_and_p']['NM_002491.2']['t_hgvs'] == 'NM_002491.2:c.208G>T'
        assert results['2-201950249-G-T']['hgvs_t_and_p']['NM_002491.2']['p_hgvs_tlc'] == 'NP_002482.1:p.(Gly70Ter)'
        assert results['2-201950249-G-T']['hgvs_t_and_p']['NM_002491.2']['p_hgvs_slc'] == 'NP_002482.1:p.(G70*)'
        assert results['2-201950249-G-T']['hgvs_t_and_p']['NM_002491.2']['transcript_variant_error'] is None
        assert 'NM_001257102.1' in results['2-201950249-G-T']['hgvs_t_and_p'].keys()
        assert results['2-201950249-G-T']['hgvs_t_and_p']['NM_001257102.1']['t_hgvs'] == 'NM_001257102.1:c.208G>T'
        assert results['2-201950249-G-T']['hgvs_t_and_p']['NM_001257102.1']['p_hgvs_tlc'] == 'NP_001244031.1:p.(Gly70Ter)'
        assert results['2-201950249-G-T']['hgvs_t_and_p']['NM_001257102.1']['p_hgvs_slc'] == 'NP_001244031.1:p.(G70*)'
        assert results['2-201950249-G-T']['hgvs_t_and_p']['NM_001257102.1']['transcript_variant_error'] is None

    def test_variant149(self):
        variant = '2-238268730-C-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '2-238268730-C-A' in results.keys()
        assert results['2-238268730-C-A']['p_vcf'] == '2-238268730-C-A'
        assert results['2-238268730-C-A']['g_hgvs'] == 'NC_000002.11:g.238268730C>A'
        assert results['2-238268730-C-A']['genomic_variant_error'] is None
        assert 'NM_004369.3' in results['2-238268730-C-A']['hgvs_t_and_p'].keys()
        assert results['2-238268730-C-A']['hgvs_t_and_p']['NM_004369.3']['t_hgvs'] == 'NM_004369.3:c.6282+1G>T'
        assert results['2-238268730-C-A']['hgvs_t_and_p']['NM_004369.3']['p_hgvs_tlc'] == 'NP_004360.2:p.?'
        assert results['2-238268730-C-A']['hgvs_t_and_p']['NM_004369.3']['p_hgvs_slc'] == 'NP_004360.2:p.?'
        assert results['2-238268730-C-A']['hgvs_t_and_p']['NM_004369.3']['transcript_variant_error'] is None
        assert 'NM_057167.3' in results['2-238268730-C-A']['hgvs_t_and_p'].keys()
        assert results['2-238268730-C-A']['hgvs_t_and_p']['NM_057167.3']['t_hgvs'] == 'NM_057167.3:c.5664+1G>T'
        assert results['2-238268730-C-A']['hgvs_t_and_p']['NM_057167.3']['p_hgvs_tlc'] == 'NP_476508.2:p.?'
        assert results['2-238268730-C-A']['hgvs_t_and_p']['NM_057167.3']['p_hgvs_slc'] == 'NP_476508.2:p.?'
        assert results['2-238268730-C-A']['hgvs_t_and_p']['NM_057167.3']['transcript_variant_error'] is None
        assert 'NM_057166.4' in results['2-238268730-C-A']['hgvs_t_and_p'].keys()
        assert results['2-238268730-C-A']['hgvs_t_and_p']['NM_057166.4']['t_hgvs'] == 'NM_057166.4:c.4461+1G>T'
        assert results['2-238268730-C-A']['hgvs_t_and_p']['NM_057166.4']['p_hgvs_tlc'] == 'NP_476507.3:p.?'
        assert results['2-238268730-C-A']['hgvs_t_and_p']['NM_057166.4']['p_hgvs_slc'] == 'NP_476507.3:p.?'
        assert results['2-238268730-C-A']['hgvs_t_and_p']['NM_057166.4']['transcript_variant_error'] is None

    def test_variant150(self):
        variant = '21-43897396-C-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '21-43897396-C-T' in results.keys()
        assert results['21-43897396-C-T']['p_vcf'] == '21-43897396-C-T'
        assert results['21-43897396-C-T']['g_hgvs'] == 'NC_000021.8:g.43897396C>T'
        assert results['21-43897396-C-T']['genomic_variant_error'] is None
        assert 'NM_080860.3' in results['21-43897396-C-T']['hgvs_t_and_p'].keys()
        assert results['21-43897396-C-T']['hgvs_t_and_p']['NM_080860.3']['t_hgvs'] == 'NM_080860.3:c.727+5G>A'
        assert results['21-43897396-C-T']['hgvs_t_and_p']['NM_080860.3']['p_hgvs_tlc'] == 'NP_543136.1:p.?'
        assert results['21-43897396-C-T']['hgvs_t_and_p']['NM_080860.3']['p_hgvs_slc'] == 'NP_543136.1:p.?'
        assert results['21-43897396-C-T']['hgvs_t_and_p']['NM_080860.3']['transcript_variant_error'] is None
        assert 'NM_080860.2' in results['21-43897396-C-T']['hgvs_t_and_p'].keys()
        assert results['21-43897396-C-T']['hgvs_t_and_p']['NM_080860.2']['t_hgvs'] == 'NM_080860.2:c.727+5G>A'
        assert results['21-43897396-C-T']['hgvs_t_and_p']['NM_080860.2']['p_hgvs_tlc'] == 'NP_543136.1:p.?'
        assert results['21-43897396-C-T']['hgvs_t_and_p']['NM_080860.2']['p_hgvs_slc'] == 'NP_543136.1:p.?'
        assert results['21-43897396-C-T']['hgvs_t_and_p']['NM_080860.2']['transcript_variant_error'] is None
        assert 'NM_001286506.1' in results['21-43897396-C-T']['hgvs_t_and_p'].keys()
        assert results['21-43897396-C-T']['hgvs_t_and_p']['NM_001286506.1']['t_hgvs'] == 'NM_001286506.1:c.613+5G>A'
        assert results['21-43897396-C-T']['hgvs_t_and_p']['NM_001286506.1']['p_hgvs_tlc'] == 'NP_001273435.1:p.?'
        assert results['21-43897396-C-T']['hgvs_t_and_p']['NM_001286506.1']['p_hgvs_slc'] == 'NP_001273435.1:p.?'
        assert results['21-43897396-C-T']['hgvs_t_and_p']['NM_001286506.1']['transcript_variant_error'] is None

    def test_variant151(self):
        variant = '22-30064360-G-GCGACGC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '22-30064360-G-GCGACGC' in results.keys()
        assert results['22-30064360-G-GCGACGC']['p_vcf'] == '22-30064360-G-GCGACGC'
        assert results['22-30064360-G-GCGACGC']['g_hgvs'] == 'NC_000022.10:g.30064360_30064361insCGACGC'
        assert results['22-30064360-G-GCGACGC']['genomic_variant_error'] is None
        assert 'NM_181830.2' in results['22-30064360-G-GCGACGC']['hgvs_t_and_p'].keys()
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181830.2']['t_hgvs'] == 'NM_181830.2:c.675_676insCGACGC'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181830.2']['p_hgvs_tlc'] == 'NP_861968.1:p.(Arg227_Arg228dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181830.2']['p_hgvs_slc'] == 'NP_861968.1:p.(R227_R228dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181830.2']['transcript_variant_error'] is None
        assert 'NM_181831.2' in results['22-30064360-G-GCGACGC']['hgvs_t_and_p'].keys()
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181831.2']['t_hgvs'] == 'NM_181831.2:c.675_676insCGACGC'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181831.2']['p_hgvs_tlc'] == 'NP_861969.1:p.(Arg227_Arg228dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181831.2']['p_hgvs_slc'] == 'NP_861969.1:p.(R227_R228dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181831.2']['transcript_variant_error'] is None
        assert 'NR_156186.1' in results['22-30064360-G-GCGACGC']['hgvs_t_and_p'].keys()
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NR_156186.1']['t_hgvs'] == 'NR_156186.1:n.1483_1484insCGACGC'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NR_156186.1']['p_hgvs_tlc'] is None
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NR_156186.1']['p_hgvs_slc'] is None
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NR_156186.1']['transcript_variant_error'] is None
        assert 'NM_181829.2' in results['22-30064360-G-GCGACGC']['hgvs_t_and_p'].keys()
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181829.2']['t_hgvs'] == 'NM_181829.2:c.801_802insCGACGC'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181829.2']['p_hgvs_tlc'] == 'NP_861967.1:p.(Arg269_Arg270dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181829.2']['p_hgvs_slc'] == 'NP_861967.1:p.(R269_R270dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181829.2']['transcript_variant_error'] is None
        assert 'NM_016418.5' in results['22-30064360-G-GCGACGC']['hgvs_t_and_p'].keys()
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_016418.5']['t_hgvs'] == 'NM_016418.5:c.924_925insCGACGC'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_016418.5']['p_hgvs_tlc'] == 'NP_057502.2:p.(Arg310_Arg311dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_016418.5']['p_hgvs_slc'] == 'NP_057502.2:p.(R310_R311dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_016418.5']['transcript_variant_error'] is None
        assert 'NM_181828.2' in results['22-30064360-G-GCGACGC']['hgvs_t_and_p'].keys()
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181828.2']['t_hgvs'] == 'NM_181828.2:c.798_799insCGACGC'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181828.2']['p_hgvs_tlc'] == 'NP_861966.1:p.(Arg268_Arg269dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181828.2']['p_hgvs_slc'] == 'NP_861966.1:p.(R268_R269dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181828.2']['transcript_variant_error'] is None
        assert 'NM_000268.3' in results['22-30064360-G-GCGACGC']['hgvs_t_and_p'].keys()
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_000268.3']['t_hgvs'] == 'NM_000268.3:c.924_925insCGACGC'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_000268.3']['p_hgvs_tlc'] == 'NP_000259.1:p.(Arg310_Arg311dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_000268.3']['p_hgvs_slc'] == 'NP_000259.1:p.(R310_R311dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_000268.3']['transcript_variant_error'] is None
        assert 'NM_181832.2' in results['22-30064360-G-GCGACGC']['hgvs_t_and_p'].keys()
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181832.2']['t_hgvs'] == 'NM_181832.2:c.924_925insCGACGC'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181832.2']['p_hgvs_tlc'] == 'NP_861970.1:p.(Arg310_Arg311dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181832.2']['p_hgvs_slc'] == 'NP_861970.1:p.(R310_R311dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181832.2']['transcript_variant_error'] is None
        assert 'NM_181825.2' in results['22-30064360-G-GCGACGC']['hgvs_t_and_p'].keys()
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181825.2']['t_hgvs'] == 'NM_181825.2:c.924_925insCGACGC'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181825.2']['p_hgvs_tlc'] == 'NP_861546.1:p.(Arg310_Arg311dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181825.2']['p_hgvs_slc'] == 'NP_861546.1:p.(R310_R311dup)'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181825.2']['transcript_variant_error'] is None
        assert 'NM_181833.2' in results['22-30064360-G-GCGACGC']['hgvs_t_and_p'].keys()
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181833.2']['t_hgvs'] == 'NM_181833.2:c.447+26086_447+26087insCGACGC'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181833.2']['p_hgvs_tlc'] == 'NP_861971.1:p.?'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181833.2']['p_hgvs_slc'] == 'NP_861971.1:p.?'
        assert results['22-30064360-G-GCGACGC']['hgvs_t_and_p']['NM_181833.2']['transcript_variant_error'] is None

    def test_variant152(self):
        variant = '3-10188187-TGTCCCGATAG-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '3-10188187-TGTCCCGATAG-T' in results.keys()
        assert results['3-10188187-TGTCCCGATAG-T']['p_vcf'] == '3-10188187-TGTCCCGATAG-T'
        assert results['3-10188187-TGTCCCGATAG-T']['g_hgvs'] == 'NC_000003.11:g.10188191_10188200del'
        assert results['3-10188187-TGTCCCGATAG-T']['genomic_variant_error'] is None
        assert 'NM_001354723.1' in results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p'].keys()
        assert results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p']['NM_001354723.1']['t_hgvs'] == 'NM_001354723.1:c.*18-3280_*18-3271del'
        assert results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p']['NM_001354723.1']['p_hgvs_tlc'] == 'NP_001341652.1:p.?'
        assert results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p']['NM_001354723.1']['p_hgvs_slc'] == 'NP_001341652.1:p.?'
        assert results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p']['NM_001354723.1']['transcript_variant_error'] is None
        assert 'NM_198156.2' in results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p'].keys()
        assert results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p']['NM_198156.2']['t_hgvs'] == 'NM_198156.2:c.341-3280_341-3271del'
        assert results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p']['NM_198156.2']['p_hgvs_tlc'] == 'NP_937799.1:p.?'
        assert results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p']['NM_198156.2']['p_hgvs_slc'] == 'NP_937799.1:p.?'
        assert results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p']['NM_198156.2']['transcript_variant_error'] is None
        assert 'NM_000551.3' in results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p'].keys()
        assert results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p']['NM_000551.3']['t_hgvs'] == 'NM_000551.3:c.341-7_343del'
        assert results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p']['NM_000551.3']['p_hgvs_tlc'] == 'NP_000542.1:p.?'
        assert results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p']['NM_000551.3']['p_hgvs_slc'] == 'NP_000542.1:p.?'
        assert results['3-10188187-TGTCCCGATAG-T']['hgvs_t_and_p']['NM_000551.3']['transcript_variant_error'] is None

    def test_variant153(self):
        variant = '3-50402127-T-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '3-50402127-T-G' in results.keys()
        assert results['3-50402127-T-G']['p_vcf'] == '3-50402127-T-G'
        assert results['3-50402127-T-G']['g_hgvs'] == 'NC_000003.11:g.50402127T>G'
        assert results['3-50402127-T-G']['genomic_variant_error'] is None
        assert 'NM_001291101.1' in results['3-50402127-T-G']['hgvs_t_and_p'].keys()
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001291101.1']['t_hgvs'] == 'NM_001291101.1:c.3201A>C'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001291101.1']['p_hgvs_tlc'] == 'NP_001278030.1:p.(Gln1067His)'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001291101.1']['p_hgvs_slc'] == 'NP_001278030.1:p.(Q1067H)'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001291101.1']['transcript_variant_error'] is None
        assert 'NR_111912.1' in results['3-50402127-T-G']['hgvs_t_and_p'].keys()
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NR_111912.1']['t_hgvs'] == 'NR_111912.1:n.443-1601T>G'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NR_111912.1']['p_hgvs_tlc'] is None
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NR_111912.1']['p_hgvs_slc'] is None
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NR_111912.1']['transcript_variant_error'] is None
        assert 'NM_001005505.2' in results['3-50402127-T-G']['hgvs_t_and_p'].keys()
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001005505.2']['t_hgvs'] == 'NM_001005505.2:c.3408A>C'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001005505.2']['p_hgvs_tlc'] == 'NP_001005505.1:p.(Gln1136His)'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001005505.2']['p_hgvs_slc'] == 'NP_001005505.1:p.(Q1136H)'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001005505.2']['transcript_variant_error'] is None
        assert 'NM_001174051.1' in results['3-50402127-T-G']['hgvs_t_and_p'].keys()
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001174051.1']['t_hgvs'] == 'NM_001174051.1:c.3423A>C'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001174051.1']['p_hgvs_tlc'] == 'NP_001167522.1:p.(Gln1141His)'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001174051.1']['p_hgvs_slc'] == 'NP_001167522.1:p.(Q1141H)'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001174051.1']['transcript_variant_error'] is None
        assert 'NM_001174051.2' in results['3-50402127-T-G']['hgvs_t_and_p'].keys()
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001174051.2']['t_hgvs'] == 'NM_001174051.2:c.3423A>C'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001174051.2']['p_hgvs_tlc'] == 'NP_001167522.1:p.(Gln1141His)'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001174051.2']['p_hgvs_slc'] == 'NP_001167522.1:p.(Q1141H)'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001174051.2']['transcript_variant_error'] is None
        assert 'NM_001005505.1' in results['3-50402127-T-G']['hgvs_t_and_p'].keys()
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001005505.1']['t_hgvs'] == 'NM_001005505.1:c.3408A>C'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001005505.1']['p_hgvs_tlc'] == 'NP_001005505.1:p.(Gln1136His)'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001005505.1']['p_hgvs_slc'] == 'NP_001005505.1:p.(Q1136H)'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_001005505.1']['transcript_variant_error'] is None
        assert 'NM_006030.2' in results['3-50402127-T-G']['hgvs_t_and_p'].keys()
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_006030.2']['t_hgvs'] == 'NM_006030.2:c.3402A>C'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_006030.2']['p_hgvs_tlc'] == 'NP_006021.2:p.(Gln1134His)'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_006030.2']['p_hgvs_slc'] == 'NP_006021.2:p.(Q1134H)'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_006030.2']['transcript_variant_error'] is None
        assert 'NM_006030.3' in results['3-50402127-T-G']['hgvs_t_and_p'].keys()
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_006030.3']['t_hgvs'] == 'NM_006030.3:c.3402A>C'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_006030.3']['p_hgvs_tlc'] == 'NP_006021.2:p.(Gln1134His)'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_006030.3']['p_hgvs_slc'] == 'NP_006021.2:p.(Q1134H)'
        assert results['3-50402127-T-G']['hgvs_t_and_p']['NM_006030.3']['transcript_variant_error'] is None

    def test_variant154(self):
        variant = '3-50402890-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '3-50402890-G-A' in results.keys()
        assert results['3-50402890-G-A']['p_vcf'] == '3-50402890-G-A'
        assert results['3-50402890-G-A']['g_hgvs'] == 'NC_000003.11:g.50402890G>A'
        assert results['3-50402890-G-A']['genomic_variant_error'] is None
        assert 'NM_001174051.2' in results['3-50402890-G-A']['hgvs_t_and_p'].keys()
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001174051.2']['t_hgvs'] == 'NM_001174051.2:c.3016C>T'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001174051.2']['p_hgvs_tlc'] == 'NP_001167522.1:p.(Pro1006Ser)'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001174051.2']['p_hgvs_slc'] == 'NP_001167522.1:p.(P1006S)'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001174051.2']['transcript_variant_error'] is None
        assert 'NR_111914.1' in results['3-50402890-G-A']['hgvs_t_and_p'].keys()
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NR_111914.1']['t_hgvs'] == 'NR_111914.1:n.126G>A'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NR_111914.1']['p_hgvs_tlc'] is None
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NR_111914.1']['p_hgvs_slc'] is None
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NR_111914.1']['transcript_variant_error'] is None
        assert 'NM_001291101.1' in results['3-50402890-G-A']['hgvs_t_and_p'].keys()
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001291101.1']['t_hgvs'] == 'NM_001291101.1:c.2788C>T'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001291101.1']['p_hgvs_tlc'] == 'NP_001278030.1:p.(Pro930Ser)'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001291101.1']['p_hgvs_slc'] == 'NP_001278030.1:p.(P930S)'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001291101.1']['transcript_variant_error'] is None
        assert 'NR_111912.1' in results['3-50402890-G-A']['hgvs_t_and_p'].keys()
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NR_111912.1']['t_hgvs'] == 'NR_111912.1:n.443-838G>A'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NR_111912.1']['p_hgvs_tlc'] is None
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NR_111912.1']['p_hgvs_slc'] is None
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NR_111912.1']['transcript_variant_error'] is None
        assert 'NM_001005505.2' in results['3-50402890-G-A']['hgvs_t_and_p'].keys()
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001005505.2']['t_hgvs'] == 'NM_001005505.2:c.2995C>T'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001005505.2']['p_hgvs_tlc'] == 'NP_001005505.1:p.(Pro999Ser)'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001005505.2']['p_hgvs_slc'] == 'NP_001005505.1:p.(P999S)'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001005505.2']['transcript_variant_error'] is None
        assert 'NM_001174051.1' in results['3-50402890-G-A']['hgvs_t_and_p'].keys()
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001174051.1']['t_hgvs'] == 'NM_001174051.1:c.3016C>T'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001174051.1']['p_hgvs_tlc'] == 'NP_001167522.1:p.(Pro1006Ser)'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001174051.1']['p_hgvs_slc'] == 'NP_001167522.1:p.(P1006S)'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001174051.1']['transcript_variant_error'] is None
        assert 'NR_111913.1' in results['3-50402890-G-A']['hgvs_t_and_p'].keys()
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NR_111913.1']['t_hgvs'] == 'NR_111913.1:n.126G>A'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NR_111913.1']['p_hgvs_tlc'] is None
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NR_111913.1']['p_hgvs_slc'] is None
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NR_111913.1']['transcript_variant_error'] is None
        assert 'NM_001005505.1' in results['3-50402890-G-A']['hgvs_t_and_p'].keys()
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001005505.1']['t_hgvs'] == 'NM_001005505.1:c.2995C>T'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001005505.1']['p_hgvs_tlc'] == 'NP_001005505.1:p.(Pro999Ser)'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001005505.1']['p_hgvs_slc'] == 'NP_001005505.1:p.(P999S)'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_001005505.1']['transcript_variant_error'] is None
        assert 'NM_006030.2' in results['3-50402890-G-A']['hgvs_t_and_p'].keys()
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_006030.2']['t_hgvs'] == 'NM_006030.2:c.2995C>T'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_006030.2']['p_hgvs_tlc'] == 'NP_006021.2:p.(Pro999Ser)'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_006030.2']['p_hgvs_slc'] == 'NP_006021.2:p.(P999S)'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_006030.2']['transcript_variant_error'] is None
        assert 'NM_006030.3' in results['3-50402890-G-A']['hgvs_t_and_p'].keys()
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_006030.3']['t_hgvs'] == 'NM_006030.3:c.2995C>T'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_006030.3']['p_hgvs_tlc'] == 'NP_006021.2:p.(Pro999Ser)'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_006030.3']['p_hgvs_slc'] == 'NP_006021.2:p.(P999S)'
        assert results['3-50402890-G-A']['hgvs_t_and_p']['NM_006030.3']['transcript_variant_error'] is None

    def test_variant155(self):
        variant = '3-57851007-AG-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '3-57851007-AG-A' in results.keys()
        assert results['3-57851007-AG-A']['p_vcf'] == '3-57851007-AG-A'
        assert results['3-57851007-AG-A']['g_hgvs'] == 'NC_000003.11:g.57851008del'
        assert results['3-57851007-AG-A']['genomic_variant_error'] is None
        assert 'NM_007159.4' in results['3-57851007-AG-A']['hgvs_t_and_p'].keys()
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_007159.4']['t_hgvs'] == 'NM_007159.4:c.1135+565del'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_007159.4']['p_hgvs_tlc'] == 'NP_009090.2:p.?'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_007159.4']['p_hgvs_slc'] == 'NP_009090.2:p.?'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_007159.4']['transcript_variant_error'] is None
        assert 'NM_001304420.1' in results['3-57851007-AG-A']['hgvs_t_and_p'].keys()
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304420.1']['t_hgvs'] == 'NM_001304420.1:c.1186+424del'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304420.1']['p_hgvs_tlc'] == 'NP_001291349.1:p.?'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304420.1']['p_hgvs_slc'] == 'NP_001291349.1:p.?'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304420.1']['transcript_variant_error'] is None
        assert 'NM_007159.3' in results['3-57851007-AG-A']['hgvs_t_and_p'].keys()
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_007159.3']['t_hgvs'] == 'NM_007159.3:c.1135+565del'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_007159.3']['p_hgvs_tlc'] == 'NP_009090.2:p.?'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_007159.3']['p_hgvs_slc'] == 'NP_009090.2:p.?'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_007159.3']['transcript_variant_error'] is None
        assert 'NM_007159.2' in results['3-57851007-AG-A']['hgvs_t_and_p'].keys()
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_007159.2']['t_hgvs'] == 'NM_007159.2:c.1135+565del'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_007159.2']['p_hgvs_tlc'] == 'NP_009090.2:p.?'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_007159.2']['p_hgvs_slc'] == 'NP_009090.2:p.?'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_007159.2']['transcript_variant_error'] is None
        assert 'NM_001304421.1' in results['3-57851007-AG-A']['hgvs_t_and_p'].keys()
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304421.1']['t_hgvs'] == 'NM_001304421.1:c.1135+565del'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304421.1']['p_hgvs_tlc'] == 'NP_001291350.1:p.?'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304421.1']['p_hgvs_slc'] == 'NP_001291350.1:p.?'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304421.1']['transcript_variant_error'] is None
        assert 'NM_001304420.2' in results['3-57851007-AG-A']['hgvs_t_and_p'].keys()
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304420.2']['t_hgvs'] == 'NM_001304420.2:c.1186+424del'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304420.2']['p_hgvs_tlc'] == 'NP_001291349.1:p.?'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304420.2']['p_hgvs_slc'] == 'NP_001291349.1:p.?'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304420.2']['transcript_variant_error'] is None
        assert 'NM_001304421.2' in results['3-57851007-AG-A']['hgvs_t_and_p'].keys()
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304421.2']['t_hgvs'] == 'NM_001304421.2:c.1135+565del'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304421.2']['p_hgvs_tlc'] == 'NP_001291350.1:p.?'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304421.2']['p_hgvs_slc'] == 'NP_001291350.1:p.?'
        assert results['3-57851007-AG-A']['hgvs_t_and_p']['NM_001304421.2']['transcript_variant_error'] is None

    def test_variant156(self):
        variant = '3-122003832-G-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '3-122003832-G-C' in results.keys()
        assert results['3-122003832-G-C']['p_vcf'] == '3-122003832-G-C'
        assert results['3-122003832-G-C']['g_hgvs'] == 'NC_000003.11:g.122003832G>C'
        assert results['3-122003832-G-C']['genomic_variant_error'] is None
        assert 'NM_001178065.1' in results['3-122003832-G-C']['hgvs_t_and_p'].keys()
        assert results['3-122003832-G-C']['hgvs_t_and_p']['NM_001178065.1']['t_hgvs'] == 'NM_001178065.1:c.3061='
        assert results['3-122003832-G-C']['hgvs_t_and_p']['NM_001178065.1']['p_hgvs_tlc'] == 'NP_001171536.1:p.(Gln1021=)'
        assert results['3-122003832-G-C']['hgvs_t_and_p']['NM_001178065.1']['p_hgvs_slc'] == 'NP_001171536.1:p.(Q1021=)'
        assert results['3-122003832-G-C']['hgvs_t_and_p']['NM_001178065.1']['transcript_variant_error'] is None
        assert 'NM_000388.3' in results['3-122003832-G-C']['hgvs_t_and_p'].keys()
        assert results['3-122003832-G-C']['hgvs_t_and_p']['NM_000388.3']['t_hgvs'] == 'NM_000388.3:c.3031='
        assert results['3-122003832-G-C']['hgvs_t_and_p']['NM_000388.3']['p_hgvs_tlc'] == 'NP_000379.2:p.(Gln1011=)'
        assert results['3-122003832-G-C']['hgvs_t_and_p']['NM_000388.3']['p_hgvs_slc'] == 'NP_000379.2:p.(Q1011=)'
        assert results['3-122003832-G-C']['hgvs_t_and_p']['NM_000388.3']['transcript_variant_error'] is None

    def test_variant157(self):
        variant = '4-153332910-C-CAGG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '4-153332910-C-CAGG' in results.keys()
        assert results['4-153332910-C-CAGG']['p_vcf'] == '4-153332910-C-CAGG'
        assert results['4-153332910-C-CAGG']['g_hgvs'] == 'NC_000004.11:g.153332912_153332913insGAG'
        assert results['4-153332910-C-CAGG']['genomic_variant_error'] is None
        assert 'NM_001349798.1' in results['4-153332910-C-CAGG']['hgvs_t_and_p'].keys()
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_001349798.1']['t_hgvs'] == 'NM_001349798.1:c.45_46insCCT'
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_001349798.1']['p_hgvs_tlc'] == 'NP_001336727.1:p.(Thr15_Gly16insPro)'
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_001349798.1']['p_hgvs_slc'] == 'NP_001336727.1:p.(T15_G16insP)'
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_001349798.1']['transcript_variant_error'] is None
        assert 'NM_033632.3' in results['4-153332910-C-CAGG']['hgvs_t_and_p'].keys()
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_033632.3']['t_hgvs'] == 'NM_033632.3:c.45_46insCCT'
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_033632.3']['p_hgvs_tlc'] == 'NP_361014.1:p.(Thr15_Gly16insPro)'
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_033632.3']['p_hgvs_slc'] == 'NP_361014.1:p.(T15_G16insP)'
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_033632.3']['transcript_variant_error'] is None
        assert 'NM_001349798.2' in results['4-153332910-C-CAGG']['hgvs_t_and_p'].keys()
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_001349798.2']['t_hgvs'] == 'NM_001349798.2:c.45_46insCCT'
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_001349798.2']['p_hgvs_tlc'] == 'NP_001336727.1:p.(Thr15_Gly16insPro)'
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_001349798.2']['p_hgvs_slc'] == 'NP_001336727.1:p.(T15_G16insP)'
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_001349798.2']['transcript_variant_error'] is None
        assert 'NM_001257069.1' in results['4-153332910-C-CAGG']['hgvs_t_and_p'].keys()
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_001257069.1']['t_hgvs'] == 'NM_001257069.1:c.45_46insCCT'
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_001257069.1']['p_hgvs_tlc'] == 'NP_001243998.1:p.(Thr15_Gly16insPro)'
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_001257069.1']['p_hgvs_slc'] == 'NP_001243998.1:p.(T15_G16insP)'
        assert results['4-153332910-C-CAGG']['hgvs_t_and_p']['NM_001257069.1']['transcript_variant_error'] is None
    # Test removed. No longer intergenic in VVTA
    # def test_variant158(self):
    # 	variant = '5-1295183-G-A'
    # 	results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
    # 	results = results.stucture_data()
    # 	print(results)
    #
    # 	assert '5-1295183-G-A' in results.keys()
    # 	assert results['5-1295183-G-A']['p_vcf'] == '5-1295183-G-A'
    # 	assert results['5-1295183-G-A']['g_hgvs'] == 'NC_000005.9:g.1295183G>A'
    # 	assert results['5-1295183-G-A']['genomic_variant_error'] is None
    # 	assert results['5-1295183-G-A']['hgvs_t_and_p'] == {'intergenic': {'alt_genomic_loci': None}}

    def test_variant159(self):
        variant = '5-77396835-TTTC-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '5-77396835-TTTC-T' in results.keys()
        assert results['5-77396835-TTTC-T']['p_vcf'] == '5-77396835-TTTC-T'
        assert results['5-77396835-TTTC-T']['g_hgvs'] == 'NC_000005.9:g.77396838_77396840del'
        assert results['5-77396835-TTTC-T']['genomic_variant_error'] is None
        assert 'NM_001271769.1' in results['5-77396835-TTTC-T']['hgvs_t_and_p'].keys()
        assert results['5-77396835-TTTC-T']['hgvs_t_and_p']['NM_001271769.1']['t_hgvs'] == 'NM_001271769.1:c.2262_2264del'
        assert results['5-77396835-TTTC-T']['hgvs_t_and_p']['NM_001271769.1']['p_hgvs_tlc'] == 'NP_001258698.1:p.(Lys755del)'
        assert results['5-77396835-TTTC-T']['hgvs_t_and_p']['NM_001271769.1']['p_hgvs_slc'] == 'NP_001258698.1:p.(K755del)'
        assert results['5-77396835-TTTC-T']['hgvs_t_and_p']['NM_001271769.1']['transcript_variant_error'] is None
        assert 'NM_003664.4' in results['5-77396835-TTTC-T']['hgvs_t_and_p'].keys()
        assert results['5-77396835-TTTC-T']['hgvs_t_and_p']['NM_003664.4']['t_hgvs'] == 'NM_003664.4:c.2409_2411del'
        assert results['5-77396835-TTTC-T']['hgvs_t_and_p']['NM_003664.4']['p_hgvs_tlc'] == 'NP_003655.3:p.(Lys804del)'
        assert results['5-77396835-TTTC-T']['hgvs_t_and_p']['NM_003664.4']['p_hgvs_slc'] == 'NP_003655.3:p.(K804del)'
        assert results['5-77396835-TTTC-T']['hgvs_t_and_p']['NM_003664.4']['transcript_variant_error'] is None
        assert 'NM_003664.3' in results['5-77396835-TTTC-T']['hgvs_t_and_p'].keys()
        assert results['5-77396835-TTTC-T']['hgvs_t_and_p']['NM_003664.3']['t_hgvs'] == 'NM_003664.3:c.2409_2411del'
        assert results['5-77396835-TTTC-T']['hgvs_t_and_p']['NM_003664.3']['p_hgvs_tlc'] == 'NP_003655.3:p.(Lys804del)'
        assert results['5-77396835-TTTC-T']['hgvs_t_and_p']['NM_003664.3']['p_hgvs_slc'] == 'NP_003655.3:p.(K804del)'
        assert results['5-77396835-TTTC-T']['hgvs_t_and_p']['NM_003664.3']['transcript_variant_error'] is None

    def test_variant160(self):
        variant = '5-118811422-GGTGA-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '5-118811422-GGTGA-G' in results.keys()
        assert results['5-118811422-GGTGA-G']['p_vcf'] == '5-118811422-GGTGA-G'
        assert results['5-118811422-GGTGA-G']['g_hgvs'] == 'NC_000005.9:g.118811425_118811428del'
        assert results['5-118811422-GGTGA-G']['genomic_variant_error'] is None
        assert 'NM_001292028.1' in results['5-118811422-GGTGA-G']['hgvs_t_and_p'].keys()
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001292028.1']['t_hgvs'] == 'NM_001292028.1:c.-110+3_-110+6del'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001292028.1']['p_hgvs_tlc'] == 'NP_001278957.1:p.?'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001292028.1']['p_hgvs_slc'] == 'NP_001278957.1:p.?'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001292028.1']['transcript_variant_error'] is None
        assert 'NM_001199291.2' in results['5-118811422-GGTGA-G']['hgvs_t_and_p'].keys()
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001199291.2']['t_hgvs'] == 'NM_001199291.2:c.377+3_377+6del'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001199291.2']['p_hgvs_tlc'] == 'NP_001186220.1:p.?'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001199291.2']['p_hgvs_slc'] == 'NP_001186220.1:p.?'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001199291.2']['transcript_variant_error'] is None
        assert 'NM_001199291.1' in results['5-118811422-GGTGA-G']['hgvs_t_and_p'].keys()
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001199291.1']['t_hgvs'] == 'NM_001199291.1:c.377+3_377+6del'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001199291.1']['p_hgvs_tlc'] == 'NP_001186220.1:p.?'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001199291.1']['p_hgvs_slc'] == 'NP_001186220.1:p.?'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001199291.1']['transcript_variant_error'] is None
        assert 'NM_001292027.1' in results['5-118811422-GGTGA-G']['hgvs_t_and_p'].keys()
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001292027.1']['t_hgvs'] == 'NM_001292027.1:c.230+3_230+6del'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001292027.1']['p_hgvs_tlc'] == 'NP_001278956.1:p.?'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001292027.1']['p_hgvs_slc'] == 'NP_001278956.1:p.?'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001292027.1']['transcript_variant_error'] is None
        assert 'NM_001199292.1' in results['5-118811422-GGTGA-G']['hgvs_t_and_p'].keys()
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001199292.1']['t_hgvs'] == 'NM_001199292.1:c.248+3_248+6del'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001199292.1']['p_hgvs_tlc'] == 'NP_001186221.1:p.?'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001199292.1']['p_hgvs_slc'] == 'NP_001186221.1:p.?'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_001199292.1']['transcript_variant_error'] is None
        assert 'NM_000414.3' in results['5-118811422-GGTGA-G']['hgvs_t_and_p'].keys()
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_000414.3']['t_hgvs'] == 'NM_000414.3:c.302+3_302+6del'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_000414.3']['p_hgvs_tlc'] == 'NP_000405.1:p.?'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_000414.3']['p_hgvs_slc'] == 'NP_000405.1:p.?'
        assert results['5-118811422-GGTGA-G']['hgvs_t_and_p']['NM_000414.3']['transcript_variant_error'] is None

    def test_variant161(self):
        variant = '5-118811422-GGTGAG-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '5-118811422-GGTGAG-G' in results.keys()
        assert results['5-118811422-GGTGAG-G']['p_vcf'] == '5-118811421-GGGTGA-G'
        assert results['5-118811422-GGTGAG-G']['g_hgvs'] == 'NC_000005.9:g.118811423_118811427del'
        assert results['5-118811422-GGTGAG-G']['genomic_variant_error'] is None
        assert 'NM_001292028.1' in results['5-118811422-GGTGAG-G']['hgvs_t_and_p'].keys()
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001292028.1']['t_hgvs'] == 'NM_001292028.1:c.-110+1_-110+5del'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001292028.1']['p_hgvs_tlc'] == 'NP_001278957.1:p.?'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001292028.1']['p_hgvs_slc'] == 'NP_001278957.1:p.?'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001292028.1']['transcript_variant_error'] is None
        assert 'NM_001199291.2' in results['5-118811422-GGTGAG-G']['hgvs_t_and_p'].keys()
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001199291.2']['t_hgvs'] == 'NM_001199291.2:c.377+1_377+5del'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001199291.2']['p_hgvs_tlc'] == 'NP_001186220.1:p.?'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001199291.2']['p_hgvs_slc'] == 'NP_001186220.1:p.?'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001199291.2']['transcript_variant_error'] is None
        assert 'NM_001199291.1' in results['5-118811422-GGTGAG-G']['hgvs_t_and_p'].keys()
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001199291.1']['t_hgvs'] == 'NM_001199291.1:c.377+1_377+5del'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001199291.1']['p_hgvs_tlc'] == 'NP_001186220.1:p.?'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001199291.1']['p_hgvs_slc'] == 'NP_001186220.1:p.?'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001199291.1']['transcript_variant_error'] is None
        assert 'NM_001292027.1' in results['5-118811422-GGTGAG-G']['hgvs_t_and_p'].keys()
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001292027.1']['t_hgvs'] == 'NM_001292027.1:c.230+1_230+5del'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001292027.1']['p_hgvs_tlc'] == 'NP_001278956.1:p.?'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001292027.1']['p_hgvs_slc'] == 'NP_001278956.1:p.?'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001292027.1']['transcript_variant_error'] is None
        assert 'NM_001199292.1' in results['5-118811422-GGTGAG-G']['hgvs_t_and_p'].keys()
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001199292.1']['t_hgvs'] == 'NM_001199292.1:c.248+1_248+5del'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001199292.1']['p_hgvs_tlc'] == 'NP_001186221.1:p.?'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001199292.1']['p_hgvs_slc'] == 'NP_001186221.1:p.?'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_001199292.1']['transcript_variant_error'] is None
        assert 'NM_000414.3' in results['5-118811422-GGTGAG-G']['hgvs_t_and_p'].keys()
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_000414.3']['t_hgvs'] == 'NM_000414.3:c.302+1_302+5del'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_000414.3']['p_hgvs_tlc'] == 'NP_000405.1:p.?'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_000414.3']['p_hgvs_slc'] == 'NP_000405.1:p.?'
        assert results['5-118811422-GGTGAG-G']['hgvs_t_and_p']['NM_000414.3']['transcript_variant_error'] is None

    def test_variant162(self):
        variant = '5-131705587-CG-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '5-131705587-CG-C' in results.keys()
        assert results['5-131705587-CG-C']['p_vcf'] == '5-131705587-CG-C'
        assert results['5-131705587-CG-C']['g_hgvs'] == 'NC_000005.9:g.131705590del'
        assert results['5-131705587-CG-C']['genomic_variant_error'] is None
        assert 'NR_110997.1' in results['5-131705587-CG-C']['hgvs_t_and_p'].keys()
        assert results['5-131705587-CG-C']['hgvs_t_and_p']['NR_110997.1']['t_hgvs'] == 'NR_110997.1:n.21del'
        assert results['5-131705587-CG-C']['hgvs_t_and_p']['NR_110997.1']['p_hgvs_tlc'] is None
        assert results['5-131705587-CG-C']['hgvs_t_and_p']['NR_110997.1']['p_hgvs_slc'] is None
        assert results['5-131705587-CG-C']['hgvs_t_and_p']['NR_110997.1']['transcript_variant_error'] is None
        assert 'NM_001308122.1' in results['5-131705587-CG-C']['hgvs_t_and_p'].keys()
        assert results['5-131705587-CG-C']['hgvs_t_and_p']['NM_001308122.1']['t_hgvs'] == 'NM_001308122.1:c.-75del'
        assert results['5-131705587-CG-C']['hgvs_t_and_p']['NM_001308122.1']['p_hgvs_tlc'] == 'NP_001295051.1:p.?'
        assert results['5-131705587-CG-C']['hgvs_t_and_p']['NM_001308122.1']['p_hgvs_slc'] == 'NP_001295051.1:p.?'
        assert results['5-131705587-CG-C']['hgvs_t_and_p']['NM_001308122.1']['transcript_variant_error'] is None
        # Transcript currently deprecated in VVTA
        # assert 'NM_003060.3' in results['5-131705587-CG-C']['hgvs_t_and_p'].keys()
        # assert results['5-131705587-CG-C']['hgvs_t_and_p']['NM_003060.3']['t_hgvs'] == 'NM_003060.3:c.-75del'
        # assert results['5-131705587-CG-C']['hgvs_t_and_p']['NM_003060.3']['p_hgvs_tlc'] == 'NP_003051.1:p.?'
        # assert results['5-131705587-CG-C']['hgvs_t_and_p']['NM_003060.3']['p_hgvs_slc'] == 'NP_003051.1:p.?'
        # assert results['5-131705587-CG-C']['hgvs_t_and_p']['NM_003060.3']['transcript_variant_error'] is None

    def test_variant163(self):
        variant = '5-148406482-T-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '5-148406482-T-C' in results.keys()
        assert results['5-148406482-T-C']['p_vcf'] == '5-148406482-T-C'
        assert results['5-148406482-T-C']['g_hgvs'] == 'NC_000005.9:g.148406482T>C'
        assert results['5-148406482-T-C']['genomic_variant_error'] is None
        assert 'NM_024577.3' in results['5-148406482-T-C']['hgvs_t_and_p'].keys()
        assert results['5-148406482-T-C']['hgvs_t_and_p']['NM_024577.3']['t_hgvs'] == 'NM_024577.3:c.2813A>G'
        assert results['5-148406482-T-C']['hgvs_t_and_p']['NM_024577.3']['p_hgvs_tlc'] == 'NP_078853.2:p.(His938Arg)'
        assert results['5-148406482-T-C']['hgvs_t_and_p']['NM_024577.3']['p_hgvs_slc'] == 'NP_078853.2:p.(H938R)'
        assert results['5-148406482-T-C']['hgvs_t_and_p']['NM_024577.3']['transcript_variant_error'] is None

    def test_variant164(self):
        variant = '6-110036337-T-TCAG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '6-110036337-T-TCAG' in results.keys()
        assert results['6-110036337-T-TCAG']['p_vcf'] == '6-110036337-T-TCAG'
        assert results['6-110036337-T-TCAG']['g_hgvs'] == 'NC_000006.11:g.110036337_110036338insCAG'
        assert results['6-110036337-T-TCAG']['genomic_variant_error'] is None
        assert 'NM_014845.5' in results['6-110036337-T-TCAG']['hgvs_t_and_p'].keys()
        assert results['6-110036337-T-TCAG']['hgvs_t_and_p']['NM_014845.5']['t_hgvs'] == 'NM_014845.5:c.123_124insCAG'
        assert results['6-110036337-T-TCAG']['hgvs_t_and_p']['NM_014845.5']['p_hgvs_tlc'] == 'NP_055660.1:p.(Ile41_Asp42insGln)'
        assert results['6-110036337-T-TCAG']['hgvs_t_and_p']['NM_014845.5']['p_hgvs_slc'] == 'NP_055660.1:p.(I41_D42insQ)'
        assert results['6-110036337-T-TCAG']['hgvs_t_and_p']['NM_014845.5']['transcript_variant_error'] is None

    def test_variant165(self):
        variant = '6-110036337-TGAT-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '6-110036337-TGAT-T' in results.keys()
        assert results['6-110036337-TGAT-T']['p_vcf'] == '6-110036336-TTGA-T'
        assert results['6-110036337-TGAT-T']['g_hgvs'] == 'NC_000006.11:g.110036338_110036340del'
        assert results['6-110036337-TGAT-T']['genomic_variant_error'] is None
        assert 'NM_014845.5' in results['6-110036337-TGAT-T']['hgvs_t_and_p'].keys()
        assert results['6-110036337-TGAT-T']['hgvs_t_and_p']['NM_014845.5']['t_hgvs'] == 'NM_014845.5:c.124_126del'
        assert results['6-110036337-TGAT-T']['hgvs_t_and_p']['NM_014845.5']['p_hgvs_tlc'] == 'NP_055660.1:p.(Asp42del)'
        assert results['6-110036337-TGAT-T']['hgvs_t_and_p']['NM_014845.5']['p_hgvs_slc'] == 'NP_055660.1:p.(D42del)'
        assert results['6-110036337-TGAT-T']['hgvs_t_and_p']['NM_014845.5']['transcript_variant_error'] is None

    def test_variant166(self):
        variant = '6-152651802-C-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '6-152651802-C-A' in results.keys()
        assert results['6-152651802-C-A']['p_vcf'] == '6-152651802-C-A'
        assert results['6-152651802-C-A']['g_hgvs'] == 'NC_000006.11:g.152651802C>A'
        assert results['6-152651802-C-A']['genomic_variant_error'] is None
        assert 'NM_182961.3' in results['6-152651802-C-A']['hgvs_t_and_p'].keys()
        assert results['6-152651802-C-A']['hgvs_t_and_p']['NM_182961.3']['t_hgvs'] == 'NM_182961.3:c.14018G>T'
        assert results['6-152651802-C-A']['hgvs_t_and_p']['NM_182961.3']['p_hgvs_tlc'] == 'NP_892006.3:p.(Arg4673Leu)'
        assert results['6-152651802-C-A']['hgvs_t_and_p']['NM_182961.3']['p_hgvs_slc'] == 'NP_892006.3:p.(R4673L)'
        assert results['6-152651802-C-A']['hgvs_t_and_p']['NM_182961.3']['transcript_variant_error'] is None
        assert 'NM_033071.3' in results['6-152651802-C-A']['hgvs_t_and_p'].keys()
        assert results['6-152651802-C-A']['hgvs_t_and_p']['NM_033071.3']['t_hgvs'] == 'NM_033071.3:c.13805G>T'
        assert results['6-152651802-C-A']['hgvs_t_and_p']['NM_033071.3']['p_hgvs_tlc'] == 'NP_149062.1:p.(Arg4602Leu)'
        assert results['6-152651802-C-A']['hgvs_t_and_p']['NM_033071.3']['p_hgvs_slc'] == 'NP_149062.1:p.(R4602L)'
        assert results['6-152651802-C-A']['hgvs_t_and_p']['NM_033071.3']['transcript_variant_error'] is None

    def test_variant167(self):
        variant = '6-152737643-C-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '6-152737643-C-G' in results.keys()
        assert results['6-152737643-C-G']['p_vcf'] == '6-152737643-C-G'
        assert results['6-152737643-C-G']['g_hgvs'] == 'NC_000006.11:g.152737643C>G'
        assert results['6-152737643-C-G']['genomic_variant_error'] is None
        assert 'NM_182961.3' in results['6-152737643-C-G']['hgvs_t_and_p'].keys()
        assert results['6-152737643-C-G']['hgvs_t_and_p']['NM_182961.3']['t_hgvs'] == 'NM_182961.3:c.5929G>C'
        assert results['6-152737643-C-G']['hgvs_t_and_p']['NM_182961.3']['p_hgvs_tlc'] == 'NP_892006.3:p.(Ala1977Pro)'
        assert results['6-152737643-C-G']['hgvs_t_and_p']['NM_182961.3']['p_hgvs_slc'] == 'NP_892006.3:p.(A1977P)'
        assert results['6-152737643-C-G']['hgvs_t_and_p']['NM_182961.3']['transcript_variant_error'] is None
        assert 'NM_033071.3' in results['6-152737643-C-G']['hgvs_t_and_p'].keys()
        assert results['6-152737643-C-G']['hgvs_t_and_p']['NM_033071.3']['t_hgvs'] == 'NM_033071.3:c.5950G>C'
        assert results['6-152737643-C-G']['hgvs_t_and_p']['NM_033071.3']['p_hgvs_tlc'] == 'NP_149062.1:p.(Ala1984Pro)'
        assert results['6-152737643-C-G']['hgvs_t_and_p']['NM_033071.3']['p_hgvs_slc'] == 'NP_149062.1:p.(A1984P)'
        assert results['6-152737643-C-G']['hgvs_t_and_p']['NM_033071.3']['transcript_variant_error'] is None

    def test_variant168(self):
        variant = '7-6026775-T-C'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '7-6026775-T-C' in results.keys()
        assert results['7-6026775-T-C']['p_vcf'] == '7-6026775-T-C'
        assert results['7-6026775-T-C']['g_hgvs'] == 'NC_000007.13:g.6026775T>C'
        assert results['7-6026775-T-C']['genomic_variant_error'] is None
        assert 'NM_001322007.1' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322007.1']['t_hgvs'] == 'NM_001322007.1:c.1303A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322007.1']['p_hgvs_tlc'] == 'NP_001308936.1:p.(Lys435Glu)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322007.1']['p_hgvs_slc'] == 'NP_001308936.1:p.(K435E)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322007.1']['transcript_variant_error'] is None
        assert 'NM_001322008.1' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322008.1']['t_hgvs'] == 'NM_001322008.1:c.1303A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322008.1']['p_hgvs_tlc'] == 'NP_001308937.1:p.(Lys435Glu)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322008.1']['p_hgvs_slc'] == 'NP_001308937.1:p.(K435E)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322008.1']['transcript_variant_error'] is None
        assert 'NM_001322003.1' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322003.1']['t_hgvs'] == 'NM_001322003.1:c.1216A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322003.1']['p_hgvs_tlc'] == 'NP_001308932.1:p.(Lys406Glu)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322003.1']['p_hgvs_slc'] == 'NP_001308932.1:p.(K406E)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322003.1']['transcript_variant_error'] is None
        assert 'NM_001322011.1' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322011.1']['t_hgvs'] == 'NM_001322011.1:c.688A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322011.1']['p_hgvs_tlc'] == 'NP_001308940.1:p.(Lys230Glu)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322011.1']['p_hgvs_slc'] == 'NP_001308940.1:p.(K230E)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322011.1']['transcript_variant_error'] is None
        assert 'NR_136154.1' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NR_136154.1']['t_hgvs'] == 'NR_136154.1:n.1708A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NR_136154.1']['p_hgvs_tlc'] is None
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NR_136154.1']['p_hgvs_slc'] is None
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NR_136154.1']['transcript_variant_error'] is None
        assert 'NR_003085.2' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NR_003085.2']['t_hgvs'] == 'NR_003085.2:n.1703='
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NR_003085.2']['p_hgvs_tlc'] is None
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NR_003085.2']['p_hgvs_slc'] is None
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NR_003085.2']['transcript_variant_error'] is None
        assert 'NM_001322015.1' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322015.1']['t_hgvs'] == 'NM_001322015.1:c.1312A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322015.1']['p_hgvs_tlc'] == 'NP_001308944.1:p.(Lys438Glu)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322015.1']['p_hgvs_slc'] == 'NP_001308944.1:p.(K438E)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322015.1']['transcript_variant_error'] is None
        assert 'NM_000535.5' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_000535.5']['t_hgvs'] == 'NM_000535.5:c.1621='
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_000535.5']['p_hgvs_tlc'] == 'NP_000526.1:p.(Glu541=)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_000535.5']['p_hgvs_slc'] == 'NP_000526.1:p.(E541=)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_000535.5']['transcript_variant_error'] is None
        assert 'NM_001322014.1' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322014.1']['t_hgvs'] == 'NM_001322014.1:c.1621A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322014.1']['p_hgvs_tlc'] == 'NP_001308943.1:p.(Lys541Glu)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322014.1']['p_hgvs_slc'] == 'NP_001308943.1:p.(K541E)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322014.1']['transcript_variant_error'] is None
        assert 'NM_001322009.1' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322009.1']['t_hgvs'] == 'NM_001322009.1:c.1216A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322009.1']['p_hgvs_tlc'] == 'NP_001308938.1:p.(Lys406Glu)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322009.1']['p_hgvs_slc'] == 'NP_001308938.1:p.(K406E)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322009.1']['transcript_variant_error'] is None
        assert 'NM_001322004.1' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322004.1']['t_hgvs'] == 'NM_001322004.1:c.1216A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322004.1']['p_hgvs_tlc'] == 'NP_001308933.1:p.(Lys406Glu)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322004.1']['p_hgvs_slc'] == 'NP_001308933.1:p.(K406E)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322004.1']['transcript_variant_error'] is None
        assert 'NM_001322010.1' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322010.1']['t_hgvs'] == 'NM_001322010.1:c.1060A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322010.1']['p_hgvs_tlc'] == 'NP_001308939.1:p.(Lys354Glu)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322010.1']['p_hgvs_slc'] == 'NP_001308939.1:p.(K354E)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322010.1']['transcript_variant_error'] is None
        assert 'NM_001322005.1' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322005.1']['t_hgvs'] == 'NM_001322005.1:c.1216A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322005.1']['p_hgvs_tlc'] == 'NP_001308934.1:p.(Lys406Glu)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322005.1']['p_hgvs_slc'] == 'NP_001308934.1:p.(K406E)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322005.1']['transcript_variant_error'] is None
        assert 'NM_000535.6' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_000535.6']['t_hgvs'] == 'NM_000535.6:c.1621A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_000535.6']['p_hgvs_tlc'] == 'NP_000526.2:p.(Lys541Glu)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_000535.6']['p_hgvs_slc'] == 'NP_000526.2:p.(K541E)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_000535.6']['transcript_variant_error'] is None
        assert 'NM_001322013.1' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322013.1']['t_hgvs'] == 'NM_001322013.1:c.1048A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322013.1']['p_hgvs_tlc'] == 'NP_001308942.1:p.(Lys350Glu)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322013.1']['p_hgvs_slc'] == 'NP_001308942.1:p.(K350E)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322013.1']['transcript_variant_error'] is None
        assert 'NM_001322012.1' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322012.1']['t_hgvs'] == 'NM_001322012.1:c.688A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322012.1']['p_hgvs_tlc'] == 'NP_001308941.1:p.(Lys230Glu)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322012.1']['p_hgvs_slc'] == 'NP_001308941.1:p.(K230E)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322012.1']['transcript_variant_error'] is None
        assert 'NM_001322006.1' in results['7-6026775-T-C']['hgvs_t_and_p'].keys()
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322006.1']['t_hgvs'] == 'NM_001322006.1:c.1465A>G'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322006.1']['p_hgvs_tlc'] == 'NP_001308935.1:p.(Lys489Glu)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322006.1']['p_hgvs_slc'] == 'NP_001308935.1:p.(K489E)'
        assert results['7-6026775-T-C']['hgvs_t_and_p']['NM_001322006.1']['transcript_variant_error'] is None

    def test_variant169(self):
        variant = '7-55242465-GGAATTAAGAGAAGCA-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '7-55242465-GGAATTAAGAGAAGCA-G' in results.keys()
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['p_vcf'] == '7-55242465-GGAATTAAGAGAAGCA-G'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['g_hgvs'] == 'NC_000007.13:g.55242466_55242480del'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['genomic_variant_error'] is None
        assert 'NM_005228.4' in results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p'].keys()
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_005228.4']['t_hgvs'] == 'NM_005228.4:c.2236_2250del'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_005228.4']['p_hgvs_tlc'] == 'NP_005219.2:p.(Glu746_Ala750del)'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_005228.4']['p_hgvs_slc'] == 'NP_005219.2:p.(E746_A750del)'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_005228.4']['transcript_variant_error'] is None
        assert 'NM_001346899.1' in results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p'].keys()
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346899.1']['t_hgvs'] == 'NM_001346899.1:c.2101_2115del'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346899.1']['p_hgvs_tlc'] == 'NP_001333828.1:p.(Glu701_Ala705del)'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346899.1']['p_hgvs_slc'] == 'NP_001333828.1:p.(E701_A705del)'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346899.1']['transcript_variant_error'] is None
        assert 'NM_001346897.1' in results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p'].keys()
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346897.1']['t_hgvs'] == 'NM_001346897.1:c.2101_2115del'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346897.1']['p_hgvs_tlc'] == 'NP_001333826.1:p.(Glu701_Ala705del)'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346897.1']['p_hgvs_slc'] == 'NP_001333826.1:p.(E701_A705del)'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346897.1']['transcript_variant_error'] is None
        assert 'NM_001346898.1' in results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p'].keys()
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346898.1']['t_hgvs'] == 'NM_001346898.1:c.2236_2250del'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346898.1']['p_hgvs_tlc'] == 'NP_001333827.1:p.(Glu746_Ala750del)'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346898.1']['p_hgvs_slc'] == 'NP_001333827.1:p.(E746_A750del)'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346898.1']['transcript_variant_error'] is None
        assert 'NM_005228.3' in results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p'].keys()
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_005228.3']['t_hgvs'] == 'NM_005228.3:c.2236_2250del'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_005228.3']['p_hgvs_tlc'] == 'NP_005219.2:p.(Glu746_Ala750del)'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_005228.3']['p_hgvs_slc'] == 'NP_005219.2:p.(E746_A750del)'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_005228.3']['transcript_variant_error'] is None
        assert 'NM_001346941.1' in results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p'].keys()
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346941.1']['t_hgvs'] == 'NM_001346941.1:c.1435_1449del'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346941.1']['p_hgvs_tlc'] == 'NP_001333870.1:p.(Glu479_Ala483del)'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346941.1']['p_hgvs_slc'] == 'NP_001333870.1:p.(E479_A483del)'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346941.1']['transcript_variant_error'] is None
        assert 'NM_001346900.1' in results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p'].keys()
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346900.1']['t_hgvs'] == 'NM_001346900.1:c.2077_2091del'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346900.1']['p_hgvs_tlc'] == 'NP_001333829.1:p.(Glu693_Ala697del)'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346900.1']['p_hgvs_slc'] == 'NP_001333829.1:p.(E693_A697del)'
        assert results['7-55242465-GGAATTAAGAGAAGCA-G']['hgvs_t_and_p']['NM_001346900.1']['transcript_variant_error'] is None

    def test_variant170(self):
        variant = '7-55248992-T-TTCCAGGAAGCCT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '7-55248992-T-TTCCAGGAAGCCT' in results.keys()
        assert results['7-55248992-T-TTCCAGGAAGCCT']['p_vcf'] == '7-55248980-C-CTCCAGGAAGCCT'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['g_hgvs'] == 'NC_000007.13:g.55248981_55248992dup'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['genomic_variant_error'] is None
        assert 'NM_005228.4' in results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p'].keys()
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_005228.4']['t_hgvs'] == 'NM_005228.4:c.2284-5_2290dup'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_005228.4']['p_hgvs_tlc'] == 'NP_005219.2:p.?'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_005228.4']['p_hgvs_slc'] == 'NP_005219.2:p.?'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_005228.4']['transcript_variant_error'] is None
        assert 'NM_001346899.1' in results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p'].keys()
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346899.1']['t_hgvs'] == 'NM_001346899.1:c.2149-5_2155dup'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346899.1']['p_hgvs_tlc'] == 'NP_001333828.1:p.?'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346899.1']['p_hgvs_slc'] == 'NP_001333828.1:p.?'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346899.1']['transcript_variant_error'] is None
        assert 'NM_001346897.1' in results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p'].keys()
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346897.1']['t_hgvs'] == 'NM_001346897.1:c.2149-5_2155dup'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346897.1']['p_hgvs_tlc'] == 'NP_001333826.1:p.?'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346897.1']['p_hgvs_slc'] == 'NP_001333826.1:p.?'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346897.1']['transcript_variant_error'] is None
        assert 'NM_001346898.1' in results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p'].keys()
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346898.1']['t_hgvs'] == 'NM_001346898.1:c.2284-5_2290dup'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346898.1']['p_hgvs_tlc'] == 'NP_001333827.1:p.?'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346898.1']['p_hgvs_slc'] == 'NP_001333827.1:p.?'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346898.1']['transcript_variant_error'] is None
        assert 'NM_005228.3' in results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p'].keys()
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_005228.3']['t_hgvs'] == 'NM_005228.3:c.2284-5_2290dup'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_005228.3']['p_hgvs_tlc'] == 'NP_005219.2:p.?'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_005228.3']['p_hgvs_slc'] == 'NP_005219.2:p.?'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_005228.3']['transcript_variant_error'] is None
        assert 'NM_001346941.1' in results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p'].keys()
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346941.1']['t_hgvs'] == 'NM_001346941.1:c.1483-5_1489dup'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346941.1']['p_hgvs_tlc'] == 'NP_001333870.1:p.?'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346941.1']['p_hgvs_slc'] == 'NP_001333870.1:p.?'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346941.1']['transcript_variant_error'] is None
        assert 'NM_001346900.1' in results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p'].keys()
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346900.1']['t_hgvs'] == 'NM_001346900.1:c.2125-5_2131dup'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346900.1']['p_hgvs_tlc'] == 'NP_001333829.1:p.?'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346900.1']['p_hgvs_slc'] == 'NP_001333829.1:p.?'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NM_001346900.1']['transcript_variant_error'] is None
        assert 'NR_047551.1' in results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p'].keys()
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NR_047551.1']['t_hgvs'] == 'NR_047551.1:n.1272_1283dup'
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NR_047551.1']['p_hgvs_tlc'] is None
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NR_047551.1']['p_hgvs_slc'] is None
        assert results['7-55248992-T-TTCCAGGAAGCCT']['hgvs_t_and_p']['NR_047551.1']['transcript_variant_error'] is None

    def test_variant171(self):
        variant = '7-75932111-C-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '7-75932111-C-A' in results.keys()
        assert results['7-75932111-C-A']['p_vcf'] == '7-75932111-C-A'
        assert results['7-75932111-C-A']['g_hgvs'] == 'NC_000007.13:g.75932111C>A'
        assert results['7-75932111-C-A']['genomic_variant_error'] is None
        assert 'NM_001540.4' in results['7-75932111-C-A']['hgvs_t_and_p'].keys()
        assert results['7-75932111-C-A']['hgvs_t_and_p']['NM_001540.4']['t_hgvs'] == 'NM_001540.4:c.82C>A'
        assert results['7-75932111-C-A']['hgvs_t_and_p']['NM_001540.4']['p_hgvs_tlc'] == 'NP_001531.1:p.(Leu28Ile)'
        assert results['7-75932111-C-A']['hgvs_t_and_p']['NM_001540.4']['p_hgvs_slc'] == 'NP_001531.1:p.(L28I)'
        assert results['7-75932111-C-A']['hgvs_t_and_p']['NM_001540.4']['transcript_variant_error'] is None
        assert 'NM_001540.3' in results['7-75932111-C-A']['hgvs_t_and_p'].keys()
        assert results['7-75932111-C-A']['hgvs_t_and_p']['NM_001540.3']['t_hgvs'] == 'NM_001540.3:c.82C>A'
        assert results['7-75932111-C-A']['hgvs_t_and_p']['NM_001540.3']['p_hgvs_tlc'] == 'NP_001531.1:p.(Leu28Ile)'
        assert results['7-75932111-C-A']['hgvs_t_and_p']['NM_001540.3']['p_hgvs_slc'] == 'NP_001531.1:p.(L28I)'
        assert results['7-75932111-C-A']['hgvs_t_and_p']['NM_001540.3']['transcript_variant_error'] is None

    def test_variant172(self):
        variant = '7-91652178-A-AAAC'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '7-91652178-A-AAAC' in results.keys()
        assert results['7-91652178-A-AAAC']['p_vcf'] == '7-91652178-A-AAAC'
        assert results['7-91652178-A-AAAC']['g_hgvs'] == 'NC_000007.13:g.91652179_91652181dup'
        assert results['7-91652178-A-AAAC']['genomic_variant_error'] is None
        assert 'NM_005751.4' in results['7-91652178-A-AAAC']['hgvs_t_and_p'].keys()
        assert results['7-91652178-A-AAAC']['hgvs_t_and_p']['NM_005751.4']['t_hgvs'] == 'NM_005751.4:c.4004_4006dup'
        assert results['7-91652178-A-AAAC']['hgvs_t_and_p']['NM_005751.4']['p_hgvs_tlc'] == 'NP_005742.4:p.(Lys1335_Leu1336insGln)'
        assert results['7-91652178-A-AAAC']['hgvs_t_and_p']['NM_005751.4']['p_hgvs_slc'] == 'NP_005742.4:p.(K1335_L1336insQ)'
        assert results['7-91652178-A-AAAC']['hgvs_t_and_p']['NM_005751.4']['transcript_variant_error'] is None
        assert 'NM_147185.2' in results['7-91652178-A-AAAC']['hgvs_t_and_p'].keys()
        assert results['7-91652178-A-AAAC']['hgvs_t_and_p']['NM_147185.2']['t_hgvs'] == 'NM_147185.2:c.4004_4006dup'
        assert results['7-91652178-A-AAAC']['hgvs_t_and_p']['NM_147185.2']['p_hgvs_tlc'] == 'NP_671714.1:p.(Lys1335_Leu1336insGln)'
        assert results['7-91652178-A-AAAC']['hgvs_t_and_p']['NM_147185.2']['p_hgvs_slc'] == 'NP_671714.1:p.(K1335_L1336insQ)'
        assert results['7-91652178-A-AAAC']['hgvs_t_and_p']['NM_147185.2']['transcript_variant_error'] is None

    def test_variant173(self):
        variant = '7-117199644-ATCT-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '7-117199644-ATCT-A' in results.keys()
        assert results['7-117199644-ATCT-A']['p_vcf'] == '7-117199644-ATCT-A'
        assert results['7-117199644-ATCT-A']['g_hgvs'] == 'NC_000007.13:g.117199646_117199648del'
        assert results['7-117199644-ATCT-A']['genomic_variant_error'] is None
        assert 'NR_149084.1' in results['7-117199644-ATCT-A']['hgvs_t_and_p'].keys()
        assert results['7-117199644-ATCT-A']['hgvs_t_and_p']['NR_149084.1']['t_hgvs'] == 'NR_149084.1:n.221+1140_221+1142del'
        assert results['7-117199644-ATCT-A']['hgvs_t_and_p']['NR_149084.1']['p_hgvs_tlc'] is None
        assert results['7-117199644-ATCT-A']['hgvs_t_and_p']['NR_149084.1']['p_hgvs_slc'] is None
        assert results['7-117199644-ATCT-A']['hgvs_t_and_p']['NR_149084.1']['transcript_variant_error'] is None
        assert 'NM_000492.3' in results['7-117199644-ATCT-A']['hgvs_t_and_p'].keys()
        assert results['7-117199644-ATCT-A']['hgvs_t_and_p']['NM_000492.3']['t_hgvs'] == 'NM_000492.3:c.1521_1523del'
        assert results['7-117199644-ATCT-A']['hgvs_t_and_p']['NM_000492.3']['p_hgvs_tlc'] == 'NP_000483.3:p.(Phe508del)'
        assert results['7-117199644-ATCT-A']['hgvs_t_and_p']['NM_000492.3']['p_hgvs_slc'] == 'NP_000483.3:p.(F508del)'
        assert results['7-117199644-ATCT-A']['hgvs_t_and_p']['NM_000492.3']['transcript_variant_error'] is None

    def test_variant174(self):
        variant = '7-140453136-AC-CT'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '7-140453136-AC-CT' in results.keys()
        assert results['7-140453136-AC-CT']['p_vcf'] == '7-140453136-AC-CT'
        assert results['7-140453136-AC-CT']['g_hgvs'] == 'NC_000007.13:g.140453136_140453137delinsCT'
        assert results['7-140453136-AC-CT']['genomic_variant_error'] is None
        assert 'NM_001354609.1' in results['7-140453136-AC-CT']['hgvs_t_and_p'].keys()
        assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NM_001354609.1']['t_hgvs'] == 'NM_001354609.1:c.1798_1799delinsAG'
        assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NM_001354609.1']['p_hgvs_tlc'] == 'NP_001341538.1:p.(Val600Arg)'
        assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NM_001354609.1']['p_hgvs_slc'] == 'NP_001341538.1:p.(V600R)'
        assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NM_001354609.1']['transcript_variant_error'] is None
        assert 'NM_004333.5' in results['7-140453136-AC-CT']['hgvs_t_and_p'].keys()
        assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NM_004333.5']['t_hgvs'] == 'NM_004333.5:c.1798_1799delinsAG'
        assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NM_004333.5']['p_hgvs_tlc'] == 'NP_004324.2:p.(Val600Arg)'
        assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NM_004333.5']['p_hgvs_slc'] == 'NP_004324.2:p.(V600R)'
        assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NM_004333.5']['transcript_variant_error'] is None
        assert 'NM_004333.4' in results['7-140453136-AC-CT']['hgvs_t_and_p'].keys()
        assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NM_004333.4']['t_hgvs'] == 'NM_004333.4:c.1798_1799delinsAG'
        assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NM_004333.4']['p_hgvs_tlc'] == 'NP_004324.2:p.(Val600Arg)'
        assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NM_004333.4']['p_hgvs_slc'] == 'NP_004324.2:p.(V600R)'
        assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NM_004333.4']['transcript_variant_error'] is None
        # Test removed. Currently Deprecated in VVTA
        # assert 'NR_148928.1' in results['7-140453136-AC-CT']['hgvs_t_and_p'].keys()
        # assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NR_148928.1']['t_hgvs'] == 'NR_148928.1:n.2896_2897delinsAG'
        # assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NR_148928.1']['p_hgvs_tlc'] is None
        # assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NR_148928.1']['p_hgvs_slc'] is None
        # assert results['7-140453136-AC-CT']['hgvs_t_and_p']['NR_148928.1']['transcript_variant_error'] is None

    def test_variant175(self):
        variant = '7-140453136-A-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '7-140453136-A-T' in results.keys()
        assert results['7-140453136-A-T']['p_vcf'] == '7-140453136-A-T'
        assert results['7-140453136-A-T']['g_hgvs'] == 'NC_000007.13:g.140453136A>T'
        assert results['7-140453136-A-T']['genomic_variant_error'] is None
        assert 'NM_001354609.1' in results['7-140453136-A-T']['hgvs_t_and_p'].keys()
        assert results['7-140453136-A-T']['hgvs_t_and_p']['NM_001354609.1']['t_hgvs'] == 'NM_001354609.1:c.1799T>A'
        assert results['7-140453136-A-T']['hgvs_t_and_p']['NM_001354609.1']['p_hgvs_tlc'] == 'NP_001341538.1:p.(Val600Glu)'
        assert results['7-140453136-A-T']['hgvs_t_and_p']['NM_001354609.1']['p_hgvs_slc'] == 'NP_001341538.1:p.(V600E)'
        assert results['7-140453136-A-T']['hgvs_t_and_p']['NM_001354609.1']['transcript_variant_error'] is None
        assert 'NM_004333.5' in results['7-140453136-A-T']['hgvs_t_and_p'].keys()
        assert results['7-140453136-A-T']['hgvs_t_and_p']['NM_004333.5']['t_hgvs'] == 'NM_004333.5:c.1799T>A'
        assert results['7-140453136-A-T']['hgvs_t_and_p']['NM_004333.5']['p_hgvs_tlc'] == 'NP_004324.2:p.(Val600Glu)'
        assert results['7-140453136-A-T']['hgvs_t_and_p']['NM_004333.5']['p_hgvs_slc'] == 'NP_004324.2:p.(V600E)'
        assert results['7-140453136-A-T']['hgvs_t_and_p']['NM_004333.5']['transcript_variant_error'] is None
        assert 'NM_004333.4' in results['7-140453136-A-T']['hgvs_t_and_p'].keys()
        assert results['7-140453136-A-T']['hgvs_t_and_p']['NM_004333.4']['t_hgvs'] == 'NM_004333.4:c.1799T>A'
        assert results['7-140453136-A-T']['hgvs_t_and_p']['NM_004333.4']['p_hgvs_tlc'] == 'NP_004324.2:p.(Val600Glu)'
        assert results['7-140453136-A-T']['hgvs_t_and_p']['NM_004333.4']['p_hgvs_slc'] == 'NP_004324.2:p.(V600E)'
        assert results['7-140453136-A-T']['hgvs_t_and_p']['NM_004333.4']['transcript_variant_error'] is None
        # Test removed. Currently Deprecated in VVTA
        # assert 'NR_148928.1' in results['7-140453136-A-T']['hgvs_t_and_p'].keys()
        # assert results['7-140453136-A-T']['hgvs_t_and_p']['NR_148928.1']['t_hgvs'] == 'NR_148928.1:n.2897T>A'
        # assert results['7-140453136-A-T']['hgvs_t_and_p']['NR_148928.1']['p_hgvs_tlc'] is None
        # assert results['7-140453136-A-T']['hgvs_t_and_p']['NR_148928.1']['p_hgvs_slc'] is None
        # assert results['7-140453136-A-T']['hgvs_t_and_p']['NR_148928.1']['transcript_variant_error'] is None

    def test_variant176(self):
        variant = '7-140453137-C-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '7-140453137-C-T' in results.keys()
        assert results['7-140453137-C-T']['p_vcf'] == '7-140453137-C-T'
        assert results['7-140453137-C-T']['g_hgvs'] == 'NC_000007.13:g.140453137C>T'
        assert results['7-140453137-C-T']['genomic_variant_error'] is None
        assert 'NM_001354609.1' in results['7-140453137-C-T']['hgvs_t_and_p'].keys()
        assert results['7-140453137-C-T']['hgvs_t_and_p']['NM_001354609.1']['t_hgvs'] == 'NM_001354609.1:c.1798G>A'
        assert results['7-140453137-C-T']['hgvs_t_and_p']['NM_001354609.1']['p_hgvs_tlc'] == 'NP_001341538.1:p.(Val600Met)'
        assert results['7-140453137-C-T']['hgvs_t_and_p']['NM_001354609.1']['p_hgvs_slc'] == 'NP_001341538.1:p.(V600M)'
        assert results['7-140453137-C-T']['hgvs_t_and_p']['NM_001354609.1']['transcript_variant_error'] is None
        assert 'NM_004333.5' in results['7-140453137-C-T']['hgvs_t_and_p'].keys()
        assert results['7-140453137-C-T']['hgvs_t_and_p']['NM_004333.5']['t_hgvs'] == 'NM_004333.5:c.1798G>A'
        assert results['7-140453137-C-T']['hgvs_t_and_p']['NM_004333.5']['p_hgvs_tlc'] == 'NP_004324.2:p.(Val600Met)'
        assert results['7-140453137-C-T']['hgvs_t_and_p']['NM_004333.5']['p_hgvs_slc'] == 'NP_004324.2:p.(V600M)'
        assert results['7-140453137-C-T']['hgvs_t_and_p']['NM_004333.5']['transcript_variant_error'] is None
        assert 'NM_004333.4' in results['7-140453137-C-T']['hgvs_t_and_p'].keys()
        assert results['7-140453137-C-T']['hgvs_t_and_p']['NM_004333.4']['t_hgvs'] == 'NM_004333.4:c.1798G>A'
        assert results['7-140453137-C-T']['hgvs_t_and_p']['NM_004333.4']['p_hgvs_tlc'] == 'NP_004324.2:p.(Val600Met)'
        assert results['7-140453137-C-T']['hgvs_t_and_p']['NM_004333.4']['p_hgvs_slc'] == 'NP_004324.2:p.(V600M)'
        assert results['7-140453137-C-T']['hgvs_t_and_p']['NM_004333.4']['transcript_variant_error'] is None
        # Test removed. Currently Deprecated in VVTA
        # assert 'NR_148928.1' in results['7-140453137-C-T']['hgvs_t_and_p'].keys()
        # assert results['7-140453137-C-T']['hgvs_t_and_p']['NR_148928.1']['t_hgvs'] == 'NR_148928.1:n.2896G>A'
        # assert results['7-140453137-C-T']['hgvs_t_and_p']['NR_148928.1']['p_hgvs_tlc'] is None
        # assert results['7-140453137-C-T']['hgvs_t_and_p']['NR_148928.1']['p_hgvs_slc'] is None
        # assert results['7-140453137-C-T']['hgvs_t_and_p']['NR_148928.1']['transcript_variant_error'] is None

    def test_variant177(self):
        variant = '7-143013488-A-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '7-143013488-A-T' in results.keys()
        assert results['7-143013488-A-T']['p_vcf'] == '7-143013488-A-T'
        assert results['7-143013488-A-T']['g_hgvs'] == 'NC_000007.13:g.143013488A>T'
        assert results['7-143013488-A-T']['genomic_variant_error'] is None
        assert 'NR_046453.1' in results['7-143013488-A-T']['hgvs_t_and_p'].keys()
        assert results['7-143013488-A-T']['hgvs_t_and_p']['NR_046453.1']['t_hgvs'] == 'NR_046453.1:n.267+3A>T'
        assert results['7-143013488-A-T']['hgvs_t_and_p']['NR_046453.1']['p_hgvs_tlc'] is None
        assert results['7-143013488-A-T']['hgvs_t_and_p']['NR_046453.1']['p_hgvs_slc'] is None
        assert results['7-143013488-A-T']['hgvs_t_and_p']['NR_046453.1']['transcript_variant_error'] is None
        assert 'NM_000083.2' in results['7-143013488-A-T']['hgvs_t_and_p'].keys()
        assert results['7-143013488-A-T']['hgvs_t_and_p']['NM_000083.2']['t_hgvs'] == 'NM_000083.2:c.180+3A>T'
        assert results['7-143013488-A-T']['hgvs_t_and_p']['NM_000083.2']['p_hgvs_tlc'] == 'NP_000074.2:p.?'
        assert results['7-143013488-A-T']['hgvs_t_and_p']['NM_000083.2']['p_hgvs_slc'] == 'NP_000074.2:p.?'
        assert results['7-143013488-A-T']['hgvs_t_and_p']['NM_000083.2']['transcript_variant_error'] is None

    def test_variant178(self):
        variant = '7-143018934-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '7-143018934-G-A' in results.keys()
        assert results['7-143018934-G-A']['p_vcf'] == '7-143018934-G-A'
        assert results['7-143018934-G-A']['g_hgvs'] == 'NC_000007.13:g.143018934G>A'
        assert results['7-143018934-G-A']['genomic_variant_error'] is None
        assert 'NR_046453.1' in results['7-143018934-G-A']['hgvs_t_and_p'].keys()
        assert results['7-143018934-G-A']['hgvs_t_and_p']['NR_046453.1']['t_hgvs'] == 'NR_046453.1:n.776G>A'
        assert results['7-143018934-G-A']['hgvs_t_and_p']['NR_046453.1']['p_hgvs_tlc'] is None
        assert results['7-143018934-G-A']['hgvs_t_and_p']['NR_046453.1']['p_hgvs_slc'] is None
        assert results['7-143018934-G-A']['hgvs_t_and_p']['NR_046453.1']['transcript_variant_error'] is None
        assert 'NM_000083.2' in results['7-143018934-G-A']['hgvs_t_and_p'].keys()
        assert results['7-143018934-G-A']['hgvs_t_and_p']['NM_000083.2']['t_hgvs'] == 'NM_000083.2:c.689G>A'
        assert results['7-143018934-G-A']['hgvs_t_and_p']['NM_000083.2']['p_hgvs_tlc'] == 'NP_000074.2:p.(Gly230Glu)'
        assert results['7-143018934-G-A']['hgvs_t_and_p']['NM_000083.2']['p_hgvs_slc'] == 'NP_000074.2:p.(G230E)'
        assert results['7-143018934-G-A']['hgvs_t_and_p']['NM_000083.2']['transcript_variant_error'] is None

    def test_variant179(self):
        variant = '7-143048771-C-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '7-143048771-C-T' in results.keys()
        assert results['7-143048771-C-T']['p_vcf'] == '7-143048771-C-T'
        assert results['7-143048771-C-T']['g_hgvs'] == 'NC_000007.13:g.143048771C>T'
        assert results['7-143048771-C-T']['genomic_variant_error'] is None
        assert 'NR_046453.1' in results['7-143048771-C-T']['hgvs_t_and_p'].keys()
        assert results['7-143048771-C-T']['hgvs_t_and_p']['NR_046453.1']['t_hgvs'] == 'NR_046453.1:n.2620C>T'
        assert results['7-143048771-C-T']['hgvs_t_and_p']['NR_046453.1']['p_hgvs_tlc'] is None
        assert results['7-143048771-C-T']['hgvs_t_and_p']['NR_046453.1']['p_hgvs_slc'] is None
        assert results['7-143048771-C-T']['hgvs_t_and_p']['NR_046453.1']['transcript_variant_error'] is None
        assert 'NM_000083.2' in results['7-143048771-C-T']['hgvs_t_and_p'].keys()
        assert results['7-143048771-C-T']['hgvs_t_and_p']['NM_000083.2']['t_hgvs'] == 'NM_000083.2:c.2680C>T'
        assert results['7-143048771-C-T']['hgvs_t_and_p']['NM_000083.2']['p_hgvs_tlc'] == 'NP_000074.2:p.(Arg894Ter)'
        assert results['7-143048771-C-T']['hgvs_t_and_p']['NM_000083.2']['p_hgvs_slc'] == 'NP_000074.2:p.(R894*)'
        assert results['7-143048771-C-T']['hgvs_t_and_p']['NM_000083.2']['transcript_variant_error'] is None

    def test_variant180(self):
        variant = '8-1871951-C-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '8-1871951-C-T' in results.keys()
        assert results['8-1871951-C-T']['p_vcf'] == '8-1871951-C-T'
        assert results['8-1871951-C-T']['g_hgvs'] == 'NC_000008.10:g.1871951C>T'
        assert results['8-1871951-C-T']['genomic_variant_error'] is None
        assert 'NM_014629.3' in results['8-1871951-C-T']['hgvs_t_and_p'].keys()
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_014629.3']['t_hgvs'] == 'NM_014629.3:c.2399C>T'
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_014629.3']['p_hgvs_tlc'] == 'NP_055444.2:p.(Pro800Leu)'
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_014629.3']['p_hgvs_slc'] == 'NP_055444.2:p.(P800L)'
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_014629.3']['transcript_variant_error'] is None
        assert 'NM_014629.2' in results['8-1871951-C-T']['hgvs_t_and_p'].keys()
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_014629.2']['t_hgvs'] == 'NM_014629.2:c.2399C>T'
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_014629.2']['p_hgvs_tlc'] == 'NP_055444.2:p.(Pro800Leu)'
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_014629.2']['p_hgvs_slc'] == 'NP_055444.2:p.(P800L)'
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_014629.2']['transcript_variant_error'] is None
        assert 'NM_001308152.1' in results['8-1871951-C-T']['hgvs_t_and_p'].keys()
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_001308152.1']['t_hgvs'] == 'NM_001308152.1:c.2285C>T'
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_001308152.1']['p_hgvs_tlc'] == 'NP_001295081.1:p.(Pro762Leu)'
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_001308152.1']['p_hgvs_slc'] == 'NP_001295081.1:p.(P762L)'
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_001308152.1']['transcript_variant_error'] is None
        assert 'NM_001308153.1' in results['8-1871951-C-T']['hgvs_t_and_p'].keys()
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_001308153.1']['t_hgvs'] == 'NM_001308153.1:c.2471C>T'
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_001308153.1']['p_hgvs_tlc'] == 'NP_001295082.1:p.(Pro824Leu)'
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_001308153.1']['p_hgvs_slc'] == 'NP_001295082.1:p.(P824L)'
        assert results['8-1871951-C-T']['hgvs_t_and_p']['NM_001308153.1']['transcript_variant_error'] is None

    def test_variant181(self):
        variant = '9-13112056-T-TG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '9-13112056-T-TG' in results.keys()
        assert results['9-13112056-T-TG']['p_vcf'] == '9-13112056-T-TG'
        assert results['9-13112056-T-TG']['g_hgvs'] == 'NC_000009.11:g.13112059dup'
        assert results['9-13112056-T-TG']['genomic_variant_error'] is None
        assert 'NM_003829.4' in results['9-13112056-T-TG']['hgvs_t_and_p'].keys()
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_003829.4']['t_hgvs'] == 'NM_003829.4:c.5603dup'
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_003829.4']['p_hgvs_tlc'] == 'NP_003820.2:p.(Thr1869AsnfsTer15)'
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_003829.4']['p_hgvs_slc'] == 'NP_003820.2:p.(T1869Nfs*15)'
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_003829.4']['transcript_variant_error'] is None
        assert 'NM_001261407.1' in results['9-13112056-T-TG']['hgvs_t_and_p'].keys()
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_001261407.1']['t_hgvs'] == 'NM_001261407.1:c.5504dup'
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_001261407.1']['p_hgvs_tlc'] == 'NP_001248336.1:p.(Thr1836AsnfsTer15)'
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_001261407.1']['p_hgvs_slc'] == 'NP_001248336.1:p.(T1836Nfs*15)'
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_001261407.1']['transcript_variant_error'] is None
        assert 'NM_001330637.1' in results['9-13112056-T-TG']['hgvs_t_and_p'].keys()
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_001330637.1']['t_hgvs'] == 'NM_001330637.1:c.5690dup'
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_001330637.1']['p_hgvs_tlc'] == 'NP_001317566.1:p.(Thr1898AsnfsTer15)'
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_001330637.1']['p_hgvs_slc'] == 'NP_001317566.1:p.(T1898Nfs*15)'
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_001330637.1']['transcript_variant_error'] is None
        assert 'NM_001261406.1' in results['9-13112056-T-TG']['hgvs_t_and_p'].keys()
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_001261406.1']['t_hgvs'] == 'NM_001261406.1:c.5591dup'
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_001261406.1']['p_hgvs_tlc'] == 'NP_001248335.1:p.(Thr1865AsnfsTer15)'
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_001261406.1']['p_hgvs_slc'] == 'NP_001248335.1:p.(T1865Nfs*15)'
        assert results['9-13112056-T-TG']['hgvs_t_and_p']['NM_001261406.1']['transcript_variant_error'] is None

    def test_variant182(self):
        variant = '9-21971208-C-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '9-21971208-C-A' in results.keys()
        assert results['9-21971208-C-A']['p_vcf'] == '9-21971208-C-A'
        assert results['9-21971208-C-A']['g_hgvs'] == 'NC_000009.11:g.21971208C>A'
        assert results['9-21971208-C-A']['genomic_variant_error'] is None
        assert 'NM_000077.4' in results['9-21971208-C-A']['hgvs_t_and_p'].keys()
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_000077.4']['t_hgvs'] == 'NM_000077.4:c.151-1G>T'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_000077.4']['p_hgvs_tlc'] == 'NP_000068.1:p.?'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_000077.4']['p_hgvs_slc'] == 'NP_000068.1:p.?'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_000077.4']['transcript_variant_error'] is None
        assert 'NM_001363763.1' in results['9-21971208-C-A']['hgvs_t_and_p'].keys()
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_001363763.1']['t_hgvs'] == 'NM_001363763.1:c.-3-1G>T'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_001363763.1']['p_hgvs_tlc'] == 'NP_001350692.1:p.?'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_001363763.1']['p_hgvs_slc'] == 'NP_001350692.1:p.?'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_001363763.1']['transcript_variant_error'] is None
        assert 'NM_058197.4' in results['9-21971208-C-A']['hgvs_t_and_p'].keys()
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_058197.4']['t_hgvs'] == 'NM_058197.4:c.*74-1G>T'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_058197.4']['p_hgvs_tlc'] == 'NP_478104.2:p.?'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_058197.4']['p_hgvs_slc'] == 'NP_478104.2:p.?'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_058197.4']['transcript_variant_error'] is None
        assert 'NM_001195132.1' in results['9-21971208-C-A']['hgvs_t_and_p'].keys()
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_001195132.1']['t_hgvs'] == 'NM_001195132.1:c.151-1G>T'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_001195132.1']['p_hgvs_tlc'] == 'NP_001182061.1:p.?'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_001195132.1']['p_hgvs_slc'] == 'NP_001182061.1:p.?'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_001195132.1']['transcript_variant_error'] is None
        assert 'NM_058195.3' in results['9-21971208-C-A']['hgvs_t_and_p'].keys()
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_058195.3']['t_hgvs'] == 'NM_058195.3:c.194-1G>T'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_058195.3']['p_hgvs_tlc'] == 'NP_478102.2:p.?'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_058195.3']['p_hgvs_slc'] == 'NP_478102.2:p.?'
        assert results['9-21971208-C-A']['hgvs_t_and_p']['NM_058195.3']['transcript_variant_error'] is None

    def test_variant183(self):
        variant = '9-35683240-T-TG'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '9-35683240-T-TG' in results.keys()
        assert results['9-35683240-T-TG']['p_vcf'] == '9-35683240-T-TG'
        assert results['9-35683240-T-TG']['g_hgvs'] == 'NC_000009.11:g.35683248dup'
        assert results['9-35683240-T-TG']['genomic_variant_error'] is None
        assert 'NM_001301226.1' in results['9-35683240-T-TG']['hgvs_t_and_p'].keys()
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_001301226.1']['t_hgvs'] == 'NM_001301226.1:c.772+1002dup'
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_001301226.1']['p_hgvs_tlc'] == 'NP_001288155.1:p.?'
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_001301226.1']['p_hgvs_slc'] == 'NP_001288155.1:p.?'
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_001301226.1']['transcript_variant_error'] is None
        assert 'NM_213674.1' in results['9-35683240-T-TG']['hgvs_t_and_p'].keys()
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_213674.1']['t_hgvs'] == 'NM_213674.1:c.772+1002dup'
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_213674.1']['p_hgvs_tlc'] == 'NP_998839.1:p.?'
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_213674.1']['p_hgvs_slc'] == 'NP_998839.1:p.?'
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_213674.1']['transcript_variant_error'] is None
        assert 'NM_003289.3' in results['9-35683240-T-TG']['hgvs_t_and_p'].keys()
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_003289.3']['t_hgvs'] == 'NM_003289.3:c.773-3dup'
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_003289.3']['p_hgvs_tlc'] == 'NP_003280.2:p.?'
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_003289.3']['p_hgvs_slc'] == 'NP_003280.2:p.?'
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_003289.3']['transcript_variant_error'] is None
        assert 'NM_001301227.1' in results['9-35683240-T-TG']['hgvs_t_and_p'].keys()
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_001301227.1']['t_hgvs'] == 'NM_001301227.1:c.773-3dup'
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_001301227.1']['p_hgvs_tlc'] == 'NP_001288156.1:p.?'
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_001301227.1']['p_hgvs_slc'] == 'NP_001288156.1:p.?'
        assert results['9-35683240-T-TG']['hgvs_t_and_p']['NM_001301227.1']['transcript_variant_error'] is None

    def test_variant184(self):
        variant = '9-135796754-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '9-135796754-G-A' in results.keys()
        assert results['9-135796754-G-A']['p_vcf'] == '9-135796754-G-A'
        assert results['9-135796754-G-A']['g_hgvs'] == 'NC_000009.11:g.135796754G>A'
        assert results['9-135796754-G-A']['genomic_variant_error'] is None
        assert 'NM_001162427.1' in results['9-135796754-G-A']['hgvs_t_and_p'].keys()
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_001162427.1']['t_hgvs'] == 'NM_001162427.1:c.580C>T'
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_001162427.1']['p_hgvs_tlc'] == 'NP_001155899.1:p.(Arg194Ter)'
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_001162427.1']['p_hgvs_slc'] == 'NP_001155899.1:p.(R194*)'
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_001162427.1']['transcript_variant_error'] is None
        assert 'NM_001362177.1' in results['9-135796754-G-A']['hgvs_t_and_p'].keys()
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_001362177.1']['t_hgvs'] == 'NM_001362177.1:c.370C>T'
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_001362177.1']['p_hgvs_tlc'] == 'NP_001349106.1:p.(Arg124Ter)'
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_001362177.1']['p_hgvs_slc'] == 'NP_001349106.1:p.(R124*)'
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_001362177.1']['transcript_variant_error'] is None
        assert 'NM_001162426.1' in results['9-135796754-G-A']['hgvs_t_and_p'].keys()
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_001162426.1']['t_hgvs'] == 'NM_001162426.1:c.733C>T'
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_001162426.1']['p_hgvs_tlc'] == 'NP_001155898.1:p.(Arg245Ter)'
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_001162426.1']['p_hgvs_slc'] == 'NP_001155898.1:p.(R245*)'
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_001162426.1']['transcript_variant_error'] is None
        assert 'NM_000368.4' in results['9-135796754-G-A']['hgvs_t_and_p'].keys()
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_000368.4']['t_hgvs'] == 'NM_000368.4:c.733C>T'
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_000368.4']['p_hgvs_tlc'] == 'NP_000359.1:p.(Arg245Ter)'
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_000368.4']['p_hgvs_slc'] == 'NP_000359.1:p.(R245*)'
        assert results['9-135796754-G-A']['hgvs_t_and_p']['NM_000368.4']['transcript_variant_error'] is None

    def test_variant185(self):
        variant = 'HG536_PATCH-10391-AC-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'HG536_PATCH-10391-AC-A' in results.keys()
        assert results['HG536_PATCH-10391-AC-A']['p_vcf'] == 'HG536_PATCH-10391-AC-A'
        assert results['HG536_PATCH-10391-AC-A']['g_hgvs'] == 'NW_003571046.1:g.10396del'
        assert results['HG536_PATCH-10391-AC-A']['genomic_variant_error'] is None
        assert 'NM_005247.2' in results['HG536_PATCH-10391-AC-A']['hgvs_t_and_p'].keys()
        assert results['HG536_PATCH-10391-AC-A']['hgvs_t_and_p']['NM_005247.2']['t_hgvs'] == 'NM_005247.2:c.616del'
        assert results['HG536_PATCH-10391-AC-A']['hgvs_t_and_p']['NM_005247.2']['p_hgvs_tlc'] == 'NP_005238.1:p.(Val206SerfsTer117)'
        assert results['HG536_PATCH-10391-AC-A']['hgvs_t_and_p']['NM_005247.2']['p_hgvs_slc'] == 'NP_005238.1:p.(V206Sfs*117)'
        assert results['HG536_PATCH-10391-AC-A']['hgvs_t_and_p']['NM_005247.2']['transcript_variant_error'] is None

    def test_variant186(self):
        variant = 'HG865_PATCH-33547-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'HG865_PATCH-33547-G-A' in results.keys()
        assert results['HG865_PATCH-33547-G-A']['p_vcf'] == 'HG865_PATCH-33547-G-A'
        assert results['HG865_PATCH-33547-G-A']['g_hgvs'] == 'NW_004070871.1:g.33547G>A'
        assert results['HG865_PATCH-33547-G-A']['genomic_variant_error'] is None
        assert 'NM_133266.4' in results['HG865_PATCH-33547-G-A']['hgvs_t_and_p'].keys()
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NM_133266.4']['t_hgvs'] == 'NM_133266.4:c.802C>T'
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NM_133266.4']['p_hgvs_tlc'] == 'NP_573573.2:p.(Leu268=)'
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NM_133266.4']['transcript_variant_error'] is None
        assert 'NM_012309.4' in results['HG865_PATCH-33547-G-A']['hgvs_t_and_p'].keys()
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NM_012309.4']['t_hgvs'] == 'NM_012309.4:c.2566C>T'
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NM_012309.4']['p_hgvs_tlc'] == 'NP_036441.2:p.(Leu856=)'
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NM_012309.4']['p_hgvs_slc'] == 'NP_036441.2:p.(L856=)'
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NM_012309.4']['transcript_variant_error'] is None
        assert 'NR_110766.1' in results['HG865_PATCH-33547-G-A']['hgvs_t_and_p'].keys()
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NR_110766.1']['t_hgvs'] == 'NR_110766.1:n.833+969C>T'
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NR_110766.1']['p_hgvs_tlc'] is None
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NR_110766.1']['p_hgvs_slc'] is None
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NR_110766.1']['transcript_variant_error'] is None
        assert 'NM_133266.3' in results['HG865_PATCH-33547-G-A']['hgvs_t_and_p'].keys()
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NM_133266.3']['t_hgvs'] == 'NM_133266.3:c.802C>T'
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NM_133266.3']['p_hgvs_tlc'] == 'NP_573573.2:p.(Leu268=)'
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NM_133266.3']['p_hgvs_slc'] == 'NP_573573.2:p.(L268=)'
        assert results['HG865_PATCH-33547-G-A']['hgvs_t_and_p']['NM_133266.3']['transcript_variant_error'] is None

    def test_variant187(self):
        variant = 'HG865_PATCH-569441-G-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'HG865_PATCH-569441-G-T' in results.keys()
        assert results['HG865_PATCH-569441-G-T']['p_vcf'] == 'HG865_PATCH-569441-G-T'
        assert results['HG865_PATCH-569441-G-T']['g_hgvs'] == 'NW_004070871.1:g.569441G>T'
        assert results['HG865_PATCH-569441-G-T']['genomic_variant_error'] is None
        assert 'NM_012309.4' in results['HG865_PATCH-569441-G-T']['hgvs_t_and_p'].keys()
        assert results['HG865_PATCH-569441-G-T']['hgvs_t_and_p']['NM_012309.4']['t_hgvs'] == 'NM_012309.4:c.960C>A'
        assert results['HG865_PATCH-569441-G-T']['hgvs_t_and_p']['NM_012309.4']['p_hgvs_tlc'] == 'NP_036441.2:p.(Tyr320Ter)'
        assert results['HG865_PATCH-569441-G-T']['hgvs_t_and_p']['NM_012309.4']['p_hgvs_slc'] == 'NP_036441.2:p.(Y320*)'
        assert results['HG865_PATCH-569441-G-T']['hgvs_t_and_p']['NM_012309.4']['transcript_variant_error'] is None
        # Test removed. Currently Deprecated in VVTA
        # assert 'NM_012309.3' in results['HG865_PATCH-569441-G-T']['hgvs_t_and_p'].keys()
        # assert results['HG865_PATCH-569441-G-T']['hgvs_t_and_p']['NM_012309.3']['t_hgvs'] is None
        # assert results['HG865_PATCH-569441-G-T']['hgvs_t_and_p']['NM_012309.3']['p_hgvs_tlc'] is None
        # assert results['HG865_PATCH-569441-G-T']['hgvs_t_and_p']['NM_012309.3']['p_hgvs_slc'] is None
        # assert results['HG865_PATCH-569441-G-T']['hgvs_t_and_p']['NM_012309.3']['transcript_variant_error'] == 'start or end or both are beyond the bounds of transcript record'

    def test_variant188(self):
        variant = 'HG865_PATCH-574546-C-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'HG865_PATCH-574546-C-T' in results.keys()
        assert results['HG865_PATCH-574546-C-T']['p_vcf'] == 'HG865_PATCH-574546-C-T'
        assert results['HG865_PATCH-574546-C-T']['g_hgvs'] == 'NW_004070871.1:g.574546C>T'
        assert results['HG865_PATCH-574546-C-T']['genomic_variant_error'] is None
        assert 'NM_012309.4' in results['HG865_PATCH-574546-C-T']['hgvs_t_and_p'].keys()
        assert results['HG865_PATCH-574546-C-T']['hgvs_t_and_p']['NM_012309.4']['t_hgvs'] == 'NM_012309.4:c.913-5058G>A'
        assert results['HG865_PATCH-574546-C-T']['hgvs_t_and_p']['NM_012309.4']['p_hgvs_tlc'] == 'NP_036441.2:p.?'
        assert results['HG865_PATCH-574546-C-T']['hgvs_t_and_p']['NM_012309.4']['p_hgvs_slc'] == 'NP_036441.2:p.?'
        assert results['HG865_PATCH-574546-C-T']['hgvs_t_and_p']['NM_012309.4']['transcript_variant_error'] is None
        # Test removed. Currently Deprecated in VVTA
        # assert 'NM_012309.3' in results['HG865_PATCH-574546-C-T']['hgvs_t_and_p'].keys()
        # assert results['HG865_PATCH-574546-C-T']['hgvs_t_and_p']['NM_012309.3']['t_hgvs'] is None
        # assert results['HG865_PATCH-574546-C-T']['hgvs_t_and_p']['NM_012309.3']['p_hgvs_tlc'] is None
        # assert results['HG865_PATCH-574546-C-T']['hgvs_t_and_p']['NM_012309.3']['p_hgvs_slc'] is None
        # assert results['HG865_PATCH-574546-C-T']['hgvs_t_and_p']['NM_012309.3']['transcript_variant_error'] == 'start or end or both are beyond the bounds of transcript record'

    def test_variant189(self):
        variant = 'HSCHR1_1_CTG31-133178-TAG-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'HSCHR1_1_CTG31-133178-TAG-T' in results.keys()
        assert results['HSCHR1_1_CTG31-133178-TAG-T']['p_vcf'] == 'HSCHR1_1_CTG31-133178-TAG-T'
        assert results['HSCHR1_1_CTG31-133178-TAG-T']['g_hgvs'] == 'NW_003315905.1:g.133179_133180del'
        assert results['HSCHR1_1_CTG31-133178-TAG-T']['genomic_variant_error'] is None
        assert 'NM_020699.4' in results['HSCHR1_1_CTG31-133178-TAG-T']['hgvs_t_and_p'].keys()
        assert results['HSCHR1_1_CTG31-133178-TAG-T']['hgvs_t_and_p']['NM_020699.4']['t_hgvs'] == 'NM_020699.4:c.774_775del'
        assert results['HSCHR1_1_CTG31-133178-TAG-T']['hgvs_t_and_p']['NM_020699.4']['p_hgvs_tlc'] == 'NP_065750.1:p.(Met259ValfsTer22)'
        assert results['HSCHR1_1_CTG31-133178-TAG-T']['hgvs_t_and_p']['NM_020699.4']['p_hgvs_slc'] == 'NP_065750.1:p.(M259Vfs*22)'
        assert results['HSCHR1_1_CTG31-133178-TAG-T']['hgvs_t_and_p']['NM_020699.4']['transcript_variant_error'] is None

    def test_variant190(self):
        variant = 'HSCHR6_MHC_MANN_CTG1-3848158-T-G'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'HSCHR6_MHC_MANN_CTG1-3848158-T-G' in results.keys()
        assert results['HSCHR6_MHC_MANN_CTG1-3848158-T-G']['p_vcf'] == 'HSCHR6_MHC_MANN_CTG1-3848158-T-G'
        assert results['HSCHR6_MHC_MANN_CTG1-3848158-T-G']['g_hgvs'] == 'NT_167246.1:g.3848158T>G'
        assert results['HSCHR6_MHC_MANN_CTG1-3848158-T-G']['genomic_variant_error'] is None
        assert 'NM_021983.4' in results['HSCHR6_MHC_MANN_CTG1-3848158-T-G']['hgvs_t_and_p'].keys()
        assert results['HSCHR6_MHC_MANN_CTG1-3848158-T-G']['hgvs_t_and_p']['NM_021983.4']['t_hgvs'] == 'NM_021983.4:c.490G>C'
        assert results['HSCHR6_MHC_MANN_CTG1-3848158-T-G']['hgvs_t_and_p']['NM_021983.4']['p_hgvs_tlc'] == 'NP_068818.4:p.(Gly164Arg)'
        assert results['HSCHR6_MHC_MANN_CTG1-3848158-T-G']['hgvs_t_and_p']['NM_021983.4']['p_hgvs_slc'] == 'NP_068818.4:p.(G164R)'
        assert results['HSCHR6_MHC_MANN_CTG1-3848158-T-G']['hgvs_t_and_p']['NM_021983.4']['transcript_variant_error'] is None

    def test_variant191(self):
        variant = 'HSCHR6_MHC_MANN_CTG1-3851043-C-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'HSCHR6_MHC_MANN_CTG1-3851043-C-A' in results.keys()
        assert results['HSCHR6_MHC_MANN_CTG1-3851043-C-A']['p_vcf'] == 'HSCHR6_MHC_MANN_CTG1-3851043-C-A'
        assert results['HSCHR6_MHC_MANN_CTG1-3851043-C-A']['g_hgvs'] == 'NT_167246.1:g.3851043C>A'
        assert results['HSCHR6_MHC_MANN_CTG1-3851043-C-A']['genomic_variant_error'] is None
        assert 'NM_021983.4' in results['HSCHR6_MHC_MANN_CTG1-3851043-C-A']['hgvs_t_and_p'].keys()
        assert results['HSCHR6_MHC_MANN_CTG1-3851043-C-A']['hgvs_t_and_p']['NM_021983.4']['t_hgvs'] == 'NM_021983.4:c.346G>T'
        assert results['HSCHR6_MHC_MANN_CTG1-3851043-C-A']['hgvs_t_and_p']['NM_021983.4']['p_hgvs_tlc'] == 'NP_068818.4:p.(Glu116Ter)'
        assert results['HSCHR6_MHC_MANN_CTG1-3851043-C-A']['hgvs_t_and_p']['NM_021983.4']['p_hgvs_slc'] == 'NP_068818.4:p.(E116*)'
        assert results['HSCHR6_MHC_MANN_CTG1-3851043-C-A']['hgvs_t_and_p']['NM_021983.4']['transcript_variant_error'] is None

    def test_variant192(self):
        variant = 'X-70443101-C-T'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'X-70443101-C-T' in results.keys()
        assert results['X-70443101-C-T']['p_vcf'] == 'X-70443101-C-T'
        assert results['X-70443101-C-T']['g_hgvs'] == 'NC_000023.10:g.70443101C>T'
        assert results['X-70443101-C-T']['genomic_variant_error'] is None
        assert 'NM_000166.5' in results['X-70443101-C-T']['hgvs_t_and_p'].keys()
        assert results['X-70443101-C-T']['hgvs_t_and_p']['NM_000166.5']['t_hgvs'] == 'NM_000166.5:c.-101C>T'
        assert results['X-70443101-C-T']['hgvs_t_and_p']['NM_000166.5']['p_hgvs_tlc'] == 'NP_000157.1:p.?'
        assert results['X-70443101-C-T']['hgvs_t_and_p']['NM_000166.5']['p_hgvs_slc'] == 'NP_000157.1:p.?'
        assert results['X-70443101-C-T']['hgvs_t_and_p']['NM_000166.5']['transcript_variant_error'] is None
        assert 'NM_001097642.2' in results['X-70443101-C-T']['hgvs_t_and_p'].keys()
        assert results['X-70443101-C-T']['hgvs_t_and_p']['NM_001097642.2']['t_hgvs'] == 'NM_001097642.2:c.-16-441C>T'
        assert results['X-70443101-C-T']['hgvs_t_and_p']['NM_001097642.2']['p_hgvs_tlc'] == 'NP_001091111.1:p.?'
        assert results['X-70443101-C-T']['hgvs_t_and_p']['NM_001097642.2']['p_hgvs_slc'] == 'NP_001091111.1:p.?'
        assert results['X-70443101-C-T']['hgvs_t_and_p']['NM_001097642.2']['transcript_variant_error'] is None

    def test_variant193(self):
        variant = 'X-153296777-G-A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'X-153296777-G-A' in results.keys()
        assert results['X-153296777-G-A']['p_vcf'] == 'X-153296777-G-A'
        assert results['X-153296777-G-A']['g_hgvs'] == 'NC_000023.10:g.153296777G>A'
        assert results['X-153296777-G-A']['genomic_variant_error'] is None
        assert 'NM_001316337.1' in results['X-153296777-G-A']['hgvs_t_and_p'].keys()
        assert results['X-153296777-G-A']['hgvs_t_and_p']['NM_001316337.1']['t_hgvs'] == 'NM_001316337.1:c.223C>T'
        assert results['X-153296777-G-A']['hgvs_t_and_p']['NM_001316337.1']['p_hgvs_tlc'] == 'NP_001303266.1:p.(Arg75Ter)'
        assert results['X-153296777-G-A']['hgvs_t_and_p']['NM_001316337.1']['p_hgvs_slc'] == 'NP_001303266.1:p.(R75*)'
        assert results['X-153296777-G-A']['hgvs_t_and_p']['NM_001316337.1']['transcript_variant_error'] is None
        assert 'NM_004992.3' in results['X-153296777-G-A']['hgvs_t_and_p'].keys()
        assert results['X-153296777-G-A']['hgvs_t_and_p']['NM_004992.3']['t_hgvs'] == 'NM_004992.3:c.502C>T'
        assert results['X-153296777-G-A']['hgvs_t_and_p']['NM_004992.3']['p_hgvs_tlc'] == 'NP_004983.1:p.(Arg168Ter)'
        assert results['X-153296777-G-A']['hgvs_t_and_p']['NM_004992.3']['p_hgvs_slc'] == 'NP_004983.1:p.(R168*)'
        assert results['X-153296777-G-A']['hgvs_t_and_p']['NM_004992.3']['transcript_variant_error'] is None
        assert 'NM_001110792.1' in results['X-153296777-G-A']['hgvs_t_and_p'].keys()
        assert results['X-153296777-G-A']['hgvs_t_and_p']['NM_001110792.1']['t_hgvs'] == 'NM_001110792.1:c.538C>T'
        assert results['X-153296777-G-A']['hgvs_t_and_p']['NM_001110792.1']['p_hgvs_tlc'] == 'NP_001104262.1:p.(Arg180Ter)'
        assert results['X-153296777-G-A']['hgvs_t_and_p']['NM_001110792.1']['p_hgvs_slc'] == 'NP_001104262.1:p.(R180*)'
        assert results['X-153296777-G-A']['hgvs_t_and_p']['NM_001110792.1']['transcript_variant_error'] is None

    def test_variant194(self):
        variant = 'NC_000023.11:g.33215693T>G'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.33215693T>G' in results.keys()
        assert results['NC_000023.11:g.33215693T>G']['p_vcf'] == 'X:33215693:T:G'
        assert results['NC_000023.11:g.33215693T>G']['g_hgvs'] == 'NC_000023.11:g.33215693T>G'
        assert results['NC_000023.11:g.33215693T>G']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.33215693T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.33215693T>G']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.7+123566A>C'
        assert results['NC_000023.11:g.33215693T>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.33215693T>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.33215693T>G']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None

    def test_variant195(self):
        variant = 'NC_000023.11:g.33211557T>C'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.33211557T>C' in results.keys()
        assert results['NC_000023.11:g.33211557T>C']['p_vcf'] == 'X:33211557:T:C'
        assert results['NC_000023.11:g.33211557T>C']['g_hgvs'] == 'NC_000023.11:g.33211557T>C'
        assert results['NC_000023.11:g.33211557T>C']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.33211557T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.33211557T>C']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.7+127702A>G'
        assert results['NC_000023.11:g.33211557T>C']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.33211557T>C']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.33211557T>C']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None

    def test_variant196(self):
        variant = 'NC_000023.11:g.33211556A>T'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.33211556A>T' in results.keys()
        assert results['NC_000023.11:g.33211556A>T']['p_vcf'] == 'X:33211556:A:T'
        assert results['NC_000023.11:g.33211556A>T']['g_hgvs'] == 'NC_000023.11:g.33211556A>T'
        assert results['NC_000023.11:g.33211556A>T']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.33211556A>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.33211556A>T']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.7+127703T>A'
        assert results['NC_000023.11:g.33211556A>T']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.33211556A>T']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.33211556A>T']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.33211556A>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.33211556A>T']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.-244T>A'
        assert results['NC_000023.11:g.33211556A>T']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.33211556A>T']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.33211556A>T']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None

    def test_variant197(self):
        variant = 'NC_000023.11:g.33211450C>G'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.33211450C>G' in results.keys()
        assert results['NC_000023.11:g.33211450C>G']['p_vcf'] == 'X:33211450:C:G'
        assert results['NC_000023.11:g.33211450C>G']['g_hgvs'] == 'NC_000023.11:g.33211450C>G'
        assert results['NC_000023.11:g.33211450C>G']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.33211450C>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.33211450C>G']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.7+127809G>C'
        assert results['NC_000023.11:g.33211450C>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.33211450C>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.33211450C>G']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.33211450C>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.33211450C>G']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.-138G>C'
        assert results['NC_000023.11:g.33211450C>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.33211450C>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.33211450C>G']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None

    def test_variant198(self):
        variant = 'NC_000006.12:g.152637185C>A'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000006.12:g.152637185C>A' in results.keys()
        assert results['NC_000006.12:g.152637185C>A']['p_vcf'] == '6:152637185:C:A'
        assert results['NC_000006.12:g.152637185C>A']['g_hgvs'] == 'NC_000006.12:g.152637185C>A'
        assert results['NC_000006.12:g.152637185C>A']['genomic_variant_error'] is None
        assert 'NM_182961.3' in results['NC_000006.12:g.152637185C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000006.12:g.152637185C>A']['hgvs_t_and_p']['NM_182961.3']['t_hgvs'] == 'NM_182961.3:c.-392+4G>T'
        assert results['NC_000006.12:g.152637185C>A']['hgvs_t_and_p']['NM_182961.3']['p_hgvs_tlc'] == 'NP_892006.3:p.?'
        assert results['NC_000006.12:g.152637185C>A']['hgvs_t_and_p']['NM_182961.3']['p_hgvs_slc'] == 'NP_892006.3:p.?'
        assert results['NC_000006.12:g.152637185C>A']['hgvs_t_and_p']['NM_182961.3']['transcript_variant_error'] is None

    def test_variant199(self):
        variant = 'NC_000006.12:g.152636926C>A'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000006.12:g.152636926C>A' in results.keys()
        assert results['NC_000006.12:g.152636926C>A']['p_vcf'] == '6:152636926:C:A'
        assert results['NC_000006.12:g.152636926C>A']['g_hgvs'] == 'NC_000006.12:g.152636926C>A'
        assert results['NC_000006.12:g.152636926C>A']['genomic_variant_error'] is None
        assert 'NM_182961.3' in results['NC_000006.12:g.152636926C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000006.12:g.152636926C>A']['hgvs_t_and_p']['NM_182961.3']['t_hgvs'] == 'NM_182961.3:c.-391-121G>T'
        assert results['NC_000006.12:g.152636926C>A']['hgvs_t_and_p']['NM_182961.3']['p_hgvs_tlc'] == 'NP_892006.3:p.?'
        assert results['NC_000006.12:g.152636926C>A']['hgvs_t_and_p']['NM_182961.3']['p_hgvs_slc'] == 'NP_892006.3:p.?'
        assert results['NC_000006.12:g.152636926C>A']['hgvs_t_and_p']['NM_182961.3']['transcript_variant_error'] is None

    def test_variant200(self):
        variant = 'NC_000023.11:g.33211312T>C'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.33211312T>C' in results.keys()
        assert results['NC_000023.11:g.33211312T>C']['p_vcf'] == 'X:33211312:T:C'
        assert results['NC_000023.11:g.33211312T>C']['g_hgvs'] == 'NC_000023.11:g.33211312T>C'
        assert results['NC_000023.11:g.33211312T>C']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.33211312T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.33211312T>C']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.7+127947A>G'
        assert results['NC_000023.11:g.33211312T>C']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.33211312T>C']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.33211312T>C']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.33211312T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.33211312T>C']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.1A>G'
        assert results['NC_000023.11:g.33211312T>C']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.(Met1?)'
        assert results['NC_000023.11:g.33211312T>C']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.(M1?)'
        assert results['NC_000023.11:g.33211312T>C']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None

    def test_variant201(self):
        variant = 'NC_000023.11:g.33211311A>T'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.33211311A>T' in results.keys()
        assert results['NC_000023.11:g.33211311A>T']['p_vcf'] == 'X:33211311:A:T'
        assert results['NC_000023.11:g.33211311A>T']['g_hgvs'] == 'NC_000023.11:g.33211311A>T'
        assert results['NC_000023.11:g.33211311A>T']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.33211311A>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.33211311A>T']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.7+127948T>A'
        assert results['NC_000023.11:g.33211311A>T']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.33211311A>T']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.33211311A>T']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.33211311A>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.33211311A>T']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.2T>A'
        assert results['NC_000023.11:g.33211311A>T']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.(Met1?)'
        assert results['NC_000023.11:g.33211311A>T']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.(M1?)'
        assert results['NC_000023.11:g.33211311A>T']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None

    def test_variant202(self):
        variant = 'NC_000023.11:g.33211310C>G'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.33211310C>G' in results.keys()
        assert results['NC_000023.11:g.33211310C>G']['p_vcf'] == 'X:33211310:C:G'
        assert results['NC_000023.11:g.33211310C>G']['g_hgvs'] == 'NC_000023.11:g.33211310C>G'
        assert results['NC_000023.11:g.33211310C>G']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.33211310C>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.33211310C>G']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.7+127949G>C'
        assert results['NC_000023.11:g.33211310C>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.33211310C>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.33211310C>G']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.33211310C>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.33211310C>G']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.3G>C'
        assert results['NC_000023.11:g.33211310C>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.(Met1?)'
        assert results['NC_000023.11:g.33211310C>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.(M1?)'
        assert results['NC_000023.11:g.33211310C>G']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None

    def test_variant203(self):
        variant = 'NC_000023.11:g.32849793A>T'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32849793A>T' in results.keys()
        assert results['NC_000023.11:g.32849793A>T']['p_vcf'] == 'X:32849793:A:T'
        assert results['NC_000023.11:g.32849793A>T']['g_hgvs'] == 'NC_000023.11:g.32849793A>T'
        assert results['NC_000023.11:g.32849793A>T']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.97T>A'
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.(Phe33Ile)'
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.(F33I)'
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.109T>A'
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.(Phe37Ile)'
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.(F37I)'
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.121T>A'
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.(Phe41Ile)'
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.(F41I)'
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-249T>A'
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32849793A>T']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant204(self):
        variant = 'NC_000023.11:g.32823295C>G'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32823295C>G' in results.keys()
        assert results['NC_000023.11:g.32823295C>G']['p_vcf'] == 'X:32823295:C:G'
        assert results['NC_000023.11:g.32823295C>G']['g_hgvs'] == 'NC_000023.11:g.32823295C>G'
        assert results['NC_000023.11:g.32823295C>G']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.333G>C'
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.(Gln111His)'
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.(Q111H)'
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.345G>C'
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.(Gln115His)'
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.(Q115H)'
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.357G>C'
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.(Gln119His)'
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.(Q119H)'
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-13G>C'
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823295C>G']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant205(self):
        variant = 'NC_000023.11:g.32823294C>A'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32823294C>A' in results.keys()
        assert results['NC_000023.11:g.32823294C>A']['p_vcf'] == 'X:32823294:C:A'
        assert results['NC_000023.11:g.32823294C>A']['g_hgvs'] == 'NC_000023.11:g.32823294C>A'
        assert results['NC_000023.11:g.32823294C>A']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.333+1G>T'
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.345+1G>T'
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.357+1G>T'
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-13+1G>T'
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823294C>A']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant206(self):
        variant = 'NC_000023.11:g.32823293A>G'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32823293A>G' in results.keys()
        assert results['NC_000023.11:g.32823293A>G']['p_vcf'] == 'X:32823293:A:G'
        assert results['NC_000023.11:g.32823293A>G']['g_hgvs'] == 'NC_000023.11:g.32823293A>G'
        assert results['NC_000023.11:g.32823293A>G']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.333+2T>C'
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.345+2T>C'
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.357+2T>C'
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-13+2T>C'
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823293A>G']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant207(self):
        variant = 'NC_000023.11:g.32823292T>G'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32823292T>G' in results.keys()
        assert results['NC_000023.11:g.32823292T>G']['p_vcf'] == 'X:32823292:T:G'
        assert results['NC_000023.11:g.32823292T>G']['g_hgvs'] == 'NC_000023.11:g.32823292T>G'
        assert results['NC_000023.11:g.32823292T>G']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.333+3A>C'
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.345+3A>C'
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.357+3A>C'
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-13+3A>C'
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823292T>G']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant208(self):
        variant = 'NC_000023.11:g.32823291T>G'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32823291T>G' in results.keys()
        assert results['NC_000023.11:g.32823291T>G']['p_vcf'] == 'X:32823291:T:G'
        assert results['NC_000023.11:g.32823291T>G']['g_hgvs'] == 'NC_000023.11:g.32823291T>G'
        assert results['NC_000023.11:g.32823291T>G']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.333+4A>C'
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.345+4A>C'
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.357+4A>C'
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-13+4A>C'
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823291T>G']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant209(self):
        variant = 'NC_000023.11:g.32823290C>G'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32823290C>G' in results.keys()
        assert results['NC_000023.11:g.32823290C>G']['p_vcf'] == 'X:32823290:C:G'
        assert results['NC_000023.11:g.32823290C>G']['g_hgvs'] == 'NC_000023.11:g.32823290C>G'
        assert results['NC_000023.11:g.32823290C>G']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.333+5G>C'
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.345+5G>C'
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.357+5G>C'
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-13+5G>C'
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823290C>G']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant210(self):
        variant = 'NC_000023.11:g.32823289T>A'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32823289T>A' in results.keys()
        assert results['NC_000023.11:g.32823289T>A']['p_vcf'] == 'X:32823289:T:A'
        assert results['NC_000023.11:g.32823289T>A']['g_hgvs'] == 'NC_000023.11:g.32823289T>A'
        assert results['NC_000023.11:g.32823289T>A']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.333+6A>T'
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.345+6A>T'
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.357+6A>T'
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-13+6A>T'
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823289T>A']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant211(self):
        variant = 'NC_000023.11:g.32823288T>G'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32823288T>G' in results.keys()
        assert results['NC_000023.11:g.32823288T>G']['p_vcf'] == 'X:32823288:T:G'
        assert results['NC_000023.11:g.32823288T>G']['g_hgvs'] == 'NC_000023.11:g.32823288T>G'
        assert results['NC_000023.11:g.32823288T>G']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.333+7A>C'
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.345+7A>C'
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.357+7A>C'
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-13+7A>C'
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32823288T>G']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant212(self):
        variant = 'NC_000023.11:g.32820903T>A'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32820903T>A' in results.keys()
        assert results['NC_000023.11:g.32820903T>A']['p_vcf'] == 'X:32820903:T:A'
        assert results['NC_000023.11:g.32820903T>A']['g_hgvs'] == 'NC_000023.11:g.32820903T>A'
        assert results['NC_000023.11:g.32820903T>A']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.333+2392A>T'
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.345+2392A>T'
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.357+2392A>T'
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-13+2392A>T'
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32820903T>A']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant213(self):
        variant = 'NC_000023.11:g.32819968A>G'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32819968A>G' in results.keys()
        assert results['NC_000023.11:g.32819968A>G']['p_vcf'] == 'X:32819968:A:G'
        assert results['NC_000023.11:g.32819968A>G']['g_hgvs'] == 'NC_000023.11:g.32819968A>G'
        assert results['NC_000023.11:g.32819968A>G']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.333+3327T>C'
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.345+3327T>C'
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.357+3327T>C'
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-13+3327T>C'
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32819968A>G']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant214(self):
        variant = 'NC_000023.11:g.32698556G>C'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32698556G>C' in results.keys()
        assert results['NC_000023.11:g.32698556G>C']['p_vcf'] == 'X:32698556:G:C'
        assert results['NC_000023.11:g.32698556G>C']['g_hgvs'] == 'NC_000023.11:g.32698556G>C'
        assert results['NC_000023.11:g.32698556G>C']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.807+556C>G'
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.819+556C>G'
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.831+556C>G'
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.462+556C>G'
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32698556G>C']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    # Center of the intron.
    def test_variant215(self):
        variant = 'NC_000023.11:g.32698555T>G'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32698555T>G' in results.keys()
        assert results['NC_000023.11:g.32698555T>G']['p_vcf'] == 'X:32698555:T:G'
        assert results['NC_000023.11:g.32698555T>G']['g_hgvs'] == 'NC_000023.11:g.32698555T>G'
        assert results['NC_000023.11:g.32698555T>G']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.807+557A>C'
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.819+557A>C'
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.831+557A>C'
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.462+557A>C'
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32698555T>G']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant216(self):
        variant = 'NC_000023.11:g.32698554A>C'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32698554A>C' in results.keys()
        assert results['NC_000023.11:g.32698554A>C']['p_vcf'] == 'X:32698554:A:C'
        assert results['NC_000023.11:g.32698554A>C']['g_hgvs'] == 'NC_000023.11:g.32698554A>C'
        assert results['NC_000023.11:g.32698554A>C']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.808-556T>G'
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.820-556T>G'
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.832-556T>G'
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.463-556T>G'
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32698554A>C']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant217(self):
        variant = 'NC_000023.11:g.32819967G>T'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32819967G>T' in results.keys()
        assert results['NC_000023.11:g.32819967G>T']['p_vcf'] == 'X:32819967:G:T'
        assert results['NC_000023.11:g.32819967G>T']['g_hgvs'] == 'NC_000023.11:g.32819967G>T'
        assert results['NC_000023.11:g.32819967G>T']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.334-3327C>A'
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.346-3327C>A'
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.358-3327C>A'
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-12-3327C>A'
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32819967G>T']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant218(self):
        variant = 'NC_000023.11:g.32819475T>A'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32819475T>A' in results.keys()
        assert results['NC_000023.11:g.32819475T>A']['p_vcf'] == 'X:32819475:T:A'
        assert results['NC_000023.11:g.32819475T>A']['g_hgvs'] == 'NC_000023.11:g.32819475T>A'
        assert results['NC_000023.11:g.32819475T>A']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.334-2835A>T'
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.346-2835A>T'
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.358-2835A>T'
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-12-2835A>T'
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32819475T>A']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant219(self):
        variant = 'NC_000023.11:g.32816643A>C'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32816643A>C' in results.keys()
        assert results['NC_000023.11:g.32816643A>C']['p_vcf'] == 'X:32816643:A:C'
        assert results['NC_000023.11:g.32816643A>C']['g_hgvs'] == 'NC_000023.11:g.32816643A>C'
        assert results['NC_000023.11:g.32816643A>C']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.334-3T>G'
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.346-3T>G'
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.358-3T>G'
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-12-3T>G'
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32816643A>C']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant220(self):
        variant = 'NC_000023.11:g.32816642T>A'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32816642T>A' in results.keys()
        assert results['NC_000023.11:g.32816642T>A']['p_vcf'] == 'X:32816642:T:A'
        assert results['NC_000023.11:g.32816642T>A']['g_hgvs'] == 'NC_000023.11:g.32816642T>A'
        assert results['NC_000023.11:g.32816642T>A']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.334-2A>T'
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.346-2A>T'
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.358-2A>T'
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-12-2A>T'
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32816642T>A']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant221(self):
        variant = 'NC_000023.11:g.32816641C>T'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32816641C>T' in results.keys()
        assert results['NC_000023.11:g.32816641C>T']['p_vcf'] == 'X:32816641:C:T'
        assert results['NC_000023.11:g.32816641C>T']['g_hgvs'] == 'NC_000023.11:g.32816641C>T'
        assert results['NC_000023.11:g.32816641C>T']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.334-1G>A'
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.346-1G>A'
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.358-1G>A'
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-12-1G>A'
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32816641C>T']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant222(self):
        variant = 'NC_000023.11:g.32816640C>T'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32816640C>T' in results.keys()
        assert results['NC_000023.11:g.32816640C>T']['p_vcf'] == 'X:32816640:C:T'
        assert results['NC_000023.11:g.32816640C>T']['g_hgvs'] == 'NC_000023.11:g.32816640C>T'
        assert results['NC_000023.11:g.32816640C>T']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.334G>A'
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.(Val112Ile)'
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.(V112I)'
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.346G>A'
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.(Val116Ile)'
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.(V116I)'
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.358G>A'
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.(Val120Ile)'
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.(V120I)'
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.-12G>A'
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32816640C>T']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant223(self):
        variant = 'NC_000023.11:g.31121921A>G'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.31121921A>G' in results.keys()
        assert results['NC_000023.11:g.31121921A>G']['p_vcf'] == 'X:31121921:A:G'
        assert results['NC_000023.11:g.31121921A>G']['g_hgvs'] == 'NC_000023.11:g.31121921A>G'
        assert results['NC_000023.11:g.31121921A>G']['genomic_variant_error'] is None
        assert 'NM_004011.3' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004011.3']['t_hgvs'] == 'NM_004011.3:c.7033T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004011.3']['p_hgvs_tlc'] == 'NP_004002.2:p.(Ter2345GlnextTer17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004011.3']['p_hgvs_slc'] == 'NP_004002.2:p.(*2345Qext*17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004011.3']['transcript_variant_error'] is None
        assert 'NM_004021.2' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004021.2']['t_hgvs'] == 'NM_004021.2:c.3644T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004021.2']['p_hgvs_tlc'] == 'NP_004012.1:p.(Val1215Ala)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004021.2']['p_hgvs_slc'] == 'NP_004012.1:p.(V1215A)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004021.2']['transcript_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.11032T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.(Ter3678GlnextTer17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.(*3678Qext*17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.11044T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.(Ter3682GlnextTer17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.(*3682Qext*17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004014.2' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004014.2']['t_hgvs'] == 'NM_004014.2:c.2869T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004014.2']['p_hgvs_tlc'] == 'NP_004005.1:p.(Ter957GlnextTer17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004014.2']['p_hgvs_slc'] == 'NP_004005.1:p.(*957Qext*17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004014.2']['transcript_variant_error'] is None
        assert 'NM_004015.2' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004015.2']['t_hgvs'] == 'NM_004015.2:c.1852T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004015.2']['p_hgvs_tlc'] == 'NP_004006.1:p.(Ter618GlnextTer17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004015.2']['p_hgvs_slc'] == 'NP_004006.1:p.(*618Qext*17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004015.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.10687T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.(Ter3563GlnextTer17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.(*3563Qext*17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None
        assert 'NM_004013.2' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004013.2']['t_hgvs'] == 'NM_004013.2:c.3676T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004013.2']['p_hgvs_tlc'] == 'NP_004004.1:p.(Ter1226GlnextTer17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004013.2']['p_hgvs_slc'] == 'NP_004004.1:p.(*1226Qext*17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004013.2']['transcript_variant_error'] is None
        assert 'NM_004022.2' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004022.2']['t_hgvs'] == 'NM_004022.2:c.3605T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004022.2']['p_hgvs_tlc'] == 'NP_004013.1:p.(Val1202Ala)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004022.2']['p_hgvs_slc'] == 'NP_004013.1:p.(V1202A)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004022.2']['transcript_variant_error'] is None
        assert 'NM_004023.2' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004023.2']['t_hgvs'] == 'NM_004023.2:c.3314T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004023.2']['p_hgvs_tlc'] == 'NP_004014.1:p.(Val1105Ala)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004023.2']['p_hgvs_slc'] == 'NP_004014.1:p.(V1105A)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004023.2']['transcript_variant_error'] is None
        assert 'NM_004018.2' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004018.2']['t_hgvs'] == 'NM_004018.2:c.1781T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004018.2']['p_hgvs_tlc'] == 'NP_004009.1:p.(Val594Ala)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004018.2']['p_hgvs_slc'] == 'NP_004009.1:p.(V594A)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004018.2']['transcript_variant_error'] is None
        assert 'NM_004017.2' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004017.2']['t_hgvs'] == 'NM_004017.2:c.1813T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004017.2']['p_hgvs_tlc'] == 'NP_004008.1:p.(Ter605GlnextTer17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004017.2']['p_hgvs_slc'] == 'NP_004008.1:p.(*605Qext*17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004017.2']['transcript_variant_error'] is None
        assert 'NM_004016.2' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004016.2']['t_hgvs'] == 'NM_004016.2:c.1820T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004016.2']['p_hgvs_tlc'] == 'NP_004007.1:p.(Val607Ala)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004016.2']['p_hgvs_slc'] == 'NP_004007.1:p.(V607A)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004016.2']['transcript_variant_error'] is None
        assert 'NM_004020.3' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004020.3']['t_hgvs'] == 'NM_004020.3:c.3346T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004020.3']['p_hgvs_tlc'] == 'NP_004011.2:p.(Ter1116GlnextTer17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004020.3']['p_hgvs_slc'] == 'NP_004011.2:p.(*1116Qext*17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004020.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.11056T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.(Ter3686GlnextTer17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.(*3686Qext*17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004012.3' in results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004012.3']['t_hgvs'] == 'NM_004012.3:c.7024T>C'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004012.3']['p_hgvs_tlc'] == 'NP_004003.1:p.(Ter2342GlnextTer17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004012.3']['p_hgvs_slc'] == 'NP_004003.1:p.(*2342Qext*17)'
        assert results['NC_000023.11:g.31121921A>G']['hgvs_t_and_p']['NM_004012.3']['transcript_variant_error'] is None

    def test_variant224(self):
        variant = 'NC_000023.11:g.31121920T>C'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.31121920T>C' in results.keys()
        assert results['NC_000023.11:g.31121920T>C']['p_vcf'] == 'X:31121920:T:C'
        assert results['NC_000023.11:g.31121920T>C']['g_hgvs'] == 'NC_000023.11:g.31121920T>C'
        assert results['NC_000023.11:g.31121920T>C']['genomic_variant_error'] is None
        assert 'NM_004011.3' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004011.3']['t_hgvs'] == 'NM_004011.3:c.7034A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004011.3']['p_hgvs_tlc'] == 'NP_004002.2:p.(Ter2345TrpextTer17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004011.3']['p_hgvs_slc'] == 'NP_004002.2:p.(*2345Wext*17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004011.3']['transcript_variant_error'] is None
        assert 'NM_004021.2' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004021.2']['t_hgvs'] == 'NM_004021.2:c.3645A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004021.2']['p_hgvs_tlc'] == 'NP_004012.1:p.(Val1215=)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004021.2']['p_hgvs_slc'] == 'NP_004012.1:p.(V1215=)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004021.2']['transcript_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.11033A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.(Ter3678TrpextTer17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.(*3678Wext*17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.11045A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.(Ter3682TrpextTer17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.(*3682Wext*17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004014.2' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004014.2']['t_hgvs'] == 'NM_004014.2:c.2870A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004014.2']['p_hgvs_tlc'] == 'NP_004005.1:p.(Ter957TrpextTer17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004014.2']['p_hgvs_slc'] == 'NP_004005.1:p.(*957Wext*17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004014.2']['transcript_variant_error'] is None
        assert 'NM_004015.2' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004015.2']['t_hgvs'] == 'NM_004015.2:c.1853A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004015.2']['p_hgvs_tlc'] == 'NP_004006.1:p.(Ter618TrpextTer17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004015.2']['p_hgvs_slc'] == 'NP_004006.1:p.(*618Wext*17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004015.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.10688A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.(Ter3563TrpextTer17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.(*3563Wext*17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None
        assert 'NM_004013.2' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004013.2']['t_hgvs'] == 'NM_004013.2:c.3677A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004013.2']['p_hgvs_tlc'] == 'NP_004004.1:p.(Ter1226TrpextTer17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004013.2']['p_hgvs_slc'] == 'NP_004004.1:p.(*1226Wext*17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004013.2']['transcript_variant_error'] is None
        assert 'NM_004022.2' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004022.2']['t_hgvs'] == 'NM_004022.2:c.3606A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004022.2']['p_hgvs_tlc'] == 'NP_004013.1:p.(Val1202=)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004022.2']['p_hgvs_slc'] == 'NP_004013.1:p.(V1202=)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004022.2']['transcript_variant_error'] is None
        assert 'NM_004023.2' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004023.2']['t_hgvs'] == 'NM_004023.2:c.3315A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004023.2']['p_hgvs_tlc'] == 'NP_004014.1:p.(Val1105=)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004023.2']['p_hgvs_slc'] == 'NP_004014.1:p.(V1105=)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004023.2']['transcript_variant_error'] is None
        assert 'NM_004018.2' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004018.2']['t_hgvs'] == 'NM_004018.2:c.1782A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004018.2']['p_hgvs_tlc'] == 'NP_004009.1:p.(Val594=)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004018.2']['p_hgvs_slc'] == 'NP_004009.1:p.(V594=)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004018.2']['transcript_variant_error'] is None
        assert 'NM_004017.2' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004017.2']['t_hgvs'] == 'NM_004017.2:c.1814A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004017.2']['p_hgvs_tlc'] == 'NP_004008.1:p.(Ter605TrpextTer17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004017.2']['p_hgvs_slc'] == 'NP_004008.1:p.(*605Wext*17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004017.2']['transcript_variant_error'] is None
        assert 'NM_004016.2' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004016.2']['t_hgvs'] == 'NM_004016.2:c.1821A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004016.2']['p_hgvs_tlc'] == 'NP_004007.1:p.(Val607=)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004016.2']['p_hgvs_slc'] == 'NP_004007.1:p.(V607=)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004016.2']['transcript_variant_error'] is None
        assert 'NM_004020.3' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004020.3']['t_hgvs'] == 'NM_004020.3:c.3347A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004020.3']['p_hgvs_tlc'] == 'NP_004011.2:p.(Ter1116TrpextTer17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004020.3']['p_hgvs_slc'] == 'NP_004011.2:p.(*1116Wext*17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004020.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.11057A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.(Ter3686TrpextTer17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.(*3686Wext*17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004012.3' in results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004012.3']['t_hgvs'] == 'NM_004012.3:c.7025A>G'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004012.3']['p_hgvs_tlc'] == 'NP_004003.1:p.(Ter2342TrpextTer17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004012.3']['p_hgvs_slc'] == 'NP_004003.1:p.(*2342Wext*17)'
        assert results['NC_000023.11:g.31121920T>C']['hgvs_t_and_p']['NM_004012.3']['transcript_variant_error'] is None

    def test_variant225(self):
        variant = 'NC_000023.11:g.31121919C>A'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.31121919C>A' in results.keys()
        assert results['NC_000023.11:g.31121919C>A']['p_vcf'] == 'X:31121919:C:A'
        assert results['NC_000023.11:g.31121919C>A']['g_hgvs'] == 'NC_000023.11:g.31121919C>A'
        assert results['NC_000023.11:g.31121919C>A']['genomic_variant_error'] is None
        assert 'NM_004011.3' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004011.3']['t_hgvs'] == 'NM_004011.3:c.7035G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004011.3']['p_hgvs_tlc'] == 'NP_004002.2:p.(Ter2345TyrextTer17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004011.3']['p_hgvs_slc'] == 'NP_004002.2:p.(*2345Yext*17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004011.3']['transcript_variant_error'] is None
        assert 'NM_004021.2' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004021.2']['t_hgvs'] == 'NM_004021.2:c.3646G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004021.2']['p_hgvs_tlc'] == 'NP_004012.1:p.(Gly1216Ter)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004021.2']['p_hgvs_slc'] == 'NP_004012.1:p.(G1216*)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004021.2']['transcript_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.11034G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.(Ter3678TyrextTer17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.(*3678Yext*17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.11046G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.(Ter3682TyrextTer17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.(*3682Yext*17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004014.2' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004014.2']['t_hgvs'] == 'NM_004014.2:c.2871G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004014.2']['p_hgvs_tlc'] == 'NP_004005.1:p.(Ter957TyrextTer17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004014.2']['p_hgvs_slc'] == 'NP_004005.1:p.(*957Yext*17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004014.2']['transcript_variant_error'] is None
        assert 'NM_004015.2' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004015.2']['t_hgvs'] == 'NM_004015.2:c.1854G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004015.2']['p_hgvs_tlc'] == 'NP_004006.1:p.(Ter618TyrextTer17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004015.2']['p_hgvs_slc'] == 'NP_004006.1:p.(*618Yext*17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004015.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.10689G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.(Ter3563TyrextTer17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.(*3563Yext*17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None
        assert 'NM_004013.2' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004013.2']['t_hgvs'] == 'NM_004013.2:c.3678G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004013.2']['p_hgvs_tlc'] == 'NP_004004.1:p.(Ter1226TyrextTer17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004013.2']['p_hgvs_slc'] == 'NP_004004.1:p.(*1226Yext*17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004013.2']['transcript_variant_error'] is None
        assert 'NM_004022.2' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004022.2']['t_hgvs'] == 'NM_004022.2:c.3607G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004022.2']['p_hgvs_tlc'] == 'NP_004013.1:p.(Gly1203Ter)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004022.2']['p_hgvs_slc'] == 'NP_004013.1:p.(G1203*)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004022.2']['transcript_variant_error'] is None
        assert 'NM_004023.2' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004023.2']['t_hgvs'] == 'NM_004023.2:c.3316G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004023.2']['p_hgvs_tlc'] == 'NP_004014.1:p.(Gly1106Ter)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004023.2']['p_hgvs_slc'] == 'NP_004014.1:p.(G1106*)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004023.2']['transcript_variant_error'] is None
        assert 'NM_004018.2' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004018.2']['t_hgvs'] == 'NM_004018.2:c.1783G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004018.2']['p_hgvs_tlc'] == 'NP_004009.1:p.(Gly595Ter)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004018.2']['p_hgvs_slc'] == 'NP_004009.1:p.(G595*)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004018.2']['transcript_variant_error'] is None
        assert 'NM_004017.2' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004017.2']['t_hgvs'] == 'NM_004017.2:c.1815G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004017.2']['p_hgvs_tlc'] == 'NP_004008.1:p.(Ter605TyrextTer17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004017.2']['p_hgvs_slc'] == 'NP_004008.1:p.(*605Yext*17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004017.2']['transcript_variant_error'] is None
        assert 'NM_004016.2' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004016.2']['t_hgvs'] == 'NM_004016.2:c.1822G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004016.2']['p_hgvs_tlc'] == 'NP_004007.1:p.(Gly608Ter)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004016.2']['p_hgvs_slc'] == 'NP_004007.1:p.(G608*)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004016.2']['transcript_variant_error'] is None
        assert 'NM_004020.3' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004020.3']['t_hgvs'] == 'NM_004020.3:c.3348G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004020.3']['p_hgvs_tlc'] == 'NP_004011.2:p.(Ter1116TyrextTer17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004020.3']['p_hgvs_slc'] == 'NP_004011.2:p.(*1116Yext*17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004020.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.11058G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.(Ter3686TyrextTer17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.(*3686Yext*17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004012.3' in results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004012.3']['t_hgvs'] == 'NM_004012.3:c.7026G>T'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004012.3']['p_hgvs_tlc'] == 'NP_004003.1:p.(Ter2342TyrextTer17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004012.3']['p_hgvs_slc'] == 'NP_004003.1:p.(*2342Yext*17)'
        assert results['NC_000023.11:g.31121919C>A']['hgvs_t_and_p']['NM_004012.3']['transcript_variant_error'] is None

    def test_variant226(self):
        variant = 'NC_000023.11:g.31121131T>G'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.31121131T>G' in results.keys()
        assert results['NC_000023.11:g.31121131T>G']['p_vcf'] == 'X:31121131:T:G'
        assert results['NC_000023.11:g.31121131T>G']['g_hgvs'] == 'NC_000023.11:g.31121131T>G'
        assert results['NC_000023.11:g.31121131T>G']['genomic_variant_error'] is None
        assert 'NM_004011.3' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004011.3']['t_hgvs'] == 'NM_004011.3:c.*788A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004011.3']['p_hgvs_tlc'] == 'NP_004002.2:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004011.3']['p_hgvs_slc'] == 'NP_004002.2:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004011.3']['transcript_variant_error'] is None
        assert 'NM_004021.2' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004021.2']['t_hgvs'] == 'NM_004021.2:c.*702A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004021.2']['p_hgvs_tlc'] == 'NP_004012.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004021.2']['p_hgvs_slc'] == 'NP_004012.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004021.2']['transcript_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.*788A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.*788A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004014.2' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004014.2']['t_hgvs'] == 'NM_004014.2:c.*788A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004014.2']['p_hgvs_tlc'] == 'NP_004005.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004014.2']['p_hgvs_slc'] == 'NP_004005.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004014.2']['transcript_variant_error'] is None
        assert 'NM_004015.2' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004015.2']['t_hgvs'] == 'NM_004015.2:c.*788A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004015.2']['p_hgvs_tlc'] == 'NP_004006.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004015.2']['p_hgvs_slc'] == 'NP_004006.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004015.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.*788A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None
        assert 'NM_004013.2' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004013.2']['t_hgvs'] == 'NM_004013.2:c.*788A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004013.2']['p_hgvs_tlc'] == 'NP_004004.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004013.2']['p_hgvs_slc'] == 'NP_004004.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004013.2']['transcript_variant_error'] is None
        assert 'NM_004022.2' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004022.2']['t_hgvs'] == 'NM_004022.2:c.*702A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004022.2']['p_hgvs_tlc'] == 'NP_004013.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004022.2']['p_hgvs_slc'] == 'NP_004013.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004022.2']['transcript_variant_error'] is None
        assert 'NM_004023.2' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004023.2']['t_hgvs'] == 'NM_004023.2:c.*702A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004023.2']['p_hgvs_tlc'] == 'NP_004014.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004023.2']['p_hgvs_slc'] == 'NP_004014.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004023.2']['transcript_variant_error'] is None
        assert 'NM_004018.2' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004018.2']['t_hgvs'] == 'NM_004018.2:c.*702A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004018.2']['p_hgvs_tlc'] == 'NP_004009.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004018.2']['p_hgvs_slc'] == 'NP_004009.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004018.2']['transcript_variant_error'] is None
        assert 'NM_004017.2' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004017.2']['t_hgvs'] == 'NM_004017.2:c.*788A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004017.2']['p_hgvs_tlc'] == 'NP_004008.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004017.2']['p_hgvs_slc'] == 'NP_004008.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004017.2']['transcript_variant_error'] is None
        assert 'NM_004016.2' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004016.2']['t_hgvs'] == 'NM_004016.2:c.*702A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004016.2']['p_hgvs_tlc'] == 'NP_004007.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004016.2']['p_hgvs_slc'] == 'NP_004007.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004016.2']['transcript_variant_error'] is None
        assert 'NM_004020.3' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004020.3']['t_hgvs'] == 'NM_004020.3:c.*788A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004020.3']['p_hgvs_tlc'] == 'NP_004011.2:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004020.3']['p_hgvs_slc'] == 'NP_004011.2:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004020.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.*788A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004012.3' in results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004012.3']['t_hgvs'] == 'NM_004012.3:c.*788A>C'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004012.3']['p_hgvs_tlc'] == 'NP_004003.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004012.3']['p_hgvs_slc'] == 'NP_004003.1:p.?'
        assert results['NC_000023.11:g.31121131T>G']['hgvs_t_and_p']['NM_004012.3']['transcript_variant_error'] is None

    # Test removed due to major changes to the transcript versions an QC requirments between UTA and VVTA.
    def test_variant227(self):
        variant = 'NC_000023.11:g.31119228G>T'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.31119228G>T' in results.keys()
        assert results['NC_000023.11:g.31119228G>T']['p_vcf'] == 'X:31119228:G:T'
        assert results['NC_000023.11:g.31119228G>T']['g_hgvs'] == 'NC_000023.11:g.31119228G>T'
        assert results['NC_000023.11:g.31119228G>T']['genomic_variant_error'] is None

    def test_variant229(self):
        variant = 'NC_000023.11:g.31118748C>A'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.31118748C>A' in results.keys()
        assert results['NC_000023.11:g.31118748C>A']['p_vcf'] == 'X:31118748:C:A'
        assert results['NC_000023.11:g.31118748C>A']['g_hgvs'] == 'NC_000023.11:g.31118748C>A'
        assert results['NC_000023.11:g.31118748C>A']['genomic_variant_error'] is None
        assert results['NC_000023.11:g.31118748C>A']['hgvs_t_and_p'] == {'intergenic': {'alt_genomic_loci': None}}

    def test_variant230(self):
        variant = 'NC_000011.10:g.70487682T>A'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000011.10:g.70487682T>A' in results.keys()
        assert results['NC_000011.10:g.70487682T>A']['p_vcf'] == '11:70487682:T:A'
        assert results['NC_000011.10:g.70487682T>A']['g_hgvs'] == 'NC_000011.10:g.70487682T>A'
        assert results['NC_000011.10:g.70487682T>A']['genomic_variant_error'] is None
        assert 'NM_133266.4' in results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_133266.4']['t_hgvs'] == 'NM_133266.4:c.847A>T'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_133266.4']['p_hgvs_tlc'] == 'NP_573573.2:p.(Met283Leu)'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_133266.4']['p_hgvs_slc'] == 'NP_573573.2:p.(M283L)'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_133266.4']['transcript_variant_error'] is None
        assert 'NM_012309.4' in results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_012309.4']['t_hgvs'] == 'NM_012309.4:c.2611A>T'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_012309.4']['p_hgvs_tlc'] == 'NP_036441.2:p.(Met871Leu)'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_012309.4']['p_hgvs_slc'] == 'NP_036441.2:p.(M871L)'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_012309.4']['transcript_variant_error'] is None
        assert 'NR_110766.1' in results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NR_110766.1']['t_hgvs'] == 'NR_110766.1:n.833+2621A>T'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NR_110766.1']['p_hgvs_tlc'] is None
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NR_110766.1']['p_hgvs_slc'] is None
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NR_110766.1']['transcript_variant_error'] is None

    def test_variant231(self):
        variant = 'NC_000011.10:g.70487682T>A'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000011.10:g.70487682T>A' in results.keys()
        assert results['NC_000011.10:g.70487682T>A']['p_vcf'] == '11:70487682:T:A'
        assert results['NC_000011.10:g.70487682T>A']['g_hgvs'] == 'NC_000011.10:g.70487682T>A'
        assert results['NC_000011.10:g.70487682T>A']['genomic_variant_error'] is None
        assert 'NM_133266.4' in results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_133266.4']['t_hgvs'] == 'NM_133266.4:c.847A>T'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_133266.4']['p_hgvs_tlc'] == 'NP_573573.2:p.(Met283Leu)'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_133266.4']['p_hgvs_slc'] == 'NP_573573.2:p.(M283L)'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_133266.4']['transcript_variant_error'] is None
        assert 'NM_012309.4' in results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_012309.4']['t_hgvs'] == 'NM_012309.4:c.2611A>T'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_012309.4']['p_hgvs_tlc'] == 'NP_036441.2:p.(Met871Leu)'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_012309.4']['p_hgvs_slc'] == 'NP_036441.2:p.(M871L)'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_012309.4']['transcript_variant_error'] is None
        assert 'NR_110766.1' in results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NR_110766.1']['t_hgvs'] == 'NR_110766.1:n.833+2621A>T'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NR_110766.1']['p_hgvs_tlc'] is None
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NR_110766.1']['p_hgvs_slc'] is None
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NR_110766.1']['transcript_variant_error'] is None

    def test_variant232(self):
        variant = 'NC_000023.11:g.32485077T>C'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32485077T>C' in results.keys()
        assert results['NC_000023.11:g.32485077T>C']['p_vcf'] == 'X:32485077:T:C'
        assert results['NC_000023.11:g.32485077T>C']['g_hgvs'] == 'NC_000023.11:g.32485077T>C'
        assert results['NC_000023.11:g.32485077T>C']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.2621='
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.(Gly874=)'
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.(G874=)'
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.2633='
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.(Gly878=)'
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.(G878=)'
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.2645='
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.(Gly882=)'
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.(G882=)'
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.2276='
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.(Gly759=)'
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.(G759=)'
        assert results['NC_000023.11:g.32485077T>C']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant233(self):
        variant = 'NC_000011.10:g.70487682T>A'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000011.10:g.70487682T>A' in results.keys()
        assert results['NC_000011.10:g.70487682T>A']['p_vcf'] == '11:70487682:T:A'
        assert results['NC_000011.10:g.70487682T>A']['g_hgvs'] == 'NC_000011.10:g.70487682T>A'
        assert results['NC_000011.10:g.70487682T>A']['genomic_variant_error'] is None
        assert 'NM_133266.4' in results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_133266.4']['t_hgvs'] == 'NM_133266.4:c.847A>T'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_133266.4']['p_hgvs_tlc'] == 'NP_573573.2:p.(Met283Leu)'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_133266.4']['p_hgvs_slc'] == 'NP_573573.2:p.(M283L)'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_133266.4']['transcript_variant_error'] is None
        assert 'NM_012309.4' in results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_012309.4']['t_hgvs'] == 'NM_012309.4:c.2611A>T'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_012309.4']['p_hgvs_tlc'] == 'NP_036441.2:p.(Met871Leu)'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_012309.4']['p_hgvs_slc'] == 'NP_036441.2:p.(M871L)'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NM_012309.4']['transcript_variant_error'] is None
        assert 'NR_110766.1' in results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p'].keys()
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NR_110766.1']['t_hgvs'] == 'NR_110766.1:n.833+2621A>T'
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NR_110766.1']['p_hgvs_tlc'] is None
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NR_110766.1']['p_hgvs_slc'] is None
        assert results['NC_000011.10:g.70487682T>A']['hgvs_t_and_p']['NR_110766.1']['transcript_variant_error'] is None

    def test_variant234(self):
        variant = 'NC_000023.11:g.32503194A>N'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32503194A>N' in results.keys()
        assert results['NC_000023.11:g.32503194A>N']['p_vcf'] is None
        assert results['NC_000023.11:g.32503194A>N']['g_hgvs'] is None
        assert results['NC_000023.11:g.32503194A>N']['genomic_variant_error'] == 'NC_000023.11:g.32503194A>N: Variant reference (A) does not agree with reference sequence (C)'
        assert results['NC_000023.11:g.32503194A>N']['hgvs_t_and_p'] is None

    def test_variant235(self):
        variant = 'NC_000023.11:g.32503194C>N'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32503194C>N' in results.keys()
        assert results['NC_000023.11:g.32503194C>N']['p_vcf'] == 'X:32503194:C:N'
        assert results['NC_000023.11:g.32503194C>N']['g_hgvs'] == 'NC_000023.11:g.32503194C>N'
        assert results['NC_000023.11:g.32503194C>N']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.2269-1352G>N'
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.?'
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.2281-1352G>N'
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.?'
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.2293-1352G>N'
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.?'
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.1924-1352G>N'
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.?'
        assert results['NC_000023.11:g.32503194C>N']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant236(self):
        variant = 'NC_000023.11:g.32644229='
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32644229=' in results.keys()
        assert results['NC_000023.11:g.32644229=']['p_vcf'] == 'X:32644229:A:A'
        assert results['NC_000023.11:g.32644229=']['g_hgvs'] == 'NC_000023.11:g.32644229='
        assert results['NC_000023.11:g.32644229=']['genomic_variant_error'] is None
        assert 'NM_000109.3' in results['NC_000023.11:g.32644229=']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_000109.3']['t_hgvs'] == 'NM_000109.3:c.1210='
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_tlc'] == 'NP_000100.2:p.(Leu404=)'
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_000109.3']['p_hgvs_slc'] == 'NP_000100.2:p.(L404=)'
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_000109.3']['transcript_variant_error'] is None
        assert 'NM_004009.3' in results['NC_000023.11:g.32644229=']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_004009.3']['t_hgvs'] == 'NM_004009.3:c.1222='
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_tlc'] == 'NP_004000.1:p.(Leu408=)'
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_004009.3']['p_hgvs_slc'] == 'NP_004000.1:p.(L408=)'
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_004009.3']['transcript_variant_error'] is None
        assert 'NM_004006.2' in results['NC_000023.11:g.32644229=']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_004006.2']['t_hgvs'] == 'NM_004006.2:c.1234='
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_tlc'] == 'NP_003997.1:p.(Leu412=)'
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_004006.2']['p_hgvs_slc'] == 'NP_003997.1:p.(L412=)'
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_004006.2']['transcript_variant_error'] is None
        assert 'NM_004010.3' in results['NC_000023.11:g.32644229=']['hgvs_t_and_p'].keys()
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_004010.3']['t_hgvs'] == 'NM_004010.3:c.865='
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_tlc'] == 'NP_004001.1:p.(Leu289=)'
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_004010.3']['p_hgvs_slc'] == 'NP_004001.1:p.(L289=)'
        assert results['NC_000023.11:g.32644229=']['hgvs_t_and_p']['NM_004010.3']['transcript_variant_error'] is None

    def test_variant237(self):
        variant = 'NC_000023.11:g.32849761C>A'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000023.11:g.32849761C>A' in results.keys()
        assert results['NC_000023.11:g.32849761C>A']['p_vcf'] is None
        assert results['NC_000023.11:g.32849761C>A']['g_hgvs'] is None
        assert results['NC_000023.11:g.32849761C>A']['genomic_variant_error'] == 'NC_000023.11:g.32849761C>A: Variant reference (C) does not agree with reference sequence (T)'
        assert results['NC_000023.11:g.32849761C>A']['hgvs_t_and_p'] is None

    def test_variant238(self):
        variant = '14-105246588-TCT-T'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '14-105246588-TCT-T' in results.keys()
        assert results['14-105246588-TCT-T']['p_vcf'] is None
        assert results['14-105246588-TCT-T']['g_hgvs'] is None
        assert results['14-105246588-TCT-T']['genomic_variant_error'] == 'NC_000014.9:g.105246588_105246590delTCTinsT: Variant reference (TCT) does not agree with reference sequence (GCC)'
        assert results['14-105246588-TCT-T']['hgvs_t_and_p'] is None

    def test_variant239(self):
        variant = '11-108218120-C-CT'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '11-108218120-C-CT' in results.keys()
        assert results['11-108218120-C-CT']['p_vcf'] == '11-108218120-C-CT'
        assert results['11-108218120-C-CT']['g_hgvs'] == 'NC_000011.10:g.108218122dup'
        assert results['11-108218120-C-CT']['genomic_variant_error'] is None
        assert 'NM_002519.2' in results['11-108218120-C-CT']['hgvs_t_and_p'].keys()
        assert results['11-108218120-C-CT']['hgvs_t_and_p']['NM_002519.2']['t_hgvs'] == 'NM_002519.2:c.37+4379dup'
        assert results['11-108218120-C-CT']['hgvs_t_and_p']['NM_002519.2']['p_hgvs_tlc'] == 'NP_002510.2:p.?'
        assert results['11-108218120-C-CT']['hgvs_t_and_p']['NM_002519.2']['p_hgvs_slc'] == 'NP_002510.2:p.?'
        assert results['11-108218120-C-CT']['hgvs_t_and_p']['NM_002519.2']['transcript_variant_error'] is None
        assert 'NM_001321307.1' in results['11-108218120-C-CT']['hgvs_t_and_p'].keys()
        assert results['11-108218120-C-CT']['hgvs_t_and_p']['NM_001321307.1']['t_hgvs'] == 'NM_001321307.1:c.37+4379dup'
        assert results['11-108218120-C-CT']['hgvs_t_and_p']['NM_001321307.1']['p_hgvs_tlc'] == 'NP_001308236.1:p.?'
        assert results['11-108218120-C-CT']['hgvs_t_and_p']['NM_001321307.1']['p_hgvs_slc'] == 'NP_001308236.1:p.?'
        assert results['11-108218120-C-CT']['hgvs_t_and_p']['NM_001321307.1']['transcript_variant_error'] is None

    def test_variant240(self):
        variant = '11-108218120-CT-C'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert '11-108218120-CT-C' in results.keys()
        assert results['11-108218120-CT-C']['p_vcf'] == '11-108218120-CT-C'
        assert results['11-108218120-CT-C']['g_hgvs'] == 'NC_000011.10:g.108218122del'
        assert results['11-108218120-CT-C']['genomic_variant_error'] is None
        assert 'NM_002519.2' in results['11-108218120-CT-C']['hgvs_t_and_p'].keys()
        assert results['11-108218120-CT-C']['hgvs_t_and_p']['NM_002519.2']['t_hgvs'] == 'NM_002519.2:c.37+4379del'
        assert results['11-108218120-CT-C']['hgvs_t_and_p']['NM_002519.2']['p_hgvs_tlc'] == 'NP_002510.2:p.?'
        assert results['11-108218120-CT-C']['hgvs_t_and_p']['NM_002519.2']['p_hgvs_slc'] == 'NP_002510.2:p.?'
        assert results['11-108218120-CT-C']['hgvs_t_and_p']['NM_002519.2']['transcript_variant_error'] is None
        assert 'NM_001321307.1' in results['11-108218120-CT-C']['hgvs_t_and_p'].keys()
        assert results['11-108218120-CT-C']['hgvs_t_and_p']['NM_001321307.1']['t_hgvs'] == 'NM_001321307.1:c.37+4379del'
        assert results['11-108218120-CT-C']['hgvs_t_and_p']['NM_001321307.1']['p_hgvs_tlc'] == 'NP_001308236.1:p.?'
        assert results['11-108218120-CT-C']['hgvs_t_and_p']['NM_001321307.1']['p_hgvs_slc'] == 'NP_001308236.1:p.?'
        assert results['11-108218120-CT-C']['hgvs_t_and_p']['NM_001321307.1']['transcript_variant_error'] is None

    def test_variant241(self):
        variant = 'X-76813150-T-ATA'
        results = vf.FormatVariant(variant, 'GRCh38', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'X-76813150-T-ATA' in results.keys()
        assert results['X-76813150-T-ATA']['p_vcf'] == 'X-76813150-T-ATA'
        assert results['X-76813150-T-ATA']['g_hgvs'] == 'NC_000023.11:g.76813150delinsATA'
        assert results['X-76813150-T-ATA']['genomic_variant_error'] is None
        assert 'NR_110406.2' in results['X-76813150-T-ATA']['hgvs_t_and_p'].keys()
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110406.2']['t_hgvs'] == 'NR_110406.2:n.412-154842delinsTAT'
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110406.2']['p_hgvs_tlc'] is None
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110406.2']['p_hgvs_slc'] is None
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110406.2']['transcript_variant_error'] is None
        assert 'NR_110402.2' in results['X-76813150-T-ATA']['hgvs_t_and_p'].keys()
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110402.2']['t_hgvs'] == 'NR_110402.2:n.412-120063delinsTAT'
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110402.2']['p_hgvs_tlc'] is None
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110402.2']['p_hgvs_slc'] is None
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110402.2']['transcript_variant_error'] is None
        assert 'NR_110401.2' in results['X-76813150-T-ATA']['hgvs_t_and_p'].keys()
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110401.2']['t_hgvs'] == 'NR_110401.2:n.412-154482delinsTAT'
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110401.2']['p_hgvs_tlc'] is None
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110401.2']['p_hgvs_slc'] is None
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110401.2']['transcript_variant_error'] is None
        assert 'NR_110404.2' in results['X-76813150-T-ATA']['hgvs_t_and_p'].keys()
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110404.2']['t_hgvs'] == 'NR_110404.2:n.307-154482delinsTAT'
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110404.2']['p_hgvs_tlc'] is None
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110404.2']['p_hgvs_slc'] is None
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110404.2']['transcript_variant_error'] is None
        assert 'NR_110400.2' in results['X-76813150-T-ATA']['hgvs_t_and_p'].keys()
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110400.2']['t_hgvs'] == 'NR_110400.2:n.307-120063delinsTAT'
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110400.2']['p_hgvs_tlc'] is None
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110400.2']['p_hgvs_slc'] is None
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110400.2']['transcript_variant_error'] is None
        assert 'NR_110405.2' in results['X-76813150-T-ATA']['hgvs_t_and_p'].keys()
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110405.2']['t_hgvs'] == 'NR_110405.2:n.307-31858delinsTAT'
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110405.2']['p_hgvs_tlc'] is None
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110405.2']['p_hgvs_slc'] is None
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110405.2']['transcript_variant_error'] is None
        assert 'NR_110403.2' in results['X-76813150-T-ATA']['hgvs_t_and_p'].keys()
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110403.2']['t_hgvs'] == 'NR_110403.2:n.389+24660delinsTAT'
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110403.2']['p_hgvs_tlc'] is None
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110403.2']['p_hgvs_slc'] is None
        assert results['X-76813150-T-ATA']['hgvs_t_and_p']['NR_110403.2']['transcript_variant_error'] is None

    def test_issue_322a(self):
        variant = 'NC_000017.10:g.48275363CC>A'
        results = vf.FormatVariant(variant, 'GRCh37', vfo,  'refseq', None)
        results = results.stucture_data()
        print(results)
        assert 'NC_000017.10:g.48275363CC>A' in results.keys()
        assert results['NC_000017.10:g.48275363CC>A']['g_hgvs'] == 'NC_000017.10:g.48275363_48275364delinsA'

    def test_issue_360a(self):
        results = VariantFormatter.simpleVariantFormatter.format('NC_012920.1:m.1011C>T',
                                                                 'GRCh38', 'refseq', None, False, True)
        print(results)
        assert 'NC_012920.1:m.1011C>T' in results.keys()
        assert "NC_001807.4" in results['NC_012920.1:m.1011C>T']['NC_012920.1:m.1011C>T']['hgvs_t_and_p'][
            'intergenic']['primary_assembly_loci']['hg19'].keys()
        assert "NC_012920.1" in results['NC_012920.1:m.1011C>T']['NC_012920.1:m.1011C>T']['hgvs_t_and_p'][
            'intergenic']['primary_assembly_loci']['hg38'].keys()
        assert "NC_012920.1" in results['NC_012920.1:m.1011C>T']['NC_012920.1:m.1011C>T']['hgvs_t_and_p'][
            'intergenic']['primary_assembly_loci']['grch37'].keys()
        assert "NC_012920.1" in results['NC_012920.1:m.1011C>T']['NC_012920.1:m.1011C>T']['hgvs_t_and_p'][
            'intergenic']['primary_assembly_loci']['grch38'].keys()

    def test_issue_360b(self):
        results = VariantFormatter.simpleVariantFormatter.format('NC_001807.4:m.1013C>T',
                                                                 'hg19', 'refseq', None, False, True)
        print(results)
        assert 'NC_001807.4:m.1013C>T' in results.keys()
        assert "NC_001807.4" in results['NC_001807.4:m.1013C>T']['NC_001807.4:m.1013C>T']['hgvs_t_and_p'][
            'intergenic']['primary_assembly_loci']['hg19'].keys()
        assert "NC_012920.1" in results['NC_001807.4:m.1013C>T']['NC_001807.4:m.1013C>T']['hgvs_t_and_p'][
            'intergenic']['primary_assembly_loci']['hg38'].keys()
        assert "NC_012920.1" in results['NC_001807.4:m.1013C>T']['NC_001807.4:m.1013C>T']['hgvs_t_and_p'][
            'intergenic']['primary_assembly_loci']['grch37'].keys()
        assert "NC_012920.1" in results['NC_001807.4:m.1013C>T']['NC_001807.4:m.1013C>T']['hgvs_t_and_p'][
            'intergenic']['primary_assembly_loci']['grch38'].keys()

    # test select_transcripts 'all' vs 'raw'
    def test_issue_392a(self):
        vfo.testing = False
        results = VariantFormatter.simpleVariantFormatter.format('NC_000008.10:g.6673379del',
                                                                 'GRCh37', 'refseq', 'all', False, True, testing=False)
        print(results)
        assert 'NC_000008.10:g.6673379del' in results.keys()
        assert 'NC_000008.10:g.6673379del' in results['NC_000008.10:g.6673379del'][
            'NC_000008.10:g.6673379del']['g_hgvs']
        assert 'NM_207411.4' not in results['NC_000008.10:g.6673379del']['NC_000008.10:g.6673379del']\
            ['hgvs_t_and_p'].keys()
        assert 'NM_001289973.1' not in results['NC_000008.10:g.6673379del']['NC_000008.10:g.6673379del']\
            ['hgvs_t_and_p'].keys()

        results = VariantFormatter.simpleVariantFormatter.format('NC_000008.10:g.6673379del',
                                                                 'GRCh37', 'refseq', 'raw', False, True, testing=False)
        print(results)
        assert 'NC_000008.10:g.6673379del' in results.keys()
        assert 'NC_000008.10:g.6673379del' in results['NC_000008.10:g.6673379del'][
            'NC_000008.10:g.6673379del']['g_hgvs']
        assert 'NM_207411.4' in results['NC_000008.10:g.6673379del']['NC_000008.10:g.6673379del']\
            ['hgvs_t_and_p'].keys()
        assert 'NM_001289973.1' in results['NC_000008.10:g.6673379del']['NC_000008.10:g.6673379del']\
            ['hgvs_t_and_p'].keys()

        vfo.testing = True

    def test_issue_392b(self):
        results = VariantFormatter.simpleVariantFormatter.format('NC_000008.11:g.6815857G>A',
                                                                 'GRCh38', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000008.11:g.6815857G>A' in results.keys()
        assert 'NC_000008.11:g.6815857G>A' in results['NC_000008.11:g.6815857G>A'][
            'NC_000008.11:g.6815857G>A']['g_hgvs']
        assert 'NC_000008.10:g.6673379del' in results['NC_000008.11:g.6815857G>A'][
            'NC_000008.11:g.6815857G>A']['hgvs_t_and_p']['NM_001289973.1']['primary_assembly_loci'][
            'grch37']['NC_000008.10']['hgvs_genomic_description']

    def test_issue_370(self):
        results = VariantFormatter.simpleVariantFormatter.format('NC_000008.10:g.24811072C>T',
                                                                 'GRCh37', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000008.10:g.24811072C>T' in results.keys()
        assert 'NM_006158.3:c.1407delinsAC' in results['NC_000008.10:g.24811072C>T'][
            'NC_000008.10:g.24811072C>T']['hgvs_t_and_p']['NM_006158.3']['t_hgvs']

    def test_issue_370(self):
        results = VariantFormatter.simpleVariantFormatter.format('NC_000008.10:g.24811072C>T',
                                                                 'GRCh37', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000008.10:g.24811072C>T' in results.keys()
        assert 'NM_006158.3:c.1407delinsAC' in results['NC_000008.10:g.24811072C>T'][
            'NC_000008.10:g.24811072C>T']['hgvs_t_and_p']['NM_006158.3']['t_hgvs']

    def test_multiple_variants(self):
        results = VariantFormatter.simpleVariantFormatter.format('["NC_000008.10:g.24811072C>T", '
                                                                 '"NC_000008.10:g.24811072C>A"]',
                                                                 'GRCh37', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000008.10:g.24811072C>T' in results.keys()
        assert 'NC_000008.10:g.24811072C>A' in results.keys()

# <LICENSE>
# Copyright (C) 2016-2024 VariantValidator Contributors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>