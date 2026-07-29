from unittest import TestCase

from VariantFormatter import simpleVariantFormatter


class TestWarnings(TestCase):

    def test_issue_360(self):
        # Mitochondrial references must use m., not g.
        result = simpleVariantFormatter.format(
            'NC_012920.1:g.100del',
            'GRCh37',
            'refseq',
            None,
            False,
            True
        )
        assert (
                "The given reference sequence (NC_012920.1) does not match the DNA type (g). "
                "For NC_012920.1, please use (m). For g. variants, please use a linear genomic reference sequence"
                in result["NC_012920.1:g.100del"]["NC_012920.1:g.100del"]["genomic_variant_error"]
        )

        # The DNA type error is more fundamental than the genome-build error, so it
        # should be reported first regardless of the supplied genome build.
        result = simpleVariantFormatter.format(
            'NC_012920.1:g.100del',
            'hg19',
            'refseq',
            None,
            False,
            True
        )
        assert (
                "does not match the DNA type (g)" in
                result["NC_012920.1:g.100del"]["NC_012920.1:g.100del"]["genomic_variant_error"]
        )

        result = simpleVariantFormatter.format(
            'NC_012920.1:m.1011C>T',
            'GRCh37',
            'refseq',
            None,
            False,
            True
        )
        assert "grch37" in result["NC_012920.1:m.1011C>T"]["NC_012920.1:m.1011C>T"][
            "hgvs_t_and_p"]["intergenic"]["primary_assembly_loci"].keys()
        assert "grch38" in result["NC_012920.1:m.1011C>T"]["NC_012920.1:m.1011C>T"][
            "hgvs_t_and_p"]["intergenic"]["primary_assembly_loci"].keys()
        assert "hg38" in result["NC_012920.1:m.1011C>T"]["NC_012920.1:m.1011C>T"][
            "hgvs_t_and_p"]["intergenic"]["primary_assembly_loci"].keys()
        assert "hg19" not in result["NC_012920.1:m.1011C>T"]["NC_012920.1:m.1011C>T"][
            "hgvs_t_and_p"]["intergenic"]["primary_assembly_loci"].keys()

        result = simpleVariantFormatter.format(
            'NC_012920.1:m.1011C>T',
            'GRCh38',
            'refseq',
            None,
            False,
            True
        )
        assert "grch37" in result["NC_012920.1:m.1011C>T"]["NC_012920.1:m.1011C>T"][
            "hgvs_t_and_p"]["intergenic"]["primary_assembly_loci"].keys()
        assert "grch38" in result["NC_012920.1:m.1011C>T"]["NC_012920.1:m.1011C>T"][
            "hgvs_t_and_p"]["intergenic"]["primary_assembly_loci"].keys()
        assert "hg38" in result["NC_012920.1:m.1011C>T"]["NC_012920.1:m.1011C>T"][
            "hgvs_t_and_p"]["intergenic"]["primary_assembly_loci"].keys()
        assert "hg19" in result["NC_012920.1:m.1011C>T"]["NC_012920.1:m.1011C>T"][
            "hgvs_t_and_p"]["intergenic"]["primary_assembly_loci"].keys()

        result = simpleVariantFormatter.format(
            'NC_012920.1:m.1011C>T',
            'hg38',
            'refseq',
            None,
            False,
            True
        )
        assert "grch37" in result["NC_012920.1:m.1011C>T"]["NC_012920.1:m.1011C>T"][
            "hgvs_t_and_p"]["intergenic"]["primary_assembly_loci"].keys()
        assert "grch38" in result["NC_012920.1:m.1011C>T"]["NC_012920.1:m.1011C>T"][
            "hgvs_t_and_p"]["intergenic"]["primary_assembly_loci"].keys()
        assert "hg38" in result["NC_012920.1:m.1011C>T"]["NC_012920.1:m.1011C>T"][
            "hgvs_t_and_p"]["intergenic"]["primary_assembly_loci"].keys()
        assert "hg19" in result["NC_012920.1:m.1011C>T"]["NC_012920.1:m.1011C>T"][
            "hgvs_t_and_p"]["intergenic"]["primary_assembly_loci"].keys()


    def test_issue_322b(self):
        results = simpleVariantFormatter.format('NC_000017.10:g.48275363CG>A',
                                                'GRCh37', 'refseq', None, False, True)

        print(results)
        assert 'NC_000017.10:g.48275363CG>A' in results.keys()
        assert "Variant reference (CG) does not agree with reference sequence (CC)" in results[
            'NC_000017.10:g.48275363CG>A']['NC_000017.10:g.48275363CG>A']['genomic_variant_error']


class TestVFGapWarnings(TestCase):

    def test_vf_series_1(self):
        results = simpleVariantFormatter.format('NC_000004.11:g.140811117C>A', 'GRCh37', 'refseq', None, False, True,
                                                testing=True)
        print(results)
        assert 'NC_000004.11:g.140811117C>A' in results.keys()
        assert 'NM_018717.4 contains 3 fewer bases between c.2276_2277, and 12 fewer bases between c.1467_1468 ' \
               'than NC_000004.11' in results['NC_000004.11:g.140811117C>A'][
            'NC_000004.11:g.140811117C>A']['hgvs_t_and_p']['NM_018717.4']['gap_statement']


    def test_vf_series_2(self):
        results = simpleVariantFormatter.format('NC_000008.10:g.24811072C>T',
                                                                 'GRCh37', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000008.10:g.24811072C>T' in results.keys()
        assert 'NM_006158.3 contains 1 fewer bases between c.1407_1408 than NC_000008.10' in results[
            'NC_000008.10:g.24811072C>T']['NC_000008.10:g.24811072C>T']['hgvs_t_and_p'][
            'NM_006158.3']['gap_statement']
        assert 'NM_006158.4 contains 1 fewer bases between c.1407_1408 than NC_000008.10' in results[
            'NC_000008.10:g.24811072C>T']['NC_000008.10:g.24811072C>T']['hgvs_t_and_p'][
            'NM_006158.4']['gap_statement']
        assert 'NM_006158.5 contains 1 fewer bases between c.1413_1414 than NC_000008.10' in results[
            'NC_000008.10:g.24811072C>T']['NC_000008.10:g.24811072C>T']['hgvs_t_and_p'][
            'NM_006158.5']['gap_statement']


    def test_vf_series_3(self):
        results = simpleVariantFormatter.format('NC_000015.9:g.72105933del',
                                                                 'GRCh37', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000015.9:g.72105933del' in results.keys()
        assert 'NM_014249.2 contains 1 fewer bases between c.947_948 than NC_000015.9' in results[
            'NC_000015.9:g.72105933del']['NC_000015.9:g.72105933del']['hgvs_t_and_p'][
            'NM_014249.2']['gap_statement']
        assert 'NM_014249.3 contains 1 fewer bases between c.947_948 than NC_000015.9' in results[
            'NC_000015.9:g.72105933del']['NC_000015.9:g.72105933del']['hgvs_t_and_p'][
            'NM_014249.3']['gap_statement']
        assert 'NM_014249.4 contains 1 fewer bases between c.951_952 than NC_000015.9' in results[
            'NC_000015.9:g.72105933del']['NC_000015.9:g.72105933del']['hgvs_t_and_p'][
            'NM_014249.4']['gap_statement']
        assert 'NM_016346.2 contains 1 fewer bases between c.947_948 than NC_000015.9' in results[
            'NC_000015.9:g.72105933del']['NC_000015.9:g.72105933del']['hgvs_t_and_p'][
            'NM_016346.2']['gap_statement']
        assert 'NM_016346.3 contains 1 fewer bases between c.947_948 than NC_000015.9' in results[
            'NC_000015.9:g.72105933del']['NC_000015.9:g.72105933del']['hgvs_t_and_p'][
            'NM_016346.3']['gap_statement']
        assert 'NM_016346.4 contains 1 fewer bases between c.951_952 than NC_000015.9' in results[
            'NC_000015.9:g.72105933del']['NC_000015.9:g.72105933del']['hgvs_t_and_p'][
            'NM_016346.4']['gap_statement']


    def test_vf_series_4(self):
        results = simpleVariantFormatter.format('NC_000019.9:g.41123095dup',
                                                                 'GRCh37', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000019.9:g.41123095dup' in results.keys()
        assert 'NM_001042544.1 contains 1 extra bases between c.3233_3235 than NC_000019.9' in results[
            'NC_000019.9:g.41123095dup']['NC_000019.9:g.41123095dup']['hgvs_t_and_p'][
            'NM_001042544.1']['gap_statement']
        assert 'NM_001042545.1 contains 1 extra bases between c.3032_3034 than NC_000019.9' in results[
            'NC_000019.9:g.41123095dup']['NC_000019.9:g.41123095dup']['hgvs_t_and_p'][
            'NM_001042545.1']['gap_statement']
        assert 'NM_001042545.2 contains 1 extra bases between c.3034_3036 than NC_000019.9' in results[
            'NC_000019.9:g.41123095dup']['NC_000019.9:g.41123095dup']['hgvs_t_and_p'][
            'NM_001042545.2']['gap_statement']
        assert 'NM_003573.2 contains 1 extra bases between c.3122_3124 than NC_000019.9' in results[
            'NC_000019.9:g.41123095dup']['NC_000019.9:g.41123095dup']['hgvs_t_and_p'][
            'NM_003573.2']['gap_statement']


    def test_vf_series_5(self):
        results = simpleVariantFormatter.format('NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=',
                                                                 'GRCh37', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=' in results.keys()
        assert 'GappedAlignmentWarning: NM_001083585.1 contains 25 fewer bases between c.*344_*345 than NC_000017.10' in results[
            'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
            'NM_001083585.1']['gap_statement']
        assert 'GappedAlignmentWarning: NM_001083585.2 contains 25 fewer bases between c.*344_*345 than NC_000017.10' in results[
            'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
            'NM_001083585.2']['gap_statement']
        assert 'GappedAlignmentWarning: NM_001083585.3 contains 25 fewer bases between c.*369_*370 than NC_000017.10' in results[
            'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
            'NM_001083585.3']['gap_statement']
        assert 'GappedAlignmentWarning: NM_001291581.1 contains 25 fewer bases between c.*344_*345 than NC_000017.10' in results[
            'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
            'NM_001291581.1']['gap_statement']
        assert 'GappedAlignmentWarning: NM_001291581.2 contains 25 fewer bases between c.*369_*370 than NC_000017.10' in results[
            'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
            'NM_001291581.2']['gap_statement']
        assert 'GappedAlignmentWarning: NM_004703.4 contains 25 fewer bases between c.*344_*345 than NC_000017.10' in results[
            'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
            'NM_004703.4']['gap_statement']
        assert 'GappedAlignmentWarning: NM_004703.5 contains 25 fewer bases between c.*344_*345 than NC_000017.10' in results[
            'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
            'NM_004703.5']['gap_statement']
        assert 'GappedAlignmentWarning: NM_004703.6 contains 25 fewer bases between c.*369_*370 than NC_000017.10' in results[
            'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
            'NM_004703.6']['gap_statement']


    def test_vf_series_6(self):
        results = simpleVariantFormatter.format('NC_000012.11:g.122064777C>A',
                                                                 'GRCh37', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000012.11:g.122064777C>A' in results.keys()
        assert 'GappedAlignmentWarning: NM_032790.3 contains 6 fewer bases between c.126_127 than NC_000012.11' in results[
            'NC_000012.11:g.122064777C>A']['NC_000012.11:g.122064777C>A']['hgvs_t_and_p'][
            'NM_032790.3']['gap_statement']


    def test_vf_series_6a_regression(self):
        results = simpleVariantFormatter.format('NC_000012.11:g.122064700del',
                                                                 'GRCh37', 'refseq', None, False, True, testing=True)
        print(results)
        assert results['NC_000012.11:g.122064700del']["NC_000012.11:g.122064700del"]['hgvs_t_and_p'][
            'NM_032790.4']['gap_statement'] == "GappedAlignmentWarning: NM_032790.4 contains 6 fewer bases between c.137_138 than NC_000012.11"
        assert results['NC_000012.11:g.122064700del']["NC_000012.11:g.122064700del"]['hgvs_t_and_p'][
            'NM_032790.4']['gapped_alignment_warning'] == "GappedAlignmentWarning: Variation described in the context of an imperfect alignment of NM_032790.4 with NC_000012.11 (genome build GRCh37)"


    def test_vf_series_7(self):
        results = simpleVariantFormatter.format('NC_000002.11:g.95847041_95847043GCG=',
                                                                 'GRCh37', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000002.11:g.95847041_95847043GCG=' in results.keys()
        assert 'GappedAlignmentWarning: NM_001017396.1 contains 3 fewer bases between c.341_342 than NC_000002.11' in results[
            'NC_000002.11:g.95847041_95847043GCG=']['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'][
            'NM_001017396.1']['gap_statement']
        assert 'GappedAlignmentWarning: NM_001017396.2 contains 3 fewer bases between c.341_342 than NC_000002.11' in results[
            'NC_000002.11:g.95847041_95847043GCG=']['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'][
            'NM_001017396.2']['gap_statement']
        assert 'GappedAlignmentWarning: NM_001282398.1 contains 3 fewer bases between c.353_354 than NC_000002.11' in results[
            'NC_000002.11:g.95847041_95847043GCG=']['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'][
            'NM_001282398.1']['gap_statement']
        assert 'GappedAlignmentWarning: NM_001291604.1 contains 3 fewer bases between c.227_228 than NC_000002.11' in results[
            'NC_000002.11:g.95847041_95847043GCG=']['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'][
            'NM_001291604.1']['gap_statement']
        assert 'GappedAlignmentWarning: NM_001291605.1 contains 3 fewer bases between c.506_507 than NC_000002.11' in results[
            'NC_000002.11:g.95847041_95847043GCG=']['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'][
            'NM_001291605.1']['gap_statement']
        assert 'GappedAlignmentWarning: NM_021088.2 contains 3 fewer bases between c.467_468 than NC_000002.11' in results[
            'NC_000002.11:g.95847041_95847043GCG=']['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'][
            'NM_021088.2']['gap_statement']
        assert 'GappedAlignmentWarning: NM_021088.3 contains 3 fewer bases between c.467_468 than NC_000002.11' in results[
            'NC_000002.11:g.95847041_95847043GCG=']['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'][
            'NM_021088.3']['gap_statement']


    def test_vf_series_8(self):
        results = simpleVariantFormatter.format('NC_000003.11:g.14561629_14561630insG',
                                                'GRCh37', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'GappedAlignmentWarning: NM_001080423.2 contains 1 extra bases between c.1308_1310 than NC_000003.11' in results[
            'NC_000003.11:g.14561629_14561630insG']['NC_000003.11:g.14561629_14561630insG']['hgvs_t_and_p'][
            'NM_001080423.2']['gap_statement']
        assert 'GappedAlignmentWarning: NM_001080423.3 contains 1 extra bases between c.1017_1019 than NC_000003.11' in results[
            'NC_000003.11:g.14561629_14561630insG']['NC_000003.11:g.14561629_14561630insG']['hgvs_t_and_p'][
            'NM_001080423.3']['gap_statement']
        assert 'GappedAlignmentWarning: NM_001080423.4 contains 1 extra bases between c.1019_1021 than NC_000003.11' in results[
            'NC_000003.11:g.14561629_14561630insG']['NC_000003.11:g.14561629_14561630insG']['hgvs_t_and_p'][
            'NM_001080423.4']['gap_statement']


    def test_vf_series_9(self):
        results = simpleVariantFormatter.format('NC_000004.11:g.140811117C>A',
                                                'GRCh37', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000004.11:g.140811117C>A' in results.keys()
        assert 'GappedAlignmentWarning: NM_018717.4 contains 3 fewer bases between c.2276_2277, and 12 fewer bases between c.1467_1468 ' \
               'than NC_000004.11' in results[
            'NC_000004.11:g.140811117C>A']['NC_000004.11:g.140811117C>A']['hgvs_t_and_p'][
            'NM_018717.4']['gap_statement']


    def test_vf_series_10(self):
        results = simpleVariantFormatter.format('NC_000009.11:g.136132908_136132909TA=',
                                                'GRCh37', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000009.11:g.136132908_136132909TA=' in results.keys()
        assert 'GappedAlignmentWarning: NM_020469.2 contains 1 extra bases between c.260_262 than NC_000009.11' in results[
            'NC_000009.11:g.136132908_136132909TA=']['NC_000009.11:g.136132908_136132909TA=']['hgvs_t_and_p'][
            'NM_020469.2']['gap_statement']
        assert 'GappedAlignmentWarning: NM_020469.3 contains 22 extra bases between c.*756_*757, and 2 extra bases between c.*797_*798, ' \
               'and 110 extra bases between c.*840_*841, and 2 extra bases between c.*4648_*4649, and 1 extra ' \
               'bases between c.260_262 than NC_000009.11' in results[
            'NC_000009.11:g.136132908_136132909TA=']['NC_000009.11:g.136132908_136132909TA=']['hgvs_t_and_p'][
            'NM_020469.3']['gap_statement']


    def test_vf_series_11(self):
        results = simpleVariantFormatter.format('NC_000019.10:g.50378563_50378564insTAC',
                                                'GRCh38', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000019.10:g.50378563_50378564insTAC' in results.keys()
        assert 'NM_001256647.1 contains 3 extra bases between c.223_227 than NC_000019.10' in results[
            'NC_000019.10:g.50378563_50378564insTAC']['NC_000019.10:g.50378563_50378564insTAC']['hgvs_t_and_p'][
            'NM_001256647.1']['gap_statement']
        assert 'NM_007121.5 contains 3 extra bases between c.514_518 than NC_000019.10' in results[
            'NC_000019.10:g.50378563_50378564insTAC']['NC_000019.10:g.50378563_50378564insTAC']['hgvs_t_and_p'][
            'NM_007121.5']['gap_statement']


    def test_vf_series_12(self):
        results = simpleVariantFormatter.format('NC_000007.13:g.149476664_149476666delinsTC',
                                                'GRCh37', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000007.13:g.149476664_149476666delinsTC' in results.keys()
        assert 'NR_163594.1 contains 1 extra bases between n.1129_1131, and 1 fewer bases between n.11675_11676 ' \
               'than NC_000007.13' in results[
            'NC_000007.13:g.149476664_149476666delinsTC']['NC_000007.13:g.149476664_149476666delinsTC'][
            'hgvs_t_and_p']['NR_163594.1']['gap_statement']


    def test_vf_series_13(self):
        results = simpleVariantFormatter.format('NC_000004.12:g.139889957_139889968del',
                                                'GRCh38', 'refseq', None, False, True, testing=True)
        print(results)
        assert 'NC_000004.12:g.139889957_139889968del' in results.keys()
        assert 'NM_018717.4 contains 3 fewer bases between c.2276_2277, and 12 fewer bases between c.1467_1468 ' \
               'than NC_000004.12' in results[
            'NC_000004.12:g.139889957_139889968del']['NC_000004.12:g.139889957_139889968del']['hgvs_t_and_p'][
            'NM_018717.4']['gap_statement']


    def test_aligned_transcript_versions_vf(self):
        results = simpleVariantFormatter.format('NC_000009.12:g.134642190_134642191inv',
                                                                 'GRCh38', 'all', "raw", False, False, testing=True)
        print(results)
        assert 'NC_000009.12:g.134642190_134642191inv' in results.keys()
        assert 'TranscriptVersionWarning: A more recent version of the selected reference sequence NM_000093.4 is available for genome build GRCh38 (NM_000093.5)' in results[
            'NC_000009.12:g.134642190_134642191inv']['NC_000009.12:g.134642190_134642191inv']['hgvs_t_and_p'][
            'NM_000093.4']['transcript_version_warning']
        assert 'TranscriptVersionWarning: A more recent version of the selected reference sequence NM_000093.3 is available for genome build GRCh38 (NM_000093.5)' in results[
            'NC_000009.12:g.134642190_134642191inv']['NC_000009.12:g.134642190_134642191inv']['hgvs_t_and_p'][
            'NM_000093.3']['transcript_version_warning']

        results = simpleVariantFormatter.format('NC_000008.11:g.10623201T>A',
                                                                 'GRCh38', 'all', "raw", False, False, testing=True)

        assert 'NC_000008.11:g.10623201T>A' in results.keys()
        assert 'TranscriptVersionWarning: A more recent version of the selected reference sequence ENST00000382483.3 is available for genome build GRCh38 (ENST00000382483.4)' in results[
            'NC_000008.11:g.10623201T>A']['NC_000008.11:g.10623201T>A']['hgvs_t_and_p'][
            'ENST00000382483.3']['transcript_version_warning']

    def test_vf_series_11_wrong_genome_build(self):
        variant = 'NC_000019.10:g.50378563_50378564insTAC'

        results = simpleVariantFormatter.format(
            variant,
            'GRCh37',
            'refseq',
            None,
            False,
            True,
            testing=True
        )

        print(results)

        assert variant in results.keys()

        entry = results[variant][variant]

        assert (
            "GenomeBuildError: chromosome ID NC_000019.10 is not associated "
            "with genome build GRCh37"
            in entry["genomic_variant_error"]
        )


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
