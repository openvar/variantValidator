import pytest
from VariantValidator import variantValidator as vv


class TestVariants(object):

    @classmethod
    def setup_class(cls):
        vv.my_config()

    def test_variant1(self):
        """Test for first variant in file, all other tests in this class will follow the same format."""
        variant = "NM_015120.4:c.35T>C"
        results = vv.validator(variant, 'GRCh37', 'all')
        print results

        assert variant in results.keys()
        assert results[variant]['hgvs_lrg_transcript_variant'] == 'LRG_741t1:c.35T>C'
        assert results[variant]['refseqgene_context_intronic_sequence'] == ''
        assert results[variant]['alt_genomic_loci'] == []
        assert results[variant]['transcript_description'] == 'Homo sapiens ALMS1, centrosome and basal body associated protein (ALMS1), mRNA'
        assert results[variant]['gene_symbol'] == 'ALMS1'
        assert results[variant]['hgvs_predicted_protein_consequence']['tlr'] == 'NP_055935.4(LRG_741p1):p.(Leu12Pro)'
        assert results[variant]['hgvs_predicted_protein_consequence']['slr'] == 'NP_055935.4:p.(L12P)'
        assert results[variant]['genome_context_intronic_sequence'] == ''
        assert results[variant]['hgvs_lrg_variant'] == 'LRG_741:g.5146T>C'
        assert results[variant]['hgvs_transcript_variant'] == 'NM_015120.4:c.35T>C'
        assert results[variant]['hgvs_refseqgene_variant'] == 'NG_011690.1:g.5146T>C'
        assert results[variant]['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000002.11:g.73613031delinsCGGA', 'vcf': {'chr': 'chr2', 'ref': 'T', 'pos': '73613031', 'alt': 'CGGA'}}
        assert results[variant]['primary_assembly_loci']['hg38'] == {'hgvs_genomic_description': 'NC_000002.12:g.73385903delinsCGGA', 'vcf': {'chr': 'chr2', 'ref': 'T', 'pos': '73385903', 'alt': 'CGGA'}}
        assert results[variant]['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000002.11:g.73613031delinsCGGA', 'vcf': {'chr': '2', 'ref': 'T', 'pos': '73613031', 'alt': 'CGGA'}}
        assert results[variant]['primary_assembly_loci']['grch38'] == {'hgvs_genomic_description': 'NC_000002.12:g.73385903delinsCGGA', 'vcf': {'chr': '2', 'ref': 'T', 'pos': '73385903', 'alt': 'CGGA'}}
        assert results[variant]['reference_sequence_records']['protein'] == 'https://www.ncbi.nlm.nih.gov/nuccore/NP_055935.4'
        assert results[variant]['reference_sequence_records']['transcript'] == 'https://www.ncbi.nlm.nih.gov/nuccore/NM_015120.4'
        assert results['flag'] == 'gene_variant'

    def test_variant2(self):
        variant = "NM_015120.4:c.39G>C"
        results = vv.validator(variant, 'GRCh37', 'all')
        print results

        assert variant in results.keys()
        assert results[variant]['hgvs_lrg_transcript_variant'] == 'LRG_741t1:c.39G>C'
        assert results[variant]['refseqgene_context_intronic_sequence'] == ''
        assert results[variant]['alt_genomic_loci'] == []
        assert results[variant]['transcript_description'] == 'Homo sapiens ALMS1, centrosome and basal body associated protein (ALMS1), mRNA'
        assert results[variant]['gene_symbol'] == 'ALMS1'
        assert results[variant]['hgvs_predicted_protein_consequence']['tlr'] == 'NP_055935.4(LRG_741p1):p.(Glu13Asp)'
        assert results[variant]['hgvs_predicted_protein_consequence']['slr'] == 'NP_055935.4:p.(E13D)'
        assert results[variant]['genome_context_intronic_sequence'] == ''
        assert results[variant]['hgvs_lrg_variant'] == 'LRG_741:g.5150G>C'
        assert results[variant]['hgvs_transcript_variant'] == 'NM_015120.4:c.39G>C'
        assert results[variant]['hgvs_refseqgene_variant'] == 'NG_011690.1:g.5150G>C'
        assert results[variant]['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000002.11:g.73613034_73613035insCGA', 'vcf': {'chr': 'chr2', 'ref': 'G', 'pos': '73613032', 'alt': 'GGAC'}}
        assert results[variant]['primary_assembly_loci']['hg38'] == {'hgvs_genomic_description': 'NC_000002.12:g.73385906_73385907insCGA', 'vcf': {'chr': 'chr2', 'ref': 'G', 'pos': '73385904', 'alt': 'GGAC'}}
        assert results[variant]['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000002.11:g.73613034_73613035insCGA', 'vcf': {'chr': '2', 'ref': 'G', 'pos': '73613032', 'alt': 'GGAC'}}
        assert results[variant]['primary_assembly_loci']['grch38'] == {'hgvs_genomic_description': 'NC_000002.12:g.73385906_73385907insCGA', 'vcf': {'chr': '2', 'ref': 'G', 'pos': '73385904', 'alt': 'GGAC'}}
        assert results['flag'] == 'gene_variant'

    def test_variant3(self):
        variant = "NM_015120.4:c.34C>T"
        results = vv.validator(variant, 'GRCh37', 'all')
        print results

        assert variant in results.keys()
        assert results[variant]['hgvs_lrg_transcript_variant'] == 'LRG_741t1:c.34C>T'
        assert results[variant]['refseqgene_context_intronic_sequence'] == ''
        assert results[variant]['alt_genomic_loci'] == []
        assert results[variant]['transcript_description'] == 'Homo sapiens ALMS1, centrosome and basal body associated protein (ALMS1), mRNA'
        assert results[variant]['gene_symbol'] == 'ALMS1'
        assert results[variant]['hgvs_predicted_protein_consequence']['tlr'] == 'NP_055935.4(LRG_741p1):p.(Leu12=)'
        assert results[variant]['hgvs_predicted_protein_consequence']['slr'] == 'NP_055935.4:p.(L12=)'
        assert results[variant]['genome_context_intronic_sequence'] == ''
        assert results[variant]['hgvs_lrg_variant'] == 'LRG_741:g.5145C>T'
        assert results[variant]['hgvs_transcript_variant'] == 'NM_015120.4:c.34C>T'
        assert results[variant]['hgvs_refseqgene_variant'] == 'NG_011690.1:g.5145C>T'
        assert results[variant]['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000002.11:g.73613030C>T', 'vcf': {'chr': 'chr2', 'ref': u'C', 'pos': '73613030', 'alt': u'T'}}
        assert results[variant]['primary_assembly_loci']['hg38'] == {'hgvs_genomic_description': 'NC_000002.12:g.73385902C>T', 'vcf': {'chr': 'chr2', 'ref': u'C', 'pos': '73385902', 'alt': u'T'}}
        assert results[variant]['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000002.11:g.73613030C>T', 'vcf': {'chr': '2', 'ref': u'C', 'pos': '73613030', 'alt': u'T'}}
        assert results[variant]['primary_assembly_loci']['grch38'] == {'hgvs_genomic_description': 'NC_000002.12:g.73385902C>T', 'vcf': {'chr': '2', 'ref': u'C', 'pos': '73385902', 'alt': u'T'}}
        assert results['flag'] == 'gene_variant'

    def test_variant4(self):
        variant = "NC_000002.11:g.73613030C>T"
        results = vv.validator(variant, 'GRCh37', 'all')
        print results

        assert 'NM_015120.4:c.34C>T' in results.keys()
        assert results['NM_015120.4:c.34C>T']['hgvs_lrg_transcript_variant'] == 'LRG_741t1:c.34C>T'
        assert results['NM_015120.4:c.34C>T']['refseqgene_context_intronic_sequence'] == ''
        assert results['NM_015120.4:c.34C>T']['alt_genomic_loci'] == []
        assert results['NM_015120.4:c.34C>T']['transcript_description'] == 'Homo sapiens ALMS1, centrosome and basal body associated protein (ALMS1), mRNA'
        assert results['NM_015120.4:c.34C>T']['gene_symbol'] == 'ALMS1'
        assert results['NM_015120.4:c.34C>T']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_055935.4(LRG_741p1):p.(Leu12=)'
        assert results['NM_015120.4:c.34C>T']['hgvs_predicted_protein_consequence']['slr'] == 'NP_055935.4:p.(L12=)'
        assert results['NM_015120.4:c.34C>T']['genome_context_intronic_sequence'] == ''
        assert results['NM_015120.4:c.34C>T']['hgvs_lrg_variant'] == 'LRG_741:g.5145C>T'
        assert results['NM_015120.4:c.34C>T']['hgvs_transcript_variant'] == 'NM_015120.4:c.34C>T'
        assert results['NM_015120.4:c.34C>T']['hgvs_refseqgene_variant'] == 'NG_011690.1:g.5145C>T'
        assert results['NM_015120.4:c.34C>T']['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000002.11:g.73613030C>T', 'vcf': {'chr': 'chr2', 'ref': u'C', 'pos': '73613030', 'alt': u'T'}}
        assert results['NM_015120.4:c.34C>T']['primary_assembly_loci']['hg38'] == {'hgvs_genomic_description': 'NC_000002.12:g.73385902C>T', 'vcf': {'chr': 'chr2', 'ref': u'C', 'pos': '73385902', 'alt': u'T'}}
        assert results['NM_015120.4:c.34C>T']['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000002.11:g.73613030C>T', 'vcf': {'chr': '2', 'ref': u'C', 'pos': '73613030', 'alt': u'T'}}
        assert results['NM_015120.4:c.34C>T']['primary_assembly_loci']['grch38'] == {'hgvs_genomic_description': 'NC_000002.12:g.73385902C>T', 'vcf': {'chr': '2', 'ref': u'C', 'pos': '73385902', 'alt': u'T'}}
        assert results['flag'] == 'gene_variant'

    def test_variant5(self):
        variant = "NC_000023.10:g.33229673A>T"
        results = vv.validator(variant, 'GRCh37', 'all')
        print results

        assert 'NM_000109.3:c.7+127703T>A' in results.keys()
        assert results['NM_000109.3:c.7+127703T>A']['hgvs_lrg_transcript_variant'] == ''
        assert results['NM_000109.3:c.7+127703T>A']['refseqgene_context_intronic_sequence'] == ''
        assert results['NM_000109.3:c.7+127703T>A']['alt_genomic_loci'] == []
        assert results['NM_000109.3:c.7+127703T>A']['transcript_description'] == 'Homo sapiens dystrophin (DMD), transcript variant Dp427c, mRNA'
        assert results['NM_000109.3:c.7+127703T>A']['gene_symbol'] == 'DMD'
        assert results['NM_000109.3:c.7+127703T>A']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_000100.2:p.?'
        assert results['NM_000109.3:c.7+127703T>A']['hgvs_predicted_protein_consequence']['slr'] == 'NP_000100.2:p.?'
        assert results['NM_000109.3:c.7+127703T>A']['genome_context_intronic_sequence'] == 'NC_000023.10(NM_000109.3):c.7+127703T>A'
        assert results['NM_000109.3:c.7+127703T>A']['hgvs_lrg_variant'] == ''
        assert results['NM_000109.3:c.7+127703T>A']['hgvs_transcript_variant'] == 'NM_000109.3:c.7+127703T>A'
        assert results['NM_000109.3:c.7+127703T>A']['hgvs_refseqgene_variant'] == ''
        assert results['NM_000109.3:c.7+127703T>A']['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000023.10:g.33229673A>T', 'vcf': {'chr': 'chrX', 'ref': u'A', 'pos': '33229673', 'alt': u'T'}}
        assert results['NM_000109.3:c.7+127703T>A']['primary_assembly_loci']['hg38'] == {'hgvs_genomic_description': 'NC_000023.11:g.33211556A>T', 'vcf': {'chr': 'chrX', 'ref': u'A', 'pos': '33211556', 'alt': u'T'}}
        assert results['NM_000109.3:c.7+127703T>A']['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000023.10:g.33229673A>T', 'vcf': {'chr': 'X', 'ref': u'A', 'pos': '33229673', 'alt': u'T'}}
        assert results['NM_000109.3:c.7+127703T>A']['primary_assembly_loci']['grch38'] == {'hgvs_genomic_description': 'NC_000023.11:g.33211556A>T', 'vcf': {'chr': 'X', 'ref': u'A', 'pos': '33211556', 'alt': u'T'}}

        assert 'NM_004006.2:c.-244T>A' in results.keys()
        assert results['NM_004006.2:c.-244T>A']['hgvs_lrg_transcript_variant'] == 'LRG_199t1:c.-244T>A'
        assert results['NM_004006.2:c.-244T>A']['refseqgene_context_intronic_sequence'] == ''
        assert results['NM_004006.2:c.-244T>A']['alt_genomic_loci'] == []
        assert results['NM_004006.2:c.-244T>A'][
                   'transcript_description'] == 'Homo sapiens dystrophin (DMD), transcript variant Dp427m, mRNA'
        assert results['NM_004006.2:c.-244T>A']['gene_symbol'] == 'DMD'
        assert results['NM_004006.2:c.-244T>A']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_003997.1(LRG_199p1):p.?'
        assert results['NM_004006.2:c.-244T>A']['hgvs_predicted_protein_consequence']['slr'] == 'NP_003997.1:p.?'
        assert results['NM_004006.2:c.-244T>A']['genome_context_intronic_sequence'] == ''
        assert results['NM_004006.2:c.-244T>A']['hgvs_lrg_variant'] == 'LRG_199:g.133054T>A'
        assert results['NM_004006.2:c.-244T>A']['hgvs_transcript_variant'] == 'NM_004006.2:c.-244T>A'
        assert results['NM_004006.2:c.-244T>A']['hgvs_refseqgene_variant'] == 'NG_012232.1:g.133054T>A'
        assert results['NM_004006.2:c.-244T>A']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000023.10:g.33229673A>T', 'vcf': {'chr': 'chrX', 'ref': u'A', 'pos': '33229673', 'alt': u'T'}}
        assert results['NM_004006.2:c.-244T>A']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000023.11:g.33211556A>T', 'vcf': {'chr': 'chrX', 'ref': u'A', 'pos': '33211556', 'alt': u'T'}}
        assert results['NM_004006.2:c.-244T>A']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000023.10:g.33229673A>T', 'vcf': {'chr': 'X', 'ref': u'A', 'pos': '33229673', 'alt': u'T'}}
        assert results['NM_004006.2:c.-244T>A']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000023.11:g.33211556A>T', 'vcf': {'chr': 'X', 'ref': u'A', 'pos': '33211556', 'alt': u'T'}}
        assert results['flag'] == 'gene_variant'

    def test_variant6(self):
        variant = "NM_001145026.1:c.715A>G"
        results = vv.validator(variant, 'GRCh37', 'all')
        print results

        assert variant in results.keys()
        assert results[variant]['hgvs_lrg_transcript_variant'] == ''
        assert results[variant]['refseqgene_context_intronic_sequence'] == ''
        assert results[variant]['alt_genomic_loci'] == []
        assert results[variant]['transcript_description'] == 'Homo sapiens protein tyrosine phosphatase, receptor type Q (PTPRQ), mRNA'
        assert results[variant]['gene_symbol'] == 'PTPRQ'
        assert results[variant]['hgvs_predicted_protein_consequence']['tlr'] == 'NP_001138498.1:p.(Arg239Gly)'
        assert results[variant]['hgvs_predicted_protein_consequence']['slr'] == 'NP_001138498.1:p.(R239G)'
        assert results[variant]['genome_context_intronic_sequence'] == ''
        assert results[variant]['hgvs_lrg_variant'] == ''
        assert results[variant]['hgvs_transcript_variant'] == 'NM_001145026.1:c.715A>G'
        assert results[variant]['hgvs_refseqgene_variant'] == ''
        assert results[variant]['primary_assembly_loci']['hg38'] == {'hgvs_genomic_description': 'NC_000012.12:g.80460707A>G', 'vcf': {'chr': 'chr12', 'ref': u'A', 'pos': '80460707', 'alt': u'G'}}
        assert results[variant]['primary_assembly_loci']['grch38'] == {'hgvs_genomic_description': 'NC_000012.12:g.80460707A>G', 'vcf': {'chr': '12', 'ref': u'A', 'pos': '80460707', 'alt': u'G'}}
        assert results['flag'] == 'gene_variant'

    def test_variant7(self):
        variant = "NC_000016.9:g.2099572TC>T"
        results = vv.validator(variant, 'GRCh37', 'all')
        print results
        assert results['flag'] == 'gene_variant'

        assert 'NM_001077183.2:c.138+821del' in results.keys()
        assert results['NM_001077183.2:c.138+821del']['hgvs_lrg_transcript_variant'] == ''
        assert results['NM_001077183.2:c.138+821del']['refseqgene_context_intronic_sequence'] == ''
        assert results['NM_001077183.2:c.138+821del']['alt_genomic_loci'] == []
        assert results['NM_001077183.2:c.138+821del']['transcript_description'] == 'Homo sapiens TSC complex subunit 2 (TSC2), transcript variant 4, mRNA'
        assert results['NM_001077183.2:c.138+821del']['gene_symbol'] == 'TSC2'
        assert results['NM_001077183.2:c.138+821del']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_001070651.1:p.?'
        assert results['NM_001077183.2:c.138+821del']['hgvs_predicted_protein_consequence']['slr'] == 'NP_001070651.1:p.?'
        assert results['NM_001077183.2:c.138+821del']['genome_context_intronic_sequence'] == 'NC_000016.9(NM_001077183.2):c.138+821del'
        assert results['NM_001077183.2:c.138+821del']['hgvs_lrg_variant'] == ''
        assert results['NM_001077183.2:c.138+821del']['hgvs_transcript_variant'] == 'NM_001077183.2:c.138+821del'
        assert results['NM_001077183.2:c.138+821del']['hgvs_refseqgene_variant'] == ''
        assert results['NM_001077183.2:c.138+821del']['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001077183.2:c.138+821del']['primary_assembly_loci']['hg38'] == {'hgvs_genomic_description': 'NC_000016.10:g.2049574del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2049571', 'alt': 'T'}}
        assert results['NM_001077183.2:c.138+821del']['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001077183.2:c.138+821del']['primary_assembly_loci']['grch38'] == {'hgvs_genomic_description': 'NC_000016.10:g.2049574del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2049571', 'alt': 'T'}}

        assert 'NM_001318831.1:c.-89+821del' in results.keys()
        assert results['NM_001318831.1:c.-89+821del']['hgvs_lrg_transcript_variant'] == ''
        assert results['NM_001318831.1:c.-89+821del']['refseqgene_context_intronic_sequence'] == ''
        assert results['NM_001318831.1:c.-89+821del']['alt_genomic_loci'] == []
        assert results['NM_001318831.1:c.-89+821del']['transcript_description'] == 'Homo sapiens TSC complex subunit 2 (TSC2), transcript variant 8, mRNA'
        assert results['NM_001318831.1:c.-89+821del']['gene_symbol'] == 'TSC2'
        assert results['NM_001318831.1:c.-89+821del']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_001305760.1:p.?'
        assert results['NM_001318831.1:c.-89+821del']['hgvs_predicted_protein_consequence']['slr'] == 'NP_001305760.1:p.?'
        assert results['NM_001318831.1:c.-89+821del']['genome_context_intronic_sequence'] == 'NC_000016.9(NM_001318831.1):c.-89+821del'
        assert results['NM_001318831.1:c.-89+821del']['hgvs_lrg_variant'] == ''
        assert results['NM_001318831.1:c.-89+821del']['hgvs_transcript_variant'] == 'NM_001318831.1:c.-89+821del'
        assert results['NM_001318831.1:c.-89+821del']['hgvs_refseqgene_variant'] == ''
        assert results['NM_001318831.1:c.-89+821del']['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001318831.1:c.-89+821del']['primary_assembly_loci']['hg38'] == {'hgvs_genomic_description': 'NC_000016.10:g.2049574del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2049571', 'alt': 'T'}}
        assert results['NM_001318831.1:c.-89+821del']['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001318831.1:c.-89+821del']['primary_assembly_loci']['grch38'] == {'hgvs_genomic_description': 'NC_000016.10:g.2049574del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2049571', 'alt': 'T'}}

        assert 'NM_001318832.1:c.171+821del' in results.keys()
        assert results['NM_001318832.1:c.171+821del']['hgvs_lrg_transcript_variant'] == ''
        assert results['NM_001318832.1:c.171+821del']['refseqgene_context_intronic_sequence'] == ''
        assert results['NM_001318832.1:c.171+821del']['alt_genomic_loci'] == []
        assert results['NM_001318832.1:c.171+821del']['transcript_description'] == 'Homo sapiens TSC complex subunit 2 (TSC2), transcript variant 9, mRNA'
        assert results['NM_001318832.1:c.171+821del']['gene_symbol'] == 'TSC2'
        assert results['NM_001318832.1:c.171+821del']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_001305761.1:p.?'
        assert results['NM_001318832.1:c.171+821del']['hgvs_predicted_protein_consequence']['slr'] == 'NP_001305761.1:p.?'
        assert results['NM_001318832.1:c.171+821del']['genome_context_intronic_sequence'] == 'NC_000016.9(NM_001318832.1):c.171+821del'
        assert results['NM_001318832.1:c.171+821del']['hgvs_lrg_variant'] == ''
        assert results['NM_001318832.1:c.171+821del']['hgvs_transcript_variant'] == 'NM_001318832.1:c.171+821del'
        assert results['NM_001318832.1:c.171+821del']['hgvs_refseqgene_variant'] == ''
        assert results['NM_001318832.1:c.171+821del']['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001318832.1:c.171+821del']['primary_assembly_loci']['hg38'] == {'hgvs_genomic_description': 'NC_000016.10:g.2049574del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2049571', 'alt': 'T'}}
        assert results['NM_001318832.1:c.171+821del']['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001318832.1:c.171+821del']['primary_assembly_loci']['grch38'] == {'hgvs_genomic_description': 'NC_000016.10:g.2049574del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2049571', 'alt': 'T'}}

        assert 'NM_001114382.1:c.138+821del' in results.keys()
        assert results['NM_001114382.1:c.138+821del']['hgvs_lrg_transcript_variant'] == ''
        assert results['NM_001114382.1:c.138+821del']['refseqgene_context_intronic_sequence'] == ''
        assert results['NM_001114382.1:c.138+821del']['alt_genomic_loci'] == []
        assert results['NM_001114382.1:c.138+821del']['transcript_description'] == 'Homo sapiens tuberous sclerosis 2 (TSC2), transcript variant 5, mRNA'
        assert results['NM_001114382.1:c.138+821del']['gene_symbol'] == 'TSC2'
        assert results['NM_001114382.1:c.138+821del']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_001107854.1:p.?'
        assert results['NM_001114382.1:c.138+821del']['hgvs_predicted_protein_consequence']['slr'] == 'NP_001107854.1:p.?'
        assert results['NM_001114382.1:c.138+821del']['genome_context_intronic_sequence'] == 'NC_000016.9(NM_001114382.1):c.138+821del'
        assert results['NM_001114382.1:c.138+821del']['hgvs_lrg_variant'] == ''
        assert results['NM_001114382.1:c.138+821del']['hgvs_transcript_variant'] == 'NM_001114382.1:c.138+821del'
        assert results['NM_001114382.1:c.138+821del']['hgvs_refseqgene_variant'] == ''
        assert results['NM_001114382.1:c.138+821del']['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001114382.1:c.138+821del']['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}

        assert 'NM_000548.4:c.138+821del' in results.keys()
        assert results['NM_000548.4:c.138+821del']['hgvs_lrg_transcript_variant'] == ''
        assert results['NM_000548.4:c.138+821del']['refseqgene_context_intronic_sequence'] == ''
        assert results['NM_000548.4:c.138+821del']['alt_genomic_loci'] == []
        assert results['NM_000548.4:c.138+821del']['transcript_description'] == 'Homo sapiens TSC complex subunit 2 (TSC2), transcript variant 1, mRNA'
        assert results['NM_000548.4:c.138+821del']['gene_symbol'] == 'TSC2'
        assert results['NM_000548.4:c.138+821del']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_000539.2(LRG_487p1):p.?'
        assert results['NM_000548.4:c.138+821del']['hgvs_predicted_protein_consequence']['slr'] == 'NP_000539.2:p.?'
        assert results['NM_000548.4:c.138+821del']['genome_context_intronic_sequence'] == 'NC_000016.9(NM_000548.4):c.138+821del'
        assert results['NM_000548.4:c.138+821del']['hgvs_lrg_variant'] == ''
        assert results['NM_000548.4:c.138+821del']['hgvs_transcript_variant'] == 'NM_000548.4:c.138+821del'
        assert results['NM_000548.4:c.138+821del']['hgvs_refseqgene_variant'] == ''
        assert results['NM_000548.4:c.138+821del']['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_000548.4:c.138+821del']['primary_assembly_loci']['hg38'] == {'hgvs_genomic_description': 'NC_000016.10:g.2049574del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2049571', 'alt': 'T'}}
        assert results['NM_000548.4:c.138+821del']['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_000548.4:c.138+821del']['primary_assembly_loci']['grch38'] == {'hgvs_genomic_description': 'NC_000016.10:g.2049574del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2049571', 'alt': 'T'}}

        assert 'NM_001363528.1:c.138+821del' in results.keys()
        assert results['NM_001363528.1:c.138+821del']['hgvs_lrg_transcript_variant'] == ''
        assert results['NM_001363528.1:c.138+821del']['refseqgene_context_intronic_sequence'] == ''
        assert results['NM_001363528.1:c.138+821del']['alt_genomic_loci'] == []
        assert results['NM_001363528.1:c.138+821del']['transcript_description'] == 'Homo sapiens TSC complex subunit 2 (TSC2), transcript variant 10, mRNA'
        assert results['NM_001363528.1:c.138+821del']['gene_symbol'] == 'TSC2'
        assert results['NM_001363528.1:c.138+821del']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_001350457.1:p.?'
        assert results['NM_001363528.1:c.138+821del']['hgvs_predicted_protein_consequence']['slr'] == 'NP_001350457.1:p.?'
        assert results['NM_001363528.1:c.138+821del']['genome_context_intronic_sequence'] == 'NC_000016.9(NM_001363528.1):c.138+821del'
        assert results['NM_001363528.1:c.138+821del']['hgvs_lrg_variant'] == ''
        assert results['NM_001363528.1:c.138+821del']['hgvs_transcript_variant'] == 'NM_001363528.1:c.138+821del'
        assert results['NM_001363528.1:c.138+821del']['hgvs_refseqgene_variant'] == ''
        assert results['NM_001363528.1:c.138+821del']['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001363528.1:c.138+821del']['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}

        assert 'NM_000548.3:c.138+821del' in results.keys()
        assert results['NM_000548.3:c.138+821del']['hgvs_lrg_transcript_variant'] == 'LRG_487t1:c.138+821del'
        assert results['NM_000548.3:c.138+821del']['refseqgene_context_intronic_sequence'] == 'NG_005895.1(NM_000548.3):c.138+821del'
        assert results['NM_000548.3:c.138+821del']['alt_genomic_loci'] == []
        assert results['NM_000548.3:c.138+821del']['transcript_description'] == 'Homo sapiens tuberous sclerosis 2 (TSC2), transcript variant 1, mRNA'
        assert results['NM_000548.3:c.138+821del']['gene_symbol'] == 'TSC2'
        assert results['NM_000548.3:c.138+821del']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_000539.2(LRG_487p1):p.?'
        assert results['NM_000548.3:c.138+821del']['hgvs_predicted_protein_consequence']['slr'] == 'NP_000539.2:p.?'
        assert results['NM_000548.3:c.138+821del']['genome_context_intronic_sequence'] == 'NC_000016.9(NM_000548.3):c.138+821del'
        assert results['NM_000548.3:c.138+821del']['hgvs_lrg_variant'] == 'LRG_487:g.5269del'
        assert results['NM_000548.3:c.138+821del']['hgvs_transcript_variant'] == 'NM_000548.3:c.138+821del'
        assert results['NM_000548.3:c.138+821del']['hgvs_refseqgene_variant'] == 'NG_005895.1:g.5269del'
        assert results['NM_000548.3:c.138+821del']['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_000548.3:c.138+821del']['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}

        assert 'NM_021055.2:c.138+821del' in results.keys()
        assert results['NM_021055.2:c.138+821del']['hgvs_lrg_transcript_variant'] == ''
        assert results['NM_021055.2:c.138+821del']['refseqgene_context_intronic_sequence'] == ''
        assert results['NM_021055.2:c.138+821del']['alt_genomic_loci'] == []
        assert results['NM_021055.2:c.138+821del']['transcript_description'] == 'Homo sapiens TSC complex subunit 2 (TSC2), transcript variant 2, mRNA'
        assert results['NM_021055.2:c.138+821del']['gene_symbol'] == 'TSC2'
        assert results['NM_021055.2:c.138+821del']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_066399.2:p.?'
        assert results['NM_021055.2:c.138+821del']['hgvs_predicted_protein_consequence']['slr'] == 'NP_066399.2:p.?'
        assert results['NM_021055.2:c.138+821del']['genome_context_intronic_sequence'] == 'NC_000016.9(NM_021055.2):c.138+821del'
        assert results['NM_021055.2:c.138+821del']['hgvs_lrg_variant'] == ''
        assert results['NM_021055.2:c.138+821del']['hgvs_transcript_variant'] == 'NM_021055.2:c.138+821del'
        assert results['NM_021055.2:c.138+821del']['hgvs_refseqgene_variant'] == ''
        assert results['NM_021055.2:c.138+821del']['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_021055.2:c.138+821del']['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}

        assert 'NM_001077183.1:c.138+821del' in results.keys()
        assert results['NM_001077183.1:c.138+821del']['hgvs_lrg_transcript_variant'] == ''
        assert results['NM_001077183.1:c.138+821del']['refseqgene_context_intronic_sequence'] == ''
        assert results['NM_001077183.1:c.138+821del']['alt_genomic_loci'] == []
        assert results['NM_001077183.1:c.138+821del']['transcript_description'] == 'Homo sapiens tuberous sclerosis 2 (TSC2), transcript variant 4, mRNA'
        assert results['NM_001077183.1:c.138+821del']['gene_symbol'] == 'TSC2'
        assert results['NM_001077183.1:c.138+821del']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_001070651.1:p.?'
        assert results['NM_001077183.1:c.138+821del']['hgvs_predicted_protein_consequence']['slr'] == 'NP_001070651.1:p.?'
        assert results['NM_001077183.1:c.138+821del']['genome_context_intronic_sequence'] == 'NC_000016.9(NM_001077183.1):c.138+821del'
        assert results['NM_001077183.1:c.138+821del']['hgvs_lrg_variant'] == ''
        assert results['NM_001077183.1:c.138+821del']['hgvs_transcript_variant'] == 'NM_001077183.1:c.138+821del'
        assert results['NM_001077183.1:c.138+821del']['hgvs_refseqgene_variant'] == ''
        assert results['NM_001077183.1:c.138+821del']['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001077183.1:c.138+821del']['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}

        assert 'NM_001318827.1:c.138+821del' in results.keys()
        assert results['NM_001318827.1:c.138+821del']['hgvs_lrg_transcript_variant'] == ''
        assert results['NM_001318827.1:c.138+821del']['refseqgene_context_intronic_sequence'] == ''
        assert results['NM_001318827.1:c.138+821del']['alt_genomic_loci'] == []
        assert results['NM_001318827.1:c.138+821del']['transcript_description'] == 'Homo sapiens TSC complex subunit 2 (TSC2), transcript variant 6, mRNA'
        assert results['NM_001318827.1:c.138+821del']['gene_symbol'] == 'TSC2'
        assert results['NM_001318827.1:c.138+821del']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_001305756.1:p.?'
        assert results['NM_001318827.1:c.138+821del']['hgvs_predicted_protein_consequence']['slr'] == 'NP_001305756.1:p.?'
        assert results['NM_001318827.1:c.138+821del']['genome_context_intronic_sequence'] == 'NC_000016.9(NM_001318827.1):c.138+821del'
        assert results['NM_001318827.1:c.138+821del']['hgvs_lrg_variant'] == ''
        assert results['NM_001318827.1:c.138+821del']['hgvs_transcript_variant'] == 'NM_001318827.1:c.138+821del'
        assert results['NM_001318827.1:c.138+821del']['hgvs_refseqgene_variant'] == ''
        assert results['NM_001318827.1:c.138+821del']['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001318827.1:c.138+821del']['primary_assembly_loci']['hg38'] == {'hgvs_genomic_description': 'NC_000016.10:g.2049574del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2049571', 'alt': 'T'}}
        assert results['NM_001318827.1:c.138+821del']['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001318827.1:c.138+821del']['primary_assembly_loci']['grch38'] == {'hgvs_genomic_description': 'NC_000016.10:g.2049574del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2049571', 'alt': 'T'}}

        assert 'NM_001114382.2:c.138+821del' in results.keys()
        assert results['NM_001114382.2:c.138+821del']['hgvs_lrg_transcript_variant'] == ''
        assert results['NM_001114382.2:c.138+821del']['refseqgene_context_intronic_sequence'] == ''
        assert results['NM_001114382.2:c.138+821del']['alt_genomic_loci'] == []
        assert results['NM_001114382.2:c.138+821del']['transcript_description'] == 'Homo sapiens TSC complex subunit 2 (TSC2), transcript variant 5, mRNA'
        assert results['NM_001114382.2:c.138+821del']['gene_symbol'] == 'TSC2'
        assert results['NM_001114382.2:c.138+821del']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_001107854.1:p.?'
        assert results['NM_001114382.2:c.138+821del']['hgvs_predicted_protein_consequence']['slr'] == 'NP_001107854.1:p.?'
        assert results['NM_001114382.2:c.138+821del']['genome_context_intronic_sequence'] == 'NC_000016.9(NM_001114382.2):c.138+821del'
        assert results['NM_001114382.2:c.138+821del']['hgvs_lrg_variant'] == ''
        assert results['NM_001114382.2:c.138+821del']['hgvs_transcript_variant'] == 'NM_001114382.2:c.138+821del'
        assert results['NM_001114382.2:c.138+821del']['hgvs_refseqgene_variant'] == ''
        assert results['NM_001114382.2:c.138+821del']['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001114382.2:c.138+821del']['primary_assembly_loci']['hg38'] == {'hgvs_genomic_description': 'NC_000016.10:g.2049574del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2049571', 'alt': 'T'}}
        assert results['NM_001114382.2:c.138+821del']['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001114382.2:c.138+821del']['primary_assembly_loci']['grch38'] == {'hgvs_genomic_description': 'NC_000016.10:g.2049574del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2049571', 'alt': 'T'}}

        assert 'NM_001318829.1:c.-9-826del' in results.keys()
        assert results['NM_001318829.1:c.-9-826del']['hgvs_lrg_transcript_variant'] == ''
        assert results['NM_001318829.1:c.-9-826del']['refseqgene_context_intronic_sequence'] == ''
        assert results['NM_001318829.1:c.-9-826del']['alt_genomic_loci'] == []
        assert results['NM_001318829.1:c.-9-826del']['transcript_description'] == 'Homo sapiens TSC complex subunit 2 (TSC2), transcript variant 7, mRNA'
        assert results['NM_001318829.1:c.-9-826del']['gene_symbol'] == 'TSC2'
        assert results['NM_001318829.1:c.-9-826del']['hgvs_predicted_protein_consequence']['tlr'] == 'NP_001305758.1:p.?'
        assert results['NM_001318829.1:c.-9-826del']['hgvs_predicted_protein_consequence']['slr'] == 'NP_001305758.1:p.?'
        assert results['NM_001318829.1:c.-9-826del']['genome_context_intronic_sequence'] == 'NC_000016.9(NM_001318829.1):c.-9-826del'
        assert results['NM_001318829.1:c.-9-826del']['hgvs_lrg_variant'] == ''
        assert results['NM_001318829.1:c.-9-826del']['hgvs_transcript_variant'] == 'NM_001318829.1:c.-9-826del'
        assert results['NM_001318829.1:c.-9-826del']['hgvs_refseqgene_variant'] == ''
        assert results['NM_001318829.1:c.-9-826del']['primary_assembly_loci']['hg19'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001318829.1:c.-9-826del']['primary_assembly_loci']['hg38'] == {'hgvs_genomic_description': 'NC_000016.10:g.2049574del', 'vcf': {'chr': 'chr16', 'ref': 'TC', 'pos': '2049571', 'alt': 'T'}}
        assert results['NM_001318829.1:c.-9-826del']['primary_assembly_loci']['grch37'] == {'hgvs_genomic_description': 'NC_000016.9:g.2099575del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2099572', 'alt': 'T'}}
        assert results['NM_001318829.1:c.-9-826del']['primary_assembly_loci']['grch38'] == {'hgvs_genomic_description': 'NC_000016.10:g.2049574del', 'vcf': {'chr': '16', 'ref': 'TC', 'pos': '2049571', 'alt': 'T'}}

    @pytest.mark.skip(reason="template for copy/paste")
    def test_variant_blank(self):
        variant = ""
        results = vv.validator(variant, 'GRCh37', 'all')
        print results

        assert variant in results.keys()
        assert results[variant]['hgvs_lrg_transcript_variant'] == ''
        assert results[variant]['refseqgene_context_intronic_sequence'] == ''
        assert results[variant]['alt_genomic_loci'] == []
        assert results[variant]['transcript_description'] == ''
        assert results[variant]['gene_symbol'] == ''
        assert results[variant]['hgvs_predicted_protein_consequence']['tlr'] == ''
        assert results[variant]['hgvs_predicted_protein_consequence']['slr'] == ''
        assert results[variant]['genome_context_intronic_sequence'] == ''
        assert results[variant]['hgvs_lrg_variant'] == ''
        assert results[variant]['hgvs_transcript_variant'] == ''
        assert results[variant]['hgvs_refseqgene_variant'] == ''
        assert results[variant]['primary_assembly_loci']['hg19'] == {}
        assert results[variant]['primary_assembly_loci']['hg38'] == {}
        assert results[variant]['primary_assembly_loci']['grch37'] == {}
        assert results[variant]['primary_assembly_loci']['grch38'] == {}
        assert results['flag'] == ''