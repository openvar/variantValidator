from VariantValidator import Validator
from unittest import TestCase


class TestVariantsEnsembl(TestCase):

    @classmethod
    def setup_class(cls):
        cls.vv = Validator()

    # COL1A1
    def test_variant1(self):
        variant = 'ENST00000225964.10:c.589-1GG>G'
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000225964.10:c.590del' in list(results.keys())
        assert results['ENST00000225964.10:c.590del']['submitted_variant'] == 'ENST00000225964.10:c.589-1GG>G'
        assert results['ENST00000225964.10:c.590del']['gene_symbol'] == 'COL1A1'
        assert results['ENST00000225964.10:c.590del']['gene_ids'] == {'hgnc_id': 'HGNC:2197', 'entrez_gene_id': '1277',
                                                               'ucsc_id': 'uc002iqm.4', 'omim_id': ['120150']}
        assert results['ENST00000225964.10:c.590del']['hgvs_transcript_variant'] == 'ENST00000225964.10:c.590del'
        assert results['ENST00000225964.10:c.590del']['genome_context_intronic_sequence'] == ''
        # assert results['ENST00000225964.10:c.590del']['refseqgene_context_intronic_sequence'] == ''
        # assert results['ENST00000225964.10:c.590del']['hgvs_refseqgene_variant'] == 'NG_007400.1:g.8639del'
        assert results['ENST00000225964.10:c.590del']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000225964.6:p.(Gly197ValfsTer68)', 'slr': 'ENSP00000225964.6:p.(G197Vfs*68)'}
        # assert results['ENST00000225964.10:c.590del']['hgvs_lrg_transcript_variant'] == 'LRG_1t1:c.590del'
        # assert results['ENST00000225964.10:c.590del']['hgvs_lrg_variant'] == 'LRG_1:g.8639del'
        self.assertCountEqual(results['ENST00000225964.10:c.590del']['alt_genomic_loci'], [])
        assert results['ENST00000225964.10:c.590del']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000017.10:g.48275364del',
            'vcf': {'chr': 'chr17', 'pos': '48275361', 'ref': 'AC', 'alt': 'A'}}
        assert results['ENST00000225964.10:c.590del']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000017.11:g.50198003del',
            'vcf': {'chr': 'chr17', 'pos': '50198000', 'ref': 'AC', 'alt': 'A'}}
        assert results['ENST00000225964.10:c.590del']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000017.10:g.48275364del',
            'vcf': {'chr': '17', 'pos': '48275361', 'ref': 'AC', 'alt': 'A'}}
        assert results['ENST00000225964.10:c.590del']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000017.11:g.50198003del',
            'vcf': {'chr': '17', 'pos': '50198000', 'ref': 'AC', 'alt': 'A'}}
        assert results['ENST00000225964.10:c.590del']['reference_sequence_records'] == {
             'transcript': 'https://www.ncbi.nlm.nih.gov/nuccore/NM_000088.4',
             'protein': 'https://www.ncbi.nlm.nih.gov/nuccore/NP_000079.2',
            # 'refseqgene': 'https://www.ncbi.nlm.nih.gov/nuccore/NG_007400.1',
            # 'lrg': 'http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_1.xml'
            }

    # COL5A1
    def test_variant2(self):
        variant = 'ENST00000371817.8:c.5071A>T'
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000371817.8:c.5071A>T' in list(results.keys())
        assert results['ENST00000371817.8:c.5071A>T']['submitted_variant'] == 'ENST00000371817.8:c.5071A>T'
        assert results['ENST00000371817.8:c.5071A>T']['gene_symbol'] == 'COL5A1'
        assert results['ENST00000371817.8:c.5071A>T']['gene_ids'] == {'hgnc_id': 'HGNC:2209', 'entrez_gene_id': '1289',
                                                               'ucsc_id': 'uc004cfe.5', 'omim_id': ['120215']}
        assert results['ENST00000371817.8:c.5071A>T']['hgvs_transcript_variant'] == 'ENST00000371817.8:c.5071A>T'
        assert results['ENST00000371817.8:c.5071A>T']['genome_context_intronic_sequence'] == ''
        # assert results['ENST00000371817.8:c.5071A>T']['refseqgene_context_intronic_sequence'] == ''
        # assert results['ENST00000371817.8:c.5071A>T']['hgvs_refseqgene_variant'] == ''
        assert results['ENST00000371817.8:c.5071A>T']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000360882.3:p.(Arg1691Ter)', 'slr': 'ENSP00000360882.3p.(R1691*)'}
        # assert results['ENST00000371817.8:c.5071A>T']['hgvs_lrg_transcript_variant'] == ''
        # assert results['ENST00000371817.8:c.5071A>T']['hgvs_lrg_variant'] == ''
        self.assertCountEqual(results['ENST00000371817.8:c.5071A>T']['alt_genomic_loci'], [])
        assert results['ENST00000371817.8:c.5071A>T']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000009.11:g.137721825A>T', 
            'vcf': {'chr': 'chr9', 'pos': '137721825', 'ref': 'A', 'alt': 'T'}}
        assert results['ENST00000371817.8:c.5071A>T']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000009.12:g.134829979A>T', 
            'vcf': {'chr': 'chr9', 'pos': '134829979', 'ref': 'A', 'alt': 'T'}}
        assert results['ENST00000371817.8:c.5071A>T']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000009.11:g.137721825A>T', 
            'vcf': {'chr': '9', 'pos': '137721825', 'ref': 'A', 'alt': 'T'}}
        assert results['ENST00000371817.8:c.5071A>T']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000009.12:g.134829979A>T', 
            'vcf': {'chr': '9', 'pos': '134829979', 'ref': 'A', 'alt': 'T'}}
        assert results['ENST00000371817.8:c.5071A>T']['reference_sequence_records'] == {
             'transcript': 'https://www.ncbi.nlm.nih.gov/nuccore/NM_000093.5',
             'protein': 'https://www.ncbi.nlm.nih.gov/nuccore/NP_000084.3'}

    # TP53
    def test_variant3(self):
        variant = 'ENST00000269305.9:c.652_654del'
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000269305.9:c.652_654del' in list(results.keys())
        assert results['ENST00000269305.9:c.652_654del']['submitted_variant'] == 'ENST00000269305.9:c.652_654del'
        assert results['ENST00000269305.9:c.652_654del']['gene_symbol'] == 'TP53'
        assert results['ENST00000269305.9:c.652_654del']['gene_ids'] == {'hgnc_id': 'HGNC:11998', 'entrez_gene_id': '7157',
                                                                   'ucsc_id': 'uc060aur.1', 'omim_id': ['191170']}
        assert results['ENST00000269305.9:c.652_654del']['hgvs_transcript_variant'] == 'ENST00000269305.9:c.652_654del'
        assert results['ENST00000269305.9:c.652_654del']['genome_context_intronic_sequence'] == ''
        # assert results['ENST00000269305.9:c.652_654del']['refseqgene_context_intronic_sequence'] == ''
        # assert results['ENST00000269305.9:c.652_654del']['hgvs_refseqgene_variant'] == 'NG_017013.2:g.17672_17674del'
        assert results['ENST00000269305.9:c.652_654del']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000269305.4:p.(Val218del)', 'slr': 'ENSP00000269305.4:p.(V218del)'}
        # assert results['ENST00000269305.9:c.652_654del']['hgvs_lrg_transcript_variant'] == 'LRG_321t1:c.652_654del'
        # assert results['ENST00000269305.9:c.652_654del']['hgvs_lrg_variant'] == 'LRG_321:g.17672_17674del'
        self.assertCountEqual(results['ENST00000269305.9:c.652_654del']['alt_genomic_loci'], [])
        assert results['ENST00000269305.9:c.652_654del']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000017.10:g.7578201_7578203del',
            'vcf': {'chr': 'chr17', 'pos': '7578194', 'ref': 'GCAC', 'alt': 'G'}}
        assert results['ENST00000269305.9:c.652_654del']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000017.11:g.7674883_7674885del',
            'vcf': {'chr': 'chr17', 'pos': '7674876', 'ref': 'GCAC', 'alt': 'G'}}
        assert results['ENST00000269305.9:c.652_654del']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000017.10:g.7578201_7578203del',
            'vcf': {'chr': '17', 'pos': '7578194', 'ref': 'GCAC', 'alt': 'G'}}
        assert results['ENST00000269305.9:c.652_654del']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000017.11:g.7674883_7674885del',
            'vcf': {'chr': '17', 'pos': '7674876', 'ref': 'GCAC', 'alt': 'G'}}
        assert results['ENST00000269305.9:c.652_654del']['reference_sequence_records'] == {
            'transcript': 'https://www.ncbi.nlm.nih.gov/nuccore/NM_000546.6',
            'protein': 'https://www.ncbi.nlm.nih.gov/nuccore/NP_000537.3'
            #'refseqgene': 'https://www.ncbi.nlm.nih.gov/nuccore/NG_017013.2',
            #'lrg': 'http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_321.xml'
            }

    # P3H1
    def test_variant4(self):
        variant = 'ENST00000296388.10:c.2055+18G>A'
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000296388.10:c.2055+18G>A' in list(results.keys())
        assert results['ENST00000296388.10:c.2055+18G>A']['submitted_variant'] == 'ENST00000296388.10:c.2055+18G>A'
        assert results['ENST00000296388.10:c.2055+18G>A']['gene_symbol'] == 'P3H1'
        assert results['ENST00000296388.10:c.2055+18G>A']['gene_ids'] == {'hgnc_id': 'HGNC:19316', 'entrez_gene_id': '64175',
                                                                   'ucsc_id': '', 'omim_id': ['610339']}
        assert results['ENST00000296388.10:c.2055+18G>A']['hgvs_transcript_variant'] == 'ENST00000296388.10:c.2055+18G>A'
        assert results['ENST00000296388.10:c.2055+18G>A'][
                   'genome_context_intronic_sequence'] == 'NC_000001.10(ENST00000296388.10):c.2055+18G>A'
        # assert results['ENST00000296388.10:c.2055+18G>A'][
        #     'refseqgene_context_intronic_sequence'] == 'NG_008123.1(ENST00000296388.10):c.2055+18G>A'
        # assert results['ENST00000296388.10:c.2055+18G>A']['hgvs_refseqgene_variant'] == 'NG_008123.1:g.24831G>A'
        assert results['ENST00000296388.10:c.2055+18G>A']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000296388.5:p.?', 'slr': 'ENSP00000296388.5:p.?'}
        # assert results['ENST00000296388.10:c.2055+18G>A']['hgvs_lrg_transcript_variant'] == 'LRG_5t1:c.2055+18G>A'
        # assert results['ENST00000296388.10:c.2055+18G>A']['hgvs_lrg_variant'] == 'LRG_5:g.24831G>A'
        self.assertCountEqual(results['ENST00000296388.10:c.2055+18G>A']['alt_genomic_loci'], [])
        assert results['ENST00000296388.10:c.2055+18G>A']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000001.10:g.43212925C>T',
            'vcf': {'chr': 'chr1', 'pos': '43212925', 'ref': 'C', 'alt': 'T'}}
        assert results['ENST00000296388.10:c.2055+18G>A']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000001.11:g.42747254C>T',
            'vcf': {'chr': 'chr1', 'pos': '42747254', 'ref': 'C', 'alt': 'T'}}
        assert results['ENST00000296388.10:c.2055+18G>A']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000001.10:g.43212925C>T',
            'vcf': {'chr': '1', 'pos': '43212925', 'ref': 'C', 'alt': 'T'}}
        assert results['ENST00000296388.10:c.2055+18G>A']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000001.11:g.42747254C>T',
            'vcf': {'chr': '1', 'pos': '42747254', 'ref': 'C', 'alt': 'T'}}
        assert results['ENST00000296388.10:c.2055+18G>A']['reference_sequence_records'] == {
            'transcript': 'https://www.ncbi.nlm.nih.gov/nuccore/NM_022356.4',
            'protein': 'https://www.ncbi.nlm.nih.gov/nuccore/NP_071751.3',
            #'refseqgene': 'https://www.ncbi.nlm.nih.gov/nuccore/NG_008123.1',
            #'lrg': 'http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_5.xml'
            }

    # BRCA1
    def test_variant5(self):
        variant = 'ENST00000357654.9:c.301+1G>C'
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        
        assert 'ENST00000357654.9:c.301+1G>C' in list(results.keys())
        assert results['ENST00000357654.9:c.301+1G>C']['submitted_variant'] == 'ENST00000357654.9:c.301+1G>C'
        assert results['ENST00000357654.9:c.301+1G>C']['gene_symbol'] == 'BRCA1'
        assert results['ENST00000357654.9:c.301+1G>C']['gene_ids'] == {'hgnc_id': 'HGNC:1100', 'entrez_gene_id': '672',
                                                                 'ucsc_id': 'uc002ict.4', 'omim_id': ['113705']}
        assert results['ENST00000357654.9:c.301+1G>C']['hgvs_transcript_variant'] == 'ENST00000357654.9:c.301+1G>C'
        assert results['ENST00000357654.9:c.301+1G>C'][
                   'genome_context_intronic_sequence'] == 'NC_000017.10(ENST00000357654.9):c.301+1G>C'
        # assert results['ENST00000357654.9:c.301+1G>C'][
        #           'refseqgene_context_intronic_sequence'] == 'NG_005905.2(ENST00000357654.9):c.301+1G>C'
        # assert results['ENST00000357654.9:c.301+1G>C']['hgvs_refseqgene_variant'] == 'NG_005905.2:g.113117G>C'
        assert results['ENST00000357654.9:c.301+1G>C']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'NP_009225.1:p.?', 'slr': 'NP_009225.1:p.?'}
        # assert results['ENST00000357654.9:c.301+1G>C']['hgvs_lrg_transcript_variant'] == 'LRG_292t1:c.301+1G>C'
        # assert results['ENST00000357654.9:c.301+1G>C']['hgvs_lrg_variant'] == 'LRG_292:g.113117G>C'
        self.assertCountEqual(results['ENST00000357654.9:c.301+1G>C']['alt_genomic_loci'], [])
        assert results['ENST00000357654.9:c.301+1G>C']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000017.10:g.41256884C>G',
            'vcf': {'chr': 'chr17', 'pos': '41256884', 'ref': 'C', 'alt': 'G'}}
        assert results['ENST00000357654.9:c.301+1G>C']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000017.11:g.43104867C>G',
            'vcf': {'chr': 'chr17', 'pos': '43104867', 'ref': 'C', 'alt': 'G'}}
        assert results['ENST00000357654.9:c.301+1G>C']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000017.10:g.41256884C>G',
            'vcf': {'chr': '17', 'pos': '41256884', 'ref': 'C', 'alt': 'G'}}
        assert results['ENST00000357654.9:c.301+1G>C']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000017.11:g.43104867C>G',
            'vcf': {'chr': '17', 'pos': '43104867', 'ref': 'C', 'alt': 'G'}}
        assert results['ENST00000357654.9:c.301+1G>C']['reference_sequence_records'] == {
            'transcript': 'https://www.ncbi.nlm.nih.gov/nuccore/NM_007294.4',
            'protein': 'https://www.ncbi.nlm.nih.gov/nuccore/NP_009225.1',
            #'refseqgene': 'https://www.ncbi.nlm.nih.gov/nuccore/NG_005905.2',
            #'lrg': 'http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_292.xml'
            }
    
    # BRCA2
    def test_variant6(self):
        variant = 'NC_000013.10:g.32929387T>C'
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000380152.8:c.7397=' in list(results.keys())
        assert results['ENST00000380152.8:c.7397=']['submitted_variant'] == 'NC_000013.10:g.32929387T>C'
        assert results['ENST00000380152.8:c.7397=']['gene_symbol'] == 'BRCA2'
        assert results['ENST00000380152.8:c.7397=']['gene_ids'] == {'hgnc_id': 'HGNC:1101', 'entrez_gene_id': '675',
                                                               'ucsc_id': 'uc001uub.2', 'omim_id': ['600185']}
        assert results['ENST00000380152.8:c.7397=']['hgvs_transcript_variant'] == 'ENST00000380152.8:c.7397='
        assert results['ENST00000380152.8:c.7397=']['genome_context_intronic_sequence'] == ''
        # assert results['ENST00000380152.8:c.7397=']['refseqgene_context_intronic_sequence'] == ''
        # assert results['ENST00000380152.8:c.7397=']['hgvs_refseqgene_variant'] == 'NG_012772.3:g.44771='
        assert results['ENST00000380152.8:c.7397=']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000369497.3:p.(Ala2466=)', 'slr': 'ENSP00000369497.3:p.(A2466=)'}
        # assert results['ENST00000380152.8:c.7397=']['hgvs_lrg_transcript_variant'] == 'LRG_293t1:c.7397='
        # assert results['ENST00000380152.8:c.7397=']['hgvs_lrg_variant'] == 'LRG_293:g.44771='
        self.assertCountEqual(results['ENST00000380152.8:c.7397=']['alt_genomic_loci'], [])
        assert results['ENST00000380152.8:c.7397=']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000013.10:g.32929387T>C',
            'vcf': {'chr': 'chr13', 'pos': '32929387', 'ref': 'T', 'alt': 'C'}}
        assert results['ENST00000380152.8:c.7397=']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000013.11:g.32355250T>C',
            'vcf': {'chr': 'chr13', 'pos': '32355250', 'ref': 'T', 'alt': 'C'}}
        assert results['ENST00000380152.8:c.7397=']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000013.10:g.32929387T>C',
            'vcf': {'chr': '13', 'pos': '32929387', 'ref': 'T', 'alt': 'C'}}
        assert results['ENST00000380152.8:c.7397=']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000013.11:g.32355250T>C',
            'vcf': {'chr': '13', 'pos': '32355250', 'ref': 'T', 'alt': 'C'}}
        assert results['ENST00000380152.8:c.7397=']['reference_sequence_records'] == {
            'transcript': 'https://www.ncbi.nlm.nih.gov/nuccore/NM_000059.4',
            'protein': 'https://www.ncbi.nlm.nih.gov/nuccore/NP_000050.3',
            #'refseqgene': 'https://www.ncbi.nlm.nih.gov/nuccore/NG_012772.3',
            #'lrg': 'http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_293.xml'
            }
    
    # HBB
    def test_variant7(self):
        variant = '11-5248232-T-A'
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000335295.4:c.20A>T' in list(results.keys())
        assert results['ENST00000335295.4:c.20A>T']['submitted_variant'] == '11-5248232-T-A'
        assert results['ENST00000335295.4:c.20A>T']['gene_symbol'] == 'HBB'
        assert results['ENST00000335295.4:c.20A>T']['gene_ids'] == {'hgnc_id': 'HGNC:4827', 'entrez_gene_id': '3043',
                                                              'ucsc_id': 'uc001mae.2', 'omim_id': ['141900']}
        assert results['ENST00000335295.4:c.20A>T']['hgvs_transcript_variant'] == 'ENST00000335295.4:c.20A>T'
        assert results['ENST00000335295.4:c.20A>T']['genome_context_intronic_sequence'] == ''
        # assert results['ENST00000335295.4:c.20A>T']['refseqgene_context_intronic_sequence'] == ''
        # assert results['ENST00000335295.4:c.20A>T']['hgvs_refseqgene_variant'] == 'NG_059281.1:g.5070A>T'
        assert results['ENST00000335295.4:c.20A>T']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000333994.3:p.(Glu7Val)', 'slr': 'ENSP00000333994.3:p.(E7V)'}
        # assert results['ENST00000335295.4:c.20A>T']['hgvs_lrg_transcript_variant'] == 'LRG_1232t1:c.20A>T'
        # assert results['ENST00000335295.4:c.20A>T']['hgvs_lrg_variant'] == 'LRG_1232:g.5070A>T'
        self.assertCountEqual(results['ENST00000335295.4:c.20A>T']['alt_genomic_loci'], [])
        assert results['ENST00000335295.4:c.20A>T']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000011.9:g.5248232T>A',
            'vcf': {'chr': 'chr11', 'pos': '5248232', 'ref': 'T', 'alt': 'A'}}
        assert results['ENST00000335295.4:c.20A>T']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000011.10:g.5227002T>A', 
            'vcf': {'chr': 'chr11', 'pos': '5227002', 'ref': 'T', 'alt': 'A'}}
        assert results['ENST00000335295.4:c.20A>T']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000011.9:g.5248232T>A',
            'vcf': {'chr': '11', 'pos': '5248232', 'ref': 'T', 'alt': 'A'}}
        assert results['ENST00000335295.4:c.20A>T']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000011.10:g.5227002T>A', 
            'vcf': {'chr': '11', 'pos': '5227002', 'ref': 'T', 'alt': 'A'}}
        assert results['ENST00000335295.4:c.20A>T']['reference_sequence_records'] == {
            'transcript': 'https://www.ncbi.nlm.nih.gov/nuccore/NM_000518.5',
            'protein': 'https://www.ncbi.nlm.nih.gov/nuccore/NP_000509.1',
            #'refseqgene': 'https://www.ncbi.nlm.nih.gov/nuccore/NG_059281.1',
            #'lrg': 'http://ftp.ebi.ac.uk/pub/databases/lrgex/pending/LRG_1232.xml'
            }
    
    # MCL1
    def test_variant8(self):
        variant = '1:150550916G>A'
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert 'ENST00000369026.3:c.740C>T' in list(results.keys())
        assert results['ENST00000369026.3:c.740C>T']['submitted_variant'] == '1:150550916G>A'
        assert results['ENST00000369026.3:c.740C>T']['gene_symbol'] == 'MCL1'
        assert results['ENST00000369026.3:c.740C>T']['gene_ids'] == {'hgnc_id': 'HGNC:6943', 'entrez_gene_id': '4170',
                                                               'ucsc_id': 'uc001euz.4', 'omim_id': ['159552']}
        assert results['ENST00000369026.3:c.740C>T']['hgvs_transcript_variant'] == 'ENST00000369026.3:c.740C>T'
        assert results['ENST00000369026.3:c.740C>T']['genome_context_intronic_sequence'] == ''
        # assert results['ENST00000369026.3:c.740C>T']['refseqgene_context_intronic_sequence'] == ''
        # assert results['ENST00000369026.3:c.740C>T']['hgvs_refseqgene_variant'] == 'NG_029146.1:g.6299C>T'
        assert results['ENST00000369026.3:c.740C>T']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000358022.2:p.(Ser247Phe)', 'slr': 'ENSP00000358022.2:p.(S247F)'}
        # assert results['ENST00000369026.3:c.740C>T']['hgvs_lrg_transcript_variant'] == ''
        # assert results['ENST00000369026.3:c.740C>T']['hgvs_lrg_variant'] == ''
        self.assertCountEqual(results['ENST00000369026.3:c.740C>T']['alt_genomic_loci'], [])
        assert results['ENST00000369026.3:c.740C>T']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000001.10:g.150550916G>A',
            'vcf': {'chr': 'chr1', 'pos': '150550916', 'ref': 'G', 'alt': 'A'}}
        assert results['ENST00000369026.3:c.740C>T']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000001.11:g.150578440G>A',
            'vcf': {'chr': 'chr1', 'pos': '150578440', 'ref': 'G', 'alt': 'A'}}
        assert results['ENST00000369026.3:c.740C>T']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000001.10:g.150550916G>A',
            'vcf': {'chr': '1', 'pos': '150550916', 'ref': 'G', 'alt': 'A'}}
        assert results['ENST00000369026.3:c.740C>T']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000001.11:g.150578440G>A',
            'vcf': {'chr': '1', 'pos': '150578440', 'ref': 'G', 'alt': 'A'}}
        assert results['ENST00000369026.3:c.740C>T']['reference_sequence_records'] == {
            'transcript': 'https://www.ncbi.nlm.nih.gov/nuccore/NM_021960.5',
            'protein': 'https://www.ncbi.nlm.nih.gov/nuccore/NP_068779.1',
            #'refseqgene': 'https://www.ncbi.nlm.nih.gov/nuccore/NG_029146.1'
            }

    # TNXB
    def test_variant9(self):
        variant = '6-32012992-CG-C'
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert 'ENST00000644971.2:c.10717del' in list(results.keys())
        assert results['ENST00000644971.2:c.10717del']['submitted_variant'] == '6-32012992-CG-C'
        assert results['ENST00000644971.2:c.10717del']['gene_symbol'] == 'TNXB'
        assert results['ENST00000644971.2:c.10717del']['gene_ids'] == {'hgnc_id': 'HGNC:11976', 'entrez_gene_id': '7148',
                                                                    'ucsc_id': 'uc063nnw.1', 'omim_id': ['600985']}
        assert results['ENST00000644971.2:c.10717del']['hgvs_transcript_variant'] == 'ENST00000644971.2:c.10717del'
        assert results['ENST00000644971.2:c.10717del']['genome_context_intronic_sequence'] == ''
        # assert results['ENST00000644971.2:c.10717del']['refseqgene_context_intronic_sequence'] == ''
        # assert results['ENST00000644971.2:c.10717del']['hgvs_refseqgene_variant'] == ''
        assert results['ENST00000644971.2:c.10717del']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000496448.1:p.(Arg3573AlafsTer91)', 'slr': 'ENSP00000496448.1:p.(R3573Afs*91)'}
        # assert results['ENST00000644971.2:c.10717del']['hgvs_lrg_transcript_variant'] == ''
        # assert results['ENST00000644971.2:c.10717del']['hgvs_lrg_variant'] == ''
        self.assertCountEqual(results['ENST00000644971.2:c.10717del']['alt_genomic_loci'], [{'grch37': {
            'hgvs_genomic_description': 'NT_167245.1:g.3292210del',
            'vcf': {'chr': 'HSCHR6_MHC_DBB_CTG1', 'pos': '3292209', 'ref': 'CG', 'alt': 'C'}}}, {'hg19': {
            'hgvs_genomic_description': 'NT_167245.1:g.3292210del',
            'vcf': {'chr': 'chr6_dbb_hap3', 'pos': '3292209', 'ref': 'CG', 'alt': 'C'}}}, {'grch37': {
            'hgvs_genomic_description': 'NT_113891.2:g.3483644del',
            'vcf': {'chr': 'HSCHR6_MHC_COX_CTG1', 'pos': '3483643', 'ref': 'CG', 'alt': 'C'}}}, {'hg19': {
            'hgvs_genomic_description': 'NT_113891.2:g.3483644del',
            'vcf': {'chr': 'chr6_cox_hap2', 'pos': '3483643', 'ref': 'CG', 'alt': 'C'}}}, {'grch37': {
            'hgvs_genomic_description': 'NT_167247.1:g.3392834del',
            'vcf': {'chr': 'HSCHR6_MHC_MCF_CTG1', 'pos': '3392833', 'ref': 'CG', 'alt': 'C'}}}, {'hg19': {
            'hgvs_genomic_description': 'NT_167247.1:g.3392834del',
            'vcf': {'chr': 'chr6_mcf_hap5', 'pos': '3392833', 'ref': 'CG', 'alt': 'C'}}}])
        assert results['ENST00000644971.2:c.10717del']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000006.11:g.32012993del',
            'vcf': {'chr': '6', 'pos': '32012992', 'ref': 'CG', 'alt': 'C'}}
        assert results['ENST00000644971.2:c.10717del']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000006.12:g.32045216del', 
            'vcf': {'chr': 'chr6', 'pos': '32045215', 'ref': 'CG', 'alt': 'C'}}
        assert results['ENST00000644971.2:c.10717del']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000006.11:g.32012993del',
            'vcf': {'chr': '6', 'pos': '32012992', 'ref': 'CG', 'alt': 'C'}}
        assert results['ENST00000644971.2:c.10717del']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000006.12:g.32045216del', 
            'vcf': {'chr': '6', 'pos': '32045215', 'ref': 'CG', 'alt': 'C'}}
        assert results['ENST00000644971.2.1:c.10717del']['reference_sequence_records'] == {
            'transcript': 'https://www.ncbi.nlm.nih.gov/nuccore/NM_001365276.2',
            'protein': 'https://www.ncbi.nlm.nih.gov/nuccore/NP_001352205.1'}

    # TSC1
    def test_variant10(self):
        variant = 'ENST00000298552.9:c.363+1dupG'
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert 'ENST00000298552.9:c.363+1dup' in list(results.keys())
        assert results['ENST00000298552.9:c.363+1dup']['submitted_variant'] == 'ENST00000298552.9:c.363+1dupG'
        assert results['ENST00000298552.9:c.363+1dup']['gene_symbol'] == 'TSC1'
        assert results['ENST00000298552.9:c.363+1dup']['gene_ids'] == {'hgnc_id': 'HGNC:12362', 'entrez_gene_id': '7248',
                                                                 'ucsc_id': 'uc004cca.3', 'omim_id': ['605284']}
        assert results['ENST00000298552.9:c.363+1dup']['hgvs_transcript_variant'] == 'ENST00000298552.9:c.363+1dup'
        assert results['ENST00000298552.9:c.363+1dup'][
                   'genome_context_intronic_sequence'] == 'NC_000009.11(ENST00000298552.9):c.363+1dup'
        # assert results['ENST00000298552.9:c.363+1dup'][
        #            'refseqgene_context_intronic_sequence'] == 'NG_012386.1(ENST00000298552.9):c.363+1dup'
        # assert results['ENST00000298552.9:c.363+1dup']['hgvs_refseqgene_variant'] == 'NG_012386.1:g.24048dup'
        assert results['ENST00000298552.9:c.363+1dup']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000298552.3:p.?', 'slr': 'ENSP00000298552.3:p.?'}
        # assert results['ENST00000298552.9:c.363+1dup']['hgvs_lrg_transcript_variant'] == 'LRG_486t1:c.363+1dup'
        # assert results['ENST00000298552.9:c.363+1dup']['hgvs_lrg_variant'] == 'LRG_486:g.24048dup'
        self.assertCountEqual(results['ENST00000298552.9:c.363+1dup']['alt_genomic_loci'], [])
        assert results['ENST00000298552.9:c.363+1dup']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000009.11:g.135800974dup',
            'vcf': {'chr': 'chr9', 'pos': '135800972', 'ref': 'A', 'alt': 'AC'}}
        assert results['ENST00000298552.9:c.363+1dup']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000009.12:g.132925587dup',
            'vcf': {'chr': 'chr9', 'pos': '132925585', 'ref': 'A', 'alt': 'AC'}}
        assert results['ENST00000298552.9:c.363+1dup']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000009.11:g.135800974dup',
            'vcf': {'chr': '9', 'pos': '135800972', 'ref': 'A', 'alt': 'AC'}}
        assert results['ENST00000298552.9:c.363+1dup']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000009.12:g.132925587dup',
            'vcf': {'chr': '9', 'pos': '132925585', 'ref': 'A', 'alt': 'AC'}}
        assert results['ENST00000298552.9:c.363+1dup']['reference_sequence_records'] == {
            'transcript': 'https://www.ncbi.nlm.nih.gov/nuccore/NM_000368.5',
            'protein': 'https://www.ncbi.nlm.nih.gov/nuccore/NP_000359.1',
            #'refseqgene': 'https://www.ncbi.nlm.nih.gov/nuccore/NG_012386.1',
            #'lrg': 'http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_486.xml'
            }

    # TSC2
    def test_variant11(self):
        variant = 'NC_000016.9:g.2099572TC>T'
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000219476.9:c.138+821del' in list(results.keys())
        assert results['ENST00000219476.9:c.138+821del']['submitted_variant'] == 'NC_000016.9:g.2099572TC>T'
        assert results['ENST00000219476.9:c.138+821del']['gene_symbol'] == 'TSC2'
        assert results['ENST00000219476.9:c.138+821del']['gene_ids'] == {'hgnc_id': 'HGNC:12363', 'entrez_gene_id': '7249',
                                                                   'ucsc_id': 'uc002con.4', 'omim_id': ['191092']}
        assert results['ENST00000219476.9:c.138+821del']['hgvs_transcript_variant'] == 'ENST00000219476.9:c.138+821del'
        assert results['ENST00000219476.9:c.138+821del'][
                   'genome_context_intronic_sequence'] == 'NC_000016.9(ENST00000219476.9):c.138+821del'
        # assert results['ENST00000219476.9:c.138+821del']['refseqgene_context_intronic_sequence'] == ''
        # assert results['ENST00000219476.9:c.138+821del']['hgvs_refseqgene_variant'] == ''
        assert results['ENST00000219476.9:c.138+821del']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000219476.3:p.?', 'slr': 'ENSP00000219476.3:p.?'}
        # assert results['ENST00000219476.9:c.138+821del']['hgvs_lrg_transcript_variant'] == ''
        # assert results['ENST00000219476.9:c.138+821del']['hgvs_lrg_variant'] == ''
        self.assertCountEqual(results['ENST00000219476.9:c.138+821del']['alt_genomic_loci'], [])
        assert results['ENST00000219476.9:c.138+821del']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000016.9:g.2099575del',
            'vcf': {'chr': 'chr16', 'pos': '2099572', 'ref': 'TC', 'alt': 'T'}}
        assert results['ENST00000219476.9:c.138+821del']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000016.10:g.2049574del',
            'vcf': {'chr': 'chr16', 'pos': '2049571', 'ref': 'TC', 'alt': 'T'}}
        assert results['ENST00000219476.9:c.138+821del']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000016.9:g.2099575del',
            'vcf': {'chr': '16', 'pos': '2099572', 'ref': 'TC', 'alt': 'T'}}
        assert results['ENST00000219476.9:c.138+821del']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000016.10:g.2049574del',
            'vcf': {'chr': '16', 'pos': '2049571', 'ref': 'TC', 'alt': 'T'}}
        assert results['ENST00000219476.9:c.138+821del']['reference_sequence_records'] == {
            'transcript': 'https://www.ncbi.nlm.nih.gov/nuccore/NM_000548.5',
            'protein': 'https://www.ncbi.nlm.nih.gov/nuccore/NP_000539.2'}



    

   