from VariantValidator import Validator
from unittest import TestCase


class TestVariantsEnsembl(TestCase):

    @classmethod
    def setup_class(cls):
        cls.vv = Validator()

    # COL1A1
    def test_variant1(self):
        variant = 'ENST00000225964.10:c.589-1GG>G'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
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
        # assert results['ENST00000225964.10:c.590del']['primary_assembly_loci']['hg19'] == {
        #     'hgvs_genomic_description': 'NC_000017.10:g.48275364del',
        #     'vcf': {'chr': 'chr17', 'pos': '48275361', 'ref': 'AC', 'alt': 'A'}}
        assert results['ENST00000225964.10:c.590del']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000017.11:g.50198003del',
            'vcf': {'chr': 'chr17', 'pos': '50198000', 'ref': 'AC', 'alt': 'A'}}
        # assert results['ENST00000225964.10:c.590del']['primary_assembly_loci']['grch37'] == {
        #     'hgvs_genomic_description': 'NC_000017.10:g.48275364del',
        #     'vcf': {'chr': '17', 'pos': '48275361', 'ref': 'AC', 'alt': 'A'}}
        assert results['ENST00000225964.10:c.590del']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000017.11:g.50198003del',
            'vcf': {'chr': '17', 'pos': '50198000', 'ref': 'AC', 'alt': 'A'}}
        assert results['ENST00000225964.10:c.590del']['reference_sequence_records'] == {
             'transcript': 'https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000225964.10',
             'protein': 'https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000225964.6',
            }

    # COL5A1
    def test_variant2(self):
        variant = 'ENST00000371817.8:c.5071A>T'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
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
            'tlr': 'ENSP00000360882.3:p.(Arg1691Ter)', 'slr': 'ENSP00000360882.3:p.(R1691*)'}
        # assert results['ENST00000371817.8:c.5071A>T']['hgvs_lrg_transcript_variant'] == ''
        # assert results['ENST00000371817.8:c.5071A>T']['hgvs_lrg_variant'] == ''
        self.assertCountEqual(results['ENST00000371817.8:c.5071A>T']['alt_genomic_loci'], [])
        # assert results['ENST00000371817.8:c.5071A>T']['primary_assembly_loci']['hg19'] == {
        #     'hgvs_genomic_description': 'NC_000009.11:g.137721825A>T',
        #     'vcf': {'chr': 'chr9', 'pos': '137721825', 'ref': 'A', 'alt': 'T'}}
        assert results['ENST00000371817.8:c.5071A>T']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000009.12:g.134829979A>T', 
            'vcf': {'chr': 'chr9', 'pos': '134829979', 'ref': 'A', 'alt': 'T'}}
        # assert results['ENST00000371817.8:c.5071A>T']['primary_assembly_loci']['grch37'] == {
        #     'hgvs_genomic_description': 'NC_000009.11:g.137721825A>T',
        #     'vcf': {'chr': '9', 'pos': '137721825', 'ref': 'A', 'alt': 'T'}}
        assert results['ENST00000371817.8:c.5071A>T']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000009.12:g.134829979A>T', 
            'vcf': {'chr': '9', 'pos': '134829979', 'ref': 'A', 'alt': 'T'}}
        assert results['ENST00000371817.8:c.5071A>T']['reference_sequence_records'] == {
             'transcript': 'https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000371817.8',
             'protein': 'https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000360882.3'}

    # TP53
    def test_variant3(self):
        variant = 'ENST00000269305.9:c.652_654del'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
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
        # assert results['ENST00000269305.9:c.652_654del']['primary_assembly_loci']['hg19'] == {
        #     'hgvs_genomic_description': 'NC_000017.10:g.7578201_7578203del',
        #     'vcf': {'chr': 'chr17', 'pos': '7578194', 'ref': 'GCAC', 'alt': 'G'}}
        assert results['ENST00000269305.9:c.652_654del']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000017.11:g.7674883_7674885del',
            'vcf': {'chr': 'chr17', 'pos': '7674876', 'ref': 'GCAC', 'alt': 'G'}}
        # assert results['ENST00000269305.9:c.652_654del']['primary_assembly_loci']['grch37'] == {
        #     'hgvs_genomic_description': 'NC_000017.10:g.7578201_7578203del',
        #     'vcf': {'chr': '17', 'pos': '7578194', 'ref': 'GCAC', 'alt': 'G'}}
        assert results['ENST00000269305.9:c.652_654del']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000017.11:g.7674883_7674885del',
            'vcf': {'chr': '17', 'pos': '7674876', 'ref': 'GCAC', 'alt': 'G'}}
        assert results['ENST00000269305.9:c.652_654del']['reference_sequence_records'] == {
            'transcript': 'https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000269305.9',
            'protein': 'https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000269305.4'
            #'refseqgene': 'https://www.ncbi.nlm.nih.gov/nuccore/NG_017013.2',
            #'lrg': 'http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_321.xml'
            }

    # P3H1
    def test_variant4(self):
        variant = 'ENST00000296388.10:c.2055+18G>A'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000296388.10:c.2055+18G>A' in list(results.keys())
        assert results['ENST00000296388.10:c.2055+18G>A']['submitted_variant'] == 'ENST00000296388.10:c.2055+18G>A'
        assert results['ENST00000296388.10:c.2055+18G>A']['gene_symbol'] == 'P3H1'
        assert results['ENST00000296388.10:c.2055+18G>A']['gene_ids'] == {'hgnc_id': 'HGNC:19316', 'entrez_gene_id': '64175',
                                                                   'ucsc_id': '', 'omim_id': ['610339']}
        assert results['ENST00000296388.10:c.2055+18G>A']['hgvs_transcript_variant'] == 'ENST00000296388.10:c.2055+18G>A'
        assert results['ENST00000296388.10:c.2055+18G>A']['genome_context_intronic_sequence'] == 'NC_000001.11(ENST00000296388.10):c.2055+18G>A'
        # assert results['ENST00000296388.10:c.2055+18G>A'][
        #     'refseqgene_context_intronic_sequence'] == 'NG_008123.1(ENST00000296388.10):c.2055+18G>A'
        # assert results['ENST00000296388.10:c.2055+18G>A']['hgvs_refseqgene_variant'] == 'NG_008123.1:g.24831G>A'
        assert results['ENST00000296388.10:c.2055+18G>A']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000296388.5:p.?', 'slr': 'ENSP00000296388.5:p.?'}
        # assert results['ENST00000296388.10:c.2055+18G>A']['hgvs_lrg_transcript_variant'] == 'LRG_5t1:c.2055+18G>A'
        # assert results['ENST00000296388.10:c.2055+18G>A']['hgvs_lrg_variant'] == 'LRG_5:g.24831G>A'
        self.assertCountEqual(results['ENST00000296388.10:c.2055+18G>A']['alt_genomic_loci'], [])
        # assert results['ENST00000296388.10:c.2055+18G>A']['primary_assembly_loci']['hg19'] == {
        #     'hgvs_genomic_description': 'NC_000001.10:g.43212925C>T',
        #     'vcf': {'chr': 'chr1', 'pos': '43212925', 'ref': 'C', 'alt': 'T'}}
        assert results['ENST00000296388.10:c.2055+18G>A']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000001.11:g.42747254C>T',
            'vcf': {'chr': 'chr1', 'pos': '42747254', 'ref': 'C', 'alt': 'T'}}
        # assert results['ENST00000296388.10:c.2055+18G>A']['primary_assembly_loci']['grch37'] == {
        #     'hgvs_genomic_description': 'NC_000001.10:g.43212925C>T',
        #     'vcf': {'chr': '1', 'pos': '43212925', 'ref': 'C', 'alt': 'T'}}
        assert results['ENST00000296388.10:c.2055+18G>A']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000001.11:g.42747254C>T',
            'vcf': {'chr': '1', 'pos': '42747254', 'ref': 'C', 'alt': 'T'}}
        assert results['ENST00000296388.10:c.2055+18G>A']['reference_sequence_records'] == {
            'transcript': 'https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000296388.10',
            'protein': 'https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000296388.5'
            }

    # BRCA1
    def test_variant5(self):
        variant = 'ENST00000357654.9:c.301+1G>C'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        
        assert 'ENST00000357654.9:c.301+1G>C' in list(results.keys())
        assert results['ENST00000357654.9:c.301+1G>C']['submitted_variant'] == 'ENST00000357654.9:c.301+1G>C'
        assert results['ENST00000357654.9:c.301+1G>C']['gene_symbol'] == 'BRCA1'
        assert results['ENST00000357654.9:c.301+1G>C']['gene_ids'] == {'hgnc_id': 'HGNC:1100', 'entrez_gene_id': '672',
                                                                 'ucsc_id': 'uc002ict.4', 'omim_id': ['113705']}
        assert results['ENST00000357654.9:c.301+1G>C']['hgvs_transcript_variant'] == 'ENST00000357654.9:c.301+1G>C'
        assert results['ENST00000357654.9:c.301+1G>C'][
                   'genome_context_intronic_sequence'] == 'NC_000017.11(ENST00000357654.9):c.301+1G>C'
        # assert results['ENST00000357654.9:c.301+1G>C'][
        #           'refseqgene_context_intronic_sequence'] == 'NG_005905.2(ENST00000357654.9):c.301+1G>C'
        # assert results['ENST00000357654.9:c.301+1G>C']['hgvs_refseqgene_variant'] == 'NG_005905.2:g.113117G>C'
        assert results['ENST00000357654.9:c.301+1G>C']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000350283.3:p.?', 'slr': 'ENSP00000350283.3:p.?'}
        # assert results['ENST00000357654.9:c.301+1G>C']['hgvs_lrg_transcript_variant'] == 'LRG_292t1:c.301+1G>C'
        # assert results['ENST00000357654.9:c.301+1G>C']['hgvs_lrg_variant'] == 'LRG_292:g.113117G>C'
        self.assertCountEqual(results['ENST00000357654.9:c.301+1G>C']['alt_genomic_loci'], [])
        # assert results['ENST00000357654.9:c.301+1G>C']['primary_assembly_loci']['hg19'] == {
        #     'hgvs_genomic_description': 'NC_000017.10:g.41256884C>G',
        #     'vcf': {'chr': 'chr17', 'pos': '41256884', 'ref': 'C', 'alt': 'G'}}
        assert results['ENST00000357654.9:c.301+1G>C']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000017.11:g.43104867C>G',
            'vcf': {'chr': 'chr17', 'pos': '43104867', 'ref': 'C', 'alt': 'G'}}
        # assert results['ENST00000357654.9:c.301+1G>C']['primary_assembly_loci']['grch37'] == {
        #     'hgvs_genomic_description': 'NC_000017.10:g.41256884C>G',
        #     'vcf': {'chr': '17', 'pos': '41256884', 'ref': 'C', 'alt': 'G'}}
        assert results['ENST00000357654.9:c.301+1G>C']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000017.11:g.43104867C>G',
            'vcf': {'chr': '17', 'pos': '43104867', 'ref': 'C', 'alt': 'G'}}
        assert results['ENST00000357654.9:c.301+1G>C']['reference_sequence_records'] == {
            'transcript': 'https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000357654.9',
            'protein': 'https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000350283.3'
            }
    
    # BRCA2
    def test_variant6(self):
        variant = 'NC_000013.10:g.32929387T>C'
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000380152.3:c.7397T>C' in list(results.keys())
        assert results['ENST00000380152.3:c.7397T>C']['submitted_variant'] == 'NC_000013.10:g.32929387T>C'
        assert results['ENST00000380152.3:c.7397T>C']['gene_symbol'] == 'BRCA2'
        assert results['ENST00000380152.3:c.7397T>C']['gene_ids'] == {'hgnc_id': 'HGNC:1101', 'entrez_gene_id': '675',
                                                               'ucsc_id': 'uc001uub.2', 'omim_id': ['600185']}
        assert results['ENST00000380152.3:c.7397T>C']['hgvs_transcript_variant'] == 'ENST00000380152.3:c.7397T>C'
        assert results['ENST00000380152.3:c.7397T>C']['genome_context_intronic_sequence'] == ''
        # assert results['ENST00000380152.3:c.7397T>C']['refseqgene_context_intronic_sequence'] == ''
        # assert results['ENST00000380152.3:c.7397T>C']['hgvs_refseqgene_variant'] == 'NG_012772.3:g.44771='
        assert results['ENST00000380152.3:c.7397T>C']['hgvs_predicted_protein_consequence'] == {
            "tlr": "ENSP00000369497.3:p.(Val2466Ala)",
            "slr": "ENSP00000369497.3:p.(V2466A)"
        }
        # assert results['ENST00000380152.3:c.7397T>C']['hgvs_lrg_transcript_variant'] == 'LRG_293t1:c.7397T>C'
        # assert results['ENST00000380152.3:c.7397T>C']['hgvs_lrg_variant'] == 'LRG_293:g.44771='
        self.assertCountEqual(results['ENST00000380152.3:c.7397T>C']['alt_genomic_loci'], [])
        assert results['ENST00000380152.3:c.7397T>C']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000013.10:g.32929387T>C',
            'vcf': {'chr': 'chr13', 'pos': '32929387', 'ref': 'T', 'alt': 'C'}}
        # assert results['ENST00000380152.3:c.7397T>C']['primary_assembly_loci']['hg38'] == {
        #     'hgvs_genomic_description': 'NC_000013.11:g.32355250T>C',
        #     'vcf': {'chr': 'chr13', 'pos': '32355250', 'ref': 'T', 'alt': 'C'}}
        assert results['ENST00000380152.3:c.7397T>C']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000013.10:g.32929387T>C',
            'vcf': {'chr': '13', 'pos': '32929387', 'ref': 'T', 'alt': 'C'}}
        # assert results['ENST00000380152.3:c.7397T>C']['primary_assembly_loci']['grch38'] == {
        #     'hgvs_genomic_description': 'NC_000013.11:g.32355250T>C',
        #     'vcf': {'chr': '13', 'pos': '32355250', 'ref': 'T', 'alt': 'C'}}

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000544455.1:c.7397T>C' in list(results.keys())
        assert results['ENST00000544455.1:c.7397T>C']['submitted_variant'] == 'NC_000013.10:g.32929387T>C'
        assert results['ENST00000544455.1:c.7397T>C']['gene_symbol'] == 'BRCA2'
        assert results['ENST00000544455.1:c.7397T>C']['gene_ids'] == {'hgnc_id': 'HGNC:1101', 'entrez_gene_id': '675',
                                                               'ucsc_id': 'uc001uub.2', 'omim_id': ['600185']}
        assert results['ENST00000544455.1:c.7397T>C']['hgvs_transcript_variant'] == 'ENST00000544455.1:c.7397T>C'
        assert results['ENST00000544455.1:c.7397T>C']['genome_context_intronic_sequence'] == ''
        # assert results['ENST00000544455.1:c.7397T>C']['refseqgene_context_intronic_sequence'] == ''
        # assert results['ENST00000544455.1:c.7397T>C']['hgvs_refseqgene_variant'] == 'NG_012772.3:g.44771='
        assert results['ENST00000544455.1:c.7397T>C']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000439902.1:p.(Val2466Ala)', 'slr': 'ENSP00000439902.1:p.(V2466A)'}
        # assert results['ENST00000544455.1:c.7397T>C']['hgvs_lrg_transcript_variant'] == 'LRG_293t1:c.7397T>C'
        # assert results['ENST00000544455.1:c.7397T>C']['hgvs_lrg_variant'] == 'LRG_293:g.44771='
        self.assertCountEqual(results['ENST00000544455.1:c.7397T>C']['alt_genomic_loci'], [])
        assert results['ENST00000544455.1:c.7397T>C']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000013.10:g.32929387T>C',
            'vcf': {'chr': 'chr13', 'pos': '32929387', 'ref': 'T', 'alt': 'C'}}
        # assert results['ENST00000544455.1:c.7397T>C']['primary_assembly_loci']['hg38'] == {
        #     'hgvs_genomic_description': 'NC_000013.11:g.32355250T>C',
        #     'vcf': {'chr': 'chr13', 'pos': '32355250', 'ref': 'T', 'alt': 'C'}}
        assert results['ENST00000544455.1:c.7397T>C']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000013.10:g.32929387T>C',
            'vcf': {'chr': '13', 'pos': '32929387', 'ref': 'T', 'alt': 'C'}}
        # assert results['ENST00000544455.1:c.7397T>C']['primary_assembly_loci']['grch38'] == {
        #     'hgvs_genomic_description': 'NC_000013.11:g.32355250T>C',
        #     'vcf': {'chr': '13', 'pos': '32355250', 'ref': 'T', 'alt': 'C'}}

    # HBG1
    def test_variant7(self):
        variant = '11-5248232-T-A' # Pseudo-VCF format
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000330597.3:c.*127A>T' in list(results.keys())
        assert results['ENST00000330597.3:c.*127A>T']['submitted_variant'] == '11-5248232-T-A'
        assert results['ENST00000330597.3:c.*127A>T']['gene_symbol'] == 'HBG1'
        assert results['ENST00000330597.3:c.*127A>T']['gene_ids'] == {'hgnc_id': 'HGNC:4831', 'entrez_gene_id': '3047',
                                                              'ucsc_id': 'uc001mah.2', 'omim_id': ['142200']}
        assert results['ENST00000330597.3:c.*127A>T']['hgvs_transcript_variant'] == 'ENST00000330597.3:c.*127A>T'
        assert results['ENST00000330597.3:c.*127A>T']['genome_context_intronic_sequence'] == ''
        # assert results['ENST00000330597.3:c.*127A>T']['refseqgene_context_intronic_sequence'] == ''
        # assert results['ENST00000330597.3:c.*127A>T']['hgvs_refseqgene_variant'] == 'NG_059281.1:g.5070A>T'
        assert results['ENST00000330597.3:c.*127A>T']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000327431.3:p.?', 'slr': 'ENSP00000327431.3:p.?'}
        # assert results['ENST00000330597.3:c.*127A>T']['hgvs_lrg_transcript_variant'] == 'LRG_1232t1:c.20A>T'
        # assert results['ENST00000330597.3:c.*127A>T']['hgvs_lrg_variant'] == 'LRG_1232:g.5070A>T'
        self.assertCountEqual(results['ENST00000330597.3:c.*127A>T']['alt_genomic_loci'], [])
        # assert results['ENST00000330597.3:c.*127A>T']['primary_assembly_loci']['hg19'] == {
        #     'hgvs_genomic_description': 'NC_000011.9:g.5269462T>A',
        #     'vcf': {'chr': 'chr11', 'pos': '5269462', 'ref': 'T', 'alt': 'A'}}
        assert results['ENST00000330597.3:c.*127A>T']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000011.10:g.5248232T>A', 
            'vcf': {'chr': 'chr11', 'pos': '5248232', 'ref': 'T', 'alt': 'A'}}
        # assert results['ENST00000330597.3:c.*127A>T']['primary_assembly_loci']['grch37'] == {
        #     'hgvs_genomic_description': 'NC_000011.9:g.5269462T>A',
        #     'vcf': {'chr': '11', 'pos': '5269462', 'ref': 'T', 'alt': 'A'}}
        assert results['ENST00000330597.3:c.*127A>T']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000011.10:g.5248232T>A', 
            'vcf': {'chr': '11', 'pos': '5248232', 'ref': 'T', 'alt': 'A'}}
        assert results['ENST00000330597.3:c.*127A>T']['reference_sequence_records'] == {
            'transcript': 'https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000330597.3',
            'protein': 'https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000327431.3'
            }

    # TNXB
    def test_variant9(self):
        variant = '6-32012992-CG-C'
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000375247.2:c.10711del' in list(results.keys())
        assert results['ENST00000375247.2:c.10711del']['submitted_variant'] == '6-32012992-CG-C'
        assert results['ENST00000375247.2:c.10711del']['gene_symbol'] == 'TNXB'
        assert results['ENST00000375247.2:c.10711del']['gene_ids'] == {'hgnc_id': 'HGNC:11976', 'entrez_gene_id': '7148',
                                                                 'ucsc_id': 'uc063nnw.1', 'omim_id': ['600985']}
        assert results['ENST00000375247.2:c.10711del']['hgvs_transcript_variant'] == 'ENST00000375247.2:c.10711del'
        assert results['ENST00000375247.2:c.10711del']['genome_context_intronic_sequence'] == ''
        #assert results['ENST00000375247.2:c.10711del']['refseqgene_context_intronic_sequence'] == ''
        #assert results['ENST00000375247.2:c.10711del']['hgvs_refseqgene_variant'] == 'NG_008337.2:g.69159del'
        assert results['ENST00000375247.2:c.10711del']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000364396.2:p.(Arg3571AlafsTer91)', 'slr': 'ENSP00000364396.2:p.(R3571Afs*91)'}
        #assert results['ENST00000375247.2:c.10711del']['hgvs_lrg_transcript_variant'] == ''
        #assert results['ENST00000375247.2:c.10711del']['hgvs_lrg_variant'] == ''
        assert results['ENST00000375247.2:c.10711del']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000006.11:g.32012993del',
            'vcf': {'chr': 'chr6', 'pos': '32012992', 'ref': 'CG', 'alt': 'C'}}
        # assert results['ENST00000375247.2:c.10711del']['primary_assembly_loci']['hg38'] == {
        #     'hgvs_genomic_description': 'NC_000006.12:g.32045216del',
        #     'vcf': {'chr': 'chr6', 'pos': '32045215', 'ref': 'CG', 'alt': 'C'}}
        assert results['ENST00000375247.2:c.10711del']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000006.11:g.32012993del',
            'vcf': {'chr': '6', 'pos': '32012992', 'ref': 'CG', 'alt': 'C'}}
        # assert results['ENST00000375247.2:c.10711del']['primary_assembly_loci']['grch38'] == {
        #     'hgvs_genomic_description': 'NC_000006.12:g.32045216del',
        #     'vcf': {'chr': '6', 'pos': '32045215', 'ref': 'CG', 'alt': 'C'}}
        assert results['ENST00000375247.2:c.10711del']['reference_sequence_records'] == {
            'transcript': 'https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000375247.2',
            'protein': 'https://grch37.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000364396.2'
            #'refseqgene': 'https://www.ncbi.nlm.nih.gov/nuccore/NG_008337.2'
            }

        assert 'ENST00000451343.1:c.4del' in list(results.keys())
        assert results['ENST00000451343.1:c.4del']['submitted_variant'] == '6-32012992-CG-C'
        assert results['ENST00000451343.1:c.4del']['gene_symbol'] == 'TNXB'
        assert results['ENST00000451343.1:c.4del']['gene_ids'] == {'hgnc_id': 'HGNC:11976', 'entrez_gene_id': '7148',
                                                             'ucsc_id': 'uc063nnw.1', 'omim_id': ['600985']}
        assert results['ENST00000451343.1:c.4del']['hgvs_transcript_variant'] == 'ENST00000451343.1:c.4del'
        assert results['ENST00000451343.1:c.4del']['genome_context_intronic_sequence'] == ''
        #assert results['ENST00000451343.1:c.4del']['refseqgene_context_intronic_sequence'] == ''
        #assert results['ENST00000451343.1:c.4del']['hgvs_refseqgene_variant'] == 'NG_008337.2:g.69159del'
        assert results['ENST00000451343.1:c.4del']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000407685.1:p.(Arg2AlafsTer91)', 'slr': 'ENSP00000407685.1:p.(R2Afs*91)'}
        #assert results['ENST00000451343.1:c.4del']['hgvs_lrg_transcript_variant'] == ''
        #assert results['ENST00000451343.1:c.4del']['hgvs_lrg_variant'] == ''
        assert results['ENST00000451343.1:c.4del']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000006.11:g.32012993del',
            'vcf': {'chr': 'chr6', 'pos': '32012992', 'ref': 'CG', 'alt': 'C'}}
        # assert results['ENST00000451343.1:c.4del']['primary_assembly_loci']['hg38'] == {
        #     'hgvs_genomic_description': 'NC_000006.12:g.32045216del',
        #     'vcf': {'chr': 'chr6', 'pos': '32045215', 'ref': 'CG', 'alt': 'C'}}
        assert results['ENST00000451343.1:c.4del']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000006.11:g.32012993del',
            'vcf': {'chr': '6', 'pos': '32012992', 'ref': 'CG', 'alt': 'C'}}
        # assert results['ENST00000451343.1:c.4del']['primary_assembly_loci']['grch38'] == {
        #     'hgvs_genomic_description': 'NC_000006.12:g.32045216del',
        #     'vcf': {'chr': '6', 'pos': '32045215', 'ref': 'CG', 'alt': 'C'}}
        assert results['ENST00000451343.1:c.4del']['reference_sequence_records'] == {
            'transcript': 'https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000451343.1',
            'protein': 'https://grch37.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000407685.1'
            }

    # TSC1
    def test_variant10(self):
        variant = 'ENST00000298552.9:c.363+1dupG'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert 'ENST00000298552.9:c.363+1dup' in list(results.keys())
        assert results['ENST00000298552.9:c.363+1dup']['submitted_variant'] == 'ENST00000298552.9:c.363+1dupG'
        assert results['ENST00000298552.9:c.363+1dup']['gene_symbol'] == 'TSC1'
        assert results['ENST00000298552.9:c.363+1dup']['gene_ids'] == {'hgnc_id': 'HGNC:12362', 'entrez_gene_id': '7248',
                                                                 'ucsc_id': 'uc004cca.3', 'omim_id': ['605284']}
        assert results['ENST00000298552.9:c.363+1dup']['hgvs_transcript_variant'] == 'ENST00000298552.9:c.363+1dup'
        assert results['ENST00000298552.9:c.363+1dup'][
                   'genome_context_intronic_sequence'] == 'NC_000009.12(ENST00000298552.9):c.363+1dup'
        # assert results['ENST00000298552.9:c.363+1dup'][
        #            'refseqgene_context_intronic_sequence'] == 'NG_012386.1(ENST00000298552.9):c.363+1dup'
        # assert results['ENST00000298552.9:c.363+1dup']['hgvs_refseqgene_variant'] == 'NG_012386.1:g.24048dup'
        assert results['ENST00000298552.9:c.363+1dup']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000298552.3:p.?', 'slr': 'ENSP00000298552.3:p.?'}
        # assert results['ENST00000298552.9:c.363+1dup']['hgvs_lrg_transcript_variant'] == 'LRG_486t1:c.363+1dup'
        # assert results['ENST00000298552.9:c.363+1dup']['hgvs_lrg_variant'] == 'LRG_486:g.24048dup'
        self.assertCountEqual(results['ENST00000298552.9:c.363+1dup']['alt_genomic_loci'], [])
        # assert results['ENST00000298552.9:c.363+1dup']['primary_assembly_loci']['hg19'] == {
        #     'hgvs_genomic_description': 'NC_000009.11:g.135800974dup',
        #     'vcf': {'chr': 'chr9', 'pos': '135800972', 'ref': 'A', 'alt': 'AC'}}
        assert results['ENST00000298552.9:c.363+1dup']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000009.12:g.132925587dup',
            'vcf': {'chr': 'chr9', 'pos': '132925585', 'ref': 'A', 'alt': 'AC'}}
        # assert results['ENST00000298552.9:c.363+1dup']['primary_assembly_loci']['grch37'] == {
        #     'hgvs_genomic_description': 'NC_000009.11:g.135800974dup',
        #     'vcf': {'chr': '9', 'pos': '135800972', 'ref': 'A', 'alt': 'AC'}}
        assert results['ENST00000298552.9:c.363+1dup']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000009.12:g.132925587dup',
            'vcf': {'chr': '9', 'pos': '132925585', 'ref': 'A', 'alt': 'AC'}}
        assert results['ENST00000298552.9:c.363+1dup']['reference_sequence_records'] == {
            'transcript': 'https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000298552.9',
            'protein': 'https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000298552.3'
            }

    # TSC2
    def test_variant11(self):
        variant = 'NC_000016.10:g.2099572TG>T'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000262304.9:c.10050+71del' in list(results.keys())
        assert results['ENST00000262304.9:c.10050+71del']['submitted_variant'] == 'NC_000016.10:g.2099572TG>T'
        assert results['ENST00000262304.9:c.10050+71del']['gene_symbol'] == 'PKD1'
        assert results['ENST00000262304.9:c.10050+71del']['hgvs_transcript_variant'
               ] == ('ENST00000262304.9:c.10050+71del')
        assert results['ENST00000262304.9:c.10050+71del'][
                   'genome_context_intronic_sequence'] == 'NC_000016.10(ENST00000262304.9):c.10050+71del'
        assert results['ENST00000262304.9:c.10050+71del']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000262304.4:p.?', 'slr': 'ENSP00000262304.4:p.?'}
        self.assertCountEqual(results['ENST00000262304.9:c.10050+71del']['alt_genomic_loci'], [])
        # assert results['ENST00000262304.9:c.10050+71del']['primary_assembly_loci']['hg19'] == {
        #     'hgvs_genomic_description': 'NC_000016.9:g.2149575del', 'vcf': {
        #         'chr': 'chr16', 'pos': '2149573', 'ref': 'TG', 'alt': 'T'}}
        assert results['ENST00000262304.9:c.10050+71del']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000016.10:g.2099574del', 'vcf': {
                'chr': 'chr16', 'pos': '2099572', 'ref': 'TG', 'alt': 'T'}}

    """
    Gapped alignment test variants
    """

    def test_variant12(self):
        variant = '19-41123094-G-GG'  # ENST00000396819.3 contains 1 extra bases between c.3233_3235 than NC_000019.9
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000396819.3:c.3033_3034insGGT' in list(results.keys())
        assert results['ENST00000396819.3:c.3033_3034insGGT']['submitted_variant'] == '19-41123094-G-GG'
        assert results['ENST00000396819.3:c.3033_3034insGGT']['gene_symbol'] == 'LTBP4'
        assert results['ENST00000396819.3:c.3033_3034insGGT']['gene_ids'] == {'hgnc_id': 'HGNC:6717', 'entrez_gene_id': '8425',
                                                                      'ucsc_id': 'uc032hxp.2', 'omim_id': ['604710']}
        assert results['ENST00000396819.3:c.3033_3034insGGT']['hgvs_transcript_variant'] == 'ENST00000396819.3:c.3033_3034insGGT'
        assert results['ENST00000396819.3:c.3033_3034insGGT']['genome_context_intronic_sequence'] == ''
        assert results['ENST00000396819.3:c.3033_3034insGGT']['refseqgene_context_intronic_sequence'] == ''
        assert results['ENST00000396819.3:c.3033_3034insGGT']['hgvs_predicted_protein_consequence'] == \
               {
                   "slr": "ENSP00000380031.3:p.(Q1011_Y1012insG)",
                   "tlr": "ENSP00000380031.3:p.(Gln1011_Tyr1012insGly)"
               }
        assert results['ENST00000396819.3:c.3033_3034insGGT']['hgvs_lrg_transcript_variant'] == ''
        assert results['ENST00000396819.3:c.3033_3034insGGT']['hgvs_lrg_variant'] == ''
        self.assertCountEqual(results['ENST00000396819.3:c.3033_3034insGGT']['alt_genomic_loci'], [])
        assert results['ENST00000396819.3:c.3033_3034insGGT']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000019.9:g.41123095dup',
            'vcf': {'chr': 'chr19', 'pos': '41123093', 'ref': 'A', 'alt': 'AG'}}
        # assert results['ENST00000396819.3:c.3033_3034insGGT']['primary_assembly_loci']['hg38'] == {
        #     'hgvs_genomic_description': 'NC_000019.10:g.40617187_40617189=',
        #     'vcf': {'chr': 'chr19', 'pos': '40617187', 'ref': 'AGG', 'alt': 'AGG'}}
        assert results['ENST00000396819.3:c.3033_3034insGGT']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000019.9:g.41123095dup',
            'vcf': {'chr': '19', 'pos': '41123093', 'ref': 'A', 'alt': 'AG'}}
        # assert results['ENST00000396819.3:c.3033_3034insGGT']['primary_assembly_loci']['grch38'] == {
        #     'hgvs_genomic_description': 'NC_000019.10:g.40617187_40617189=',
        #     'vcf': {'chr': '19', 'pos': '40617187', 'ref': 'AGG', 'alt': 'AGG'}}

    def test_variant12b(self):
        variant = 'ENST00000396819.3:c.3033_3034insGGT'  # ENST00000396819.3 contains 1 extra bases between c.3233_3235 than NC_000019.9
        results = self.vv.validate(variant, 'GRCh37', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000396819.3:c.3033_3034insGGT' in list(results.keys())
        assert results['ENST00000396819.3:c.3033_3034insGGT']['submitted_variant'] == 'ENST00000396819.3:c.3033_3034insGGT'
        assert results['ENST00000396819.3:c.3033_3034insGGT']['gene_symbol'] == 'LTBP4'
        assert results['ENST00000396819.3:c.3033_3034insGGT']['gene_ids'] == {'hgnc_id': 'HGNC:6717', 'entrez_gene_id': '8425',
                                                                      'ucsc_id': 'uc032hxp.2', 'omim_id': ['604710']}
        assert results['ENST00000396819.3:c.3033_3034insGGT']['hgvs_transcript_variant'] == 'ENST00000396819.3:c.3033_3034insGGT'
        assert results['ENST00000396819.3:c.3033_3034insGGT']['genome_context_intronic_sequence'] == ''
        assert results['ENST00000396819.3:c.3033_3034insGGT']['refseqgene_context_intronic_sequence'] == ''
        assert results['ENST00000396819.3:c.3033_3034insGGT']['hgvs_predicted_protein_consequence'] == \
               {
                   "slr": "ENSP00000380031.3:p.(Q1011_Y1012insG)",
                   "tlr": "ENSP00000380031.3:p.(Gln1011_Tyr1012insGly)"
               }
        assert results['ENST00000396819.3:c.3033_3034insGGT']['hgvs_lrg_transcript_variant'] == ''
        assert results['ENST00000396819.3:c.3033_3034insGGT']['hgvs_lrg_variant'] == ''
        self.assertCountEqual(results['ENST00000396819.3:c.3033_3034insGGT']['alt_genomic_loci'], [])
        assert results['ENST00000396819.3:c.3033_3034insGGT']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000019.9:g.41123095dup',
            'vcf': {'chr': 'chr19', 'pos': '41123093', 'ref': 'A', 'alt': 'AG'}}
        # assert results['ENST00000396819.3:c.3033_3034insGGT']['primary_assembly_loci']['hg38'] == {
        #     'hgvs_genomic_description': 'NC_000019.10:g.40617187_40617189=',
        #     'vcf': {'chr': 'chr19', 'pos': '40617187', 'ref': 'AGG', 'alt': 'AGG'}}
        assert results['ENST00000396819.3:c.3033_3034insGGT']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000019.9:g.41123095dup',
            'vcf': {'chr': '19', 'pos': '41123093', 'ref': 'A', 'alt': 'AG'}}
        # assert results['ENST00000396819.3:c.3033_3034insGGT']['primary_assembly_loci']['grch38'] == {
        #     'hgvs_genomic_description': 'NC_000019.10:g.40617187_40617189=',
        #     'vcf': {'chr': '19', 'pos': '40617187', 'ref': 'AGG', 'alt': 'AGG'}}

    # 
    def test_variant13(self):
        # No GRCh37 coding transcripts for NR2E3
        variant = '15-72105928-AC-A'  # ENST00000398840.2 contains 1 fewer bases between NC_000015.9
        results = self.vv.validate(variant, 'GRCh37', 'ENST00000398840.2', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000398840.2:n.1133_1141=' in list(results.keys())
        assert results['ENST00000398840.2:n.1133_1141=']['submitted_variant'] == '15-72105928-AC-A'
        assert results['ENST00000398840.2:n.1133_1141=']['gene_symbol'] == 'NR2E3'
        assert results['ENST00000398840.2:n.1133_1141=']['hgvs_transcript_variant'] == 'ENST00000398840.2:n.1133_1141='
        assert results['ENST00000398840.2:n.1133_1141='][
                   'genome_context_intronic_sequence'] == ''
        # assert results['ENST00000398840.2:n.1133_1141=']['refseqgene_context_intronic_sequence'] == ''
        # assert results['ENST00000398840.2:n.1133_1141=']['hgvs_refseqgene_variant'] == ''
        assert results['ENST00000398840.2:n.1133_1141=']['hgvs_predicted_protein_consequence'] == {
            'tlr': '', 'slr': ''}
        # assert results['ENST00000398840.2:n.1133_1141=']['hgvs_lrg_transcript_variant'] == ''
        # assert results['ENST00000398840.2:n.1133_1141=']['hgvs_lrg_variant'] == ''
        self.assertCountEqual(results['ENST00000398840.2:n.1133_1141=']['alt_genomic_loci'], [])
        assert results['ENST00000398840.2:n.1133_1141=']['primary_assembly_loci']['hg19'] == {
                "hgvs_genomic_description": "NC_000015.9:g.72105933del",
                "vcf": {
                    "alt": "A",
                    "chr": "chr15",
                    "pos": "72105928",
                    "ref": "AC"
                }}
        # assert results['ENST00000398840.2:n.1133_1141=']['primary_assembly_loci']['hg38'] == {
        #     'hgvs_genomic_description': 'NC_000016.10:g.2049574del',
        #     'vcf': {'chr': 'chr16', 'pos': '2049571', 'ref': 'TC', 'alt': 'T'}}
        assert results['ENST00000398840.2:n.1133_1141=']['primary_assembly_loci']['grch37'] == {
                "hgvs_genomic_description": "NC_000015.9:g.72105933del",
                "vcf": {
                    "alt": "A",
                    "chr": "15",
                    "pos": "72105928",
                    "ref": "AC"
                }}
        # assert results['ENST00000398840.2:n.1133_1141=']['primary_assembly_loci']['grch38'] == {
        #     'hgvs_genomic_description': 'NC_000016.10:g.2049574del',
        #     'vcf': {'chr': '16', 'pos': '2049571', 'ref': 'TC', 'alt': 'T'}}

    def test_variant13b(self):
        # No GRCh37 coding transcripts for NR2E3
        variant = 'ENST00000398840.2:n.1133_1141='  # ENST00000398840.2 contains 1 fewer bases between NC_000015.9
        results = self.vv.validate(variant, 'GRCh37', 'ENST00000398840.2', transcript_set="ensembl").format_as_dict(
            test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000398840.2:n.1133_1141=' in list(results.keys())
        assert results['ENST00000398840.2:n.1133_1141=']['submitted_variant'] == 'ENST00000398840.2:n.1133_1141='
        assert results['ENST00000398840.2:n.1133_1141=']['gene_symbol'] == 'NR2E3'
        assert results['ENST00000398840.2:n.1133_1141=']['hgvs_transcript_variant'] == 'ENST00000398840.2:n.1133_1141='
        assert results['ENST00000398840.2:n.1133_1141='][
                   'genome_context_intronic_sequence'] == ''
        # assert results['ENST00000398840.2:n.1133_1141=']['refseqgene_context_intronic_sequence'] == ''
        # assert results['ENST00000398840.2:n.1133_1141=']['hgvs_refseqgene_variant'] == ''
        assert results['ENST00000398840.2:n.1133_1141=']['hgvs_predicted_protein_consequence'] == {
            'tlr': '', 'slr': ''}
        # assert results['ENST00000398840.2:n.1133_1141=']['hgvs_lrg_transcript_variant'] == ''
        # assert results['ENST00000398840.2:n.1133_1141=']['hgvs_lrg_variant'] == ''
        self.assertCountEqual(results['ENST00000398840.2:n.1133_1141=']['alt_genomic_loci'], [])
        assert results['ENST00000398840.2:n.1133_1141=']['primary_assembly_loci']['hg19'] == {
            "hgvs_genomic_description": "NC_000015.9:g.72105933del",
            "vcf": {
                "alt": "A",
                "chr": "chr15",
                "pos": "72105928",
                "ref": "AC"
            }}
        # assert results['ENST00000398840.2:n.1133_1141=']['primary_assembly_loci']['hg38'] == {
        #     'hgvs_genomic_description': 'NC_000016.10:g.2049574del',
        #     'vcf': {'chr': 'chr16', 'pos': '2049571', 'ref': 'TC', 'alt': 'T'}}
        assert results['ENST00000398840.2:n.1133_1141=']['primary_assembly_loci']['grch37'] == {
            "hgvs_genomic_description": "NC_000015.9:g.72105933del",
            "vcf": {
                "alt": "A",
                "chr": "15",
                "pos": "72105928",
                "ref": "AC"
            }}
        # assert results['ENST00000398840.2:n.1133_1141=']['primary_assembly_loci']['grch38'] == {
        #     'hgvs_genomic_description': 'NC_000016.10:g.2049574del',
        #     'vcf': {'chr': '16', 'pos': '2049571', 'ref': 'TC', 'alt': 'T'}}

    def test_variant14(self):
        # Because this is a 3 nt deletion, there is no gap in the transcript becuse the transcript matches the genome
        variant = 'NC_000002.11:g.95847041_95847043GCG='  # ENST00000340539.5 contains 3 fewer bases between NC_000002.11
        results = self.vv.validate(variant, 'GRCh37', 'ENST00000340539.5', transcript_set="ensembl").format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000340539.5:c.468_470=' in list(results.keys())
        assert results['ENST00000340539.5:c.468_470=']['submitted_variant'] == 'NC_000002.11:g.95847041_95847043GCG='
        assert results['ENST00000340539.5:c.468_470=']['gene_symbol'] == 'ZNF2'
        assert results['ENST00000340539.5:c.468_470=']['hgvs_transcript_variant'] == 'ENST00000340539.5:c.468_470='
        assert results['ENST00000340539.5:c.468_470='][
                   'genome_context_intronic_sequence'] == ''
        # assert results['ENST00000340539.5:c.468_470=']['refseqgene_context_intronic_sequence'] == ''
        # assert results['ENST00000340539.5:c.468_470=']['hgvs_refseqgene_variant'] == ''
        assert results['ENST00000340539.5:c.468_470=']['hgvs_predicted_protein_consequence'] == {
            "slr": "ENSP00000345392.5:p.(L156_R157=)",
            "tlr": "ENSP00000345392.5:p.(Leu156_Arg157=)"
        }
        # assert results['ENST00000340539.5:c.468_470=']['hgvs_lrg_transcript_variant'] == ''
        # assert results['ENST00000340539.5:c.468_470=']['hgvs_lrg_variant'] == ''
        self.assertCountEqual(results['ENST00000340539.5:c.468_470=']['alt_genomic_loci'], [])
        assert results['ENST00000340539.5:c.468_470=']['primary_assembly_loci'][
                   'hg19']["hgvs_genomic_description"] == "NC_000002.11:g.95847041_95847043="
        # assert results['ENST00000340539.5:c.468_470=']['primary_assembly_loci']['hg38'] == {
        #     'hgvs_genomic_description': 'NC_000016.10:g.2049574del',
        #     'vcf': {'chr': 'chr16', 'pos': '2049571', 'ref': 'TC', 'alt': 'T'}}
        assert results['ENST00000340539.5:c.468_470=']['primary_assembly_loci'][
                   'grch37']["hgvs_genomic_description"] == "NC_000002.11:g.95847041_95847043="
        # assert results['ENST00000340539.5:c.468_470=']['primary_assembly_loci']['grch38'] == {
        #     'hgvs_genomic_description': 'NC_000016.10:g.2049574del',
        #     'vcf': {'chr': '16', 'pos': '2049571', 'ref': 'TC', 'alt': 'T'}}


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
