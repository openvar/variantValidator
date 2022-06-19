from VariantValidator import Validator
from unittest import TestCase


class TestVariantsEnsembl(TestCase):

    @classmethod
    def setup_class(cls):
        cls.vv = Validator()

    # COL1A1
    def test_variant1(self):
        variant = 'ENST00000225964.10:c.589-1GG>G'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000225964.10:c.590del' in list(results.keys())
        assert results['ENST00000225964.10:c.590del']['submitted_variant'] == 'ENST00000225964.10:c.589-1GG>G'
        assert results['ENST00000225964.10:c.590del']['gene_symbol'] == 'COL1A1'
        assert results['ENST00000225964.10:c.590del']['gene_ids'] == {'hgnc_id': 'HGNC:2197', 'entrez_gene_id': '1277',
                                                               'ucsc_id': 'uc002iqm.4', 'omim_id': ['120150']}
        assert results['ENST00000225964.10:c.590del']['hgvs_transcript_variant'] == 'ENST00000225964.10:c.590del'
        assert results['ENST00000225964.10:c.590del']['genome_context_intronic_sequence'] == ''
        # assert results['NM_000088.3:c.590del']['refseqgene_context_intronic_sequence'] == ''
        # assert results['NM_000088.3:c.590del']['hgvs_refseqgene_variant'] == 'NG_007400.1:g.8639del'
        assert results['ENST00000225964.10:c.590del']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000225964.6:p.(Gly197ValfsTer68)', 'slr': 'ENSP00000225964.6:p.(G197Vfs*68)'}
        # assert results['NM_000088.3:c.590del']['hgvs_lrg_transcript_variant'] == 'LRG_1t1:c.590del'
        # assert results['NM_000088.3:c.590del']['hgvs_lrg_variant'] == 'LRG_1:g.8639del'
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
        variant = 'COL5A1:c.5071A>T'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'warning'
        assert 'validation_warning_1' in list(results.keys())
        assert results['validation_warning_1']['submitted_variant'] == 'COL5A1:c.5071A>T'
        assert results['validation_warning_1']['gene_symbol'] == ''
        assert results['validation_warning_1']['gene_ids'] == {}
        assert results['validation_warning_1']['hgvs_transcript_variant'] == ''
        assert results['validation_warning_1']['genome_context_intronic_sequence'] == ''
        assert results['validation_warning_1']['refseqgene_context_intronic_sequence'] == ''
        assert results['validation_warning_1']['hgvs_refseqgene_variant'] == ''
        assert results['validation_warning_1']['hgvs_predicted_protein_consequence'] == {'tlr': '', 'slr': ''}
        assert results['validation_warning_1']['hgvs_lrg_transcript_variant'] == ''
        assert results['validation_warning_1']['hgvs_lrg_variant'] == ''
        self.assertCountEqual(results['validation_warning_1']['alt_genomic_loci'], [])
        # assert 'hg19' not in list(results['validation_warning_1']['primary_assembly_loci'].keys())
        # # assert 'hg38' not in list(results['validation_warning_1']['primary_assembly_loci'].keys())
        # assert 'grch37' not in list(results['validation_warning_1']['primary_assembly_loci'].keys())
        # assert 'grch38' not in list(results['validation_warning_1']['primary_assembly_loci'].keys())
        assert results['validation_warning_1']['reference_sequence_records'] == ''

    # TP53
    def test_variant3(self):
        variant = 'ENST00000269305.9'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'ENST00000269305.9:c.652_654del' in list(results.keys())
        assert results['ENST00000269305.9:c.652_654del']['submitted_variant'] == 'ENST00000269305.9'
        assert results['ENST00000269305.9:c.652_654del']['gene_symbol'] == 'TP53'
        assert results['ENST00000269305.9:c.652_654del']['gene_ids'] == {'hgnc_id': 'HGNC:11998', 'entrez_gene_id': '7157',
                                                                   'ucsc_id': 'uc060aur.1', 'omim_id': ['191170']}
        assert results['ENST00000269305.9:c.652_654del']['hgvs_transcript_variant'] == 'ENST00000269305.9:c.652_654del'
        assert results['ENST00000269305.9:c.652_654del']['genome_context_intronic_sequence'] == ''
        assert results['ENST00000269305.9:c.652_654del']['refseqgene_context_intronic_sequence'] == ''
        assert results['ENST00000269305.9:c.652_654del']['hgvs_refseqgene_variant'] == 'NG_017013.2:g.17672_17674del'
        assert results['ENST00000269305.9:c.652_654del']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'ENSP00000269305.4:p.(Val218del)', 'slr': 'ENSP00000269305.4:p.(V218del)'}
        assert results['ENST00000269305.9:c.652_654del']['hgvs_lrg_transcript_variant'] == 'LRG_321t1:c.652_654del'
        assert results['ENST00000269305.9:c.652_654del']['hgvs_lrg_variant'] == 'LRG_321:g.17672_17674del'
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

    

   