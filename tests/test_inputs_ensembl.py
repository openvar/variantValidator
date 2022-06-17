from VariantValidator import Validator
from unittest import TestCase

class TestVariantsEnsembl(TestCase):

    @classmethod
    def setup_class(cls):
        cls.vv = Validator()


    # COL1A1
    def test_variant1(self):
        variant = 'NM_000088.3:c.589-1GG>G'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
        print(results)

        assert results['flag'] == 'gene_variant'
        assert 'NM_000088.3:c.590del' in list(results.keys())
        assert results['NM_000088.3:c.590del']['submitted_variant'] == 'NM_000088.3:c.589-1GG>G'
        assert results['NM_000088.3:c.590del']['gene_symbol'] == 'COL1A1'
        assert results['NM_000088.3:c.590del']['gene_ids'] == {'hgnc_id': 'HGNC:2197', 'entrez_gene_id': '1277',
                                                               'ucsc_id': 'uc002iqm.4', 'omim_id': ['120150']}
        assert results['NM_000088.3:c.590del']['hgvs_transcript_variant'] == 'NM_000088.3:c.590del'
        assert results['NM_000088.3:c.590del']['genome_context_intronic_sequence'] == ''
        # assert results['NM_000088.3:c.590del']['refseqgene_context_intronic_sequence'] == ''
        # assert results['NM_000088.3:c.590del']['hgvs_refseqgene_variant'] == 'NG_007400.1:g.8639del'
        assert results['NM_000088.3:c.590del']['hgvs_predicted_protein_consequence'] == {
            'tlr': 'NP_000079.2:p.(Gly197ValfsTer68)', 'slr': 'NP_000079.2:p.(G197Vfs*68)'}
        # assert results['NM_000088.3:c.590del']['hgvs_lrg_transcript_variant'] == 'LRG_1t1:c.590del'
        # assert results['NM_000088.3:c.590del']['hgvs_lrg_variant'] == 'LRG_1:g.8639del'
        self.assertCountEqual(results['NM_000088.3:c.590del']['alt_genomic_loci'], [])
        assert results['NM_000088.3:c.590del']['primary_assembly_loci']['hg19'] == {
            'hgvs_genomic_description': 'NC_000017.10:g.48275364del',
            'vcf': {'chr': 'chr17', 'pos': '48275361', 'ref': 'AC', 'alt': 'A'}}
        assert results['NM_000088.3:c.590del']['primary_assembly_loci']['hg38'] == {
            'hgvs_genomic_description': 'NC_000017.11:g.50198003del',
            'vcf': {'chr': 'chr17', 'pos': '50198000', 'ref': 'AC', 'alt': 'A'}}
        assert results['NM_000088.3:c.590del']['primary_assembly_loci']['grch37'] == {
            'hgvs_genomic_description': 'NC_000017.10:g.48275364del',
            'vcf': {'chr': '17', 'pos': '48275361', 'ref': 'AC', 'alt': 'A'}}
        assert results['NM_000088.3:c.590del']['primary_assembly_loci']['grch38'] == {
            'hgvs_genomic_description': 'NC_000017.11:g.50198003del',
            'vcf': {'chr': '17', 'pos': '50198000', 'ref': 'AC', 'alt': 'A'}}
        assert results['NM_000088.3:c.590del']['reference_sequence_records'] == {
             'transcript': 'https://www.ncbi.nlm.nih.gov/nuccore/NM_000088.3',
             'protein': 'https://www.ncbi.nlm.nih.gov/nuccore/NP_000079.2',
            # 'refseqgene': 'https://www.ncbi.nlm.nih.gov/nuccore/NG_007400.1',
            # 'lrg': 'http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_1.xml'
            }
