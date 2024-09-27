import VariantValidator
from VariantValidator import Validator
from unittest import TestCase
vfo = VariantValidator.Validator()


class TestWarnings(TestCase):

    @classmethod
    def setup_class(cls):
        cls.vv = Validator()
        cls.vv.testing = True

    def test_t_in_rna_string(self):
        variant = 'ENST00000376372.9:r.235_236insGCCCACCCACCTGCCAG'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'The IUPAC RNA alphabet dictates that RNA variants must use the character u in place of t' in \
               results['validation_warning_1']['validation_warnings']

    def test_issue_169(self):
        variant = 'NC_000017.11(NC_000017.11(ENST00000357654.9):c.4422-63C>G'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'ExonBoundaryError: Position c.4422-63 does not correspond with an exon boundary for transcript ENST00000357654.9' in \
               results['validation_warning_1']['validation_warnings']

    def test_issue_176(self):
        variant = 'NC_000023.11(ENST00000357033.9):c.8810A>G' # Unlike the RefSeq test ENST seq always matches the genome (in theory)
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'ENST00000357033.9:c.8810A>G: Variant reference (A) does not agree with reference sequence (G)' in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_issue_180a(self):
        variant = 'NC_000017.10:g.41232400_41236235del383'
        results = self.vv.validate(variant, 'hg19', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'Length implied by coordinates must equal sequence deletion length' in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_issue_180b(self):
        variant = 'NC_000017.11(ENST00000357654.9):c.4186-1642_4358-983del10'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'Length implied by coordinates must equal sequence deletion length' in \
               results['validation_warning_1']['validation_warnings']

    def test_issue_180b1(self):
        variant = 'NC_000017.11(ENST00000357654.9):c.4176-1642_4368-983del10'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert ('ExonBoundaryError: Position c.4176-1642 does not correspond with an exon boundary for transcript '
                'ENST00000357654.9') in \
               results['validation_warning_1']['validation_warnings']

    def test_issue_180c(self):
        variant = 'NC_000017.11(ENST00000225964.10):c.589-1del2'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'Length implied by coordinates must equal sequence deletion length' in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_issue_195a(self):
        variant = 'ENST00000225964.10(COL1A1):c.590delG'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'Removing redundant gene symbol COL1A1 from variant description' in \
               results['ENST00000225964.10:c.590del']['validation_warnings'][0]
        assert 'Removing redundant reference bases from variant description' in \
               results['ENST00000225964.10:c.590del']['validation_warnings'][1]

    def test_issue_216a(self):
        variant = 'ENST00000396884.8:c.850_877dup27'
        results = self.vv.validate(variant, 'hg38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'Length implied by coordinates must equal sequence duplication length' in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_issue_216b(self):
        variant = 'ENST00000396884.8:c.850_877dup28'
        results = self.vv.validate(variant, 'hg38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'Trailing digits are not permitted in HGVS variant descriptions' in \
               results['ENST00000396884.8:c.850_877dup']['validation_warnings'][0]

    def test_issue_239(self):
        variant = 'ENST00000396884.8:c.1047dupT'
        results = self.vv.validate(variant, 'hg38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'Removing redundant reference bases from variant description' in \
               results['ENST00000396884.8:c.1047dup']['validation_warnings'][0]

    def test_issue_338(self):
        # Also issue 357
        variant = 'ENST00000225964.10:C.589G>T'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert ('Reference type incorrectly stated in the variant description ENST00000225964.10:C.589G>T Valid '
                'types are g,c,n,r, or p') in \
               results['ENST00000225964.10:c.589G>T']['validation_warnings'][0]

        variant = 'enst00000225964.10:c.589G>T'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'characters being in the wrong case' in \
               results['ENST00000225964.10:c.589G>T']['validation_warnings'][0]


    def test_issue_359(self):
        variant = 'ENST00000396884.8:c.483ins'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'The inserted sequence must be provided for insertions or deletion-insertions' in \
               results['validation_warning_1']['validation_warnings'][0]
        assert 'An insertion must be provided with the two positions between which the insertion has taken place' in \
               results['validation_warning_1']['validation_warnings'][1]

    def test_issue_359_b(self):
        variant = 'ENST00000396884.8:c.483ins(10)'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'The length of the variant is not formatted following the HGVS guidelines. Please rewrite e.g. (10) ' \
               'to N[10]' in \
               results['validation_warning_1']['validation_warnings'][0]
        assert 'An insertion must be provided with the two positions between which the insertion has taken place' in \
               results['validation_warning_1']['validation_warnings'][1]

    def test_issue_359_c(self):
        variant = 'ENST00000396884.8:c.483ins10'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'The length of the variant is not formatted following the HGVS guidelines. Please rewrite e.g. 10 ' \
               'to N[10]' in \
               results['validation_warning_1']['validation_warnings'][0]
        assert 'An insertion must be provided with the two positions between which the insertion has taken place' in \
               results['validation_warning_1']['validation_warnings'][1]

    def test_issue_359_d(self):
        variant = 'ENST00000396884.8:c.483insA[10]'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'ENST00000396884.8:c.483insA[10] may also be written as ENST00000396884.8:c.483insAAAAAAAAAA' in \
               results['validation_warning_1']['validation_warnings'][0]
        assert 'An insertion must be provided with the two positions between which the insertion has taken place' in \
               results['validation_warning_1']['validation_warnings'][1]
        assert 'Insertion length must be 1 e.g. 483_484insAAAAAAAAAA' in \
               results['validation_warning_1']['validation_warnings'][2]

    def test_issue_359_e(self):
        variant = 'ENST00000396884.8:c.483delinsA[10]'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'ENST00000396884.8:c.483delinsA[10] may also be written as ENST00000396884.8:c.483delinsAAAAAAAAAA' in \
               results['ENST00000396884.8:c.483delinsAAAAAAAAAA']['validation_warnings'][0]

    def test_issue_359_f(self):
        variant = 'ENST00000396884.8:c.483_484insA[10]'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'ENST00000396884.8:c.483_484insA[10] may also be written as ENST00000396884.8:c.483_484insAAAAAAAAAA' in \
               results['ENST00000396884.8:c.484_485insAAAAAAAAAA']['validation_warnings'][0]
    def test_issue_359_g(self):
        variant = 'ENST00000396884.8:c.483_484ins[A[10];T]'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'ENST00000396884.8:c.483_484ins[A[10];T] may also be written as ENST00000396884.8:c.483_484insAAAAAAAAAAT' in \
               results['ENST00000396884.8:c.484_485insAAAAAAAAATA']['validation_warnings'][0]

    def test_issue_359_h(self):
        variant = 'ENST00000396884.8:c.483_484delins[A[10];T]'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'ENST00000396884.8:c.483_484delins[A[10];T] may also be written as ' \
               'ENST00000396884.8:c.483_484delinsAAAAAAAAAAT' in \
               results['ENST00000396884.8:c.483_484delinsAAAAAAAAAAT']['validation_warnings'][0]

    def test_issue_359_i(self):
        variant = 'ENST00000396884.8:c.483ins(10_20)'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'The length of the variant is not formatted following the HGVS guidelines. Please rewrite e.g. (10_20) '\
               'to N[(10_20)](where N is an unknown nucleotide and [(10_20)] is an uncertain number of N nucleotides ' \
               'ranging from 10 to 20)' in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_issue_359_j(self):
        variant = 'ENST00000396884.8:c.483ins[(20_10)]'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'The length of the variant is not formatted following the HGVS guidelines. Please rewrite (20_10) to ' \
               'N[(10_20)]' in \
               results['validation_warning_1']['validation_warnings'][0]
        assert 'An insertion must be provided with the two positions between which the insertion has taken place' in \
               results['validation_warning_1']['validation_warnings'][1]

    def test_issue_359_k(self):
        variant = 'ENST00000396884.8:c.483ins[(20_20)]'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'The length of the variant is not formatted following the HGVS guidelines. Please rewrite ' \
               '(20_20) to N[(20)]' in \
               results['validation_warning_1']['validation_warnings'][0]
        assert 'An insertion must be provided with the two positions between which the insertion has taken place' in \
               results['validation_warning_1']['validation_warnings'][1]

    def test_issue_359_l(self):
        variant = 'ENST00000396884.8:c.483_484ins[(10_20)]'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'The variant description is syntactically correct but no further validation is possible because the ' \
               'description contains uncertainty' in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_issue_359_m(self):
        variant = 'ENST00000396884.8:c.483ins[(10_20)]'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'An insertion must be provided with the two positions between which the insertion has taken place' in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_issue_360(self):
        variant = 'NC_012920.1:g.100del'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert "The given reference sequence (NC_012920.1) does not match the DNA type (g). For NC_012920.1, " \
               "please use (m). For g. variants, please use a linear genomic reference sequence" in \
               results['mitochondrial_variant_1']['validation_warnings'][0]

    def test_issue_360a(self):
        variant = 'NC_012920.1:g.100del'
        results = self.vv.validate(variant, 'hg19', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert "NC_012920.1 is not associated with genome build hg19, instead use genome build GRCh37" in \
               results['mitochondrial_variant_1']['validation_warnings'][0]

    def test_issue_360b(self):
        variant = 'NC_001807.4:g.100del'
        results = self.vv.validate(variant, 'GRCh37', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert "NC_001807.4 is not associated with genome build GRCh37, instead use genome build hg19" in \
               results['mitochondrial_variant_1']['validation_warnings'][0]

    def test_issue_351(self):
        variant = 'M:m.1000_100del'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'This is not a valid HGVS variant description, because no reference sequence ID has been provided, ' \
               'instead use NC_012920.1:m.1000_100del' in \
               results['validation_warning_1']['validation_warnings'][0]
        assert 'The variant positions are valid but we cannot normalize variants spanning the origin of ' \
               'circular reference sequences' in \
               results['validation_warning_1']['validation_warnings'][1]

        variant = 'chr1:g.100000del'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'This is not a valid HGVS variant description, because no reference sequence ID has been provided' in \
               results['intergenic_variant_1']['validation_warnings'][0]

    def test_issue_352(self):
        variant = 'NC_000001.10:o.100_1000del'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'Reference sequence type o. should only be used for circular reference sequences that are ' \
               'not mitochondrial. Instead use m.' in \
               results['validation_warning_1']['validation_warnings'][0]

        variant = 'NC_012920.1:o.100_1000del'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'Reference sequence type o. should only be used for circular reference sequences that are not ' \
               'mitochondrial. Instead use m.' in \
               results['mitochondrial_variant_1']['validation_warnings'][0]

    def test_issue_365(self):
        variant = 'ENST00000553106.6:c.1315+5_1315+6insGTGTAACAG'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results['ENST00000553106.6:c.1315+5_1315+6insGTGTAACAG']['validation_warnings'] == []

    def test_issue_46(self):
        variant = 'ENSP00000269305.4:p.H175_V178delinsX'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results['validation_warning_1']['validation_warnings'] == [
            "The amino acid at position 175 of ENSP00000269305.4 is R not H",
            "The amino acid at position 178 of ENSP00000269305.4 is H not V"
        ]

        variant = 'ENSP00000269305.4:p.H175delinsX'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results['validation_warning_1']['validation_warnings'] == [
            "The amino acid at position 175 of ENSP00000269305.4 is R not H"
        ]

        variant = 'ENSP00000269305.4:p.R175_H178delinsX'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results['validation_warning_1']['validation_warnings'] == [
            "Protein level variant descriptions are not fully supported due to redundancy in the genetic code",
            "ENSP00000269305.4:p.Arg175_His178delinsTer is HGVS compliant and contains a valid reference amino acid description"
        ]
        assert results['validation_warning_1'][
                   'hgvs_predicted_protein_consequence']["tlr"] == "ENSP00000269305.4:p.Arg175_His178delinsXaa"

        variant = 'ENSP00000269305.4:p.R175delinsX'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results['validation_warning_1']['validation_warnings'] == [
            "Protein level variant descriptions are not fully supported due to redundancy in the genetic code",
            "ENSP00000269305.4:p.Arg175delinsTer is HGVS compliant and contains a valid reference amino acid description"
        ]
        assert results['validation_warning_1'][
                   'hgvs_predicted_protein_consequence']["tlr"] == "ENSP00000269305.4:p.Arg175delinsXaa"

    def test_issue_432(self):
        variant = 'ENST00000318312.12:c.1779+7A>G'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)

        print(results)
        assert results['validation_warning_1']['validation_warnings'] == [
            "ENST00000318312.12:c.1779+7A>G auto-mapped to ENST00000318312.12:c.*4A>G",
            "ENST00000318312.12:c.*4A>G: Variant reference (A) does not agree with reference sequence (C)"
        ]

    def test_issue_455(self):
        variant = 'ENSP00000269305.4:p.?'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)

        print(results)
        assert results['validation_warning_1']['validation_warnings'] == [
            "Protein level variant descriptions are not fully supported due to redundancy in the genetic code",
            "ENSP00000269305.4:p.? is HGVS compliant and contains a valid reference amino acid description"
        ]

    def test_issue_518a(self):
        variant = 'ENST00000636147.2(CLN3):c.791-802_1056+1445del'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)

        print(results)
        assert results['ENST00000636147.2:c.790+532_1056+1445del'][
                   'validation_warnings'] == [
            "Removing redundant gene symbol CLN3 from variant description",
            "ExonBoundaryError: Position c.791-802 has been updated to position to 790+532 ensuring correct HGVS "
            "numbering for transcript ENST00000636147.2"
        ]

    def missing_dot(self):
        variant = 'chr11:g,108121787G>A'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert "Unable to identify a dot (.) in the variant description chr11:g,108121787G>A following the reference " \
               "sequence type (g,c,n,r, or p). A dot is required in HGVS variant descriptions to separate the " \
               "reference type from the variant position i.e. <accession>:<type>. e.g. :g." in results[
            'validation_warning_1']["validation_warnings"]

    def missing_colon(self):
        variant = 'chr11g.108121787G>A'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert "Unable to identify a colon (:) in the variant description chr11g.108121787G>A. A colon is required in " \
               "HGVS variant descriptions to separate the reference accession from the reference type i.e. " \
               "<accession>:<type>. e.g. :c." in results[
            'validation_warning_1']["validation_warnings"]

    def p1_a(self):
        variant = 'LRG_199p1:p.(Met1Ala)'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert "Variant NP_003997.1:p.(Met1Ala) affects the initiation amino acid so is better " \
               "described as NP_003997.1:p.(Met1?)" in results[
                'validation_warning_1']["validation_warnings"]

    def p1_b(self):
        variant = 'LRG_199p1:p.Met1Ala'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert "Variant NP_003997.1:p.Met1Ala affects the initiation amino acid so is better " \
               "described as NP_003997.1:p.(Met1?)" in results[
                'validation_warning_1']["validation_warnings"]

    def uppercase_ref_type(self):
        variant = 'DPYD:C.1905+1G>A'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert ("Reference type incorrectly stated in the variant description DPYD:C.1905+1G>A Valid types are "
                "g,c,n,r, or p") in results['validation_warning_1']["validation_warnings"]

    def test_invalid_aa(self):
        variant = 'ENST00000003084.11:p.Z1335P'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert "Invalid amino acid Z stated in description ENST00000003084.11:p.Z1335P" in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_g_with_tc_ref(self):
        variant = 'ENST00000225964.10:g.2559del'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert "Transcript reference sequence input as genomic (g.) reference sequence. " \
               "Did you mean ENST00000225964.10:c.2559del?" in \
               results['validation_warning_1']['validation_warnings']

        variant = 'ENST00000225964.10:g.2559+54_2560del'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert "Transcript reference sequence input as genomic (g.) reference sequence. " \
               "Did you mean ENST00000225964.10:c.2559+54_2560del?" in \
               results['validation_warning_1']['validation_warnings']

    def test_p_with_tc_ref(self):
        variant = 'ENST00000225964.10:p.(Gly197Cys)'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert "Using a nucleotide reference sequence (NM_ NR_ NG_ NC_) to specify protein-level (p.) " \
               "variation is not HGVS compliant. Please select an appropriate protein reference sequence (NP_)" in \
               results['validation_warning_1']['validation_warnings']

    def test_invalid_reference_set_refseq(self):
        variant = 'ENST00000225964.10:g.2559+54_2560del'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="refseq").format_as_dict(test=True)
        print(results)
        assert "InvalidFieldError: The transcript ENST00000225964.10 is not in the RefSeq data set. Please select Ensembl" in \
               results['validation_warning_1']['validation_warnings']

    def test_invalid_reference_set_ensembl_NM(self):
        variant = 'NM_000828.4:c.-3_-2delinsTT'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert "InvalidFieldError: The transcript NM_000828.4 is not in the Ensembl data set. Please select RefSeq" in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_invalid_reference_set_ensembl_NR(self):
        variant = 'NR_015120.4'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert "InvalidFieldError: The transcript NR_015120.4 is not in the Ensembl data set. Please select RefSeq" in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_ensembl_nmd_type(self):
        variant = '1:25808717:C:T'
        results = self.vv.validate(variant, 'GRCh38', 'ENST00000494537.2', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert ("TranscriptTypeError: Cannot identify an in-frame Termination codon in the reference mRNA sequence. "
                "ENST00000494537.2 may not be a valid coding sequence") in \
               results['validation_warning_1']['validation_warnings'][0]


"""
Series of tests that assess uncertain positions using the Ensembl transcript set
"""

#     def test_uncertain_1(self):
#         variant = 'NC_000005.9:g.(90136803_90144453)_(90159675_90261231)dup'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "Uncertain positions are not fully supported, however the syntax is valid" in \
#                results['NM_032119.4:c.(17019+1_17020-1)_(17856+1_17857-1)dup']['validation_warnings']
#         assert "Only a single transcript can be processed, updating to Select" in \
#                results['NM_032119.4:c.(17019+1_17020-1)_(17856+1_17857-1)dup']['validation_warnings']
#         assert results['NM_032119.4:c.(17019+1_17020-1)_(17856+1_17857-1)dup']['primary_assembly_loci'] == {
#             "grch38": {
#                 "hgvs_genomic_description": "NC_000005.9:g.(90136803_90144453)_(90159675_90261231)dup"
#             }}
#         assert results['NM_032119.4:c.(17019+1_17020-1)_(17856+1_17857-1)dup'][
#                    'hgvs_transcript_variant'] == "NM_032119.4:c.(17019+1_17020-1)_(17856+1_17857-1)dup"
#
#     def test_uncertain_2(self):
#         variant = 'NM_006138.4:n.(1_20)_(30_36)del'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "Coding transcript reference sequence input as non-coding transcript (n.) reference sequence. " \
#                "Did you mean NM_006138.4:c.(1_20)_(30_36)del?" in \
#                results['validation_warning_1']['validation_warnings']
#
#     def test_uncertain_3(self):
#         variant = 'NM_006138.4:c.(1_20)_(30_36)del'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "Uncertain positions are not fully supported, however the syntax is valid" in \
#                results['NM_006138.4:c.(1_20)_(30_36)del']['validation_warnings']
#         assert results['NM_006138.4:c.(1_20)_(30_36)del']['hgvs_transcript_variant'] == "NM_006138.4:c.(1_20)_(30_36)del"
#         assert results['NM_006138.4:c.(1_20)_(30_36)del']['primary_assembly_loci'] == {
#             "grch38": {
#                 "hgvs_genomic_description": "NC_000011.10:g.(60061161_60061180)_(60061190_60061196)del"
#             }}
#
#     def test_uncertain_4(self):
#         variant = 'NM_032119.3:c.(17019+1_17020-1)_(17856+1_17857-1)dup'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "Uncertain positions are not fully supported, however the syntax is valid" in \
#                results['NM_032119.3:c.(17019+1_17020-1)_(17856+1_17857-1)dup']['validation_warnings']
#         assert results['NM_032119.3:c.(17019+1_17020-1)_(17856+1_17857-1)dup'][
#                    'hgvs_transcript_variant'] == "NM_032119.3:c.(17019+1_17020-1)_(17856+1_17857-1)dup"
#         assert results['NM_032119.3:c.(17019+1_17020-1)_(17856+1_17857-1)dup']['primary_assembly_loci'] == {
#             "grch38": {
#                 "hgvs_genomic_description": "NC_000005.10:g.(90840986_90848636)_(90863858_90965414)dup"
#             }}
#
#     def test_uncertain_5(self):
#         variant = 'NC_000005.9:g.(90159675_90261231)_(90136803_90144453)dup'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "Position 90159675_90261231 is > or overlaps 90136803_90144453" in results[
#             'validation_warning_1']["validation_warnings"]
#
#     def test_uncertain_6(self):
#         variant = 'NC_000005.9:g.(90144453_90136803)_(90159675_90261231)dup'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "base start position must be <= end position in position 90144453_90136803" in results[
#             'validation_warning_1']["validation_warnings"]
#
#     def test_uncertain_7(self):
#         variant = 'NC_000003.12:g.(63912602_63912844)insN[15]'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NC_000003.12:g.(63912602_63912844)insN[15] may also be written as " \
#                "NC_000003.12:g.(63912602_63912844)insNNNNNNNNNNNNNNN" in results[
#             'NM_001377405.1:c.(4_246)insNNNNNNNNNNNNNNN']["validation_warnings"]
#         assert results['NM_001377405.1:c.(4_246)insNNNNNNNNNNNNNNN'][
#                    'hgvs_transcript_variant'] == "NM_001377405.1:c.(4_246)insNNNNNNNNNNNNNNN"
#         assert results['NM_001377405.1:c.(4_246)insNNNNNNNNNNNNNNN']['primary_assembly_loci'] == {
#             "grch38": {
#                 "hgvs_genomic_description": "NC_000003.12:g.(63912602_63912844)insNNNNNNNNNNNNNNN"
#             }}
#
#     def test_uncertain_8(self):
#         variant = 'NC_000003.12:g.(63912602_63912844)delN[15]'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NC_000003.12:g.(63912602_63912844)delN[15] may also be written as " \
#                "NC_000003.12:g.(63912602_63912844)delNNNNNNNNNNNNNNN" in results[
#             'NM_001377405.1:c.(4_246)delNNNNNNNNNNNNNNN']["validation_warnings"]
#         assert results['NM_001377405.1:c.(4_246)delNNNNNNNNNNNNNNN'][
#                    'hgvs_transcript_variant'] == "NM_001377405.1:c.(4_246)delNNNNNNNNNNNNNNN"
#         assert results['NM_001377405.1:c.(4_246)delNNNNNNNNNNNNNNN']['primary_assembly_loci'] == {
#             "grch38": {
#                 "hgvs_genomic_description": "NC_000003.12:g.(63912602_63912844)delNNNNNNNNNNNNNNN"
#             }}
#
#     def test_uncertain_9(self):
#         variant = 'NM_001377405.1:c.(4_246)delN[15]'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NM_001377405.1:c.(4_246)delN[15] may also be written as " \
#                "NM_001377405.1:c.(4_246)delNNNNNNNNNNNNNNN" in results[
#             'NM_001377405.1:c.(4_246)delNNNNNNNNNNNNNNN']["validation_warnings"]
#         assert results['NM_001377405.1:c.(4_246)delNNNNNNNNNNNNNNN'][
#                    'hgvs_transcript_variant'] == "NM_001377405.1:c.(4_246)delNNNNNNNNNNNNNNN"
#         assert results['NM_001377405.1:c.(4_246)delNNNNNNNNNNNNNNN']['primary_assembly_loci'] == {
#             "grch38": {
#                 "hgvs_genomic_description": "NC_000003.12:g.(63912602_63912844)delNNNNNNNNNNNNNNN"
#             }}
#
#     def test_uncertain_10(self):
#         variant = 'NM_001377405.1:c.(4_246)insN[15]'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NM_001377405.1:c.(4_246)insN[15] may also be written as " \
#                "NM_001377405.1:c.(4_246)insNNNNNNNNNNNNNNN" in results[
#             'NM_001377405.1:c.(4_246)insNNNNNNNNNNNNNNN']["validation_warnings"]
#         assert results['NM_001377405.1:c.(4_246)insNNNNNNNNNNNNNNN'][
#                    'hgvs_transcript_variant'] == "NM_001377405.1:c.(4_246)insNNNNNNNNNNNNNNN"
#         assert results['NM_001377405.1:c.(4_246)insNNNNNNNNNNNNNNN']['primary_assembly_loci'] == {
#             "grch38": {
#                 "hgvs_genomic_description": "NC_000003.12:g.(63912602_63912844)insNNNNNNNNNNNNNNN"
#             }}

"""
Series of tests fot the Allele syntax using the Ensembl transcript set
"""

#     def test_alleles_1(self):
#         variant = 'NM_000093.5:c.[14del;17G>A]'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "AlleleSyntaxError: Variants [14del;17G>A] should be merged into NM_000093.5:c.16_17delinsA" in results[
#             'validation_warning_1']["validation_warnings"]
#
#     def test_alleles_2(self):
#         variant = 'NM_000088.4:c.[4del;6C>G]'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "AlleleSyntaxError: Variants [4del;6C>G] should be merged into NM_000088.4:c.5_6delinsG" in results[
#             'validation_warning_1']["validation_warnings"]
#
#     def test_alleles_3(self):
#         variant = 'NM_000088.4:c.[589-1del;591T>A]'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "AlleleSyntaxError: Intronic variants can only be validated if a genomic/gene reference sequence" \
#                " is also provided " \
#                "e.g. NC_000017.11(NM_000088.3):c.589-1G>T" in results[
#             'validation_warning_1']["validation_warnings"]
#
#     def test_alleles_4(self):
#         variant = 'NC_000017.11(NM_000088.4):c.[589-1del;591T>A]'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "AlleleSyntaxError: Variants [589-1del;591T>A] should be merged into " \
#                "NM_000088.4:c.590_591delinsA" in results[
#             'validation_warning_1']["validation_warnings"]
#
#     def test_alleles_5(self):
#         variant = 'NC_000009.12(NM_000093.5):c.[277del;277+2T>A]'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "AlleleSyntaxError: Variants [277del;277+2T>A] should be merged into " \
#                "NM_000093.5:c.277+1_277+2delinsA" in results[
#             'validation_warning_1']["validation_warnings"]


"""
Series of gap warning tests using the Ensembl transcript set using VariantFormatter
"""
# class TestVFGapWarnings(TestCase):
#
#     def test_vf_series_1(self):
#         results = simpleVariantFormatter.format('NC_000004.11:g.140811117C>A', 'GRCh37', 'refseq', None, False, True,
#                                                 testing=True)
#         print(results)
#         assert 'NC_000004.11:g.140811117C>A' in results.keys()
#         assert 'NM_018717.4 contains 3 fewer bases between c.2276_2277, and 12 fewer bases between c.1467_1468 ' \
#                'than NC_000004.11' in results['NC_000004.11:g.140811117C>A'][
#             'NC_000004.11:g.140811117C>A']['hgvs_t_and_p']['NM_018717.4']['gap_statement']
#
#     def test_vf_series_2(self):
#         results = simpleVariantFormatter.format('NC_000008.10:g.24811072C>T',
#                                                                  'GRCh37', 'refseq', None, False, True, testing=True)
#         print(results)
#         assert 'NC_000008.10:g.24811072C>T' in results.keys()
#         assert 'NM_006158.3 contains 1 fewer bases between c.1407_1408 than NC_000008.10' in results[
#             'NC_000008.10:g.24811072C>T']['NC_000008.10:g.24811072C>T']['hgvs_t_and_p'][
#             'NM_006158.3']['gap_statement']
#         assert 'NM_006158.4 contains 1 fewer bases between c.1407_1408 than NC_000008.10' in results[
#             'NC_000008.10:g.24811072C>T']['NC_000008.10:g.24811072C>T']['hgvs_t_and_p'][
#             'NM_006158.4']['gap_statement']
#         assert 'NM_006158.5 contains 1 fewer bases between c.1413_1414 than NC_000008.10' in results[
#             'NC_000008.10:g.24811072C>T']['NC_000008.10:g.24811072C>T']['hgvs_t_and_p'][
#             'NM_006158.5']['gap_statement']
#
#     def test_vf_series_3(self):
#         results = simpleVariantFormatter.format('NC_000015.9:g.72105933del',
#                                                                  'GRCh37', 'refseq', None, False, True, testing=True)
#         print(results)
#         assert 'NC_000015.9:g.72105933del' in results.keys()
#         assert 'NM_014249.2 contains 1 fewer bases between c.947_948 than NC_000015.9' in results[
#             'NC_000015.9:g.72105933del']['NC_000015.9:g.72105933del']['hgvs_t_and_p'][
#             'NM_014249.2']['gap_statement']
#         assert 'NM_014249.3 contains 1 fewer bases between c.947_948 than NC_000015.9' in results[
#             'NC_000015.9:g.72105933del']['NC_000015.9:g.72105933del']['hgvs_t_and_p'][
#             'NM_014249.3']['gap_statement']
#         assert 'NM_014249.4 contains 1 fewer bases between c.951_952 than NC_000015.9' in results[
#             'NC_000015.9:g.72105933del']['NC_000015.9:g.72105933del']['hgvs_t_and_p'][
#             'NM_014249.4']['gap_statement']
#         assert 'NM_016346.2 contains 1 fewer bases between c.947_948 than NC_000015.9' in results[
#             'NC_000015.9:g.72105933del']['NC_000015.9:g.72105933del']['hgvs_t_and_p'][
#             'NM_016346.2']['gap_statement']
#         assert 'NM_016346.3 contains 1 fewer bases between c.947_948 than NC_000015.9' in results[
#             'NC_000015.9:g.72105933del']['NC_000015.9:g.72105933del']['hgvs_t_and_p'][
#             'NM_016346.3']['gap_statement']
#         assert 'NM_016346.4 contains 1 fewer bases between c.951_952 than NC_000015.9' in results[
#             'NC_000015.9:g.72105933del']['NC_000015.9:g.72105933del']['hgvs_t_and_p'][
#             'NM_016346.4']['gap_statement']
#
#     def test_vf_series_4(self):
#         results = simpleVariantFormatter.format('NC_000019.9:g.41123095dup',
#                                                                  'GRCh37', 'refseq', None, False, True, testing=True)
#         print(results)
#         assert 'NC_000019.9:g.41123095dup' in results.keys()
#         assert 'NM_001042544.1 contains 1 extra bases between c.3233_3235 than NC_000019.9' in results[
#             'NC_000019.9:g.41123095dup']['NC_000019.9:g.41123095dup']['hgvs_t_and_p'][
#             'NM_001042544.1']['gap_statement']
#         assert 'NM_001042545.1 contains 1 extra bases between c.3032_3034 than NC_000019.9' in results[
#             'NC_000019.9:g.41123095dup']['NC_000019.9:g.41123095dup']['hgvs_t_and_p'][
#             'NM_001042545.1']['gap_statement']
#         assert 'NM_001042545.2 contains 1 extra bases between c.3034_3036 than NC_000019.9' in results[
#             'NC_000019.9:g.41123095dup']['NC_000019.9:g.41123095dup']['hgvs_t_and_p'][
#             'NM_001042545.2']['gap_statement']
#         assert 'NM_003573.2 contains 1 extra bases between c.3122_3124 than NC_000019.9' in results[
#             'NC_000019.9:g.41123095dup']['NC_000019.9:g.41123095dup']['hgvs_t_and_p'][
#             'NM_003573.2']['gap_statement']
#
#     def test_vf_series_5(self):
#         results = simpleVariantFormatter.format('NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=',
#                                                                  'GRCh37', 'refseq', None, False, True, testing=True)
#         print(results)
#         assert 'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=' in results.keys()
#         assert 'NM_001083585.1 contains 25 fewer bases between c.*344_*345 than NC_000017.10' in results[
#             'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
#             'NM_001083585.1']['gap_statement']
#         assert 'NM_001083585.2 contains 25 fewer bases between c.*344_*345 than NC_000017.10' in results[
#             'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
#             'NM_001083585.2']['gap_statement']
#         assert 'NM_001083585.3 contains 25 fewer bases between c.*369_*370 than NC_000017.10' in results[
#             'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
#             'NM_001083585.3']['gap_statement']
#         assert 'NM_001291581.1 contains 25 fewer bases between c.*344_*345 than NC_000017.10' in results[
#             'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
#             'NM_001291581.1']['gap_statement']
#         assert 'NM_001291581.2 contains 25 fewer bases between c.*369_*370 than NC_000017.10' in results[
#             'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
#             'NM_001291581.2']['gap_statement']
#         assert 'NM_004703.4 contains 25 fewer bases between c.*344_*345 than NC_000017.10' in results[
#             'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
#             'NM_004703.4']['gap_statement']
#         assert 'NM_004703.5 contains 25 fewer bases between c.*344_*345 than NC_000017.10' in results[
#             'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
#             'NM_004703.5']['gap_statement']
#         assert 'NM_004703.6 contains 25 fewer bases between c.*369_*370 than NC_000017.10' in results[
#             'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=']['hgvs_t_and_p'][
#             'NM_004703.6']['gap_statement']
#
#     def test_vf_series_6(self):
#         results = simpleVariantFormatter.format('NC_000012.11:g.122064777C>A',
#                                                                  'GRCh37', 'refseq', None, False, True, testing=True)
#         print(results)
#         assert 'NC_000012.11:g.122064777C>A' in results.keys()
#         assert 'NM_032790.3 contains 6 fewer bases between c.126_127 than NC_000012.11' in results[
#             'NC_000012.11:g.122064777C>A']['NC_000012.11:g.122064777C>A']['hgvs_t_and_p'][
#             'NM_032790.3']['gap_statement']
#
#     def test_vf_series_7(self):
#         results = simpleVariantFormatter.format('NC_000002.11:g.95847041_95847043GCG=',
#                                                                  'GRCh37', 'refseq', None, False, True, testing=True)
#         print(results)
#         assert 'NC_000002.11:g.95847041_95847043GCG=' in results.keys()
#         assert 'NM_001017396.1 contains 3 fewer bases between c.341_342 than NC_000002.11' in results[
#             'NC_000002.11:g.95847041_95847043GCG=']['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'][
#             'NM_001017396.1']['gap_statement']
#         assert 'NM_001017396.2 contains 3 fewer bases between c.341_342 than NC_000002.11' in results[
#             'NC_000002.11:g.95847041_95847043GCG=']['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'][
#             'NM_001017396.2']['gap_statement']
#         assert 'NM_001282398.1 contains 3 fewer bases between c.353_354 than NC_000002.11' in results[
#             'NC_000002.11:g.95847041_95847043GCG=']['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'][
#             'NM_001282398.1']['gap_statement']
#         assert 'NM_001291604.1 contains 3 fewer bases between c.227_228 than NC_000002.11' in results[
#             'NC_000002.11:g.95847041_95847043GCG=']['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'][
#             'NM_001291604.1']['gap_statement']
#         assert 'NM_001291605.1 contains 3 fewer bases between c.506_507 than NC_000002.11' in results[
#             'NC_000002.11:g.95847041_95847043GCG=']['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'][
#             'NM_001291605.1']['gap_statement']
#         assert 'NM_021088.2 contains 3 fewer bases between c.467_468 than NC_000002.11' in results[
#             'NC_000002.11:g.95847041_95847043GCG=']['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'][
#             'NM_021088.2']['gap_statement']
#         assert 'NM_021088.3 contains 3 fewer bases between c.467_468 than NC_000002.11' in results[
#             'NC_000002.11:g.95847041_95847043GCG=']['NC_000002.11:g.95847041_95847043GCG=']['hgvs_t_and_p'][
#             'NM_021088.3']['gap_statement']
#
#     def test_vf_series_8(self):
#         results = simpleVariantFormatter.format('NC_000003.11:g.14561629_14561630insG',
#                                                 'GRCh37', 'refseq', None, False, True, testing=True)
#         print(results)
#         assert 'NM_001080423.2 contains 1 extra bases between c.1308_1310 than NC_000003.11' in results[
#             'NC_000003.11:g.14561629_14561630insG']['NC_000003.11:g.14561629_14561630insG']['hgvs_t_and_p'][
#             'NM_001080423.2']['gap_statement']
#         assert 'NM_001080423.3 contains 1 extra bases between c.1017_1019 than NC_000003.11' in results[
#             'NC_000003.11:g.14561629_14561630insG']['NC_000003.11:g.14561629_14561630insG']['hgvs_t_and_p'][
#             'NM_001080423.3']['gap_statement']
#         assert 'NM_001080423.4 contains 1 extra bases between c.1019_1021 than NC_000003.11' in results[
#             'NC_000003.11:g.14561629_14561630insG']['NC_000003.11:g.14561629_14561630insG']['hgvs_t_and_p'][
#             'NM_001080423.4']['gap_statement']
#
#     def test_vf_series_9(self):
#         results = simpleVariantFormatter.format('NC_000004.11:g.140811117C>A',
#                                                 'GRCh37', 'refseq', None, False, True, testing=True)
#         print(results)
#         assert 'NC_000004.11:g.140811117C>A' in results.keys()
#         assert 'NM_018717.4 contains 3 fewer bases between c.2276_2277, and 12 fewer bases between c.1467_1468 ' \
#                'than NC_000004.11' in results[
#             'NC_000004.11:g.140811117C>A']['NC_000004.11:g.140811117C>A']['hgvs_t_and_p'][
#             'NM_018717.4']['gap_statement']
#
#     def test_vf_series_10(self):
#         results = simpleVariantFormatter.format('NC_000009.11:g.136132908_136132909TA=',
#                                                 'GRCh37', 'refseq', None, False, True, testing=True)
#         print(results)
#         assert 'NC_000009.11:g.136132908_136132909TA=' in results.keys()
#         assert 'NM_020469.2 contains 1 extra bases between c.260_262 than NC_000009.11' in results[
#             'NC_000009.11:g.136132908_136132909TA=']['NC_000009.11:g.136132908_136132909TA=']['hgvs_t_and_p'][
#             'NM_020469.2']['gap_statement']
#         assert 'NM_020469.3 contains 22 extra bases between c.*756_*757, and 2 extra bases between c.*797_*798, ' \
#                'and 110 extra bases between c.*840_*841, and 2 extra bases between c.*4648_*4649, and 1 extra ' \
#                'bases between c.260_262 than NC_000009.11' in results[
#             'NC_000009.11:g.136132908_136132909TA=']['NC_000009.11:g.136132908_136132909TA=']['hgvs_t_and_p'][
#             'NM_020469.3']['gap_statement']
#
#     def test_vf_series_11(self):
#         results = simpleVariantFormatter.format('NC_000019.10:g.50378563_50378564insTAC',
#                                                 'GRCh38', 'refseq', None, False, True, testing=True)
#         print(results)
#         assert 'NC_000019.10:g.50378563_50378564insTAC' in results.keys()
#         assert 'NM_001256647.1 contains 3 extra bases between c.223_227 than NC_000019.10' in results[
#             'NC_000019.10:g.50378563_50378564insTAC']['NC_000019.10:g.50378563_50378564insTAC']['hgvs_t_and_p'][
#             'NM_001256647.1']['gap_statement']
#         assert 'NM_007121.5 contains 3 extra bases between c.514_518 than NC_000019.10' in results[
#             'NC_000019.10:g.50378563_50378564insTAC']['NC_000019.10:g.50378563_50378564insTAC']['hgvs_t_and_p'][
#             'NM_007121.5']['gap_statement']
#
#     def test_vf_series_12(self):
#         results = simpleVariantFormatter.format('NC_000007.13:g.149476664_149476666delinsTC',
#                                                 'GRCh37', 'refseq', None, False, True, testing=True)
#         print(results)
#         assert 'NC_000007.13:g.149476664_149476666delinsTC' in results.keys()
#         assert 'NR_163594.1 contains 1 extra bases between n.1129_1131, and 1 fewer bases between n.11675_11676 ' \
#                'than NC_000007.13' in results[
#             'NC_000007.13:g.149476664_149476666delinsTC']['NC_000007.13:g.149476664_149476666delinsTC'][
#             'hgvs_t_and_p']['NR_163594.1']['gap_statement']
#
#     def test_vf_series_13(self):
#         results = simpleVariantFormatter.format('NC_000004.12:g.139889957_139889968del',
#                                                 'GRCh38', 'refseq', None, False, True, testing=True)
#         print(results)
#         assert 'NC_000004.12:g.139889957_139889968del' in results.keys()
#         assert 'NM_018717.4 contains 3 fewer bases between c.2276_2277, and 12 fewer bases between c.1467_1468 ' \
#                'than NC_000004.12' in results[
#             'NC_000004.12:g.139889957_139889968del']['NC_000004.12:g.139889957_139889968del']['hgvs_t_and_p'][
#             'NM_018717.4']['gap_statement']
#

"""
Series of gap warning tests using the Ensembl transcript set using VariantValidator
"""

# class TestVVGapWarnings(TestCase):
#
#     @classmethod
#     def setup_class(cls):
#         cls.vv = Validator()
#         cls.vv.testing = True
#
#     def test_vv_series_1(self):
#         variant = 'NC_000004.11:g.140811117C>A'
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NM_018717.4 contains 3 fewer bases between c.2276_2277, and 12 fewer bases between c.1467_1468 than NC_000004.11" in \
#                results['NM_018717.4:c.1472_1473insTCAGCAGCAGCA']['validation_warnings']
#
#     def test_vv_series_2(self):
#         variant = 'NC_000008.10:g.24811072C>T'
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NM_006158.5 contains 1 fewer bases between c.1413_1414 than NC_000008.10" in \
#                results['NM_006158.5:c.1407delinsAC']['validation_warnings']
#         assert "NM_006158.4 contains 1 fewer bases between c.1407_1408 than NC_000008.10" in \
#                results['NM_006158.4:c.1407delinsAC']['validation_warnings']
#         assert "NM_006158.3 contains 1 fewer bases between c.1407_1408 than NC_000008.10" in \
#                results['NM_006158.3:c.1407delinsAC']['validation_warnings']
#
#     def test_vv_series_3(self):
#         variant = 'NC_000015.9:g.72105933del'
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NM_016346.4 contains 1 fewer bases between c.951_952 than NC_000015.9" in \
#                results['NM_016346.4:c.951_952=']['validation_warnings']
#         assert "NM_016346.3 contains 1 fewer bases between c.947_948 than NC_000015.9" in \
#                results['NM_016346.3:c.947_948=']['validation_warnings']
#         assert "NM_016346.2 contains 1 fewer bases between c.947_948 than NC_000015.9" in \
#                results['NM_016346.2:c.947_948=']['validation_warnings']
#         assert "NM_014249.4 contains 1 fewer bases between c.951_952 than NC_000015.9" in \
#                results['NM_014249.4:c.951_952=']['validation_warnings']
#         assert "NM_014249.3 contains 1 fewer bases between c.947_948 than NC_000015.9" in \
#                results['NM_014249.3:c.947_948=']['validation_warnings']
#         assert "NM_014249.2 contains 1 fewer bases between c.947_948 than NC_000015.9" in \
#                results['NM_014249.2:c.947_948=']['validation_warnings']
#
#     def test_vv_series_4(self):
#         variant = 'NC_000019.9:g.41123095dup'
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NM_003573.2 contains 1 extra bases between c.3122_3124 than NC_000019.9" in \
#                results['NM_003573.2:c.3122_3124=']['validation_warnings']
#         assert "NM_001042545.2 contains 1 extra bases between c.3034_3036 than NC_000019.9" in \
#                results['NM_001042545.2:c.3033_3036=']['validation_warnings']
#         assert "NM_001042545.1 contains 1 extra bases between c.3032_3034 than NC_000019.9" in \
#                results['NM_001042545.1:c.3032_3034=']['validation_warnings']
#         assert "NM_001042544.1 contains 1 extra bases between c.3233_3235 than NC_000019.9" in \
#                results['NM_001042544.1:c.3233_3235=']['validation_warnings']
#
#     def test_vv_series_5(self):
#         variant = 'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG='
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NM_004703.6 contains 25 fewer bases between c.*369_*370 than NC_000017.10" in \
#                results['NM_004703.6:c.*344_*368dup']['validation_warnings']
#         assert "NM_004703.5 contains 25 fewer bases between c.*344_*345 than NC_000017.10" in \
#                results['NM_004703.5:c.*344_*368dup']['validation_warnings']
#         assert "NM_004703.4 contains 25 fewer bases between c.*344_*345 than NC_000017.10" in \
#                results['NM_004703.4:c.*344_*368dup']['validation_warnings']
#         assert "NM_001291581.2 contains 25 fewer bases between c.*369_*370 than NC_000017.10" in \
#                results['NM_001291581.2:c.*344_*368dup']['validation_warnings']
#         assert "NM_001291581.1 contains 25 fewer bases between c.*344_*345 than NC_000017.10" in \
#                results['NM_001291581.1:c.*344_*368dup']['validation_warnings']
#         assert "NM_001083585.3 contains 25 fewer bases between c.*369_*370 than NC_000017.10" in \
#                results['NM_001083585.3:c.*344_*368dup']['validation_warnings']
#         assert "NM_001083585.2 contains 25 fewer bases between c.*344_*345 than NC_000017.10" in \
#                results['NM_001083585.2:c.*344_*368dup']['validation_warnings']
#         assert "NM_001083585.1 contains 25 fewer bases between c.*344_*345 than NC_000017.10" in \
#                results['NM_001083585.1:c.*344_*368dup']['validation_warnings']
#
#     def test_vv_series_6(self):
#         variant = 'NC_000012.11:g.122064777C>A'
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NM_032790.3 contains 6 fewer bases between c.126_127 than NC_000012.11" in \
#                results['NM_032790.3:c.129_130insACACCG']['validation_warnings']
#
#     def test_vv_series_7(self):
#         variant = 'NC_000002.11:g.95847041_95847043GCG='
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NM_021088.3 contains 3 fewer bases between c.467_468 than NC_000002.11" in \
#                results['NM_021088.3:c.471_473dup']['validation_warnings']
#         assert "NM_021088.2 contains 3 fewer bases between c.467_468 than NC_000002.11" in \
#                results['NM_021088.2:c.471_473dup']['validation_warnings']
#         assert "NM_001291605.1 contains 3 fewer bases between c.506_507 than NC_000002.11" in \
#                results['NM_001291605.1:c.510_512dup']['validation_warnings']
#         assert "NM_001291604.1 contains 3 fewer bases between c.227_228 than NC_000002.11" in \
#                results['NM_001291604.1:c.231_233dup']['validation_warnings']
#         assert "NM_001282398.1 contains 3 fewer bases between c.353_354 than NC_000002.11" in \
#                results['NM_001282398.1:c.357_359dup']['validation_warnings']
#         assert "NM_001017396.2 contains 3 fewer bases between c.341_342 than NC_000002.11" in \
#                results['NM_001017396.2:c.345_347dup']['validation_warnings']
#         assert "NM_001017396.1 contains 3 fewer bases between c.341_342 than NC_000002.11" in \
#                results['NM_001017396.1:c.345_347dup']['validation_warnings']
#
#     def test_vv_series_8(self):
#         variant = 'NC_000003.11:g.14561629_14561630insG'
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NM_001080423.4 contains 1 extra bases between c.1019_1021 than NC_000003.11" in \
#                results['NM_001080423.4:c.1019_1021=']['validation_warnings']
#         assert "NM_001080423.3 contains 1 extra bases between c.1017_1019 than NC_000003.11" in \
#                results['NM_001080423.3:c.1017_1020=']['validation_warnings']
#         assert "NM_001080423.2 contains 1 extra bases between c.1308_1310 than NC_000003.11" in \
#                results['NM_001080423.2:c.1308_1311=']['validation_warnings']
#
#     def test_vv_series_9(self):
#         variant = 'NC_000004.11:g.140811117C>A'
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NM_018717.4 contains 3 fewer bases between c.2276_2277, and 12 fewer bases between c.1467_1468 than NC_000004.11" in \
#                results['NM_018717.4:c.1472_1473insTCAGCAGCAGCA']['validation_warnings']
#
#     def test_vv_series_10(self):
#         variant = 'NC_000009.11:g.136132908_136132909TA='
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NM_020469.3 contains 22 extra bases between c.*756_*757, and 2 extra bases between c.*797_*798, and 110 extra bases between c.*840_*841, and 2 extra bases between c.*4648_*4649, and 1 extra bases between c.260_262 than NC_000009.11" in \
#                results['NM_020469.3:c.261del']['validation_warnings']
#         assert "NM_020469.2 contains 1 extra bases between c.260_262 than NC_000009.11" in \
#                results['NM_020469.2:c.261del']['validation_warnings']
#
#     def test_vv_series_11(self):
#         variant = 'NC_000019.10:g.50378563_50378564insTAC'
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NM_007121.5 contains 3 extra bases between c.514_518 than NC_000019.10" in \
#                results['NM_007121.5:c.515A>T']['validation_warnings']
#         assert "NM_001256647.1 contains 3 extra bases between c.223_227 than NC_000019.10" in \
#                results['NM_001256647.1:c.224A>T']['validation_warnings']
#
#     def test_vv_series_12(self):
#         variant = 'NC_000007.13:g.149476664_149476666delinsTC'
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NR_163594.1 contains 1 extra bases between n.1129_1131, and 1 fewer bases between n.11675_11676 than NC_000007.13" in \
#                results['NR_163594.1:n.1122_1124delinsT']['validation_warnings']
#
#     def test_vv_series_13(self):
#         variant = 'NC_000004.12:g.139889957_139889968del'
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "NM_018717.4 contains 3 fewer bases between c.2276_2277, and 12 fewer bases between c.1467_1468 than NC_000004.12" in \
#                results['NM_018717.4:c.1466_1468=']['validation_warnings']
#
#     def test_vv_series_14(self):
#         variant = 'NM_000516.7:c.2780+73C>T'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "CDSError: Variant start position and/or end position are beyond the CDS end position and likely also beyond the end of the selected reference sequence" in \
#                results['validation_warning_1']['validation_warnings']
#
#     def test_vv_series_15(self):
#         variant = 'NM_000518.5:c.89+25del'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "ExonBoundaryError: Position c.89+25 does not correspond with an exon boundary for transcript NM_000518.5" in \
#                results['validation_warning_1']['validation_warnings']
#
#     def test_vv_series_16(self):
#         variant = 'NM_207122.2:c.1174_1174+1insAT'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "ExonBoundaryError: Position c.1174+1 does not correspond with an exon boundary for transcript NM_207122.2" in \
#                results['validation_warning_1']['validation_warnings']
#
#     def test_vv_series_17(self):
#         variant = 'chr17:g.7578554_7578555delinsCC'
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "This is not a valid HGVS variant description, because no reference sequence ID has been provided" in \
#                results['NM_001276761.3:c.259T>G']['validation_warnings']
#
#     def test_vv_series_18(self):
#         variant = 'NR_033955.2:r.164c>a'
#         results = self.vv.validate(variant, 'GRCh37', 'all').format_as_dict(test=True)
#         print(results)
#         assert "Invalid variant type for non-coding transcript. Instead use n." in \
#                results['validation_warning_1']['validation_warnings']
#
#     def test_vv_series_17(self):
#         variant = 'NM_000086.2(CLN3):c.791-802_1056+1445del'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "ExonBoundaryError: Position c.791-802 has been updated to position to 790+532 ensuring correct HGVS " \
#                "numbering for transcript NM_000086.2" in \
#                results['NM_000086.2:c.790+532_1056+1445del']['validation_warnings']
#
#     def test_vv_series_17a(self):
#         variant = 'NM_000088.4:c.2559_2559+54del'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "ExonBoundaryError: Position c.2559+54 has been updated to position to 2560-35 ensuring correct HGVS " \
#                "numbering for transcript NM_000088.4" in \
#                results['NM_000088.4:c.2559_2560-35del']['validation_warnings']
#
#     def test_vv_series_17b(self):
#         variant = 'NM_000086.2:c.790_791-802del'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "ExonBoundaryError: Position c.791-802 has been updated to position to 790+532 ensuring correct HGVS " \
#                "numbering for transcript NM_000086.2" in \
#                results['NM_000086.2:c.790+1_790+533del']['validation_warnings']
#
#     def test_vv_series_17c(self):
#         variant = 'NM_000088.4:c.2559+54_2560del'
#         results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
#         print(results)
#         assert "ExonBoundaryError: Position c.2559+54 has been updated to position to 2560-35 ensuring correct HGVS " \
#                "numbering for transcript NM_000088.4" in \
#                results['NM_000088.4:c.2560-34_2561del']['validation_warnings']
#

#


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
