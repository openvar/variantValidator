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
        assert 'RnaAlphabetError: The IUPAC RNA alphabet dictates that RNA variants must use the character u in place of t' in \
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
        assert 'ReferenceMismatchError: ENST00000357033.9:c.8810A>G: Variant reference (A) does not agree with reference sequence (G)' in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_issue_180a(self):
        variant = 'NC_000017.10:g.41232400_41236235del383'
        results = self.vv.validate(variant, 'hg19', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'InvalidRangeError: Length implied by coordinates must equal sequence deletion length' in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_issue_180b(self):
        variant = 'NC_000017.11(ENST00000357654.9):c.4186-1642_4358-983del10'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'InvalidRangeError: Length implied by coordinates must equal sequence deletion length' in \
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
        assert 'InvalidRangeError: Length implied by coordinates must equal sequence deletion length' in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_issue_195a(self):
        variant = 'ENST00000225964.10(COL1A1):c.590delG'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'VariantSyntaxError: Removing redundant gene symbol COL1A1 from variant description' in \
               results['ENST00000225964.10:c.590del']['validation_warnings'][0]
        assert 'VariantSyntaxError: Removing redundant reference bases from variant description' in \
               results['ENST00000225964.10:c.590del']['validation_warnings'][1]

    def test_issue_216a(self):
        variant = 'ENST00000396884.8:c.850_877dup27'
        results = self.vv.validate(variant, 'hg38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'InvalidRangeError: Length implied by coordinates must equal sequence duplication length' in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_issue_216b(self):
        variant = 'ENST00000396884.8:c.850_877dup28'
        results = self.vv.validate(variant, 'hg38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'VariantSyntaxError: Trailing digits are not permitted in HGVS variant descriptions' in \
               results['ENST00000396884.8:c.850_877dup']['validation_warnings'][0]

    def test_issue_239(self):
        variant = 'ENST00000396884.8:c.1047dupT'
        results = self.vv.validate(variant, 'hg38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'VariantSyntaxError: Removing redundant reference bases from variant description' in \
               results['ENST00000396884.8:c.1047dup']['validation_warnings'][0]

    def test_issue_338(self):
        # Also issue 357
        variant = 'ENST00000225964.10:C.589G>T'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert ('ReferenceTypeError: Reference type incorrectly stated in the variant description ENST00000225964.10:C.589G>T Valid '
                'types are g,c,n,r, or p') in \
               results['ENST00000225964.10:c.589G>T']['validation_warnings'][0]

        variant = 'enst00000225964.10:c.589G>T'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'InvalidCaseError: This is not a valid HGVS description because characters are in the wrong case. Please check the use of upper- and lowercase characters.' in \
               results['ENST00000225964.10:c.589G>T']['validation_warnings'][0]


    def test_issue_359(self):
        variant = 'ENST00000396884.8:c.483ins'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'InsertionSequenceError: The inserted sequence must be provided for insertions or deletion-insertions' in \
               results['validation_warning_1']['validation_warnings'][0]
        assert 'InsertionLengthError: An insertion must be provided with the two positions between which the insertion has taken place' in \
               results['validation_warning_1']['validation_warnings'][1]

    def test_issue_359_b(self):
        variant = 'ENST00000396884.8:c.483ins(10)'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert 'InsertionLengthError: The length of the variant is not formatted following the HGVS guidelines. Please rewrite e.g. (10) ' \
               'to N[10]' in \
               results['validation_warning_1']['validation_warnings'][0]
        assert 'InsertionLengthError: An insertion must be provided with the two positions between which the insertion has taken place' in \
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
            "AminoMismatchError: The amino acid at position 175 of ENSP00000269305.4 is R not H",
            "AminoMismatchError: The amino acid at position 178 of ENSP00000269305.4 is H not V"
        ]

        variant = 'ENSP00000269305.4:p.H175delinsX'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results['validation_warning_1']['validation_warnings'] == [
            "AminoMismatchError: The amino acid at position 175 of ENSP00000269305.4 is R not H"
        ]

        variant = 'ENSP00000269305.4:p.R175_H178delinsX'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results['validation_warning_1']['validation_warnings'] == [
            "ProteinSupportWarning: Protein level variant descriptions are not fully supported due to redundancy in the genetic code",
            "ProteinSupportWarning: ENSP00000269305.4:p.Arg175_His178delinsTer is HGVS compliant and contains a valid reference amino acid description"
        ]
        assert results['validation_warning_1'][
                   'hgvs_predicted_protein_consequence']["tlr"] == "ENSP00000269305.4:p.Arg175_His178delinsXaa"

        variant = 'ENSP00000269305.4:p.R175delinsX'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results['validation_warning_1']['validation_warnings'] == [
            "ProteinSupportWarning: Protein level variant descriptions are not fully supported due to redundancy in the genetic code",
            "ProteinSupportWarning: ENSP00000269305.4:p.Arg175delinsTer is HGVS compliant and contains a valid reference amino acid description"
        ]
        assert results['validation_warning_1'][
                   'hgvs_predicted_protein_consequence']["tlr"] == "ENSP00000269305.4:p.Arg175delinsXaa"

    def test_issue_432(self):
        variant = 'ENST00000318312.12:c.1779+7A>G'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)

        print(results)
        assert results['validation_warning_1']['validation_warnings'] == [
            "VariantMappingWarning: ENST00000318312.12:c.1779+7A>G auto-mapped to ENST00000318312.12:c.*4A>G",
            "ReferenceMismatchError: ENST00000318312.12:c.*4A>G: Variant reference (A) does not agree with reference sequence (C)"
        ]

    def test_issue_455(self):
        variant = 'ENSP00000269305.4:p.?'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)

        print(results)
        assert results['validation_warning_1']['validation_warnings'] == [
            "ProteinSupportWarning: Protein level variant descriptions are not fully supported due to redundancy in the genetic code",
            "ProteinSupportWarning: ENSP00000269305.4:p.? is HGVS compliant and contains a valid reference amino acid description"
        ]

    def test_issue_518a(self):
        variant = 'ENST00000636147.2(CLN3):c.791-802_1056+1445del'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)

        print(results)
        assert results['ENST00000636147.2:c.790+532_1056+1445del'][
                   'validation_warnings'] == [
            "VariantSyntaxError: Removing redundant gene symbol CLN3 from variant description",
            "ExonBoundaryError: Position c.791-802 has been updated to position to 790+532 ensuring correct HGVS "
            "numbering for transcript ENST00000636147.2"
        ]

    def test_missing_dot(self):
        variant = 'chr11:g,108121787G>A'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results[
            'validation_warning_1']["validation_warnings"] == ['VariantSyntaxError: Unable to identify a dot (.) in the variant description chr11:g,108121787G>A following the reference sequence type (g,c,n,r, or p). A dot is required in HGVS variant descriptions to separate the reference type from the variant position i.e. <accession>:<type>. e.g. :g.']

    def test_missing_colon(self):
        variant = 'chr11g.108121787G>A'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results[
            'validation_warning_1']["validation_warnings"] == ['VariantSyntaxError: Unable to identify a colon (:) in the variant description chr11g.108121787G>A. A colon is required in HGVS variant descriptions to separate the reference accession from the reference type i.e. <accession>:<type>. e.g. :c.', 'ReferenceSequenceError: This is not a valid HGVS variant description, because no reference sequence ID has been provided', 'ReferenceMismatchError: NC_000011.10:g.108121787G>A: Variant reference (G) does not agree with reference sequence (T)']

    def test_p1_a(self):
        variant = 'LRG_199p1:p.(Met1Ala)'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results[
                'validation_warning_1']["validation_warnings"] == ['VariantMappingWarning: LRG_199p1:p.(Met1Ala) automapped to equivalent RefSeq record NP_003997.1:p.(Met1Ala)', 'ProteinSupportWarning: Protein level variant descriptions are not fully supported due to redundancy in the genetic code', 'InitiationCodonWarning: Variant NP_003997.1:p.(Met1Ala) affects the initiation amino acid so is better described as NP_003997.1:p.(Met1?)']

    def test_p1_b(self):
        variant = 'LRG_199p1:p.Met1Ala'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results[
                'validation_warning_1']["validation_warnings"] == ['VariantMappingWarning: LRG_199p1:p.Met1Ala automapped to equivalent RefSeq record NP_003997.1:p.Met1Ala', 'ProteinSupportWarning: Protein level variant descriptions are not fully supported due to redundancy in the genetic code', 'InitiationCodonWarning: Variant NP_003997.1:p.Met1Ala affects the initiation amino acid so is better described as NP_003997.1:p.(Met1?)']

    def test_uppercase_ref_type(self):
        variant = 'DPYD:C.1905+1G>A'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert ('ReferenceTypeError: Reference type incorrectly stated in the variant description DPYD:C.1905+1G>A Valid types are g,c,n,r, or p') in results['validation_warning_1']["validation_warnings"]

    def test_invalid_aa(self):
        variant = 'ENST00000003084.11:p.Z1335P'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert ("ReferenceTypeError: Using a nucleotide reference sequence (NM_ NR_ NG_ NC_) to specify protein-level "
                "(p.) variation is not HGVS compliant. Please select an appropriate protein reference sequence (NP_)")

        variant = 'ENSP00000003084.11:p.Z1335P'
        results = self.vv.validate(variant, 'GRCh38', 'all', liftover_level='primary', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert "Invalid amino acid Z stated in description ENSP00000003084.11:p.Z1335P" in \
               results['validation_warning_1']['validation_warnings'][0]

    def test_g_with_tc_ref(self):
        variant = 'ENST00000225964.10:g.2559del'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results['validation_warning_1']['validation_warnings'] == ['ReferenceTypeError: Transcript reference sequence input as genomic (g.) reference sequence. Did you mean ENST00000225964.10:c.2559del?']

        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="refseq").format_as_dict(test=True)
        assert "InvalidFieldError: The transcript ENST00000225964.10 is not in the RefSeq data set. Please select Ensembl" in \
               results['validation_warning_1']['validation_warnings']


    def test_g_with_tc_ref_b(self):
        variant = 'ENST00000225964.10:g.2559+54_2560del'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results['validation_warning_1']['validation_warnings'] == ['ReferenceTypeError: Transcript reference sequence input as genomic (g.) reference sequence. Did you mean ENST00000225964.10:c.2559+54_2560del?']

    def test_p_with_tc_ref(self):
        variant = 'ENST00000225964.10:p.(Gly197Cys)'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert results['validation_warning_1']['validation_warnings'] == ['ReferenceTypeError: Using a nucleotide reference sequence (NM_ NR_ NC_ NG_ NT_ NW_) to specify protein-level (p.) variation is not HGVS compliant. Please select an appropriate protein reference sequence (NP_)']

    def test_invalid_reference_set_refseq(self):
        variant = 'ENST00000225964.10:g.2559+54_2560del'
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="refseq").format_as_dict(test=True)
        print(results)
        assert "InvalidFieldError: The transcript ENST00000225964.10 is not in the RefSeq data set. Please select Ensembl" in \
               results['validation_warning_1']['validation_warnings']
        results = self.vv.validate(variant, 'GRCh38', 'all', transcript_set="ensembl").format_as_dict(test=True)
        print(results)
        assert (results['validation_warning_1']['validation_warnings'] ==
                ['ReferenceTypeError: Transcript reference sequence input as genomic (g.) reference sequence. '
                 'Did you mean ENST00000225964.10:c.2559+54_2560del?'])

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


# <LICENSE>
# Copyright (C) 2016-2026 VariantValidator Contributors
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
