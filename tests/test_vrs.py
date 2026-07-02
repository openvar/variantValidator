"""
This test set tests on HGVS to VRS output conversion module to check first the
internal functions and then the overall outputs.

The extra test data, particularly test checksums, is taken from the known good
output given in the VRS documentation.
"""
import copy
import os
import unittest
import json
from configparser import ConfigParser

import vvhgvs
from vvhgvs.location import Interval
from biocommons.seqrepo import SeqRepo

from VariantValidator.modules.vrs_utils import HGVS_to_VRS
from VariantValidator.modules.complex_descriptions import FEInterval
from VariantValidator import Validator
from VariantValidator.modules.hgvs_utils import hgvs_obj_from_existing_edit,\
        _hgvs_offset_pos_from_str_in
from VariantValidator.settings import CONFIG_DIR

vv = Validator()
vv.alt_aln_method = "splign"


class MockVVData():
    "Mock class should be == to VV wrt to VRS output production"
    def __init__(
            self,
            selected_assembly = 'GRCh38', original = 'TestOrigVarText',
            warnings = None, lovd_messages = None, lovd_corrections = None,
            stable_gene_ids = None, description = '', gene_symbol = '',
            output_type_flag = None, rna_data = None,
            hgvs_transcript_variant = '', primary_assembly_loci = None,
            alt_genomic_loci = None, hgvs_predicted_protein_consequence = None,
            hgvs_refseqgene_variant = None, hgvs_genomic = None):
        self.selected_assembly = selected_assembly
        self.original = original # submitted_variant
        self.warnings = []
        if self.warnings:
            self.warnings = warnings
        self.lovd_messages = []
        if lovd_messages:
            self.lovd_messages = lovd_messages
        self.lovd_corrections = []
        if lovd_corrections:
            self.lovd_corrections = lovd_corrections
        self.hgvs_genomic = hgvs_genomic

        ## should be set to null if intergenic
        self.stable_gene_ids = []
        if stable_gene_ids:
            self.stable_gene_ids = stable_gene_ids

        self.description = description # transcript_description
        self.gene_symbol = gene_symbol

        # one of 'gene' 'warning' 'mitochondrial' or 'intergenic'
        self.output_type_flag = output_type_flag


        # now set individual hgvs datasets

        # rna_data should be set to null if not r type variant, or else a dict
        # containing: 'rna_variant', 'translation_slr', and 'usage_warnings'
        self.rna_data = rna_data

        # non r type transcript data
        self.hgvs_transcript_variant = hgvs_transcript_variant

        # mappings for transcript variant or main genomic variant given
        self.primary_assembly_loci = primary_assembly_loci
        self.alt_genomic_loci = alt_genomic_loci

        # prot variation
        self.hgvs_predicted_protein_consequence = \
                hgvs_predicted_protein_consequence
        # RSG variation
        self.hgvs_refseqgene_variant = hgvs_refseqgene_variant
        # VM / data provider
        self.hdp = vvhgvs.dataproviders.uta.connect(pooling=True)
        self.lose_vm = vvhgvs.variantmapper.VariantMapper(
                self.hdp,
                replace_reference=True,
                prevalidation_level=None
                ) # Variant mapper
    def process_warnings(self):
        return self.warnings

class TestVRSOutputTrim(unittest.TestCase):
    """
    Tests for the internal vrs_util.py module, particularly the HGVS_to_VRS
    object's trimming function which reduces a variant to its most minimal
    sequence length of differing bases.

    Tests are for:
        minimise totally,
        no change,
        minimise from partial,
        minimise from whole region,
        minimise from left, and
        minimise from right

    """
    @classmethod
    def setup_class(self):
        self.vrs = HGVS_to_VRS()

    def test_vrs_trim_total(self):
        ref_seq, alt_seq, prefix_bases_removed, sufix_bases_removed = \
                self.vrs._trim('AACAA','AACAA')
        assert ref_seq == ''
        assert alt_seq == ''
        assert not prefix_bases_removed and not sufix_bases_removed

    def test_vrs_trim_no_change(self):
        ref_seq, alt_seq, prefix_bases_removed, sufix_bases_removed = \
                self.vrs._trim('ACGT','CTG')
        assert ref_seq == 'ACGT'
        assert alt_seq == 'CTG'
        assert not prefix_bases_removed and not sufix_bases_removed

    def test_vrs_trim_both_edge_partial(self):
        ref_seq, alt_seq, prefix_bases_removed, sufix_bases_removed = \
                self.vrs._trim('AAAA','ACA')
        assert ref_seq == 'AA'
        assert alt_seq == 'C'
        assert prefix_bases_removed == 1
        assert sufix_bases_removed == -1
        ref_seq, alt_seq, prefix_bases_removed, sufix_bases_removed = \
                self.vrs._trim('ACA','AAAA')
        assert alt_seq == 'AA'
        assert ref_seq == 'C'
        assert prefix_bases_removed == 1
        assert sufix_bases_removed == -1

    def test_vrs_trim_both_edge_total(self):
        ref_seq, alt_seq, prefix_bases_removed, sufix_bases_removed = \
                self.vrs._trim('AAAA','AACAA')
        assert ref_seq == ''
        assert alt_seq == 'C'
        assert prefix_bases_removed == 2
        assert sufix_bases_removed == -2
        ref_seq, alt_seq, prefix_bases_removed, sufix_bases_removed = \
                self.vrs._trim('AACAA','AAAA')
        assert alt_seq == ''
        assert ref_seq == 'C'
        assert prefix_bases_removed == 2
        assert sufix_bases_removed == -2

    def test_vrs_trim_from_left(self):
        ref_seq, alt_seq, prefix_bases_removed, sufix_bases_removed = \
                self.vrs._trim('AAAA','CAA')
        assert ref_seq == 'AA'
        assert alt_seq == 'C'
        assert prefix_bases_removed == 0
        assert sufix_bases_removed == -2
        ref_seq, alt_seq, prefix_bases_removed, sufix_bases_removed = \
                self.vrs._trim('CAA','AAAA')
        assert ref_seq == 'C'
        assert alt_seq == 'AA'
        assert prefix_bases_removed == 0
        assert sufix_bases_removed == -2

    def test_vrs_trim_from_right(self):
        ref_seq, alt_seq, prefix_bases_removed, sufix_bases_removed = \
                self.vrs._trim('AAAA','AAC')
        assert ref_seq == 'AA'
        assert alt_seq == 'C'
        assert prefix_bases_removed == 2
        assert sufix_bases_removed == 0
        ref_seq, alt_seq, prefix_bases_removed, sufix_bases_removed = \
                self.vrs._trim('AAC','AAAA')
        assert ref_seq == 'C'
        assert alt_seq == 'AA'
        assert prefix_bases_removed == 2
        assert sufix_bases_removed == 0

class TestVRSOutputPush(unittest.TestCase):
    """
    Tests for the internal vrs_util.py modules push left & right used in
    normalisation of VRS input, the testes are:
        normal push,
        push already done,
        push from full opposed end,
        push with flank amount less than repeat unit,
        push where flank fetch amount == repeat unit length,
        push with flank fetch amount large enough to overlap transcript end, and
        push with fetch end == end coordinate
    known quirks:
        ref_left_flank/ref_right_flank input needs to be set to something (or
            we return as if we hit the end), for now, due to loop logic
        push does not care (or know) about the difference between del an ins
    push function data/key
        returns: ref_(left/right)_flank, start, bases_non_derived
            ref_(left/right)_flank, bases to add to the flank after roll
            start, new start location
            bases_non_derived are the no. of bases in the roll seq not in flank
        input: target_ac, start,ref_flank, roll_seq
            roll_seq is the minimised ins or del in question,
            ref_(left/right)_flank is the current left/right flank,
            roll_seq minimised (i.e. trimmed) version of input variant
    """
    @classmethod
    def setup_class(self):
        self.vrs = HGVS_to_VRS()
        self.longer_vrs = copy.copy(self.vrs)
        self.longer_vrs.ref_fetch_blocksize = 10
        self.mid_vrs = copy.copy(self.vrs)
        self.mid_vrs.ref_fetch_blocksize = 5
        self.short_vrs = copy.copy(self.vrs)
        self.short_vrs.ref_fetch_blocksize = 2
        self.ex_short_vrs = copy.copy(self.vrs)
        self.ex_short_vrs.ref_fetch_blocksize = 1
        # push left data:
        # NM_003073.5
        # cds start 204, seq AGAAGAAGA start 1288:1297
        # flank 'TGAGATGG' 'AGAAGAAGA' 'TCCGCGAC'
        # this has a hanging part repeat of either CA or GC depending VRS
        # normalisation should get the whole +partial when hgvs does not
        # NM_002111.8
        # cds start 145 seq:
        # CA GCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA
        # start 198:261 (or rather 196:261)
        # with flank is:
        # CAAGTCCTTC
        # CA GCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA
        # ACAGCCGCCA
        # longer mid seq (not used atm)
        # NM_001350922.2:c.240+14TGGG[2] NC_000010.10:g.128358789_128358796=
        # start seq
        # NM_001330096.1
        # rep CT stop 26
        # seq CTCTCTCTCTCTCTCTCTCTCTCTCT TTCTTTCTTTCCGATGGGGAAGAGAGGGT
        # NM_001330112.1
        # rep GGGGCGGGTC stop 28 (or 20 for complete segments)
        # seq (with flank) GGGGCGGGTCGGGGCGGGTC GGGGCGGG CCCACGGGCGGCCGGATTTG
        # end seq
        # NM_052928.3
        # rep GTTTT start 4375 len 4405
        # seq (w flank):
        # ATAACACATATGCCTCCTTCTGAGTTGTTG GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT
        # NM_002129.4
        # rep AT start 1408 len 1441
        # seq ATTAGCTAAATTGTTCCTCAGGTGTGTG T ATATATATATACATATATATATATATATATAT
    def test_vrs_push_left_from_mid(self):
        # push left from middle
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_l('NR_110010.2',27,'C','GC')
        assert start == 25
        assert bases_non_derived == 0
        assert ref_left_flank == 'GC'
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_l('NM_003073.5',1291,'AGA','AGA')
        assert start == 1288
        assert bases_non_derived == 0
        assert ref_left_flank == 'AGA'
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_l('NM_002111.8',202,'AG','CAG')
        assert start == 196
        assert bases_non_derived == 0
        assert ref_left_flank == 'CAGCAG'

    def test_vrs_push_left_already_done(self):
        # push left already done (from far push base)
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_l('NR_110010.2',25,'A','GC')
        assert start == 25
        assert bases_non_derived == 2
        assert ref_left_flank == ''
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_l('NM_003073.5',1288,'GG','AGA')
        assert start == 1288
        assert bases_non_derived == 3
        assert ref_left_flank == ''
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_l('NM_002111.8',196,'CAAGTCCTTC','CAG')
        assert start == 196
        assert bases_non_derived == 3
        assert ref_left_flank == ''

    def test_vrs_push_left_from_full_right(self):
        # push left from full right
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_l('NR_110010.2',32,'G','CG')
        assert start == 25
        assert bases_non_derived == 0
        assert ref_left_flank == 'GCGCGCG'
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_l('NM_003073.5',1297,'A','AGA')
        assert start == 1288
        assert bases_non_derived == 0
        assert ref_left_flank == 'AGAAGAAGA'
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_l('NM_002111.8',261,'CA','GCA')
        assert start == 196
        assert bases_non_derived == 0
        assert ref_left_flank == \
            'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA'

    def test_vrs_push_left_from_full_right_w_fetch(self):
        # push left from full right (pre fetched flank less than needed, to force internal fetch)
        ref_left_flank, start, bases_non_derived = \
                self.longer_vrs._push_l('NR_110010.2',32,'G','CG')
        assert start == 25
        assert bases_non_derived == 0
        assert ref_left_flank == 'GCGCGCG'
        ref_left_flank, start, bases_non_derived = \
                self.longer_vrs._push_l('NM_003073.5',1297,'AGA','AGA')
        assert start == 1288
        assert bases_non_derived == 0
        assert ref_left_flank == 'AGAAGAAGA'

    def test_vrs_push_left_fetch_lt_repeat_size(self):
        # push left given flank to fetch amount less than the size repeat unit (roll_seq)
        ref_left_flank, start, bases_non_derived = \
                self.ex_short_vrs._push_l('NR_110010.2',32,'G','CG')
        assert start == 25
        assert bases_non_derived == 0
        assert ref_left_flank == 'GCGCGCG'
        ref_left_flank, start, bases_non_derived = \
                self.short_vrs._push_l('NM_003073.5',1297,'GA','AGA')
        assert start == 1288
        assert bases_non_derived == 0
        assert ref_left_flank == 'AGAAGAAGA'
        ref_left_flank, start, bases_non_derived = \
                self.short_vrs._push_l('NM_001330112.1',10,'GTC', 'GGGGCGGGTC')
        assert start == 0
        assert bases_non_derived == 0
        assert ref_left_flank == 'GGGGCGGGTC'

    def test_vrs_push_left_fetch_eq_repeat_size(self):
        # push left, flank to fetch amount == repeat unit length (roll_seq)
        ref_left_flank, start, bases_non_derived = \
                self.short_vrs._push_l('NR_110010.2',32,'CG','CG')
        assert start == 25
        assert bases_non_derived == 0
        assert ref_left_flank == 'GCGCGCG'

    def test_vrs_push_left_partial_roll(self):
        # test partial roll
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_l('NR_110010.2',27,'C','GCTGC')
        assert start == 25
        assert bases_non_derived == 3
        assert ref_left_flank == 'GC'
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_l('NM_003073.5',1291,'AGA','AGACAGA')
        assert start == 1288
        assert bases_non_derived == 4
        assert ref_left_flank == 'AGA'
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_l('NM_002111.8',199,'CAG','CAGTCAG')
        assert start == 196
        assert bases_non_derived == 4
        assert ref_left_flank == 'CAG'

    def test_vrs_push_left_hits_seq_start(self):
        # test hits start
        # NM_001330096.1
        # rep CT stop 26
        # seq CTCTCTCTCTCTCTCTCTCTCTCTCT TTCTTTCTTTCCGATGGGGAAGAGAGGGT
        # NM_001330112.1
        # rep GGGGCGGGTC stop 28 (or 20 for complete segments)
        # seq (with flank) GGGGCGGGTCGGGGCGGGTC GGGGCGGG CCCACGGGCGGCCGGATTTG
        # test corner case roll fetch hits 0 and rolls to 0
        ref_left_flank, start, bases_non_derived = \
                self.short_vrs._push_l('NM_001330096.1',4,'CT','CT')
        assert start == 0
        assert bases_non_derived == 0
        assert ref_left_flank == 'CTCT'
        ref_left_flank, start, bases_non_derived = self.longer_vrs._push_l(
                'NM_001330112.1',20,'GGGGCGGGTC','GGGGCGGGTC')
        assert start == 0
        assert bases_non_derived == 0
        assert ref_left_flank == 'GGGGCGGGTCGGGGCGGGTC'
        # test corner cases: roll fetch crosses 0 (truncated final fetch)
        # and rolls to 0
        ref_left_flank, start, bases_non_derived = \
                self.longer_vrs._push_l('NM_001330096.1',4,'CT','CT')
        assert start == 0
        assert bases_non_derived == 0
        assert ref_left_flank == 'CTCT'
        ref_left_flank, start, bases_non_derived = self.longer_vrs._push_l(
                'NM_001330112.1',20,'GGGTCGGGGCGGGTC','GGGGCGGGTC')
        assert start == 0
        assert bases_non_derived == 0
        assert ref_left_flank == 'GGGGCGGGTCGGGGCGGGTC'

        # Push right data
        # end seq
        # NM_052928.3
        # rep GTTTT start 4375 len 4405
        # seq (w flank):
        # ATAACACATATGCCTCCTTCTGAGTTGTTG GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT
        # NM_002129.4
        # rep AT start 1408 len 1441
        # seq ATTAGCTAAATTGTTCCTCAGGTGTGTG T ATATATATATACATATATATATATATATATAT
        # NR_110010.2
        # seq start 25:32
        # seq w flank ATGGTGGCGA GCGCGCG AGTGCAGAAG
        # NM_003073.5
        # cds start 204 start 1288:1297
        # seq w flank TGAGATGG AGAAGAAGA TCCGCGAC
        # this (NM_003073.5 rep) has a hanging partial repeat, either CA or GC
        # depending, VRS normalisation should get whole +partial, hgvs does not
        # NM_002111.8
        # cds start 145 start 198:261 or rather 196:261
        # seq with flank:
        # CAAGTCCTTC
        # CA GCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA
        # ACAGCCGCCA
        # longer mid seq (not used atm)
        # NM_001350922.2:c.240+14TGGG[2] NC_000010.10:g.128358789_128358796=
        # start seq
        # NM_001330096.1
        # rep CT stop 26
        # seq CTCTCTCTCTCTCTCTCTCTCTCTCT TTCTTTCTTTCCGATGGGGAAGAGAGGGT
        # NM_001330112.1
        # rep GGGGCGGGTC stop 28 (or 20 for complete segments)
        # seq (with flank) GGGGCGGGTCGGGGCGGGTC GGGGCGGG CCCACGGGCGGCCGGATTTG
        # end seq
        # NM_052928.3
        # rep GTTTT start 4375 len 4405
        # seq (w flank):
        # ATAACACATATGCCTCCTTCTGAGTTGTTG GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT
        # NM_002129.4
        # rep AT start 1408 len 1441
        # seq ATTAGCTAAATTGTTCCTCAGGTGTGTG T ATATATATATAC ATATATATATATATATATAT

    def test_vrs_push_right_from_middle(self):
        # push right from middle
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_r('NR_110010.2',27,'G','GC')
        assert start == 32
        assert bases_non_derived == 0
        assert ref_left_flank == 'GCGCG'
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_r('NM_003073.5',1291,'AGA','AGA')
        assert start == 1297
        assert bases_non_derived == 0
        assert ref_left_flank == 'AGAAGA'
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_r('NM_002111.8',202,'CA','CAG')
        assert start == 261
        assert bases_non_derived == 0
        assert ref_left_flank == \
                'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA'

    def test_vrs_push_right_already_done(self):
        # push right already done (from far push base)
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_r('NR_110010.2',32,'A','GC')
        assert start == 32
        assert bases_non_derived == 2
        assert ref_left_flank == ''
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_r('NM_003073.5',1297,'TCCGCGAC','AGA')
        assert start == 1297
        assert bases_non_derived == 3
        assert ref_left_flank == ''
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_r('NM_002111.8',261,'ACAGCCGCCA','CAG')
        assert start == 261
        assert bases_non_derived == 3
        assert ref_left_flank == ''

    def test_vrs_push_right_from_full_left(self):
        # push right from full left
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_r('NR_110010.2',25,'G','GC')
        assert start == 32
        assert bases_non_derived == 0
        assert ref_left_flank == 'GCGCGCG'
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_r('NM_003073.5',1288,'A','AGA')
        assert start == 1297
        assert bases_non_derived == 0
        assert ref_left_flank == 'AGAAGAAGA'
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_r('NM_002111.8',196,'CA','CAG')
        assert start == 261
        assert bases_non_derived == 0
        assert ref_left_flank == \
            'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA'

    def test_vrs_push_right_from_full_left_w_forced_fetch(self):
        # push right from full left (pre fetched flank less than needed, to force internal fetch)
        ref_left_flank, start, bases_non_derived = \
                self.longer_vrs._push_r('NR_110010.2',25,'G','GC')
        assert start == 32
        assert bases_non_derived == 0
        assert ref_left_flank == 'GCGCGCG'
        ref_left_flank, start, bases_non_derived = \
                self.longer_vrs._push_r('NM_003073.5',1288,'AGA','AGA')
        assert start == 1297
        assert bases_non_derived == 0
        assert ref_left_flank == 'AGAAGAAGA'

    def test_vrs_push_right_fetch_lt_repeat_unit(self):
        # push right given flank to fetch amount less than the size repeat unit (roll_seq)
        ref_left_flank, start, bases_non_derived = \
                self.ex_short_vrs._push_r('NR_110010.2',25,'G','GC')
        assert start == 32
        assert bases_non_derived == 0
        assert ref_left_flank == 'GCGCGCG'
        ref_left_flank, start, bases_non_derived = \
                self.short_vrs._push_r('NM_003073.5',1288,'AG','AGA')
        assert start == 1297
        assert bases_non_derived == 0
        assert ref_left_flank == 'AGAAGAAGA'
        ref_left_flank, start, bases_non_derived = \
                self.short_vrs._push_r('NM_001330112.1',0,'GG', 'GGGGCGGGTC')
        assert start == 28
        assert bases_non_derived == 0
        assert ref_left_flank == 'GGGGCGGGTCGGGGCGGGTCGGGGCGGG'

    def test_vrs_push_right_fetch_eq_repeat_unit(self):
        # push right, flank to fetch amount == repeat unit length (roll_seq)
        ref_left_flank, start, bases_non_derived = \
                self.short_vrs._push_r('NR_110010.2',25,'GC','GC')
        assert start == 32
        assert bases_non_derived == 0
        assert ref_left_flank == 'GCGCGCG'

    def test_vrs_push_right_partial_roll(self):
        # test partial roll
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_r('NR_110010.2',27,'G','GCTGC')
        assert start == 29
        assert bases_non_derived == 3
        assert ref_left_flank == 'GC'
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_r('NM_003073.5',1288,'AGA','AGACAGA')
        assert start == 1291
        assert bases_non_derived == 4
        assert ref_left_flank == 'AGA'
        ref_left_flank, start, bases_non_derived = \
                self.vrs._push_r('NM_002111.8',196,'CAG','CAGTCAG')
        assert start == 199
        assert bases_non_derived == 4
        assert ref_left_flank == 'CAG'

    def test_vrs_push_right_hits_end(self):
        # test hits end
        # NM_052928.3
        # rep GTTTT start 4375 len 4405
        # seq (w flank):
        # ATAACACATATGCCTCCTTCTGAGTTGTTG GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT
        # NM_002129.4
        # rep AT start 1408 len 1441
        # seq ATTAGCTAAATTGTTCCTCAGGTGTGTG T ATATATATATAC ATATATATATATATATATAT'
        # test corner case roll fetch hits last base and rolls to last base
        ref_left_flank, start, bases_non_derived = \
                self.mid_vrs._push_r('NM_052928.3',4375,'GTT','GTTTT')
        assert start == 4405
        assert bases_non_derived == 0
        assert ref_left_flank == 'GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT'
        ref_left_flank, start, bases_non_derived = \
                self.short_vrs._push_r('NM_002129.4',1421,'AT','AT')
        assert start == 1441
        assert bases_non_derived == 0
        assert ref_left_flank == 'ATATATATATATATATATAT'
    def test_vrs_push_right_hits_end_handle_last_fetch_gt_final_base(self):
        # test corner case roll fetch crosses last base (truncated final fetch)
        # and rolls to last base
        ref_left_flank, start, bases_non_derived = \
                self.longer_vrs._push_r('NM_052928.3',4375,'GTT','GTTTT')
        assert start == 4405
        assert bases_non_derived == 0
        assert ref_left_flank == 'GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT'
        ref_left_flank, start, bases_non_derived = \
                self.longer_vrs._push_r('NM_002129.4',1421,'AT','AT')
        assert start == 1441
        assert bases_non_derived == 0
        assert ref_left_flank == 'ATATATATATATATATATAT'

class TestVRSOutputNorm(unittest.TestCase):
    """
    Tests for the internal vrs_util.py module's HGVS_to_VRS object's
    normalise_outwards function, which should minimise and then push outwards
    (in both directions) test structure as for the push left/right above.
    known quirks:
        Because the first step is minimise, followed by rolling the seq in both
            directions the actual minimised seq can differ in roll position
            e.g. CG vs GC. But otherwise for the same underlying variant the
            results should always be the same,
    push function data/key:
        returns: full_ref, full_alt,min_ref,min_alt,vrs_start,vrs_end,derived
             full_ref: normalised (maximised) version of the ref
             full_alt: normalised (maximised) version of the alt
             min_ref & min_alt: minimised version of the same, (used along with
                derived_state for variant type detection)
             derived_state: True/False, is the variant derived from the
                immediate flank (eg cnv/tandem repeat)
        input: ref_ac, ref_seq,alt_seq, vrs_start, vrs_end
            ref_ac is the accession that the change is relevant for
            ref_seq is the original sequence of the region
            alt_seq is the sequence of the variant in question
            vrs_start & vrs_end are the 0 based VRS type coordinates of the ref
                location
    tests done:
        Done for both del and ins:
        expand left, expand right, expand both directions, expand already done,
        expand flank fetch amount less than repeat unit, expand flank amount
        same as repeat unit, expand flank fetch amount large enough to overlap
        transcript start/end, expand flank amount == start/end coordinate from
        current

        Also test early abort on == data.
   """
    @classmethod
    def setup_class(self):
        self.vrs = HGVS_to_VRS()
        self.ex_short_vrs = copy.copy(self.vrs)
        self.ex_short_vrs.ref_fetch_blocksize = 1
        # NR_110010.2
        # seq start 25:32
        # seq with flank ATGGTGGCGA GCGCGCG AGTGCAGAAG
        # NM_003073.5
        # cds start 204 start 1288:1297
        # seq with flank TGAGATGG AGAAGAAGA TCCGCGAC
        # this has a hanging part repeat of either CA or GC depending
        # VRS normalisation should get the whole +partial when hgvs does not
        # NM_002111.8
        # cds start 145 start 198:261 or rather 196:261
        # seq with flank:
        # CAAGTCCTTC
        # CA GCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA
        # ACAGCCGCCA
        # longer mid seq (not used atm)
        # NM_001350922.2:c.240+14TGGG[2] NC_000010.10:g.128358789_128358796=
        # start seq
        # NM_001330096.1
        # rep CT stop 26
        # seq CTCTCTCTCTCTCTCTCTCTCTCTCT TTCTTTCTTTCCGATGGGGAAGAGAGGGT
        # NM_001330112.1
        # rep GGGGCGGGTC stop 28 (or 20 for complete segments)
        # seq (with flank) GGGGCGGGTCGGGGCGGGTC GGGGCGGG CCCACGGGCGGCCGGATTTG
        # end seq
        # NM_052928.3
        # rep GTTTT start 4375 len 4405
        # seq (with flank):
        # ATAACACATATGCCTCCTTCTGAGTTGTTG GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT
        # NM_002129.4
        # rep AT start 1408 len 1441
        # seq ATTAGCTAAATTGTTCCTCAGGTGTGTG T ATATATATATACATATATATATATATATATAT

        # result set key:
        # full_ref, full_alt, min_ref, min_alt, vrs_start, vrs_end, derived
        self.NR_110010_2_ins_result_set = (
                'GCGCGCG','GCGCGCGCG','','CG',25,32,True)
        self.NR_110010_2_ins_result_set_end = (
                'GCGCGCG','GCGCGCGCG','','GC',25,32,True)
        self.NM_003073_5_ins_result_set = (
                'AGAAGAAGA','AGAAGAAGAAGA','','AGA',1288,1297,True)
        self.NM_002111_8_ins_result_set = (
         'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA',
         'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA',
         '','CAG',196,261,True)
        self.NM_002111_8_ins_result_set_end = (
         'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA',
         'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA',
         '','GCA',196,261,True)
        # start seq
        # NM_001330096.1
        # rep CT stop 26
        # seq CTCTCTCTCTCTCTCTCTCTCTCTCT TTCTTTCTTTCCGATGGGGAAGAGAGGGT
        # NM_001330112.1
        # rep GGGGCGGGTC stop 28 (or 20 for complete segments)
        # seq (with flank) GGGGCGGGTCGGGGCGGGTC GGGGCGGG CCCACGGGCGGCCGGATTTG
        # end seq
        # NM_052928.3
        # rep GTTTT start 4375 len 4405
        # seq (w flank):
        # ATAACACATATGCCTCCTTCTGAGTTGTTG GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT
        # NM_002129.4
        # rep AT start 1408 len 1441
        # seq ATTAGCTAAATTGTTCCTCAGGTGTGTG T ATATATATATACATATATATATATATATATAT
        self.res_set_start_NM_001330096_1 = (
                'CTCTCTCTCTCTCTCTCTCTCTCTCT',
                'CTCTCTCTCTCTCTCTCTCTCTCTCTCT',
                '','TC',0,26,True)
        self.res_set_start_NM_001330096_1_end = (
                'CTCTCTCTCTCTCTCTCTCTCTCTCT',
                'CTCTCTCTCTCTCTCTCTCTCTCTCTCT',
                '','CT',0,26,True)
        self.res_set_start_NM_001330112_1 = (
                'GGGGCGGGTCGGGGCGGGTCGGGGCGGG',
                'GGGGCGGGTCGGGGCGGGTCGGGGCGGGTCGGGGCGGG',
                '','GGGGCGGGTC',0,28,True)
        self.res_set_start_NM_001330112_1_end = (
                'GGGGCGGGTCGGGGCGGGTCGGGGCGGG',
                'GGGGCGGGTCGGGGCGGGTCGGGGCGGGTCGGGGCGGG',
                '','TCGGGGCGGG',0,28,True)

        self.res_set_end_NM_052928_3 = (
                'GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT',
                'GTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTT',
                '','GTTTT',4375,4405,True)
        self.res_set_end_NM_052928_3_end = (
                'GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT',
                'GTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTT',
                '','TTGTT' ,4375,4405,True)
        self.res_set_end_NM_002129_4 = (
                'ATATATATATATATATATAT',
                'ATATATATATATATATATATAT',
                '','AT',1421,1441,True)
        self.res_set_end_NM_002129_4_end = (
                'ATATATATATATATATATAT',
                'ATATATATATATATATATATAT',
                '','TA',1421,1441,True)

        # del
        self.NR_110010_2_del_result_set = (
                'GCGCGCG','GCGCG','CG','',25,32,True)
        self.NR_110010_2_del_result_set_end = (
                'GCGCGCG','GCGCG','GC','',25,32,True)
        self.NM_003073_5_del_result_set = (
                'AGAAGAAGA','AGAAGA','AGA','',1288,1297,True)
        self.NM_002111_8_del_result_set = (
            'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA',
            'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA',
            'CAG','',196,261,True)
        self.NM_002111_8_del_result_set_end = (
            'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA',
            'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA',
            'GCA','',196,261,True)
        # start seq
        # NM_001330096.1
        # rep CT stop 26
        # seq CTCTCTCTCTCTCTCTCTCTCTCTCT TTCTTTCTTTCCGATGGGGAAGAGAGGGT
        # NM_001330112.1
        # rep GGGGCGGGTC stop 28 (or 20 for complete segments)
        # seq (with flank) GGGGCGGGTCGGGGCGGGTC GGGGCGGG CCCACGGGCGGCCGGATTTG
        # end seq
        # NM_052928.3
        # rep GTTTT start 4375 len 4405
        # seq (with flank):
        # ATAACACATATGCCTCCTTCTGAGTTGTTG GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT
        # NM_002129.4
        # rep AT start 1408 len 1441 seq
        # ATTAGCTAAATTGTTCCTCAGGTGTGTG T ATATATATATACATATATATATATATATATAT
        self.res_set_d_start_NM_001330096_1 = (
                'CTCTCTCTCTCTCTCTCTCTCTCTCT',
                'CTCTCTCTCTCTCTCTCTCTCTCT',
                'TC','',0,26,True)
        self.res_set_d_start_NM_001330096_1_end = (
                'CTCTCTCTCTCTCTCTCTCTCTCTCT',
                'CTCTCTCTCTCTCTCTCTCTCTCT',
                'CT','',0,26,True)
        self.res_set_d_start_NM_001330112_1 = (
                'GGGGCGGGTCGGGGCGGGTCGGGGCGGG',
                'GGGGCGGGTCGGGGCGGG',
                'GGGGCGGGTC','',0,28,True)
        self.res_set_d_start_NM_001330112_1_end = (
                'GGGGCGGGTCGGGGCGGGTCGGGGCGGG',
                'GGGGCGGGTCGGGGCGGG',
                'TCGGGGCGGG','',0,28,True)

        self.res_set_d_end_NM_052928_3 = (
                'GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT',
                'GTTTTGTTTTGTTTTGTTTTGTTTT',
                'GTTTT', '', 4375,4405,True)
        self.res_set_d_end_NM_052928_3_end = (
                'GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT',
                'GTTTTGTTTTGTTTTGTTTTGTTTT',
                'TTGTT', '', 4375,4405,True)
        self.res_set_d_end_NM_002129_4 = (
                'ATATATATATATATATATAT',
                'ATATATATATATATATAT',
                'AT','',1421,1441,True)
        self.res_set_d_end_NM_002129_4_end = (
                'ATATATATATATATATATAT',
                'ATATATATATATATATAT',
                'TA','',1421,1441,True)

    def test_vrs_normalise_outwards(self):
        # == return un-trimmed
        res = self.vrs.normalise_outwards('NR_110010.2', 'CG','CG', 31,33)
        NR_110010_2_eq_set = ('CG','CG','CG','CG',31,33,0)
        assert res == NR_110010_2_eq_set

    def test_vrs_normalise_outwards_ex_left_ins(self):
        # expand left ins
        res = self.vrs.normalise_outwards('NR_110010.2', '','CG', 32,32)
        assert res == self.NR_110010_2_ins_result_set
        res = self.vrs.normalise_outwards('NR_110010.2', 'G','GCG', 31,32)
        assert res == self.NR_110010_2_ins_result_set_end
        res = self.vrs.normalise_outwards('NM_003073.5', '','AGA', 1297,1297)
        assert res == self.NM_003073_5_ins_result_set
        res = self.vrs.normalise_outwards(
                'NM_003073.5', 'AGA','AGAAGA', 1294,1297)
        assert res == self.NM_003073_5_ins_result_set
        res = self.vrs.normalise_outwards('NM_002111.8', '','GCA', 261,261)
        assert res == self.NM_002111_8_ins_result_set_end
        res = self.vrs.normalise_outwards('NM_002111.8', 'CA','CAGCA', 259,261)
        assert res == self.NM_002111_8_ins_result_set

    def test_vrs_normalise_outwards_ex_right_ins(self):
        # expand right ins
        res = self.vrs.normalise_outwards('NR_110010.2', '','GC', 25,25)
        assert res == self.NR_110010_2_ins_result_set_end
        res = self.vrs.normalise_outwards('NR_110010.2', 'G','GCG', 25,26)
        assert res == self.NR_110010_2_ins_result_set_end
        res = self.vrs.normalise_outwards('NM_003073.5', '','AGA', 1288,1288)
        assert res == self.NM_003073_5_ins_result_set
        res = self.vrs.normalise_outwards(
                'NM_003073.5', 'AGA','AGAAGA', 1288,1291)
        assert res == self.NM_003073_5_ins_result_set
        res = self.vrs.normalise_outwards('NM_002111.8', '','CAG', 196,196)
        assert res == self.NM_002111_8_ins_result_set
        res = self.vrs.normalise_outwards('NM_002111.8', 'CA','CAGCA', 196,198)
        assert res == self.NM_002111_8_ins_result_set

    def test_vrs_normalise_outwards_ex_both_ins(self):
        # expand both directions ins
        res = self.vrs.normalise_outwards('NR_110010.2', '','GC', 27,27)
        assert res == self.NR_110010_2_ins_result_set_end
        res = self.vrs.normalise_outwards('NR_110010.2', 'G','GCG', 27,28)
        assert res == self.NR_110010_2_ins_result_set_end
        res = self.vrs.normalise_outwards('NM_003073.5', '','AGA', 1291,1291)
        assert res == self.NM_003073_5_ins_result_set
        res = self.vrs.normalise_outwards(
                'NM_003073.5', 'AGA','AGAAGA', 1291,1294)
        assert res == self.NM_003073_5_ins_result_set
        res = self.vrs.normalise_outwards('NM_002111.8', '','CAG', 202,202)
        assert res == self.NM_002111_8_ins_result_set
        res = self.vrs.normalise_outwards('NM_002111.8', 'CA','CAGCA', 202,204)
        assert res == self.NM_002111_8_ins_result_set

    def test_vrs_normalise_outwards_ex_done_ins(self):
        # expand already done ins
        res = self.vrs.normalise_outwards(
                'NR_110010.2', 'GCGCGCG','GCGCGCGCG', 25,32)
        assert res == self.NR_110010_2_ins_result_set_end
        res = self.vrs.normalise_outwards(
                'NM_003073.5', 'AGAAGAAGA','AGAAGAAGAAGA', 1288,1297)
        assert res == self.NM_003073_5_ins_result_set
        res = self.vrs.normalise_outwards(
         'NM_002111.8',
         'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA',
         'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA',
         196,261)
        assert res == self.NM_002111_8_ins_result_set

    def test_vrs_normalise_outwards_ex_into_start_ins(self):
        # expand into start ins
        res = self.vrs.normalise_outwards('NM_001330096.1', '','CT', 24,24)
        assert res == self.res_set_start_NM_001330096_1_end
        res = self.vrs.normalise_outwards('NM_001330096.1', 'T','TCT', 23,24)
        assert res == self.res_set_start_NM_001330096_1
        res = self.vrs.normalise_outwards(
                'NM_001330112.1', '','GGGGCGGGTC', 10,10)
        assert res == self.res_set_start_NM_001330112_1
        res = self.vrs.normalise_outwards(
                'NM_001330112.1', 'TC','TCGGGGCGGGTC', 8,10)
        assert res == self.res_set_start_NM_001330112_1_end

    def test_vrs_normalise_outwards_ex_into_start_done_ins(self):
        # expand start already done ins
        res = self.vrs.normalise_outwards(
                'NM_001330096.1',
                'CTCTCTCTCTCTCTCTCTCTCTCTCT',
                'CTCTCTCTCTCTCTCTCTCTCTCTCTCT',
                0,26)
        assert res == self.res_set_start_NM_001330096_1_end
        res = self.vrs.normalise_outwards(
                'NM_001330112.1',
                'GGGGCGGGTCGGGGCGGGTCGGGGCGGG',
                'GGGGCGGGTCGGGGCGGGTCGGGGCGGGTCGGGGCGGG',
                0,28)
        assert res == self.res_set_start_NM_001330112_1


    def test_vrs_normalise_outwards_ex_into_end_ins(self):
        # expand into end ins
        res = self.vrs.normalise_outwards('NM_052928.3', '','GTTTT', 4380,4380)
        assert res == self.res_set_end_NM_052928_3
        res = self.vrs.normalise_outwards(
                'NM_052928.3', 'TT','TTGTTTT', 4378,4380)
        assert res == self.res_set_end_NM_052928_3_end
        res = self.vrs.normalise_outwards('NM_002129.4', '','TA', 1422,1422)
        assert res == self.res_set_end_NM_002129_4_end
        res = self.vrs.normalise_outwards('NM_002129.4', 'TA','TATA', 1422,1424)
        assert res == self.res_set_end_NM_002129_4_end

    def test_vrs_normalise_outwards_ex_into_end_done_ins(self):
        # expand end already done ins
        res = self.vrs.normalise_outwards(
                'NM_052928.3',
                'GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT',
                'GTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTT',
                4375,4405)
        assert res == self.res_set_end_NM_052928_3
        res = self.vrs.normalise_outwards(
                'NM_002129.4',
                'ATATATATATATATATATAT',
                'ATATATATATATATATATATAT',
                1421,1441)
        assert res == self.res_set_end_NM_002129_4


        ##DEL SET##
    def test_vrs_normalise_outwards_ex_left_del(self):
        # expand left del
        res = self.vrs.normalise_outwards('NR_110010.2', 'CG','', 30,32)
        assert res == self.NR_110010_2_del_result_set
        res = self.vrs.normalise_outwards('NR_110010.2', 'GCG','G', 29,32)
        assert res == self.NR_110010_2_del_result_set_end
        res = self.vrs.normalise_outwards('NM_003073.5', 'AGA','', 1294,1297)
        assert res == self.NM_003073_5_del_result_set
        res = self.vrs.normalise_outwards(
                'NM_003073.5', 'AGAAGA','AGA', 1291,1297)
        assert res == self.NM_003073_5_del_result_set
        res = self.vrs.normalise_outwards('NM_002111.8', 'GCA','', 258,261)
        assert res == self.NM_002111_8_del_result_set_end
        res = self.vrs.normalise_outwards('NM_002111.8', 'CAGCA','CA', 256,261)
        assert res == self.NM_002111_8_del_result_set

    def test_vrs_normalise_outwards_ex_left_del_flank_sepfetch(self):
        # repeat left shift with self.ex_short_vrs to test flank fetch with:
        # minimised ref length > 2*flank fetch
        # this is used as the cut-off for fetching flanks separately
        res = self.ex_short_vrs.normalise_outwards(
                'NM_003073.5', 'AGA','', 1294,1297)
        assert res == self.NM_003073_5_del_result_set
        res = self.ex_short_vrs.normalise_outwards(
                'NM_003073.5', 'AGAAGA','AGA', 1291,1297)
        assert res == self.NM_003073_5_del_result_set
        res = self.ex_short_vrs.normalise_outwards(
                'NM_002111.8', 'GCA','', 258,261)
        assert res == self.NM_002111_8_del_result_set_end
        res = self.ex_short_vrs.normalise_outwards(
                'NM_002111.8', 'CAGCA','CA', 256,261)
        assert res == self.NM_002111_8_del_result_set

    def test_vrs_normalise_outwards_ex_right_del(self):
        # expand right del
        res = self.vrs.normalise_outwards('NR_110010.2', 'GC','', 25,27)
        assert res == self.NR_110010_2_del_result_set_end
        res = self.vrs.normalise_outwards('NR_110010.2', 'GCG','G', 25,28)
        assert res == self.NR_110010_2_del_result_set_end
        res = self.vrs.normalise_outwards('NM_003073.5', 'AGA','', 1288,1291)
        assert res == self.NM_003073_5_del_result_set
        res = self.vrs.normalise_outwards(
                'NM_003073.5', 'AGAAGA','AGA', 1288,1294)
        assert res == self.NM_003073_5_del_result_set
        res = self.vrs.normalise_outwards('NM_002111.8', 'CAG','', 196,199)
        assert res == self.NM_002111_8_del_result_set
        res = self.vrs.normalise_outwards('NM_002111.8', 'CAGCA','CA', 196,201)
        assert res == self.NM_002111_8_del_result_set

    def test_vrs_normalise_outwards_ex_right_del_flank_sepfetch(self):
        # repeat left shift with self.ex_short_vrs to test flank fetch with:
        # minimised ref length > 2*flank fetch
        # this is used as the cut-off for fetching flanks separately
        res = self.ex_short_vrs.normalise_outwards(
                'NM_003073.5', 'AGA','', 1288,1291)
        assert res == self.NM_003073_5_del_result_set
        res = self.ex_short_vrs.normalise_outwards(
                'NM_003073.5', 'AGAAGA','AGA', 1288,1294)
        assert res == self.NM_003073_5_del_result_set
        res = self.ex_short_vrs.normalise_outwards(
                'NM_002111.8', 'CAG','', 196,199)
        assert res == self.NM_002111_8_del_result_set
        res = self.ex_short_vrs.normalise_outwards(
                'NM_002111.8', 'CAGCA','CA', 196,201)
        assert res == self.NM_002111_8_del_result_set

    def test_vrs_normalise_outwards_ex_both_del(self):
        # expand both directions del
        res = self.vrs.normalise_outwards('NR_110010.2', 'GC','', 27,29)
        assert res == self.NR_110010_2_del_result_set_end
        res = self.vrs.normalise_outwards('NR_110010.2', 'GCG','G', 27,30)
        assert res == self.NR_110010_2_del_result_set_end
        res = self.vrs.normalise_outwards('NM_003073.5', 'AGA','', 1291,1294)
        assert res == self.NM_003073_5_del_result_set
        res = self.vrs.normalise_outwards('NM_003073.5', 'AGAAGA','AGA', 1291,1297)
        assert res == self.NM_003073_5_del_result_set
        res = self.vrs.normalise_outwards('NM_002111.8', 'CAG','', 202,205)
        assert res == self.NM_002111_8_del_result_set
        res = self.vrs.normalise_outwards('NM_002111.8', 'CAGCA','CA', 202,207)
        assert res == self.NM_002111_8_del_result_set

    def test_vrs_normalise_outwards_ex_done_del(self):
        # expand already done del
        res = self.vrs.normalise_outwards(
                'NR_110010.2', 'GCGCGCG','GCGCG', 25,32)
        assert res == self.NR_110010_2_del_result_set_end
        res = self.vrs.normalise_outwards(
                'NM_003073.5', 'AGAAGAAGA','AGAAGA', 1288,1297)
        assert res == self.NM_003073_5_del_result_set
        res = self.vrs.normalise_outwards(
            'NM_002111.8',
            'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA',
            'CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA',
            196,261)
        assert res == self.NM_002111_8_del_result_set

    def test_vrs_normalise_outwards_ex_into_start_del(self):
        # expand into start del
        res = self.vrs.normalise_outwards('NM_001330096.1', 'CT','', 22,24)
        assert res == self.res_set_d_start_NM_001330096_1_end
        res = self.vrs.normalise_outwards('NM_001330096.1', 'TCT','T', 21,24)
        assert res == self.res_set_d_start_NM_001330096_1
        res = self.vrs.normalise_outwards(
                'NM_001330112.1', 'GGGGCGGGTC','', 10,20)
        assert res == self.res_set_d_start_NM_001330112_1
        res = self.vrs.normalise_outwards(
                'NM_001330112.1', 'TCGGGGCGGGTC', 'TC', 8,20)
        assert res == self.res_set_d_start_NM_001330112_1_end

    def test_vrs_normalise_outwards_ex_into_start_done_del(self):
        # expand start already done del
        res = self.vrs.normalise_outwards(
                'NM_001330096.1',
                'CTCTCTCTCTCTCTCTCTCTCTCTCT',
                'CTCTCTCTCTCTCTCTCTCTCTCT',
                0,26)
        assert res == self.res_set_d_start_NM_001330096_1_end
        res = self.vrs.normalise_outwards(
                'NM_001330112.1',
                'GGGGCGGGTCGGGGCGGGTCGGGGCGGG',
                'GGGGCGGGTCGGGGCGGG',
                0,28)
        assert res == self.res_set_d_start_NM_001330112_1

    def test_vrs_normalise_outwards_ex_into_end_del(self):
        # expand into end del
        res = self.vrs.normalise_outwards('NM_052928.3', 'GTTTT','', 4380,4385)
        assert res == self.res_set_d_end_NM_052928_3
        res = self.vrs.normalise_outwards(
                'NM_052928.3', 'TTGTTTT','TT', 4378,4385)
        assert res == self.res_set_d_end_NM_052928_3_end
        res = self.vrs.normalise_outwards('NM_002129.4', 'TA','', 1422,1424)
        assert res == self.res_set_d_end_NM_002129_4_end
        res = self.vrs.normalise_outwards('NM_002129.4', 'TATA','TA', 1422,1426)
        assert res == self.res_set_d_end_NM_002129_4_end

    def test_vrs_normalise_outwards_ex_into_end_done_del(self):
        # expand end already done del
        res = self.vrs.normalise_outwards(
                'NM_052928.3',
                'GTTTTGTTTTGTTTTGTTTTGTTTTGTTTT',
                'GTTTTGTTTTGTTTTGTTTTGTTTT',
                4375,4405)
        assert res == self.res_set_d_end_NM_052928_3
        res = self.vrs.normalise_outwards(
                'NM_002129.4',
                'ATATATATATATATATATAT',
                'ATATATATATATATATAT',
                1421,1441)
        assert res == self.res_set_d_end_NM_002129_4

class TestVRSOutputVRSObjectIDs(unittest.TestCase):
    """
    Tests for the internal vrs_util.py module, particularly the HGVS_to_VRS
    object functions required for fetching and making VRS IDs.

    A good number of these tests take their required output from the VRS docs
    (current as of VRS 2.1).
    tests:
        Test VRS IDs for all data types with specified id's in the VRS docs,
            using stored VRS data as input
        Test that the VRS id generation code produces the expected output for
            the given input.
    """
    @classmethod
    def setup_class(self):
        self.vrs = HGVS_to_VRS()

    def test_vrs_seq_ids_bad_id(self):
        # test 1 bad ID to check that it asserts in some way
        with self.assertRaises(KeyError):
            self.vrs.id_fetch('LIVE_PARROT')

    def test_vrs_seq_ids(self):
        """
        Test the basic fetching of VRS sequence IDs using given hgvs ID.
        The VRS id_fetch() calls get_vrs_id_for_seq() to get the base id
        via SeqRepo, this might change in the future, but we want to test the
        id_fetch in use for this, not the base function, as the id_fetch
        function in use is the target of these tests.
        """
        # the internal set of sequences used:
        # NR_110010.2, NM_003073.5,  NM_002111.8, NM_001350922.2,
        # NM_001330096.1, NM_001330112.1, NM_052928.3, and NM_002129.4
        vrs_seq_id = self.vrs.id_fetch('NR_110010.2')
        assert vrs_seq_id == 'SQ.Ksaa229gzGC4Uc3ytCixai4vziud-MKi'
        vrs_seq_id = self.vrs.id_fetch('NM_003073.5')
        assert vrs_seq_id == 'SQ.pL2lXe9Qjg9ME_5vP0o85MgTUIqTAklc'
        vrs_seq_id = self.vrs.id_fetch('NM_002111.8')
        assert vrs_seq_id == 'SQ.4qBhnA470l_6xfLWJyi8WkOKeW6u5KJJ'
        vrs_seq_id = self.vrs.id_fetch('NM_001350922.2')
        assert vrs_seq_id == 'SQ.OgH_wLQV1lYs_Z107fts3xI9CuarHG7_'
        vrs_seq_id = self.vrs.id_fetch('NM_001330096.1')
        assert vrs_seq_id == 'SQ.Xzqbjy7R3hDF81IInpgNugs7N3c0qk2M'
        vrs_seq_id = self.vrs.id_fetch('NM_001330112.1')
        assert vrs_seq_id == 'SQ.L_1dwsmKTbVRGrbiF4p6FrgaWkWKWycu'
        vrs_seq_id = self.vrs.id_fetch('NM_052928.3')
        assert vrs_seq_id == 'SQ.pGqCl6h3ASKqxNA6o0aLmx_6NIoQio9l'
        vrs_seq_id = self.vrs.id_fetch('NM_002129.4')
        assert vrs_seq_id == 'SQ.8rOzThqGEimTT1vXu8KEb73YqgPBMfWJ'
        # sequence IDs found in the vrs documentation
        # allele docs
        vrs_seq_id = self.vrs.id_fetch('NC_000001.11')
        assert vrs_seq_id == 'SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO'
        # sequence docs and from:
        # https://github.com/ga4gh/vrs/blob/2.0/validation/models.yaml
        vrs_seq_id = self.vrs.id_fetch('NC_000007.14')
        assert vrs_seq_id == "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"

    def test_vrs_flatten(self):
        """ Test the VRS flattening logic for basic functionality and currently
        unused but, in the future needed, behaviour """
        # Test flattening start stop spans
        flat = self.vrs.vrs_flatten({
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"
                },
            "start": [44908820,44908821],
            "end": [44908822,44908823]
            })
        assert flat == \
            '{"end":"[44908822,44908823]","sequenceReference":{'\
            '"refgetAccession":"SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","type":"'\
            'SequenceReference"},"start":"[44908820,44908821]","type":'\
            '"SequenceLocation"}'
        flat = self.vrs.vrs_flatten({
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"
                },
            "start": [None,44908821],
            "end": [44908822,None]
            })
        assert flat == \
            '{"end":"[44908822,null]","sequenceReference":{"refgetAccession":'\
            '"SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","type":"SequenceReference"}'\
            ',"start":"[null,44908821]","type":"SequenceLocation"}'
        # Test specified digest_keys
        flat = self.vrs.vrs_flatten({
            "type": "SequenceReference",
            "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            "test":"bla"
            }, dict_keys=["refgetAccession","test","type"])
        assert flat == \
            '{"refgetAccession":"SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","test":'\
            '"bla","type":"SequenceReference"}'


    def test_vrs_object_ids_assert_invalid_type(self):
        # test type assert on flatten & hash on non valid data type
        with self.assertRaises(TypeError):
            vrs_eg_seq_loc = {
                    "type": "SequenceReference",
                    "refgetAccession": self.vrs
                    }
            self.vrs.vrs_flatten(vrs_eg_seq_loc,digest=True)

    def test_vrs_object_ids_assert_fail_no_vrs_id_prefix(self):
        # test value error on VRS type without VRS id prefix
        with self.assertRaises(ValueError):
            vrs_state = {
                    "type": "ReferenceLengthExpression",
                    "length": 11,
                    "repeatSubunitLength": 3
                    }
            self.vrs.vrs_flatten(vrs_state,digest=True)

    # now test working ids
    def test_vrs_object_ids_seq_loc(self):
        vrs_eg_seq_loc = {
                "id": "ga4gh:SL.4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT",
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"
                    },
                "start": 44908821,
                "end": 44908822
                }
        vrs_id = self.vrs.vrs_flatten(vrs_eg_seq_loc,digest=True)
        assert vrs_id == vrs_eg_seq_loc['id']

    def test_vrs_object_ids_full_allele(self):
        vrs_eg_allele = {
            "id": "ga4gh:VA.Oop4kjdTtKcg1kiZjIJAAR3bp7qi4aNT",
            "type": "Allele",
            "expressions": [
                {
                    "syntax": "spdi",
                    "value": "NC_000001.11:40819438:CTCCTCCT:CTCCTCCTCCT"
                    }
                ],
            "location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                    "residueAlphabet": "na",
                    "id": "NC_000001.11"
                    },
                "start": 40819438,
                "end": 40819446
                },
            "state": {
                "type": "ReferenceLengthExpression",
                "length": 11,
                "repeatSubunitLength": 3
                }
            }
        vrs_id = self.vrs.vrs_flatten(vrs_eg_allele,digest=True)
        assert vrs_id == vrs_eg_allele['id']

    def test_vrs_object_ids_allele_from_eg_page(self):
        vrs_eg_allele_from_id_eg_page = {
          "location": {
              "type": "SequenceLocation",
              "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
              },
              "start": 44908821,
              "end": 44908822
          },
          "state": {
            "type": "LiteralSequenceExpression",
            "sequence": "T"
          },
          "type": "Allele"
        }
        #{"location":"wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz","state":
        # {"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}
        # ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt_
        # note that the _ on the end is not genarated even when running the
        # sample code given by
        # https://vrs.ga4gh.org/en/latest/conventions/computed_identifiers.html#digest-serialization
        # *exactly* as-is, the rest of the ID is identical
        vrs_id = self.vrs.vrs_flatten(vrs_eg_allele_from_id_eg_page,digest=True)
        assert vrs_id == 'ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt'

    def test_vrs_object_ids_b_moddels_notebook(self):
        # test the ids from the basic models notebook
        # https://github.com/ga4gh/vrs-python/blob/main/notebooks/getting_started/3_Basic_Models.ipynb
        #'NM_002439.5'# SQ.Pw3Ch0x3XWD6ljsnIfmk_NERcZCI9sNM
        vrs_eg_allele_from_basic_mod_nb = {
         'id': 'ga4gh:VA.5C67OBmCLuHPgDkCQj7EOMih58BS2Eor',
         'type': 'Allele',
         'location': {'type': 'SequenceLocation',
          'sequenceReference': {'type': 'SequenceReference',
            'refgetAccession': 'SQ.Pw3Ch0x3XWD6ljsnIfmk_NERcZCI9sNM'},
          'start': 80656509,
          'end': 80656510},
         'state': {'type': 'LiteralSequenceExpression', 'sequence': 'TT'}}
        vrs_id_loc = self.vrs.vrs_flatten(
                vrs_eg_allele_from_basic_mod_nb['location'],digest=True)
        assert vrs_id_loc == 'ga4gh:SL.lGxOP1JRd4dysmrOVaskO5P_35DyCLnx'
        vrs_id = self.vrs.vrs_flatten(
                vrs_eg_allele_from_basic_mod_nb,digest=True)
        assert vrs_id == 'ga4gh:VA.5C67OBmCLuHPgDkCQj7EOMih58BS2Eor'

    def test_vrs_object_ids_2_0_validation_models1(self):
        # test id creation for Alleles and sequence loacations found in
        # https://github.com/ga4gh/vrs/blob/2.0/validation/models.yaml
        vrs_allele_valid_model_test_set_lit= {'type': 'Allele',
         'location': {'type': 'SequenceLocation',
          'sequenceReference': {'type': 'SequenceReference',
            'refgetAccession': 'SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl'},
          'start': 44908821,
          'end': 44908822},
         'state': {'type': 'LiteralSequenceExpression', 'sequence': 'T'}}
        vrs_id = self.vrs.vrs_flatten(
                vrs_allele_valid_model_test_set_lit['location'],digest=True)
        assert vrs_id == 'ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz'
        vrs_id = self.vrs.vrs_flatten(
                vrs_allele_valid_model_test_set_lit,digest=True)
        assert vrs_id == 'ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt'

    def test_vrs_object_ids_2_0_validation_models2(self):
        # test id creation for Alleles and sequence loacations found in
        # https://github.com/ga4gh/vrs/blob/2.0/validation/models.yaml
        vrs_allele_valid_model_test_set_len= {'type': 'Allele',
         'location': {'type': 'SequenceLocation',
          'sequenceReference': {'type': 'SequenceReference',
            'refgetAccession': 'SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO'},
          'start': 40819438,
          'end': 40819446},
         'state': {
          'type': 'ReferenceLengthExpression',
          'length': 11,
          'repeatSubunitLength': 3}}
        vrs_id = self.vrs.vrs_flatten(
                vrs_allele_valid_model_test_set_len['location'],digest=True)
        assert vrs_id == 'ga4gh:SL.nQGBuvRQOLEboA5TYtcz975fp_GulxbZ'
        vrs_id = self.vrs.vrs_flatten(
                vrs_allele_valid_model_test_set_len,digest=True)
        assert  vrs_id == 'ga4gh:VA.Oop4kjdTtKcg1kiZjIJAAR3bp7qi4aNT'

    def test_vrs_object_ids_notebook_ex_allele_translator1(self):
        # https://github.com/ga4gh/vrs-python/blob/main/notebooks/getting_started/4_Exploring_the_AlleleTranslator.ipynb
        hgvs_out_lit = {'type': 'Allele',
         'location': {'type': 'SequenceLocation',
          'sequenceReference': {'type': 'SequenceReference',
           'refgetAccession': 'SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI'},
          'start': 80656488,
          'end': 80656489},
         'state': {'type': 'LiteralSequenceExpression', 'sequence': 'T'}}
        vrs_id = self.vrs.vrs_flatten(hgvs_out_lit['location'],digest=True)
        assert vrs_id == 'ga4gh:SL.JiLRuuyS5wefF_6-Vw7m3Yoqqb2YFkss'
        vrs_id = self.vrs.vrs_flatten(hgvs_out_lit,digest=True)
        assert vrs_id == 'ga4gh:VA.ebezGL6HoAhtGJyVnB_mE5BH18ntKev4'

    def test_vrs_object_ids_notebook_ex_allele_translator2(self):
        # https://github.com/ga4gh/vrs-python/blob/main/notebooks/getting_started/4_Exploring_the_AlleleTranslator.ipynb
        hgvs_out_lit2 = {'type': 'Allele',
         'location': {'type': 'SequenceLocation',
          'sequenceReference': {'type': 'SequenceReference',
           'refgetAccession': 'SQ.vbjOdMfHJvTjK_nqvFvpaSKhZillW0SX'},
          'start': 79952307,
          'end': 79952308},
         'state': {'type': 'LiteralSequenceExpression', 'sequence': 'T'}}
        vrs_id = self.vrs.vrs_flatten(hgvs_out_lit2['location'],digest=True)
        assert vrs_id == 'ga4gh:SL.Y-itBtqe9IwbxyL4EVZ4T_X9TUsdbJ22'
        vrs_id = self.vrs.vrs_flatten(hgvs_out_lit2,digest=True)
        assert vrs_id == 'ga4gh:VA.hEyB1sGiQrdrPFIq4u4CF17uAuUs2Wvx'

class TestVRSOutputVRSObjects(unittest.TestCase):
    """
    Tests for the internal vrs_util.py module, particularly the HGVS_to_VRS
    object function hgvs_single_var_to_vrs, which is the the basic core
    hgvs->vrs function.

    Starts with standard input as from
    https://github.com/ga4gh/vrs/blob/2.0/validation/models.yaml
    then moves to extra tests on shortcuts and error/null return states.

    """
    @classmethod
    def setup_class(self):
        self.vrs = HGVS_to_VRS()

    def test_vrs_object_output_NC_000019_10_g_44908822C_T_set(self):

        # allele test
        # hgvs in NC_000019.10:g44908822_44908821C>T
        hgvs = vv.hp.parse('NC_000019.10:g.44908822C>T')
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(hgvs)
        assert hgvs_allele['type'] == 'Allele'
        assert hgvs_allele['state']['type'] == 'LiteralSequenceExpression'
        assert hgvs_allele['location']['type'] == 'SequenceLocation'
        assert hgvs_allele['location']['sequenceReference'][
                'type'] == 'SequenceReference'
        assert hgvs_allele['location']['end'] == 44908822
        assert hgvs_allele['location']['start'] == 44908821
        assert hgvs_allele['location']['sequenceReference'][
                'refgetAccession'] == 'SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl'
        assert hgvs_allele['location'][
                'id'] == 'ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz'
        assert hgvs_allele['state']['sequence'] == 'T'
        assert hgvs_allele['id'] == 'ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt'

    def test_vrs_object_output_NC_000001_11_g_40819439_40819446delCTCCTCCTinsCTCCTCCTCCT(self):
        # No need to test the 'canonical' JSON serilised form as this, is not available as output,
        # and should be redundant with the checksum ID test
        # hgvs in  "NC_000001.11:40819439_40819446delCTCCTCCTinsCTCCTCCTCCT
        hgvs = vv.hp.parse(
            "NC_000001.11:g.40819439_40819446delCTCCTCCTinsCTCCTCCTCCT")
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(hgvs)
        assert hgvs_allele['type'] == 'Allele'
        assert hgvs_allele['state']['type'] == 'ReferenceLengthExpression'
        assert hgvs_allele['location']['type'] == 'SequenceLocation'
        assert hgvs_allele['location']['sequenceReference'][
                'type'] == 'SequenceReference'
        assert hgvs_allele['location']['end'] == 40819446
        assert hgvs_allele['location']['start'] == 40819438
        assert hgvs_allele['location']['sequenceReference'][
                'refgetAccession'] == 'SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO'
        assert hgvs_allele['location'][
                'id'] == 'ga4gh:SL.nQGBuvRQOLEboA5TYtcz975fp_GulxbZ'
        assert hgvs_allele['state']['length'] == 11
        assert hgvs_allele['state']['repeatSubunitLength'] == 3
        assert hgvs_allele['id'] == 'ga4gh:VA.Oop4kjdTtKcg1kiZjIJAAR3bp7qi4aNT'

    def test_vrs_object_output_NC_000005_10_g_80656489C_T(self):
        # test for good output from hgvs input matching
        # https://github.com/ga4gh/vrs-python/blob/main/notebooks/getting_started/4_Exploring_the_AlleleTranslator.ipynb
        #From HGVS
        hgvs = vv.hp.parse("NC_000005.10:g.80656489C>T")
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(hgvs)
        assert hgvs_allele['type'] == 'Allele'
        assert hgvs_allele['state']['type'] == 'LiteralSequenceExpression'
        assert hgvs_allele['location']['type'] == 'SequenceLocation'
        assert hgvs_allele['location']['sequenceReference'][
                'type'] == 'SequenceReference'
        assert hgvs_allele['location']['end'] == 80656489
        assert hgvs_allele['location']['start'] == 80656488
        assert hgvs_allele['location']['sequenceReference'][
                'refgetAccession'] == 'SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI'
        assert hgvs_allele['location'][
                'id'] == 'ga4gh:SL.JiLRuuyS5wefF_6-Vw7m3Yoqqb2YFkss'
        assert hgvs_allele['state']['sequence'] == 'T'
        assert hgvs_allele['id'] == 'ga4gh:VA.ebezGL6HoAhtGJyVnB_mE5BH18ntKev4'

    def test_vrs_object_output_NC_000005_9_g_79952308C_T(self):
        hgvs = vv.hp.parse("NC_000005.9:g.79952308C>T")
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(hgvs)
        assert hgvs_allele['type'] == 'Allele'
        assert hgvs_allele['state']['type'] == 'LiteralSequenceExpression'
        assert hgvs_allele['location']['type'] == 'SequenceLocation'
        assert hgvs_allele['location']['sequenceReference'][
                'type'] == 'SequenceReference'
        assert hgvs_allele['location']['end'] == 79952308
        assert hgvs_allele['location']['start'] == 79952307
        assert hgvs_allele['location']['sequenceReference'][
                'refgetAccession'] == 'SQ.vbjOdMfHJvTjK_nqvFvpaSKhZillW0SX'
        assert hgvs_allele['location'][
                'id'] == 'ga4gh:SL.Y-itBtqe9IwbxyL4EVZ4T_X9TUsdbJ22'
        assert hgvs_allele['state']['sequence'] == 'T'
        assert hgvs_allele['id'] == 'ga4gh:VA.hEyB1sGiQrdrPFIq4u4CF17uAuUs2Wvx'


    def test_vrs_obj_out_norm_identity_diff_input_start(self):
        """
        test VRS IDs are the same for multiple hgvs inputs that should
        normalise to the same output VRS start hgvs input.
        """
        # hgvs in  "NC_000001.11:40819439_40819446delCTCCTCCTinsCTCCTCCTCCT
        hgvs = vv.hp.parse("NC_000001.11:g.40819438_40819439insCTC")
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(hgvs)
        assert hgvs_allele['type'] == 'Allele'
        assert hgvs_allele['state']['type'] == 'ReferenceLengthExpression'
        assert hgvs_allele['location']['type'] == 'SequenceLocation'
        assert hgvs_allele['location']['sequenceReference'][
                'type'] == 'SequenceReference'
        assert hgvs_allele['location']['end'] == 40819446
        assert hgvs_allele['location']['start'] == 40819438
        assert hgvs_allele['location']['sequenceReference'][
                'refgetAccession'] == 'SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO'
        assert hgvs_allele['location'][
                'id'] == 'ga4gh:SL.nQGBuvRQOLEboA5TYtcz975fp_GulxbZ'
        assert hgvs_allele['state']['length'] == 11
        assert hgvs_allele['state']['repeatSubunitLength'] == 3
        assert hgvs_allele['id'] == 'ga4gh:VA.Oop4kjdTtKcg1kiZjIJAAR3bp7qi4aNT'

    def test_vrs_obj_out_norm_identity_diff_input_end(self):
        """
        test VRS IDs are the same for multiple hgvs inputs that should
        normalise to the same output VRS end hgvs input.
        """
        # hgvs in  "NC_000001.11:40819439_40819446delCTCCTCCTinsCTCCTCCTCCT
        hgvs = vv.hp.parse("NC_000001.11:g.40819446_40819447insCCT")
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(hgvs)
        assert hgvs_allele['type'] == 'Allele'
        assert hgvs_allele['state']['type'] == 'ReferenceLengthExpression'
        assert hgvs_allele['location']['type'] == 'SequenceLocation'
        assert hgvs_allele['location']['sequenceReference'][
                'type'] == 'SequenceReference'
        assert hgvs_allele['location']['end'] == 40819446
        assert hgvs_allele['location']['start'] == 40819438
        assert hgvs_allele['location']['sequenceReference'][
                'refgetAccession'] == 'SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO'
        assert hgvs_allele['location'][
                'id'] == 'ga4gh:SL.nQGBuvRQOLEboA5TYtcz975fp_GulxbZ'
        assert hgvs_allele['state']['length'] == 11
        assert hgvs_allele['state']['repeatSubunitLength'] == 3
        assert hgvs_allele['id'] == 'ga4gh:VA.Oop4kjdTtKcg1kiZjIJAAR3bp7qi4aNT'

    def test_vrs_obj_out_norm_identity_diff_input_mid(self):
        """
        test VRS IDs are the same for multiple hgvs inputs that should
        normalise to the same output VRS middle hgvs input.
        """
        # hgvs in  "NC_000001.11:40819439_40819446delCTCCTCCTinsCTCCTCCTCCT
        hgvs = vv.hp.parse("NC_000001.11:g.40819443_40819444insCCT")
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(hgvs)
        assert hgvs_allele['type'] == 'Allele'
        assert hgvs_allele['state']['type'] == 'ReferenceLengthExpression'
        assert hgvs_allele['location']['type'] == 'SequenceLocation'
        assert hgvs_allele['location']['sequenceReference'][
                'type'] == 'SequenceReference'
        assert hgvs_allele['location']['end'] == 40819446
        assert hgvs_allele['location']['start'] == 40819438
        assert hgvs_allele['location']['sequenceReference'][
                'refgetAccession'] == 'SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO'
        assert hgvs_allele['location'][
                'id'] == 'ga4gh:SL.nQGBuvRQOLEboA5TYtcz975fp_GulxbZ'
        assert hgvs_allele['state']['length'] == 11
        assert hgvs_allele['state']['repeatSubunitLength'] == 3
        assert hgvs_allele['id'] == 'ga4gh:VA.Oop4kjdTtKcg1kiZjIJAAR3bp7qi4aNT'

    def test_vrs_obj_out_norm_identity_diff_set_dup(self):
        """
        test VRS IDs are the same for multiple hgvs inputs that should
        normalise to the same output VRS middle hgvs input.
        Dup, like ins also has a separate code path from del/delins/subtitue
        test this as well.
        """
        # hgvs in  "NC_000001.11:40819439_40819446delCTCCTCCTinsCTCCTCCTCCT
        hgvs = vv.hp.parse("NC_000001.11:g.40819442_40819444dup")
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(hgvs)
        assert hgvs_allele['type'] == 'Allele'
        assert hgvs_allele['state']['type'] == 'ReferenceLengthExpression'
        assert hgvs_allele['location']['type'] == 'SequenceLocation'
        assert hgvs_allele['location']['sequenceReference'][
                'type'] == 'SequenceReference'
        assert hgvs_allele['location']['end'] == 40819446
        assert hgvs_allele['location']['start'] == 40819438
        assert hgvs_allele['location']['sequenceReference'][
                'refgetAccession'] == 'SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO'
        assert hgvs_allele['location'][
                'id'] == 'ga4gh:SL.nQGBuvRQOLEboA5TYtcz975fp_GulxbZ'
        assert hgvs_allele['state']['length'] == 11
        assert hgvs_allele['state']['repeatSubunitLength'] == 3
        assert hgvs_allele['id'] == 'ga4gh:VA.Oop4kjdTtKcg1kiZjIJAAR3bp7qi4aNT'

    def test_vrs_object_output_part_ambig_pos(self):
        """
        test variant with partially ambiguous positions
        eg ACCA to ACGCCA could be an ins either before (of CG) or after
        (of GC) the fist C
        so VRS expands to cover the first C
        hgvs inputs should be normalised 3' in this case, but this also should
        work either way
        done for TCCT >TCGCCT and TCCT >TCCGCT in fully expanded delins + both
        directions
        """
        hgvs = vv.hp.parse("NC_000001.11:g.40819444_40819445delCCinsCCGC")
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(hgvs)
        assert hgvs_allele['type'] == 'Allele'
        assert hgvs_allele['state']['type'] == 'LiteralSequenceExpression'
        assert hgvs_allele['location']['type'] == 'SequenceLocation'
        assert hgvs_allele['location']['sequenceReference'][
                'type'] == 'SequenceReference'
        assert hgvs_allele['location']['end'] == 40819445
        assert hgvs_allele['location']['start'] == 40819444
        assert hgvs_allele['location']['sequenceReference'][
                'refgetAccession'] == 'SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO'
        assert hgvs_allele['location'][
                'id'] == 'ga4gh:SL.UhfMqgh_jnDrVNWaFx-_ksYssV1-9Hr3'
        assert hgvs_allele['state']['sequence'] == 'CGC'
        assert hgvs_allele['state']['type'] == 'LiteralSequenceExpression'
        assert hgvs_allele['id'] == 'ga4gh:VA.r-IMwTjeSjYhwqkaIX7tloTeeMnuEGYy'

        hgvs = vv.hp.parse("NC_000001.11:g.40819445delCinsCGC")
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(hgvs)
        assert hgvs_allele['type'] == 'Allele'
        assert hgvs_allele['state']['type'] == 'LiteralSequenceExpression'
        assert hgvs_allele['location']['type'] == 'SequenceLocation'
        assert hgvs_allele['location']['sequenceReference'][
                'type'] == 'SequenceReference'
        assert hgvs_allele['location']['end'] == 40819445
        assert hgvs_allele['location']['start'] == 40819444
        assert hgvs_allele['location']['sequenceReference'][
                'refgetAccession'] == 'SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO'
        assert hgvs_allele['location'][
                'id'] == 'ga4gh:SL.UhfMqgh_jnDrVNWaFx-_ksYssV1-9Hr3'
        assert hgvs_allele['state']['sequence'] == 'CGC'
        assert hgvs_allele['state']['type'] == 'LiteralSequenceExpression'
        assert hgvs_allele['id'] == 'ga4gh:VA.r-IMwTjeSjYhwqkaIX7tloTeeMnuEGYy'


        # extra tests VRS varified, outside of the examples?

    def test_vrs_object_output_eq_shortcut(self):
        "hgvs_single_var_to_vrs on == variant to test internal shortcut"
        hgvs = vv.hp.parse("NC_000001.11:g.40819444_40819445=")
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(hgvs)
        assert hgvs_allele['type'] == 'Allele'
        assert hgvs_allele['state']['type'] == 'ReferenceLengthExpression'
        assert hgvs_allele['location']['type'] == 'SequenceLocation'
        assert hgvs_allele['location']['sequenceReference'][
                'type'] == 'SequenceReference'
        assert hgvs_allele['location']['end'] == 40819445
        assert hgvs_allele['location']['start'] == 40819443
        assert hgvs_allele['location']['sequenceReference'][
                'refgetAccession'] == 'SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO'
        assert hgvs_allele['location'][
                'id'] == 'ga4gh:SL.zxbyl3-hn1MHmcZYF0N8WWaAzNZcryGk'
        assert hgvs_allele['state']['length'] == 2
        assert hgvs_allele['state']['repeatSubunitLength'] == 2
        assert hgvs_allele['state']['type'] == 'ReferenceLengthExpression'

    def test_vrs_object_output_fail_with_no_mapper_for_c_input(self):
        # Test C type works, and failes propperly with mapper missing
        with self.assertRaises(ValueError):
            self.vrs.hgvs_single_var_to_vrs(
                vv.hp.parse('NM_015120.4:c.35T>C'))

    def test_vrs_object_output_for_c_input(self):
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(
                vv.hp.parse('NM_015120.4:c.35T>C'),
                variant = MockVVData())
        assert hgvs_allele == {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.mq1Saz29534GzZ6T_jCfUjsqNlifC7pZ',
                    'label': 'NM_015120.4'
                    },
                'start': 145,
                'end': 146,
                'id': 'ga4gh:SL.x_rdZXUcg7V8Ae5pPVv0OS-faxKyzNQb'
                },
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'C'
                },
            'id': 'ga4gh:VA.x8gKEBS-s638D8uuuQXVJKIgCDcuQJCW'
            }

    def test_vrs_object_output_no_norm_shortcut_ins(self):
        # test that ins and del that don't at all normalise outwards
        # use the shortcut correctly
        hgvs = vv.hp.parse("NC_000001.11:g.40819444_40819445insG")
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(hgvs)
        assert hgvs_allele['type'] == 'Allele'
        assert hgvs_allele['state']['type'] == 'LiteralSequenceExpression'
        assert hgvs_allele['location']['type'] == 'SequenceLocation'
        assert hgvs_allele['location']['sequenceReference'][
                'type'] == 'SequenceReference'
        assert hgvs_allele['location']['end'] == 40819444
        assert hgvs_allele['location']['start'] == 40819444
        assert hgvs_allele['location']['sequenceReference'][
                'refgetAccession'] == 'SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO'
        assert hgvs_allele['location'][
                'id'] == 'ga4gh:SL.aBJ1_py6od8pI4OWHQegaTxdAsszkbBy'
        assert hgvs_allele['state']['sequence'] == 'G'
        assert hgvs_allele['state']['type'] == 'LiteralSequenceExpression'
        assert hgvs_allele['id'] == 'ga4gh:VA.i_1aVFr-IQvl8bXoIX_UjwAU6_IwAmCD'

    def test_vrs_object_output_no_norm_shortcut_del(self):
        hgvs = vv.hp.parse("NC_000001.11:g.40819444_40819445delCC")
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(hgvs)
        assert hgvs_allele['type'] == 'Allele'
        assert hgvs_allele['state']['type'] == 'ReferenceLengthExpression'
        assert hgvs_allele['location']['type'] == 'SequenceLocation'
        assert hgvs_allele['location']['sequenceReference'][
                'type'] == 'SequenceReference'
        assert hgvs_allele['location']['end'] == 40819445
        assert hgvs_allele['location']['start'] == 40819443
        assert hgvs_allele['location']['sequenceReference'][
                'refgetAccession'] == 'SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO'
        assert hgvs_allele['location'][
                'id'] == 'ga4gh:SL.zxbyl3-hn1MHmcZYF0N8WWaAzNZcryGk'
        assert hgvs_allele['state']['repeatSubunitLength'] == 2
        assert hgvs_allele['state']['length'] == 0
        assert hgvs_allele['state']['type'] == 'ReferenceLengthExpression'
        assert hgvs_allele['id'] == 'ga4gh:VA.gpJtRDCma3r-7snQ90wFU7i-hRPuJKGh'

    def test_vrs_object_output_intronic_null_result(self):
        # test null result on intronic input
        hgvs_allele = self.vrs.hgvs_single_var_to_vrs(
                # not valid but works for test
                vv.hp.parse('NM_015120.4:c.35-6T>C'),
                variant = MockVVData())
        assert hgvs_allele is None

    def test_vrs_object_output_cache(self):
        # test cached = type
        vrs_cached = copy.copy(self.vrs)
        vrs_cached.cache = {}
        hgvs = vv.hp.parse("NC_000001.11:g.40819444_40819445=")
        hgvs_allele = vrs_cached.hgvs_single_var_to_vrs(hgvs)
        assert hgvs_allele['type'] == 'Allele'
        assert hgvs_allele['state']['type'] == 'ReferenceLengthExpression'
        assert hgvs_allele['location']['type'] == 'SequenceLocation'
        assert hgvs_allele['location']['sequenceReference'][
                'type'] == 'SequenceReference'
        assert hgvs_allele['location']['end'] == 40819445
        assert hgvs_allele['location']['start'] == 40819443
        assert hgvs_allele['location']['sequenceReference'][
                'refgetAccession'] == 'SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO'
        assert hgvs_allele['location'][
                'id'] == 'ga4gh:SL.zxbyl3-hn1MHmcZYF0N8WWaAzNZcryGk'
        assert hgvs_allele['state']['length'] == 2
        assert hgvs_allele['state']['repeatSubunitLength'] == 2
        assert hgvs_allele['state']['type'] == 'ReferenceLengthExpression'
        assert vrs_cached.cache["NC_000001.11:g.40819444_40819445="] == \
                hgvs_allele

class TestVRSOutputFromDummyVV(unittest.TestCase):
    """
    Test the generation of VRS output from dummy VV result objects, unlike raw
    VV this relies objects this does not rely on VV working consistently with
    previous results. Unfortunately this also means that this test set may need
    to be changed in the case of VV alterations/improvements.

    See also the direct VV using tests below in the case of failures.
    """

    @classmethod
    def setup_class(self):
        self.vrs = HGVS_to_VRS()
        # core set of dummy input
        self.dummy_vv_output = MockVVData(
                selected_assembly='GRCh38',
                original = 'test_input',
                warnings=[],
                lovd_messages = [],
                lovd_corrections = [],
                stable_gene_ids = {
                    'hgnc_id': 'HGNC:TST',
                    'entrez_gene_id': '0000',
                    'ucsc_id': 'uc0tst.1',
                    'omim_id': ['0000000']},
                gene_symbol = 'Dummy_S',
                description = "Text transcript description",)

    def test_vv_output_to_vrs_output_val_err(self):
        # test set 1: test that we trap and return variants with corrections
        # without VRS data
        # 1a test variant.warnings == ['Validation error'] response
        err_fail = copy.copy(self.dummy_vv_output)
        err_fail.warnings=['Validation error']
        output = self.vrs.variant_validator_output_set_to_vrs(err_fail)
        assert output['selected_assembly'] == err_fail.selected_assembly
        assert output['submitted_variant'] == err_fail.original
        assert output['warnings_and_messages'] == {
                'validation_warnings': ['Validation error'],
                'lovd_messages': [],
                'lovd_corrections': []
                }
        assert output['gene_ids']['hgnc_id'] == 'HGNC:TST'
        assert output['gene_ids']['entrez_gene_id'] == '0000'
        assert output['gene_ids']['ucsc_id'] == 'uc0tst.1'
        assert output['gene_ids']['omim_id'] == ['0000000']
        assert output['gene_ids']['current_symbol'] == 'Dummy_S'
        assert output['transcript_description'] == err_fail.description

    def test_vv_output_to_vrs_output_warn_type(self):
        # 1b test variant.output_type_flag = 'warning'
        err_warn = copy.copy(self.dummy_vv_output)
        err_warn.output_type_flag = 'warning'
        output = self.vrs.variant_validator_output_set_to_vrs(err_warn)
        assert output['selected_assembly'] == err_warn.selected_assembly
        assert output['submitted_variant'] == err_warn.original
        assert output['warnings_and_messages'] == {
                'validation_warnings': [],
                'lovd_messages': [],
                'lovd_corrections': []
                }
        assert output['gene_ids']['hgnc_id'] == 'HGNC:TST'
        assert output['gene_ids']['entrez_gene_id'] == '0000'
        assert output['gene_ids']['ucsc_id'] == 'uc0tst.1'
        assert output['gene_ids']['omim_id'] == ['0000000']
        assert output['gene_ids']['current_symbol'] == 'Dummy_S'
        assert output['transcript_description'] ==  err_warn.description

    def test_vv_output_to_vrs_output_warn_for_working_type(self):
        # test warnings in output for variants with output_type_flag and w/o
        err_warn = copy.copy(self.dummy_vv_output)
        err_warn.output_type_flag = 'warning'
        err_warn.warnings=['VV_Warn']
        err_warn.lovd_messages = ['LO_Warn']
        err_warn.lovd_corrections = ['LO_corr']
        output = self.vrs.variant_validator_output_set_to_vrs(err_warn)
        assert output['selected_assembly'] == err_warn.selected_assembly
        assert output['submitted_variant'] == err_warn.original
        assert output['warnings_and_messages'] == {
                'validation_warnings': ['VV_Warn'],
                'lovd_messages': ['LO_Warn'],
                'lovd_corrections': ['LO_corr']
                }
        assert output['gene_ids']['hgnc_id'] == 'HGNC:TST'
        assert output['gene_ids']['entrez_gene_id'] == '0000'
        assert output['gene_ids']['ucsc_id'] == 'uc0tst.1'
        assert output['gene_ids']['omim_id'] == ['0000000']
        assert output['gene_ids']['current_symbol'] == 'Dummy_S'
        assert output['transcript_description'] == err_warn.description

    def test_vv_output_to_vrs_output_gene_symbol(self):
        # test that gene symbols get handled regardless of stable gene id state
        no_gene_id = copy.copy(self.dummy_vv_output)
        no_gene_id.stable_gene_ids = None
        output = self.vrs.variant_validator_output_set_to_vrs(no_gene_id)
        assert output['gene_ids'] == {'current_symbol': 'Dummy_S'}
        no_gene_id.gene_symbol = None
        output = self.vrs.variant_validator_output_set_to_vrs(no_gene_id)
        assert output['gene_ids'] is None

    # seq tests done using data from input test 1
    def test_vv_output_to_vrs_output_r_type(self):
        # test 2: test that we catch and return the simpler r type input
        # 2a test for working output to valid input
        rna_test = copy.copy(self.dummy_vv_output)
        rna_test.rna_data = {
                'rna_variant':vv.hp.parse('NM_015120.4:c.35T>C'),
                'translation_slr':vv.hp.parse('NP_055935.4:p.(L12P)'),
                'usage_warnings':"RNA warning"}
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(rna_test)
        assert vrs_out['vrs_transcript_variant'] == {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.mq1Saz29534GzZ6T_jCfUjsqNlifC7pZ',
                    'label': 'NM_015120.4'
                    },
                'start': 145,
                'end': 146,
                'id': 'ga4gh:SL.x_rdZXUcg7V8Ae5pPVv0OS-faxKyzNQb'
                },
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'C'
                },
            'id': 'ga4gh:VA.x8gKEBS-s638D8uuuQXVJKIgCDcuQJCW'
            }
        assert vrs_out['vrs_predicted_protein_consequence'] == {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.ggxH8oENu3WDprE6Z_Qg1HvOiT0d8UNa',
                    'label': 'NP_055935.4'
                    },
                'start': 11,
                'end': 12,
                'id': 'ga4gh:SL.qOpPs9dwgqK6MuH3CTQUDoRYbi9IWwSh'
                },
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'P'
                },
            'id': 'ga4gh:VA.cJPulDXE1G73sZS4dHaZX_27-c4pNjYo'
            }

    def test_vv_output_to_vrs_output_r_type_missmatch1(self):
        # 2b test that r type variant causes assertions on
        # hgvs_transcript_variant or  primary_assembly_loci being set
        # (we don't check object type)
        rna_test_bad_1 = copy.copy(self.dummy_vv_output)
        rna_test_bad_1.rna_data = {
                'rna_variant':vv.hp.parse('NM_015120.4:c.35T>C'),
                'translation_slr':vv.hp.parse('NP_055935.4:p.(L12P)'),
                'usage_warnings':"RNA warning"}
        rna_test_bad_1.hgvs_transcript_variant = 1
        with self.assertRaises(Exception):
            self.vrs.variant_validator_output_set_to_vrs(rna_test_bad_1)

    def test_vv_output_to_vrs_output_r_type_missmatch2(self):
        rna_test_bad_2 = copy.copy(self.dummy_vv_output)
        rna_test_bad_2.rna_data = {
                'rna_variant':vv.hp.parse('NM_015120.4:c.35T>C'),
                'translation_slr':vv.hp.parse('NP_055935.4:p.(L12P)'),
                'usage_warnings':"RNA warning"}
        rna_test_bad_2.primary_assembly_loci = 2
        with self.assertRaises(Exception):
            self.vrs.variant_validator_output_set_to_vrs(rna_test_bad_2)

    def test_vv_output_to_vrs_output_working_exonic(self):
        # test 3: test for working responses to c/n type input that is exonic
        c_test = copy.copy(self.dummy_vv_output)
        c_test.hgvs_transcript_variant = vv.hp.parse('NM_015120.4:c.35T>C')
        hgvs_hg37 = vv.hp.parse('NC_000002.11:g.73613031delinsCGGA')
        hgvs_hg38 = vv.hp.parse('NC_000002.12:g.73385903delinsCGGA')
        c_test.primary_assembly_loci = {
            'hg19': {
                'hgvs_genomic_description': hgvs_hg37,
                'vcf': {
                    'chr': 'chr2', 'pos': '73613031', 'ref': 'T', 'alt': 'CGGA'}
                },
            'hg38': {
                'hgvs_genomic_description': hgvs_hg37,
                'vcf': {
                    'chr': 'chr2', 'pos': '73385903', 'ref': 'T', 'alt': 'CGGA'}
                },
            'grch37' : {
                'hgvs_genomic_description': hgvs_hg37,
                'vcf': {
                    'chr': '2', 'pos': '73613031', 'ref': 'T', 'alt': 'CGGA'}
                },
            'grch38': {
                'hgvs_genomic_description': hgvs_hg38,
                'vcf': {
                    'chr': '2', 'pos': '73385903', 'ref': 'T', 'alt': 'CGGA'}
                }
            }
        c_test.hgvs_predicted_protein_consequence = {
                'prot':vv.hp.parse('NP_055935.4:p.(L12P)')}
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(c_test)
        assert vrs_out['vrs_transcript_variant'] == {
                'type': 'Allele',
                'location': {
                    'type': 'SequenceLocation',
                    'sequenceReference': {
                        'type': 'SequenceReference',
                        'refgetAccession': 'SQ.mq1Saz29534GzZ6T_jCfUjsqNlifC7pZ',
                        'label': 'NM_015120.4'
                        },
                    'start': 145,
                    'end': 146,
                    'id': 'ga4gh:SL.x_rdZXUcg7V8Ae5pPVv0OS-faxKyzNQb'
                    },
                'state': {
                    'type': 'LiteralSequenceExpression',
                    'sequence': 'C'
                    },
                'id': 'ga4gh:VA.x8gKEBS-s638D8uuuQXVJKIgCDcuQJCW'
                }
        assert vrs_out['vrs_predicted_protein_consequence'] == {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.ggxH8oENu3WDprE6Z_Qg1HvOiT0d8UNa',
                    'label': 'NP_055935.4'
                    },
                'start': 11,
                'end': 12,
                'id': 'ga4gh:SL.qOpPs9dwgqK6MuH3CTQUDoRYbi9IWwSh'
                },
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'P'
                },
            'id': 'ga4gh:VA.cJPulDXE1G73sZS4dHaZX_27-c4pNjYo'
            }
        assert vrs_out['VRS_mappings_for_primary_assemblies'] == {
         'grch37': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.9KdcA9ZpY1Cpvxvg8bMSLYDUpsX6GDLO',
            'label': 'NC_000002.11'
           },
           'start': 73613030,
           'end': 73613031,
           'id': 'ga4gh:SL.IzA6GY6mvDiJ1EzMbe1Cxv6jyhg4u8t7'
           },
          'state': {
           'type': 'LiteralSequenceExpression',
           'sequence': 'CGGA'
          },
          'id': 'ga4gh:VA.rKc_HxeylEKYQzYo_e_kY1lqYLF6Ch5m'
         },
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g',
            'label': 'NC_000002.12'
           },
           'start': 73385902,
           'end': 73385903,
           'id': 'ga4gh:SL.ndPBndyntJTTRHuU6tYpgbN6TUsRv94B'
           },
          'state': {
           'type': 'LiteralSequenceExpression',
           'sequence': 'CGGA'
           },
          'id': 'ga4gh:VA.eDiJCZTUpLf1UQHaUFFP1UaPsHfbOjyC'
         }
        }
    def test_vv_output_to_vrs_output_working_exonic_w_rsg(self):
        c_test = copy.copy(self.dummy_vv_output)
        c_test.hgvs_transcript_variant = vv.hp.parse('NM_015120.4:c.35T>C')
        c_test.hgvs_refseqgene_variant = vv.hp.parse('NG_011690.1:g.5146T>C')
        hgvs_hg37 = vv.hp.parse('NC_000002.11:g.73613031delinsCGGA')
        hgvs_hg38 = vv.hp.parse('NC_000002.12:g.73385903delinsCGGA')
        c_test.primary_assembly_loci = {
                'hg19': {'hgvs_genomic_description': hgvs_hg37},
                'hg38': {'hgvs_genomic_description': hgvs_hg37},
                'grch37' : {'hgvs_genomic_description': hgvs_hg37},
                'grch38': {'hgvs_genomic_description': hgvs_hg38}
                }
        c_test.hgvs_predicted_protein_consequence = {
                'prot':vv.hp.parse('NP_055935.4:p.(L12P)')}
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(c_test)
        assert vrs_out['vrs_transcript_variant'] == {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.mq1Saz29534GzZ6T_jCfUjsqNlifC7pZ',
                    'label': 'NM_015120.4'
                    },
                'start': 145,
                'end': 146,
                'id': 'ga4gh:SL.x_rdZXUcg7V8Ae5pPVv0OS-faxKyzNQb'
                },
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'C'
                },
            'id': 'ga4gh:VA.x8gKEBS-s638D8uuuQXVJKIgCDcuQJCW'
            }
        assert vrs_out['vrs_predicted_protein_consequence'] == {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.ggxH8oENu3WDprE6Z_Qg1HvOiT0d8UNa',
                    'label': 'NP_055935.4'
                    },
                'start': 11,
                'end': 12,
                'id': 'ga4gh:SL.qOpPs9dwgqK6MuH3CTQUDoRYbi9IWwSh'
                },
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'P'
                },
            'id': 'ga4gh:VA.cJPulDXE1G73sZS4dHaZX_27-c4pNjYo'
            }
        assert vrs_out['VRS_mappings_for_primary_assemblies'] == {
         'grch37': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.9KdcA9ZpY1Cpvxvg8bMSLYDUpsX6GDLO',
            'label': 'NC_000002.11'
            },
           'start': 73613030,
           'end': 73613031,
           'id': 'ga4gh:SL.IzA6GY6mvDiJ1EzMbe1Cxv6jyhg4u8t7'
           },
          'state': {
           'type': 'LiteralSequenceExpression',
           'sequence': 'CGGA'
           },
          'id': 'ga4gh:VA.rKc_HxeylEKYQzYo_e_kY1lqYLF6Ch5m'
         },
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g',
            'label': 'NC_000002.12'
           },
           'start': 73385902,
           'end': 73385903,
           'id': 'ga4gh:SL.ndPBndyntJTTRHuU6tYpgbN6TUsRv94B'
           },
          'state': {
           'type': 'LiteralSequenceExpression',
           'sequence': 'CGGA'
           },
          'id': 'ga4gh:VA.eDiJCZTUpLf1UQHaUFFP1UaPsHfbOjyC'
         }
        }
        assert vrs_out['vrs_refseqgene_variant']  == {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.r47Jw6neAxxIQgeyZLmdPrqFN4S6h7hy',
                    'label': 'NG_011690.1'},
                'start': 5145,
                'end': 5146,
                'id': 'ga4gh:SL.3CEnQ3JJo2jQbraJ0UHIO-lqBbR4xs15'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'C'},
            'id': 'ga4gh:VA.TjhuWhUG-PJapHphK3r0HWyU5b9CgGG_'}

    def test_vv_output_to_vrs_output_intronic_main_chr(self):
        # test 4: test for working responses to c/n type input that is intronic
        # 4a intronic matches main genomic ref, taken from input test 7
        c_int_test = copy.copy(self.dummy_vv_output)
        c_int_test.hgvs_transcript_variant = vv.hp.parse(
                'NM_000548.4:c.138+821del')
        c_int_test.genome_context_intronic_sequence = copy.copy(
                c_int_test.hgvs_transcript_variant)
        c_int_test.genome_context_intronic_sequence.rel_ac = 'NC_000016.9'
        c_int_test.hgvs_predicted_protein_consequence = {
                'prot': vv.hp.parse('NP_000539.2:p.?')}

        hgvs_hg37 = vv.hp.parse('NC_000016.9:g.2099575del')
        hgvs_hg38 = vv.hp.parse('NC_000016.10:g.2049574del')
        c_int_test.primary_assembly_loci = {
                'hg19':{'hgvs_genomic_description': hgvs_hg37},
                'hg38':{'hgvs_genomic_description': hgvs_hg38},
                'grch37':{'hgvs_genomic_description': hgvs_hg37},
                'grch38':{'hgvs_genomic_description': hgvs_hg38}}
        c_int_test.hgvs_genomic = copy.copy(hgvs_hg37)

        vrs_out = self.vrs.variant_validator_output_set_to_vrs(c_int_test)
        warnings_and_messages = {
            'validation_warnings': [],
            'lovd_messages': [],
            'lovd_corrections': [],
            'vrs_output_warnings': [
                'VRSIntronWarning: VRS does not handle mappings shared between '
                'genomic and transcript reference sequences. As such only the '
                'genomic mappings for this hgvs intronic transcript variant are'
                ' preserved.']}
        vrs_intronic_genomic_variant = {
           'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.W6wLoIFOn4G7cjopxPxYNk2lcEqhLQFb',
                    'label': 'NC_000016.9'},
                'start': 2099572,
                'end': 2099575,
                'id': 'ga4gh:SL.OSRCAmPHLKeVCIccBJw2PE2fjBBZVO5g'},
            'state': {
                'type': 'ReferenceLengthExpression',
                'length': 2, 'repeatSubunitLength': 1},
            'id': 'ga4gh:VA.R4sVDvRvhoB_h0TTNcxbRkEkqL8Iwzww'}
        VRS_mappings_for_primary_assemblies = {
         'grch37': vrs_intronic_genomic_variant,
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.yC_0RBj3fgBlvgyAuycbzdubtLxq-rE0',
            'label': 'NC_000016.10'},
           'start': 2049571,
           'end': 2049574,
           'id': 'ga4gh:SL.LBYz9IIcNSR00VhmUb06RRpf9WQKh_xN'},
          'state': {
           'type': 'ReferenceLengthExpression',
           'length': 2,
           'repeatSubunitLength': 1},
          'id': 'ga4gh:VA.swtaZ9h8m_gc6bDXPtNAALWSRBquhDcs'}}
        assert vrs_out['warnings_and_messages']['vrs_output_warnings'] == \
                warnings_and_messages['vrs_output_warnings']
        assert vrs_out['warnings_and_messages'] == warnings_and_messages
        assert vrs_out['vrs_intronic_genomic_variant'] == \
                vrs_intronic_genomic_variant
        assert vrs_out['VRS_mappings_for_primary_assemblies'] == \
                VRS_mappings_for_primary_assemblies

    def test_vv_output_to_vrs_output_intronic_alt_chr_match(self):
        # 4b intronic matches alt genomic ref, adapted from input test 189
        # also tests RSG+ intronic RSG
        c_int_alt = copy.copy(self.dummy_vv_output)
        c_int_alt.hgvs_transcript_variant = vv.hp.parse(
                'NM_012309.4:c.913-5058G>A')
        # orig 'genome_context_intronic_sequence' = \
        #        'NC_000011.10(NM_012309.4):c.913-5058G>A'
        # we switch to NW_004070871.1 which would normally only turn up in
        # output if used in input to test such input handling
        c_int_alt.hgvs_transcript_variant.rel_ac = 'NW_004070871.1'
        c_int_alt.refseqgene_context_intronic_sequence =  vv.hp.parse(
                'NM_012309.4:c.913-5058G>A')
        c_int_alt.refseqgene_context_intronic_sequence.rel_ac = 'NG_042866.1'
        c_int_alt.hgvs_refseqgene_variant = vv.hp.parse(
                'NG_042866.1:g.149464G>A')
        c_int_alt.hgvs_predicted_protein_consequence = {
                'prot': vv.hp.parse('NP_036441.2:p.?')}
        hgvs_hg37_alt = vv.hp.parse('NW_004070871.1:g.574546C>T')
        c_int_alt.alt_genomic_loci = [{
                'grch37': {
                    'hgvs_genomic_description': hgvs_hg37_alt,
            }}, {
                'hg19': {
                    'hgvs_genomic_description': hgvs_hg37_alt,
            }}]
        hgvs_hg37 = vv.hp.parse('NC_000011.10:g.71080333C>T')
        c_int_alt.primary_assembly_loci = {
                'hg38': {'hgvs_genomic_description': hgvs_hg37},
                'grch38': {'hgvs_genomic_description': hgvs_hg37}}
        c_int_alt.hgvs_genomic = hgvs_hg37_alt
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(c_int_alt)
        warnings_and_messages = {
                'validation_warnings': [],
                'lovd_messages': [],
                'lovd_corrections': [],
                'vrs_output_warnings': [
                    'VRSIntronWarning: VRS does not handle mappings shared '
                    'between genomic and transcript reference sequences. As '
                    'such only the genomic mappings for this hgvs intronic '
                    'transcript variant are preserved.']}
        VRS_mappings_for_primary_assemblies = {
         'grch38': {
          'type': 'Allele',
          'location': {
            'type': 'SequenceLocation',
            'sequenceReference': {
             'type': 'SequenceReference',
             'refgetAccession': 'SQ.2NkFm8HK88MqeNkCgj78KidCAXgnsfV1',
             'label': 'NC_000011.10'},
            'start': 71080332,
            'end': 71080333,
            'id': 'ga4gh:SL.s01bDxQQkHA1YafVKjPKqmPdy_BCgwWp'},
           'state': {
            'type': 'LiteralSequenceExpression',
            'sequence': 'T'},
           'id': 'ga4gh:VA.Exskx3X6aUlO_msp7n5jyBYUDkiHG2Xw'}}
        VRS_mappings_for_alt_genomic_loci = [{
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.EZG3Q7xvM2XQ-HQTX5w1BEEmWo-AGqtB',
                    'label': 'NW_004070871.1'},
                'start': 574545,
                'end': 574546,
                'id': 'ga4gh:SL.c9V2DUbRs4P15CWnthnRRgptVohYcy__'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'T'},
            'id': 'ga4gh:VA.BzVOTYucxGEVNsslkihaJ_kgH9GN5Tro'}]
        vrs_intronic_genomic_variant = {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.EZG3Q7xvM2XQ-HQTX5w1BEEmWo-AGqtB',
                    'label': 'NW_004070871.1'},
                'start': 574545,
                'end': 574546,
                'id': 'ga4gh:SL.c9V2DUbRs4P15CWnthnRRgptVohYcy__'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'T'},
            'id': 'ga4gh:VA.BzVOTYucxGEVNsslkihaJ_kgH9GN5Tro'}
        assert vrs_out['warnings_and_messages'] == warnings_and_messages
        assert vrs_out['VRS_mappings_for_primary_assemblies'] == \
                VRS_mappings_for_primary_assemblies
        assert vrs_out['VRS_mappings_for_alt_genomic_loci'] == \
                VRS_mappings_for_alt_genomic_loci
        assert vrs_out['vrs_intronic_genomic_variant'] == \
                vrs_intronic_genomic_variant

    def test_vv_output_to_vrs_output_intergenic_type_data(self):
        # test 5 g type input (intergenic/psudo-intergenic variant input)
        # taken from input test 1
        intergenic_vv_res = copy.copy(self.dummy_vv_output)
        # submitted_variant NC_000017.10:g.48279242G>T
        intergenic_vv_res.hgvs_refseqgene_variant = vv.hp.parse(
                'NG_007400.1:g.4759C>A')
        hgvs_hg37 = vv.hp.parse('NC_000017.10:g.48279242G>T')
        hgvs_hg38 = vv.hp.parse('NC_000017.11:g.50201881G>T')
        intergenic_vv_res.hgvs_genomic = hgvs_hg38
        intergenic_vv_res.primary_assembly_loci = {
                'hg19': {'hgvs_genomic_description': hgvs_hg37},
                'hg38': {'hgvs_genomic_description': hgvs_hg38},
                'grch37': {'hgvs_genomic_description': hgvs_hg37},
                'grch38': {'hgvs_genomic_description': hgvs_hg38}}
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(intergenic_vv_res)
        warnings_and_messages = {
                'validation_warnings': [],
                'lovd_messages': [],
                'lovd_corrections': []}
        vrs_refseqgene_variant = {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.QHD1sq0MO0ekfVC3dyTMD76jTTuJdHzC',
                    'label': 'NG_007400.1'},
                'start': 4758,
                'end': 4759,
                'id': 'ga4gh:SL.Zo3m-4l8X4ec1R9WJoFNnRZJpsf1vxD-'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'A'},
            'id': 'ga4gh:VA.T6RVcuFrzTqxCS16D3stnBB8uaf86LoT'}
        VRS_mapping_for_grch37 = {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.AjWXsI7AkTK35XW9pgd3UbjpC3MAevlz',
                    'label': 'NC_000017.10'},
                'start': 48279241,
                'end': 48279242,
                'id': 'ga4gh:SL.BTi2Q9MAcw-eIZjsascBF4kPH7JFHiwK'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'T'},
            'id': 'ga4gh:VA.Mu21lzIQgKp1Bzq0-E-DiVNNVS1GeF9l'}
        VRS_mapping_for_grch38 = {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7',
                    'label': 'NC_000017.11'},
                'start': 50201880,
                'end': 50201881,
                'id': 'ga4gh:SL.4sKmwc3lN8mn6DxNKFmV_piFt6zW4nkD'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'T'},
            'id': 'ga4gh:VA.VMA149-dG25GHnAHd76yLApQ_5lMS2-4'}
        VRS_mappings_for_primary_assemblies = {
            'grch37': VRS_mapping_for_grch37,
            'grch38': VRS_mapping_for_grch38}
        assert vrs_out['vrs_intergenic_genomic_variant'] == \
                VRS_mapping_for_grch38
        assert vrs_out['VRS_mappings_for_alt_genomic_loci'] == []
        assert vrs_out['warnings_and_messages'] == warnings_and_messages
        assert vrs_out['VRS_mappings_for_primary_assemblies'] == \
                VRS_mappings_for_primary_assemblies
        assert vrs_out['vrs_refseqgene_variant'] == vrs_refseqgene_variant

    def test_vv_output_to_vrs_output_intronic(self):
        # test variant.hgvs_refseqgene_variant used for output with transcript data
        ex_vv_res = copy.copy(self.dummy_vv_output)
        ex_vv_res.hgvs_transcript_variant = vv.hp.parse(
                'NM_000548.3:c.138+821del')
        ex_vv_res.hgvs_refseqgene_variant = vv.hp.parse(
                'NG_005895.1:g.5269del')
        hgvs_hg37 = vv.hp.parse('NC_000016.9:g.2099575del')
        ex_vv_res.selected_assembly = 'GRCh37'
        ex_vv_res.primary_assembly_loci = {
                'hg19': {'hgvs_genomic_description': hgvs_hg37},
                'grch37': {'hgvs_genomic_description': hgvs_hg37}}
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(ex_vv_res)
        print(vrs_out)
        assert vrs_out['VRS_mappings_for_primary_assemblies'] == { 'grch37': {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.W6wLoIFOn4G7cjopxPxYNk2lcEqhLQFb',
                    'label': 'NC_000016.9'},
                'start': 2099572,
                'end': 2099575,
                'id': 'ga4gh:SL.OSRCAmPHLKeVCIccBJw2PE2fjBBZVO5g'},
            'state': {
                'type': 'ReferenceLengthExpression',
                'length': 2,
                'repeatSubunitLength': 1},
            'id': 'ga4gh:VA.R4sVDvRvhoB_h0TTNcxbRkEkqL8Iwzww'}}

        assert vrs_out['vrs_refseqgene_variant'] == {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.Zuvuz5j8xEuvGT-iyjeCPW96y9eJ-2Rt',
                    'label': 'NG_005895.1'},
                'start': 5266,
                'end': 5269,
                'id': 'ga4gh:SL.7lSDh7A9L3fdqd91RzMX5JcLXzKoRL4m'},
            'state': {
                'type': 'ReferenceLengthExpression',
                'length': 2, 'repeatSubunitLength': 1},
            'id': 'ga4gh:VA.3-TnA4NFzRb1rPNumnbdwe6vzFK93Dy3'}

class TestVRSOutputRangeTypeFromVV(unittest.TestCase):
    """
    Test VRS range response to VV uncertain location input types.
    """
    @classmethod
    def setup_class(self):
        self.vrs = HGVS_to_VRS()
        # core set of dummy input
        self.dummy_vv_output = MockVVData(
                selected_assembly='GRCh38',
                original = 'test_input',
                warnings=[],
                lovd_messages = [],
                lovd_corrections = [],
                stable_gene_ids = {
                    'hgnc_id': 'HGNC:TST',
                    'entrez_gene_id': '0000',
                    'ucsc_id': 'uc0tst.1',
                    'omim_id': ['0000000']},
                gene_symbol = 'Dummy_S',
                description = "Text transcript description",)
    def _two_span_var(self, ac, var_type, s1, e1, s2, e2, edit):
        if var_type in ['g','m','o']:
            edit = vv.hp.parse_dna_edit(edit)
        else:
            edit = vv.hp.parse_rna_edit(edit)
        o_start_pos, start_end = _hgvs_offset_pos_from_str_in(
                s1,
                None,
                ref_type=var_type,
                end=e1)
        start = Interval(start=o_start_pos,end=start_end)
        end_start, o_end_pos = _hgvs_offset_pos_from_str_in(
                s2,
                None,
                ref_type=var_type,
                end=e2)
        end = Interval(start=end_start,end=o_end_pos)
        full_obj = hgvs_obj_from_existing_edit(
                ac,
                var_type,
                FEInterval(
                    start=start,
                    end=end),
                edit)
        return full_obj

        # possible data for additional tests?
        #SequenceLocation (redundant data == to below skipped)
        #    in:
        #      end: [44908822, 44908922]
        #      start: [44908721, 44908821]
        #     out:
        #      ga4gh_digest: 8-sGv9AY7GJT6QVgqbxhMXFNamnWcFJu
        #      ga4gh_identify: ga4gh:SL.8-sGv9AY7GJT6QVgqbxhMXFNamnWcFJu
        #      ga4gh_serialize: '{"end":[44908822,44908922],"sequenceReference":{"refgetAccession":"SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","type":"SequenceReference"},"start":[44908721,44908821],"type":"SequenceLocation"}'
        #  - name: "SequenceLocation w/Definite and Indefinite Ranges"
        #    in:
        #      end: [44908822, null]
        #      start: [44908721, 44908821]
        #    out:
        #      ga4gh_digest: XQAXpesghmuDHziAcDCAmESBOPKTBhwD
        #      ga4gh_identify: ga4gh:SL.XQAXpesghmuDHziAcDCAmESBOPKTBhwD
        #      ga4gh_1_3_identify: ga4gh:VSL.zGh9Zy42Zu9R0sbyB1rXsxd33BIyiORk
        #      ga4gh_serialize: '{"end":[44908822,null],"sequenceReference":{"refgetAccession":"SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","type":"SequenceReference"},"start":[44908721,44908821],"type":"SequenceLocation"}'
        #      ga4gh_1_3_serialize: '{"interval":{"end":{"comparator":">=","type":"IndefiniteRange","value":44908822},"start":{"max":44908821,"min":44908721,"type":"DefiniteRange"},"type":"SequenceInterval"},"sequence_id":"F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","type":"SequenceLocation"}'
        #  - name: "SequenceLocation w/more Definite and Indefinite Ranges"
        #    in:
        #      end: [null, 44908822]
        #      start: [44908721, 44908821]
        #    out:
        #      ga4gh_digest: OYplG0vkUojmK2hDejylSykx-np3HPFP
        #      ga4gh_identify: ga4gh:SL.OYplG0vkUojmK2hDejylSykx-np3HPFP
        #      ga4gh_serialize: '{"end":[null,44908822],"sequenceReference":{"refgetAccession":"SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","type":"SequenceReference"},"start":[44908721,44908821],"type":"SequenceLocation"}'

    def test_vrs_range_unc_del(self):
        # test unknown range input for del ins and dup
        # 'NM_001377405.1:c.(4_246)delN[15]'
        # "NC_000003.12:g.(63912602_63912844)delNNNNNNNNNNNNNNN"
        tx_unk_var = vv.hp.parse('NM_001377405.1:c.(4_246)delNNNNNNNNNNNNNNN')
        gen_unk_var = vv.hp.parse(
                "NC_000003.12:g.(63912602_63912844)delNNNNNNNNNNNNNNN")
        unc_test = copy.copy(self.dummy_vv_output)
        unc_test.hgvs_transcript_variant = tx_unk_var
        unc_test.primary_assembly_loci = {
                'hg38': {
                    'hgvs_genomic_description': gen_unk_var,
                    },
                'grch38': {
                    'hgvs_genomic_description': gen_unk_var,
                    }
                }
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(unc_test)
        warnings_and_messages = {
                'validation_warnings': [],
                'lovd_messages': [],
                'lovd_corrections': []}
        assert vrs_out['warnings_and_messages'] == warnings_and_messages
        vrs_transcript_variant = {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.xVHtmFKoGeLwmbDXLl_nxw0q4c2FRMYg',
                    'label': 'NM_001377405.1'},
                'start': [405, 648],
                'end': [405, 648],
                'id': 'ga4gh:SL.wV0h-QOnWPUSj3T9ImdcphfFnqveZIMq'},
            'state': {
                'type': 'LengthExpression',
                'length': [0,228]},# this would better be something like -15
                                   # but does not translate well into VRS
            'id': 'ga4gh:VA.EpZ3bKdZlmN7K4mmRaLPtXavbriZJCT5'}
        assert vrs_out['vrs_transcript_variant'] == vrs_transcript_variant
        VRS_mappings_for_primary_assemblies = {
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX',
            'label': 'NC_000003.12'},
           'start': [63912601,63912844],
           'end': [63912601,63912844],
           'id': 'ga4gh:SL.u2JpyMzpWPiqASuoO6yvXH_615hlTpRQ'},
          'state': {
           'type': 'LengthExpression',
           'length': [0,228]},# as before this is the post-change state
                              # 15 bases del from unc pos maps badly
          'id': 'ga4gh:VA.Nt3i0Fc2JYErd93LznLg1G-l3hYTSN80'}}
        assert vrs_out['VRS_mappings_for_primary_assemblies'] == \
                VRS_mappings_for_primary_assemblies
        VRS_mappings_for_alt_genomic_loci = []
        assert vrs_out['VRS_mappings_for_alt_genomic_loci'] ==\
                VRS_mappings_for_alt_genomic_loci

    def test_vrs_range_unc_ins(self):
        # now test ins delins dup and sub
        # 'NM_001377405.1:c.(4_246)insTTTT'
        # "NC_000003.12:g.(63912602_63912844)insTTT"
        tx_unk_var = vv.hp.parse('NM_001377405.1:c.(4_246)insTTTT')
        gen_unk_var = vv.hp.parse("NC_000003.12:g.(63912602_63912844)insTTTT")
        unc_test = copy.copy(self.dummy_vv_output)
        unc_test.hgvs_transcript_variant = tx_unk_var
        unc_test.primary_assembly_loci = {
                'hg38': {
                    'hgvs_genomic_description': gen_unk_var,
                    },
                'grch38': {
                    'hgvs_genomic_description': gen_unk_var,
                    }
                }
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(unc_test)
        vrs_transcript_variant = {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.xVHtmFKoGeLwmbDXLl_nxw0q4c2FRMYg',
                    'label': 'NM_001377405.1'},
                'start': [405, 648],
                'end': [405, 648],
                'id': 'ga4gh:SL.wV0h-QOnWPUSj3T9ImdcphfFnqveZIMq'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'TTTT'},
            'id': 'ga4gh:VA.HLa3ZXHcJk2LPUFAX8d7qYE-NJcIVl0D'}
        VRS_mappings_for_primary_assemblies = {
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX',
            'label': 'NC_000003.12'},
           'start': [63912601, 63912844],
           'end': [63912601, 63912844],
           'id': 'ga4gh:SL.u2JpyMzpWPiqASuoO6yvXH_615hlTpRQ'},
          'state': {
           'type': 'LiteralSequenceExpression',
           'sequence': 'TTTT'},
          'id': 'ga4gh:VA.KsKPm1k9vo1xoWj2AsKbzWrXzcKZ_rdG'}}
        assert vrs_out['vrs_transcript_variant'] == vrs_transcript_variant
        assert vrs_out['VRS_mappings_for_primary_assemblies'] ==\
                VRS_mappings_for_primary_assemblies
        assert vrs_out['VRS_mappings_for_alt_genomic_loci'] == []

    def test_vrs_range_unc_delins(self):
        # 'NM_001377405.1:c.(4_246)delinsTTTT'
        # "NC_000003.12:g.(63912602_63912844)delinsTTT"
        tx_unk_var = vv.hp.parse('NM_001377405.1:c.(4_246)delinsTTTT')
        gen_unk_var = vv.hp.parse(
                "NC_000003.12:g.(63912602_63912844)delinsTTTT")
        unc_test = copy.copy(self.dummy_vv_output)
        unc_test.hgvs_transcript_variant = tx_unk_var
        unc_test.primary_assembly_loci = {
                'hg38': {
                    'hgvs_genomic_description': gen_unk_var,
                    },
                'grch38': {
                    'hgvs_genomic_description': gen_unk_var,
                    }
                }
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(unc_test)
        vrs_transcript_variant = {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.xVHtmFKoGeLwmbDXLl_nxw0q4c2FRMYg',
                    'label': 'NM_001377405.1'},
                'start': [405, 648],
                'end': [405, 648],
                'id': 'ga4gh:SL.wV0h-QOnWPUSj3T9ImdcphfFnqveZIMq'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'TTTT'},
            'id': 'ga4gh:VA.HLa3ZXHcJk2LPUFAX8d7qYE-NJcIVl0D'}
        VRS_mappings_for_primary_assemblies = {
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX',
            'label': 'NC_000003.12'},
           'start': [63912601, 63912844],
           'end': [63912601, 63912844],
           'id': 'ga4gh:SL.u2JpyMzpWPiqASuoO6yvXH_615hlTpRQ'},
          'state': {
           'type': 'LiteralSequenceExpression',
           'sequence': 'TTTT'},
          'id': 'ga4gh:VA.KsKPm1k9vo1xoWj2AsKbzWrXzcKZ_rdG'}}
        assert vrs_out['vrs_transcript_variant'] == vrs_transcript_variant
        assert vrs_out['VRS_mappings_for_primary_assemblies'] == \
                VRS_mappings_for_primary_assemblies
        assert vrs_out['VRS_mappings_for_alt_genomic_loci'] == []

    def test_vrs_range_unc_dup(self):
        # 'NM_001377405.1:c.(4_246)dup'
        # "NC_000003.12:g.(63912602_63912844)dup"
        tx_unk_var = vv.hp.parse('NM_001377405.1:c.(4_246)dup')
        gen_unk_var = vv.hp.parse("NC_000003.12:g.(63912602_63912844)dup")
        unc_test = copy.copy(self.dummy_vv_output)
        unc_test.hgvs_transcript_variant = tx_unk_var
        unc_test.primary_assembly_loci = {
                'hg38': {
                    'hgvs_genomic_description': gen_unk_var,
                    },
                'grch38': {
                    'hgvs_genomic_description': gen_unk_var,
                    }
                }
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(unc_test)
        vrs_transcript_variant = {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.xVHtmFKoGeLwmbDXLl_nxw0q4c2FRMYg',
                    'label': 'NM_001377405.1'},
                'start': [405, 648],
                'end': [405, 648],
                'id': 'ga4gh:SL.wV0h-QOnWPUSj3T9ImdcphfFnqveZIMq'},
            'state': {
                'type': 'LengthExpression',
                'length': [0, 486]},
            'id': 'ga4gh:VA.oqDJt3ITorI0hrvsS7LmBUx-mBtQZRS7'}
        VRS_mappings_for_primary_assemblies = {
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX',
            'label': 'NC_000003.12'},
           'start': [63912601, 63912844],
           'end': [63912601, 63912844],
           'id': 'ga4gh:SL.u2JpyMzpWPiqASuoO6yvXH_615hlTpRQ'},
          'state': {
           'type': 'LengthExpression',
           'length': [0, 486]},
          'id': 'ga4gh:VA.gtya_8q_3eubU0Hijsw_bR1DNU3KoPo4'}}
        assert vrs_out['vrs_transcript_variant'] == vrs_transcript_variant
        assert vrs_out['VRS_mappings_for_primary_assemblies'] == \
                VRS_mappings_for_primary_assemblies
        assert vrs_out['VRS_mappings_for_alt_genomic_loci'] == []

    def test_vrs_range_unc_basechanage(self):
        # 'NM_001377405.1:c.(4_246)C>G' # probably nonsense but we don't
        # (can't) normalise so it works to test the code
        # "NC_000003.12:g.(63912602_63912844)C>G"
        tx_unk_var = vv.hp.parse('NM_001377405.1:c.(4_246)C>G')
        gen_unk_var = vv.hp.parse("NC_000003.12:g.(63912602_63912844)C>G")
        unc_test = copy.copy(self.dummy_vv_output)
        unc_test.hgvs_transcript_variant = tx_unk_var
        unc_test.primary_assembly_loci = {
                'hg38': {
                    'hgvs_genomic_description': gen_unk_var,
                    },
                'grch38': {
                    'hgvs_genomic_description': gen_unk_var,
                    }
                }
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(unc_test)
        vrs_transcript_variant = {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.xVHtmFKoGeLwmbDXLl_nxw0q4c2FRMYg',
                    'label': 'NM_001377405.1'},
                'start': [405, 648],
                'end': [405, 648],
                'id': 'ga4gh:SL.wV0h-QOnWPUSj3T9ImdcphfFnqveZIMq'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'G'},
            'id': 'ga4gh:VA.6r6uUT5Xopvs9mzMgnoLGnxn0UozdQ69'}
        VRS_mappings_for_primary_assemblies = {
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX',
            'label': 'NC_000003.12'},
           'start': [63912601, 63912844],
           'end': [63912601, 63912844],
           'id': 'ga4gh:SL.u2JpyMzpWPiqASuoO6yvXH_615hlTpRQ'},
          'state': {
           'type': 'LiteralSequenceExpression',
           'sequence': 'G'},
          'id': 'ga4gh:VA.p4PEk1oAj25G_7C81VjYTNUmzxvgA1pu'}}
        assert vrs_out['vrs_transcript_variant'] == vrs_transcript_variant
        assert vrs_out['VRS_mappings_for_primary_assemblies'] == \
                VRS_mappings_for_primary_assemblies
        assert vrs_out['VRS_mappings_for_alt_genomic_loci'] == []

    # test range input with spans for each end
    #  taken from test_uncertain_3
    # 'NM_006138.4:c.(1_20)_(30_36)del'
    # "grch38"  'NC_000011.10:g.(60061161_60061180)_(60061190_60061196)del'

    def test_vrs_range_unc_span_intronic(self):
        # test that intronic detection works on span of span types
        tx_unk_var_intr = self._two_span_var(
                'NM_006138.4','c','1','20','30','36','del')
        tx_unk_var_intr.posedit.pos.start.start.offset = True
        assert self.vrs._intronic(tx_unk_var_intr)

    def test_vrs_range_unc_span_del(self):
        #'NM_006138.4:c.(1_20)_(30_36)del'
        tx_unk_var = self._two_span_var(
                'NM_006138.4','c','1','20','30','36','del')
        gen_unk_var =  self._two_span_var(
                'NC_000011.10','g',
                '60061161','60061180',
                '60061190','60061196',
                'del')
        unc_multi_end_test = copy.copy(self.dummy_vv_output)
        unc_multi_end_test.hgvs_transcript_variant = tx_unk_var
        unc_multi_end_test.primary_assembly_loci = {
                'hg38': {
                    'hgvs_genomic_description': gen_unk_var,
                    },
                'grch38': {
                    'hgvs_genomic_description': gen_unk_var,
                    }
                }
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(unc_multi_end_test)
        vrs_transcript_variant = {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.kDpc3QS9BP_MPB7-8c32tLm_ofkpKuki',
                    'label': 'NM_006138.4'},
                'start': [128, 148],
                'end': [157, 164],
                'id': 'ga4gh:SL.VKym0ZocRoqcMIFALbtocM9G-JXI_NPq'},
            'state': {
                'type': 'LengthExpression',
                'length': [9, 36]},
            'id': 'ga4gh:VA._QbcvpURg3S57MKEBF5kNXJagBt8LD3Q'}
        VRS_mappings_for_primary_assemblies = {
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.2NkFm8HK88MqeNkCgj78KidCAXgnsfV1',
            'label': 'NC_000011.10'},
           'start': [60061160, 60061180],
           'end': [60061189, 60061196],
           'id': 'ga4gh:SL.3SYax-eEvXhMFCS3jWUrqpOQmBwbt5Xd'},
          'state': {
           'type': 'LengthExpression',
           'length': [9, 36]},
          'id': 'ga4gh:VA.BimsbJi2CRmY5zw8dTShvosXgdW4KQGf'}}
        assert vrs_out['vrs_transcript_variant'] == vrs_transcript_variant
        assert vrs_out['VRS_mappings_for_primary_assemblies'] == VRS_mappings_for_primary_assemblies
        assert vrs_out['VRS_mappings_for_alt_genomic_loci'] == []

    def test_vrs_range_unc_span_delins(self):
        #'NM_006138.4:c.(1_30)_(20_36)delinsCCC'
        tx_unk_var = self._two_span_var(
                'NM_006138.4','c','1','20','30','36','delinsCCC')
        gen_unk_var =  self._two_span_var(
                'NC_000011.10','g',
                '60061161','60061180',
                '60061190','60061196',
                'delinsCCC')
        unc_multi_end_test = copy.copy(self.dummy_vv_output)
        unc_multi_end_test.hgvs_transcript_variant = tx_unk_var
        unc_multi_end_test.primary_assembly_loci = {
                'hg38': {
                    'hgvs_genomic_description': gen_unk_var,
                    },
                'grch38': {
                    'hgvs_genomic_description': gen_unk_var,
                    }
                }
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(
                unc_multi_end_test)
        vrs_transcript_variant = {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.kDpc3QS9BP_MPB7-8c32tLm_ofkpKuki',
                    'label': 'NM_006138.4'},
                'start': [128, 148],
                'end': [157, 164],
                'id': 'ga4gh:SL.VKym0ZocRoqcMIFALbtocM9G-JXI_NPq'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'CCC'},
            'id': 'ga4gh:VA.hJy6_vMrmQUVQQzu3ZBmPboypcESYKmk'}
        VRS_mappings_for_primary_assemblies = {
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.2NkFm8HK88MqeNkCgj78KidCAXgnsfV1',
            'label': 'NC_000011.10'},
           'start': [60061160, 60061180],
           'end': [60061189, 60061196],
           'id': 'ga4gh:SL.3SYax-eEvXhMFCS3jWUrqpOQmBwbt5Xd'},
          'state': {
           'type': 'LiteralSequenceExpression',
           'sequence': 'CCC'},
          'id': 'ga4gh:VA.W4I8QnB7DOZzL5LFvdVLwWeEI0BdEWvX'}}
        assert vrs_out['vrs_transcript_variant'] == vrs_transcript_variant
        assert vrs_out['VRS_mappings_for_primary_assemblies'] == \
                VRS_mappings_for_primary_assemblies
        assert vrs_out['VRS_mappings_for_alt_genomic_loci'] == []

    def test_vrs_range_unc_span_ins(self):
        #ins
        #'NM_006138.4:c.(1_30)_(20_36)insCCC'
        tx_unk_var = self._two_span_var(
                'NM_006138.4','c','1','30','20','36','insCCC')
        gen_unk_var =  self._two_span_var(
                'NC_000011.10','g',
                '60061161','60061190',
                '60061180','60061196',
                'insCCC')
        unc_multi_end_test = copy.copy(self.dummy_vv_output)
        unc_multi_end_test.hgvs_transcript_variant = tx_unk_var
        unc_multi_end_test.primary_assembly_loci = {
                'hg38': {
                    'hgvs_genomic_description': gen_unk_var,
                    },
                'grch38': {
                    'hgvs_genomic_description': gen_unk_var,
                    }
                }
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(
                unc_multi_end_test)
        vrs_transcript_variant = {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.kDpc3QS9BP_MPB7-8c32tLm_ofkpKuki',
                    'label': 'NM_006138.4'},
                'start': [128, 158],
                'end': [147, 164],
                'id': 'ga4gh:SL.aXpOvtXB05WvJoNnaJEw8XCUJFKJfSRT'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'CCC'},
            'id': 'ga4gh:VA.mTaxUtyrZx27pR7KgWFH392KzesV6PNm'}
        VRS_mappings_for_primary_assemblies =  {
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.2NkFm8HK88MqeNkCgj78KidCAXgnsfV1',
            'label': 'NC_000011.10'},
           'start': [60061160, 60061190],
           'end': [60061179, 60061196],
           'id': 'ga4gh:SL.HJElyJSuBdeXx6zYpvI6kl2gRrmYYYAw'},
          'state': {
           'type': 'LiteralSequenceExpression',
           'sequence': 'CCC'},
          'id': 'ga4gh:VA.zqcJ7xJM5ucY5wAEEBf-VWcqnyRinBL8'}}

        assert vrs_out['vrs_transcript_variant'] == vrs_transcript_variant
        assert vrs_out['VRS_mappings_for_primary_assemblies'] == \
                VRS_mappings_for_primary_assemblies


    def test_vrs_range_unc_span_dup(self):
        #'NM_006138.4:c.(1_30)_(20_36)dup'
        tx_unk_var = self._two_span_var(
                'NM_006138.4','c','1','20','30','36','dup')
        gen_unk_var =  self._two_span_var(
                'NC_000011.10','g',
                '60061161','60061180',
                '60061190','60061196',
                'dup')
        unc_multi_end_test = copy.copy(self.dummy_vv_output)
        unc_multi_end_test.hgvs_transcript_variant = tx_unk_var
        unc_multi_end_test.primary_assembly_loci = {
                'hg38': {
                    'hgvs_genomic_description': gen_unk_var,
                    },
                'grch38': {
                    'hgvs_genomic_description': gen_unk_var,
                    }
                }
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(
                unc_multi_end_test)
        vrs_transcript_variant = {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.kDpc3QS9BP_MPB7-8c32tLm_ofkpKuki',
                    'label': 'NM_006138.4'},
                'start': [128, 148],
                'end': [157, 164],
                'id': 'ga4gh:SL.VKym0ZocRoqcMIFALbtocM9G-JXI_NPq'},
            'state': {
                'type': 'LengthExpression',
                'length': [18, 72]},
            'id': 'ga4gh:VA.QldDqxFkU7wNQbB9QJvfiR3apKiWd25L'}
        VRS_mappings_for_primary_assemblies = {
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.2NkFm8HK88MqeNkCgj78KidCAXgnsfV1',
            'label': 'NC_000011.10'},
           'start': [60061160, 60061180],
           'end': [60061189, 60061196],
           'id': 'ga4gh:SL.3SYax-eEvXhMFCS3jWUrqpOQmBwbt5Xd'},
          'state': {
           'type': 'LengthExpression',
           'length': [18, 72]},
          'id': 'ga4gh:VA.6_mydc1dT4GgT4lH7hsj1kzL7CuyHi41'}}
        assert vrs_out['vrs_transcript_variant'] == vrs_transcript_variant
        assert vrs_out['VRS_mappings_for_primary_assemblies'] == \
                VRS_mappings_for_primary_assemblies
        assert vrs_out['VRS_mappings_for_alt_genomic_loci'] == []

    def test_vrs_range_unc_span_r_type_rna_var(self):
        # test same for r. type variants
        # prot is dummy but same type as expected
        # note that r type are actually C or n type as far as biocommons.hgvs
        # is concerned, r is used as similar to ":n.", but without caring about
        # coding status
        rna_test = copy.copy(self.dummy_vv_output)
        tx_unk_var = self._two_span_var(
                'NM_006138.4','c','1','20','30','36','dup')
        rna_test.rna_data = {
                'rna_variant':tx_unk_var,
                'translation_slr':vv.hp.parse('NP_055935.4:p.?'),
                'usage_warnings':"RNA warning"}
        vrs_out = self.vrs.variant_validator_output_set_to_vrs(rna_test)
        vrs_transcript_variant = {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.kDpc3QS9BP_MPB7-8c32tLm_ofkpKuki',
                    'label': 'NM_006138.4'},
                'start': [128, 148],
                'end': [157, 164],
                'id': 'ga4gh:SL.VKym0ZocRoqcMIFALbtocM9G-JXI_NPq'},
            'state': {
                'type': 'LengthExpression',
                'length': [18, 72]},
            'id': 'ga4gh:VA.QldDqxFkU7wNQbB9QJvfiR3apKiWd25L'}
        print(vrs_out)
        assert vrs_out['vrs_transcript_variant'] == vrs_transcript_variant
        assert 'VRS_mappings_for_primary_assemblies' not in vrs_out
        assert 'VRS_mappings_for_alt_genomic_loci' not in vrs_out

class TestVRSObjectInit(unittest.TestCase):
    """
    Test all remaining untested object init features, these are not normally
    exposed by VV, so may be more subject to change than might otherwise be
    expected.
    """

    def test_vrs_object_seqrepo(self):
        # test setting SR loc
        config = ConfigParser()
        config.read(CONFIG_DIR)
        seq_repo_path = os.path.join(config["seqrepo"]["location"],
                               config["seqrepo"]["version"])
        vrs = HGVS_to_VRS(seq_repo=SeqRepo(seq_repo_path))
        #do ID fetch based on set seqrepo as above
        vrs_seq_id = vrs.id_fetch('NR_110010.2')
        assert vrs_seq_id == 'SQ.Ksaa229gzGC4Uc3ytCixai4vziud-MKi'

    # there is an additional (for now untested fallback for when
    # VariantValidator.settings import CONFIG_DIR
    # fails, but testing this here is hard, and since VV can't work in this case
    # is redundant

    def test_vrs_object_init_id_fetch(self):
        # test setting ID fetch method
        def _test_id_fetch(input_id):
            if input_id == 'test_true':
                return 'True'
            return False
        vrs = HGVS_to_VRS(id_fetch=_test_id_fetch)
        vrs_seq_id = vrs.id_fetch('NR_110010.2')
        assert vrs_seq_id is False
        vrs_seq_id = vrs.id_fetch('test_true')
        assert vrs_seq_id == 'True'

    def test_vrs_object_init_preset_cache(self):
        # dummy set a pre-filled cache with normally used and filled with
        # previous hgvs sourced VRS output, due to how common repeats VV output
        # can be(multi transcript but same genomic), but with weird contents
        # for testing. Test via hgvs_single_var_to_vrs().
        vrs = HGVS_to_VRS(cache = {'not_hgvs':{'not_vrs'}})
        vrs_out = vrs.hgvs_single_var_to_vrs('not_hgvs')
        assert vrs_out == {'not_vrs'}
        assert vrs.cache == {'not_hgvs':{'not_vrs'}}

    def test_vrs_object_init_normal_cache(self):
        # test normal cache behaviour
        vrs = HGVS_to_VRS(cache=True)
        assert vrs.cache == {}

class TestVRSObjectVV(unittest.TestCase):
    """
    Tests for VV usage of the code working as intended, should mostly be a
    simple shim on vv_output_to_vrs_output internally so not a complex test set
    at the moment.

    Testing transcript, intronic transcript, genomic, and error type results
    """
    @classmethod
    def setup_class(self):
        self.vv = Validator()
        self.vv.testing = True

    def test_format_as_vrs_tx(self):
        variant = 'NM_015120.4:c.35T>C'
        results = self.vv.validate(variant, 'GRCh37', 'all')
        resultsf = results.format_as_vrs()
        print(results.format_as_dict())
        print(resultsf)
        result = resultsf['ga4gh:VA.x8gKEBS-s638D8uuuQXVJKIgCDcuQJCW']
        assert result['selected_assembly'] == 'GRCh37'
        assert result['submitted_variant'] == 'NM_015120.4:c.35T>C'
        assert result['warnings_and_messages'] == {
                'validation_warnings': [],
                'lovd_messages': None,
                'lovd_corrections': None}
        assert result['transcript_description'] == \
                'Homo sapiens ALMS1 centrosome and basal body associated '\
                'protein (ALMS1), transcript variant 1, mRNA'
        assert result['gene_ids'] == {
                'hgnc_id': 'HGNC:428',
                'entrez_gene_id': '7840',
                'ensembl_gene_id': 'ENSG00000116127',
                'ucsc_id': 'uc032nrd.1',
                'omim_id': ['606844'],
                'ccds_ids': ['CCDS42697'],
                'current_symbol': 'ALMS1'}
        assert result['vrs_transcript_variant'] == {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.mq1Saz29534GzZ6T_jCfUjsqNlifC7pZ',
                    'label': 'NM_015120.4'
                    },
                'start': 145,
                'end': 146,
                'id': 'ga4gh:SL.x_rdZXUcg7V8Ae5pPVv0OS-faxKyzNQb'
                },
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'C'
                },
            'id': 'ga4gh:VA.x8gKEBS-s638D8uuuQXVJKIgCDcuQJCW'
            }
        assert result['vrs_refseqgene_variant']  == {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.r47Jw6neAxxIQgeyZLmdPrqFN4S6h7hy',
                    'label': 'NG_011690.1'},
                'start': 5145,
                'end': 5146,
                'id': 'ga4gh:SL.3CEnQ3JJo2jQbraJ0UHIO-lqBbR4xs15'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'C'},
            'id': 'ga4gh:VA.TjhuWhUG-PJapHphK3r0HWyU5b9CgGG_'}

        assert result['vrs_predicted_protein_consequence'] == {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.ggxH8oENu3WDprE6Z_Qg1HvOiT0d8UNa',
                    'label': 'NP_055935.4'
                    },
                'start': 11,
                'end': 12,
                'id': 'ga4gh:SL.qOpPs9dwgqK6MuH3CTQUDoRYbi9IWwSh'
                },
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'P'
                },
            'id': 'ga4gh:VA.cJPulDXE1G73sZS4dHaZX_27-c4pNjYo'
            }
        assert result['VRS_mappings_for_primary_assemblies'] == {
         'grch37': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.9KdcA9ZpY1Cpvxvg8bMSLYDUpsX6GDLO',
            'label': 'NC_000002.11'
           },
           'start': 73613030,
           'end': 73613031,
           'id': 'ga4gh:SL.IzA6GY6mvDiJ1EzMbe1Cxv6jyhg4u8t7'
          },
         'state': {
          'type': 'LiteralSequenceExpression',
          'sequence': 'CGGA'
         },
         'id': 'ga4gh:VA.rKc_HxeylEKYQzYo_e_kY1lqYLF6Ch5m'
         },
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g',
            'label': 'NC_000002.12'
            },
           'start': 73385902,
           'end': 73385903,
           'id': 'ga4gh:SL.ndPBndyntJTTRHuU6tYpgbN6TUsRv94B'
           },
          'state': {
           'type': 'LiteralSequenceExpression',
           'sequence': 'CGGA'
           },
          'id': 'ga4gh:VA.eDiJCZTUpLf1UQHaUFFP1UaPsHfbOjyC'
          }
         }
        assert result['VRS_mappings_for_alt_genomic_loci'] == [{
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.8wk0aF1q6kFjIccOvc8bq3AJVuqWersd',
                    'label': 'NW_025791766.1'},
                'start': 55397,
                'end': 55398,
                'id': 'ga4gh:SL.vknogVomLiGh3q7IfPcZDAU9Kg6nSvyY'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'C'},
            'id': 'ga4gh:VA.pdjaBvFWh8JgChnGQzOTt9s5zxhNudPS'}]
    def test_format_as_vrs_tx_intronic(self):
        # from test_variant5
        variant = 'NC_000023.10:g.33229673A>T'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_vrs()
        print(results)
        #'NM_000109.3:c.7+127703T>A'
        result = results['ga4gh:VA.6LfJDlOTAm1QeCDWfbE0986oJwQ25NSd']
        assert result['selected_assembly'] == 'GRCh37'
        assert result['submitted_variant'] == 'NC_000023.10:g.33229673A>T'
        assert result['warnings_and_messages'] == {
            'validation_warnings': [],
            'lovd_messages': None, 'lovd_corrections': None,
            'vrs_output_warnings': [
                'VRSIntronWarning: VRS does not handle mappings shared between '
                'genomic and transcript reference sequences. As such only the '
                'genomic mappings for this hgvs intronic transcript variant are'
                ' preserved.']}
        assert result['gene_ids']['hgnc_id'] == 'HGNC:2928'
        assert result['gene_ids']['entrez_gene_id'] == '1756'
        assert result['gene_ids']['ensembl_gene_id'] =='ENSG00000198947'
        assert result['gene_ids']['ucsc_id'] == 'uc004dda.2'
        assert result['transcript_description'] == \
                'Homo sapiens dystrophin (DMD), transcript variant Dp427c, mRNA'
        assert result['vrs_intronic_genomic_variant'] == {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm',
                    'label': 'NC_000023.10'},
                'start': 33229672,
                'end': 33229673,
                'id': 'ga4gh:SL.WqDN0NCy9-Y7PVGwqvI0VtcogTVZCsM4'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'T'},
            'id': 'ga4gh:VA.6LfJDlOTAm1QeCDWfbE0986oJwQ25NSd'}
        assert result['VRS_mappings_for_primary_assemblies'] == {
         'grch37':  {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm',
            'label': 'NC_000023.10'},
           'start': 33229672,
           'end': 33229673,
           'id': 'ga4gh:SL.WqDN0NCy9-Y7PVGwqvI0VtcogTVZCsM4'},
         'state': {
          'type': 'LiteralSequenceExpression',
          'sequence': 'T'},
         'id': 'ga4gh:VA.6LfJDlOTAm1QeCDWfbE0986oJwQ25NSd'},
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP',
            'label': 'NC_000023.11'},
           'start': 33211555,
           'end': 33211556,
           'id': 'ga4gh:SL.wdORKcdcBA4pikHhob-AdKP-T-IlgKbV'},
          'state': {
           'type': 'LiteralSequenceExpression',
           'sequence': 'T'},
          'id': 'ga4gh:VA.I70ooeAvjS4QEfhx9ZYZI3ews44FjqiH'}}
    def test_format_as_vrs_genomic(self):
        # from test_variant16
        variant = 'NC_000017.10:g.48279242G>T'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_vrs()
        print(results)
        result = results['ga4gh:VA.Mu21lzIQgKp1Bzq0-E-DiVNNVS1GeF9l']
        assert result['selected_assembly'] == 'GRCh37'
        assert result['submitted_variant'] == 'NC_000017.10:g.48279242G>T'
        assert result['warnings_and_messages'] == {
            'validation_warnings': [
                'TranscriptIdentificationWarning: No individual transcripts '
                'have been identified that fully overlap the described '
                'variation in the genomic sequence. Large variants might span '
                'one or more genes and are currently only described at the '
                'genome (g.) level.'],
            'lovd_messages': None,
            'lovd_corrections': None}
        assert result['gene_ids'] == {
                'hgnc_id': 'HGNC:2197',
                'entrez_gene_id': '1277',
                'ensembl_gene_id': 'ENSG00000108821',
                'ucsc_id': 'uc002iqm.4',
                'omim_id': ['120150'],
                'ccds_ids': ['CCDS11561'],
                'current_symbol': 'COL1A1'}
        assert result['transcript_description'] == ''
        assert result['vrs_intergenic_genomic_variant'] == {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.AjWXsI7AkTK35XW9pgd3UbjpC3MAevlz',
                    'label': 'NC_000017.10'},
                'start': 48279241,
                'end': 48279242,
                'id': 'ga4gh:SL.BTi2Q9MAcw-eIZjsascBF4kPH7JFHiwK'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'T'},
            'id': 'ga4gh:VA.Mu21lzIQgKp1Bzq0-E-DiVNNVS1GeF9l'}

        assert result['vrs_refseqgene_variant'] == {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequenceReference': {
                    'type': 'SequenceReference',
                    'refgetAccession': 'SQ.QHD1sq0MO0ekfVC3dyTMD76jTTuJdHzC',
                    'label': 'NG_007400.1'},
                'start': 4758,
                'end': 4759,
                'id': 'ga4gh:SL.Zo3m-4l8X4ec1R9WJoFNnRZJpsf1vxD-'},
            'state': {
                'type': 'LiteralSequenceExpression',
                'sequence': 'A'},
            'id': 'ga4gh:VA.T6RVcuFrzTqxCS16D3stnBB8uaf86LoT'}
        assert result['VRS_mappings_for_primary_assemblies'] == {
         'grch37': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.AjWXsI7AkTK35XW9pgd3UbjpC3MAevlz',
            'label': 'NC_000017.10'},
           'start': 48279241,
           'end': 48279242,
           'id': 'ga4gh:SL.BTi2Q9MAcw-eIZjsascBF4kPH7JFHiwK'},
          'state': {
           'type': 'LiteralSequenceExpression',
           'sequence': 'T'},
          'id': 'ga4gh:VA.Mu21lzIQgKp1Bzq0-E-DiVNNVS1GeF9l'},
         'grch38': {
          'type': 'Allele',
          'location': {
           'type': 'SequenceLocation',
           'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7',
            'label': 'NC_000017.11'},
           'start': 50201880,
           'end': 50201881,
           'id': 'ga4gh:SL.4sKmwc3lN8mn6DxNKFmV_piFt6zW4nkD'},
          'state': {
           'type': 'LiteralSequenceExpression',
           'sequence': 'T'},
          'id': 'ga4gh:VA.VMA149-dG25GHnAHd76yLApQ_5lMS2-4'}}

    def test_format_as_vrs_err_result(self):
        variant = 'NM_007075.3:r.235_236insGCCCACCCACCTGCCAG'
        results = self.vv.validate(variant, 'GRCh37', 'all').format_as_vrs()
        assert 'error_1' in results
        assert results['error_1'] == {
            'selected_assembly': 'GRCh37',
            'submitted_variant': 'NM_007075.3:r.235_236insGCCCACCCACCTGCCAG',
            'warnings_and_messages': {
                'validation_warnings': [
                    'The IUPAC RNA alphabet dictates that RNA variants must '
                    'use the character u in place of t'],
                'lovd_messages': None,
                'lovd_corrections': None},
            'gene_ids': {},
            'transcript_description': ''}

    def test_format_as_vrs_multi_err_result(self):
        variants = ['NM_007075.3:r.235_236insGCCCACCCACCTGCCAG',
                    'NC_000017.10:g.41232400_41236235del383']
        variants = json.dumps(variants)
        results = self.vv.validate(variants, 'GRCh37', 'all').format_as_vrs()
        assert 'error_1' in results
        assert 'error_2' in results
        assert results['error_1'] == {
            'selected_assembly': 'GRCh37',
            'submitted_variant': 'NM_007075.3:r.235_236insGCCCACCCACCTGCCAG',
            'warnings_and_messages': {
                'validation_warnings':
                ['The IUPAC RNA alphabet dictates that RNA variants must '
                 'use the character u in place of t'],
                'lovd_messages': None,
                'lovd_corrections': None},
            'gene_ids': {},
            'transcript_description': ''}
        assert results['error_2'] == {
         'selected_assembly': 'GRCh37',
         'submitted_variant': 'NC_000017.10:g.41232400_41236235del383',
         'warnings_and_messages': {
          'validation_warnings': [
           'Length implied by coordinates must equal sequence deletion length',
           'Trailing digits are not permitted in HGVS variant descriptions',
           'Refer to http://varnomen.hgvs.org/recommendations/DNA/variant/'],
          'lovd_messages': None,
          'lovd_corrections': None},
         'gene_ids': {},
         'transcript_description': ''}
