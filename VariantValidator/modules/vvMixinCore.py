'''
This module contains the main function for variant validator. It's added to the Validator object in the vvObjects file.
'''

import hgvs
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.dataproviders.seqfetcher
import hgvs.assemblymapper
import hgvs.variantmapper
import hgvs.sequencevariant
import hgvs.validator
import hgvs.exceptions
import hgvs.location
import hgvs.posedit
import hgvs.edit
import hgvs.normalizer

# IMPORT PYTHON MODULES
import re
import time
#import datetime
import copy
import os
import sys
from operator import itemgetter
#from pyliftover import LiftOver
import traceback
#from configparser import ConfigParser

#from Bio.Seq import Seq

# Import variantanalyser and peripheral VV modules
#import ref_seq_type
#import external
#import output_formatter
#import variantanalyser
from .vvLogging import logger
import hgvs
from . import vvHGVS
#from variantanalyser import functions as va_func
#from variantanalyser import dbControls as va_dbCrl
#from variantanalyser import hgvs2vcf as vvHGVS
#from variantanalyser import batch as va_btch
#from variantanalyser import g_to_g as va_g2g
#from variantanalyser import supported_chromosome_builds as va_scb
#from variantanalyser.liftover import liftover as lift_over
from .vvLiftover import liftover as lift_over #???

from . import vvFunctions as fn
from . import vvDatabase
from . import vvChromosomes
from . import vvMixinConverters
from .vvFunctions import VariantValidatorError
from .variant import Variant
from . import format_converters
from . import use_checking
from . import collect_info
from . import mappers


class Mixin(vvMixinConverters.Mixin):
    def validate(self, batch_variant, selected_assembly, select_transcripts, transcriptSet = "refseq"):
        '''
        This is the main validator function.
        :param batch_variant: A string containing the variant to be validated
        :param selected_assembly: The version of the genome assembly to use.
        :param select_transcripts: Can be an array of different transcripts, or 'all'
        Selecting multiple transcripts will lead to a multiple variant outputs.
        :param transcriptSet: 'refseq' or 'ensembl'. Currently only 'refseq' is supported
        :return:
        '''
        logger.info(batch_variant + ' : ' + selected_assembly)

        if transcriptSet == "refseq":
            alt_aln_method = 'splign'
        elif transcriptSet == "ensembl":
            alt_aln_method = 'genebuild'
            logger.warning("Ensembl is currently not supported")
            raise Exception("Ensembl is currently not supported")
        else:
            raise Exception("The transcriptSet variable '%s' is invalid, it must be 'refseq' or 'ensembl'" % transcriptSet)

        # Take start time
        start_time = time.time()

        # Set pre defined variables
        # SeqFetcher
        # sf = hgvs.dataproviders.seqfetcher.SeqFetcher()
        primary_assembly=None

        self.selected_assembly = selected_assembly
        self.select_transcripts = select_transcripts
        self.alt_aln_method = alt_aln_method

        try:
            # Validation
            ############

            # Create a dictionary of transcript ID : ''
            select_transcripts_dict = {}
            select_transcripts_dict_plus_version = {}
            if select_transcripts != 'all':
                select_transcripts_list = select_transcripts.split('|')
                for id in select_transcripts_list:
                    id = id.strip()
                    if re.match('LRG', id):
                        id = self.db.get_RefSeqTranscriptID_from_lrgTranscriptID(id)
                        if id == 'none':
                            continue
                    select_transcripts_dict_plus_version[id] = ''
                    id = id.split('.')[0]
                    select_transcripts_dict[id] = ''
            # Set up gene list dictionary
            input_genes = {}

            # Remove genes if transcripts selected
            # if select_transcripts != 'all':

            # split the batch queries into a list
            batch_queries = batch_variant.split('|')

            # Turn each variant into a dictionary. The dictionary will be compiled during validation
            self.batch_list = []
            for queries in batch_queries:
                queries = queries.strip()
                query = Variant(queries)
                self.batch_list.append(query)

            # Create List to carry batch data
            batch_out = []

            # Enter the validation loop
            ###########################
            # Allow order by input
            ordering = 0

            """
            Set a flag to mark the final output type
            flag : warning
            flag : error
            flag : intragenic
            flag : gene
            """

            logger.debug("Batch list length " + str(len(self.batch_list)))
            for my_variant in self.batch_list:
                # Start timing
                logger.traceStart(my_variant)

                # Create Normalizers
                my_variant.hn = hgvs.normalizer.Normalizer(self.hdp,
                                                cross_boundaries=False,
                                                shuffle_direction=3,
                                                alt_aln_method=alt_aln_method
                                                )
                hn = hgvs.normalizer.Normalizer(self.hdp,
                                                cross_boundaries=False,
                                                shuffle_direction=3,
                                                alt_aln_method=alt_aln_method
                                                )
                my_variant.reverse_normalizer = hgvs.normalizer.Normalizer(self.hdp,
                                                                cross_boundaries=False,
                                                                shuffle_direction=5,
                                                                alt_aln_method=alt_aln_method
                                                                )
                reverse_normalizer = hgvs.normalizer.Normalizer(self.hdp,
                                                                cross_boundaries=False,
                                                                shuffle_direction=5,
                                                                alt_aln_method=alt_aln_method
                                                                )
                # This will be used to order the final output
                if not my_variant.order:
                    ordering = ordering + 1
                    my_variant.order = ordering

                # Bug catcher
                try:
                    # Note, ID is not touched. It is always the input variant description.
                    # Quibble will be altered but id will not if type = g.
                    input = my_variant.quibble
                    logger.trace("Commenced validation of " + str(my_variant.quibble), my_variant)

                    if not my_variant.is_ascii():
                        chars, positions = my_variant.get_non_ascii()
                        error = 'Submitted variant description contains an invalid character(s) %s at position(s) %s: '\
                                'Please remove this character and re-submit: A useful search function for ' \
                                'Unicode characters can be found at https://unicode-search.net/' % (chars, positions)
                        my_variant.warnings += ': ' + error
                        logger.warning(error)
                        continue

                    # Remove whitespace
                    my_variant.remove_whitespace()
                    if my_variant.quibble != my_variant.original:
                        caution = 'Whitespace removed from variant description %s' % my_variant.original
                        my_variant.warnings += ': ' + caution
                        logger.info(caution)

                    stash_input = copy.copy(input)
                    # Set the primary_assembly
                    if not my_variant.primary_assembly:
                        if selected_assembly == 'hg19':
                            primary_assembly = 'GRCh37'
                        elif selected_assembly == 'hg38':
                            primary_assembly = 'GRCh38'
                        # Ensure genome build is correctly formatted
                        elif re.search('GRC', selected_assembly, re.IGNORECASE):
                            selected_assembly = selected_assembly.replace('g', 'G')
                            selected_assembly = selected_assembly.replace('r', 'R')
                            selected_assembly = selected_assembly.replace('c', 'C')
                            selected_assembly = selected_assembly.replace('H', 'h')
                            primary_assembly = selected_assembly
                        # Catch invalid genome build
                        if primary_assembly in self.genome_builds:
                            my_variant.primary_assembly = primary_assembly
                        else:
                            my_variant.primary_assembly = 'GRCh38'
                            primary_assembly = 'GRCh38'
                            my_variant.warnings += ': Invalid genome build has been specified. Automap has selected the default build (GRCh38)'
                            logger.warning(
                                'Invalid genome build has been specified. Automap has selected the default build ' + my_variant.primary_assembly)
                    else:
                        primary_assembly = my_variant.primary_assembly
                    logger.trace("Completed string formatting", my_variant)

                    # Set variables that batch will not use but are required
                    crossing = 'false'
                    boundary = 'false'

                    # VCF type 1
                    toskip = format_converters.vcf2hgvs_stage1(my_variant, self)
                    if toskip:
                        continue

                    # API type non-HGVS
                    # e.g. Chr16:2099572TC>T
                    toskip = format_converters.vcf2hgvs_stage2(my_variant, self)
                    if toskip:
                        continue

                    toskip = format_converters.vcf2hgvs_stage3(my_variant, self)
                    if toskip:
                        continue

                    toskip = format_converters.gene_symbol_catch(my_variant, self, select_transcripts_dict_plus_version)
                    if toskip:
                        continue

                    # NG_:c. or NC_:c.
                    toskip = format_converters.refseq_catch(my_variant, self, select_transcripts_dict_plus_version)
                    if toskip:
                        continue

                    # Find not_sub type in input e.g. GGGG>G
                    toskip = format_converters.vcf2hgvs_stage4(my_variant, self)
                    if toskip:
                        continue

                    toskip = format_converters.indel_catching(my_variant, self)
                    if toskip:
                        continue

                    # Tackle compound variant descriptions NG or NC (NM_) i.e. correctly input NG/NC_(NM_):c.
                    format_converters.intronic_converter(my_variant)

                    # Extract variants from HGVS allele descriptions
                    # http://varnomen.hgvs.org/recommendations/DNA/variant/alleles/
                    toskip = format_converters.allele_parser(my_variant, self)
                    if toskip:
                        continue

                    input = my_variant.quibble

                    print("Original: %s" % my_variant.original)
                    print("Quibble: %s" % my_variant.quibble)

                    caution = ''
                    # INITIAL USER INPUT FORMATTING
                    invalid = my_variant.format_quibble()
                    if invalid:
                        if re.search(r'\w+\:[gcnmrp]', my_variant.quibble) and not re.search(r'\w+\:[gcnmrp]\.', my_variant.quibble):
                            error = 'Variant description ' + my_variant.quibble + ' lacks the . character between <type> and <position> in the expected pattern <accession>:<type>.<position>'
                        else:
                            error = 'Variant description ' + my_variant.quibble + ' is not in an accepted format'
                        my_variant.warnings += ': ' + error
                        logger.warning(error)
                        continue

                    formatted_variant = my_variant.quibble
                    input = my_variant.quibble
                    stash_input = my_variant.quibble
                    my_variant.stashed = stash_input
                    format_type = my_variant.reftype

                    # Validator specific variables, note, not all will be necessary for batch, but keep to ensure that batch works
                    # vars = []
                    # refseq_gene = ''
                    # relevant = ''
                    warning = ''
                    automap = 'false'
                    # vmapped = 'false'
                    # coords = 'false'
                    # ensembl_gene = 'false'
                    hgnc_gene_info = 'false'
                    # issue_link = 'false'
                    # cr_available = 'false'
                    # rcmds_tab = 'false'

                    logger.trace("Variant input formatted, proceeding to validate.", my_variant)

                    # Conversions
                    # Conversions are not currently supported. The HGVS format for conversions
                    # is rarely seen wrt genomic sequencing data and needs to be re-evaluated
                    if 'con' in my_variant.quibble:
                        my_variant.warnings += ': ' + 'Gene conversions currently unsupported'
                        logger.warning('Gene conversions currently unsupported')
                        continue

                    # Primary check that hgvs will accept the variant
                    error = 'false'
                    # Change RNA bases to upper case but nothing else
                    if format_type == ":r.":
                        formatted_variant = formatted_variant.upper()
                        formatted_variant = formatted_variant.replace(':R.', ':r.')
                        # lowercase the supported variant types
                        formatted_variant = formatted_variant.replace('DEL', 'del')
                        formatted_variant = formatted_variant.replace('INS', 'ins')
                        formatted_variant = formatted_variant.replace('INV', 'inv')
                        formatted_variant = formatted_variant.replace('DUP', 'dup')

                    try:
                        input_parses = self.hp.parse_hgvs_variant(formatted_variant)
                        print(input_parses, input_parses.ac, type(input_parses.ac))
                        my_variant.hgvs_formatted = input_parses
                    except hgvs.exceptions.HGVSError as e:
                        my_variant.warnings += ': ' + str(e)
                        logger.warning(error)
                        continue

                    if 'LRG' in my_variant.hgvs_formatted.ac:
                        my_variant.hgvs_formatted.ac.replace('T', 't')
                    else:
                        my_variant.hgvs_formatted.ac = my_variant.hgvs_formatted.ac.upper()
                    if hasattr(my_variant.hgvs_formatted.posedit.edit, 'alt'):
                        if my_variant.hgvs_formatted.posedit.edit.alt is not None:
                            my_variant.hgvs_formatted.posedit.edit.alt = my_variant.hgvs_formatted.posedit.edit.alt.upper()
                    if hasattr(my_variant.hgvs_formatted.posedit.edit, 'ref'):
                        if my_variant.hgvs_formatted.posedit.edit.ref is not None:
                            my_variant.hgvs_formatted.posedit.edit.ref = my_variant.hgvs_formatted.posedit.edit.ref.upper()
                    formatted_variant = str(my_variant.hgvs_formatted)
                    input = str(my_variant.hgvs_formatted)

                    assert formatted_variant == str(my_variant.hgvs_formatted)

                    my_variant.set_quibble(str(my_variant.hgvs_formatted))

                    # ENST support needs to be re-evaluated, but is very low priority
                    # ENST not supported by ACMG and is under review by HGVS
                    if my_variant.refsource == 'ENS':
                        trap_ens_in = str(my_variant.hgvs_formatted)
                        sim_tx = self.hdp.get_similar_transcripts(my_variant.hgvs_formatted.ac)
                        for line in sim_tx:
                            print(line)
                            if line[2] and line[3] and line[4] and line[5] and line[6]:
                                print("RESET")
                                my_variant.hgvs_formatted.ac = line[1]
                                my_variant.set_quibble(str(my_variant.hgvs_formatted))
                                formatted_variant = my_variant.quibble
                                break
                        if my_variant.refsource == 'ENS':
                            error = 'Unable to map ' + my_variant.hgvs_formatted.ac + ' to an equivalent RefSeq transcript'
                            my_variant.warnings += ': ' + error
                            logger.warning(error)
                            continue
                        else:
                            my_variant.warnings += ': ' + str(trap_ens_in) + ' automapped to equivalent RefSeq transcript ' + my_variant.quibble
                            logger.warning(str(trap_ens_in) + ' automapped to equivalent RefSeq transcript ' + my_variant.quibble)
                    logger.trace("HVGS acceptance test passed", my_variant)

                    # Check whether supported genome build is requested for non g. descriptions
                    historic_assembly = 'false'
                    mapable_assemblies = {
                        'GRCh37': True,
                        'GRCh38': True,
                        'NCBI36': False
                    }
                    is_mapable = mapable_assemblies.get(primary_assembly)
                    if is_mapable:

                        # These objects cannot be moved outside of the main function because they gather data from the
                        # iuser input e.g. alignment method and genome build
                        # They initiate quickly, so no need to move them unnecessarily

                        # Create easy variant mapper (over variant mapper) and splign locked evm
                        my_variant.evm = hgvs.assemblymapper.AssemblyMapper(self.hdp,
                                                                 assembly_name=primary_assembly,
                                                                 alt_aln_method=alt_aln_method,
                                                                 normalize=True,
                                                                 replace_reference=True
                                                                 )

                        evm = hgvs.assemblymapper.AssemblyMapper(self.hdp,
                                                                 assembly_name=primary_assembly,
                                                                 alt_aln_method=alt_aln_method,
                                                                 normalize=True,
                                                                 replace_reference=True
                                                                 )

                        # Setup a reverse normalize instance and non-normalize evm
                        my_variant.no_norm_evm = hgvs.assemblymapper.AssemblyMapper(self.hdp,
                                                                         assembly_name=primary_assembly,
                                                                         alt_aln_method=alt_aln_method,
                                                                         normalize=False,
                                                                         replace_reference=True
                                                                         )
                        no_norm_evm = hgvs.assemblymapper.AssemblyMapper(self.hdp,
                                                                         assembly_name=primary_assembly,
                                                                         alt_aln_method=alt_aln_method,
                                                                         normalize=False,
                                                                         replace_reference=True
                                                                         )

                        # Create a specific minimal evm with no normalizer and no replace_reference
                        my_variant.min_evm = hgvs.assemblymapper.AssemblyMapper(self.hdp,
                                                                     assembly_name=primary_assembly,
                                                                     alt_aln_method=alt_aln_method,
                                                                     normalize=False,
                                                                     replace_reference=False
                                                                     )

                    else:
                        error = 'Mapping of ' + formatted_variant + ' to genome assembly ' + primary_assembly + ' is not supported'
                        my_variant.warnings += ': ' + error
                        logger.warning(error)
                        continue

                    # Catch interval end > interval start
                    # hgvs did/does not handle 3' UTR position ordering well. This function
                    # ensures that end pos is not > start pos wrt 3' UTRs.
                    # Also identifies some variants which span into the downstream sequence
                    # i.e. out of bounds
                    astr = re.compile(r'\*')
                    if '*' in str(my_variant.hgvs_formatted.posedit):
                        input_parses_copy = copy.deepcopy(my_variant.hgvs_formatted)
                        input_parses_copy.type = "c"
                        # Map to n. position
                        # Create easy variant mapper (over variant mapper) and splign locked evm
                        try:
                            to_n = evm.c_to_n(input_parses_copy)
                        except hgvs.exceptions.HGVSError as e:
                            fn.exceptPass()
                        else:
                            if to_n.posedit.pos.end.base < to_n.posedit.pos.start.base:
                                error = 'Interval end position < interval start position '
                                my_variant.warnings += ': ' + error
                                logger.warning(error)
                                continue
                    elif my_variant.hgvs_formatted.posedit.pos.end.base < my_variant.hgvs_formatted.posedit.pos.start.base:
                        error = 'Interval end position ' + str(my_variant.hgvs_formatted.posedit.pos.end.base) + \
                                ' < interval start position ' + str(my_variant.hgvs_formatted.posedit.pos.start.base)
                        my_variant.warnings += ': ' + error
                        logger.warning(error)
                        continue

                    assert formatted_variant == str(my_variant.hgvs_formatted)

                    # Catch missing version number in refseq
                    ref_type = re.compile(r"^N\w\w\d")
                    is_version = re.compile(r"\d\.\d")
                    en_type = re.compile(r'^ENS')
                    lrg_type = re.compile(r'LRG')
                    if my_variant.refsource == 'RefSeq' and not is_version.search(str(my_variant.hgvs_formatted)):
                        error = 'RefSeq variant accession numbers MUST include a version number'
                        my_variant.warnings += ': ' + str(error)
                        continue
                    logger.trace("HVGS interval/version mapping complete", my_variant)

                    # handle LRG inputs

                    if my_variant.refsource == 'LRG':
                        format_converters.lrg_to_refseq(my_variant, self)
                        formatted_variant = my_variant.quibble
                        input = str(my_variant.hgvs_formatted)
                        stash_input = input
                        logger.trace("LRG check for conversion to refseq completed", my_variant)

                    # Additional Incorrectly input variant capture training
                    if my_variant.refsource == 'RefSeq':
                        toskip = use_checking.refseq_common_mistakes(my_variant)
                        if toskip:
                            continue
                        logger.trace("Passed 'common mistakes' catcher", my_variant)

                    assert formatted_variant == str(my_variant.hgvs_formatted)

                    # Primary validation of the input
                    toskip = use_checking.structure_checks(my_variant, self)
                    print(toskip, my_variant.hgvs_formatted, my_variant.quibble)
                    if toskip:
                        continue
                    logger.trace("Variant structure and contents searches passed", my_variant)

                    assert formatted_variant == str(my_variant.hgvs_formatted)

                    # Mitochondrial variants
                    toskip = format_converters.mitochondrial(my_variant, self)
                    if toskip:
                        continue

                    toskip = format_converters.proteins(my_variant, self)
                    if toskip:
                        continue

                    trapped_input = str(my_variant.hgvs_formatted)
                    my_variant.trapped = trapped_input
                    toskip = format_converters.rna(my_variant, self)
                    if toskip:
                        continue

                    assert formatted_variant == str(my_variant.hgvs_formatted)

                    # COLLECT gene symbol, name and ACCESSION INFORMATION
                    # Gene symbol
                    if my_variant.reftype != ':g.':
                        toskip = collect_info.get_transcript_info(my_variant, self)
                        print(toskip, my_variant.hgvs_formatted, my_variant.hgvs_genomic)
                        if toskip:
                            continue

                    assert formatted_variant == str(my_variant.hgvs_formatted)
                    # Now start mapping from genome to transcripts

                    if my_variant.reftype == ':g.':
                        toskip = mappers.gene_to_transcripts(my_variant, self)
                        print(toskip, my_variant.hgvs_formatted, my_variant.hgvs_genomic)
                        if toskip:
                            continue

                    assert formatted_variant == str(my_variant.hgvs_formatted)
                    # TYPE = :c.

                    if format_type == ':c.' or format_type == ':n.':
                        print('hgvs_formatted:', my_variant.hgvs_formatted)
                        print('input:', input)
                        print('trapped:', my_variant.trapped)
                        print('quibble:', my_variant.quibble)
                        print('formatted_variant', formatted_variant)
                        #print(my_variant.hgvs_formatted, my_variant.trapped, input)
                        toskip = mappers.transcripts_to_gene(my_variant, self)
                        print(toskip, my_variant.hgvs_formatted)
                        if toskip:
                            print("CARRYING ON")
                            continue

                    # Set the data
                    my_variant.output_type_flag = 'gene'
                    my_variant.description = hgnc_gene_info
                    # my_variant.coding = str(hgvs_coding)
                    # my_variant.genomic_r = str(hgvs_refseq)
                    # my_variant.genomic_g = str(hgvs_genomic)
                    # my_variant.protein = str(hgvs_protein)
                    my_variant.primary_assembly = primary_assembly
                    # if gap_compensation is True:
                    #     my_variant.test_stash_tx_left = test_stash_tx_left
                    #     my_variant.test_stash_tx_right = test_stash_tx_right
                    # finish timing
                    logger.traceEnd(my_variant)
                # Report errors to User and VV admin
                except KeyboardInterrupt:
                    raise
                except:
                    my_variant.output_type_flag = 'error'
                    error = 'Validation error'
                    my_variant.warnings = str(error)
                    exc_type, exc_value, last_traceback = sys.exc_info()
                    te = traceback.format_exc()
                    tbk = [str(exc_type), str(exc_value), str(te)]
                    er = str('\n'.join(tbk))
                    logger.error(str(exc_type) + " " + str(exc_value))
                    logger.debug(er)
                    #debug
                    raise

            # Outside the for loop
            ######################
            logger.trace("End of for loop")
            # order the rows
            by_order = sorted(self.batch_list, key=lambda x: x.order)

            for variant in by_order:
                if not variant.write:
                    continue

                # Blank VCF
                #                     chr = ''
                #                     pos = ''
                #                     ref = ''
                #                     alt = ''

                # Fromulate a json type response
                dict_out = {}

                # Set gap compensation bool
                gap_compensation = True

                # warngins
                warnings = variant.warnings
                warnings = re.sub('del[GATC][GATC][GATC][GATC]+', 'del', warnings)
                warnings = re.sub('^: ', '', warnings)
                warnings = re.sub('::', ':', warnings)

                # Submitted variant
                submitted = variant.original

                # Genomic sequence variation
                genomic_variant = variant.genomic_g

                # genomic accession
                if genomic_variant != '':
                    hgvs_genomic_variant = self.hp.parse_hgvs_variant(genomic_variant)
                    genomic_variant = fn.valstr(hgvs_genomic_variant)
                    genomic_accession = hgvs_genomic_variant.ac
                else:
                    genomic_accession = ''

                # RefSeqGene variation
                refseqgene_variant = variant.genomic_r
                refseqgene_variant = refseqgene_variant.strip()
                if 'RefSeqGene' in refseqgene_variant or refseqgene_variant == '':
                    warnings = warnings + ': ' + refseqgene_variant
                    refseqgene_variant = ''
                    lrg_variant = ''
                    hgvs_refseqgene_variant = 'false'
                else:
                    hgvs_refseqgene_variant = self.hp.parse_hgvs_variant(refseqgene_variant)
                    rsg_ac = self.db.get_lrgID_from_RefSeqGeneID(str(hgvs_refseqgene_variant.ac))
                    if rsg_ac[0] == 'none':
                        lrg_variant = ''
                    else:
                        hgvs_lrg = copy.deepcopy(hgvs_refseqgene_variant)
                        hgvs_lrg.ac = rsg_ac[0]
                        lrg_variant = fn.valstr(hgvs_lrg)
                        if rsg_ac[1] == 'public':
                            pass
                        else:
                            warnings = warnings + ': The current status of ' + str(
                                hgvs_lrg.ac) + ' is pending therefore changes may be made to the LRG reference sequence'

                # Transcript sequence variation
                tx_variant = variant.coding
                if tx_variant != '':
                    if '(' in tx_variant and ')' in tx_variant:
                        tx_variant = tx_variant.split('(')[1]
                        tx_variant = tx_variant.replace(')', '')

                    # transcript accession
                    hgvs_tx_variant = self.hp.parse_hgvs_variant(tx_variant)
                    tx_variant = fn.valstr(hgvs_tx_variant)
                    hgvs_transcript_variant = self.hp.parse_hgvs_variant(tx_variant)
                    transcript_accession = hgvs_transcript_variant.ac

                    # Handle LRG
                    lrg_status = 'public'
                    lrg_transcript = self.db.get_lrgTranscriptID_from_RefSeqTranscriptID(transcript_accession)
                    if lrg_transcript == 'none':
                        lrg_transcript_variant = ''
                    else:
                        # Note - LRG availability is dependant on UTA containing the data. In some
                        # instances we will be able to display the LRG_tx without being able to
                        # display the LRG gene data

                        # if not re.search('RefSeqGene', refseqgene_variant) or refseqgene_variant != '':
                        # if hgvs_refseqgene_variant != 'RefSeqGene record not available' and hgvs_refseqgene_variant != 'false':
                        try:
                            hgvs_lrg_t = self.vm.g_to_t(hgvs_refseqgene_variant, transcript_accession)
                            hgvs_lrg_t.ac = lrg_transcript
                            lrg_transcript_variant = fn.valstr(hgvs_lrg_t)
                        except:
                            if hgvs_transcript_variant.posedit.pos.start.offset == 0 and hgvs_transcript_variant.posedit.pos.end.offset == 0:
                                hgvs_lrg_t = copy.copy(hgvs_transcript_variant)
                                hgvs_lrg_t.ac = lrg_transcript
                                lrg_transcript_variant = fn.valstr(hgvs_lrg_t)
                            else:
                                lrg_transcript_variant = ''
                else:
                    transcript_accession = ''
                    lrg_transcript_variant = ''

                # Look for intronic variants
                if transcript_accession != '' and genomic_accession != '':
                    # Remove del bases
                    str_transcript = fn.valstr(hgvs_transcript_variant)
                    hgvs_transcript_variant = self.hp.parse_hgvs_variant(str_transcript)
                    try:
                        self.vr.validate(hgvs_transcript_variant)
                    except hgvs.exceptions.HGVSError as e:
                        error = str(e)
                        if 'intronic variant' in error:
                            genome_context_transcript_variant = genomic_accession + '(' + transcript_accession + '):c.' + str(
                                hgvs_transcript_variant.posedit)
                            if refseqgene_variant != '':
                                hgvs_refseqgene_variant = self.hp.parse_hgvs_variant(refseqgene_variant)
                                refseqgene_accession = hgvs_refseqgene_variant.ac
                                hgvs_coding_from_refseqgene = self.vm.g_to_t(hgvs_refseqgene_variant,
                                                                        hgvs_transcript_variant.ac)
                                hgvs_coding_from_refseqgene = fn.valstr(hgvs_coding_from_refseqgene)
                                hgvs_coding_from_refseqgene = self.hp.parse_hgvs_variant(hgvs_coding_from_refseqgene)
                                RefSeqGene_context_transcript_variant = refseqgene_accession + '(' + transcript_accession + '):c.' + str(
                                    hgvs_coding_from_refseqgene.posedit.pos) + str(
                                    hgvs_coding_from_refseqgene.posedit.edit)
                            else:
                                RefSeqGene_context_transcript_variant = ''
                        else:
                            genome_context_transcript_variant = ''  # transcript_variant
                            RefSeqGene_context_transcript_variant = ''
                    else:
                        genome_context_transcript_variant = ''  # transcript_variant
                        RefSeqGene_context_transcript_variant = ''
                else:
                    genome_context_transcript_variant = ''
                    RefSeqGene_context_transcript_variant = ''

                # Protein description
                predicted_protein_variant = variant.protein
                if 'NP_' in predicted_protein_variant:
                    rs_p, pred_prot_posedit = predicted_protein_variant.split(':')
                    lrg_p = self.db.get_lrgProteinID_from_RefSeqProteinID(rs_p)
                    if 'LRG' in lrg_p:
                        predicted_protein_variant = rs_p + '(' + lrg_p + '):' + pred_prot_posedit

                # Gene
                if transcript_accession != '':
                    try:
                        gene_symbol = self.db.get_gene_symbol_from_transcriptID(transcript_accession)
                    except:
                        gene_symbol = 'Unable to verify gene symbol for ' + str(transcript_accession)
                else:
                    gene_symbol = ''

                # Transcript description
                transcript_description = variant.description

                # Stashed variants
                # if valid.test_stash_tx_left:
                #     test_stash_tx_left = valid.test_stash_tx_left
                # if valid.test_stash_tx_right:
                #     test_stash_tx_right = valid.test_stash_tx_right

                # Multiple genomic variants
                # multi_gen_vars = []
                if tx_variant != '':
                    multi_gen_vars, hgvs_coding = mappers.final_tx_to_multiple_genomic(variant, self, tx_variant)

                else:
                    # HGVS genomic in the absence of a transcript variant
                    if genomic_variant != '':
                        multi_gen_vars = [hgvs_genomic_variant]
                    else:
                        multi_gen_vars = []

                # Dictionaries of genomic loci
                alt_genomic_dicts = []
                primary_genomic_dicts = {}

                if len(multi_gen_vars) != 0:
                    for alt_gen_var in multi_gen_vars:
                        for build in self.genome_builds:
                            test = vvChromosomes.supported_for_mapping(alt_gen_var.ac, build)
                            if test == 'true':
                                try:
                                    vcf_dict = vvHGVS.report_hgvs2vcf(alt_gen_var, build, reverse_normalizer, self.sf)
                                except hgvs.exceptions.HGVSInvalidVariantError as e:
                                    continue
                                # Identify primary assembly positions
                                if re.match('NC_', alt_gen_var.ac):
                                    if re.match('GRC', build):
                                        primary_genomic_dicts[build.lower()] = {
                                            'hgvs_genomic_description': fn.valstr(alt_gen_var),
                                            'vcf': {'chr': vcf_dict['grc_chr'],
                                                    'pos': vcf_dict['pos'],
                                                    'ref': vcf_dict['ref'],
                                                    'alt': vcf_dict['alt']
                                                    }
                                        }

                                    else:
                                        primary_genomic_dicts[build.lower()] = {
                                            'hgvs_genomic_description': fn.valstr(alt_gen_var),
                                            'vcf': {'chr': vcf_dict['ucsc_chr'],
                                                    'pos': vcf_dict['pos'],
                                                    'ref': vcf_dict['ref'],
                                                    'alt': vcf_dict['alt']
                                                    }
                                        }
                                    if build == 'GRCh38':
                                        vcf_dict = vvHGVS.report_hgvs2vcf(alt_gen_var, 'hg38', reverse_normalizer,
                                                                          self.sf)
                                        primary_genomic_dicts['hg38'] = {
                                            'hgvs_genomic_description': fn.valstr(alt_gen_var),
                                            'vcf': {'chr': vcf_dict['ucsc_chr'],
                                                    'pos': vcf_dict['pos'],
                                                    'ref': vcf_dict['ref'],
                                                    'alt': vcf_dict['alt']
                                                    }
                                        }

                                    continue

                                else:
                                    if re.match('GRC', build):
                                        dict = {build.lower(): {'hgvs_genomic_description': fn.valstr(alt_gen_var),
                                                                'vcf': {'chr': vcf_dict['grc_chr'],
                                                                        'pos': vcf_dict['pos'],
                                                                        'ref': vcf_dict['ref'],
                                                                        'alt': vcf_dict['alt']
                                                                        }
                                                                }
                                                }
                                    else:
                                        dict = {build.lower(): {'hgvs_genomic_description': fn.valstr(alt_gen_var),
                                                                'vcf': {'chr': vcf_dict['ucsc_chr'],
                                                                        'pos': vcf_dict['pos'],
                                                                        'ref': vcf_dict['ref'],
                                                                        'alt': vcf_dict['alt']
                                                                        }
                                                                }
                                                }
                                    # Append
                                    alt_genomic_dicts.append(dict)

                                    if build == 'GRCh38':
                                        vcf_dict = vvHGVS.report_hgvs2vcf(alt_gen_var, 'hg38', reverse_normalizer,
                                                                          self.sf)
                                        dict = {'hg38': {'hgvs_genomic_description': fn.valstr(alt_gen_var),
                                                         'vcf': {'chr': vcf_dict['ucsc_chr'],
                                                                 'pos': vcf_dict['pos'],
                                                                 'ref': vcf_dict['ref'],
                                                                 'alt': vcf_dict['alt']
                                                                 }
                                                         }
                                                }
                                        # Append
                                        alt_genomic_dicts.append(dict)
                                        continue
                            else:
                                # May need to account for ALT NC_
                                pass

                # Warn not directly mapped to specified genome build
                if genomic_accession != '':
                    caution = ''
                    if primary_assembly.lower() not in list(primary_genomic_dicts.keys()):
                        warnings = warnings + ': ' + str(
                            hgvs_coding) + ' cannot be mapped directly to genome build ' + primary_assembly + ': See alternative genomic loci or alternative genome builds for aligned genomic positions'

                warn_list = warnings.split(': ')
                warnings_out = []
                for warning in warn_list:
                    warning.strip()
                    warning = warning.replace("'", "")
                    if warning == '':
                        continue
                    warnings_out.append(warning)
                # Remove duplicate elements but maintain the order
                seen = {}
                no_rep_list = [seen.setdefault(x, x) for x in warnings_out if x not in seen]
                warnings_out = no_rep_list

                # Ensure Variants have had the refs removed.
                # if not hasattr(posedit, refseqgene_variant):
                if refseqgene_variant != '':
                    try:
                        refseqgene_variant = fn.valstr(hgvs_refseqgene_variant)
                    except:
                        fn.exceptPass()

                # Add single letter AA code to protein descriptions
                predicted_protein_variant_dict = {"tlr": str(predicted_protein_variant), "slr": ''}
                if predicted_protein_variant != '':
                    if not 'Non-coding :n.' in predicted_protein_variant:
                        try:
                            format_p = predicted_protein_variant
                            format_p = re.sub(r'\(LRG_.+?\)', '', format_p)
                            re_parse_protein = self.hp.parse_hgvs_variant(format_p)
                            re_parse_protein_singleAA = fn.single_letter_protein(re_parse_protein)
                            predicted_protein_variant_dict["slr"] = str(re_parse_protein_singleAA)
                        except hgvs.exceptions.HGVSParseError:
                            fn.exceptPass()
                    else:
                        predicted_protein_variant_dict["slr"] = str(predicted_protein_variant)

                # Populate the dictionary
                dict_out['submitted_variant'] = submitted
                dict_out['gene_symbol'] = gene_symbol
                dict_out['transcript_description'] = transcript_description
                dict_out['hgvs_transcript_variant'] = tx_variant
                dict_out['genome_context_intronic_sequence'] = genome_context_transcript_variant
                dict_out['refseqgene_context_intronic_sequence'] = RefSeqGene_context_transcript_variant
                dict_out['hgvs_refseqgene_variant'] = refseqgene_variant
                dict_out['hgvs_predicted_protein_consequence'] = predicted_protein_variant_dict
                dict_out['validation_warnings'] = warnings_out
                dict_out['hgvs_lrg_transcript_variant'] = lrg_transcript_variant
                dict_out['hgvs_lrg_variant'] = lrg_variant
                dict_out['alt_genomic_loci'] = alt_genomic_dicts
                dict_out['primary_assembly_loci'] = primary_genomic_dicts
                dict_out['reference_sequence_records'] = ''

                # Add links to reference_sequence_records
                ref_records = self.db.get_urls(dict_out)
                if ref_records != {}:
                    dict_out['reference_sequence_records'] = ref_records

                # Append to a list for return
                batch_out.append(dict_out)


            """
            Structure the output into dictionaries rather than a list with descriptive keys
            and a validation type flag
            """
            logger.trace("Populating output dictionary")
            # Create output dictionary
            validation_output = {'flag': None}

            # For gene outputs, i.e. those that hit transcripts
            # dotter = ''
            if my_variant.output_type_flag == 'gene':
                validation_output['flag'] = 'gene_variant'
                validation_error_counter = 0
                validation_obsolete_counter = 0
                for valid_v in batch_out:
                    if valid_v['validation_warnings'] == ['Validation error']:
                        validation_error_counter = validation_error_counter + 1
                        identification_key = 'Validation_Error_%s' % (str(validation_error_counter))
                    else:
                        obs_obs = False
                        for ob_rec in valid_v['validation_warnings']:
                            if 'obsolete' in ob_rec:
                                validation_obsolete_counter = validation_obsolete_counter + 1
                                obs_obs = True
                                break
                        if obs_obs is True:
                            identification_key = 'obsolete_record_%s' % (str(validation_obsolete_counter))
                        else:
                            identification_key = '%s' % (str(valid_v['hgvs_transcript_variant']))

                    # if identification_key not in validation_output.keys():
                    validation_output[identification_key] = valid_v
                    # else:
                    # dotter = dotter + ' '
                    # validation_output[identification_key + dotter] = valid_v

            # For warning only outputs
            # Should only ever be 1 output as an error or a warning of the following types
            # Gene symbol as reference sequence
            # Gene as transcript reference sequence
            if my_variant.output_type_flag == 'warning':
                validation_output['flag'] = 'warning'
                validation_error_counter = 0
                validation_warning_counter = 0
                if len(batch_out) == 0:
                    validation_output['flag'] = 'empty_result'
                for valid_v in batch_out:
                    if valid_v['validation_warnings'] == ['Validation error']:
                        validation_error_counter = validation_error_counter + 1
                        identification_key = 'validation_error_%s' % (str(validation_error_counter))
                    else:
                        validation_warning_counter = validation_warning_counter + 1
                        identification_key = 'validation_warning_%s' % (str(validation_warning_counter))
                    validation_output[identification_key] = valid_v

            # Intergenic variants
            validation_intergenic_counter = 0
            if my_variant.output_type_flag == 'intergenic':
                validation_output['flag'] = 'intergenic'
                for valid_v in batch_out:
                    validation_intergenic_counter = validation_intergenic_counter + 1
                    identification_key = 'Intergenic_Variant_%s' % (str(validation_intergenic_counter))

                    # Attempt to liftover between genome builds
                    # Note: pyliftover uses the UCSC liftOver tool.
                    # https://pypi.org/project/pyliftover/
                    genomic_position_info = valid_v['primary_assembly_loci']
                    for g_p_key in list(genomic_position_info.keys()):

                        # Identify the current build and hgvs_genomic descripsion
                        if re.match('hg', g_p_key):
                            # incoming_vcf = genomic_position_info[g_p_key]['vcf']
                            # set builds
                            if g_p_key == 'hg38':
                                build_to = 'hg19'
                                build_from = 'hg38'
                            if g_p_key == 'hg19':
                                build_to = 'hg38'
                                build_from = 'hg19'
                        elif re.match('grc', g_p_key):
                            # incoming_vcf = genomic_position_info[g_p_key]['vcf']
                            # set builds
                            if g_p_key == 'grch38':
                                build_to = 'GRCh37'
                                build_from = 'GRCh38'
                            if g_p_key == 'grch37':
                                build_to = 'GRCh38'
                                build_from = 'GRCh37'

                        # Liftover
                        lifted_response = lift_over(genomic_position_info[g_p_key]['hgvs_genomic_description'], build_from, build_to, hn, self.vm, self.vr, self.hdp, self.hp, reverse_normalizer, self.sf, evm)

                        # Sort the respomse into primary assembly and ALT
                        primary_assembly_loci = {}
                        alt_genomic_loci = []
                        for build_key, accession_dict in list(lifted_response.items()):
                            try:
                                accession_key = list(accession_dict.keys())[0]
                                if re.match('NC_', accession_dict[accession_key]['hgvs_genomic_description']):
                                    primary_assembly_loci[build_key.lower()] = accession_dict[accession_key]
                                else:
                                    alt_genomic_loci.append({build_key.lower(): accession_dict[accession_key]})

                            # KeyError if the dicts are empty
                            except KeyError:
                                continue

                        # Add the dictionaries from lifted response to the output
                        if primary_assembly_loci != {}:
                            valid_v['primary_assembly_loci'] = primary_assembly_loci
                        if alt_genomic_loci != []:
                            valid_v['alt_genomic_loci'] = alt_genomic_loci

                    # Finalise the output dictionary
                    validation_output[identification_key] = valid_v

            # Add error strings to validation output
            # '''
            metadata = {}
            logger.info("Variant successfully validated")
            logs = []
            logString = logger.getString()
            for l in logger.getString().split("\n"):
                logs.append(l)

            if os.environ.get("ADD_LOGS")=="True":
                metadata["logs"] = logs
            metadata["variant"] = batch_variant
            #metadata["assembly"] = selected_assembly
            #metadata["transcripts"] = select_transcripts
            #metadata['seqrepo_directory'] = self.seqrepoPath
            #metadata['uta_url'] = self.utaPath
            #metadata['py_liftover_directory'] = self.liftoverPath
            #metadata['variantvalidator_data_url'] = self.db.path
            #metadata['entrez_id'] = self.entrezID
            metadata['variantvalidator_version'] = self.version
            metadata['variantvalidator_hgvs_version'] = self.hgvsVersion
            metadata['uta_schema'] = self.utaSchema
            metadata['seqrepo_db'] = self.seqrepoVersion
            validation_output["metadata"] = metadata
            # '''
            # Measure time elapsed
            time_now = time.time()
            elapsed_time = time_now - start_time
            logger.debug('validation time = ' + str(elapsed_time))

            # return batch_out
            return validation_output

        # Bug catcher
        except KeyboardInterrupt:
            raise
        except BaseException as e:
            # Debug mode
            exc_type, exc_value, last_traceback = sys.exc_info()
            te = traceback.format_exc()
            # tr = ''.join(traceback.format_stack())
            tbk = [str(exc_type), str(exc_value), str(te)]
            er = '\n'.join(tbk)
            # raise VariantValidatorError('Validation error')
            # Return
            # return
            logger.critical(str(exc_type) + " " + str(exc_value))
            raise
            logger.debug(str(er))
