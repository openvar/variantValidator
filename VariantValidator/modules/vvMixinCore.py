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
from . import variant
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
                query = variant.Variant(queries)
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

                        # TODO: Need to check this as it's only being using outside of this loop as well as inside!
                        rec_var = ''

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

            for valid in by_order:
                if not valid.write:
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
                warnings = valid.warnings
                warnings = re.sub('del[GATC][GATC][GATC][GATC]+', 'del', warnings)
                warnings = re.sub('^: ', '', warnings)
                warnings = re.sub('::', ':', warnings)

                # Submitted variant
                submitted = valid.original

                # Genomic sequence variation
                genomic_variant = valid.genomic_g

                # genomic accession
                if genomic_variant != '':
                    hgvs_genomic_variant = self.hp.parse_hgvs_variant(genomic_variant)
                    genomic_variant = fn.valstr(hgvs_genomic_variant)
                    genomic_accession = hgvs_genomic_variant.ac
                else:
                    genomic_accession = ''

                # RefSeqGene variation
                refseqgene_variant = valid.genomic_r
                refseqgene_variant = refseqgene_variant.strip()
                if re.search('RefSeqGene', refseqgene_variant) or refseqgene_variant == '':
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
                tx_variant = valid.coding
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
                        if re.search('intronic variant', error):
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
                predicted_protein_variant = valid.protein
                if re.match('NP_', predicted_protein_variant):
                    rs_p, pred_prot_posedit = predicted_protein_variant.split(':')
                    lrg_p = self.db.get_lrgProteinID_from_RefSeqProteinID(rs_p)
                    if re.match('LRG', lrg_p):
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
                transcript_description = valid.description

                # Stashed variants
                # if valid.test_stash_tx_left:
                #     test_stash_tx_left = valid.test_stash_tx_left
                # if valid.test_stash_tx_right:
                #     test_stash_tx_right = valid.test_stash_tx_right

                # Multiple genomic variants
                # multi_gen_vars = []
                if tx_variant != '':
                    hgvs_coding = self.hp.parse_hgvs_variant(str(tx_variant))
                    # Gap gene black list
                    try:
                        gene_symbol = self.db.get_gene_symbol_from_transcriptID(hgvs_coding.ac)
                    except Exception:
                        fn.exceptPass()
                    else:
                        # If the gene symbol is not in the list, the value False will be returned
                        gap_compensation = vvChromosomes.gap_black_list(gene_symbol)

                    # Look for variants spanning introns
                    try:
                        hgvs_coding = hn.normalize(hgvs_coding)
                    except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                        error = str(e)
                        if re.search('boundary', str(error)) or re.search('spanning', str(error)):
                            gap_compensation = False
                        else:
                            pass
                    except hgvs.exceptions.HGVSError:
                        fn.exceptPass()

                    # Warn gap code status
                    logger.warning("gap_compensation_3 = " + str(gap_compensation))
                    multi_g = []
                    multi_list = []
                    mapping_options = self.hdp.get_tx_mapping_options(hgvs_coding.ac)
                    for alt_chr in mapping_options:
                        if (re.match('NC_', alt_chr[1]) or re.match('NT_', alt_chr[1]) or re.match('NW_',
                                                                                                   alt_chr[1])) and \
                                alt_chr[2] == alt_aln_method:
                            multi_list.append(alt_chr[1])

                    for alt_chr in multi_list:
                        try:
                            # Re set ori
                            ori = self.tx_exons(tx_ac=hgvs_coding.ac, alt_ac=alt_chr,
                                                   alt_aln_method=alt_aln_method)
                            orientation = int(ori[0]['alt_strand'])
                            hgvs_alt_genomic = self.myvm_t_to_g(hgvs_coding, alt_chr, no_norm_evm, hn)
                            # Set hgvs_genomic accordingly
                            hgvs_genomic = copy.deepcopy(hgvs_alt_genomic)

                            # genomic_possibilities
                            # 1. take the simple 3 pr normalized hgvs_genomic
                            # 2. Lock in hgvs_genomic at its most 5 prime position wrt genome
                            hgvs_genomic_possibilities = []

                            # Loop out gap code under these circumstances!
                            if gap_compensation is True:
                                logger.warning('g_to_t gap code 3 active')
                                rn_hgvs_genomic = reverse_normalizer.normalize(hgvs_alt_genomic)
                                hgvs_genomic_possibilities.append(rn_hgvs_genomic)
                                if orientation != -1:
                                    try:
                                        chromosome_normalized_hgvs_coding = reverse_normalizer.normalize(
                                            hgvs_coding)
                                    except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                                        chromosome_normalized_hgvs_coding = hgvs_coding
                                else:
                                    try:
                                        chromosome_normalized_hgvs_coding = hn.normalize(hgvs_coding)
                                    except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                                        error = str(e)
                                        chromosome_normalized_hgvs_coding = hgvs_coding

                                most_3pr_hgvs_genomic = self.myvm_t_to_g(chromosome_normalized_hgvs_coding,
                                                                            alt_chr,
                                                                            no_norm_evm, hn)
                                hgvs_genomic_possibilities.append(most_3pr_hgvs_genomic)

                                # First to the right
                                hgvs_stash = copy.deepcopy(hgvs_coding)
                                try:
                                    hgvs_stash = no_norm_evm.c_to_n(hgvs_stash)
                                except:
                                    fn.exceptPass()
                                try:
                                    stash_ac = hgvs_stash.ac
                                    stash_dict = vvHGVS.hard_right_hgvs2vcf(hgvs_stash, primary_assembly, hn, self.sf)
                                    stash_pos = int(stash_dict['pos'])
                                    stash_ref = stash_dict['ref']
                                    stash_alt = stash_dict['alt']
                                    # Generate an end position
                                    stash_end = str(stash_pos + len(stash_ref) - 1)
                                    # make a not real deletion insertion
                                    stash_hgvs_not_delins = self.hp.parse_hgvs_variant(
                                        stash_ac + ':' + hgvs_stash.type + '.' + str(
                                            stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)
                                    try:
                                        stash_hgvs_not_delins = no_norm_evm.n_to_c(stash_hgvs_not_delins)
                                    except:
                                        fn.exceptPass()
                                        # Store a tx copy for later use
                                    test_stash_tx_right = copy.deepcopy(stash_hgvs_not_delins)
                                    stash_genomic = self.myvm_t_to_g(test_stash_tx_right, hgvs_alt_genomic.ac,
                                                                        no_norm_evm, hn)
                                    # Stash the outputs if required
                                    # test variants = NC_000006.11:g.90403795G= (causes double identity)
                                    #                 NC_000002.11:g.73675227_73675228insCTC (? incorrect assumed insertion position)
                                    #                 NC_000003.11:g.14561629_14561630GC= NC_000003.11:g.14561629_14561630insG (Odd gap position)
                                    # if test_stash_tx_right.posedit.edit.type == 'identity' and stash_genomic.posedit.edit.type == 'identity':
                                    # pass
                                    if len(test_stash_tx_right.posedit.edit.ref) == ((
                                                                                             stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                                        stash_tx_right = test_stash_tx_right
                                        if hasattr(test_stash_tx_right.posedit.edit,
                                                   'alt') and test_stash_tx_right.posedit.edit.alt is not None:
                                            alt = test_stash_tx_right.posedit.edit.alt
                                        else:
                                            alt = ''
                                        if hasattr(stash_genomic.posedit.edit,
                                                   'alt') and stash_genomic.posedit.edit.alt is not None:
                                            g_alt = stash_genomic.posedit.edit.alt
                                        else:
                                            g_alt = ''
                                        if (len(alt) - (
                                                test_stash_tx_right.posedit.pos.end.base - test_stash_tx_right.posedit.pos.start.base) + 1) != (
                                                len(g_alt) - (
                                                stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                                            hgvs_genomic_possibilities.append(stash_genomic)
                                        else:
                                            hgvs_genomic_possibilities.append('')
                                    elif test_stash_tx_right.posedit.edit.type == 'identity':
                                        reform_ident = str(test_stash_tx_right).split(':')[0]
                                        reform_ident = reform_ident + ':c.' + str(test_stash_tx_right.posedit.pos) + 'del' + str(
                                            test_stash_tx_right.posedit.edit.ref)  # + 'ins' + str(test_stash_tx_right.posedit.edit.alt)
                                        hgvs_reform_ident = self.hp.parse_hgvs_variant(reform_ident)
                                        try:
                                            hn.normalize(hgvs_reform_ident)
                                        except hgvs.exceptions.HGVSError as e:
                                            error = str(e)
                                            if re.search('spanning the exon-intron boundary', error):
                                                stash_tx_right = test_stash_tx_right
                                                hgvs_genomic_possibilities.append('')
                                        else:
                                            stash_tx_right = test_stash_tx_right
                                            hgvs_genomic_possibilities.append(stash_genomic)
                                    else:
                                        try:
                                            hn.normalize(test_stash_tx_right)
                                        except hgvs.exceptions.HGVSUnsupportedOperationError:
                                            hgvs_genomic_possibilities.append('')
                                        else:
                                            stash_tx_right = test_stash_tx_right
                                            hgvs_genomic_possibilities.append(stash_genomic)
                                except hgvs.exceptions.HGVSError as e:
                                    fn.exceptPass()
                                except ValueError:
                                    fn.exceptPass()

                                # Then to the left
                                hgvs_stash = copy.deepcopy(hgvs_coding)
                                try:
                                    hgvs_stash = no_norm_evm.c_to_n(hgvs_stash)
                                except:
                                    fn.exceptPass()
                                try:
                                    stash_ac = hgvs_stash.ac
                                    stash_dict = vvHGVS.hard_left_hgvs2vcf(hgvs_stash, primary_assembly,
                                                                           reverse_normalizer, self.sf)
                                    stash_pos = int(stash_dict['pos'])
                                    stash_ref = stash_dict['ref']
                                    stash_alt = stash_dict['alt']
                                    # Generate an end position
                                    stash_end = str(stash_pos + len(stash_ref) - 1)
                                    # make a not real deletion insertion
                                    stash_hgvs_not_delins = self.hp.parse_hgvs_variant(
                                        stash_ac + ':' + hgvs_stash.type + '.' + str(
                                            stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)
                                    try:
                                        stash_hgvs_not_delins = no_norm_evm.n_to_c(stash_hgvs_not_delins)
                                    except:
                                        fn.exceptPass()
                                        # Store a tx copy for later use
                                    test_stash_tx_left = copy.deepcopy(stash_hgvs_not_delins)
                                    stash_genomic = self.myvm_t_to_g(test_stash_tx_left, hgvs_alt_genomic.ac,
                                                                        no_norm_evm, hn)
                                    # Stash the outputs if required
                                    # test variants = NC_000006.11:g.90403795G= (causes double identity)
                                    #                 NC_000002.11:g.73675227_73675228insCTC
                                    #                 NC_000003.11:g.14561629_14561630GC= NC_000003.11:g.14561629_14561630insG (Odd gap position)
                                    # if test_stash_tx_left.posedit.edit.type == 'identity' and stash_genomic.posedit.edit.type == 'identity':
                                    # pass
                                    if len(test_stash_tx_left.posedit.edit.ref) == ((
                                                                                            stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):  # len(stash_genomic.posedit.edit.ref):
                                        stash_tx_left = test_stash_tx_left
                                        if hasattr(test_stash_tx_left.posedit.edit,
                                                   'alt') and test_stash_tx_left.posedit.edit.alt is not None:
                                            alt = test_stash_tx_left.posedit.edit.alt
                                        else:
                                            alt = ''
                                        if hasattr(stash_genomic.posedit.edit,
                                                   'alt') and stash_genomic.posedit.edit.alt is not None:
                                            g_alt = stash_genomic.posedit.edit.alt
                                        else:
                                            g_alt = ''
                                        if (len(alt) - (
                                                test_stash_tx_left.posedit.pos.end.base - test_stash_tx_left.posedit.pos.start.base) + 1) != (
                                                len(g_alt) - (
                                                stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                                            hgvs_genomic_possibilities.append(stash_genomic)
                                        else:
                                            hgvs_genomic_possibilities.append('')
                                    elif test_stash_tx_left.posedit.edit.type == 'identity':
                                        reform_ident = str(test_stash_tx_left).split(':')[0]
                                        reform_ident = reform_ident + ':c.' + str(test_stash_tx_left.posedit.pos) + 'del' + str(
                                            test_stash_tx_left.posedit.edit.ref)  # + 'ins' + str(test_stash_tx_left.posedit.edit.alt)
                                        hgvs_reform_ident = self.hp.parse_hgvs_variant(reform_ident)
                                        try:
                                            hn.normalize(hgvs_reform_ident)
                                        except hgvs.exceptions.HGVSError as e:
                                            error = str(e)
                                            if re.search('spanning the exon-intron boundary', error):
                                                stash_tx_left = test_stash_tx_left
                                                hgvs_genomic_possibilities.append('')
                                        else:
                                            stash_tx_left = test_stash_tx_left
                                            hgvs_genomic_possibilities.append(stash_genomic)
                                    else:
                                        try:
                                            hn.normalize(test_stash_tx_left)
                                        except hgvs.exceptions.HGVSUnsupportedOperationError:
                                            hgvs_genomic_possibilities.append('')
                                        else:
                                            stash_tx_left = test_stash_tx_left
                                            hgvs_genomic_possibilities.append(stash_genomic)
                                except hgvs.exceptions.HGVSError as e:
                                    fn.exceptPass()
                                except ValueError:
                                    fn.exceptPass()

                                # direct mapping from reverse_normalized transcript insertions in the delins format
                                try:
                                    if hgvs_coding.posedit.edit.type == 'ins':
                                        most_5pr_hgvs_transcript_variant = copy.deepcopy(hgvs_coding)
                                        most_3pr_hgvs_transcript_variant = reverse_normalizer.normalize(hgvs_coding)
                                        try:
                                            n_3pr = self.vm.c_to_n(most_3pr_hgvs_transcript_variant)
                                            n_5pr = self.vm.c_to_n(most_5pr_hgvs_transcript_variant)
                                        except:
                                            n_3pr = most_3pr_hgvs_transcript_variant
                                            n_5pr = most_5pr_hgvs_transcript_variant
                                        # Make into a delins by adding the ref bases to the variant ref and alt
                                        pr3_ref = self.sf.fetch_seq(hgvs_coding.ac, n_3pr.posedit.pos.start.base - 1,
                                                               n_3pr.posedit.pos.end.base)
                                        pr5_ref = self.sf.fetch_seq(hgvs_coding.ac, n_5pr.posedit.pos.start.base - 1,
                                                               n_5pr.posedit.pos.end.base)
                                        most_3pr_hgvs_transcript_variant.posedit.edit.ref = pr3_ref
                                        most_5pr_hgvs_transcript_variant.posedit.edit.ref = pr5_ref
                                        most_3pr_hgvs_transcript_variant.posedit.edit.alt = pr3_ref[
                                                                                                0] + most_3pr_hgvs_transcript_variant.posedit.edit.alt + \
                                                                                            pr3_ref[1]
                                        most_5pr_hgvs_transcript_variant.posedit.edit.alt = pr5_ref[
                                                                                                0] + most_5pr_hgvs_transcript_variant.posedit.edit.alt + \
                                                                                            pr5_ref[1]
                                        # Map to the genome
                                        genomic_from_most_3pr_hgvs_transcript_variant = self.vm.t_to_g(
                                            most_3pr_hgvs_transcript_variant, hgvs_genomic.ac)
                                        genomic_from_most_5pr_hgvs_transcript_variant = self.vm.t_to_g(
                                            most_5pr_hgvs_transcript_variant, hgvs_genomic.ac)

                                        # Normalize - If the variant spans a gap it should then form a static genomic variant
                                        try:
                                            genomic_from_most_3pr_hgvs_transcript_variant = hn.normalize(
                                                genomic_from_most_3pr_hgvs_transcript_variant)
                                        except hgvs.exceptions.HGVSInvalidVariantError as e:
                                            error = str(e)
                                            if error == 'base start position must be <= end position':
                                                start = genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.start.base
                                                end = genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.end.base
                                                genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.start.base = end
                                                genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.end.base = start
                                                genomic_from_most_3pr_hgvs_transcript_variant = hn.normalize(
                                                    genomic_from_most_3pr_hgvs_transcript_variant)
                                        try:
                                            genomic_from_most_5pr_hgvs_transcript_variant = hn.normalize(
                                                genomic_from_most_5pr_hgvs_transcript_variant)
                                        except hgvs.exceptions.HGVSInvalidVariantError as e:
                                            error = str(e)
                                            if error == 'base start position must be <= end position':
                                                start = genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.start.base
                                                end = genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.end.base
                                                genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.start.base = end
                                                genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.end.base = start
                                                genomic_from_most_5pr_hgvs_transcript_variant = hn.normalize(
                                                    genomic_from_most_5pr_hgvs_transcript_variant)

                                        try:
                                            if genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.alt is None:
                                                genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.alt = ''
                                        except Exception as e:
                                            if str(e) == "'Dup' object has no attribute 'alt'":
                                                genomic_from_most_3pr_hgvs_transcript_variant_delins_from_dup = genomic_from_most_3pr_hgvs_transcript_variant.ac + ':' + genomic_from_most_3pr_hgvs_transcript_variant.type + '.' + str(
                                                    genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.start.base) + '_' + str(
                                                    genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.end.base) + 'del' + genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.ref + 'ins' + genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.ref + genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.ref
                                                genomic_from_most_3pr_hgvs_transcript_variant = self.hp.parse_hgvs_variant(
                                                    genomic_from_most_3pr_hgvs_transcript_variant_delins_from_dup)

                                        try:
                                            if most_3pr_hgvs_transcript_variant.posedit.edit.alt is None:
                                                most_3pr_hgvs_transcript_variant.posedit.edit.alt = ''
                                        except Exception as e:
                                            if str(e) == "'Dup' object has no attribute 'alt'":
                                                most_3pr_hgvs_transcript_variant_delins_from_dup = most_3pr_hgvs_transcript_variant.ac + ':' + most_3pr_hgvs_transcript_variant.type + '.' + str(
                                                    most_3pr_hgvs_transcript_variant.posedit.pos.start.base) + '_' + str(
                                                    most_3pr_hgvs_transcript_variant.posedit.pos.end.base) + 'del' + most_3pr_hgvs_transcript_variant.posedit.edit.ref + 'ins' + most_3pr_hgvs_transcript_variant.posedit.edit.ref + most_3pr_hgvs_transcript_variant.posedit.edit.ref
                                                most_3pr_hgvs_transcript_variant = self.hp.parse_hgvs_variant(
                                                    most_3pr_hgvs_transcript_variant_delins_from_dup)

                                        try:
                                            if genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.alt is None:
                                                genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.alt = ''
                                        except Exception as e:
                                            if str(e) == "'Dup' object has no attribute 'alt'":
                                                genomic_from_most_5pr_hgvs_transcript_variant_delins_from_dup = genomic_from_most_5pr_hgvs_transcript_variant.ac + ':' + genomic_from_most_5pr_hgvs_transcript_variant.type + '.' + str(
                                                    genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.start.base) + '_' + str(
                                                    genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.end.base) + 'del' + genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.ref + 'ins' + genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.ref + genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.ref
                                                genomic_from_most_5pr_hgvs_transcript_variant = self.hp.parse_hgvs_variant(
                                                    genomic_from_most_5pr_hgvs_transcript_variant_delins_from_dup)

                                        try:
                                            if most_5pr_hgvs_transcript_variant.posedit.edit.alt is None:
                                                most_5pr_hgvs_transcript_variant.posedit.edit.alt = ''
                                        except Exception as e:
                                            if str(e) == "'Dup' object has no attribute 'alt'":
                                                most_5pr_hgvs_transcript_variant_delins_from_dup = most_5pr_hgvs_transcript_variant.ac + ':' + most_5pr_hgvs_transcript_variant.type + '.' + str(
                                                    most_5pr_hgvs_transcript_variant.posedit.pos.start.base) + '_' + str(
                                                    most_5pr_hgvs_transcript_variant.posedit.pos.end.base) + 'del' + most_5pr_hgvs_transcript_variant.posedit.edit.ref + 'ins' + most_5pr_hgvs_transcript_variant.posedit.edit.ref + most_5pr_hgvs_transcript_variant.posedit.edit.ref
                                                most_5pr_hgvs_transcript_variant = self.hp.parse_hgvs_variant(
                                                    most_5pr_hgvs_transcript_variant_delins_from_dup)

                                        if len(
                                                genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.alt) < len(
                                            most_3pr_hgvs_transcript_variant.posedit.edit.alt):
                                            hgvs_genomic_possibilities.append(
                                                genomic_from_most_3pr_hgvs_transcript_variant)
                                        if len(
                                                genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.alt) < len(
                                            most_5pr_hgvs_transcript_variant.posedit.edit.alt):
                                            hgvs_genomic_possibilities.append(
                                                genomic_from_most_5pr_hgvs_transcript_variant)

                                except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                                    error = str(e)
                                    if re.match('Normalization of intronic variants is not supported',
                                                error) or re.match(
                                        'Unsupported normalization of variants spanning the exon-intron boundary',
                                        error):
                                        pass
                                    fn.exceptPass()

                                # Set variables for problem specific warnings
                                gapped_alignment_warning = ''
                                corrective_action_taken = ''
                                gapped_transcripts = ''
                                auto_info = ''

                                # Mark as not disparity detected
                                disparity_deletion_in = ['false', 'false']
                                # Loop through to see if a gap can be located
                                possibility_counter = 0
                                for possibility in hgvs_genomic_possibilities:
                                    possibility_counter = possibility_counter + 1
                                    # Loop out stash possibilities which will not spot gaps so are empty
                                    if possibility == '':
                                        continue

                                    # Use VCF generation code to push hgvs_genomic as for 5 prime as possible to uncover gaps
                                    hgvs_genomic_variant = possibility
                                    stored_hgvs_genomic_variant = copy.deepcopy(hgvs_genomic_variant)

                                    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
                                    try:
                                        reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(
                                            hgvs_genomic_variant)
                                    except hgvs.exceptions.HGVSError as e:
                                        # Strange error caused by gap in genomic
                                        error = str(e)
                                        if re.search('base start position must be <= end position', error):
                                            if hgvs_genomic.posedit.edit.type == 'delins':
                                                start = hgvs_genomic.posedit.pos.start.base
                                                end = hgvs_genomic.posedit.pos.end.base
                                                lhb = self.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                                                rhb = self.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                                                hgvs_genomic.posedit.edit.ref = lhb + rhb
                                                hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                                                hgvs_genomic.posedit.pos.start.base = end
                                                hgvs_genomic.posedit.pos.end.base = start
                                                reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(
                                                    hgvs_genomic)
                                            if hgvs_genomic.posedit.edit.type == 'del':
                                                start = hgvs_genomic.posedit.pos.start.base
                                                end = hgvs_genomic.posedit.pos.end.base
                                                lhb = self.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                                                rhb = self.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                                                hgvs_genomic.posedit.edit.ref = lhb + rhb
                                                hgvs_genomic.posedit.edit.alt = lhb + rhb
                                                hgvs_genomic.posedit.pos.start.base = end
                                                hgvs_genomic.posedit.pos.end.base = start
                                                reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(
                                                    hgvs_genomic)
                                        if re.search('insertion length must be 1', error):
                                            if hgvs_genomic.posedit.edit.type == 'ins':
                                                start = hgvs_genomic.posedit.pos.start.base
                                                end = hgvs_genomic.posedit.pos.end.base
                                                ref_bases = self.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, end)
                                                lhb = self.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                                                rhb = self.sf.fetch_seq(str(hgvs_genomic.ac), start, end)
                                                hgvs_genomic.posedit.edit.ref = lhb + rhb
                                                hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                                                reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(
                                                    hgvs_genomic)

                                    hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)
                                    # Store a copy for later use
                                    stored_hgvs_genomic_5pr = copy.deepcopy(hgvs_genomic_5pr)

                                    # Make VCF
                                    vcf_dict = vvHGVS.hgvs2vcf(reverse_normalized_hgvs_genomic, primary_assembly,
                                                               reverse_normalizer, self.sf)
                                    chr = vcf_dict['chr']
                                    pos = vcf_dict['pos']
                                    ref = vcf_dict['ref']
                                    alt = vcf_dict['alt']

                                    # Look for exonic gaps within transcript or chromosome
                                    no_normalized_c = 'false'  # Mark true to not produce an additional normalization of c.

                                    # Generate an end position
                                    end = str(int(pos) + len(ref) - 1)
                                    pos = str(pos)

                                    # Store a not real deletion insertion to test for gapping
                                    stored_hgvs_not_delins = self.hp.parse_hgvs_variant(str(
                                        hgvs_genomic_5pr.ac) + ':' + hgvs_genomic_5pr.type + '.' + pos + '_' + end + 'del' + ref + 'ins' + alt)
                                    v = [chr, pos, ref, alt]

                                    # Save a copy of current hgvs_coding
                                    try:
                                        saved_hgvs_coding = no_norm_evm.g_to_t(stored_hgvs_not_delins,
                                                                               hgvs_coding.ac)
                                    except Exception as e:
                                        if str(
                                                e) == 'start or end or both are beyond the bounds of transcript record':
                                            saved_hgvs_coding = hgvs_coding
                                            continue

                                    # Detect intronic variation using normalization
                                    intronic_variant = 'false'
                                    # Look for normalized variant options that do not match hgvs_coding
                                    if orientation == -1:
                                        # position genomic at its most 5 prime position
                                        try:
                                            query_genomic = reverse_normalizer.normalize(hgvs_genomic)
                                        except:
                                            query_genomic = hgvs_genomic
                                        # Map to the transcript ant test for movement
                                        try:
                                            hgvs_seek_var = evm.g_to_t(query_genomic, hgvs_coding.ac)
                                        except hgvs.exceptions.HGVSError as e:
                                            hgvs_seek_var = saved_hgvs_coding
                                        else:
                                            seek_var = fn.valstr(hgvs_seek_var)
                                            seek_ac = str(hgvs_seek_var.ac)
                                        if (
                                                hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                                                hgvs_coding.posedit.pos.start.base + hgvs_coding.posedit.pos.start.offset) and (
                                                hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                                                hgvs_coding.posedit.pos.end.base + hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
                                            pass
                                        else:
                                            hgvs_seek_var = saved_hgvs_coding

                                    elif orientation != -1:
                                        # position genomic at its most 3 prime position
                                        try:
                                            query_genomic = hn.normalize(hgvs_genomic)
                                        except:
                                            query_genomic = hgvs_genomic
                                        # Map to the transcript and test for movement
                                        try:
                                            hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
                                        except hgvs.exceptions.HGVSError as e:
                                            hgvs_seek_var = saved_hgvs_coding
                                        seek_var = fn.valstr(hgvs_seek_var)
                                        seek_ac = str(hgvs_seek_var.ac)
                                        if (
                                                hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                                                hgvs_coding.posedit.pos.start.base + hgvs_coding.posedit.pos.start.offset) and (
                                                hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                                                hgvs_coding.posedit.pos.end.base + hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
                                            pass
                                        else:
                                            hgvs_seek_var = saved_hgvs_coding

                                    try:
                                        intron_test = hn.normalize(hgvs_seek_var)
                                    except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                                        error = str(e)
                                        if re.match('Normalization of intronic variants is not supported',
                                                    error) or re.match(
                                            'Unsupported normalization of variants spanning the exon-intron boundary',
                                            error):
                                            if re.match(
                                                    'Unsupported normalization of variants spanning the exon-intron boundary',
                                                    error):
                                                intronic_variant = 'hard_fail'
                                            else:
                                                # Double check to see whether the variant is actually intronic?
                                                for exon in ori:
                                                    genomic_start = int(exon['alt_start_i'])
                                                    genomic_end = int(exon['alt_end_i'])
                                                    if (
                                                            hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                                                            hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                                                        intronic_variant = 'false'
                                                        break
                                                    else:
                                                        intronic_variant = 'true'

                                    if intronic_variant != 'hard_fail':
                                        if re.search(r'\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search(r'\d+\-',
                                                                                                           str(
                                                                                                               hgvs_seek_var.posedit.pos)) or re.search(
                                            r'\*\d+\+', str(
                                                hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+\-', str(
                                            hgvs_seek_var.posedit.pos)):
                                            # Double check to see whether the variant is actually intronic?
                                            for exon in ori:
                                                genomic_start = int(exon['alt_start_i'])
                                                genomic_end = int(exon['alt_end_i'])
                                                if (
                                                        hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                                                        hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                                                    intronic_variant = 'false'
                                                    break
                                                else:
                                                    intronic_variant = 'true'

                                    if intronic_variant != 'true':
                                        # Flag RefSeqGene for ammendment
                                        # amend_RefSeqGene = 'false'
                                        # Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping transcript lengths
                                        if stored_hgvs_not_delins != '':
                                            # Refresh hgvs_not_delins from stored_hgvs_not_delins
                                            hgvs_not_delins = copy.deepcopy(stored_hgvs_not_delins)
                                            # This test will only occur in dup of single base, insertion or substitution
                                            if not re.search('_', str(hgvs_not_delins.posedit.pos)):
                                                if re.search('dup',
                                                             hgvs_genomic_5pr.posedit.edit.type) or re.search(
                                                    'ins', hgvs_genomic_5pr.posedit.edit.type):
                                                    # For gap in chr, map to t. - but becaouse we have pushed to 5 prime by norm, add 1 to end pos
                                                    plussed_hgvs_not_delins = copy.deepcopy(hgvs_not_delins)
                                                    plussed_hgvs_not_delins.posedit.pos.end.base = plussed_hgvs_not_delins.posedit.pos.end.base + 1
                                                    plussed_hgvs_not_delins.posedit.edit.ref = ''
                                                    transcript_variant = no_norm_evm.g_to_t(plussed_hgvs_not_delins,
                                                                                            str(
                                                                                                saved_hgvs_coding.ac))
                                                    if ((
                                                            transcript_variant.posedit.pos.end.base - transcript_variant.posedit.pos.start.base) > (
                                                            hgvs_genomic_5pr.posedit.pos.end.base - hgvs_genomic_5pr.posedit.pos.start.base)):
                                                        if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                                                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                            start = hgvs_not_delins.posedit.pos.start.base - 1
                                                            end = hgvs_not_delins.posedit.pos.end.base
                                                            ref_bases = self.sf.fetch_seq(str(hgvs_not_delins.ac), start,
                                                                                     end)
                                                            hgvs_not_delins.posedit.edit.ref = ref_bases
                                                            hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                               :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                                     1:] + ref_bases[
                                                                                                           1:]
                                                        elif re.search('ins', str(
                                                                hgvs_genomic_5pr.posedit.edit)) and re.search('del',
                                                                                                              str(
                                                                                                                  hgvs_genomic_5pr.posedit.edit)):
                                                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                        elif re.search('ins', str(
                                                                hgvs_genomic_5pr.posedit.edit)) and not re.search(
                                                            'del',
                                                            str(
                                                                hgvs_genomic_5pr.posedit.edit)):
                                                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                            start = hgvs_not_delins.posedit.pos.start.base - 1
                                                            end = hgvs_not_delins.posedit.pos.end.base
                                                            ref_bases = self.sf.fetch_seq(str(hgvs_not_delins.ac), start,
                                                                                     end)
                                                            hgvs_not_delins.posedit.edit.ref = ref_bases
                                                            hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                               :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                                     1:] + ref_bases[
                                                                                                           1:]
                                                    else:
                                                        if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                                                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                            start = hgvs_not_delins.posedit.pos.start.base - 1
                                                            end = hgvs_not_delins.posedit.pos.end.base
                                                            ref_bases = self.sf.fetch_seq(str(hgvs_not_delins.ac), start,
                                                                                     end)
                                                            hgvs_not_delins.posedit.edit.ref = ref_bases
                                                            hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                               :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                                     1:] + ref_bases[
                                                                                                           1:]
                                                        elif re.search('ins', str(
                                                                hgvs_genomic_5pr.posedit.edit)) and re.search('del',
                                                                                                              str(
                                                                                                                  hgvs_genomic_5pr.posedit.edit)):
                                                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                        elif re.search('ins', str(
                                                                hgvs_genomic_5pr.posedit.edit)) and not re.search(
                                                            'del',
                                                            str(
                                                                hgvs_genomic_5pr.posedit.edit)):
                                                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                            start = hgvs_not_delins.posedit.pos.start.base - 1
                                                            end = hgvs_not_delins.posedit.pos.end.base
                                                            ref_bases = self.sf.fetch_seq(str(hgvs_not_delins.ac), start,
                                                                                     end)
                                                            hgvs_not_delins.posedit.edit.ref = ref_bases
                                                            hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                               :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                                     1:] + ref_bases[
                                                                                                           1:]
                                                else:
                                                    pass
                                            else:
                                                pass
                                            tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                    saved_hgvs_coding.ac)
                                            # Create normalized version of tx_hgvs_not_delins
                                            rn_tx_hgvs_not_delins = copy.deepcopy(tx_hgvs_not_delins)
                                            # Check for +1 base and adjust
                                            if re.search(r'\+',
                                                         str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search(
                                                r'\+',
                                                str(
                                                    rn_tx_hgvs_not_delins.posedit.pos.start)):
                                                # Remove offsetting to span the gap
                                                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                                rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                                                rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
                                                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                                try:
                                                    rn_tx_hgvs_not_delins.posedit.edit.alt = ''
                                                except:
                                                    fn.exceptPass()

                                            elif re.search(r'\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                                                # move tx end base to next available non-offset base
                                                rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
                                                rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                                                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                                if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                                    test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                                else:
                                                    test_tx_var = rn_tx_hgvs_not_delins
                                                # re-make genomic and tx
                                                hgvs_not_delins = self.myvm_t_to_g(test_tx_var, alt_chr,
                                                                        no_norm_evm, hn)
                                                rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                           str(
                                                                                               saved_hgvs_coding.ac))
                                            elif re.search(r'\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                                                # move tx start base to previous available non-offset base
                                                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                                if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                                    test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                                else:
                                                    test_tx_var = rn_tx_hgvs_not_delins
                                                # re-make genomic and tx
                                                hgvs_not_delins = self.myvm_t_to_g(test_tx_var, alt_chr,
                                                                        no_norm_evm, hn)
                                                rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                           str(
                                                                                               saved_hgvs_coding.ac))
                                                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                            #                                                 else:
                                            #                                                     pass

                                            # Check for -ve base and adjust
                                            elif re.search(r'\-',
                                                           str(
                                                               rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search(
                                                r'\-',
                                                str(
                                                    rn_tx_hgvs_not_delins.posedit.pos.start)):
                                                # Remove offsetting to span the gap
                                                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                                rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                                                rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
                                                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                                try:
                                                    rn_tx_hgvs_not_delins.posedit.edit.alt = ''
                                                except:
                                                    fn.exceptPass()
                                            elif re.search(r'\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                                                # move tx end base back to next available non-offset base
                                                rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                                                # Delete the ref
                                                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                                # Add the additional base to the ALT
                                                start = rn_tx_hgvs_not_delins.posedit.pos.end.base - 1
                                                end = rn_tx_hgvs_not_delins.posedit.pos.end.base
                                                ref_bases = self.sf.fetch_seq(str(tx_hgvs_not_delins.ac), start, end)
                                                rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + ref_bases
                                                if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                                    test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                                else:
                                                    test_tx_var = rn_tx_hgvs_not_delins
                                                # re-make genomic and tx
                                                hgvs_not_delins = self.myvm_t_to_g(test_tx_var, alt_chr,
                                                                        no_norm_evm, hn)
                                                rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                           str(
                                                                                               saved_hgvs_coding.ac))
                                            elif re.search(r'\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                                                # move tx start base to previous available non-offset base
                                                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                                rn_tx_hgvs_not_delins.posedit.pos.start.base = rn_tx_hgvs_not_delins.posedit.pos.start.base - 1
                                                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                                if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                                    test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                                else:
                                                    test_tx_var = rn_tx_hgvs_not_delins
                                                # re-make genomic and tx
                                                hgvs_not_delins = self.myvm_t_to_g(test_tx_var, alt_chr,
                                                                        no_norm_evm, hn)
                                                rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                           str(
                                                                                               saved_hgvs_coding.ac))
                                                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                            else:
                                                fn.exceptPass()

                                            # Logic
                                            if len(hgvs_not_delins.posedit.edit.ref) < len(
                                                    rn_tx_hgvs_not_delins.posedit.edit.ref):
                                                gap_length = len(rn_tx_hgvs_not_delins.posedit.edit.ref) - len(
                                                    hgvs_not_delins.posedit.edit.ref)
                                                disparity_deletion_in = ['chromosome', gap_length]
                                            elif len(hgvs_not_delins.posedit.edit.ref) > len(
                                                    rn_tx_hgvs_not_delins.posedit.edit.ref):
                                                gap_length = len(hgvs_not_delins.posedit.edit.ref) - len(
                                                    rn_tx_hgvs_not_delins.posedit.edit.ref)
                                                disparity_deletion_in = ['transcript', gap_length]
                                            else:
                                                re_capture_tx_variant = []
                                                for possibility in hgvs_genomic_possibilities:
                                                    if possibility == '':
                                                        continue
                                                    hgvs_t_possibility = self.vm.g_to_t(possibility, hgvs_coding.ac)
                                                    if hgvs_t_possibility.posedit.edit.type == 'ins':
                                                        try:
                                                            hgvs_t_possibility = self.vm.c_to_n(hgvs_t_possibility)
                                                        except:
                                                            continue
                                                        if hgvs_t_possibility.posedit.pos.start.offset != 0 or hgvs_t_possibility.posedit.pos.end.offset != 0:
                                                            continue
                                                        ins_ref = self.sf.fetch_seq(hgvs_t_possibility.ac,
                                                                               hgvs_t_possibility.posedit.pos.start.base - 1,
                                                                               hgvs_t_possibility.posedit.pos.start.base + 1)
                                                        try:
                                                            hgvs_t_possibility = self.vm.n_to_c(hgvs_t_possibility)
                                                        except:
                                                            continue
                                                        hgvs_t_possibility.posedit.edit.ref = ins_ref
                                                        hgvs_t_possibility.posedit.edit.alt = ins_ref[
                                                                                                  0] + hgvs_t_possibility.posedit.edit.alt + \
                                                                                              ins_ref[1]
                                                    if possibility.posedit.edit.type == 'ins':
                                                        ins_ref = self.sf.fetch_seq(possibility.ac,
                                                                               possibility.posedit.pos.start.base - 1,
                                                                               possibility.posedit.pos.end.base)
                                                        possibility.posedit.edit.ref = ins_ref
                                                        possibility.posedit.edit.alt = ins_ref[
                                                                                           0] + possibility.posedit.edit.alt + \
                                                                                       ins_ref[1]
                                                    if len(hgvs_t_possibility.posedit.edit.ref) < len(
                                                            possibility.posedit.edit.ref):
                                                        gap_length = len(possibility.posedit.edit.ref) - len(
                                                            hgvs_t_possibility.posedit.edit.ref)
                                                        re_capture_tx_variant = ['transcript', gap_length,
                                                                                 hgvs_t_possibility]
                                                        hgvs_not_delins = possibility
                                                        hgvs_genomic_5pr = possibility
                                                        break

                                                if re_capture_tx_variant != []:
                                                    try:
                                                        tx_hgvs_not_delins = self.vm.c_to_n(re_capture_tx_variant[2])
                                                    except:
                                                        tx_hgvs_not_delins = re_capture_tx_variant[2]
                                                    disparity_deletion_in = re_capture_tx_variant[0:-1]
                                                else:
                                                    pass

                                        # Final sanity checks
                                        try:
                                            self.vm.g_to_t(hgvs_not_delins, tx_hgvs_not_delins.ac)
                                        except Exception as e:
                                            if str(
                                                    e) == 'start or end or both are beyond the bounds of transcript record':
                                                continue
                                        try:
                                            hn.normalize(tx_hgvs_not_delins)
                                        except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                                            error = str(e)

                                            if re.match('Normalization of intronic variants is not supported',
                                                        error) or re.match(
                                                'Unsupported normalization of variants spanning the exon-intron boundary',
                                                error):
                                                if re.match(
                                                        'Unsupported normalization of variants spanning the exon-intron boundary',
                                                        error):
                                                    continue
                                                elif re.match('Normalization of intronic variants is not supported',
                                                              error):
                                                    # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                                                    disparity_deletion_in = ['transcript', 'Requires Analysis']

                                        # Recreate hgvs_genomic
                                        if disparity_deletion_in[0] == 'transcript':
                                            hgvs_genomic = hgvs_not_delins

                                        # Find oddly placed gaps where the tx variant is encompassed in the gap
                                        if disparity_deletion_in[0] == 'false' and (
                                                possibility_counter == 3 or possibility_counter == 4):
                                            rg = reverse_normalizer.normalize(hgvs_not_delins)
                                            rtx = self.vm.g_to_t(rg, tx_hgvs_not_delins.ac)
                                            fg = hn.normalize(hgvs_not_delins)
                                            ftx = self.vm.g_to_t(fg, tx_hgvs_not_delins.ac)
                                            if (
                                                    rtx.posedit.pos.start.offset == 0 and rtx.posedit.pos.end.offset == 0) and (
                                                    ftx.posedit.pos.start.offset != 0 and ftx.posedit.pos.end.offset != 0):
                                                exons = self.hdp.get_tx_exons(ftx.ac, hgvs_not_delins.ac, alt_aln_method)
                                                exonic = False
                                                for ex_test in exons:
                                                    if ftx.posedit.pos.start.base in range(ex_test[6], ex_test[
                                                        7]) and ftx.posedit.pos.end.base in range(ex_test[6],
                                                                                                  ex_test[7]):
                                                        exonic = True
                                                if exonic is True:
                                                    hgvs_not_delins = fg
                                                    hgvs_genomic = fg
                                                    hgvs_genomic_5pr = fg
                                                    try:
                                                        tx_hgvs_not_delins = self.vm.c_to_n(ftx)
                                                    except Exception:
                                                        tx_hgvs_not_delins = ftx
                                                    disparity_deletion_in = ['transcript', 'Requires Analysis']

                                        # Pre-processing of tx_hgvs_not_delins
                                        try:
                                            if tx_hgvs_not_delins.posedit.edit.alt is None:
                                                tx_hgvs_not_delins.posedit.edit.alt = ''
                                        except Exception as e:
                                            if str(e) == "'Dup' object has no attribute 'alt'":
                                                tx_hgvs_not_delins_delins_from_dup = tx_hgvs_not_delins.ac + ':' + tx_hgvs_not_delins.type + '.' + str(
                                                    tx_hgvs_not_delins.posedit.pos.start) + '_' + str(
                                                    tx_hgvs_not_delins.posedit.pos.end) + 'del' + tx_hgvs_not_delins.posedit.edit.ref + 'ins' + tx_hgvs_not_delins.posedit.edit.ref + tx_hgvs_not_delins.posedit.edit.ref
                                                tx_hgvs_not_delins = self.hp.parse_hgvs_variant(
                                                    tx_hgvs_not_delins_delins_from_dup)

                                        if disparity_deletion_in[0] == 'transcript':
                                            # amend_RefSeqGene = 'true'
                                            # ANY VARIANT WHOLLY WITHIN THE GAP
                                            if (re.search(r'\+',
                                                          str(tx_hgvs_not_delins.posedit.pos.start)) or re.search(
                                                r'\-', str(tx_hgvs_not_delins.posedit.pos.start))) and (
                                                    re.search(r'\+',
                                                              str(tx_hgvs_not_delins.posedit.pos.end)) or re.search(
                                                r'\-', str(tx_hgvs_not_delins.posedit.pos.end))):
                                                gapped_transcripts = gapped_transcripts + ' ' + str(
                                                    tx_hgvs_not_delins.ac)

                                                # Copy the current variant
                                                tx_gap_fill_variant = copy.deepcopy(tx_hgvs_not_delins)
                                                try:
                                                    if tx_gap_fill_variant.posedit.edit.alt is None:
                                                        tx_gap_fill_variant.posedit.edit.alt = ''
                                                except Exception as e:
                                                    if str(e) == "'Dup' object has no attribute 'alt'":
                                                        tx_gap_fill_variant_delins_from_dup = tx_gap_fill_variant.ac + ':' + tx_gap_fill_variant.type + '.' + str(
                                                            tx_gap_fill_variant.posedit.pos.start) + '_' + str(
                                                            tx_gap_fill_variant.posedit.pos.end) + 'del' + tx_gap_fill_variant.posedit.edit.ref + 'ins' + tx_gap_fill_variant.posedit.edit.ref + tx_gap_fill_variant.posedit.edit.ref
                                                        tx_gap_fill_variant = self.hp.parse_hgvs_variant(
                                                            tx_gap_fill_variant_delins_from_dup)

                                                # Identify which half of the NOT-intron the start position of the variant is in
                                                if re.search(r'\-', str(tx_gap_fill_variant.posedit.pos.start)):
                                                    tx_gap_fill_variant.posedit.pos.start.base = tx_gap_fill_variant.posedit.pos.start.base - 1
                                                    tx_gap_fill_variant.posedit.pos.start.offset = int(
                                                        '0')  # int('+1')
                                                    tx_gap_fill_variant.posedit.pos.end.offset = int(
                                                        '0')  # int('-1')
                                                    tx_gap_fill_variant.posedit.edit.alt = ''
                                                    tx_gap_fill_variant.posedit.edit.ref = ''
                                                elif re.search(r'\+', str(tx_gap_fill_variant.posedit.pos.start)):
                                                    tx_gap_fill_variant.posedit.pos.start.offset = int(
                                                        '0')  # int('+1')
                                                    tx_gap_fill_variant.posedit.pos.end.base = tx_gap_fill_variant.posedit.pos.end.base + 1
                                                    tx_gap_fill_variant.posedit.pos.end.offset = int(
                                                        '0')  # int('-1')
                                                    tx_gap_fill_variant.posedit.edit.alt = ''
                                                    tx_gap_fill_variant.posedit.edit.ref = ''

                                                try:
                                                    tx_gap_fill_variant = self.vm.n_to_c(tx_gap_fill_variant)
                                                except:
                                                    fn.exceptPass()
                                                genomic_gap_fill_variant = self.vm.t_to_g(tx_gap_fill_variant,
                                                                                     reverse_normalized_hgvs_genomic.ac)
                                                genomic_gap_fill_variant.posedit.edit.alt = genomic_gap_fill_variant.posedit.edit.ref

                                                try:
                                                    c_tx_hgvs_not_delins = self.vm.n_to_c(tx_hgvs_not_delins)
                                                except Exception:
                                                    c_tx_hgvs_not_delins = copy.copy(tx_hgvs_not_delins)
                                                genomic_gap_fill_variant_alt = self.vm.t_to_g(c_tx_hgvs_not_delins,
                                                                                         hgvs_genomic_5pr.ac)

                                                # Ensure an ALT exists
                                                try:
                                                    if genomic_gap_fill_variant_alt.posedit.edit.alt is None:
                                                        genomic_gap_fill_variant_alt.posedit.edit.alt = 'X'
                                                except Exception as e:
                                                    if str(e) == "'Dup' object has no attribute 'alt'":
                                                        genomic_gap_fill_variant_delins_from_dup = genomic_gap_fill_variant.ac + ':' + genomic_gap_fill_variant.type + '.' + str(
                                                            genomic_gap_fill_variant.posedit.pos.start.base) + '_' + str(
                                                            genomic_gap_fill_variant.posedit.pos.end.base) + 'del' + genomic_gap_fill_variant.posedit.edit.ref + 'ins' + genomic_gap_fill_variant.posedit.edit.ref + genomic_gap_fill_variant.posedit.edit.ref
                                                        genomic_gap_fill_variant = self.hp.parse_hgvs_variant(
                                                            genomic_gap_fill_variant_delins_from_dup)
                                                        genomic_gap_fill_variant_alt_delins_from_dup = genomic_gap_fill_variant_alt.ac + ':' + genomic_gap_fill_variant_alt.type + '.' + str(
                                                            genomic_gap_fill_variant_alt.posedit.pos.start.base) + '_' + str(
                                                            genomic_gap_fill_variant_alt.posedit.pos.end.base) + 'del' + genomic_gap_fill_variant_alt.posedit.edit.ref + 'ins' + genomic_gap_fill_variant_alt.posedit.edit.ref + genomic_gap_fill_variant_alt.posedit.edit.ref
                                                        genomic_gap_fill_variant_alt = self.hp.parse_hgvs_variant(
                                                            genomic_gap_fill_variant_alt_delins_from_dup)

                                                # Correct insertion alts
                                                if genomic_gap_fill_variant_alt.posedit.edit.type == 'ins':
                                                    append_ref = self.sf.fetch_seq(genomic_gap_fill_variant_alt.ac,
                                                                              genomic_gap_fill_variant_alt.posedit.pos.start.base - 1,
                                                                              genomic_gap_fill_variant_alt.posedit.pos.end.base)
                                                    genomic_gap_fill_variant_alt.posedit.edit.alt = append_ref[
                                                                                                        0] + genomic_gap_fill_variant_alt.posedit.edit.alt + \
                                                                                                    append_ref[1]

                                                # Split the reference and replacing alt sequence into a dictionary
                                                reference_bases = list(genomic_gap_fill_variant.posedit.edit.ref)
                                                if genomic_gap_fill_variant_alt.posedit.edit.alt is not None:
                                                    alternate_bases = list(
                                                        genomic_gap_fill_variant_alt.posedit.edit.alt)
                                                else:
                                                    # Deletions with no ins
                                                    pre_alternate_bases = list(
                                                        genomic_gap_fill_variant_alt.posedit.edit.ref)
                                                    alternate_bases = []
                                                    for base in pre_alternate_bases:
                                                        alternate_bases.append('X')

                                                # Create the dictionaries
                                                ref_start = genomic_gap_fill_variant.posedit.pos.start.base
                                                alt_start = genomic_gap_fill_variant_alt.posedit.pos.start.base
                                                ref_base_dict = {}
                                                for base in reference_bases:
                                                    ref_base_dict[ref_start] = str(base)
                                                    ref_start = ref_start + 1

                                                alt_base_dict = {}

                                                # NEED TO SEARCH FOR RANGE = and replace with interval_range
                                                # Need to search for int and replace with integer

                                                # Note, all variants will be forced into the format delete insert
                                                # Deleted bases in the ALT will be substituted for X
                                                for integer in range(
                                                        genomic_gap_fill_variant_alt.posedit.pos.start.base,
                                                        genomic_gap_fill_variant_alt.posedit.pos.end.base + 1, 1):
                                                    if integer == alt_start:
                                                        alt_base_dict[integer] = str(''.join(alternate_bases))
                                                    else:
                                                        alt_base_dict[integer] = 'X'

                                                # Generate the alt sequence
                                                alternate_sequence_bases = []
                                                for integer in range(
                                                        genomic_gap_fill_variant.posedit.pos.start.base,
                                                        genomic_gap_fill_variant.posedit.pos.end.base + 1,
                                                        1):
                                                    if integer in list(alt_base_dict.keys()):
                                                        alternate_sequence_bases.append(alt_base_dict[integer])
                                                    else:
                                                        alternate_sequence_bases.append(ref_base_dict[integer])
                                                alternate_sequence = ''.join(alternate_sequence_bases)
                                                alternate_sequence = alternate_sequence.replace('X', '')

                                                # Add the new alt to the gap fill variant and generate transcript variant
                                                genomic_gap_fill_variant.posedit.edit.alt = alternate_sequence
                                                hgvs_refreshed_variant = self.vm.g_to_t(genomic_gap_fill_variant,
                                                                                   tx_gap_fill_variant.ac)

                                                # Set warning
                                                gap_size = str(len(genomic_gap_fill_variant.posedit.edit.ref) - 2)
                                                disparity_deletion_in[1] = [gap_size]
                                                auto_info = auto_info + str(
                                                    stored_hgvs_not_delins.ac) + ':g.' + str(
                                                    stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + gap_size + ' genomic base(s) that fail to align to transcript ' + str(
                                                    tx_hgvs_not_delins.ac)
                                                non_valid_caution = 'true'

                                                # Alignment position
                                                for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                                if re.match('NM_', str(for_location_c)):
                                                    for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                                if re.match(r'\-', str(for_location_c.posedit.pos.start.offset)):
                                                    gps = for_location_c.posedit.pos.start.base - 1
                                                    gpe = for_location_c.posedit.pos.start.base
                                                else:
                                                    gps = for_location_c.posedit.pos.start.base
                                                    gpe = for_location_c.posedit.pos.start.base + 1
                                                gap_position = ' between positions c.' + str(gps) + '_' + str(
                                                    gpe) + '\n'
                                                auto_info = auto_info + '%s' % (gap_position)

                                            else:
                                                if tx_hgvs_not_delins.posedit.pos.start.offset == 0 and tx_hgvs_not_delins.posedit.pos.end.offset == 0:
                                                    # In this instance, we have identified a transcript gap but the n. version of
                                                    # the transcript variant but do not have a position which actually hits the gap,
                                                    # so the variant likely spans the gap, and is not picked up by an offset.
                                                    try:
                                                        c1 = self.vm.n_to_c(tx_hgvs_not_delins)
                                                    except:
                                                        c1 = tx_hgvs_not_delins
                                                    g1 = self.nr_vm.t_to_g(c1, hgvs_genomic.ac)
                                                    g3 = self.nr_vm.t_to_g(c1, hgvs_genomic.ac)
                                                    g2 = self.vm.t_to_g(c1, hgvs_genomic.ac)
                                                    ng2 = hn.normalize(g2)
                                                    g3.posedit.pos.end.base = g3.posedit.pos.start.base + (
                                                            len(g3.posedit.edit.ref) - 1)
                                                    try:
                                                        c2 = self.vm.g_to_t(g3, c1.ac)
                                                        if c2.posedit.pos.start.offset == 0 and c2.posedit.pos.end.offset == 0:
                                                            pass
                                                        else:
                                                            tx_hgvs_not_delins = c2
                                                            try:
                                                                tx_hgvs_not_delins = self.vm.c_to_n(tx_hgvs_not_delins)
                                                            except hgvs.exceptions.HGVSError:
                                                                fn.exceptPass()
                                                    except hgvs.exceptions.HGVSInvalidVariantError:
                                                        fn.exceptPass()

                                                if re.search(r'\+', str(
                                                        tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                                                    r'\+',
                                                    str(
                                                        tx_hgvs_not_delins.posedit.pos.end)):
                                                    auto_info = auto_info + str(
                                                        stored_hgvs_not_delins.ac) + ':g.' + str(
                                                        stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                                                        disparity_deletion_in[
                                                            1]) + ' genomic base(s) that fail to align to transcript ' + str(
                                                        tx_hgvs_not_delins.ac)
                                                    non_valid_caution = 'true'
                                                    try:
                                                        c2 = self.vm.n_to_c(tx_hgvs_not_delins)
                                                    except:
                                                        c2 = tx_hgvs_not_delins
                                                    c1 = copy.deepcopy(c2)
                                                    c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                                                    c1.posedit.pos.start.offset = 0
                                                    c1.posedit.pos.end = c2.posedit.pos.start
                                                    c1.posedit.edit.ref = ''
                                                    c1.posedit.edit.alt = ''
                                                    if orientation != -1:
                                                        g1 = self.vm.t_to_g(c1, hgvs_genomic.ac)
                                                        g2 = self.vm.t_to_g(c2, hgvs_genomic.ac)
                                                        g1.posedit.edit.alt = g1.posedit.edit.ref
                                                    else:
                                                        g1 = self.vm.t_to_g(c2, hgvs_genomic.ac)
                                                        g2 = self.vm.t_to_g(c1, hgvs_genomic.ac)
                                                        g2.posedit.edit.alt = g2.posedit.edit.ref
                                                    reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                                    alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                                    g3 = copy.deepcopy(g1)
                                                    g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                                    g3.posedit.edit.ref = reference
                                                    g3.posedit.edit.alt = alternate
                                                    c3 = self.vm.g_to_t(g3, c1.ac)
                                                    hgvs_refreshed_variant = c3
                                                    # Alignment position
                                                    for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                                    if re.match('NM_', str(for_location_c)):
                                                        for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                                        gps = for_location_c.posedit.pos.start.base
                                                        gpe = for_location_c.posedit.pos.start.base + 1
                                                    gap_position = ' between positions c.' + str(gps) + '_' + str(
                                                        gpe) + '\n'
                                                    # Warn update
                                                    auto_info = auto_info + '%s' % (gap_position)
                                                elif re.search(r'\+', str(
                                                        tx_hgvs_not_delins.posedit.pos.end)) and not re.search(r'\+',
                                                                                                               str(
                                                                                                                   tx_hgvs_not_delins.posedit.pos.start)):
                                                    auto_info = auto_info + 'Genome position ' + str(
                                                        stored_hgvs_not_delins.ac) + ':g.' + str(
                                                        stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(
                                                        disparity_deletion_in[1]) + '-bp gap in transcript ' + str(
                                                        tx_hgvs_not_delins.ac)
                                                    gapped_transcripts = gapped_transcripts + ' ' + str(
                                                        tx_hgvs_not_delins.ac)
                                                    non_valid_caution = 'true'
                                                    try:
                                                        c1 = self.vm.n_to_c(tx_hgvs_not_delins)
                                                    except:
                                                        c1 = tx_hgvs_not_delins
                                                    c2 = copy.deepcopy(c1)
                                                    c2.posedit.pos.start = c1.posedit.pos.end
                                                    c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                                                    c2.posedit.pos.end.offset = 0
                                                    c2.posedit.edit.ref = ''
                                                    c2.posedit.edit.alt = ''
                                                    if orientation != -1:
                                                        g1 = self.vm.t_to_g(c1, hgvs_genomic.ac)
                                                        g2 = self.vm.t_to_g(c2, hgvs_genomic.ac)
                                                        g2.posedit.edit.alt = g2.posedit.edit.ref
                                                    else:
                                                        g1 = self.vm.t_to_g(c2, hgvs_genomic.ac)
                                                        g2 = self.vm.t_to_g(c1, hgvs_genomic.ac)
                                                        g1.posedit.edit.alt = g1.posedit.edit.ref
                                                    reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                                    alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                                    g3 = copy.deepcopy(g1)
                                                    g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                                    g3.posedit.edit.ref = reference
                                                    g3.posedit.edit.alt = alternate
                                                    c3 = self.vm.g_to_t(g3, c1.ac)
                                                    hgvs_refreshed_variant = c3
                                                    # Alignment position
                                                    for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                                    if re.match('NM_', str(for_location_c)):
                                                        for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                                    gps = for_location_c.posedit.pos.end.base
                                                    gpe = for_location_c.posedit.pos.end.base + 1
                                                    gap_position = ' between positions c.' + str(gps) + '_' + str(
                                                        gpe) + '\n'
                                                    # Warn update
                                                    auto_info = auto_info + '%s' % (gap_position)
                                                elif re.search(r'\-', str(
                                                        tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                                                    r'\-',
                                                    str(
                                                        tx_hgvs_not_delins.posedit.pos.end)):
                                                    auto_info = auto_info + str(
                                                        stored_hgvs_not_delins.ac) + ':g.' + str(
                                                        stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                                                        disparity_deletion_in[
                                                            1]) + ' genomic base(s) that fail to align to transcript ' + str(
                                                        tx_hgvs_not_delins.ac)
                                                    non_valid_caution = 'true'
                                                    try:
                                                        c2 = self.vm.n_to_c(tx_hgvs_not_delins)
                                                    except:
                                                        c2 = tx_hgvs_not_delins
                                                    c1 = copy.deepcopy(c2)
                                                    c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                                                    c1.posedit.pos.start.offset = 0
                                                    c1.posedit.pos.end = c2.posedit.pos.start
                                                    c1.posedit.edit.ref = ''
                                                    c1.posedit.edit.alt = ''
                                                    if orientation != -1:
                                                        g1 = self.vm.t_to_g(c1, hgvs_genomic.ac)
                                                        g2 = self.vm.t_to_g(c2, hgvs_genomic.ac)
                                                        g1.posedit.edit.alt = g1.posedit.edit.ref
                                                    else:
                                                        g1 = self.vm.t_to_g(c2, hgvs_genomic.ac)
                                                        g2 = self.vm.t_to_g(c1, hgvs_genomic.ac)
                                                        g2.posedit.edit.alt = g2.posedit.edit.ref
                                                    reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                                    alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                                    g3 = copy.deepcopy(g1)
                                                    g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                                    g3.posedit.edit.ref = reference
                                                    g3.posedit.edit.alt = alternate
                                                    c3 = self.vm.g_to_t(g3, c1.ac)
                                                    hgvs_refreshed_variant = c3
                                                    # Alignment position
                                                    for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                                    if re.match('NM_', str(for_location_c)):
                                                        for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                                    gps = for_location_c.posedit.pos.start.base - 1
                                                    gpe = for_location_c.posedit.pos.start.base
                                                    gap_position = ' between positions c.' + str(gps) + '_' + str(
                                                        gpe) + '\n'
                                                    # Warn update
                                                    auto_info = auto_info + '%s' % (gap_position)
                                                elif re.search(r'\-', str(
                                                        tx_hgvs_not_delins.posedit.pos.end)) and not re.search(r'\-',
                                                                                                               str(
                                                                                                                   tx_hgvs_not_delins.posedit.pos.start)):
                                                    auto_info = auto_info + 'Genome position ' + str(
                                                        stored_hgvs_not_delins.ac) + ':g.' + str(
                                                        stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(
                                                        disparity_deletion_in[1]) + '-bp gap in transcript ' + str(
                                                        tx_hgvs_not_delins.ac)
                                                    gapped_transcripts = gapped_transcripts + ' ' + str(
                                                        tx_hgvs_not_delins.ac)
                                                    non_valid_caution = 'true'
                                                    try:
                                                        c1 = self.vm.n_to_c(tx_hgvs_not_delins)
                                                    except:
                                                        c1 = tx_hgvs_not_delins
                                                    c2 = copy.deepcopy(c1)
                                                    c2.posedit.pos.start = c1.posedit.pos.end
                                                    c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                                                    c2.posedit.pos.end.offset = 0
                                                    c2.posedit.edit.ref = ''
                                                    c2.posedit.edit.alt = ''
                                                    if orientation != -1:
                                                        g1 = self.vm.t_to_g(c1, hgvs_genomic.ac)
                                                        g2 = self.vm.t_to_g(c2, hgvs_genomic.ac)
                                                        g2.posedit.edit.alt = g2.posedit.edit.ref
                                                    else:
                                                        g1 = self.vm.t_to_g(c2, hgvs_genomic.ac)
                                                        g2 = self.vm.t_to_g(c1, hgvs_genomic.ac)
                                                        g1.posedit.edit.alt = g1.posedit.edit.ref
                                                    reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                                    alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                                    g3 = copy.deepcopy(g1)
                                                    g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                                    g3.posedit.edit.ref = reference
                                                    g3.posedit.edit.alt = alternate
                                                    c3 = self.vm.g_to_t(g3, c1.ac)
                                                    hgvs_refreshed_variant = c3
                                                    # Alignment position
                                                    for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                                    if re.match('NM_', str(for_location_c)):
                                                        for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                                    gps = for_location_c.posedit.pos.end.base - 1
                                                    gpe = for_location_c.posedit.pos.end.base
                                                    gap_position = ' between positions c.' + str(gps) + '_' + str(
                                                        gpe) + '\n'
                                                    # Warn update
                                                    auto_info = auto_info + '%s' % (gap_position)
                                                else:
                                                    auto_info = auto_info + str(
                                                        stored_hgvs_not_delins.ac) + ':g.' + str(
                                                        stored_hgvs_not_delins.posedit.pos) + ' contains ' + str(
                                                        disparity_deletion_in[
                                                            1]) + ' genomic base(s) that fail to align to transcript ' + str(
                                                        tx_hgvs_not_delins.ac) + '\n'
                                                    hgvs_refreshed_variant = tx_hgvs_not_delins

                                        # GAP IN THE CHROMOSOME
                                        elif disparity_deletion_in[0] == 'chromosome':
                                            # amend_RefSeqGene = 'true'
                                            if possibility_counter == 3:
                                                hgvs_refreshed_variant = stash_tx_right
                                            elif possibility_counter == 4:
                                                hgvs_refreshed_variant = stash_tx_left
                                            else:
                                                hgvs_refreshed_variant = chromosome_normalized_hgvs_coding
                                            # Warn
                                            auto_info = auto_info + str(hgvs_refreshed_variant.ac) + ':c.' + str(
                                                hgvs_refreshed_variant.posedit.pos) + ' contains ' + str(
                                                disparity_deletion_in[
                                                    1]) + ' transcript base(s) that fail to align to chromosome ' + str(
                                                hgvs_genomic.ac) + '\n'
                                        else:
                                            # Keep the same by re-setting rel_var
                                            hgvs_refreshed_variant = hgvs_coding
                                        # amend_RefSeqGene = 'false'

                                        # Edit the output
                                        if re.match('NM_', str(hgvs_refreshed_variant.ac)) and not re.search('c',
                                                                                                             str(
                                                                                                                 hgvs_refreshed_variant.type)):
                                            hgvs_refreshed_variant = no_norm_evm.n_to_c(hgvs_refreshed_variant)
                                        else:
                                            pass

                                        try:
                                            hn.normalize(hgvs_refreshed_variant)
                                        except Exception as e:
                                            error = str(e)
                                            # Ensure the final variant is not intronic nor does it cross exon boundaries
                                            if re.match('Normalization of intronic variants is not supported',
                                                        error) or re.match(
                                                'Unsupported normalization of variants spanning the exon-intron boundary',
                                                error):
                                                hgvs_refreshed_variant = saved_hgvs_coding
                                            else:
                                                continue

                                        # Quick check to make sure the coding variant has not changed
                                        try:
                                            to_test = hn.normalize(hgvs_refreshed_variant)
                                        except:
                                            to_test = hgvs_refreshed_variant
                                        if str(to_test.posedit.edit) != str(hgvs_coding.posedit.edit):
                                            # Try the next available genomic option
                                            if hgvs_coding.posedit.edit.type == 'identity' and to_test.posedit.edit.type == 'identity':
                                                hgvs_coding = to_test
                                            else:
                                                continue

                                        # Update hgvs_genomic
                                        hgvs_alt_genomic = self.myvm_t_to_g(hgvs_refreshed_variant, alt_chr,
                                                                        no_norm_evm, hn)
                                        if hgvs_alt_genomic.posedit.edit.type == 'identity':
                                            re_c = self.vm.g_to_t(hgvs_alt_genomic, hgvs_refreshed_variant.ac)
                                            if (hn.normalize(re_c)) != (hn.normalize(hgvs_refreshed_variant)):
                                                shuffle_left_g = copy.copy(hgvs_alt_genomic)
                                                shuffle_left_g.posedit.edit.ref = ''
                                                shuffle_left_g.posedit.edit.alt = ''
                                                shuffle_left_g.posedit.pos.start.base = shuffle_left_g.posedit.pos.start.base - 1
                                                shuffle_left_g.posedit.pos.end.base = shuffle_left_g.posedit.pos.end.base - 1
                                                shuffle_left_g = reverse_normalizer.normalize(shuffle_left_g)
                                                re_c = self.vm.g_to_t(shuffle_left_g, hgvs_refreshed_variant.ac)
                                                if (hn.normalize(re_c)) != (hn.normalize(hgvs_refreshed_variant)):
                                                    hgvs_alt_genomic = shuffle_left_g

                                                    # If it is intronic, these vairables will not have been set
                                    else:
                                        # amend_RefSeqGene = 'false'
                                        no_normalized_c = 'false'

                                    # Break if gap has been detected
                                    if disparity_deletion_in[0] != 'false':
                                        break

                                # Normailse hgvs_genomic
                                try:
                                    hgvs_alt_genomic = hn.normalize(hgvs_alt_genomic)
                                except hgvs.exceptions.HGVSError as e:
                                    # Strange error caused by gap in genomic
                                    error = str(e)
                                    if re.search('base start position must be <= end position', error) and \
                                            disparity_deletion_in[0] == 'chromosome':
                                        if hgvs_alt_genomic.posedit.edit.type == 'delins':
                                            start = hgvs_alt_genomic.posedit.pos.start.base
                                            end = hgvs_alt_genomic.posedit.pos.end.base
                                            lhb = self.sf.fetch_seq(str(hgvs_alt_genomic.ac), end - 1, end)
                                            rhb = self.sf.fetch_seq(str(hgvs_alt_genomic.ac), start - 1, start)
                                            hgvs_alt_genomic.posedit.edit.ref = lhb + rhb
                                            hgvs_alt_genomic.posedit.edit.alt = lhb + hgvs_alt_genomic.posedit.edit.alt + rhb
                                            hgvs_alt_genomic.posedit.pos.start.base = end
                                            hgvs_alt_genomic.posedit.pos.end.base = start
                                            hgvs_alt_genomic = hn.normalize(hgvs_alt_genomic)
                                        if hgvs_alt_genomic.posedit.edit.type == 'del':
                                            start = hgvs_alt_genomic.posedit.pos.start.base
                                            end = hgvs_alt_genomic.posedit.pos.end.base
                                            lhb = self.sf.fetch_seq(str(hgvs_alt_genomic.ac), end - 1, end)
                                            rhb = self.sf.fetch_seq(str(hgvs_alt_genomic.ac), start - 1, start)
                                            hgvs_alt_genomic.posedit.edit.ref = lhb + rhb
                                            hgvs_alt_genomic.posedit.edit.alt = lhb + rhb
                                            hgvs_alt_genomic.posedit.pos.start.base = end
                                            hgvs_alt_genomic.posedit.pos.end.base = start
                                            hgvs_alt_genomic = hn.normalize(hgvs_alt_genomic)

                                # Refresh the :g. variant
                                multi_g.append(hgvs_alt_genomic)
                            else:
                                multi_g.append(hgvs_alt_genomic)
                                corrective_action_taken = 'false'

                        # In this instance, the gap code has generally found an incomplete-alignment rather than a
                        # truly gapped alignment.
                        except KeyError:
                            warnings = warnings + ': Suspected incomplete alignment between transcript %s and ' \
                                                  'genomic reference sequence %s' % (hgvs_coding.ac,
                                                                                     alt_chr)
                            continue
                        except hgvs.exceptions.HGVSError as e:
                            exc_type, exc_value, last_traceback = sys.exc_info()
                            te = traceback.format_exc()
                            error = str(te)
                            logger.error(str(exc_type) + " " + str(exc_value))
                            logger.debug(error)
                            continue

                    if multi_g != []:

                        multi_gen_vars = multi_g  # '|'.join(multi_g)
                    else:
                        multi_gen_vars = []
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
