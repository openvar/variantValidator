import hgvs
import hgvs.exceptions
import hgvs.normalizer
import re
import copy
import sys
import logging
import json
from hgvs.assemblymapper import AssemblyMapper
from . import hgvs_utils
from . import utils as fn
from . import seq_data
from . import vvMixinConverters
from .variant import Variant
from . import format_converters
from . import use_checking
from . import mappers
from . import valoutput
from .liftover import liftover

logger = logging.getLogger(__name__)


class Mixin(vvMixinConverters.Mixin):
    """
    This module contains the main function for variant validator.
    It's added to the Validator object in the vvObjects file.
    """

    def validate(self, batch_variant, selected_assembly, select_transcripts, transcript_set="refseq"):
        """
        This is the main validator function.
        :param batch_variant: A string containing the variant to be validated
        :param selected_assembly: The version of the genome assembly to use.
        :param select_transcripts: Can be an array of different transcripts, or 'all'
        Selecting multiple transcripts will lead to a multiple variant outputs.
        :param transcript_set: 'refseq' or 'ensembl'. Currently only 'refseq' is supported
        :return:
        """
        logger.debug("Running validate with inputs %s and assembly %s", batch_variant, selected_assembly)

        if transcript_set == "refseq":
            self.alt_aln_method = 'splign'
        elif transcript_set == "ensembl":
            self.alt_aln_method = 'genebuild'
            logger.warning("Ensembl is currently not supported")
            raise Exception("Ensembl is currently not supported")
        else:
            raise Exception("The transcriptSet variable '%s' is invalid, it must be 'refseq' or 'ensembl'" %
                            transcript_set)

        primary_assembly = None

        self.selected_assembly = selected_assembly
        self.select_transcripts = select_transcripts

        try:
            # Validation
            ############

            # Create a dictionary of transcript ID : ''
            select_transcripts_dict = {}
            select_transcripts_dict_plus_version = {}
            if select_transcripts != 'all':
                select_transcripts_list = select_transcripts.split('|')
                for trans_id in select_transcripts_list:
                    trans_id = trans_id.strip()
                    if 'LRG' in trans_id:
                        trans_id = self.db.get_refseq_transcript_id_from_lrg_transcript_id(trans_id)
                        if trans_id == 'none':
                            continue
                    select_transcripts_dict_plus_version[trans_id] = ''
                    trans_id = trans_id.split('.')[0]
                    select_transcripts_dict[trans_id] = ''

            # split the batch queries into a list
            batch_queries = batch_variant.split('|')

            # Turn each variant into a dictionary. The dictionary will be compiled during validation
            self.batch_list = []
            for queries in batch_queries:
                queries = queries.strip()
                query = Variant(queries)
                self.batch_list.append(query)
                logger.info("Submitting variant with format %s", queries)

            # Create List to carry batch data output
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

                # Create Normalizers
                my_variant.hn = hgvs.normalizer.Normalizer(self.hdp,
                                                           cross_boundaries=False,
                                                           shuffle_direction=3,
                                                           alt_aln_method=self.alt_aln_method
                                                           )
                my_variant.reverse_normalizer = hgvs.normalizer.Normalizer(self.hdp,
                                                                           cross_boundaries=False,
                                                                           shuffle_direction=5,
                                                                           alt_aln_method=self.alt_aln_method
                                                                           )
                # This will be used to order the final output
                if not my_variant.order:
                    ordering = ordering + 1
                    my_variant.order = ordering

                # Bug catcher
                try:
                    # Note, ID is not touched. It is always the input variant description.
                    # Quibble will be altered but id will not if type = g.
                    logger.info("Started validation of %s (originally %s)", str(my_variant.quibble),
                                 my_variant.original)

                    if not my_variant.is_ascii():
                        chars, positions = my_variant.get_non_ascii()
                        error = 'Submitted variant description contains an invalid character(s) %s at position(s) %s: '\
                                'Please remove this character and re-submit: A useful search function for ' \
                                'Unicode characters can be found at https://unicode-search.net/' % (chars, positions)
                        my_variant.warnings.append(error)
                        logger.warning(error)
                        continue

                    # Remove whitespace
                    my_variant.remove_whitespace()
                    if my_variant.quibble != my_variant.original:
                        caution = 'Whitespace removed from variant description %s' % my_variant.original
                        my_variant.warnings.append(caution)
                        logger.debug(caution)

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
                        if primary_assembly in self.genome_builds or primary_assembly == 'hg38':
                            my_variant.primary_assembly = primary_assembly
                        else:
                            my_variant.primary_assembly = 'GRCh38'
                            primary_assembly = 'GRCh38'
                            my_variant.warnings.append('Invalid genome build has been specified. Automap has selected '
                                                       'the default build (GRCh38)')
                            logger.warning(
                                'Invalid genome build has been specified. Automap has selected the '
                                'default build ' + my_variant.primary_assembly)
                    else:
                        primary_assembly = my_variant.primary_assembly
                    logger.debug("Completed string formatting")

                    toskip = format_converters.initial_format_conversions(my_variant, self,
                                                                          select_transcripts_dict_plus_version)
                    if toskip:
                        continue

                    # INITIAL USER INPUT FORMATTING
                    invalid = my_variant.format_quibble()
                    if invalid:
                        if re.search(r'\w+:[gcnmrp]', my_variant.quibble) and not \
                                re.search(r'\w+:[gcnmrp]\.', my_variant.quibble):
                            error = 'Variant description ' + my_variant.quibble + ' lacks the . character between ' \
                                    '<type> and <position> in the expected pattern <accession>:<type>.<position>'
                        else:
                            error = 'Variant description ' + my_variant.quibble + ' is not in an accepted format'
                        my_variant.warnings.append(error)
                        logger.warning(error)
                        continue

                    formatted_variant = my_variant.quibble
                    stash_input = my_variant.quibble
                    my_variant.post_format_conversion = stash_input
                    format_type = my_variant.reftype

                    hgnc_gene_info = 'false'

                    logger.debug("Variant input formatted, proceeding to validate.")

                    # Conversions
                    # Conversions are not currently supported. The HGVS format for conversions
                    # is rarely seen wrt genomic sequencing data and needs to be re-evaluated
                    if 'con' in my_variant.quibble:
                        my_variant.warnings.append('Gene conversions currently unsupported')
                        logger.warning('Gene conversions currently unsupported')
                        continue

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
                        my_variant.hgvs_formatted = input_parses
                    except hgvs.exceptions.HGVSError as e:
                        my_variant.warnings.append(str(e))
                        logger.warning(str(e))
                        continue

                    if 'LRG' in my_variant.hgvs_formatted.ac:
                        my_variant.hgvs_formatted.ac.replace('T', 't')
                    else:
                        my_variant.hgvs_formatted.ac = my_variant.hgvs_formatted.ac.upper()
                    if hasattr(my_variant.hgvs_formatted.posedit.edit, 'alt'):
                        if my_variant.hgvs_formatted.posedit.edit.alt is not None:
                            my_variant.hgvs_formatted.posedit.edit.alt = \
                                my_variant.hgvs_formatted.posedit.edit.alt.upper()
                    if hasattr(my_variant.hgvs_formatted.posedit.edit, 'ref'):
                        if my_variant.hgvs_formatted.posedit.edit.ref is not None:
                            my_variant.hgvs_formatted.posedit.edit.ref = \
                                my_variant.hgvs_formatted.posedit.edit.ref.upper()
                    formatted_variant = str(my_variant.hgvs_formatted)

                    my_variant.set_quibble(str(my_variant.hgvs_formatted))

                    # ENST support needs to be re-evaluated, but is very low priority
                    # ENST not supported by ACMG and is under review by HGVS
                    if my_variant.refsource == 'ENS':
                        trap_ens_in = str(my_variant.hgvs_formatted)
                        sim_tx = self.hdp.get_similar_transcripts(my_variant.hgvs_formatted.ac)
                        for line in sim_tx:
                            if line[2] and line[3] and line[4] and line[5] and line[6]:
                                my_variant.hgvs_formatted.ac = line[1]
                                my_variant.set_quibble(str(my_variant.hgvs_formatted))
                                formatted_variant = my_variant.quibble
                                break
                        if my_variant.refsource == 'ENS':
                            error = 'Unable to map ' + my_variant.hgvs_formatted.ac + \
                                    ' to an equivalent RefSeq transcript'
                            my_variant.warnings.append(error)
                            logger.warning(error)
                            continue
                        else:
                            my_variant.warnings.append(str(trap_ens_in) + ' automapped to equivalent RefSeq transcript '
                                                       + my_variant.quibble)
                            logger.info(str(trap_ens_in) + ' automapped to equivalent RefSeq '
                                                              'transcript ' + my_variant.quibble)
                    logger.debug("HVGS acceptance test passed")

                    # Check whether supported genome build is requested for non g. descriptions
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
                        my_variant.evm = AssemblyMapper(self.hdp,
                                                        assembly_name=primary_assembly,
                                                        alt_aln_method=self.alt_aln_method,
                                                        normalize=True,
                                                        replace_reference=True
                                                        )

                        # Setup a reverse normalize instance and non-normalize evm
                        my_variant.no_norm_evm = AssemblyMapper(self.hdp,
                                                                assembly_name=primary_assembly,
                                                                alt_aln_method=self.alt_aln_method,
                                                                normalize=False,
                                                                replace_reference=True
                                                                )

                        # Create a specific minimal evm with no normalizer and no replace_reference
                        my_variant.min_evm = AssemblyMapper(self.hdp,
                                                            assembly_name=primary_assembly,
                                                            alt_aln_method=self.alt_aln_method,
                                                            normalize=False,
                                                            replace_reference=False
                                                            )

                    else:
                        error = 'Mapping of ' + formatted_variant + ' to genome assembly ' + \
                                primary_assembly + ' is not supported'
                        my_variant.warnings.append(error)
                        logger.warning(error)
                        continue

                    # Catch interval end > interval start
                    # hgvs did/does not handle 3' UTR position ordering well. This function
                    # ensures that end pos is not > start pos wrt 3' UTRs.
                    # Also identifies some variants which span into the downstream sequence
                    # i.e. out of bounds
                    if '*' in str(my_variant.hgvs_formatted.posedit):
                        input_parses_copy = copy.deepcopy(my_variant.hgvs_formatted)
                        input_parses_copy.type = "c"
                        # Map to n. position
                        # Create easy variant mapper (over variant mapper) and splign locked evm
                        try:
                            to_n = my_variant.evm.c_to_n(input_parses_copy)
                        except hgvs.exceptions.HGVSError as e:
                            logger.debug("Except passed, %s", e)
                        else:
                            if to_n.posedit.pos.end.base < to_n.posedit.pos.start.base:
                                error = 'Interval end position < interval start position '
                                my_variant.warnings.append(error)
                                logger.warning(error)
                                continue
                    elif my_variant.hgvs_formatted.posedit.pos.end.base < my_variant.hgvs_formatted.posedit.pos.start.base:
                        error = 'Interval end position ' + str(my_variant.hgvs_formatted.posedit.pos.end.base) + \
                                ' < interval start position ' + str(my_variant.hgvs_formatted.posedit.pos.start.base)
                        my_variant.warnings.append(error)
                        logger.warning(error)
                        continue

                    # Catch missing version number in refseq
                    is_version = re.compile(r"\d\.\d")
                    if my_variant.refsource == 'RefSeq' and not is_version.search(str(my_variant.hgvs_formatted)):
                        error = 'RefSeq variant accession numbers MUST include a version number'
                        my_variant.warnings.append(error)
                        continue
                    logger.debug("HVGS interval/version mapping complete")

                    # handle LRG inputs

                    if my_variant.refsource == 'LRG':
                        format_converters.lrg_to_refseq(my_variant, self)
                        logger.debug("LRG check for conversion to refseq completed")

                    # Additional Incorrectly input variant capture training
                    if my_variant.refsource == 'RefSeq':
                        toskip = use_checking.refseq_common_mistakes(my_variant)
                        if toskip:
                            continue
                        logger.debug("Passed 'common mistakes' catcher")

                    # Primary validation of the input
                    toskip = use_checking.structure_checks(my_variant, self)
                    if toskip:
                        continue
                    logger.debug("Variant structure and contents searches passed")

                    # Mitochondrial variants
                    toskip = format_converters.mitochondrial(my_variant, self)
                    if toskip:
                        continue

                    toskip = format_converters.proteins(my_variant, self)
                    if toskip:
                        continue

                    trapped_input = str(my_variant.hgvs_formatted)
                    my_variant.pre_RNA_conversion = trapped_input
                    toskip = format_converters.rna(my_variant, self)
                    if toskip:
                        continue

                    # COLLECT gene symbol, name and ACCESSION INFORMATION
                    # Gene symbol
                    if my_variant.reftype != ':g.':
                        toskip = self._get_transcript_info(my_variant)
                        if toskip:
                            continue

                    # Now start mapping from genome to transcripts
                    if my_variant.reftype == ':g.':
                        toskip = mappers.gene_to_transcripts(my_variant, self)
                        if toskip:
                            continue

                    if format_type == ':c.' or format_type == ':n.':
                        toskip = mappers.transcripts_to_gene(my_variant, self, select_transcripts_dict_plus_version)
                        if toskip:
                            continue

                    # Set the data
                    my_variant.output_type_flag = 'gene'
                    my_variant.description = hgnc_gene_info
                    my_variant.primary_assembly = primary_assembly
                    logger.info("Completed initial validation for %s", my_variant.quibble)
                # Report errors to User and VV admin
                except KeyboardInterrupt:
                    raise
                except Exception:
                    my_variant.output_type_flag = 'error'
                    error = 'Validation error'
                    my_variant.warnings.append(error)
                    exc_type, exc_value, last_traceback = sys.exc_info()
                    logger.error(str(exc_type) + " " + str(exc_value))
                    raise

            # Outside the for loop
            ######################
            logger.debug("End of 1st for loop")
            # order the rows
            by_order = sorted(self.batch_list, key=lambda x: x.order)

            for variant in by_order:
                if not variant.write:
                    continue

                # Genomic sequence variation
                genomic_variant = variant.genomic_g
                hgvs_genomic_variant = genomic_variant

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
                    variant.warnings.append(refseqgene_variant)
                    refseqgene_variant = ''
                    lrg_variant = ''
                    hgvs_refseqgene_variant = 'false'
                else:
                    hgvs_refseqgene_variant = self.hp.parse_hgvs_variant(refseqgene_variant)
                    rsg_ac = self.db.get_lrg_id_from_refseq_gene_id(str(hgvs_refseqgene_variant.ac))
                    if rsg_ac[0] == 'none':
                        lrg_variant = ''
                    else:
                        hgvs_lrg = copy.deepcopy(hgvs_refseqgene_variant)
                        hgvs_lrg.ac = rsg_ac[0]
                        lrg_variant = fn.valstr(hgvs_lrg)
                        if rsg_ac[1] != 'public':
                            variant.warnings.append('The current status of ' + str(hgvs_lrg.ac) + ' is pending '
                                                    'therefore changes may be made to the LRG reference sequence')

                # Transcript sequence variation
                tx_variant = variant.coding
                hgvs_transcript_variant = tx_variant
                hgvs_tx_variant = None
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
                    lrg_transcript = self.db.get_lrg_transcript_id_from_refseq_transcript_id(transcript_accession)
                    if lrg_transcript == 'none':
                        lrg_transcript_variant = ''
                    else:
                        # Note - LRG availability is dependant on UTA containing the data. In some
                        # instances we will be able to display the LRG_tx without being able to
                        # display the LRG gene data

                        try:
                            hgvs_lrg_t = self.vm.g_to_t(hgvs_refseqgene_variant, transcript_accession)
                            hgvs_lrg_t.ac = lrg_transcript
                            lrg_transcript_variant = fn.valstr(hgvs_lrg_t)
                        except Exception:
                            if hgvs_transcript_variant.posedit.pos.start.offset == 0 and \
                                    hgvs_transcript_variant.posedit.pos.end.offset == 0:
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
                            genome_context_transcript_variant = genomic_accession + '(' + transcript_accession +\
                                                                '):c.' + str(hgvs_transcript_variant.posedit)
                            if refseqgene_variant != '':
                                hgvs_refseqgene_variant = self.hp.parse_hgvs_variant(refseqgene_variant)
                                refseqgene_accession = hgvs_refseqgene_variant.ac
                                hgvs_coding_from_refseqgene = self.vm.g_to_t(hgvs_refseqgene_variant,
                                                                             hgvs_transcript_variant.ac)
                                hgvs_coding_from_refseqgene = fn.valstr(hgvs_coding_from_refseqgene)
                                hgvs_coding_from_refseqgene = self.hp.parse_hgvs_variant(hgvs_coding_from_refseqgene)
                                refseqgene_context_transcript_variant = refseqgene_accession + '(' + \
                                    transcript_accession + '):c.' + str(hgvs_coding_from_refseqgene.posedit.pos) + str(
                                        hgvs_coding_from_refseqgene.posedit.edit)
                            else:
                                refseqgene_context_transcript_variant = ''
                        else:
                            genome_context_transcript_variant = ''  # transcript_variant
                            refseqgene_context_transcript_variant = ''
                    else:
                        genome_context_transcript_variant = ''  # transcript_variant
                        refseqgene_context_transcript_variant = ''
                else:
                    genome_context_transcript_variant = ''
                    refseqgene_context_transcript_variant = ''

                # Protein description
                predicted_protein_variant = variant.protein
                if 'NP_' in predicted_protein_variant:
                    rs_p, pred_prot_posedit = predicted_protein_variant.split(':')
                    lrg_p = self.db.get_lrg_protein_id_from_ref_seq_protein_id(rs_p)
                    if 'LRG' in lrg_p:
                        predicted_protein_variant = rs_p + '(' + lrg_p + '):' + pred_prot_posedit

                # Gene
                if transcript_accession == '':
                    variant.gene_symbol = ''

                if tx_variant != '':
                    multi_gen_vars = mappers.final_tx_to_multiple_genomic(variant, self, tx_variant)

                else:
                    # HGVS genomic in the absence of a transcript variant
                    if genomic_variant != '':
                        multi_gen_vars = [hgvs_genomic_variant]
                    else:
                        multi_gen_vars = []

                # Dictionaries of genomic loci
                alt_genomic_dicts = []
                primary_genomic_dicts = {}

                for alt_gen_var in multi_gen_vars:
                    try:
                        alt_gen_var = variant.hn.normalize(alt_gen_var)
                    except hgvs.exceptions.HGVSInvalidVariantError:
                        continue
                    for build in self.genome_builds:
                        test = seq_data.supported_for_mapping(alt_gen_var.ac, build)
                        if test:
                            try:
                                vcf_dict = hgvs_utils.report_hgvs2vcf(alt_gen_var, build, variant.reverse_normalizer,
                                                                      self.sf)
                            except hgvs.exceptions.HGVSInvalidVariantError:
                                continue
                            # Identify primary assembly positions
                            if 'NC_' in alt_gen_var.ac:
                                if 'GRC' in build:
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
                                    vcf_dict = hgvs_utils.report_hgvs2vcf(alt_gen_var, 'hg38', variant.reverse_normalizer,
                                                                          self.sf)
                                    primary_genomic_dicts['hg38'] = {
                                        'hgvs_genomic_description': fn.valstr(alt_gen_var),
                                        'vcf': {'chr': vcf_dict['ucsc_chr'],
                                                'pos': vcf_dict['pos'],
                                                'ref': vcf_dict['ref'],
                                                'alt': vcf_dict['alt']
                                                }
                                    }

                            else:
                                if 'GRC' in build:
                                    alt_dict = {build.lower(): {'hgvs_genomic_description': fn.valstr(alt_gen_var),
                                                                'vcf': {'chr': vcf_dict['grc_chr'],
                                                                        'pos': vcf_dict['pos'],
                                                                        'ref': vcf_dict['ref'],
                                                                        'alt': vcf_dict['alt']
                                                                        }
                                                                }
                                                }
                                else:
                                    alt_dict = {build.lower(): {'hgvs_genomic_description': fn.valstr(alt_gen_var),
                                                                'vcf': {'chr': vcf_dict['ucsc_chr'],
                                                                        'pos': vcf_dict['pos'],
                                                                        'ref': vcf_dict['ref'],
                                                                        'alt': vcf_dict['alt']
                                                                        }
                                                                }
                                                }
                                # Append
                                alt_genomic_dicts.append(alt_dict)

                                if build == 'GRCh38':
                                    vcf_dict = hgvs_utils.report_hgvs2vcf(alt_gen_var, 'hg38', variant.reverse_normalizer,
                                                                          self.sf)
                                    alt_dict = {'hg38': {'hgvs_genomic_description': fn.valstr(alt_gen_var),
                                                         'vcf': {'chr': vcf_dict['ucsc_chr'],
                                                                 'pos': vcf_dict['pos'],
                                                                 'ref': vcf_dict['ref'],
                                                                 'alt': vcf_dict['alt']
                                                                 }
                                                         }
                                                }
                                    # Append
                                    alt_genomic_dicts.append(alt_dict)

                # Warn not directly mapped to specified genome build
                if genomic_accession != '':
                    if primary_assembly.lower() not in list(primary_genomic_dicts.keys()):
                        variant.warnings.extend([
                            str(variant.hgvs_coding) + ' cannot be mapped directly to genome build ' + primary_assembly,
                            'See alternative genomic loci or alternative genome builds for aligned genomic positions'
                        ])

                # Ensure Variants have had the refs removed.
                # if not hasattr(posedit, refseqgene_variant):
                if refseqgene_variant != '':
                    try:
                        refseqgene_variant = fn.valstr(hgvs_refseqgene_variant)
                    except Exception as e:
                        logger.debug("Except passed, %s", e)

                # Add single letter AA code to protein descriptions
                predicted_protein_variant_dict = {"tlr": str(predicted_protein_variant), "slr": ''}
                if predicted_protein_variant != '':
                    if 'Non-coding :n.' not in predicted_protein_variant:
                        try:
                            format_p = predicted_protein_variant
                            format_p = re.sub(r'\(LRG_.+?\)', '', format_p)
                            re_parse_protein = self.hp.parse_hgvs_variant(format_p)
                            re_parse_protein_single_aa = fn.single_letter_protein(re_parse_protein)
                            predicted_protein_variant_dict["slr"] = str(re_parse_protein_single_aa)
                        except hgvs.exceptions.HGVSParseError as e:
                            logger.debug("Except passed, %s", e)
                    else:
                        predicted_protein_variant_dict["slr"] = str(predicted_protein_variant)

                # Add stable gene_ids
                stable_gene_ids = {}
                if variant.gene_symbol != '':
                    gene_stable_info = self.db.get_stable_gene_id_info(variant.gene_symbol)

                    # Add or update stable ID and transcript data
                    if gene_stable_info[1] == 'No data' and hgvs_tx_variant is not None:
                        self.db.update_transcript_info_record(hgvs_tx_variant.ac, self)
                        gene_stable_info = self.db.get_stable_gene_id_info(variant.gene_symbol)

                    # Update gene_symbol
                    if variant.gene_symbol != str(gene_stable_info[1]) and str(gene_stable_info[1]) != 'No data':
                        variant.gene_symbol = str(gene_stable_info[1])

                    try:
                        # Dictionary the output
                        stable_gene_ids['hgnc_id'] = gene_stable_info[2]
                        stable_gene_ids['entrez_gene_id'] = gene_stable_info[3]
                        # stable_gene_ids['ensembl_gene_id'] = gene_stable_info[4]
                        stable_gene_ids['ucsc_id'] = gene_stable_info[5]
                        stable_gene_ids['omim_id'] = json.loads(gene_stable_info[6])
                        # stable_gene_ids['vega_id'] = gene_stable_info[7]
                        # stable_gene_ids['ccds_id'] = gene_stable_info[8]
                    except IndexError as e:
                        logger.debug("Except pass, %s", e)

                variant.stable_gene_ids = stable_gene_ids
                variant.hgvs_transcript_variant = tx_variant
                variant.genome_context_intronic_sequence = genome_context_transcript_variant
                variant.refseqgene_context_intronic_sequence = refseqgene_context_transcript_variant
                variant.hgvs_refseqgene_variant = refseqgene_variant
                variant.hgvs_predicted_protein_consequence = predicted_protein_variant_dict
                variant.hgvs_lrg_transcript_variant = lrg_transcript_variant
                variant.hgvs_lrg_variant = lrg_variant
                variant.alt_genomic_loci = alt_genomic_dicts
                variant.primary_assembly_loci = primary_genomic_dicts
                variant.reference_sequence_records = ''
                variant.validated = True

                # Add links to reference_sequence_records
                ref_records = self.db.get_urls(variant.output_dict())
                if ref_records != {}:
                    variant.reference_sequence_records = ref_records

                if variant.output_type_flag == 'intergenic':
                    # Attempt to liftover between genome builds
                    # Note: pyliftover uses the UCSC liftOver tool.
                    # https://pypi.org/project/pyliftover/
                    genomic_position_info = variant.primary_assembly_loci
                    for g_p_key in list(genomic_position_info.keys()):
                        build_to = ''
                        build_from = ''

                        # Identify the current build and hgvs_genomic descripsion
                        if 'hg' in g_p_key:
                            # incoming_vcf = genomic_position_info[g_p_key]['vcf']
                            # set builds
                            if g_p_key == 'hg38':
                                build_to = 'hg19'
                                build_from = 'hg38'
                            if g_p_key == 'hg19':
                                build_to = 'hg38'
                                build_from = 'hg19'
                        elif 'grc' in g_p_key:
                            # incoming_vcf = genomic_position_info[g_p_key]['vcf']
                            # set builds
                            if g_p_key == 'grch38':
                                build_to = 'GRCh37'
                                build_from = 'GRCh38'
                            if g_p_key == 'grch37':
                                build_to = 'GRCh38'
                                build_from = 'GRCh37'

                        # Liftover
                        lifted_response = liftover(genomic_position_info[g_p_key]['hgvs_genomic_description'],
                                                   build_from,
                                                   build_to, variant.hn, variant.reverse_normalizer,
                                                   variant.evm, self)

                        # Sort the respomse into primary assembly and ALT
                        primary_assembly_loci = {}
                        alt_genomic_loci = []
                        for build_key, accession_dict in list(lifted_response.items()):
                            try:
                                accession_key = list(accession_dict.keys())[0]
                                if 'NC_' in accession_dict[accession_key]['hgvs_genomic_description']:
                                    primary_assembly_loci[build_key.lower()] = accession_dict[accession_key]
                                else:
                                    alt_genomic_loci.append({build_key.lower(): accession_dict[accession_key]})

                            # KeyError if the dicts are empty
                            except KeyError:
                                continue
                            except IndexError:
                                continue

                        # Add the dictionaries from lifted response to the output
                        if primary_assembly_loci != {}:
                            variant.primary_assembly_loci = primary_assembly_loci
                        if alt_genomic_loci:
                            variant.alt_genomic_loci = alt_genomic_loci

                # Append to a list for return
                batch_out.append(variant)

            output = valoutput.ValOutput(batch_out, self)
            return output

        # Bug catcher
        except KeyboardInterrupt:
            raise
        except BaseException:
            # Debug mode
            exc_type, exc_value, last_traceback = sys.exc_info()
            logger.critical(str(exc_type) + " " + str(exc_value))
            raise fn.VariantValidatorError('Validation error')

    def gene2transcripts(self, query):
        """
        Generates a list of transcript (UTA supported) and transcript names from a gene symbol or RefSeq transcript ID
        :param query: string gene symbol or RefSeq ID (e.g. NANOG or NM_024865.3)
        :return: dictionary of transcript information
        """
        query = query.upper()
        if re.search(r'\d+ORF\d+', query):
            query = query.replace('ORF', 'orf')

        # Quick check for LRG
        elif 'LRG' in query:
            lrg_id = query.split('T')[0]
            lrg_to_hgnc = self.db.get_lrg_data_from_lrg_id(lrg_id)
            if lrg_to_hgnc and lrg_to_hgnc[0] != 'none':
                query = lrg_to_hgnc[2]

        # Quick check for blank form
        if query == '':
            return {'error': 'Please enter HGNC gene name or transcript identifier (NM_, NR_, or ENST)'}

        hgnc = query
        if 'NM_' in hgnc or 'NR_' in hgnc:  # or re.match('ENST', hgnc):
            if '.' in hgnc:
                try:
                    tx_info = self.hdp.get_tx_identity_info(hgnc)
                    hgnc = tx_info[6]
                except hgvs.exceptions.HGVSError as e:
                    return {'error': str(e)}
            else:
                found_res = False
                for version in range(25):
                    refresh_hgnc = hgnc + '.' + str(version)
                    try:
                        tx_info = self.hdp.get_tx_identity_info(refresh_hgnc)
                        hgnc = tx_info[6]
                        found_res = True
                        break
                    except hgvs.exceptions.HGVSError as e:
                        logger.debug("Except passed, %s", e)
                if not found_res:
                    return {'error': 'No transcript definition for (tx_ac=' + hgnc + ')'}

        # First perform a search against the input gene symbol or the symbol inferred from UTA
        initial = fn.hgnc_rest(path="/fetch/symbol/" + hgnc)
        # Check for a record
        if str(initial['record']['response']['numFound']) != '0':
            current_sym = hgnc
            previous = initial
        # No record found, is it a previous symbol?
        else:
            # Look up current name
            current = fn.hgnc_rest(path="/search/prev_symbol/" + hgnc)
            # Look for historic names
            # If historic names = 0
            if str(current['record']['response']['numFound']) == '0':
                current_sym = hgnc
            else:
                current_sym = current['record']['response']['docs'][0]['symbol']
            # Look up previous symbols and gene name
            # Re-set the previous variable
            previous = fn.hgnc_rest(path="/fetch/symbol/" + current_sym)

        if len(previous['record']['response']['docs']) == 0:
            return {'error': 'Unable to recognise gene symbol %s' % current_sym}

        # Extract the relevant data
        if 'prev_symbol' in list(previous['record']['response']['docs'][0].keys()):
            previous_sym = previous['record']['response']['docs'][0]['prev_symbol'][0]
        else:
            previous_sym = current_sym

        # Get gene name
        if 'name' in list(previous['record']['response']['docs'][0].keys()):
            gene_name = previous['record']['response']['docs'][0]['name']
        else:
            # error = current_sym + ' is not a valid HGNC gene symbol'
            gene_name = 'Gene symbol %s not found in the HGNC database of human gene names www.genenames.org' % query
            return {'error': gene_name}

        # Look up previous name
        if 'prev_name' in list(previous['record']['response']['docs'][0].keys()):
            previous_name = previous['record']['response']['docs'][0]['prev_name'][0]
        else:
            previous_name = gene_name

        # Get transcripts
        tx_for_gene = self.hdp.get_tx_for_gene(current_sym)
        if len(tx_for_gene) == 0:
            tx_for_gene = self.hdp.get_tx_for_gene(previous_sym)
        if len(tx_for_gene) == 0:
            return {'error': 'Unable to retrieve data from the UTA, please contact admin'}

        # Loop through each transcript and get the relevant transcript description
        genes_and_tx = []
        recovered = []
        for line in tx_for_gene:
            if line[3].startswith('NM_') or line[3].startswith('NR_'):
                # Transcript ID
                tx = line[3]
                tx_description = self.db.get_transcript_description(tx)
                if tx_description == 'none':
                    self.db.update_transcript_info_record(tx, self)
                    tx_description = self.db.get_transcript_description(tx)
                # Check for duplicates
                if tx not in recovered:
                    recovered.append(tx)
                    if len(line) >= 3 and isinstance(line[1], int):
                        genes_and_tx.append({'reference': tx,
                                             'description': tx_description,
                                             'coding_start': line[1] + 1 + 1,
                                             'coding_end': line[2]
                                             })
                    else:
                        genes_and_tx.append({'reference': tx,
                                             'description': tx_description,
                                             'coding_start': 'non-coding',
                                             'coding_end': 'non-coding'
                                             })
                    # LRG information
                    lrg_transcript = self.db.get_lrg_transcript_id_from_refseq_transcript_id(tx)
                    if lrg_transcript != 'none':
                        genes_and_tx.append({'reference': lrg_transcript,
                                             'description': tx_description,
                                             'coding_start': line[1] + 1 + 1,
                                             'coding_end': line[2]
                                             })

        # Return data table
        g2d_data = {'current_symbol': current_sym,
                    'previous_symbol': previous_sym,
                    'current_name': gene_name,
                    'previous_name': previous_name,
                    'transcripts': genes_and_tx
                    }

        return g2d_data

    def hgvs2ref(self, query):
        """
        Fetch reference sequence from a HGVS variant description
        :param query:
        :return:
        """
        logger.debug('Fetching reference sequence for ' + query)
        # Dictionary to store the data
        reference = {'variant': query,
                     'start_position': '',
                     'end_position': '',
                     'warning': '',
                     'sequence': '',
                     'error': ''}
        # Step 1: parse the query. Dictionary the parse error if parsing fails
        try:
            input_hgvs_query = self.hp.parse_hgvs_variant(query)
        except Exception as e:
            reference['error'] = str(e)
            return reference
        # Step 2: If the variant is a c., it needs to transferred to n.
        try:
            hgvs_query = self.vm.c_to_n(input_hgvs_query)
        except:
            hgvs_query = input_hgvs_query

        # For transcript reference sequences
        if hgvs_query.type == 'c' or hgvs_query.type == 'n':
            # Step 4: Check for intronic sequence
            if hgvs_query.posedit.pos.start.offset != 0 and hgvs_query.posedit.pos.end.offset != 0:
                reference['warning'] = 'Intronic sequence variation: Use genomic reference sequence'
                return reference

            elif hgvs_query.posedit.pos.start.offset != 0 or hgvs_query.posedit.pos.end.offset != 0:
                reference['warning'] = 'Partial intronic sequence variation: Returning exonic and/or UTR sequence only'

        elif hgvs_query.type != 'g' and hgvs_query.type != 'p':
            return reference

        # Step 3: split the variant description into the parts required for seqfetching
        accession = hgvs_query.ac
        start = hgvs_query.posedit.pos.start.base - 1
        end = hgvs_query.posedit.pos.end.base

        # Step 5: try and fetch the sequence using SeqFetcher. Dictionary an error if this fails
        try:
            sequence = self.sf.fetch_seq(accession, start, end)
        except Exception as e:
            reference['error'] = str(e)
            logger.warning(str(e))
        else:
            reference['start_position'] = str(input_hgvs_query.posedit.pos.start)
            reference['end_position'] = str(input_hgvs_query.posedit.pos.end)
            reference['sequence'] = sequence

        # Return the resulting reference sequence and error message
        return reference

    def _get_transcript_info(self, variant):
        """
        Collect transcript information from a non-genomic variant.
        Should only be called during the validator process
        """

        hgvs_vt = self.hp.parse_hgvs_variant(str(variant.hgvs_formatted))
        try:
            self.hdp.get_tx_identity_info(str(hgvs_vt.ac))
        except hgvs.exceptions.HGVSError as e:
            error = 'Please inform UTA admin of the following error: ' + str(e)
            reason = "VariantValidator cannot recover information for transcript " + str(
                hgvs_vt.ac) + ' because it is not available in the Universal Transcript Archive'
            variant.warnings.append(reason)
            logger.warning(str(reason) + ": " + str(error))
            return True

        # Get accurate transcript descriptions from the relevant databases
        # RefSeq databases
        if self.alt_aln_method != 'genebuild':
            # Gene description  - requires GenBank search to get all the required info, i.e. transcript variant ID
            # accession number
            hgvs_object = self.hp.parse_hgvs_variant(str(variant.hgvs_formatted))
            accession = hgvs_object.ac
            # Look for the accession in our database
            # Connect to database and send request
            entry = self.db.in_entries(accession, 'transcript_info')

            # Analyse the returned data and take the necessary actions
            # If the error key exists
            if 'error' in entry:
                # Open a hgvs exception log file in append mode
                error = entry['description']
                variant.warnings.extend([str(error), 'A Database error occurred, please contact admin'])
                logger.warning(str(error) + ": A Database error occurred, please contact admin")
                return True

            # If the accession key is found
            elif 'accession' in entry:
                # If the current entry is too old
                if entry['expiry'] == 'true':
                    try:
                        entry = self.db.data_add(accession=accession, validator=self)
                    except hgvs.exceptions.HGVSError:
                        error = 'Transcript %s is not currently supported' % accession
                        variant.warnings.append(error)
                        logger.warning(error)
                        return True
                    except fn.ObsoleteSeqError as e:
                        error = 'Unable to assign transcript identity records to %s. %s' % (accession, str(e))
                        variant.warnings.append(error)
                        logger.info(error)
                        return True
                    except fn.DatabaseConnectionError as e:
                        error = '%s. Please try again later and if the problem persists contact admin.' % str(e)
                        variant.warnings.append(error)
                        logger.warning(error)
                        return True
                    variant.description = entry['description']
                    variant.gene_symbol = entry['hgnc_symbol']
                else:
                    variant.description = entry['description']
            # If the none key is found add the description to the database
            elif 'none' in entry:
                try:
                    entry = self.db.data_add(accession=accession, validator=self)
                except fn.ObsoleteSeqError as e:
                    error = 'Unable to assign transcript identity records to %s. %s' % (accession, str(e))
                    variant.warnings.append(error)
                    logger.info(error)
                    return True
                except fn.DatabaseConnectionError as e:
                    error = '%s. Please try again later and if the problem persists contact admin.' % str(e)
                    variant.warnings.append(error)
                    logger.warning(error)
                    return True
                variant.description = entry['description']
                variant.gene_symbol = entry['hgnc_symbol']

            # If no correct keys are found
            else:
                # Open a hgvs exception log file in append mode
                error = 'Unknown error type'
                variant.warnings.extend([error, ': A Database error occurred, please contact admin'])
                logger.warning(error)
                return True

        # Ensembl databases
        else:
            # accession number
            hgvs_object = self.hp.parse_hgvs_variant(str(variant.hgvs_formatted))
            accession = hgvs_object.ac
            # Look for the accession in our database
            # Connect to database and send request
            entry = self.db.in_entries(accession, 'transcript_info')

            # Analyse the returned data and take the necessary actions
            # If the error key exists
            if 'error' in entry:
                # Open a hgvs exception log file in append mode
                error = entry['description']
                variant.warnings.extend([str(error), ': A Database error occurred, please contact admin'])
                logger.warning(str(error))
                return True

            # If the accession key is found
            elif 'accession' in entry:
                # If the current entry is too old
                if entry['expiry'] == 'true':
                    entry = self.db.data_add(accession=accession, validator=self)
                    variant.description = entry['description']
                else:
                    variant.description = entry['description']
            # If the none key is found add the description to the database
            elif 'none' in entry:
                try:
                    entry = self.db.data_add(accession=accession, validator=self)
                except Exception as e:
                    logger.info(str(e))
                    error = 'Unable to assign transcript identity records to ' + accession + \
                            ', potentially an obsolete record or there is an issue retrieving data from NCBI. ' \
                            'Please try again later and if the problem persists contact admin'
                    variant.warnings.append(error)
                    logger.info(error)
                    return True
                variant.description = entry['description']

            # If no correct keys are found
            else:
                # Open a hgvs exception log file in append mode
                error = 'Unknown error type'
                variant.warnings.extend([error, ': A Database error occurred, please contact admin'])
                logger.warning(error)
                return True
        return False

# <LICENSE>
# Copyright (C) 2019 VariantValidator Contributors
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
