import hgvs
import hgvs.exceptions
import hgvs.normalizer
import re
import copy
import sys
import traceback
from hgvs.assemblymapper import AssemblyMapper
from .vvLogging import logger
from . import vvHGVS
from . import vvFunctions as fn
from . import vvChromosomes
from . import vvMixinConverters
from .variant import Variant
from . import format_converters
from . import use_checking
from . import collect_info
from . import mappers
from . import valoutput


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
        logger.info(batch_variant + ' : ' + selected_assembly)

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
                        trans_id = self.db.get_RefSeqTranscriptID_from_lrgTranscriptID(trans_id)
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
                # Start timing
                logger.traceStart(my_variant)

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
                            my_variant.warnings += ': Invalid genome build has been specified. ' \
                                                   'Automap has selected the default build (GRCh38)'
                            logger.warning(
                                'Invalid genome build has been specified. Automap has selected the '
                                'default build ' + my_variant.primary_assembly)
                    else:
                        primary_assembly = my_variant.primary_assembly
                    logger.trace("Completed string formatting", my_variant)

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
                        my_variant.warnings += ': ' + error
                        logger.warning(error)
                        continue

                    formatted_variant = my_variant.quibble
                    stash_input = my_variant.quibble
                    my_variant.stashed = stash_input
                    format_type = my_variant.reftype

                    hgnc_gene_info = 'false'

                    logger.trace("Variant input formatted, proceeding to validate.", my_variant)

                    # Conversions
                    # Conversions are not currently supported. The HGVS format for conversions
                    # is rarely seen wrt genomic sequencing data and needs to be re-evaluated
                    if 'con' in my_variant.quibble:
                        my_variant.warnings += ': ' + 'Gene conversions currently unsupported'
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
                        my_variant.warnings += ': ' + str(e)
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
                            my_variant.warnings += ': ' + error
                            logger.warning(error)
                            continue
                        else:
                            my_variant.warnings += ': ' + str(trap_ens_in) + ' automapped to equivalent ' \
                                                                             'RefSeq transcript ' + my_variant.quibble
                            logger.warning(str(trap_ens_in) + ' automapped to equivalent RefSeq '
                                                              'transcript ' + my_variant.quibble)
                    logger.trace("HVGS acceptance test passed", my_variant)

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
                        my_variant.warnings += ': ' + error
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
                        except hgvs.exceptions.HGVSError:
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

                    # Catch missing version number in refseq
                    is_version = re.compile(r"\d\.\d")
                    if my_variant.refsource == 'RefSeq' and not is_version.search(str(my_variant.hgvs_formatted)):
                        error = 'RefSeq variant accession numbers MUST include a version number'
                        my_variant.warnings += ': ' + str(error)
                        continue
                    logger.trace("HVGS interval/version mapping complete", my_variant)

                    # handle LRG inputs

                    if my_variant.refsource == 'LRG':
                        format_converters.lrg_to_refseq(my_variant, self)
                        logger.trace("LRG check for conversion to refseq completed", my_variant)

                    # Additional Incorrectly input variant capture training
                    if my_variant.refsource == 'RefSeq':
                        toskip = use_checking.refseq_common_mistakes(my_variant)
                        if toskip:
                            continue
                        logger.trace("Passed 'common mistakes' catcher", my_variant)

                    # Primary validation of the input
                    toskip = use_checking.structure_checks(my_variant, self)
                    if toskip:
                        continue
                    logger.trace("Variant structure and contents searches passed", my_variant)

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

                    # COLLECT gene symbol, name and ACCESSION INFORMATION
                    # Gene symbol
                    if my_variant.reftype != ':g.':
                        toskip = collect_info.get_transcript_info(my_variant, self)
                        if toskip:
                            continue

                    # Now start mapping from genome to transcripts
                    if my_variant.reftype == ':g.':
                        toskip = mappers.gene_to_transcripts(my_variant, self)
                        if toskip:
                            continue

                    if format_type == ':c.' or format_type == ':n.':
                        toskip = mappers.transcripts_to_gene(my_variant, self)
                        if toskip:
                            continue

                    # Set the data
                    my_variant.output_type_flag = 'gene'
                    my_variant.description = hgnc_gene_info
                    my_variant.primary_assembly = primary_assembly
                    logger.traceEnd(my_variant)
                # Report errors to User and VV admin
                except KeyboardInterrupt:
                    raise
                except Exception:
                    my_variant.output_type_flag = 'error'
                    error = 'Validation error'
                    my_variant.warnings = str(error)
                    exc_type, exc_value, last_traceback = sys.exc_info()
                    te = traceback.format_exc()
                    tbk = [str(exc_type), str(exc_value), str(te)]
                    er = str('\n'.join(tbk))
                    logger.error(str(exc_type) + " " + str(exc_value))
                    logger.debug(er)
                    raise

            # Outside the for loop
            ######################
            logger.trace("End of for loop")
            # order the rows
            by_order = sorted(self.batch_list, key=lambda x: x.order)

            for variant in by_order:
                if not variant.write:
                    continue

                # warngins
                warnings = variant.warnings
                warnings = re.sub('del[GATC][GATC][GATC][GATC]+', 'del', warnings)
                warnings = re.sub('^: ', '', warnings)
                warnings = re.sub('::', ':', warnings)

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
                        if rsg_ac[1] != 'public':
                            warnings = warnings + ': The current status of ' + str(
                                hgvs_lrg.ac) + ' is pending therefore changes may be made to the LRG reference sequence'

                # Transcript sequence variation
                tx_variant = variant.coding
                hgvs_transcript_variant = tx_variant
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
                    lrg_transcript = self.db.get_lrgTranscriptID_from_RefSeqTranscriptID(transcript_accession)
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
                    for build in self.genome_builds:
                        test = vvChromosomes.supported_for_mapping(alt_gen_var.ac, build)
                        if test == 'true':
                            try:
                                vcf_dict = vvHGVS.report_hgvs2vcf(alt_gen_var, build, variant.reverse_normalizer,
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
                                    vcf_dict = vvHGVS.report_hgvs2vcf(alt_gen_var, 'hg38', variant.reverse_normalizer,
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
                                    vcf_dict = vvHGVS.report_hgvs2vcf(alt_gen_var, 'hg38', variant.reverse_normalizer,
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
                        warnings = warnings + ': ' + str(
                            variant.hgvs_coding) + ' cannot be mapped directly to genome build ' + primary_assembly + \
                            ': See alternative genomic loci or alternative genome builds for aligned genomic positions'

                warn_list = warnings.split(': ')
                warnings_out = []
                for warning in warn_list:
                    warning.strip()
                    warning = warning.replace("'", "")
                    if warning == '':
                        continue
                    if warning not in warnings_out:
                        # Remove duplicate elements but maintain the order
                        warnings_out.append(warning)

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
                    if 'Non-coding :n.' not in predicted_protein_variant:
                        try:
                            format_p = predicted_protein_variant
                            format_p = re.sub(r'\(LRG_.+?\)', '', format_p)
                            re_parse_protein = self.hp.parse_hgvs_variant(format_p)
                            re_parse_protein_single_aa = fn.single_letter_protein(re_parse_protein)
                            predicted_protein_variant_dict["slr"] = str(re_parse_protein_single_aa)
                        except hgvs.exceptions.HGVSParseError:
                            fn.exceptPass()
                    else:
                        predicted_protein_variant_dict["slr"] = str(predicted_protein_variant)

                variant.gene_symbol = gene_symbol
                variant.hgvs_transcript_variant = tx_variant
                variant.genome_context_intronic_sequence = genome_context_transcript_variant
                variant.refseqgene_context_intronic_sequence = refseqgene_context_transcript_variant
                variant.hgvs_refseqgene_variant = refseqgene_variant
                variant.hgvs_predicted_protein_consequence = predicted_protein_variant_dict
                variant.validation_warnings = warnings_out
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
            raise
