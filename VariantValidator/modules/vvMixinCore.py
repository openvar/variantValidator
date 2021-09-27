import vvhgvs
import vvhgvs.exceptions
import vvhgvs.normalizer
import re
import copy
import sys
import logging
import json
from vvhgvs.assemblymapper import AssemblyMapper
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

                    # Select LRG equivalent transcripts
                    if 'LRG' in trans_id:
                        trans_id = self.db.get_refseq_transcript_id_from_lrg_transcript_id(trans_id)
                        if trans_id == 'none':
                            continue

                    # Create dictionaries
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
            flag : mitochondrial
            """

            logger.debug("Batch list length " + str(len(self.batch_list)))
            for my_variant in self.batch_list:

                # Create Normalizers
                my_variant.hn = vvhgvs.normalizer.Normalizer(self.hdp,
                                                             cross_boundaries=False,
                                                             shuffle_direction=3,
                                                             alt_aln_method=self.alt_aln_method
                                                             )
                my_variant.reverse_normalizer = vvhgvs.normalizer.Normalizer(self.hdp,
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

                    # Remove whitespace and quotes
                    my_variant.remove_whitespace()
                    my_variant.remove_quotes()

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
                            my_variant.selected_assembly = selected_assembly
                        else:
                            my_variant.primary_assembly = 'GRCh38'
                            my_variant.selected_assembly = selected_assembly
                            primary_assembly = 'GRCh38'
                            my_variant.warnings.append('Invalid genome build has been specified. Automap has selected '
                                                       'the default build (GRCh38)')
                            logger.warning(
                                'Invalid genome build has been specified. Automap has selected the '
                                'default build ' + my_variant.primary_assembly)
                    else:
                        primary_assembly = my_variant.primary_assembly
                    logger.debug("Completed string formatting")

                    try:
                        toskip = format_converters.initial_format_conversions(my_variant, self,
                                                                              select_transcripts_dict_plus_version)

                    except vvhgvs.exceptions.HGVSError as e:
                        checkref = str(e)
                        try:
                            # Test intronic variants for incorrect boundaries (see issue #169)
                            test_variant = copy.copy(my_variant)
                            test_variant.hgvs_formatted = str(my_variant.quibble)

                            # Create easy variant mapper (over variant mapper) and splign locked evm
                            test_variant.evm = AssemblyMapper(self.hdp,
                                                              assembly_name=primary_assembly,
                                                              alt_aln_method=self.alt_aln_method,
                                                              normalize=True,
                                                              replace_reference=True
                                                              )

                            # Setup a reverse normalize instance and non-normalize evm
                            test_variant.no_norm_evm = AssemblyMapper(self.hdp,
                                                                      assembly_name=primary_assembly,
                                                                      alt_aln_method=self.alt_aln_method,
                                                                      normalize=False,
                                                                      replace_reference=True
                                                                      )

                            mappers.transcripts_to_gene(test_variant, self, select_transcripts_dict_plus_version)
                        except mappers.MappersError:
                            my_variant.output_type_flag = 'warning'
                            continue

                        except vvhgvs.exceptions.HGVSParseError as e:
                            my_variant.warnings.append(str(e))
                            logger.warning(str(e))
                            continue

                        # Other issues to collect, for example, the specified position in NC_ does not agree with g.
                        # See issue #176
                        except Exception:
                            if 'does not agree with reference sequence' in checkref:
                                my_variant.warnings.append(str(e))
                                logger.warning(str(e))
                                continue

                        if 'base start position must be <= end position' in str(e):
                            toskip = None
                        else:
                            my_variant.warnings.append(str(e))
                            if "The entered coordinates do not agree with the intron/exon boundaries for the selected " \
                               "transcript" not in my_variant.warnings[0]:
                                my_variant.warnings.reverse()
                            logger.warning(str(e))
                            continue

                    if toskip:
                        continue

                    # INITIAL USER INPUT FORMATTING
                    # Requested warnings from https://github.com/openvar/variantValidator/issues/195
                    if re.search(r'\(.+?\)', my_variant.quibble):  # Pattern looks for (....)
                        gene_symbol_query = re.search(r'\(.+?\)', my_variant.quibble).group(0)
                        gene_symbol_query = gene_symbol_query.replace('(', '')
                        gene_symbol_query = gene_symbol_query.replace(')', '')
                        is_it_a_gene = self.db.get_hgnc_symbol(gene_symbol_query)
                        if is_it_a_gene != 'none':
                            warning = "Removing redundant gene symbol %s from variant description" % is_it_a_gene
                            my_variant.warnings.append(warning)
                            logger.warning(warning)
                    if re.search('del[GATC]+', my_variant.quibble) or re.search('inv[GATC]+', my_variant.quibble) or\
                            re.search('dup[GATC]+', my_variant.quibble):
                        warning = "Removing redundant reference bases from variant description"
                        my_variant.warnings.append(warning)
                        logger.warning(warning)

                    invalid = my_variant.format_quibble()
                    if invalid:
                        if re.search(r'\w+:[gcnmrp],', my_variant.quibble):
                            error = 'Variant description ' + my_variant.quibble + ' contained the , character between '\
                                    '<type> and <position> in the expected pattern <accession>:<type>.<position> and ' \
                                    'has been auto-corrected'
                            my_variant.quibble = my_variant.quibble.replace(',', '.')
                            my_variant.warnings.append(error)
                            logger.warning(error)
                            pass
                        elif re.search(r'\w+:[gcnmrp]', my_variant.quibble) and not \
                                re.search(r'\w+:[gcnmrp]\.', my_variant.quibble):
                            error = 'Variant description ' + my_variant.quibble + ' lacks the . character between ' \
                                    '<type> and <position> in the expected pattern <accession>:<type>.<position>'
                            my_variant.warnings.append(error)
                            logger.warning(error)
                            continue
                        else:
                            error = 'Variant description ' + my_variant.quibble + ' is not in an accepted format'
                            my_variant.warnings.append(error)
                            logger.warning(error)
                            continue

                    formatted_variant = my_variant.quibble
                    stash_input = my_variant.quibble
                    my_variant.post_format_conversion = stash_input

                    logger.debug("Variant input formatted, proceeding to validate.")

                    # Conversions
                    # Conversions are not currently supported. The HGVS format for conversions
                    # is rarely seen wrt genomic sequencing data and needs to be re-evaluated
                    if 'con' in my_variant.quibble:
                        my_variant.warnings.append('Gene conversions currently unsupported')
                        logger.warning('Gene conversions currently unsupported')
                        continue

                    # Change RNA bases to upper case but nothing else
                    if my_variant.reftype == ":r.":
                        formatted_variant = formatted_variant.upper()
                        formatted_variant = formatted_variant.replace(':R.', ':r.')
                        # lowercase the supported variant types
                        formatted_variant = formatted_variant.replace('DEL', 'del')
                        formatted_variant = formatted_variant.replace('INS', 'ins')
                        formatted_variant = formatted_variant.replace('INV', 'inv')
                        formatted_variant = formatted_variant.replace('DUP', 'dup')

                    # Handle <position><edit><position> style variants
                    # Refer to https://github.com/openvar/variantValidator/issues/161
                    # Example provided is NC_000017.10:g.41199848_41203626delins41207680_41207915
                    # Current theory, should apply to delins, ins.
                    # We may also need to expand to http://varnomen.hgvs.org/recommendations/DNA/variant/insertion/
                    # complex insertions

                    # Validate syntax of the now HGVS variants
                    try:
                        input_parses = self.hp.parse_hgvs_variant(str(formatted_variant))
                        my_variant.hgvs_formatted = input_parses
                    except vvhgvs.exceptions.HGVSError as e:
                        # Look for T not U!
                        posedit = formatted_variant.split(':')[-1]
                        if 'T' in posedit and "r." in posedit:
                            e = 'The IUPAC RNA alphabet dictates that RNA variants must use the character u in ' \
                                'place of t'
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
                        except vvhgvs.exceptions.HGVSError as e:
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
                    try:
                        toskip = use_checking.structure_checks(my_variant, self)
                    except vvhgvs.exceptions.HGVSDataNotAvailableError as e:
                        my_variant.warnings.append(str(e))
                        continue
                    if toskip:
                        continue
                    logger.debug("Variant structure and contents searches passed")

                    # Mitochondrial variants
                    toskip = format_converters.mitochondrial(my_variant, self)
                    if toskip:
                        continue

                    # Protein variants
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
                        toskip = mappers.gene_to_transcripts(my_variant, self, select_transcripts_dict)
                        if toskip:
                            continue

                    if my_variant.reftype == ':c.' or my_variant.reftype == ':n.':
                        try:
                            toskip = mappers.transcripts_to_gene(my_variant, self, select_transcripts_dict_plus_version)
                        except mappers.MappersError:
                            my_variant.output_type_flag = 'warning'
                            continue
                        except vvhgvs.exceptions.HGVSInvalidVariantError:
                            my_variant.output_type_flag = 'warning'
                            continue
                        if toskip:
                            continue

                    # Set the data
                    my_variant.output_type_flag = 'gene'
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
                logger.debug("Formatting variant " + variant.quibble)
                if not variant.write:
                    continue

                # Genomic sequence variation
                # Check for gapped delins<No_Alt>
                if re.search('ins$', variant.genomic_g) and 'del' in variant.genomic_g:
                    variant.genomic_g = variant.genomic_g.replace('ins', '')

                genomic_variant = variant.genomic_g
                hgvs_genomic_variant = genomic_variant

                # genomic accession
                logger.debug("genomic accession")
                if genomic_variant != '':
                    str_genomic_variant = fn.remove_reference_string(genomic_variant)
                    hgvs_genomic_variant = self.hp.parse_hgvs_variant(str_genomic_variant)
                    genomic_variant = fn.valstr(hgvs_genomic_variant)
                    genomic_accession = hgvs_genomic_variant.ac
                else:
                    genomic_accession = ''

                # RefSeqGene variation
                logger.debug("RefSeqGene variation")
                refseqgene_variant = variant.genomic_r
                refseqgene_variant = refseqgene_variant.strip()
                if 'RefSeqGene' in refseqgene_variant or refseqgene_variant == '':
                    variant.warnings.append(refseqgene_variant)
                    refseqgene_variant = ''
                    lrg_variant = ''
                    hgvs_refseqgene_variant = 'false'
                else:
                    hgvs_refseqgene_variant = self.hp.parse_hgvs_variant(fn.remove_reference_string(refseqgene_variant))
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
                logger.debug("Transcript sequence variation")
                tx_variant = variant.coding
                hgvs_transcript_variant = tx_variant
                hgvs_tx_variant = None
                if tx_variant != '':
                    if '(' in tx_variant and ')' in tx_variant:
                        tx_variant = tx_variant.split('(')[1]
                        tx_variant = tx_variant.replace(')', '')

                    # transcript accession
                    logger.debug("transcript accession")
                    hgvs_tx_variant = self.hp.parse_hgvs_variant(fn.remove_reference_string(tx_variant))
                    tx_variant = fn.valstr(hgvs_tx_variant)
                    hgvs_transcript_variant = self.hp.parse_hgvs_variant(fn.remove_reference_string(tx_variant))
                    transcript_accession = hgvs_transcript_variant.ac

                    # Handle LRG
                    logger.debug("Handle LRG")
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
                logger.debug("Look for intronic variants")
                if transcript_accession != '' and genomic_accession != '':
                    # Remove del bases
                    str_transcript = fn.valstr(hgvs_transcript_variant)
                    hgvs_transcript_variant = self.hp.parse_hgvs_variant(str_transcript)
                    try:
                        self.vr.validate(hgvs_transcript_variant)
                    except vvhgvs.exceptions.HGVSError as e:
                        error = str(e)
                        if 'intronic variant' in error:
                            genome_context_transcript_variant = genomic_accession + '(' + transcript_accession +\
                                                                '):c.' + str(hgvs_transcript_variant.posedit)
                            if refseqgene_variant != '':
                                hgvs_refseqgene_variant = self.hp.parse_hgvs_variant(
                                    fn.remove_reference_string(refseqgene_variant))
                                refseqgene_accession = hgvs_refseqgene_variant.ac
                                hgvs_coding_from_refseqgene = self.vm.g_to_t(hgvs_refseqgene_variant,
                                                                             hgvs_transcript_variant.ac)
                                hgvs_coding_from_refseqgene = fn.valstr(hgvs_coding_from_refseqgene)
                                hgvs_coding_from_refseqgene = self.hp.parse_hgvs_variant(
                                    fn.remove_reference_string(hgvs_coding_from_refseqgene))
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
                logger.debug("Protein description")
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

                # Identify Pseudo Autosomal Regions
                chrX = False
                chrY = False
                par = False
                for g_var in multi_gen_vars:
                    if 'NC_000023' in g_var.ac:
                        chrX = True
                    if 'NC_000024' in g_var.ac:
                        chrY = True
                if chrX is True and chrY is True:
                    par = True

                for alt_gen_var in multi_gen_vars:
                    if 'NC_' in alt_gen_var.ac:
                        if 'NC_000' not in alt_gen_var.ac and 'NC_012920.1' not in alt_gen_var.ac and \
                                'NC_001807.4' not in alt_gen_var.ac:
                            continue
                    try:
                        alt_gen_var = variant.hn.normalize(alt_gen_var)
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        continue
                    except vvhgvs.exceptions.HGVSDataNotAvailableError:
                        continue

                    for build in self.genome_builds:
                        test = seq_data.supported_for_mapping(alt_gen_var.ac, build)
                        if test:
                            try:
                                vcf_dict = hgvs_utils.report_hgvs2vcf(alt_gen_var, build, variant.reverse_normalizer,
                                                                      self.sf)
                            except vvhgvs.exceptions.HGVSInvalidVariantError:
                                continue
                            # Identify primary assembly positions
                            if 'NC_' in alt_gen_var.ac and par is False:
                                if 'NC_000' not in alt_gen_var.ac and 'NC_012920.1' not in alt_gen_var.ac and \
                                        'NC_001807.4' not in alt_gen_var.ac:
                                    continue
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

                            elif 'NC_' not in alt_gen_var.ac and par is False:
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
                                    vcf_dict = hgvs_utils.report_hgvs2vcf(alt_gen_var, 'hg38',
                                                                          variant.reverse_normalizer,
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

                            # Identify primary assembly positions
                            elif 'NC_000023' in alt_gen_var.ac and par is True:
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
                                    vcf_dict = hgvs_utils.report_hgvs2vcf(alt_gen_var, 'hg38',
                                                                          variant.reverse_normalizer,
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

                # Add predicted protein variant dictionary
                if predicted_protein_variant != '':
                    predicted_protein_variant_dict = {}
                    predicted_protein_variant_dict["slr"] = ''
                    predicted_protein_variant_dict["tlr"] = ''
                    predicted_protein_variant_dict["lrg_tlr"] = ''
                    predicted_protein_variant_dict["lrg_slr"] = ''
                    if 'Non-coding :n.' not in predicted_protein_variant:
                        try:
                            # Note this code is needed if we decide to come in line with Mutalyzer  - see issue #214
                            # Requires commenting out of issue #214 code in MixinInit

                            # # remove trailing aas after ter in insertions and delins
                            # if 'fs' not in predicted_protein_variant \
                            #         and 'delins' in predicted_protein_variant\
                            #         and 'ext' not in predicted_protein_variant:
                            #     predicted_protein_variant = re.sub(r'Ter\w+', 'Ter', predicted_protein_variant)
                            #     predicted_protein_variant_dict["tlr"] = predicted_protein_variant
                            #
                            #     # Remove Training Ter from delins except for ext proteins
                            #     if re.search(r'delins\w+Ter', predicted_protein_variant):
                            #         format_p = predicted_protein_variant.replace('Ter', '')
                            #         format_p = re.sub(r'\(LRG_.+?\)', '', format_p)
                            #         re_parse_protein = self.hp.parse_hgvs_variant(format_p)
                            #         p_seq = self.sf.fetch_seq(format_p.split(':')[0])
                            #
                            #         end_aa = fn.one_to_three((p_seq[-1]))
                            #         p_len = len(p_seq)
                            #         prot_st, posedit = format_p.split('delins')
                            #         prot_sta, prot_aa_st = prot_st.split(':p.')
                            #         prot_st = prot_sta + ':p.'
                            #         prot_aa_st = prot_aa_st.split('_')[0]
                            #         prot_st = prot_st + prot_aa_st
                            #
                            #         # Create edit
                            #         if re_parse_protein.posedit.pos.start.base != p_len:
                            #             pre_posedit = '_%s%sdelins' % (end_aa, str(p_len))
                            #         else:
                            #             pre_posedit = 'delins'
                            #         posedit = pre_posedit + posedit
                            #
                            #         # Assemble protein variants
                            #         predicted_protein_variant_dict["tlr"] = '%s%s%s%s' % (
                            #             predicted_protein_variant.split(':p.')[0],
                            #             ":p.",
                            #             prot_st.split(':p.')[1],
                            #             posedit)
                            #
                            #         predicted_protein_variant = prot_st + posedit

                            # Add single letter AA code to protein descriptions
                            predicted_protein_variant_dict = {"tlr": str(predicted_protein_variant), "slr": ''}
                            if re.search('p.=', predicted_protein_variant_dict['tlr']) \
                                    or re.search('p.?', predicted_protein_variant_dict['tlr']):
                                # Replace p.= with p.(=)
                                predicted_protein_variant_dict['tlr'] = predicted_protein_variant_dict['tlr'].replace(
                                    'p.=',
                                    'p.(=)')

                            # Remove LRG
                            format_p = predicted_protein_variant_dict['tlr']
                            if 'LRG' in format_p:
                                format_lrg = copy.copy(format_p)
                                format_p = re.sub(r'\(LRG_.+?\)', '', format_p)
                                format_lrg = format_lrg.split('(')[1]
                                format_lrg = format_lrg.replace(')', '')
                            else:
                                format_lrg = None
                                pass

                            re_parse_protein = self.hp.parse_hgvs_variant(format_p)

                            # Set formatted tlr
                            predicted_protein_variant_dict['tlr'] = str(copy.copy(re_parse_protein))
                            re_parse_protein_single_aa = fn.single_letter_protein(re_parse_protein)

                            # Replace p.= with p.(=)
                            if re.search('p.=', re_parse_protein_single_aa
                                         ) or re.search('p.?', re_parse_protein_single_aa):
                                re_parse_protein_single_aa = re_parse_protein_single_aa.replace('p.=',
                                                                                                'p.(=)')

                            predicted_protein_variant_dict["slr"] = str(re_parse_protein_single_aa)

                            # set LRG outputs
                            if format_lrg is not None:
                                predicted_protein_variant_dict["lrg_tlr"] = \
                                    format_lrg.split(':')[0] + ':' + predicted_protein_variant_dict["tlr"].split(':')[1]
                                predicted_protein_variant_dict["lrg_slr"] = \
                                    format_lrg.split(':')[0] + ':' + predicted_protein_variant_dict["slr"].split(':')[1]
                            else:
                                predicted_protein_variant_dict["lrg_tlr"] = ''
                                predicted_protein_variant_dict["lrg_slr"] = ''

                        except vvhgvs.exceptions.HGVSParseError as e:
                            logger.debug("Except passed, %s", e)
                else:
                    predicted_protein_variant_dict = {}
                    predicted_protein_variant_dict["slr"] = ''
                    predicted_protein_variant_dict["tlr"] = ''
                    predicted_protein_variant_dict["lrg_tlr"] = ''
                    predicted_protein_variant_dict["lrg_slr"] = ''

                # Add stable gene_ids
                stable_gene_ids = {}
                if variant.gene_symbol != '':
                    gene_stable_info = self.db.get_stable_gene_id_info(variant.gene_symbol)

                    # Add or update stable ID and transcript data
                    if gene_stable_info[1] == 'No data' and hgvs_tx_variant is not None:
                        try:
                            self.db.update_transcript_info_record(hgvs_tx_variant.ac, self)
                        except fn.DatabaseConnectionError as e:
                            error = 'Currently unable to update all gene_ids or transcript information records ' \
                                    'because ' \
                                    'VariantValidator %s' % str(e)
                            my_variant.warnings.append(error)
                            logger.warning(error)
                            try:
                                self.db.update_transcript_info_record(hgvs_tx_variant.ac, self,
                                                                      bypass_with_symbol=variant.gene_symbol)
                            except fn.DatabaseConnectionError as e:
                                error = 'Currently unable to update all gene_ids or transcript information records ' \
                                        'because ' \
                                        'VariantValidator %s' % str(e)
                                my_variant.warnings.append(error)
                                logger.warning(error)

                        gene_stable_info = self.db.get_stable_gene_id_info(variant.gene_symbol)

                    # Update gene_symbol
                    if variant.gene_symbol != str(gene_stable_info[1]) and str(gene_stable_info[1]) != 'No data':
                        variant.gene_symbol = str(gene_stable_info[1])

                    try:
                        # Dictionary the output
                        stable_gene_ids['hgnc_id'] = gene_stable_info[2]
                        stable_gene_ids['entrez_gene_id'] = gene_stable_info[3]
                        stable_gene_ids['ensembl_gene_id'] = gene_stable_info[4]
                        stable_gene_ids['ucsc_id'] = gene_stable_info[5]
                        stable_gene_ids['omim_id'] = json.loads(gene_stable_info[6])
                        # stable_gene_ids['vega_id'] = gene_stable_info[7]

                        # reformat ccds return into a Python list
                        my_ccds = gene_stable_info[8].replace('[', '')
                        my_ccds = my_ccds.replace(']', '')
                        my_ccds = my_ccds.replace('"','')
                        my_ccds = my_ccds.replace(',', '')
                        ccds_list = my_ccds.split()
                        stable_gene_ids['ccds_ids'] = ccds_list

                    except IndexError as e:
                        logger.debug("Except pass, %s", e)

                # Add Reference sequence annotations
                reference_annotations = {}
                if hgvs_tx_variant is not None:
                    annotation_info = self.db.get_transcript_annotation(hgvs_tx_variant.ac)
                    # Add or update stable ID and transcript data
                    try:
                        annotation_info = json.loads(annotation_info)
                        annotation_info.keys()
                    except Exception:
                        try:
                            self.db.update_transcript_info_record(hgvs_tx_variant.ac, self)
                        except fn.DatabaseConnectionError as e:
                            error = 'Currently unable to update all gene_ids or transcript information records ' \
                                    'because ' \
                                    'VariantValidator %s' % str(e)
                            my_variant.warnings.append(error)
                            logger.warning(error)
                    annotation_info = self.db.get_transcript_annotation(hgvs_tx_variant.ac)
                    reference_annotations = json.loads(annotation_info)

                variant.stable_gene_ids = stable_gene_ids
                variant.annotations = reference_annotations
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
                if (variant.output_type_flag == 'intergenic') or ('grch37' not in variant.primary_assembly_loci.keys())\
                        or ('grch38' not in variant.primary_assembly_loci.keys()):

                    # Simple cache
                    lo_cache = {}

                    # Attempt to liftover between genome builds
                    # Note: pyliftover uses the UCSC liftOver tool.
                    # https://pypi.org/project/pyliftover/
                    genomic_position_info = variant.primary_assembly_loci
                    for g_p_key in list(genomic_position_info.keys()):
                        build_to = ''
                        build_from = ''

                        # Identify the current build and hgvs_genomic description
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

                        # Genome to Genome liftover if tx not annotated to the build
                        g_to_g = False
                        if variant.output_type_flag != 'intergenic':
                            g_to_g = True
                        # Liftover
                        if genomic_position_info[g_p_key]['hgvs_genomic_description'] not in lo_cache.keys():
                            lifted_response = liftover(genomic_position_info[g_p_key]['hgvs_genomic_description'],
                                                       build_from,
                                                       build_to, variant.hn, variant.reverse_normalizer,
                                                       variant.evm, self, specify_tx=False, liftover_level=False,
                                                       g_to_g=g_to_g)
                            lo_cache[genomic_position_info[g_p_key]['hgvs_genomic_description']] = lifted_response
                        else:
                            lifted_response = lo_cache[genomic_position_info[g_p_key]['hgvs_genomic_description']]

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

                # Remove duplicate warnings
                variant_warnings = []
                accession = variant.hgvs_transcript_variant.split(':')[0]
                term = "(" + accession + ")"
                for vt in variant.warnings:
                    #  Do not warn a transcript update is available for the most recent transcript
                    if term in vt and "A more recent version of the selected reference sequence" in vt:
                        continue
                    else:
                        variant_warnings.append(vt)
                variant.warnings = variant_warnings

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

        # Search for gene symbol on Transcript inputs
        hgnc = query
        if 'NM_' in hgnc or 'NR_' in hgnc:  # or re.match('ENST', hgnc):

            # Remove version
            if '.' in hgnc:
                hgnc = hgnc.split('.')[0]

            # Find latest version in UTA
            found_res = False
            for version in range(25):
                refresh_hgnc = hgnc + '.' + str(version)
                try:
                    self.hdp.get_tx_identity_info(refresh_hgnc)
                    tx_found = refresh_hgnc
                    found_res = True
                    break
                except vvhgvs.exceptions.HGVSError as e:
                    logger.debug("Except passed, %s", e)
            if not found_res:
                return {'error': 'No transcript definition for (tx_ac=' + hgnc + ')'}

            # update record and correct symbol
            try:
                self.db.update_transcript_info_record(tx_found, self)
            except fn.DatabaseConnectionError as e:
                error = 'Currently unable to update gene_ids or transcript information records because ' \
                        'VariantValidator %s' % str(e)
                # my_variant.warnings.append(error)
                logger.warning(error)

            try:
                tx_info = self.hdp.get_tx_identity_info(tx_found)
            except vvhgvs.exceptions.HGVSError as e:
                return {'error': str(e)}
            hgnc = tx_info[6]
            hgnc = self.db.get_hgnc_symbol(hgnc)

        # First perform a search against the input gene symbol or the symbol inferred from UTA
        symbol_identified = False
        vvta_record = self.hdp.get_gene_info(hgnc)
        # Check for a record
        if vvta_record is not None:
            current_sym = hgnc
            gene_name = vvta_record[3]
            hgnc_id = vvta_record[0]
            previous_sym = vvta_record[5]
            symbol_identified = True

        # No record found, is it a previous symbol?
        else:
            # Look up current name
            vvta_record = self.hdp.get_gene_info_by_alias(hgnc)
            if vvta_record is not None:
                if len(vvta_record) == 1:
                    current_sym = vvta_record[0][1]
                    gene_name = vvta_record[0][3]
                    hgnc_id = vvta_record[0][0]
                    previous_sym = hgnc
                    symbol_identified = True
                if len(vvta_record) > 1:
                    return {'error': '%s is a previous symbol for %s genes. '
                            'Refer to https://www.genenames.org/' % (current_sym, str(len(vvta_record)))}

        if symbol_identified is False:
            return {'error': 'Unable to recognise gene symbol %s' % hgnc}

        # Get transcripts
        tx_for_gene = self.hdp.get_tx_for_gene(current_sym)
        if len(tx_for_gene) == 0:
            tx_for_gene = self.hdp.get_tx_for_gene(previous_sym)
        if len(tx_for_gene) == 0:
            return {'error': 'Unable to retrieve data from the VVTA, please contact admin'}

        # Loop through each transcript and get the relevant transcript description
        genes_and_tx = []
        recovered = []
        for line in tx_for_gene:
            if (line[3].startswith('NM_') or line[3].startswith('NR_')) and '..' not in line[3]:
                # Transcript ID
                tx = line[3]

                # Protein id
                prot_id = self.hdp.get_pro_ac_for_tx_ac(tx)

                # Get additional tx_ information
                tx_exons = self.hdp.get_tx_exons(tx, line[4], line[5])
                tx_orientation = tx_exons[0]['alt_strand']

                # Fetch the sequence to get the length
                tx_seq = self.sf.fetch_seq(tx)
                tx_len = len(tx_seq)

                # Collect genomic span for the transcript against known genomic/gene reference sequences
                gen_start_pos = None
                gen_end_pos = None
                exon_set = []
                # get total exons
                total_exons = len(tx_exons)
                # Set exon counter for current exon
                current_exon_number = 0
                for tx_pos in tx_exons:
                    current_exon_number = current_exon_number + 1
                    # Collect the exon_set information
                    """
                    tx_exons have the following attributes::
                                {
                                    'tes_exon_set_id' : 98390
                                    'aes_exon_set_id' : 298679
                                    'tx_ac'           : 'NM_199425.2'
                                    'alt_ac'          : 'NC_000020.10'
                                    'alt_strand'      : -1
                                    'alt_aln_method'  : 'splign'
                                    'ord'             : 2
                                    'tx_exon_id'      : 936834
                                    'alt_exon_id'     : 2999028
                                    'tx_start_i'      : 786
                                    'tx_end_i'        : 1196
                                    'alt_start_i'     : 25059178
                                    'alt_end_i'       : 25059588
                                    'cigar'           : '410='
                                }                    
                    """
                    current_exon = {"transcript_start": tx_pos['tx_start_i'] + 1,
                                    "transcript_end": tx_pos['tx_end_i'],
                                    "genomic_start": tx_pos['alt_start_i'] + 1,
                                    "genomic_end": tx_pos['alt_end_i'],
                                    "cigar": tx_pos['cigar'],
                                    "exon_number": current_exon_number
                                    #"total_exons": total_exons
                                    }
                    exon_set.append(current_exon)
                    start_pos = tx_pos['alt_start_i']
                    end_pos = tx_pos['alt_end_i']
                    if gen_start_pos is None:
                        gen_start_pos = start_pos
                    else:
                        if int(start_pos) < int(gen_start_pos):
                            gen_start_pos = int(start_pos)
                    if gen_end_pos is None:
                        gen_end_pos = end_pos
                    else:
                        if int(end_pos) > int(gen_end_pos):
                            gen_end_pos = int(end_pos)

                # reverse the exon_set to maintain gene and not genome orientation if gene is -1 orientated
                if tx_orientation == -1:
                    exon_set.reverse()

                if ('NG_' in line[4] or 'NC_000' in line[4]) and line[5] != 'blat':
                    gen_span = True
                else:
                    gen_span = False

                tx_description = self.db.get_transcript_description(tx)
                if tx_description == 'none':
                    try:
                        self.db.update_transcript_info_record(tx, self)
                    except fn.DatabaseConnectionError as e:
                        error = 'Currently unable to update gene_ids or transcript information records because ' \
                                'VariantValidator %s' % str(e)
                        # my_variant.warnings.append(error)
                        logger.warning(error)
                    tx_description = self.db.get_transcript_description(tx)

                # Get annotation
                tx_annotation = self.db.get_transcript_annotation(tx)
                tx_annotation = json.loads(tx_annotation)

                # Check for duplicates
                if tx not in recovered:
                    recovered.append(tx)
                    if len(line) >= 3 and isinstance(line[1], int):
                        genes_and_tx.append({'reference': tx,
                                             'description': tx_description,
                                             'annotations': tx_annotation,
                                             'translation': prot_id,
                                             'length': tx_len,
                                             'coding_start': line[1] + 1,
                                             'coding_end': line[2],
                                             # 'orientation': tx_orientation,
                                             'genomic_spans': {}
                                             })
                    else:
                        genes_and_tx.append({'reference': tx,
                                             'description': tx_description,
                                             'annotations': tx_annotation,
                                             'translation': prot_id,
                                             'length': tx_len,
                                             'coding_start': None,
                                             'coding_end': None,
                                             # 'orientation': tx_orientation,
                                             'genomic_spans': {}
                                             })
                    # LRG information
                    lrg_transcript = self.db.get_lrg_transcript_id_from_refseq_transcript_id(tx)
                    if lrg_transcript != 'none':
                        genes_and_tx.append({'reference': lrg_transcript,
                                             'description': tx_description,
                                             'annotations': tx_annotation,
                                             'length': tx_len,
                                             'translation': lrg_transcript.replace('t', 'p'),
                                             'coding_start': line[1] + 1,
                                             'coding_end': line[2],
                                             'genomic_spans': {}
                                             })

                # Add the genomic span information
                if gen_span is True:
                    for check_tx in genes_and_tx:
                        lrg_transcript = self.db.get_lrg_transcript_id_from_refseq_transcript_id(tx)
                        if check_tx['reference'] == tx:
                            if gen_start_pos < gen_end_pos:
                                check_tx['genomic_spans'][line[4]] = {'start_position': gen_start_pos + 1,
                                                                      'end_position': gen_end_pos,
                                                                      'orientation': tx_orientation,
                                                                      'exon_structure': exon_set,
                                                                      "total_exons": total_exons}
                            else:
                                check_tx['genomic_spans'][line[4]] = {'start_position': gen_end_pos + 1,
                                                                      'end_position': gen_start_pos,
                                                                      'orientation': tx_orientation,
                                                                      'exon_structure': exon_set,
                                                                      "total_exons": total_exons}
                        if lrg_transcript != 'none':
                            if check_tx['reference'] == lrg_transcript:
                                if 'NG_' in line[4]:
                                    lrg_id = self.db.get_lrg_id_from_refseq_gene_id(line[4])
                                    if lrg_id[0] in lrg_transcript:
                                        check_tx['genomic_spans'][line[4]] = {'start_position': gen_start_pos + 1,
                                                                              'end_position': gen_end_pos,
                                                                              'orientation': 1,
                                                                              'exon_structure': exon_set,
                                                                              "total_exons": total_exons}

                                        check_tx['genomic_spans'][lrg_id[0]] = {'start_position': gen_start_pos + 1,
                                                                                'end_position': gen_end_pos,
                                                                                'orientation': 1,
                                                                                'exon_structure': exon_set,
                                                                                "total_exons": total_exons}
        # Return data table
        g2d_data = {'current_symbol': current_sym,
                    'previous_symbol': previous_sym,
                    'current_name': gene_name,
                    'hgnc': hgnc_id,
                    # 'previous_name': previous_name,
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
        logger.debug("Looking for transcript info")
        hgvs_vt = self.hp.parse_hgvs_variant(str(variant.hgvs_formatted))
        try:
            self.hdp.get_tx_identity_info(str(hgvs_vt.ac))
        except vvhgvs.exceptions.HGVSError as e:
            error = 'Please inform admin of the following error: ' + str(e)
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
                    except vvhgvs.exceptions.HGVSError:
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
                        # If the none key is found add the description to the database
                        if 'UTA' in str(e):
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
                    # Allows bypass with current record if external databases not available
                    if 'UTA' in str(e):
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

    def update_transcript_record(self, tx_id, **kwargs):
        """
        Siplle function allowing transcript_table to be updated
        :param tx_id:
        :param genome_build (GRCh37 or GRCh38)
        :return:
        """
        self.db.update_transcript_info_record(tx_id, self, **kwargs)

# <LICENSE>
# Copyright (C) 2016-2021 VariantValidator Contributors
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
