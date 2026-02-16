import vvhgvs
import vvhgvs.exceptions
import vvhgvs.normalizer
from vvhgvs.location import Interval
from vvhgvs.sequencevariant import SequenceVariant
from vvhgvs.posedit import PosEdit
import re
import copy
import sys
import logging
import json
import time
from vvhgvs.assemblymapper import AssemblyMapper
from VariantValidator.modules import hgvs_utils
from VariantValidator.modules import utils as fn
from VariantValidator.modules import vvMixinConverters
from VariantValidator.modules.variant import Variant
from VariantValidator.modules import format_converters
from VariantValidator.modules import use_checking
from VariantValidator.modules import mappers
from VariantValidator.modules import valoutput
from VariantValidator.modules import exon_numbering
from VariantValidator.modules.liftover import liftover
from VariantValidator.modules import gene2transcripts
from VariantValidator.modules import lovd_api
from VariantValidator.modules import initial_formatting
from VariantValidator.modules import vcf_to_pvcf
from VariantValidator.modules.seq_state_to_expanded_repeat import\
        convert_seq_state_to_expanded_repeat
from VariantValidator.modules.hgvs_utils import hgvs_delins_parts_to_hgvs_obj,\
        unset_hgvs_obj_ref, to_vv_hgvs

logger = logging.getLogger(__name__)

class Mixin(vvMixinConverters.Mixin):
    """
    This module contains the main function for variant validator.
    It's added to the Validator object in the vvObjects file.
    """

    def __init__(self):
        super().__init__()
        self.lovd_syntax_check = None

    def validate(self,
                 batch_variant,
                 selected_assembly,
                 select_transcripts,
                 transcript_set=None,
                 liftover_level=False,
                 lovd_syntax_check=False):
        """
        This is the main validator function.
        :param batch_variant: A string containing the variant to be validated
        :param selected_assembly: The version of the genome assembly to use.
        :param select_transcripts: Can be an array of different transcripts, or 'all'
        :param liftover_level: True or None or primary - liftover to different gene/genome builds or not
        Selecting multiple transcripts will lead to a multiple variant outputs.
        :param transcript_set: 'refseq' or 'ensembl'
        :param lovd_syntax_check: True or False
        :return:
        """
        logger.debug("Running validate with inputs %s and assembly %s", batch_variant, selected_assembly)

        if transcript_set == "refseq" or transcript_set is None:
            self.alt_aln_method = 'splign'
        elif transcript_set == "ensembl":
            self.alt_aln_method = 'genebuild'
            # Dangerous to liftover based on Ensembl data?
            liftover_level = None
        else:
            raise Exception("The transcriptSet variable '%s' is invalid, it must be 'refseq' or 'ensembl'" %
                            transcript_set)

        # Set the primary assembly
        primary_assembly = None
        self.selected_assembly = selected_assembly
        self.select_transcripts = select_transcripts

        # Set LOVD syntax checker
        self.lovd_syntax_check = lovd_syntax_check

        # Validation
        ############
        try:
            # Create a dictionary of transcript ID : ''
            select_transcripts_dict = {}
            select_transcripts_dict_plus_version = {}
            if select_transcripts != 'all' and select_transcripts != 'raw':
                try:
                    select_transcripts_list = json.loads(select_transcripts)
                except json.decoder.JSONDecodeError:
                    select_transcripts_list = [select_transcripts]

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
            try:
                batch_queries = json.loads(batch_variant)
            except json.decoder.JSONDecodeError:
                batch_queries = [batch_variant]
            if isinstance(batch_queries, int):
                batch_queries = [str(batch_queries)]
            # Turn each variant into a dictionary. The dictionary will be compiled during validation
            self.batch_list = []
            for queries in batch_queries:
                # try:
                #     queries = vcf_to_pvcf.vcf_to_shorthand(queries)
                # except vcf_to_pvcf.VcfConversionError as e:
                #     logger.info(f"Cannot convert {queries} into PVCF format {e}")
                #     pass
                if isinstance(queries, int):
                    queries = str(queries)
                    queries = str(queries)
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
                # Add also to validator
                my_variant.reverse_normalizer = vvhgvs.normalizer.Normalizer(self.hdp,
                                                                             cross_boundaries=False,
                                                                             shuffle_direction=5,
                                                                             alt_aln_method=self.alt_aln_method
                                                                             )
                my_variant.cross_hn = vvhgvs.normalizer.Normalizer(self.hdp,
                                                                   cross_boundaries=True,
                                                                   shuffle_direction=3,
                                                                   alt_aln_method=self.alt_aln_method
                                                                   )

                # This will be used to order the final output
                if not my_variant.order:
                    ordering = ordering + 1
                    my_variant.order = ordering

                if type(my_variant.quibble) is not str:
                    # already tidied input mapped to transcripts so no need to re-validate for user input type issues
                    if not my_variant.hgvs_formatted:
                        my_variant.hgvs_formatted = my_variant.quibble
                    if not my_variant.reftype:
                        my_variant.reftype = f':{my_variant.quibble.type}.'
                    if my_variant.reftype != ':g.':
                        toskip = self._get_transcript_info(my_variant)
                        if toskip:
                            continue
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

                    if my_variant.reftype in [':c.', ':n.']:
                        my_variant.gene_symbol = self.db.get_gene_symbol_from_transcript_id(
                                my_variant.quibble.ac)
                        try:
                            toskip = mappers.transcripts_to_gene(
                                    my_variant,
                                    self,
                                    select_transcripts_dict_plus_version)
                        except mappers.MappersError:
                            my_variant.output_type_flag = 'warning'
                            continue
                        except vvhgvs.exceptions.HGVSInvalidVariantError:
                            my_variant.output_type_flag = 'warning'
                            continue
                        # For not already normalised input gene data fetching
                        # normally happens in structure_checks!
                        # (do we want to change the normal location too ?)
                        my_variant.gene_symbol = self.db.get_gene_symbol_from_transcript_id(
                                my_variant.quibble.ac)
                        if my_variant.gene_symbol == 'none':
                            my_variant.gene_symbol = ''
                        if toskip:
                            continue
                        my_variant.hgvs_transcript_variant = my_variant.quibble

                    # set output to variant type specific
                    if my_variant.reftype in [':n.',':t.',':c.'] and my_variant.hgvs_transcript_variant != '':
                        my_variant.output_type_flag = 'gene'
                    elif my_variant.reftype == ':g.':
                        my_variant.output_type_flag = 'intergenic'
                    elif my_variant.reftype == ':m.':
                        my_variant.output_type_flag = 'mitochondrial'

                    continue

                # Bug catcher
                try:
                    # Note, ID is not touched. It is always the input variant description.
                    # Quibble will be altered but id will not if type = g.
                    logger.info("Started validation of %s (originally %s)", str(my_variant.quibble),
                                my_variant.original)

                    # Find brackets and other at the beginning of the descriptions
                    if my_variant.non_alphanum_start():
                        my_variant.warnings.append(
                            "VariantSyntaxError: Variant descriptions must begin with a reference sequence "
                            "identifier or a chromosome number. Refer to "
                            "the examples provided at https://variantvalidator.org/service/validate/")
                        continue

                    if not my_variant.is_ascii():
                        chars, positions = my_variant.get_non_ascii()
                        error = 'VariantSyntaxError: Submitted ' \
                                'variant description contains an invalid character(s) %s at position(s) %s: ' \
                                'Please remove this character and re-submit: A useful search function for ' \
                                'Unicode characters can be found at https://unicode-search.net/' % (chars, positions)
                        my_variant.warnings.append(error)
                        logger.warning(error)
                        continue

                    # VCF line handling - Note: handling csv brings too many issues, so stick to tabs tsv
                    if "\t" in my_variant.quibble and not re.search(r"[gcrnmo]\.", my_variant.quibble):
                        try:
                            my_variant.quibble = vcf_to_pvcf.vcf_to_shorthand(my_variant.quibble)
                            my_variant.warnings.append(f"VcfConversionWarning: VCF line identified and converted "
                                                       f"to {my_variant.quibble}")
                        except vcf_to_pvcf.VcfConversionError as e:
                            logger.info(f"Cannot convert {my_variant.quibble} into PVCF format {e}")
                            continue
                    elif (re.search("\s+", my_variant.quibble) and not
                            re.search(r"[gcrnmo]\.", my_variant.quibble)):
                        try:
                            my_variant.quibble = vcf_to_pvcf.vcf_to_shorthand(my_variant.quibble)
                            my_variant.warnings.append(f"VcfConversionWarning: VCF line identified and converted "
                                                       f"to {my_variant.quibble}")
                        except vcf_to_pvcf.VcfConversionError as e:
                            logger.info(f"Cannot convert {my_variant.quibble} into PVCF format {e}")
                            continue

                    # Remove whitespace and quotes
                    my_variant.remove_whitespace()
                    my_variant.remove_quotes()
                    my_variant.remove_typos()

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

                    # Create the additional required normalizers which come from allele merge code and other sources
                    self.primary_assembly = primary_assembly
                    self.create_additional_normalizers_and_mappers()

                    logger.debug("Completed string formatting")

                    # Are submitted ENST transcripts coding or noncoding?
                    if "ENST" in my_variant.quibble or "NM_" in my_variant.quibble or "NR_" in my_variant.quibble:

                        match = re.search("(ENST|NM_|NR_)\d+\.\d+", my_variant.quibble)

                        if match:
                            result = match.group()

                            # Check if Ens submitted as RefSeq set and vice versa
                            if "ENST" in result and transcript_set == "refseq":
                                my_variant.warnings.append(
                                    "InvalidFieldError: The transcript " + result + " is not in the RefSeq "
                                    "data set. Please select Ensembl")
                                continue
                            elif ("NM_" in result or "NR_" in result) and transcript_set == "ensembl":
                                my_variant.warnings.append(
                                    "InvalidFieldError: The transcript " + result + " is not in the Ensembl "
                                    "data set. Please select RefSeq")
                                continue
                            try:
                                to_code_or_not_to_code = self.hdp.get_tx_identity_info(result)
                            except vvhgvs.exceptions.HGVSDataNotAvailableError as e:
                                if "No transcript definition for" in str(e):
                                    my_variant.warnings.append("The transcript " + result + " is not in "
                                                               "our database. Please check the transcript ID")
                                    versions_available = []
                                    for i in range(1, 20):
                                        accession, version = result.split(".")
                                        version = int(version) + i
                                        accession_version = accession + "." + str(version)
                                        try:
                                            to_code_or_not_to_code = self.hdp.get_tx_identity_info(accession_version)
                                        except vvhgvs.exceptions.HGVSDataNotAvailableError:
                                            continue
                                        else:
                                            versions_available.append(accession_version)
                                    if versions_available:
                                        my_variant.warnings.append("The following versions of the requested "
                                                                   "transcript are available in our database: "
                                                                   + "|".join(versions_available))
                                        continue
                                    else:
                                        continue

                            if to_code_or_not_to_code[3] is None:
                                my_variant.transcript_type = 'n'
                            else:
                                my_variant.transcript_type = 'c'

                    # Sort out o.
                    if ":o." in my_variant.quibble:
                        if "NC_012920.1" in my_variant.quibble or "NC_001807.4" in my_variant.quibble:
                            my_variant.quibble = my_variant.quibble.replace(":o.", ":m.")
                            my_variant.warnings.append("Reference sequence type o. should only be used for circular "
                                                       "reference sequences that are not mitochondrial. Instead use m.")
                        else:
                            my_variant.warnings.append("Reference sequence type o. should only be used for circular "
                                                       "reference sequences that are not mitochondrial. Instead use m.")
                    try:
                        toskip = format_converters.initial_format_conversions(my_variant, self,
                                                                              select_transcripts_dict_plus_version)

                    except vvhgvs.exceptions.HGVSError as e:
                        # import traceback
                        # traceback.print_exc()
                        checkref = str(e)
                        try:
                            # Test intronic variants for incorrect boundaries (see issue #169)
                            test_variant = copy.copy(my_variant)
                            test_variant.hgvs_formatted = my_variant.quibble
                            if type(test_variant.hgvs_formatted) is str:
                                test_variant.hgvs_formatted = self.hp.parse_hgvs_variant(
                                        test_variant.hgvs_formatted)

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
                            if re.search("ins\d+$", my_variant.quibble):
                                my_variant.warnings.append("The length of the variant is not formatted following the "
                                                           "HGVS guidelines. Please rewrite e.g. '10' to 'N[10]'"
                                                           "(where N is an unknown nucleotide)")
                                try:
                                    if "_" not in my_variant.quibble.split(":")[1] and \
                                            "del" not in my_variant.quibble.split(":")[1]:
                                        my_variant.warnings.append("An insertion must be provided with the two "
                                                                   "positions between which the insertion has taken "
                                                                   "place")
                                except IndexError:
                                    pass
                                continue
                            else:
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
                            if "The entered coordinates do not agree with the intron/exon boundaries for the selected "\
                               "transcript" not in my_variant.warnings[0]:
                                my_variant.warnings.reverse()
                            logger.warning(str(e))
                            continue

                    else:
                        if my_variant.warnings is not None and my_variant.hgvs_genomic is not None:
                            if "NC_" in str(my_variant.hgvs_genomic) and my_variant.reformat_output == "uncertain_pos":
                                my_variant.primary_assembly_loci = {my_variant.primary_assembly.lower():
                                                                    {"hgvs_genomic_description":my_variant.hgvs_genomic,
                                                                     "vcf": {"chr": None,
                                                                             "pos": None,
                                                                             "ref": None,
                                                                             "alt": None},}}
                    if type(my_variant.quibble) is str:
                        lovd_response = lovd_api.lovd_syntax_check(my_variant.original,
                                                                   do_lovd_check=self.lovd_syntax_check)
                        if "lovd_api_error" not in lovd_response.keys():
                            my_variant.output_type_flag = 'warning'
                            my_variant.lovd_syntax_check = lovd_response
                    if toskip:
                        continue

                    # INITIAL USER INPUT FORMATTING
                    initial_formatting.initial_user_formattng(my_variant, self)

                    # Set some configurations
                    formatted_variant = my_variant.quibble
                    stash_input = my_variant.quibble
                    my_variant.post_format_conversion = stash_input

                    logger.debug("Variant input formatted, proceeding to validate.")


                    # Change RNA bases to upper case but nothing else
                    if my_variant.reftype == ":r.":
                        query_r_var = copy.deepcopy(formatted_variant)
                        edit = formatted_variant.posedit.edit
                        if edit.ref:
                            edit.ref = edit.ref.lower()
                        if not edit.type in ['inv', 'dup'] and edit.alt:
                            edit.alt = edit.alt.lower()
                        formatted_variant.posedit.edit = edit
                        # do we need to limit the supported variant types?
                        # the case for the reftype needs to already have been checked at this point
                        if str(query_r_var) != str(formatted_variant):
                            e = "This not a valid HGVS description, due to characters being in the wrong case. " \
                                "Please check the use of upper- and lowercase characters."
                            my_variant.warnings.append(str(e))
                            logger.warning(str(e))


                    my_variant.hgvs_formatted = formatted_variant

                    if 'LRG' in my_variant.hgvs_formatted.ac:
                        my_variant.hgvs_formatted.ac.replace('T', 't')
                    else:
                        my_variant.hgvs_formatted.ac = my_variant.hgvs_formatted.ac.upper()

                    if my_variant.hgvs_formatted.type == "p" and my_variant.hgvs_formatted.posedit is None \
                            and ":p.?" in str(my_variant.hgvs_formatted):

                        # Protein variants needed early!
                        toskip = format_converters.proteins(my_variant, self)
                        if toskip:
                            continue

                    if hasattr(my_variant.hgvs_formatted.posedit.edit, 'alt'):
                        if my_variant.hgvs_formatted.posedit.edit.alt is not None:
                            my_variant.hgvs_formatted.posedit.edit.alt = \
                                my_variant.hgvs_formatted.posedit.edit.alt.upper()
                    if hasattr(my_variant.hgvs_formatted.posedit.edit, 'ref'):
                        if my_variant.hgvs_formatted.posedit.edit.ref is not None:
                            my_variant.hgvs_formatted.posedit.edit.ref = \
                                my_variant.hgvs_formatted.posedit.edit.ref.upper()
                    try:
                        formatted_variant = str(my_variant.hgvs_formatted)
                    except KeyError as e:
                        if "p" in my_variant.hgvs_formatted.type:
                            error = "Invalid amino acid %s stated in description %s" % (
                                    str(e),
                                    str(my_variant.quibble.format({'p_3_letter':False})))
                            my_variant.warnings.append(error)
                            continue

                    my_variant.set_quibble(my_variant.hgvs_formatted)
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
                        # user input e.g. alignment method and genome build
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

                    # skip if we have fuzzy ends
                    if isinstance(my_variant.quibble.posedit.pos.start, Interval) or \
                            isinstance(my_variant.quibble.posedit.pos.end, Interval):
                        continue

                    if '*' in str(my_variant.hgvs_formatted.posedit):
                        input_parses_copy = copy.deepcopy(my_variant.hgvs_formatted)
                        input_parses_copy.type = "c"
                        # Map to n. position
                        # Create easy variant mapper (over variant mapper) and splign locked evm
                        try:
                            to_n = my_variant.evm.c_to_n(input_parses_copy)
                        except Exception as e:
                            logger.debug("Error with to_n")
                            logger.debug(str(e))
                        except vvhgvs.exceptions.HGVSError as e:
                            logger.debug("Except passed, %s", e)
                        else:
                            if to_n.posedit.pos.end.base < to_n.posedit.pos.start.base:
                                error = 'Interval end position < interval start position '
                                my_variant.warnings.append(error)
                                logger.warning(error)
                                continue

                    elif my_variant.hgvs_formatted.posedit.pos.end.base < \
                            my_variant.hgvs_formatted.posedit.pos.start.base:
                        if my_variant.hgvs_formatted.ac not in ["NC_012920.1", "NC_001807.4"]:
                            error = 'Interval end position ' +\
                                    str(my_variant.hgvs_formatted.posedit.pos.end.base) + \
                                    ' < interval start position ' + \
                                    str(my_variant.hgvs_formatted.posedit.pos.start.base)
                            my_variant.warnings.append(error)
                            logger.warning(error)
                            continue

                    # Catch missing version number in refseq/ens
                    is_version = re.compile(r"\d\.\d")
                    if ((my_variant.refsource == 'RefSeq' or my_variant.refsource == 'ENS') and
                            not is_version.search(my_variant.hgvs_formatted.ac)):
                        error = 'RefSeq variant accession numbers MUST include a version number'
                        my_variant.warnings.append(error)
                        continue

                    logger.debug("HVGS interval/version mapping complete")

                    # handle LRG inputs
                    if my_variant.refsource == 'LRG':
                        format_converters.lrg_to_refseq(my_variant, self)
                        logger.debug("LRG check for conversion to refseq completed")

                    # Additional Incorrectly input variant capture training
                    if my_variant.refsource == 'RefSeq' or my_variant.refsource == 'ENS':
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

                    # Protein variants
                    toskip = format_converters.proteins(my_variant, self)
                    if toskip:
                        continue

                    # RNA variants
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
                        try:
                            toskip = mappers.gene_to_transcripts(my_variant, self, select_transcripts_dict)
                        except IndexError:
                            my_variant.output_type_flag = 'warning'
                            error = '%s cannot be validated in the context of genome build %s, ' \
                                    'try an alternative genome build' \
                                    % (str(my_variant.quibble), my_variant.primary_assembly)
                            my_variant.warnings.append(error)
                            toskip = True
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
                    logger.info("Completed initial validation for %s", str(my_variant.quibble))

                # Report errors to User and VV admin
                except KeyboardInterrupt:
                    raise
                except Exception:
                    lovd_response = lovd_api.lovd_syntax_check(my_variant.original.strip(),
                                                               do_lovd_check=self.lovd_syntax_check)
                    if "lovd_api_error" not in lovd_response.keys():
                        my_variant.output_type_flag = 'warning'
                        my_variant.lovd_syntax_check = lovd_response
                        continue
                    else:
                        my_variant.output_type_flag = 'warning'
                        error = (f"InvalidVariantError: Variant description {my_variant.original} could "
                                 f"not be validated by either "
                                 f"VariantValidator or the LOVD syntax checker. Please refer to the HGVS nomenclature "
                                 f"website at https://hgvs-nomenclature.org/stable/. For additional assistance "
                                 f"contact us at https://variantvalidator.org/help/contact/")
                        my_variant.warnings.append(error)
                        exc_type, exc_value, last_traceback = sys.exc_info()
                        logger.error(str(exc_type) + " " + str(exc_value))
                        continue

            # Outside the for loop
            ######################
            logger.debug("End of 1st for loop")
            # order the rows
            by_order = sorted(self.batch_list, key=lambda x: x.order)

            for variant in by_order:
                if type(variant.quibble) is str:
                    logger.debug(f"Formatting variant {variant.quibble}")
                else:
                    logger.debug("Formatting variant " + variant.quibble.format({'p_3_letter':False}))
                if not variant.write:
                    continue

                # Genomic sequence variation
                # Check for gapped delins<No_Alt>
                if variant.genomic_g and variant.genomic_g.posedit.edit.type == 'delins':
                    variant.genomic_g = hgvs_delins_parts_to_hgvs_obj(
                            variant.genomic_g.ac,
                            variant.genomic_g.type,
                            variant.genomic_g.posedit.pos,
                            '','')

                hgvs_genomic_variant = variant.genomic_g

                # genomic accession
                logger.debug("genomic accession")
                if hgvs_genomic_variant:
                    hgvs_genomic_variant = unset_hgvs_obj_ref(hgvs_genomic_variant)
                    genomic_accession = hgvs_genomic_variant.ac
                else:
                    genomic_accession = None
                # RefSeqGene variation
                logger.debug("RefSeqGene variation")
                refseqgene_variant = variant.genomic_r

                if not refseqgene_variant or type(refseqgene_variant) is str and 'RefSeqGene' in refseqgene_variant:
                    variant.warnings.append(refseqgene_variant)
                    lrg_variant = ''
                    refseqgene_variant = ''
                else:
                    rsg_ac = self.db.get_lrg_id_from_refseq_gene_id(str(refseqgene_variant.ac))
                    if rsg_ac[0] == 'none':
                        lrg_variant = ''
                    else:
                        hgvs_lrg = copy.deepcopy(refseqgene_variant)
                        hgvs_lrg.ac = rsg_ac[0]
                        lrg_variant = fn.valstr(hgvs_lrg)
                        if rsg_ac[1] != 'public':
                            variant.warnings.append('The current status of ' + str(hgvs_lrg.ac) + ' is pending '
                                                    'therefore changes may be made to the LRG reference sequence')

                # Transcript sequence variation
                logger.debug("Transcript sequence variation")
                hgvs_tx_variant = None
                if variant.coding:
                    if '(' in str(hgvs_tx_variant) and ')' in str(hgvs_tx_variant):
                        assert False

                    # transcript accession
                    logger.debug("transcript accession")
                    hgvs_tx_variant = unset_hgvs_obj_ref(variant.coding)
                    hgvs_transcript_variant = copy.deepcopy(hgvs_tx_variant)
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
                            hgvs_lrg_t = self.vm.g_to_t(refseqgene_variant, transcript_accession)
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
                    hgvs_transcript_variant = ''

                # Look for intronic variants
                logger.debug("Look for intronic variants")
                if transcript_accession != '' and genomic_accession:
                    # Remove del bases
                    hgvs_transcript_variant = unset_hgvs_obj_ref(hgvs_transcript_variant)
                    try:
                        self.vr.validate(hgvs_transcript_variant)
                    except vvhgvs.exceptions.HGVSError as e:
                        error = str(e)
                        if 'intronic variant' in error:
                            genome_context_transcript_variant = genomic_accession + '(' + transcript_accession +\
                                                                '):c.' + hgvs_transcript_variant.posedit.format({'max_ref_length': 0})
                            if refseqgene_variant:
                                refseqgene_variant = unset_hgvs_obj_ref(refseqgene_variant)
                                refseqgene_accession = refseqgene_variant.ac
                                try:
                                    hgvs_coding_from_refseqgene = self.vm.g_to_t(refseqgene_variant,
                                                                             hgvs_transcript_variant.ac)
                                except vvhgvs.exceptions.HGVSInvalidIntervalError:
                                    hgvs_coding_from_refseqgene = hgvs_transcript_variant
                                hgvs_coding_from_refseqgene = unset_hgvs_obj_ref(hgvs_coding_from_refseqgene)
                                refseqgene_context_transcript_variant = refseqgene_accession + '(' + \
                                    transcript_accession + '):c.' + str(hgvs_coding_from_refseqgene.posedit.pos) +\
                                        hgvs_coding_from_refseqgene.posedit.edit.format({'max_ref_length': 0})
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

                if type(predicted_protein_variant) is not str and 'NP_' in predicted_protein_variant.ac:
                    lrg_p = self.db.get_lrg_protein_id_from_ref_seq_protein_id(
                            predicted_protein_variant.ac)
                    if 'LRG' in lrg_p:
                        predicted_protein_variant.ac= predicted_protein_variant.ac + '(' + lrg_p + ')'

                # Gene
                if transcript_accession == '':
                    variant.gene_symbol = ''

                if hgvs_tx_variant:
                    multi_gen_vars = mappers.final_tx_to_multiple_genomic(variant,
                                                                          self,
                                                                          hgvs_tx_variant,#.format({'max_ref_length': 0}),
                                                                          liftover_level=liftover_level)

                else:
                    # HGVS genomic in the absence of a transcript variant
                    if hgvs_genomic_variant:
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
                    try:
                        vcf_dict = hgvs_utils.report_hgvs2vcf(alt_gen_var, 'All', variant.reverse_normalizer,
                                                              self.sf)
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        continue
                    # Identify primary assembly positions
                    primary = False
                    if 'NC_' in alt_gen_var.ac and par is False:
                        if 'NC_000' not in alt_gen_var.ac and 'NC_012920.1' not in alt_gen_var.ac and \
                                'NC_001807.4' not in alt_gen_var.ac:
                            continue
                        primary =True
                    elif 'NC_' not in alt_gen_var.ac and par is False:
                        pass
                    elif 'NC_000023' in alt_gen_var.ac and par is True:
                        primary =True
                    if primary:
                        for genome_build in vcf_dict['chrs_by_genome']:
                            primary_genomic_dicts[genome_build] = {
                                'hgvs_genomic_description': alt_gen_var,
                                'vcf': {'chr': vcf_dict['chrs_by_genome'][genome_build],
                                        'pos': vcf_dict['pos'],
                                        'ref': vcf_dict['ref'],
                                        'alt': vcf_dict['alt']
                                        }
                                }
                    else:
                        for genome_build in vcf_dict['chrs_by_genome']:
                            alt_dict = {genome_build: {
                                'hgvs_genomic_description': alt_gen_var,
                                'vcf': {'chr': vcf_dict['chrs_by_genome'][genome_build],
                                        'pos': vcf_dict['pos'],
                                        'ref': vcf_dict['ref'],
                                        'alt': vcf_dict['alt']
                                        }
                                }}
                            alt_genomic_dicts.append(alt_dict)

                # Clean up mito genome issues
                cp_lifted_response = copy.deepcopy(primary_genomic_dicts)
                for key, val in cp_lifted_response.items():
                    if key == "hg19" and "NC_012920.1" == val["hgvs_genomic_description"].ac:
                        primary_genomic_dicts.pop(key)
                    elif key == "grch37" and "NC_001807.4" == val["hgvs_genomic_description"].ac:
                        primary_genomic_dicts.pop(key)

                # Warn not directly mapped to specified genome build
                if genomic_accession:
                    if primary_assembly.lower() not in list(primary_genomic_dicts.keys()):
                        errors = [str(variant.hgvs_coding) + ' is not part of genome build ' + primary_assembly]

                        if self.alt_aln_method == "splign":
                            errors.append(str(variant.hgvs_coding) + ' cannot be mapped directly to genome build ' + primary_assembly)
                            errors.append('See alternative genomic loci or alternative genome builds for aligned genomic positions')

                        elif self.alt_aln_method == "genebuild":
                            # Get the alternative genome build to recommend
                            if primary_assembly == "GRCh38" or primary_assembly == "hg38":
                                alt_build = "GRCh37"
                            elif primary_assembly == "GRCh37" or primary_assembly == "hg19":
                                alt_build = "GRCh38"
                                # Shows the alternative genome build too
                            errors.append(str(variant.hgvs_coding) + ' cannot be mapped directly to genome build ' + primary_assembly
                                            + ', did you mean ' + alt_build + '?')

                        variant.warnings.extend(errors)


                # Ensure Variants have had the refs removed.
                # if not hasattr(posedit, refseqgene_variant):
                if refseqgene_variant:
                    try:
                        refseqgene_variant =  unset_hgvs_obj_ref(refseqgene_variant)
                    except Exception as e:
                        logger.debug("Except passed, %s", e)
                    if variant.gene_symbol == "" and refseqgene_variant:
                        gene_symbol = self.db.get_gene_symbol_from_refseq_id(refseqgene_variant.ac)
                        variant.gene_symbol = gene_symbol

                # Add predicted protein variant dictionary this is the output form so str for final is OK
                if predicted_protein_variant != '':
                    predicted_protein_variant_dict = {}
                    predicted_protein_variant_dict["slr"] = ''
                    predicted_protein_variant_dict["tlr"] = ''
                    predicted_protein_variant_dict["lrg_tlr"] = ''
                    predicted_protein_variant_dict["lrg_slr"] = ''
                    if type(predicted_protein_variant) is not str:
                        # add protein descriptions if not N type edit
                        add_p_descps = True
                        try:
                            if "N" in str(hgvs_tx_variant.posedit.edit):
                                add_p_descps = False
                        except AttributeError:
                            pass
                        if add_p_descps is True:
                            try:
                                # Remove LRG if present and store presence for later
                                if 'LRG' in predicted_protein_variant.ac:
                                    format_lrg = predicted_protein_variant.ac
                                    if "(" in format_lrg:
                                        format_lrg = format_lrg.split('(')[1]
                                        format_lrg = format_lrg.replace(')', '')
                                        predicted_protein_variant.ac  = re.sub(
                                                r'\(LRG_.+?\)', '', predicted_protein_variant.ac)
                                else:
                                    format_lrg = None
                                # Add single letter AA code to protein descriptions
                                predicted_protein_variant_dict = {"tlr": str(
                                    predicted_protein_variant.format({'max_ref_length': 0})
                                    ), "slr": ''}

                                if re.search("[A-Z][a-z][a-z]1[A-Z][a-z][a-z]", str(
                                    predicted_protein_variant.posedit)):
                                    cp_warnings = []
                                    for each_warning in variant.warnings:
                                        if "is HGVS compliant and contains a valid reference " \
                                           "amino acid description" not in each_warning:
                                            cp_warnings.append(each_warning)
                                        else:
                                            aa_1 = self.sf.fetch_seq(
                                                    predicted_protein_variant.ac,
                                                    start_i=0,
                                                    end_i=1)
                                            aa_1 = fn.one_to_three(aa_1)
                                            cp_format_p = f"{predicted_protein_variant.ac}:p.({aa_1}1?)"
                                            cp_warnings.append(
                                                    f"Variant {predicted_protein_variant} affects the initiation amino acid"
                                                    f" so is better described as {cp_format_p}")
                                            predicted_protein_variant = vvhgvs.sequencevariant.SequenceVariant(
                                                    ac=predicted_protein_variant.ac,
                                                    type='p',
                                                    posedit="({aa_1}1?)")
                                            variant.warnings = cp_warnings

                                # Set formatted tlr
                                predicted_protein_variant_dict['tlr'] = \
                                        predicted_protein_variant.format({
                                            'max_ref_length': 0})
                                predicted_protein_variant_dict['slr']= \
                                        predicted_protein_variant.format({
                                            'max_ref_length': 0,
                                            'p_3_letter':False})
                                # set LRG outputs
                                if format_lrg is not None:
                                    predicted_protein_variant_dict["lrg_tlr"] = \
                                        format_lrg + ':' + \
                                        predicted_protein_variant_dict["tlr"].split(':')[1]
                                    predicted_protein_variant_dict["lrg_slr"] = \
                                        format_lrg + ':' + \
                                        predicted_protein_variant_dict["slr"].split(':')[1]
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

                # Add missing gene info which should be there (May have come from uncertain positions for example)
                if variant.hgvs_transcript_variant and variant.gene_symbol == '':
                    # the last hold out for this is uncertain positions
                    variant.gene_symbol = self.db.get_gene_symbol_from_transcript_id(
                        variant.hgvs_transcript_variant.ac)
                elif variant.hgvs_refseqgene_variant and variant.gene_symbol == '':
                    variant.gene_symbol = self.db.get_gene_symbol_from_refseq_id(
                        variant.hgvs_refseqgene_variant.ac)

                # Add stable gene_ids
                stable_gene_ids = {}
                if variant.gene_symbol != '':
                    gene_stable_info = self.db.get_stable_gene_id_info(variant.gene_symbol)

                    # Add or update stable ID and transcript data
                    if gene_stable_info[1] == 'No data' and hgvs_tx_variant is not None:
                        try:
                            self.db.update_transcript_info_record(hgvs_tx_variant.ac, self,
                                                                  genome_build=self.selected_assembly)
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
                            self.db.update_transcript_info_record(hgvs_tx_variant.ac, self,
                                                                  genome_build=self.selected_assembly)
                        except fn.DatabaseConnectionError as e:
                            error = 'Currently unable to update all gene_ids or transcript information records ' \
                                    'because ' \
                                    'VariantValidator %s' % str(e)
                            my_variant.warnings.append(error)
                            logger.warning(error)
                    i = 1
                    while i in range(10):
                        annotation_info = self.db.get_transcript_annotation(hgvs_tx_variant.ac)
                        try:
                            json.loads(annotation_info)
                        except json.decoder.JSONDecodeError:
                            i += 1
                            time.sleep(2)
                        else:
                            break
                    reference_annotations = json.loads(annotation_info)

                variant.stable_gene_ids = stable_gene_ids
                variant.annotations = reference_annotations
                variant.genome_context_intronic_sequence = genome_context_transcript_variant
                variant.refseqgene_context_intronic_sequence = refseqgene_context_transcript_variant
                variant.hgvs_refseqgene_variant = refseqgene_variant
                variant.hgvs_predicted_protein_consequence = predicted_protein_variant_dict
                variant.hgvs_lrg_transcript_variant = lrg_transcript_variant
                variant.hgvs_lrg_variant = lrg_variant
                variant.alt_genomic_loci = alt_genomic_dicts
                if variant.reformat_output != "uncertain_pos":
                    variant.primary_assembly_loci = primary_genomic_dicts
                    if hgvs_tx_variant:
                        variant.hgvs_transcript_variant = hgvs_tx_variant
                else:
                     for mapping in variant.alt_genomic_loci:
                        for gennome in mapping:
                            mapping[gennome]["vcf"] = {
                                    'chr': None,
                                    'pos': None,
                                    'ref': None,
                                    'alt': None}
                if not variant.hgvs_transcript_variant and hgvs_tx_variant:
                    variant.hgvs_transcript_variant = hgvs_tx_variant
                variant.reference_sequence_records = ''
                variant.validated = True

                # Add links to reference_sequence_records
                pre_out = {
                        'hgvs_transcript_variant':'',
                        'hgvs_predicted_protein_consequence':{'slr':''},
                        'hgvs_refseqgene_variant':'',
                        'hgvs_lrg_variant':'',
                        'selected_assembly':self.selected_assembly}
                if variant.hgvs_transcript_variant:
                    pre_out['hgvs_transcript_variant'] = variant.hgvs_transcript_variant.ac
                if variant.hgvs_refseqgene_variant:
                    pre_out['hgvs_refseqgene_variant'] = variant.hgvs_refseqgene_variant.ac
                if variant.hgvs_lrg_variant:# is str
                    pre_out['hgvs_lrg_variant'] = variant.hgvs_lrg_variant
                if variant.hgvs_predicted_protein_consequence:
                    pre_out['hgvs_predicted_protein_consequence']['slr'] = \
                            variant.hgvs_predicted_protein_consequence['slr']
                ref_records = self.db.get_urls(pre_out)
                if ref_records != {}:
                    variant.reference_sequence_records = ref_records

                # Loop out uncertain position variants
                if variant.reformat_output != "uncertain_pos":

                    # Liftover intergenic positions genome to genome
                    if (variant.output_type_flag == 'intergenic' and liftover_level is not None) or \
                            (('grch37' not in variant.primary_assembly_loci.keys() or
                              'grch38' not in variant.primary_assembly_loci.keys() or
                              'hg38' not in variant.primary_assembly_loci.keys() or
                              'hg19' not in variant.primary_assembly_loci.keys())
                             and liftover_level is not None):

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
                            if variant.output_type_flag != 'intergenic' and variant.output_type_flag != "other":
                                g_to_g = True

                            # Lift-over
                            if (str(genomic_position_info[g_p_key]['hgvs_genomic_description']) not in lo_cache.keys()
                                ) or ("NC_012920.1" == genomic_position_info[g_p_key]['hgvs_genomic_description'].ac
                                    and build_from == "hg38" and build_to == "hg19"):

                                lifted_response = liftover(genomic_position_info[g_p_key]['hgvs_genomic_description'],
                                                           build_from,
                                                           build_to, variant.hn, variant.reverse_normalizer,
                                                           variant.evm,
                                                           self,
                                                           specify_tx=False,
                                                           liftover_level=liftover_level,
                                                           g_to_g=g_to_g,
                                                           genomic_data_w_vcf=genomic_position_info)
                                if "NC_012920.1" == genomic_position_info[g_p_key]['hgvs_genomic_description'].ac or \
                                        "NC_001807.4" == genomic_position_info[g_p_key]['hgvs_genomic_description'].ac:
                                    capture_corrected_response = False
                                    for key, val in lifted_response.items():
                                        if "grch38" in key:
                                            capture_corrected_response = val
                                        if val == {} and "grch37" in key:
                                            if capture_corrected_response is not False:
                                                lifted_response[key] = capture_corrected_response
                                    for key, val in lifted_response.items():
                                        if "grch38" in key:
                                            capture_corrected_response = val
                                        if val == {} and "grch37" in key:
                                            if capture_corrected_response is not False:
                                                lifted_response[key] = capture_corrected_response

                                lo_cache[str(genomic_position_info[g_p_key]['hgvs_genomic_description'])] \
                                        = lifted_response
                            else:
                                lifted_response = \
                                        lo_cache[str(genomic_position_info[g_p_key]['hgvs_genomic_description'])]

                            # Sort the respomse into primary assembly and ALT
                            primary_assembly_loci = {}
                            alt_genomic_loci = []

                            for build_key, accession_dict in list(lifted_response.items()):
                                try:
                                    accession_key = list(accession_dict.keys())[0]
                                    if accession_dict[accession_key]['hgvs_genomic_description'].ac.startswith('NC_'):
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

                # Add exon numbering information see issue #
                if variant.coding != "":
                    try:
                        exs = exon_numbering.finds_exon_number(variant, self)
                        variant.exonic_positions = exs
                    except KeyError:
                        pass

                """
                Update / Remove warnings and add error codes

                The following if statements can be used to amend warnings and add missing error codes using String find
                or RegEx searches in existing warnings
                """

                # Remove duplicate warnings
                variant_warnings = []
                accession = ''
                if variant.hgvs_transcript_variant:
                    accession = variant.hgvs_transcript_variant.ac
                term = str(accession)
                term_2 = "%s automapped to" % str(hgvs_tx_variant)
                term_3 = "%s automapped to" % str(hgvs_genomic_variant)

                for vt in variant.warnings:
                    vt = str(vt)

                    # Base mismatch errors
                    if "does not agree with reference sequence" in vt:
                        if "mapped to" in vt:
                            variant_warnings.append(f"IntronBaseMismatchWarning: {vt}")
                        else:
                            variant_warnings.append(f"ReferenceMismatchError: {vt}")
                        continue

                    # Update out of bounds error
                    if vt == "start or end or both are beyond the bounds of transcript record":
                        variant_warnings.append(f"OutOfBoundsError: {vt}")
                        continue

                    # Do not warn transcript not part of build if it's not the relevant transcript
                    if "is not part of genome build" in vt and term not in vt:
                        continue
                    # Do not warn transcript cannot be mapped to build if it's not the relevant transcript
                    elif "cannot be mapped directly to genome build" in vt and term not in vt:
                        continue
                    #  Do not warn a transcript update is available for the most recent transcript
                    elif term in vt and "A more recent version of the selected reference sequence" in vt:
                        vt = vt.split(": ")
                        vt = ": ".join([vt[0], vt[1]])
                        variant_warnings.append(vt)

                    elif "expected a letter" in vt:
                        variant_warnings.append("InvalidVariantError: Accepted formats are HGVS, pseudoVCF. "
                                                "Refer to the examples provided at https://variantvalidator."
                                                "org/service/validate/ for more information.")

                    # Remove spurious updates away form the correct output
                    elif (term_2 in vt and hgvs_tx_variant) or (term_3 in vt and hgvs_genomic_variant):
                        continue
                    # Suppress "RefSeqGene record not available"
                    elif "RefSeqGene record not available" in vt:
                        continue
                    # convert SeqRepo errors into a more user friendly form
                    elif vt.startswith('Failed to fetch') and 'SeqRepo' in vt:
                        acc = vt.split()[3]
                        vt = (
                                f"Failed to find {acc} in our sequence store: This may mean that "+
                                "the sequence has been mistyped, or it may be missing from our "+
                                "data, possibly due to being either deprecated or yet to be added.")
                        variant_warnings.append(vt)
                    elif 'NP_' in vt and 'transcript' in vt:
                        continue
                    elif "Xaa" in vt:
                        vt = vt.replace("Xaa", "Ter")
                        variant_warnings.append(vt)
                    elif "automapped to" in vt:
                        vt = re.sub(r"(del|dup)[A-Z]+", r"\1", vt)
                        variant_warnings.append(vt)
                        continue
                    else:
                        variant_warnings.append(vt)
                variant.warnings = variant_warnings

                # Reformat as required to add back variation that would/does get lost on mapping
                if variant.reformat_output is not False:
                    if "|" in variant.reformat_output and "=" in str(variant.quibble):
                        def _apply_met_variation(data):
                            if isinstance(data, dict):
                                for key in data:
                                    if isinstance(data[key],dict) or isinstance(data[key],list):
                                        data[key] = _apply_met_variation(data[key])
                                    elif isinstance(data[key],SequenceVariant) and not data[key].type == 'p':
                                        if type(data[key].posedit) is PosEdit:
                                            data[key] = to_vv_hgvs(data[key])
                                        data[key].posedit.met_variation = variant.reformat_output
                                    elif isinstance(data[key], str) and data[key].endswith('=') and not data[key].endswith('|met=') and not ':p.' in data[key]:
                                        data[key] = data[key][:-1] + variant.reformat_output
                            elif isinstance(data,list):
                                for index, value in enumerate(data):
                                    if isinstance(value,dict) or isinstance(value,list):
                                        data[index] = _apply_met_variation(value)
                                    elif isinstance(value,SequenceVariant) and not value.type == 'p':
                                        if type(value.posedit) is PosEdit:
                                            value = to_vv_hgvs(value)
                                        value.posedit.met_variation = variant.reformat_output
                                        data[index] = value
                                    elif isinstance(value, str) and value.endswith('=') and not value.endswith('|met=') and not ':p.' in value:
                                        data[index] = value[:-1] + variant.reformat_output
                            elif isinstance(data,SequenceVariant) and not data.type == 'p':
                                if type(data.posedit) is PosEdit:
                                    data = to_vv_hgvs(data)
                                data.posedit.met_variation = variant.reformat_output
                            elif isinstance(data, str) and data.endswith('=') and not data.endswith('|met=') and not ':p.' in data:
                                data = data[:-1] + variant.reformat_output
                            return data

                        attributes = dir(variant)
                        for attribute in attributes:
                            if "__" in attribute or attribute == "reformat_output":
                                continue
                            item = (getattr(variant, attribute))
                            item = _apply_met_variation(item)
                            setattr(variant, attribute, item)

                # Add expanded repeat information
                logger.info(f"expanded repeat is {variant.expanded_repeat}")
                if variant.expanded_repeat is not None:
                    starting_tx_posedit = None
                    if variant.hgvs_transcript_variant:
                        starting_tx_posedit = variant.hgvs_transcript_variant.posedit.format(
                                {'max_ref_length': 0})
                    hgd = "hgvs_genomic_description"
                    ex_rep_start = variant.expanded_repeat["variant"].ac[:3]
                    ex_rep_var = copy.copy(variant.expanded_repeat["variant"])
                    # If the ER came from a gene description, set the HGVS transcript variant
                    if ex_rep_start in ["NG_", "LRG"]:
                        try:
                            variant.hgvs_transcript_variant = convert_seq_state_to_expanded_repeat(
                                    variant.hgvs_transcript_variant, self,
                                    known_repeat_unit=variant.expanded_repeat["repeat_sequence"],
                                    genomic_reference=ex_rep_var.ac)
                        except Exception:
                            pass
                    elif ex_rep_start == 'NC_':
                        variant.primary_assembly_loci[self.primary_assembly.lower()
                            ][hgd] = ex_rep_var
                        if self.primary_assembly == "GRCh37":
                            variant.primary_assembly_loci['grch37'][hgd] = ex_rep_var
                            variant.primary_assembly_loci["hg19"][hgd] = ex_rep_var
                            try:
                                variant.primary_assembly_loci["grch38"][hgd] = \
                                    convert_seq_state_to_expanded_repeat(
                                        variant.primary_assembly_loci["grch38"][hgd], self,
                                        known_repeat_unit=variant.expanded_repeat["repeat_sequence"])
                            except Exception:
                                pass
                            variant.primary_assembly_loci["hg38"][hgd] = \
                                    variant.primary_assembly_loci["grch38"][hgd]
                        else:
                            variant.primary_assembly_loci["grch38"][hgd] = ex_rep_var
                            variant.primary_assembly_loci["hg38"][hgd] = ex_rep_var
                            try:
                                variant.primary_assembly_loci["grch37"][hgd] = \
                                    convert_seq_state_to_expanded_repeat(
                                        variant.primary_assembly_loci["grch37"][hgd],
                                        self, known_repeat_unit=variant.expanded_repeat["repeat_sequence"])
                                variant.primary_assembly_loci["hg19"][hgd] = \
                                    variant.primary_assembly_loci["grch37"][hgd]
                            except Exception:
                                pass
                        try:
                            variant.hgvs_transcript_variant = \
                                convert_seq_state_to_expanded_repeat(
                                    variant.hgvs_transcript_variant, self,
                                    genomic_reference=ex_rep_var.ac,
                                known_repeat_unit=variant.expanded_repeat["repeat_sequence"])
                        except Exception:
                            pass

                    elif ex_rep_start in ["NM_", "NR_"] or ex_rep_var.ac.startswith("ENST"):
                        variant.hgvs_transcript_variant = copy.copy(ex_rep_var)
                        try:
                            prim_a = self.primary_assembly.lower()
                            variant.primary_assembly_loci[prim_a][hgd] = \
                                convert_seq_state_to_expanded_repeat(
                                    variant.primary_assembly_loci[prim_a][hgd], self,
                                    known_repeat_unit=variant.expanded_repeat["repeat_sequence"])
                            if self.primary_assembly == "GRCh37":
                                variant.primary_assembly_loci["hg19"][hgd] = \
                                    variant.primary_assembly_loci['grch37'][hgd]
                                variant.primary_assembly_loci["grch38"][hgd] = \
                                    convert_seq_state_to_expanded_repeat(
                                        variant.primary_assembly_loci["grch38"][hgd], self,
                                        known_repeat_unit=variant.expanded_repeat["repeat_sequence"])
                                variant.primary_assembly_loci["hg38"][hgd] = \
                                    variant.primary_assembly_loci["grch38"][hgd]
                            else:
                                variant.primary_assembly_loci["hg38"][hgd] = \
                                    variant.primary_assembly_loci[prim_a][hgd]
                                variant.primary_assembly_loci["grch37"][hgd] = \
                                    convert_seq_state_to_expanded_repeat(
                                        variant.primary_assembly_loci["grch37"][hgd], self,
                                        known_repeat_unit=variant.expanded_repeat["repeat_sequence"])
                                variant.primary_assembly_loci["hg19"][hgd] = \
                                    variant.primary_assembly_loci["grch37"][hgd]
                        except Exception:
                            pass

                    if variant.hgvs_refseqgene_variant and not ex_rep_start in ["NG_", "LRG"]:
                        assert variant.hgvs_refseqgene_variant.ac.startswith('NG_')
                        try:
                            variant.hgvs_refseqgene_variant = (
                            convert_seq_state_to_expanded_repeat(
                                variant.hgvs_refseqgene_variant, self,
                            known_repeat_unit=variant.expanded_repeat["repeat_sequence"]))
                            if variant.hgvs_lrg_variant:
                                variant.hgvs_lrg_variant = (
                                    f"{variant.hgvs_lrg_variant.split(':g.')[0]}:g." +
                                    variant.hgvs_refseqgene_variant.posedit.format(
                                        {'max_ref_length': 0}))
                        except Exception:
                            pass
                        # LRG transcript data can be, but may not be, the same as the main mapped tx
                        # so we need to test for identity before using. LRG is already text and
                        # won't be reused for VRS. This skips ex_rep conversion in the corner case.
                        if variant.hgvs_lrg_transcript_variant and \
                                variant.hgvs_lrg_transcript_variant.split(':c.')[1] ==\
                                starting_tx_posedit and variant.hgvs_transcript_variant:
                            variant.hgvs_lrg_transcript_variant = \
                                f"{variant.hgvs_lrg_transcript_variant.split(':c.')[0]}:c."\
                                + variant.hgvs_transcript_variant.posedit.format(
                                    {'max_ref_length': 0})

                # Some variant objects must be strings for back-compatibility with old output
                if variant.hgvs_transcript_variant:
                    variant.hgvs_transcript_variant = \
                            variant.hgvs_transcript_variant.format({'max_ref_length': 0})
                else:
                    variant.hgvs_transcript_variant = ''
                if variant.hgvs_refseqgene_variant:
                    variant.hgvs_refseqgene_variant = \
                            variant.hgvs_refseqgene_variant.format({'max_ref_length': 0})
                else:
                    variant.hgvs_refseqgene_variant = ''
                hgd = "hgvs_genomic_description"
                for gen in  variant.primary_assembly_loci.keys():
                    variant.primary_assembly_loci[gen][hgd] = \
                        variant.primary_assembly_loci[gen][hgd].format({'max_ref_length': 0})
                for loc in variant.alt_genomic_loci:
                    for gen in loc.keys():
                        loc[gen][hgd] = loc[gen][hgd].format({'max_ref_length': 0})

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

    def gene2transcripts(self, query, validator=False, bypass_web_searches=False, select_transcripts=None,
                         transcript_set="refseq", genome_build=None, batch_output=False, bypass_genomic_spans=False,
                         lovd_syntax_check=False):

        try:
            gene_symbols = json.loads(query)
            if isinstance(gene_symbols, int):
                batch_output = False
            elif isinstance(gene_symbols, str):
                batch_output = False
            elif isinstance(gene_symbols, list):
                batch_output = True
            else:
                batch_output = False
        except json.decoder.JSONDecodeError:
            if isinstance(query, list):
                gene_symbols = query
                batch_output = True
        except TypeError:
            if isinstance(query, list):
                gene_symbols = query
                batch_output = True
            pass

        if batch_output is False:
            g2d_data = gene2transcripts.gene2transcripts(self, query, validator, bypass_web_searches,
                                                     select_transcripts, transcript_set, genome_build,
                                                     bypass_genomic_spans, lovd_syntax_check)
        else:
            g2d_data = []
            for symbol in gene_symbols:
                data_for_gene = gene2transcripts.gene2transcripts(self, symbol, validator, bypass_web_searches,
                                                              select_transcripts, transcript_set, genome_build,
                                                              bypass_genomic_spans, lovd_syntax_check)
                g2d_data.append(data_for_gene)

        # return
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
        hgvs_vt = variant.hgvs_formatted
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
            hgvs_object = variant.hgvs_formatted
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
                    else:
                        return True
                except Exception as e:
                    error = 'Unable to assign transcript identity records to %s.  Please try again later ' \
                            'and if the problem persists contact admin. error=%s.' % (accession, str(e))
                    variant.warnings.append(error)
                    logger.info(error)
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
            hgvs_object = variant.hgvs_formatted
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
                    try:
                        self.db.update_transcript_info_record(accession, validator=self,
                                                              genome_build=variant.selected_assembly)
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
                        elif "Cannot retrieve data from Ensembl REST for record" in str(e):
                            error = ('%s. Please try an alternate genome build, or try again later and if the problem '
                                     'persists contact admin.') % str(e)
                            variant.warnings.append(error)
                            logger.warning(error)
                    entry = self.db.in_entries(accession, 'transcript_info')
                    variant.description = entry['description']
                else:
                    variant.description = entry['description']

            # If the none key is found add the description to the database
            elif 'none' in entry:
                try:
                    self.db.update_transcript_info_record(accession, validator=self,
                                                          genome_build=variant.selected_assembly)
                except Exception as e:
                    logger.info(str(e))
                    error = 'Unable to assign transcript identity records to ' + accession + \
                            ', potentially an obsolete record or there is an issue retrieving data from Ensembl. ' \
                            'Please try again later and if the problem' \
                            'persists contact admin'
                    variant.warnings.append(error)
                    logger.info(error)
                    return True
                entry = self.db.in_entries(accession, 'transcript_info')
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
