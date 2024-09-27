from cmath import log
import vvhgvs
import vvhgvs.exceptions
import vvhgvs.normalizer
import re
import copy
import sys
import logging
import json
import time
from vvhgvs.assemblymapper import AssemblyMapper
from VariantValidator.modules import hgvs_utils
from VariantValidator.modules import utils as fn
from VariantValidator.modules import seq_data
from VariantValidator.modules import vvMixinConverters
from VariantValidator.modules.variant import Variant
from VariantValidator.modules import format_converters
from VariantValidator.modules import use_checking
from VariantValidator.modules import mappers
from VariantValidator.modules import valoutput
from VariantValidator.modules import exon_numbering
from VariantValidator.modules.liftover import liftover
from VariantValidator.modules import complex_descriptions
from VariantValidator.modules import gene2transcripts
from VariantValidator.modules import expanded_repeats, methyl_syntax

logger = logging.getLogger(__name__)

class Mixin(vvMixinConverters.Mixin):
    """
    This module contains the main function for variant validator.
    It's added to the Validator object in the vvObjects file.
    """

    def validate(self,
                 batch_variant,
                 selected_assembly,
                 select_transcripts,
                 transcript_set=None,
                 liftover_level=False):
        """
        This is the main validator function.
        :param batch_variant: A string containing the variant to be validated
        :param selected_assembly: The version of the genome assembly to use.
        :param select_transcripts: Can be an array of different transcripts, or 'all'
        :param liftover_level: True or False - liftover to different gene/genome builds or not
        Selecting multiple transcripts will lead to a multiple variant outputs.
        :param transcript_set: 'refseq' or 'ensembl'
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

            # Turn each variant into a dictionary. The dictionary will be compiled during validation
            self.batch_list = []
            for queries in batch_queries:
                if isinstance(batch_queries, int):
                    batch_queries = str(batch_queries)
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
                my_variant.cross_hn = vvhgvs.normalizer.Normalizer(self.hdp,
                                                                   cross_boundaries=True,
                                                                   shuffle_direction=3,
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
                            if "NC_" in my_variant.hgvs_genomic and my_variant.reformat_output == "uncertain_pos":
                                my_variant.primary_assembly_loci = {my_variant.primary_assembly.lower():
                                                                    {"hgvs_genomic_description":
                                                                     my_variant.hgvs_genomic,
                                                                     "vcf": {"chr": None,
                                                                             "pos": None,
                                                                             "ref": None,
                                                                             "alt": None},}}
                    if toskip:
                        continue

                    # INITIAL USER INPUT FORMATTING
                    """
                    In this section of the code we are compiling HGVS errors and providing improved warnings/error 
                    messages
                    """
                    # 1. Requested warnings from https://github.com/openvar/variantValidator/issues/195
                    if re.search(r'\(.+?\)', my_variant.quibble):  # Pattern looks for (....)
                        gene_symbol_query = re.search(r'\(.+?\)', my_variant.quibble).group(0)
                        gene_symbol_query = gene_symbol_query.replace('(', '')
                        gene_symbol_query = gene_symbol_query.replace(')', '')
                        is_it_a_gene = self.db.get_hgnc_symbol(gene_symbol_query)
                        if is_it_a_gene != 'none':
                            warning = "Removing redundant gene symbol %s from variant description" % is_it_a_gene
                            my_variant.warnings.append(warning)
                            logger.warning(warning)
                    if re.search('del[GATC]+', my_variant.original) or re.search('inv[GATC]+', my_variant.original) \
                            or \
                       re.search('dup[GATC]+', my_variant.original) or re.search('ins[GATC]+', my_variant.original):

                        if not re.search('ins[GATC]+', my_variant.original):
                            warning = "Removing redundant reference bases from variant description"
                            my_variant.warnings.append(warning)
                            logger.warning(warning)

                    invalid = my_variant.format_quibble()

                    # 2. expand options for issue https://github.com/openvar/variantValidator/issues/338
                    test_for_invalid_case_in_accession = my_variant.original.split(":")[0]

                    # Basically, all reference sequences must be upper case, so we make an upper-case query accession
                    # to test the input accession against and try to spot a discrepancy
                    # The exception to the rule is LTG transcripts e.g. LRG_1t1 which we handle immediately below!
                    query_for_invalid_case_in_accession = test_for_invalid_case_in_accession.upper()
                    if re.match("LRG", test_for_invalid_case_in_accession, flags=re.IGNORECASE):
                        if "LRG" not in test_for_invalid_case_in_accession:
                            e = "This not a valid HGVS description, due to characters being in the wrong case. " \
                                "Please check the use of upper- and lowercase characters."
                            my_variant.warnings.append(str(e))
                            logger.warning(str(e))
                        if "T" in test_for_invalid_case_in_accession:
                            e = "This not a valid HGVS description, due to characters being in the wrong case. " \
                                "Please check the use of upper- and lowercase characters."
                            my_variant.warnings.append(str(e))
                            logger.warning(str(e))
                            my_variant.quibble = my_variant.quibble.replace("T", "t")

                    # Reference sequence types other than LRG
                    elif (test_for_invalid_case_in_accession != query_for_invalid_case_in_accession) \
                            and "LRG" not in test_for_invalid_case_in_accession:
                        # See issue #357
                        if re.match("chr", test_for_invalid_case_in_accession, re.IGNORECASE
                                    ) or re.match("GRCh", test_for_invalid_case_in_accession, re.IGNORECASE
                                                  ) or re.match("hg", test_for_invalid_case_in_accession, re.IGNORECASE
                                                                ):
                            e = "This is not a valid HGVS variant description, because no reference sequence ID " \
                                "has been provided"
                        else:
                            e = "This not a valid HGVS description, due to characters being in the wrong case. " \
                                "Please check the use of upper- and lowercase characters."
                        my_variant.warnings.append(str(e))
                        logger.warning(str(e))

                    if invalid:
                        if re.search(r'\w+:[gcnmrp],', my_variant.quibble):
                            error = 'Variant description ' + my_variant.quibble + ' contained the , character between '\
                                    '<type> and <position> in the expected pattern <accession>:<type>.<position> and ' \
                                    'has been auto-corrected'
                            my_variant.quibble = my_variant.quibble.replace(',', '.')
                            my_variant.warnings.append(error)
                            logger.warning(error)
                            pass

                        # Upper case type see issue #338
                        elif re.search(r":[GCNMR].", str(my_variant.quibble)):
                            rs_type_upper = re.search(r":[GCNMR].", str(my_variant.quibble))
                            e = "This not a valid HGVS description, due to characters being in the wrong case. " \
                                "Please check the use of upper- and lowercase characters."
                            my_variant.warnings.append(str(e))
                            logger.warning(str(e))
                            my_variant.quibble = my_variant.quibble.replace(rs_type_upper.group(0),
                                                                            rs_type_upper.group(0).lower())

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

                    else:
                        # 3. Here we handle syntax errors in ins and delins variants
                        # https://github.com/openvar/variantValidator/issues/359
                        if re.search("ins$", my_variant.quibble):
                            my_variant.warnings.append("The inserted sequence must be provided for insertions or "
                                                       "deletion-insertions")
                            try:
                                if "_" not in my_variant.quibble.split(":")[1] and \
                                        "del" not in my_variant.quibble.split(":")[1]:
                                    my_variant.warnings.append("An insertion must be provided with the two positions "
                                                               "between which the insertion has taken place")
                            except IndexError:
                                pass
                            continue

                        elif re.search("ins\(\d+\)$", my_variant.quibble):
                            my_variant.warnings.append("The length of the variant is not formatted following the "
                                                       "HGVS guidelines. Please rewrite e.g. '(10)' to 'N[10]'"
                                                       "(where N is an unknown nucleotide)")
                            try:
                                if "_" not in my_variant.quibble.split(":")[1] and \
                                        "del" not in my_variant.quibble.split(":")[1]:
                                    my_variant.warnings.append("An insertion must be provided with the two positions "
                                                               "between which the insertion has taken place")
                            except IndexError:
                                pass
                            continue

                        elif re.search("ins\d+$", my_variant.quibble):
                            my_variant.warnings.append("The length of the variant is not formatted following the HGVS "
                                                       "guidelines. Please rewrite e.g. '10' to 'N[10]'"
                                                       "(where N is an unknown nucleotide)")
                            try:
                                if "_" not in my_variant.quibble.split(":")[1] and \
                                        "del" not in my_variant.quibble.split(":")[1]:
                                    my_variant.warnings.append("An insertion must be provided with the two positions "
                                                               "between which the insertion has taken place")
                            except IndexError:
                                pass
                            continue

                        elif re.search("ins\(\d+_\d+\)$", my_variant.quibble):
                            my_variant.warnings.append("The length of the variant is not formatted following the HGVS "
                                                       "guidelines. Please rewrite e.g. '(10_20)' to 'N[(10_20)]'"
                                                       "(where N is an unknown nucleotide and [(10_20)] is an uncertain"
                                                       " number of N nucleotides ranging from 10 to 20)")
                            continue

                        elif re.search("ins\[\(\d+_\d+\)\]$", my_variant.quibble):
                            counts = re.findall("\d+", my_variant.quibble.split("ins")[1])

                            if int(counts[1]) < int(counts[0]):
                                wrn = "The length of the variant is not formatted following the HGVS guidelines. " \
                                      "Please rewrite (%s_%s) to N[(%s_%s)]" % (counts[0], counts[1],
                                                                                counts[1], counts[0])
                                my_variant.warnings.append(wrn)
                            elif int(counts[1]) == int(counts[0]):
                                wrn = "The length of the variant is not formatted following the HGVS guidelines. " \
                                      "Please rewrite (%s_%s) to N[(%s)]" % (counts[0], counts[1], counts[1])
                                my_variant.warnings.append(wrn)

                            try:
                                if not re.search("\d_\d", my_variant.quibble.split("ins")[0]) and \
                                        "del" not in my_variant.quibble.split(":")[1]:
                                    my_variant.warnings.append("An insertion must be provided with the two positions "
                                                               "between which the insertion has taken place")
                            except IndexError:
                                pass

                            if my_variant.warnings == []:
                                wrn = "The variant description is syntactically correct " \
                                      "but no further validation is possible because the description contains " \
                                      "uncertainty"
                                my_variant.warnings.append(wrn)

                            continue

                        elif re.search("(?:delins|del|ins)[NGATC]\[\d+\]$", my_variant.quibble) or \
                                re.search("(?:delins|del|ins)\[[NGATC]\[\d+\];", my_variant.quibble):

                            match = re.search("(?:delins|del|ins)", my_variant.quibble)[0]

                            if re.search(f"{match}\[[GATCN]+\[\d+\];", my_variant.quibble):
                                sections = my_variant.quibble.split(match)[1]
                                sections = sections[1:-1]
                                sections_listed = sections.split(";")
                                sections_edited = []
                                for stn in sections_listed:
                                    if '[' in stn and ']' in stn:
                                        sections_edited.append(stn)
                                    else:
                                        sections_edited.append(stn + "[1]")

                                ins_seq_in_full = []
                                for each_stn in sections_edited:
                                    bases, count = each_stn.split("[")
                                    count = int(count.replace("]", ""))
                                    for i in range(count):
                                        ins_seq_in_full.append(bases)

                            else:
                                bases, count = my_variant.quibble.split("[")
                                bases = bases.split(match)[1]
                                count = int(count.replace("]", ""))
                                ins_seq_in_full = []
                                for i in range(count):
                                    ins_seq_in_full.append(bases)
                            ins_seq_in_full = "".join(ins_seq_in_full)
                            vt_in_full = my_variant.quibble.split(match)[0] + match + ins_seq_in_full
                            warn = "%s may also be written as %s" % (my_variant.quibble, vt_in_full)
                            my_variant.warnings.append(warn)

                            try:
                                if "_" not in my_variant.quibble.split(":")[1] and \
                                        "del" not in my_variant.quibble.split(":")[1]:
                                    my_variant.warnings.append("An insertion must be provided with the two positions "
                                                               "between which the insertion has taken place")
                            except IndexError:
                                pass

                            # Mark the variant to not be written and re-submit for validation
                            my_variant.write = False
                            query = Variant(my_variant.original, quibble=vt_in_full,
                                            warnings=my_variant.warnings, primary_assembly=my_variant.primary_assembly,
                                            order=my_variant.order)

                            self.batch_list.append(query)
                            continue

                    # Format expanded repeat syntax
                    """
                    Waiting for HGVS nomenclature changes
                    """
                    try:
                        toskip = expanded_repeats.convert_tandem(my_variant, self, my_variant.primary_assembly,
                                                                 "all")
                    except expanded_repeats.RepeatSyntaxError as e:
                        my_variant.warnings = [str(e)]
                        continue
                    except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                        my_variant.warnings = ["HgvsSyntaxError: " + str(e)]
                        continue
                    except vvhgvs.exceptions.HGVSDataNotAvailableError as e:
                        if "invalid coordinates:" in str(e):
                            my_variant.warnings = [(f"ExonBoundaryError: Stated position "
                                                    f"does not correspond with an exon boundary for "
                                                    f"transcript {my_variant.quibble.split(':')[0]}")]
                            continue
                    except Exception as e:
                        my_variant.warnings = ["ExpandedRepeatError: " + str(e)]
                        continue

                    if toskip:
                        if my_variant.quibble != my_variant.expanded_repeat["variant"]:
                            my_variant.warnings.append(f"ExpandedRepeatWarning: {my_variant.quibble} updated "
                                                       f"to {my_variant.expanded_repeat['variant']}")
                        ins_bases = (my_variant.expanded_repeat["repeat_sequence"] *
                                     int(my_variant.expanded_repeat["copy_number"]))
                        repeat_to_delins = self.hp.parse(f"{my_variant.expanded_repeat['reference']}:"
                                                         f"{my_variant.expanded_repeat['prefix']}."
                                                         f"{my_variant.expanded_repeat['position']}"
                                                         f"delins{ins_bases}")

                        try:
                            repeat_to_delins = my_variant.hn.normalize(repeat_to_delins)
                        except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                            pass
                        my_variant.quibble = fn.valstr(repeat_to_delins)
                        my_variant.warnings.append(f"ExpandedRepeatWarning: {my_variant.expanded_repeat['variant']} "
                                                   f"should only be used as an annotation for the core "
                                                   f"HGVS descriptions provided")

                    # Methylation Syntax
                    methyl_syntax.methyl_syntax(my_variant)

                    # Set some configurations
                    formatted_variant = my_variant.quibble
                    stash_input = my_variant.quibble
                    my_variant.post_format_conversion = stash_input
                    logger.debug("Variant input formatted, proceeding to validate.")

                    # Conversions
                    # are not currently supported. The HGVS format for conversions
                    # is rarely seen wrt genomic sequencing data and needs to be re-evaluated
                    if 'con' in my_variant.quibble:
                        my_variant.warnings.append('Conversions are no longer valid HGVS Sequence Variant Descriptions')
                        logger.warning('Conversions are no longer valid HGVS Sequence Variant Descriptions')
                        continue

                    # Change RNA bases to upper case but nothing else
                    if my_variant.reftype == ":r.":
                        query_r_var = formatted_variant
                        formatted_variant = formatted_variant.upper()
                        formatted_variant = formatted_variant.replace(':R.', ':r.')

                        # lowercase the supported variant types
                        formatted_variant = formatted_variant.replace('DEL', 'del')
                        formatted_variant = formatted_variant.replace('INS', 'ins')
                        formatted_variant = formatted_variant.replace('INV', 'inv')
                        formatted_variant = formatted_variant.replace('DUP', 'dup')
                        ref, edit_ori = formatted_variant.split(":r.")
                        edit = copy.copy(edit_ori)
                        edit = edit.replace("G", "g")
                        edit = edit.replace("A", "a")
                        edit = edit.replace("T", "t")
                        edit = edit.replace("C", "c")
                        edit = edit.replace("U", "u")
                        formatted_variant = ref + ":r." + edit
                        if query_r_var != formatted_variant:
                            e = "This not a valid HGVS description, due to characters being in the wrong case. " \
                                "Please check the use of upper- and lowercase characters."
                            my_variant.warnings.append(str(e))
                            logger.warning(str(e))
                        formatted_variant = formatted_variant.replace(edit, edit_ori)

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

                        # Pass over for uncertain positions
                        if my_variant.reformat_output == "uncertain_pos":
                            continue

                        # Check for common mistakes
                        toskip = use_checking.refseq_common_mistakes(my_variant)
                        if toskip:
                            continue

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
                            error = "Invalid amino acid %s stated in description %s" % (str(e),
                                                                                        my_variant.quibble)
                            my_variant.warnings.append(error)
                            continue

                    my_variant.set_quibble(str(my_variant.hgvs_formatted))
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

                    # Fuzzy ends
                    try:
                        complex_descriptions.fuzzy_ends(my_variant, self)
                    except complex_descriptions.FuzzyPositionError as e:
                        my_variant.warnings.append(str(e))
                        logger.warning(str(e))
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
                        if "NC_012920.1" not in my_variant.hgvs_formatted.ac and \
                                "NC_001807.4" not in my_variant.hgvs_formatted.ac:
                            error = 'Interval end position ' + str(my_variant.hgvs_formatted.posedit.pos.end.base) + \
                                    ' < interval start position ' + str(my_variant.hgvs_formatted.posedit.pos.start.base)
                            my_variant.warnings.append(error)
                            logger.warning(error)
                            continue

                    # Catch missing version number in refseq
                    is_version = re.compile(r"\d\.\d")
                    if ((my_variant.refsource == 'RefSeq' or my_variant.refsource == 'ENS') and
                            not is_version.search(str(my_variant.hgvs_formatted))):
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
                                    % (my_variant.quibble, my_variant.primary_assembly)
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
                                try:
                                    hgvs_coding_from_refseqgene = self.vm.g_to_t(hgvs_refseqgene_variant,
                                                                             hgvs_transcript_variant.ac)
                                except vvhgvs.exceptions.HGVSInvalidIntervalError:
                                    hgvs_coding_from_refseqgene = hgvs_transcript_variant
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
                    multi_gen_vars = mappers.final_tx_to_multiple_genomic(variant,
                                                                          self,
                                                                          tx_variant,
                                                                          liftover_level=liftover_level)

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
                                    vcf_dict = hgvs_utils.report_hgvs2vcf(alt_gen_var, 'hg38'
                                                                          , variant.reverse_normalizer,
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
                                    vcf_dict = hgvs_utils.report_hgvs2vcf(alt_gen_var, 'hg38',
                                                                          variant.reverse_normalizer,
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

                # Clean up mito genome issues
                cp_lifted_response = copy.deepcopy(primary_genomic_dicts)
                for key, val in cp_lifted_response.items():
                    if key == "hg19" and "NC_012920.1" in val["hgvs_genomic_description"]:
                        primary_genomic_dicts.pop(key)
                    elif key == "grch37" and "NC_001807.4" in val["hgvs_genomic_description"]:
                        primary_genomic_dicts.pop(key)

                # Warn not directly mapped to specified genome build
                if genomic_accession != '':
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
                if refseqgene_variant != '':
                    try:
                        refseqgene_variant = fn.valstr(hgvs_refseqgene_variant)
                    except Exception as e:
                        logger.debug("Except passed, %s", e)
                    if variant.gene_symbol == "" and refseqgene_variant != "":
                        gene_symbol = self.db.get_gene_symbol_from_refseq_id(refseqgene_variant.split(":")[0])
                        variant.gene_symbol = gene_symbol

                # Add predicted protein variant dictionary
                if predicted_protein_variant != '':
                    predicted_protein_variant_dict = {}
                    predicted_protein_variant_dict["slr"] = ''
                    predicted_protein_variant_dict["tlr"] = ''
                    predicted_protein_variant_dict["lrg_tlr"] = ''
                    predicted_protein_variant_dict["lrg_slr"] = ''
                    if 'Non-coding :n.' not in predicted_protein_variant:
                        add_p_descps = True
                        try:
                            if "N" in str(hgvs_tx_variant.posedit.edit):
                                add_p_descps = False
                        except AttributeError:
                            pass
                        if add_p_descps is True:
                            try:
                                # Add single letter AA code to protein descriptions
                                predicted_protein_variant_dict = {"tlr": str(predicted_protein_variant), "slr": ''}
                                if re.search('p.=', predicted_protein_variant_dict['tlr']) \
                                        or re.search('p.?', predicted_protein_variant_dict['tlr']):
                                    # Replace p.= with p.(=)
                                    predicted_protein_variant_dict['tlr'] = \
                                        predicted_protein_variant_dict['tlr'].replace(
                                        'p.=',
                                        'p.(=)')

                                # Remove LRG
                                format_p = predicted_protein_variant_dict['tlr']

                                if 'LRG' in format_p:
                                    format_lrg = copy.copy(format_p)
                                    format_p = re.sub(r'\(LRG_.+?\)', '', format_p)
                                    if "(" in format_lrg:
                                        format_lrg = format_lrg.split('(')[1]
                                        format_lrg = format_lrg.replace(')', '')
                                else:
                                    format_lrg = None
                                    pass

                                if re.search("[A-Z][a-z][a-z]1[A-Z][a-z][a-z]", format_p):
                                    cp_warnings = []
                                    for each_warning in variant.warnings:
                                        if "is HGVS compliant and contains a valid reference " \
                                           "amino acid description" not in each_warning:
                                            cp_warnings.append(each_warning)
                                        else:
                                            cp_format_p = copy.copy(format_p)
                                            cp_format_p = cp_format_p.split(":")[0]
                                            aa_1 = self.sf.fetch_seq(cp_format_p, start_i=0, end_i=1)
                                            aa_1 = fn.one_to_three(aa_1)
                                            cp_format_p = f"{cp_format_p}:p.({aa_1}1?)"
                                            cp_warnings.append(f"Variant {format_p} affects the initiation amino acid"
                                                               f" so is better described as {cp_format_p}")
                                            format_p = cp_format_p
                                            variant.warnings = cp_warnings

                                re_parse_protein = self.hp.parse_hgvs_variant(format_p)

                                # Set formatted tlr
                                predicted_protein_variant_dict['tlr'] = str(copy.copy(re_parse_protein))
                                re_parse_protein_single_aa = fn.single_letter_protein(re_parse_protein)

                                # Replace p.= with p.(=)
                                if re.search('p.=', re_parse_protein_single_aa
                                             ) or re.search('p.?', re_parse_protein_single_aa):
                                    re_parse_protein_single_aa = re_parse_protein_single_aa.replace('p.=',
                                                                                                    'p.(=)')
                                if re.search("p.\(Ter\d+=\)", predicted_protein_variant_dict['tlr']):
                                    predicted_protein_variant_dict['tlr'] = \
                                        predicted_protein_variant_dict['tlr'].split("p.")[0]
                                    predicted_protein_variant_dict['tlr'] = \
                                        predicted_protein_variant_dict['tlr'] + "p.(Ter=)"
                                    re_parse_protein_single_aa = re_parse_protein_single_aa.split("p.")[0]
                                    re_parse_protein_single_aa = re_parse_protein_single_aa + "p.(*=)"

                                elif re.search("p.Ter\d+=", predicted_protein_variant_dict['tlr']):
                                    predicted_protein_variant_dict['tlr'] = \
                                        predicted_protein_variant_dict['tlr'].split("p.")[0]
                                    predicted_protein_variant_dict['tlr'] = \
                                        predicted_protein_variant_dict['tlr'] + "p.Ter="
                                    re_parse_protein_single_aa = re_parse_protein_single_aa.split("p.")[0]
                                    re_parse_protein_single_aa = re_parse_protein_single_aa + "p.*="

                                # Capture instances of variation affecting p.1
                                if re.search("[A-Z][a-z[a-z]1[?]", predicted_protein_variant_dict['tlr']):
                                    match = re.search("([A-Z][a-z][a-z]1[?])", predicted_protein_variant_dict['tlr'])
                                    captured_aa = match.group(1).split("1")[0]
                                    captured_aa_sl = fn.three_to_one(captured_aa)
                                    predicted_protein_variant_dict['tlr'] = \
                                        predicted_protein_variant_dict['tlr'].split("p.")[0]
                                    predicted_protein_variant_dict['tlr'] = \
                                        predicted_protein_variant_dict['tlr'] + "p.(" + captured_aa + "1?)"
                                    re_parse_protein_single_aa = re_parse_protein_single_aa.split("p.")[0]
                                    re_parse_protein_single_aa = re_parse_protein_single_aa + "p.(" \
                                                                                              + captured_aa_sl + "1?)"
                                predicted_protein_variant_dict["slr"] = str(re_parse_protein_single_aa)

                                # set LRG outputs
                                if format_lrg is not None:
                                    predicted_protein_variant_dict["lrg_tlr"] = \
                                        format_lrg.split(':')[0] + ':' + \
                                        predicted_protein_variant_dict["tlr"].split(':')[1]
                                    predicted_protein_variant_dict["lrg_slr"] = \
                                        format_lrg.split(':')[0] + ':' + \
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
                if variant.hgvs_transcript_variant is not None and variant.gene_symbol == '':
                    variant.gene_symbol = self.db.get_gene_symbol_from_transcript_id(
                        variant.hgvs_transcript_variant.split(":")[0])
                elif variant.hgvs_refseqgene_variant is not None and variant.gene_symbol == '':
                    variant.gene_symbol = self.db.get_gene_symbol_from_refseq_id(
                        variant.hgvs_refseqgene_variant.split(":")[0])

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
                    variant.hgvs_transcript_variant = tx_variant

                if variant.hgvs_transcript_variant is None:
                    variant.hgvs_transcript_variant = tx_variant
                variant.reference_sequence_records = ''
                variant.validated = True

                # Add links to reference_sequence_records
                ref_records = self.db.get_urls(variant.output_dict())
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
                            if (genomic_position_info[g_p_key]['hgvs_genomic_description'] not in lo_cache.keys()) or (
                                    "NC_012920.1" in genomic_position_info[g_p_key]['hgvs_genomic_description']
                                    and build_from == "hg38" and build_to == "hg19"):

                                lifted_response = liftover(genomic_position_info[g_p_key]['hgvs_genomic_description'],
                                                           build_from,
                                                           build_to, variant.hn, variant.reverse_normalizer,
                                                           variant.evm,
                                                           self,
                                                           specify_tx=False,
                                                           liftover_level=liftover_level,
                                                           g_to_g=g_to_g)

                                if "NC_012920.1" in genomic_position_info[g_p_key]['hgvs_genomic_description'] or \
                                        "NC_001807.4:" in genomic_position_info[g_p_key]['hgvs_genomic_description']:
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

                # Add exon numbering information see issue #
                if variant.coding != "":
                    try:
                        exs = exon_numbering.finds_exon_number(variant, self)
                        variant.exonic_positions = exs
                    except KeyError:
                        pass

                # Remove duplicate warnings
                variant_warnings = []
                accession = variant.hgvs_transcript_variant.split(':')[0]
                term = str(accession)
                term_2 = "%s automapped to" % tx_variant
                term_3 = "%s automapped to" % genomic_variant
                for vt in variant.warnings:
                    vt = str(vt)


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
                    elif (term_2 in vt and tx_variant != "") or (term_3 in vt and genomic_variant != ""):
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

                # Reformat as required
                if variant.reformat_output is not False:
                    if "|" in variant.reformat_output and "=" in variant.quibble:
                        attributes = dir(variant)
                        for attribute in attributes:
                            if "__" in attribute:
                                continue
                            item = (getattr(variant, attribute))
                            if isinstance(item, str):
                                if ("p." not in item and "=" in item) and attribute != "reformat_output":
                                    setattr(variant, attribute, item.replace("=", variant.reformat_output))
                            if isinstance(item, dict) or isinstance(item, list):
                                try:
                                    stringy = json.dumps(item)
                                    if re.search(r":[gcnr].", stringy) and "=" in stringy:
                                        stringy = stringy.replace("=", variant.reformat_output)
                                        setattr(variant, attribute, json.loads(stringy))
                                except json.JSONDecodeError:
                                    pass

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
                         transcript_set="refseq", genome_build=None, batch_output=False):

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
                                                         select_transcripts, transcript_set, genome_build)
        else:
            g2d_data = []
            for symbol in gene_symbols:
                data_for_gene = gene2transcripts.gene2transcripts(self, symbol, validator, bypass_web_searches,
                                                                  select_transcripts, transcript_set, genome_build)
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
