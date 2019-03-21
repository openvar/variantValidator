# -*- coding: utf-8 -*-
"""
VariantValidator.py
List of top level VariantValidator functions
This API is configured by reading the configuration information in the config.ini file
located in the configuration module contained the the root variantValidator directory.
These configurations can be over-ridden by setting environment variables, see
README.txt
This version of the VariantValidator API 0.1.0 contains the following functions:
1 my_config
my_config is a simple function that allows the user to determine whether VariantValidator is
correctly configured, i.e. is the tool searching in the correct locations for its data?
The function also returns version information
# Example
my_config()
2. validator
validator is the primary VariantValidator function which validates sequence variation
descriptions. validator uses sub functions in the variantanalyser module
contained the the root variantValidator directory, and functions priovided by the
hgvs Python package (https://github.com/biocommons/hgvs/) to "manipulate biological
sequence variants according to Human Genome Variation Society recommendations"
# Example
variant = ' NM_000088.3:c.589G>T'
selected_assembly = 'GRCh37' # or GRCh37, hg19, hg38
select_transcripts = 'all' # Or a pipe delimited, white-space-less, string of transcript
IDs validation = validator(variant, selected_assembly, select_transcripts)
# Accepted input formats
NM_000088.3:c.589G>T
NC_000017.10:g.48275363C>A
NG_007400.1:g.8638G>T
LRG_1:g.8638G>T
LRG_1t1:c.589G>T
17-50198002-C-A (GRCh38)
chr17:50198002C>A (GRCh38)
3. gene2transcripts
This function is similar to the Gene to Transcripts function
https://variantvalidator.org/ref_finder/ except the data is returned within a structured
python object
# HGNC example
variantValidator.validator.gene2transcripts('HTT')
# RefSeq Transcript example
variantValidator.validator.gene2transcripts('NM_002111.8')
4. hgvs2ref
This function retuns the reference sequence with respect to HGVS variation descriptions
The function will only return REFERENCE SEQUENCE i.e. if a c. descriptions overlaps an
intron/exon boundary, only the exonic sequence will be returned
# Example
hgvs2ref('NM_000088.3:c.589_594del')
"""

# IMPORT HGVS MODULES
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
import datetime
import copy
import os
import sys
import warnings
from operator import itemgetter
from pyliftover import LiftOver
import traceback
from configparser import ConfigParser

from Bio.Seq import Seq

# Import variantanalyser and peripheral VV modules
import ref_seq_type
import external
import output_formatter
import variantanalyser
from variantanalyser.vvLogging import logger
from variantanalyser import functions as va_func
from variantanalyser import dbControls as va_dbCrl
from variantanalyser import hgvs2vcf as va_H2V
from variantanalyser import batch as va_btch
from variantanalyser import g_to_g as va_g2g
from variantanalyser import supported_chromosome_builds as va_scb
from variantanalyser import gap_genes as gapGenes
from variantanalyser.liftover import liftover as lift_over

__version__ = None

# Ensure configuration is on the OS
if os.environ.get('CONF_ROOT') is None:
    import configuration

    CONF_ROOT = os.environ.get('CONF_ROOT')
else:
    CONF_ROOT = os.environ.get('CONF_ROOT')
# Define global configuration variables
HGVS_SEQREPO_DIR = "Unspecified"
HGVS_SEQREPO_DIR = "Unspecified"
UTA_DB_URL = 'Unspecified'
VALIDATOR_DB_URL = 'Unspecified'
PYLIFTOVER_DIR = 'Unspecified'
ENTREZ_ID = 'Unspecified'
VERSION = 'Unspecified'
hgvs_version = 'Unspecified'


def exceptPass(validation=None):
    exc_type, exc_value, last_traceback = sys.exc_info()
    te = traceback.format_exc()
    tbk = [str(exc_type), str(exc_value), str(te)]
    er = str('\n'.join(tbk))
    if last_traceback:
        logger.warning(
            "Except pass for " + str(exc_type) + " " + str(exc_value) + " at line " + str(last_traceback.tb_lineno))
    else:
        logger.warning("Except pass for " + str(exc_type) + " " + str(exc_value))
    logger.debug(er)


# Config Section Mapping function
def ConfigSectionMap(section, c):
    dict1 = {}
    options = c.options(section)
    for option in options:
        try:
            dict1[option] = c.get(section, option)
            if dict1[option] == -1:
                logger.warning("skip: %s" % option)
        except:
            logger.warning("exception on %s!" % option)
            dict1[option] = None
    return dict1


def loadConfigFile():
    # Set the version from the config.ini
    global __version__
    Config = ConfigParser()
    with open(os.path.join(CONF_ROOT, 'config.ini')) as f:
        Config.read_file(f)
    __version__ = ConfigSectionMap("variantValidator",Config)['version']
    if re.match('^\d+\.\d+\.\d+$', __version__) is not None:
        _is_released_version = True
    # Load database environments from config
    levelString = ConfigSectionMap("logging", Config)['level']
    consoleString = ConfigSectionMap("logging", Config)['console']
    if consoleString.lower()=="true":
        consoleString="console"
    fileString = ConfigSectionMap("logging", Config)['file']
    if fileString.lower()=="true":
        fileString="file"
    traceString = ConfigSectionMap("logging", Config)['trace']
    if traceString.lower()=="true":
        traceString="trace"
    logString = levelString+" "+consoleString+" "+fileString+" "+traceString
    os.environ["VALIDATOR_DEBUG"] = logString
    #print "ac", os.environ["VALIDATOR_DEBUG"]
    #print("ls", logString)


# Custom Exceptions
class variantValidatorError(Exception):
    pass


# PRE COMPILE VARIABLES
hgvs.global_config.uta.pool_max = 25
hdp = hgvs.dataproviders.uta.connect(pooling=True)
# From the hgvs parser import, create an instance of hgvs.parser.Parser
hp = hgvs.parser.Parser()
# Configure hgvs package global settings
hgvs.global_config.formatting.max_ref_length = 1000000
# Validator
vr = hgvs.validator.Validator(hdp)
# Variant mapper
vm = hgvs.variantmapper.VariantMapper(hdp)
# Create a lose vm instance
lose_vm = hgvs.variantmapper.VariantMapper(hdp,
                                           replace_reference=True,
                                           prevalidation_level=None
                                           )
nr_vm = hgvs.variantmapper.VariantMapper(hdp, replace_reference=False)
# Create seqfetcher object
sf = hgvs.dataproviders.seqfetcher.SeqFetcher()

# Set current genome builds
genome_builds = ['GRCh37', 'hg19', 'GRCh38']

# Obtain environment variables needed within the top-level function
PYLIFTOVER_DIR = os.environ.get('PYLIFTOVER_DIR')


# method for final validation and stringifying parsed hgvs variants prior to printing/passing to html
def valstr(hgvs_variant):
    """
    Function to ensure the required number of reference bases are displayed in descriptions
    """
    cp_hgvs_variant = copy.deepcopy(hgvs_variant)
    if cp_hgvs_variant.posedit.edit.type == 'identity':
        if len(cp_hgvs_variant.posedit.edit.ref) > 1:
            cp_hgvs_variant = output_formatter.remove_reference(cp_hgvs_variant)
        cp_hgvs_variant = str(cp_hgvs_variant)
    else:
        cp_hgvs_variant = output_formatter.remove_reference(cp_hgvs_variant)
        cp_hgvs_variant = str(cp_hgvs_variant)
    return cp_hgvs_variant


# Check configuration variables
def my_config():
    loadConfigFile()
    global HGVS_SEQREPO_DIR
    global UTA_DB_URL
    global VALIDATOR_DB_URL
    global PYLIFTOVER_DIR
    global ENTREZ_ID
    global VERSION
    global hgvs_version
    try:
        HGVS_SEQREPO_DIR = os.environ.get('HGVS_SEQREPO_DIR')
    except Exception:
        HGVS_SEQREPO_DIR = 'Unspecified'
    try:
        UTA_DB_URL = os.environ.get('UTA_DB_URL')
    except Exception:
        UTA_DB_URL = 'Unspecified'
    VALIDATOR_DB_URL = os.environ.get('VALIDATOR_DB_URL')
    try:
        PYLIFTOVER_DIR = os.environ.get('PYLIFTOVER_DIR')
    except Exception:
        PYLIFTOVER_DIR = 'Unspecified'
    ENTREZ_ID = os.environ.get('ENTREZ_ID')
    VERSION = __version__,
    VERSION = str(VERSION[0])
    hgvs_version = hgvs.__version__,
    hgvs_version = str(hgvs_version[0])
    locate = {
        'variantvalidator_version': VERSION, #
        'variantvalidator_hgvs_version': hgvs_version, #
        'uta_schema': str(hdp.data_version()), #
        'seqrepo_db': HGVS_SEQREPO_DIR.split('/')[-1] #
    }
    return locate


# Validator code
""" 
This is the primary VariantValidator function
"""


def validator(batch_variant, selected_assembly, select_transcripts, transcriptSet="refseq"):
    logger.info(batch_variant + ' : ' + selected_assembly)
    # Take start time
    start_time = time.time()

    # Set pre defined variables
    # SeqFetcher
    # sf = hgvs.dataproviders.seqfetcher.SeqFetcher()

    try:
        # Validation
        ############

        # Create a dictionary of transcript ID : ''
        if select_transcripts != 'all':
            select_transcripts_list = select_transcripts.split('|')
            select_transcripts_dict = {}
            select_transcripts_dict_plus_version = {}
            for id in select_transcripts_list:
                id = id.strip()
                if re.match('LRG', id):
                    id = va_dbCrl.data.get_RefSeqTranscriptID_from_lrgTranscriptID(id)
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
        batch_list = []
        for queries in batch_queries:
            queries = queries.strip()
            query = {'quibble': queries, 'id': queries, 'warnings': '', 'description': '', 'coding': '', 'coding_g': '',
                     'genomic_r': '', 'genomic_g': '', 'protein': '', 'write': 'true', 'primary_assembly': 'false',
                     'order': 'false'}
            batch_list.append(query)

        # Create List to carry batch data
        batch_out = []

        # Ensure batch_list is pulled into the function so that it can be appended to
        batch_list = batch_list

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
        set_output_type_flag = 'warning'
        logger.debug("Batch list length " + str(len(batch_list)))
        for validation in batch_list:
            # Start timing
            logger.traceStart(validation)
            # Re-set cautions and automaps

            if transcriptSet == "refseq":
                alt_aln_method = 'splign'
            elif transcriptSet == "ensembl":
                alt_aln_method = 'genebuild'
                logger.warning("Ensembl is currently not supported")
                validation['warnings'] += ': ' + "Ensembl is currently not supported"
                continue
            else:
                logger.warning(
                    "The transcript set variable " + transcriptSet + " is invalid, it needs to be 'refseq' or 'ensembl'")
                validation[
                    'warnings'] += ': ' + "The transcript set variable " + transcriptSet + " is invalid, it needs to be 'refseq' or 'ensembl'"
                continue

            # Create Normalizers
            hn = hgvs.normalizer.Normalizer(hdp,
                                            cross_boundaries=False,
                                            shuffle_direction=3,
                                            alt_aln_method=alt_aln_method
                                            )
            reverse_normalizer = hgvs.normalizer.Normalizer(hdp,
                                                            cross_boundaries=False,
                                                            shuffle_direction=5,
                                                            alt_aln_method=alt_aln_method
                                                            )

            # Blank cautions
            caution = ''
            automap = ''

            # This will be used to order the final output
            if str(validation['order']) == 'false':
                ordering = ordering + 1
                validation['order'] = ordering
            else:
                pass
            # Bug catcher
            try:
                # Note, ID is not touched. It is always the input variant description. Quibble will be altered but id will not if type = g.
                input = validation['quibble']
                logger.trace("Commenced validation of " + str(input), validation)

                # Test for rich text unicode characters
                try:
                    unicode_test = u"{}".format(input)
                except UnicodeDecodeError as e:
                    # Format the trapped character into unicode for styled printing
                    my_unicode = e[1]
                    my_unicode = my_unicode.decode('utf-8')

                    # Test for rich text unicode characters
                    try:
                        str(my_unicode)
                    except UnicodeEncodeError as e:
                        # Format the trapped character into unicode for styled printing
                        unicoded_it = e[1]
                        unicoded_it_list = unicoded_it.split()
                        for try_me in unicoded_it_list:
                            try:
                                str(try_me)
                            except UnicodeEncodeError as e:
                                found_unicode = try_me
                                found_error = str(e)
                                found_at = found_unicode.encode('raw_unicode_escape')
                                break
                        # Extract character from the error
                        unicode = re.findall("u'\\\\\w+'", found_error)
                        character = unicode[0]
                        search_term = character.replace("u'", '')
                        search_term = search_term.replace("'", '')
                        found_at_decoded = found_at.decode('raw_unicode_escape')
                        found_at = found_at_decoded.encode('raw_unicode_escape')
                        string_char = str(character)
                        # Create a human readable U+ representation
                        human_code = re.sub("u'\\\\\w", 'U+', string_char)
                        human_code = human_code.replace("'", "")
                        format_human = u"{}".format(human_code)
                        format_human = format_human.upper()
                        found_at = re.sub(search_term, u'<' + format_human + u'>', found_at)
                        slasher = re.compile("\\\\")
                        found_at = re.sub(slasher, '', found_at)
                        validation['id'] = found_at
                        error = u'Submitted variant description contains an invalid character which is represented by Unicode character ' + format_human + u' at position ' + found_at + u': Please remove this character and re-submit: A useful search function for Unicode characters can be found at https://unicode-search.net/'
                        validation['warnings'] = validation['warnings'] + ': ' + error
                        logger.warning(error)
                        continue
                    else:
                        pass
                else:
                    pass

                # Remove whitespace
                ws = copy.copy(input)
                input = input.strip()
                input = ''.join(input.split())
                if input != ws:
                    caution = 'Whitespace removed from variant description ' + str(ws)
                    validation['warnings'] = validation['warnings'] + ': ' + caution
                    logger.info(caution)
                stash_input = copy.copy(input)
                # Set the primary_assembly
                if validation['primary_assembly'] == 'false':
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
                    valid_build = False
                    for genome_build in genome_builds:
                        if selected_assembly == genome_build:
                            valid_build = True
                    if valid_build is False:
                        selected_assembly = 'GRCh38'
                        primary_assembly = selected_assembly
                        if 'GRCh38' in input or 'GRCh37' in input or 'hg19' in input or 'hg38' in input:
                            pass
                        else:
                            validation['warnings'] = validation[
                                                     'warnings'] + ': Invalid genome build has been specified. Automap has selected the default build (GRCh38)'
                            logger.warning(
                            'Invalid genome build has been specified. Automap has selected the default build ' + primary_assembly)
                        validation['primary_assembly'] = primary_assembly
                    else:
                        validation['primary_assembly'] = primary_assembly
                else:
                    primary_assembly = validation['primary_assembly']
                logger.trace("Completed string formatting", validation)
                # Set variables that batch will not use but are required
                crossing = 'false'
                boundary = 'false'

                # VCF type 1
                """
                VCF2HGVS stage 1. converts chr-pos-ref-alt into chr:posRef>Alt
                The output format is a common mistake caused by inaccurate conversion of 
                VCF variants into HGVS - hence the need for conversion step 2
                """
                if re.search('[-:]\d+[-:][GATC]+[-:][GATC]+', input):
                    input = input.replace(':', '-')
                    # Extract primary_assembly if provided
                    if re.match('GRCh3\d+-', input) or re.match('hg\d+-', input):
                        in_list = input.split('-')
                        selected_assembly = in_list[0]
                        input = '-'.join(in_list[1:])
                    pre_input = copy.deepcopy(input)
                    vcf_elements = pre_input.split('-')
                    input = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], vcf_elements[3])
                elif re.search('[-:]\d+[-:][GATC]+[-:]', input):
                    input = input.replace(':', '-')
                    # Extract primary_assembly if provided
                    if re.match('GRCh3\d+-', input) or re.match('hg\d+-', input):
                        in_list = input.split('-')
                        selected_assembly = in_list[0]
                        input = '-'.join(in_list[1:])
                    pre_input = copy.deepcopy(input)
                    vcf_elements = pre_input.split('-')
                    validation[
                        'warnings'] = 'Not stating ALT bases is ambiguous because VCF specification 4.0 would treat ' + pre_input + ' as a deletion whereas VCF specification 4.1 onwards would treat ' + pre_input + ' as ALT = REF'
                    validation['warnings'] = validation['warnings'] + ': VariantValidator has output both alternatives'
                    logger.resub('Not stating ALT bases is ambiguous because VCF specification 4.0 would treat ' +
                                 pre_input + ' as a deletion whereas VCF specification 4.1 onwards would treat ' + pre_input +
                                 ' as ALT = REF. Validator will output both alternatives.')
                    validation['write'] = 'false'
                    input_A = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], 'del')
                    input_B = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], vcf_elements[2])
                    queryA = {'quibble': input_A, 'id': validation['id'], 'warnings': validation['warnings'],
                              'description': '', 'coding': '', 'coding_g': '', 'genomic_r': '', 'genomic_g': '',
                              'protein': '', 'write': 'true', 'primary_assembly': primary_assembly, 'order': ordering}
                    queryB = {'quibble': input_B, 'id': validation['id'], 'warnings': validation['warnings'],
                              'description': '', 'coding': '', 'coding_g': '', 'genomic_r': '', 'genomic_g': '',
                              'protein': '', 'write': 'true', 'primary_assembly': primary_assembly, 'order': ordering}
                    batch_list.append(queryA)
                    batch_list.append(queryB)
                    continue
                elif re.search('[-:]\d+[-:][-:][GATC]+', input) or re.search('[-:]\d+[-:][.][-:][GATC]+', input):
                    input = input.replace(':', '-')
                    if re.search('-.-', input):
                        input = input.replace('-.-', '-ins-')
                    if re.search('--', input):
                        input = input.replace('--', '-ins-')
                    # Extract primary_assembly if provided
                    if re.match('GRCh3\d+-', input) or re.match('hg\d+-', input):
                        in_list = input.split('-')
                        selected_assembly = in_list[0]
                        input = '-'.join(in_list[1:])
                    pre_input = copy.deepcopy(input)
                    vcf_elements = pre_input.split('-')
                    input = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], vcf_elements[3])
                    stash_input = input
                logger.trace("Completed VCF-HVGS step 1", validation)
                # API type non-HGVS
                # e.g. Chr16:2099572TC>T
                """
                VCF2HGVS conversion step 2 identifies the correct chromosomal reference 
                sequence based upon the non compliant identifier e.g. <Chr16>:2099572TC>T.
                The data is currently stored in variantanalyser.supported_chromosome_builds. 
                Anticipated future builds will be transferred to MySQL which can be more 
                easily updated and maintained.
                LRGs and LRG_ts also need to be assigned the correct reference sequence identifier.
                The LRG ID data ia stored in the VariantValidator MySQL database.
                The reference sequence type is also assigned. 
                """
                if re.search('\w+\:', input) and not re.search('\w+\:[gcnmrp]\.', input):
                    if re.search('\w+\:[gcnmrp]', input) and not re.search('\w+\:[gcnmrp]\.', input):
                        # Missing dot
                        pass
                    else:
                        try:
                            if re.search('GRCh37', input) or re.search('hg19', input):
                                primary_assembly = 'GRCh37'
                            elif re.search('GRCh38', input) or re.search('hg38', input):
                                primary_assembly = 'GRCh38'
                            pre_input = copy.deepcopy(input)
                            input_list = input.split(':')
                            pos_ref_alt = str(input_list[1])
                            positionAndEdit = input_list[1]
                            if not re.match('N[CGTWMRP]_', input) and not re.match('LRG_', input):
                                chr_num = str(input_list[0])
                                chr_num = chr_num.upper()
                                chr_num = chr_num.strip()
                                if re.match('CHR', chr_num):
                                    chr_num = chr_num.replace('CHR', '')
                                # Use selected assembly
                                accession = va_scb.to_accession(chr_num, selected_assembly)
                                if accession is None:
                                    validation['warnings'] = validation[
                                                                 'warnings'] + ': ' + chr_num + \
                                                             ' is not part of genome build ' + selected_assembly
                                    logger.warning(chr_num + ' is not part of genome build ' + selected_assembly)
                                    continue
                            else:
                                accession = input_list[0]
                            if re.search('>', pre_input):
                                if re.search('del', pre_input):
                                    pos = re.match('\d+', pos_ref_alt)
                                    position = pos.group(0)
                                    old_ref, old_alt = pos_ref_alt.split('>')
                                    old_ref = old_ref.replace(position, '')
                                    position = int(position) - 1
                                    required_base = sf.fetch_seq(accession, start_i=position - 1, end_i=position)
                                    ref = required_base + old_ref
                                    alt = required_base
                                    positionAndEdit = str(position) + ref + '>' + alt
                                elif re.search('ins', pre_input):
                                    pos = re.match('\d+', pos_ref_alt)
                                    position = pos.group(0)
                                    old_ref, old_alt = pos_ref_alt.split('>')
                                    # old_ref = old_ref.replace(position, '')
                                    position = int(position) - 1
                                    required_base = sf.fetch_seq(accession, start_i=position - 1, end_i=position)
                                    ref = required_base
                                    alt = required_base + old_alt
                                    positionAndEdit = str(position) + ref + '>' + alt
                            # Assign reference sequence type
                            ref_type = ref_seq_type.ref_type_assign(accession)
                            if re.match('LRG_', accession):
                                if ref_type == ':g.':
                                    accession = va_dbCrl.data.get_RefSeqGeneID_from_lrgID(accession)
                                else:
                                    accession = va_dbCrl.data.get_RefSeqTranscriptID_from_lrgTranscriptID(accession)
                            else:
                                accession = accession
                            input = str(accession) + ref_type + str(positionAndEdit)
                            stash_input = input
                        except:
                            exceptPass(validation)

                # Descriptions lacking the colon :
                if re.search('[gcnmrp]\.', input) and not re.search(':[gcnmrp]\.', input):
                    error = 'Unable to identify a colon (:) in the variant description %s. A colon is required in HGVS variant descriptions to separate the reference accession from the reference type i.e. <accession>:<type>. e.g. :c.' % (
                        input)
                    validation['warnings'] = validation['warnings'] + ': ' + error
                    logger.warning(error)
                    continue

                # Ambiguous chr reference
                logger.trace("Completed VCF-HVGS step 2", validation)
                """
                VCF2HGVS conversion step 3 is similar to step 2 but handles 
                formats like Chr16:g.2099572TC>T which are provided by Alamut and other
                software
                """
                if re.search('\w+:[gcnmrp]\.', input) and not re.match('N[CGTWMRP]_', input):
                    # Take out lowercase Accession characters
                    lower_cased_list = input.split(':')
                    if re.search('LRG', lower_cased_list[0], re.IGNORECASE):
                        lower_case_accession = lower_cased_list[0]
                        lower_case_accession = lower_case_accession.replace('l', 'L')
                        lower_case_accession = lower_case_accession.replace('r', 'R')
                        lower_case_accession = lower_case_accession.replace('g', 'G')
                    else:
                        lower_case_accession = lower_cased_list[0]
                        lower_case_accession = lower_case_accession.upper()
                    input = ''.join(lower_cased_list[1:])
                    input = lower_case_accession + ':' + input
                    if not re.match('LRG_', input) and not re.match('ENS', input) and not re.match('N[MRPC]_', input):
                        try:
                            if re.search('GRCh37', input) or re.search('hg19', input):
                                primary_assembly = 'GRCh37'
                            elif re.search('GRCh38', input) or re.search('hg38', input):
                                primary_assembly = 'GRCh38'
                            pre_input = copy.deepcopy(input)
                            input_list = input.split(':')
                            query_a_symbol = input_list[0]
                            is_it_a_gene = va_dbCrl.data.get_hgnc_symbol(query_a_symbol)
                            if is_it_a_gene == 'none':
                                pos_ref_alt = str(input_list[1])
                                positionAndEdit = input_list[1]
                                chr_num = str(input_list[0])
                                chr_num = chr_num.upper()
                                chr_num = chr_num.strip()
                                if re.match('CHR', chr_num):
                                    chr_num = chr_num.replace('CHR', '')  # Use selected assembly
                                accession = va_scb.to_accession(chr_num, selected_assembly)
                                if accession is None:
                                    validation['warnings'] = validation['warnings'] + ': ' + chr_num + \
                                                             ' is not part of genome build ' + selected_assembly
                                    continue
                                input = str(accession) + ':' + str(positionAndEdit)
                                stash_input = input
                            else:
                                pass
                        except Exception as e:
                            exc_type, exc_value, last_traceback = sys.exc_info()
                            te = traceback.format_exc()
                            tbk = [str(exc_type), str(exc_value), str(te)]
                            er = str('\n'.join(tbk))
                            logger.warning(str(exc_type) + " " + str(exc_value))
                            logger.debug(er)

                # GENE_SYMBOL:c. n. types
                logger.trace("Completed VCF-HGVS step 3", validation)
                """
                Searches for gene symbols that have been used as reference sequence
                identifiers. Provides a sufficiently repremanding warning, but also provides
                correctly formatted variant descriptions with appropriate transcript 
                reference sequence identifiers i.e. NM_ ....
                Note: the output from the function must be validated because VV has no way
                of knowing which the users intended reference sequence was, and the exon 
                boundaries etc of the alternative transcript variants may not be equivalent 
                """
                if re.search('\w+\:[cn]\.', input):
                    try:
                        pre_input = copy.deepcopy(input)
                        query_a_symbol = pre_input.split(':')[0]
                        tx_edit = pre_input.split(':')[1]
                        is_it_a_gene = va_dbCrl.data.get_hgnc_symbol(query_a_symbol)
                        if is_it_a_gene != 'none':
                            uta_symbol = va_dbCrl.data.get_uta_symbol(is_it_a_gene)
                            available_transcripts = hdp.get_tx_for_gene(uta_symbol)
                            select_from_these_transcripts = {}
                            for tx in available_transcripts:
                                if re.match('NM_', tx[3]) or re.match('NR_', tx[3]):
                                    if tx[3] not in select_from_these_transcripts.keys():
                                        select_from_these_transcripts[tx[3]] = ''
                                    else:
                                        continue
                                else:
                                    continue
                            select_from_these_transcripts = '|'.join(select_from_these_transcripts.keys())
                            if select_transcripts != 'all':
                                validation['write'] = 'false'
                                for transcript in select_transcripts_dict_plus_version.keys():
                                    validation[
                                        'warnings'] = 'HGVS variant nomenclature does not allow the use of a gene symbol (' + \
                                                      query_a_symbol + ') in place of a valid reference sequence'
                                    refreshed_description = transcript + ':' + tx_edit
                                    query = {'quibble': refreshed_description, 'id': validation['id'],
                                             'warnings': validation['warnings'], 'description': '', 'coding': '',
                                             'coding_g': '', 'genomic_r': '', 'genomic_g': '', 'protein': '',
                                             'write': 'true', 'primary_assembly': primary_assembly, 'order': ordering}
                                    batch_list.append(query)
                                    logger.resub('HGVS variant nomenclature does not allow the use of a gene symbol (' + \
                                                 query_a_symbol + ') in place of a valid reference sequence')
                            else:
                                validation['warnings'] = validation['warnings'] + \
                                                         ': ' + 'HGVS variant nomenclature does not allow the use of a gene symbol (' + \
                                                         query_a_symbol + ') in place of a valid reference sequence: Re-submit ' + input + \
                                                         ' and specify transcripts from the following: ' + 'select_transcripts=' + select_from_these_transcripts
                                logger.warning('HGVS variant nomenclature does not allow the use of a gene symbol (' + \
                                               query_a_symbol + ') in place of a valid reference sequence: Re-submit ' + input + \
                                               ' and specify transcripts from the following: ' + 'select_transcripts=' + select_from_these_transcripts)
                            continue
                        else:
                            pass
                    except:
                        exceptPass()
                logger.trace("Gene symbol reference catching complete", validation)

                # NG_:c. or NC_:c.
                """
                Similar to the GENE_SYMBOL:c. n. types function, but spots RefSeqGene or
                Chromosomal reference sequence identifiers used in the context of c. variant
                descriptions
                """
                if re.search('\w+\:[cn]', input):
                    try:
                        if re.match('^NG_', input):
                            refSeqGeneID = input.split(':')[0]
                            tx_edit = input.split(':')[1]
                            gene_symbol = va_dbCrl.data.get_gene_symbol_from_refSeqGeneID(refSeqGeneID)
                            if gene_symbol != 'none':
                                uta_symbol = va_dbCrl.data.get_uta_symbol(gene_symbol)
                                available_transcripts = hdp.get_tx_for_gene(uta_symbol)
                                select_from_these_transcripts = {}
                                for tx in available_transcripts:
                                    if re.match('NM_', tx[3]) or re.match('NR_', tx[3]):
                                        if tx[3] not in select_from_these_transcripts.keys():
                                            select_from_these_transcripts[tx[3]] = ''
                                        else:
                                            continue
                                    else:
                                        continue
                                select_from_these_transcripts = '|'.join(select_from_these_transcripts.keys())
                                if select_transcripts != 'all':
                                    validation['write'] = 'false'
                                    for transcript in select_transcripts_dict_plus_version.keys():
                                        validation[
                                            'warnings'] = 'NG_:c.PositionVariation descriptions should not be used unless a transcript reference sequence has also been provided e.g. NG_(NM_):c.PositionVariation'
                                        refreshed_description = refSeqGeneID + '(' + transcript + ')' + ':' + tx_edit
                                        query = {'quibble': refreshed_description, 'id': validation['id'],
                                                 'warnings': validation['warnings'], 'description': '', 'coding': '',
                                                 'coding_g': '', 'genomic_r': '', 'genomic_g': '', 'protein': '',
                                                 'write': 'true', 'primary_assembly': primary_assembly,
                                                 'order': ordering}
                                        logger.resub(
                                            'NG_:c.PositionVariation descriptions should not be used unless a transcript reference sequence has also been provided e.g. NG_(NM_):c.PositionVariation. Resubmitting corrected version.')
                                        batch_list.append(query)
                                else:
                                    validation['warnings'] = validation[
                                                                 'warnings'] + ': ' + 'A transcript reference sequence has not been provided e.g. NG_(NM_):c.PositionVariation. Re-submit ' + input + ' but also specify transcripts from the following: ' + 'select_transcripts=' + select_from_these_transcripts
                                    logger.warning(
                                        + 'A transcript reference sequence has not been provided e.g. NG_(NM_):c.PositionVariation. Re-submit ' +
                                        str(
                                            input) + ' but also specify transcripts from the following: ' + 'select_transcripts=' + str(
                                            select_from_these_transcripts))
                                continue
                            else:
                                validation['warnings'] = validation[
                                                             'warnings'] + ': ' + 'A transcript reference sequence has not been provided e.g. NG_(NM_):c.PositionVariation'
                                logger.warning(
                                    'A transcript reference sequence has not been provided e.g. NG_(NM_):c.PositionVariation')
                            continue
                        elif re.match('^NC_', input):
                            validation['warnings'] = validation[
                                                         'warnings'] + ': ' + 'A transcript reference sequence has not been provided e.g. NC_(NM_):c.PositionVariation. Unable to predict available transripts because chromosomal position is not specified'
                            logger.warning(
                                'A transcript reference sequence has not been provided e.g. NC_(NM_):c.PositionVariation. Unable to predict available transripts because chromosomal position is not specified')
                            continue
                        else:
                            pass
                    except:
                        exceptPass()

                logger.trace("Chromosomal/RefSeqGene reference catching complete", validation)
                # Find not_sub type in input e.g. GGGG>G
                """
                VCF2HGVS conversion step 4 has two purposes
                1. VCF is frequently inappropriately converted into HGVS like descriptions
                such as GGGG>G which is actually a delins, del or ins. The function assigns 
                the correct edit type
                2. Detects and extracts multiple ALT sequences into HGVS descriptions and
                automatically submits them for validation  
                """
                not_sub = copy.deepcopy(input)
                not_sub_find = re.compile("([GATCgatc]+)>([GATCgatc]+)")
                if not_sub_find.search(not_sub):
                    try:
                        # If the length of either side of the substitution delimer (>) is >1
                        matches = not_sub_find.search(not_sub)
                        if len(matches.group(1)) > 1 or len(matches.group(2)) > 1 or re.search(
                                "([GATCgatc]+)>([GATCgatc]+),([GATCgatc]+)", input):
                            # Search for and remove range
                            interval_range = re.compile("([0-9]+)_([0-9]+)")
                            if interval_range.search(not_sub):
                                m = not_sub_find.search(not_sub)
                                start = m.group(1)
                                delete = m.group(2)
                                beginning_string, middle_string = not_sub.split(':')
                                middle_string = middle_string.split('_')[0]
                                end_string = start + '>' + delete
                                not_sub = beginning_string + ':' + middle_string + end_string
                            # Split description
                            split_colon = not_sub.split(':')
                            ref_ac = split_colon[0]
                            remainder = split_colon[1]
                            split_dot = remainder.split('.')
                            ref_type = split_dot[0]
                            remainder = split_dot[1]
                            posedit = remainder
                            split_greater = remainder.split('>')
                            insert = split_greater[1]
                            remainder = split_greater[0]
                            # Split remainder using matches
                            r = re.compile("([0-9]+)([GATCgatc]+)")
                            try:
                                m = r.search(remainder)
                                start = m.group(1)
                                delete = m.group(2)
                                starts = posedit.split(delete)[0]
                                re_try = ref_ac + ':' + ref_type + '.' + starts + 'del' + delete[0] + 'ins' + insert
                                hgvs_re_try = hp.parse_hgvs_variant(re_try)
                                hgvs_re_try.posedit.edit.ref = delete
                                start_pos = str(hgvs_re_try.posedit.pos.start)
                                if re.search('\-', start_pos):
                                    base, offset = start_pos.split('-')
                                    new_offset = 0 - int(offset) + (len(delete))
                                    end_pos = int(base)
                                    hgvs_re_try.posedit.pos.end.base = int(end_pos)
                                    hgvs_re_try.posedit.pos.end.offset = int(new_offset) - 1
                                    not_delins = ref_ac + ':' + ref_type + '.' + start_pos + '_' + str(
                                        hgvs_re_try.posedit.pos.end) + 'del' + delete + 'ins' + insert
                                elif re.search('\+', start_pos):
                                    base, offset = start_pos.split('+')
                                    end_pos = int(base) + (len(delete) - int(offset) - 1)
                                    new_offset = 0 + int(offset) + (len(delete) - 1)
                                    hgvs_re_try.posedit.pos.end.base = int(end_pos)
                                    hgvs_re_try.posedit.pos.end.offset = int(new_offset)
                                    not_delins = ref_ac + ':' + ref_type + '.' + start_pos + '_' + str(
                                        hgvs_re_try.posedit.pos.end) + 'del' + delete + 'ins' + insert
                                else:
                                    end_pos = int(start_pos) + (len(delete) - 1)
                                    not_delins = ref_ac + ':' + ref_type + '.' + start_pos + '_' + str(
                                        end_pos) + 'del' + delete + 'ins' + insert
                            except:
                                exceptPass()
                                not_delins = not_sub
                            # Parse into hgvs object
                            try:
                                hgvs_not_delins = hp.parse_hgvs_variant(not_delins)
                            except hgvs.exceptions.HGVSError as e:
                                # Sort out multiple ALTS from VCF inputs
                                if re.search("([GATCgatc]+)>([GATCgatc]+),([GATCgatc]+)", not_delins):
                                    header, alts = not_delins.split('>')
                                    # Split up the alts into a list
                                    alt_list = alts.split(',')
                                    # Assemble and re-submit
                                    for alt in alt_list:
                                        validation[
                                            'warnings'] = 'Multiple ALT sequences detected: auto-submitting all possible combinations'
                                        validation['write'] = 'false'
                                        refreshed_description = header + '>' + alt
                                        query = {'quibble': refreshed_description, 'id': validation['id'],
                                                 'warnings': validation['warnings'], 'description': '', 'coding': '',
                                                 'coding_g': '', 'genomic_r': '', 'genomic_g': '', 'protein': '',
                                                 'write': 'true', 'primary_assembly': primary_assembly,
                                                 'order': ordering}
                                        batch_list.append(query)
                                        logger.resub(
                                            'Multiple ALT sequences detected. Auto-submitting all possible combinations.')
                                    continue
                                else:
                                    error = str(e)
                                    issue_link = ''
                                    validation['warnings'] = validation['warnings'] + ': ' + error
                                    logger.warning(str(e))
                                    continue

                            # Re-Stash the input as an HGVS
                            stash_input = copy.copy(hgvs_not_delins)
                            try:
                                not_delins = str(hn.normalize(hgvs_not_delins))
                            except hgvs.exceptions.HGVSError as e:
                                error = str(e)
                                if re.search('Normalization of intronic variants is not supported', error):
                                    not_delins = not_delins
                                else:
                                    issue_link = ''
                                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                    logger.warning(str(e))
                                    continue
                            # Create warning
                            caution = 'Variant description ' + input + ' is not HGVS compliant'
                            automap = input + ' automapped to ' + not_delins
                            validation['warnings'] = validation['warnings'] + ': ' + automap
                            # Change input to normalized variant
                            input = not_delins
                        else:
                            pass
                    except:
                        exceptPass()
                else:
                    pass
                logger.trace("Completed VCF-HVGS step 4", validation)

                # Tackle edit1234 type
                """
                Warns that descriptions such as c.ins12 or g.del69 are not HGVS compliant
                Strips the trailing numbers and tries to parse the description into an
                hgvs object.
                If parses, provides a warning including links to the VarNomen web page, but
                continues validation
                If not, an error message is generated and the loop continues
                """
                edit_pass = re.compile('_\d+$')
                edit_fail = re.compile('\d+$')
                if edit_fail.search(input):
                    if edit_pass.search(input):
                        pass
                    else:
                        error = 'false'
                        issue_link = 'false'
                        failed = copy.deepcopy(input)
                        # Catch the trailing digits
                        digits = re.search(r"(\d+$)", failed)
                        digits = digits.group(1)
                        remove = str(digits)  + 'end_anchor'
                        failed = failed + 'end_anchor'
                        failed = failed.replace(remove, '')

                        # Remove them so that the string SHOULD parse
                        try:
                            hgvs_failed = hp.parse_hgvs_variant(failed)
                        except hgvs.exceptions.HGVSError as e:
                            error = str(e)
                            error = 'The syntax of the input variant description is invalid '
                            if re.search('ins$', failed):
                                issue_link = 'http://varnomen.hgvs.org/recommendations/DNA/variant/insertion/'
                                error = error + ' please refer to ' + issue_link
                            validation['warnings'] = validation['warnings'] + error
                            logger.warning(str(error) + " " + str(e))
                            continue
                        hgvs_failed.posedit.edit = str(hgvs_failed.posedit.edit).replace(digits, '')
                        failed = str(hgvs_failed)
                        hgvs_failed = hp.parse_hgvs_variant(failed)
                        automap = 'Non HGVS compliant variant description ' + input + ' automapped to ' + failed
                        validation['warnings'] = validation['warnings'] + ': ' + automap
                        logger.warning(automap)
                        input = failed

                logger.trace("Ins/Del reference catching complete", validation)
                # Tackle compound variant descriptions NG or NC (NM_) i.e. correctly input NG/NC_(NM_):c.
                """
                Fully HGVS compliant intronic variant descriptions take the format e.g
                NG_007400.1(NM_000088.3):c.589-1G>T. However, hgvs cannot parse and map
                these variant strings. 
                This function:
                Removes the g. reference sequence
                NG_007400.1(NM_000088.3):c.589-1G>T ---> (NM_000088.3):c.589-1G>T
                Removes the parintheses
                (NM_000088.3):c.589-1G>T ---> NM_000088.3:c.589-1G>T
                hgvs can now parse the string into an hgvs variant object and manipulate it
                """
                caution = ''
                compounder = re.compile('\(NM_')
                compounder_b = re.compile('\(ENST')
                if compounder.search(input):
                    # Find pattern e.g. +0000 and assign to a variable
                    transy = re.search(r"(NM_.+)", input)
                    transy = transy.group(1)
                    transy = transy.replace(')', '')
                    input = transy
                logger.trace("HVGS typesetting complete", validation)
                # Extract variants from HGVS allele descriptions
                # http://varnomen.hgvs.org/recommendations/DNA/variant/alleles/
                """
                HGVS allele string parsing function Occurance #1
                Takes a single HGVS allele description and separates each allele into a 
                list of HGVS variants. The variants are then automatically submitted for
                validation. 
                Note: In this context, it is inappropriate to validate descriptions
                containing intronic variant descriptions. In such instances, allele
                descriptions should be re-submitted by the user at the gene or genome level 
                """
                if (re.search(':[gcnr].\[', input) and re.search('\;', input)) or (
                        re.search(':[gcrn].\d+\[', input) and re.search('\;', input)) or (re.search('\(\;\)', input)):
                    # handle LRG inputs
                    if re.match('^LRG', input):
                        if re.match('^LRG\d+', input):
                            string, remainder = input.split(':')
                            reference = string.replace('LRG', 'LRG_')
                            input = reference + ':' + remainder
                            caution = string + ' updated to ' + reference
                        if not re.match('^LRG_\d+', input):
                            pass
                        elif re.match('^LRG_\d+:g.', input) or re.match('^LRG_\d+:p.', input) or re.match('^LRG_\d+:c.',
                                                                                                          input) or re.match(
                            '^LRG_\d+:n.', input):
                            lrg_reference, variation = input.split(':')
                            refseqgene_reference = va_dbCrl.data.get_RefSeqGeneID_from_lrgID(lrg_reference)
                            if refseqgene_reference != 'none':
                                input = refseqgene_reference + ':' + variation
                                if caution == '':
                                    caution = lrg_reference + ':' + variation + ' automapped to ' + refseqgene_reference + ':' + variation
                                else:
                                    caution = caution + ': ' + lrg_reference + ':' + variation + ' automapped to ' + refseqgene_reference + ':' + variation
                                validation['warnings'] = validation['warnings'] + ': ' + str(caution)
                                logger.warning(str(caution))
                        elif re.match('^LRG_\d+t\d+:c.', input) or re.match('^LRG_\d+t\d+:n.', input) or re.match(
                                '^LRG_\d+t\d+:p.', input) or re.match('^LRG_\d+t\d+:g.', input):
                            lrg_reference, variation = input.split(':')
                            refseqtranscript_reference = va_dbCrl.data.get_RefSeqTranscriptID_from_lrgTranscriptID(
                                lrg_reference)
                            if refseqtranscript_reference != 'none':
                                input = refseqtranscript_reference + ':' + variation
                                if caution == '':
                                    caution = lrg_reference + ':' + variation + ' automapped to ' + refseqtranscript_reference + ':' + variation
                                else:
                                    caution = caution + ': ' + lrg_reference + ':' + variation + ' automapped to ' + refseqtranscript_reference + ':' + variation
                                validation['warnings'] = validation['warnings'] + ': ' + str(caution)
                                logger.warning(str(caution))
                        else:
                            pass
                    try:
                        # Submit to allele extraction function
                        alleles = va_func.hgvs_alleles(input, hp, vr, hn, vm, sf)
                        validation['warnings'] = validation[
                                                     'warnings'] + ': ' + 'Allele description syntax is valid' + ': Automap has extracted and validated all variant descriptions'
                        logger.resub('Automap has extracted possible variant descriptions, resubmitting')
                        for allele in alleles:
                            query = {'quibble': allele, 'id': validation['id'], 'warnings': validation['warnings'],
                                     'description': '', 'coding': '', 'coding_g': '', 'genomic_r': '', 'genomic_g': '',
                                     'protein': '', 'write': 'true', 'primary_assembly': primary_assembly,
                                     'order': ordering}
                            coding = 'intergenic'
                            batch_list.append(query)
                        validation['write'] = 'false'
                        continue
                    except va_func.alleleVariantError as e:
                        if re.search("Cannot validate sequence of an intronic variant", str(e)):
                            validation['warnings'] = validation[
                                                         'warnings'] + ': ' + 'Intronic positions not supported for HGVS Allele descriptions'
                            logger.warning('Intronic positions not supported for HGVS Allele descriptions')
                            continue
                        elif re.search('does not agree with reference sequence', str(e)):
                            er = str(e).replace(': ', '; ')
                            validation['warnings'] = validation[
                                                         'warnings'] + ': ' + er
                            continue
                        else:
                            raise variantValidatorError(str(e))
                logger.trace("HVGS String allele parsing pass 1 complete", validation)
                # INITIAL USER INPUT FORMATTING
                """
                Removes whitespace from the ends of the string
                Removes anything in brackets
                Identifies variant type
                Returns a dictionary containing the formatted input string and the variant type
                Accepts c, g, n, r currently
                """
                formatted = va_func.user_input(input)

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

                # Check the initial validity of the input
                if formatted == 'invalid':
                    if re.search('\w+\:[gcnmrp]', input) and not re.search('\w+\:[gcnmrp]\.', input):
                        error = 'Variant description ' + input + ' lacks the . character between <type> and <position> in the expected pattern <accession>:<type>.<position>'
                    else:
                        error = 'Variant description ' + input + ' is not in an accepted format'
                    validation['warnings'] = validation[
                                                 'warnings'] + ': ' + error
                    logger.warning(error)
                    continue
                else:
                    variant = formatted['variant']
                    input = formatted['variant']
                    stash_input = formatted['variant']
                    type = formatted['type']
                logger.trace("Variant input formatted, proceeding to validate.", validation)
                # Conversions
                """
                Conversions are not currently supported. The HGVS format for conversions
                is rarely seen wrt genomic sequencing data and needs to be re-evaluated
                """
                conversion = re.compile('con')
                if conversion.search(variant):
                    validation['warnings'] = validation['warnings'] + ': ' + 'Gene conversions currently unsupported'
                    logger.warning('Gene conversions currently unsupported')
                    continue

                # Primary check that hgvs will accept the variant
                error = 'false'
                # Change RNA bases to upper case but nothing else
                if type == ":r.":
                    variant = variant.upper()
                    variant = variant.replace(':R.', ':r.')
                    # lowercase the supported variant types
                    variant = variant.replace('DEL', 'del')
                    variant = variant.replace('INS', 'ins')
                    variant = variant.replace('INV', 'inv')
                    variant = variant.replace('DUP', 'dup')

                try:
                    input_parses = hp.parse_hgvs_variant(variant)
                except hgvs.exceptions.HGVSError as e:
                    error = str(e)
                if error == 'false':
                    if 'LRG' in input_parses.ac:
                        input_parses.ac.replace('T', 't')
                    else:
                        input_parses.ac = input_parses.ac.upper()
                    if hasattr(input_parses.posedit.edit, 'alt'):
                        if input_parses.posedit.edit.alt is not None:
                            input_parses.posedit.edit.alt = input_parses.posedit.edit.alt.upper()
                    if hasattr(input_parses.posedit.edit, 'ref'):
                        if input_parses.posedit.edit.ref is not None:
                            input_parses.posedit.edit.ref = input_parses.posedit.edit.ref.upper()
                    variant = str(input_parses)
                    input = str(input_parses)
                    pass
                else:
                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                    logger.warning(error)
                    continue

                """
                ENST support needs to be re-evaluated, but is very low priority
                ENST not supported by ACMG and is under review by HGVS
                """
                if re.match('^ENST', str(input_parses)):
                    trap_ens_in = str(input_parses)
                    sim_tx = hdp.get_similar_transcripts(input_parses.ac)
                    for line in sim_tx:
                        if str(line[2]) == 'True' and str(line[3]) == 'True' and str(line[4]) == 'True' and str(
                                line[5]) == 'True' and str(line[6]) == 'True':
                            input_parses.ac = (line[1])
                            input = str(input_parses)
                            variant = input
                            break
                    if re.match('^ENST', str(input_parses)):
                        error = 'Unable to map ' + str(input_parses.ac) + ' to an equivalent RefSeq transcript'
                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                        logger.warning(str(error))
                        continue
                    else:
                        validation['warnings'] = validation['warnings'] + ': ' + str(
                            trap_ens_in) + ' automapped to equivalent RefSeq transcript ' + variant
                        logger.warning(str(trap_ens_in) + ' automapped to equivalent RefSeq transcript ' + variant)
                logger.trace("HVGS acceptance test passed", validation)
                # Check whether supported genome build is requested for non g. descriptions
                historic_assembly = 'false'
                mapable_assemblies = {
                    'GRCh37': 'true',
                    'GRCh38': 'true',
                    'NCBI36': 'false'
                }
                is_mapable = mapable_assemblies.get(primary_assembly)
                if is_mapable == 'true':

                    # These objects cannot be moved outside of the main function because they gather data from the
                    # iuser input e.g. alignment method and genome build
                    # They initiate quickly, so no need to move them unnecessarily

                    # Create easy variant mapper (over variant mapper) and splign locked evm
                    evm = hgvs.assemblymapper.AssemblyMapper(hdp,
                                                             assembly_name=primary_assembly,
                                                             alt_aln_method=alt_aln_method,
                                                             normalize=True,
                                                             replace_reference=True
                                                             )

                    # Setup a reverse normalize instance and non-normalize evm
                    no_norm_evm = hgvs.assemblymapper.AssemblyMapper(hdp,
                                                                     assembly_name=primary_assembly,
                                                                     alt_aln_method=alt_aln_method,
                                                                     normalize=False,
                                                                     replace_reference=True
                                                                     )

                    # Create a specific minimal evm with no normalizer and no replace_reference
                    min_evm = hgvs.assemblymapper.AssemblyMapper(hdp,
                                                                 assembly_name=primary_assembly,
                                                                 alt_aln_method=alt_aln_method,
                                                                 normalize=False,
                                                                 replace_reference=False
                                                                 )

                else:
                    error = 'Mapping of ' + variant + ' to genome assembly ' + primary_assembly + ' is not supported'
                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                    logger.warning(str(error))
                    continue
                # Catch interval end > interval start
                """
                hgvs did/does not handle 3' UTR position ordering well. This function
                ensures that end pos is not > start pos wrt 3' UTRs. 
                Also identifies some variants which span into the downstream sequence
                i.e. out of bounds
                """
                astr = re.compile('\*')
                if astr.search(str(input_parses.posedit)):
                    input_parses_copy = copy.deepcopy(input_parses)
                    input_parses_copy.type = "c"
                    # Map to n. position
                    # Create easy variant mapper (over variant mapper) and splign locked evm
                    try:
                        to_n = evm.c_to_n(input_parses_copy)
                    except hgvs.exceptions.HGVSError as e:
                        exceptPass()
                    else:
                        if to_n.posedit.pos.end.base < to_n.posedit.pos.start.base:
                            error = 'Interval end position < interval start position '
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            logger.warning(str(error))
                            continue
                elif input_parses.posedit.pos.end.base < input_parses.posedit.pos.start.base:
                    error = 'Interval end position ' + str(
                        input_parses.posedit.pos.end.base) + ' < interval start position ' + str(
                        input_parses.posedit.pos.start.base)
                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                    logger.warning(str(error))
                    continue
                else:
                    pass

                # Catch missing version number in refseq
                ref_type = re.compile("^N\w\w\d")
                is_version = re.compile("\d\.\d")
                en_type = re.compile('^ENS')
                lrg_type = re.compile('LRG')
                if (ref_type.search(str(input_parses)) and is_version.search(str(input_parses))) or (
                        en_type.search(str(input_parses))):
                    pass
                else:
                    if lrg_type.search(str(input_parses)):
                        pass
                    if ref_type.search(str(input_parses)):
                        error = 'RefSeq variant accession numbers MUST include a version number'
                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                        continue
                logger.trace("HVGS interval/version mapping complete", validation)

                # handle LRG inputs
                """
                LRG and LRG_t reference sequence identifiers need to be replaced with 
                equivalent RefSeq identifiers. The lookup data is stored in the 
                VariantValidator  MySQL database
                """
                if re.match('^LRG', str(input_parses)):
                    if re.match('^LRG\d+', str(input_parses.ac)):
                        string = str(input_parses.ac)
                        reference = string.replace('LRG', 'LRG_')
                        input_parses.ac = reference

                        caution = string + ' updated to ' + reference
                    if not re.match('^LRG_\d+', str(input_parses)):
                        pass
                    elif re.match('^LRG_\d+:g.', str(input_parses)) or re.match('^LRG_\d+:p.',
                                                                                str(input_parses)) or re.match(
                        '^LRG_\d+:c.', str(input_parses)) or re.match('^LRG_\d+:n.', str(input_parses)):
                        lrg_reference, variation = str(input_parses).split(':')
                        refseqgene_reference = va_dbCrl.data.get_RefSeqGeneID_from_lrgID(lrg_reference)
                        if refseqgene_reference != 'none':
                            input_parses.ac = refseqgene_reference
                            variant = str(input_parses)
                            input = str(input_parses)
                            stash_input = input
                            if caution == '':
                                caution = lrg_reference + ':' + variation + ' automapped to ' + refseqgene_reference + ':' + variation
                            else:
                                caution = caution + ': ' + lrg_reference + ':' + variation + ' automapped to ' + refseqgene_reference + ':' + variation
                            validation['warnings'] = validation['warnings'] + ': ' + str(caution)
                            logger.warning(str(caution))
                    elif re.match('^LRG_\d+t\d+:c.', str(input_parses)) or re.match('^LRG_\d+t\d+:n.',
                                                                                    str(input_parses)) or re.match(
                        '^LRG_\d+t\d+:p.', str(input_parses)) or re.match('^LRG_\d+t\d+:g.', str(input_parses)):
                        lrg_reference, variation = str(input_parses).split(':')
                        refseqtranscript_reference = va_dbCrl.data.get_RefSeqTranscriptID_from_lrgTranscriptID(
                            lrg_reference)
                        if refseqtranscript_reference != 'none':
                            input_parses.ac = refseqtranscript_reference
                            variant = str(input_parses)
                            input = str(input_parses)
                            stash_input = input
                            if caution == '':
                                caution = lrg_reference + ':' + variation + ' automapped to ' + refseqtranscript_reference + ':' + variation
                            else:
                                caution = caution + ': ' + lrg_reference + ':' + variation + ' automapped to ' + refseqtranscript_reference + ':' + variation
                            validation['warnings'] = validation['warnings'] + ': ' + str(caution)
                            logger.warning(str(caution))
                    else:
                        pass
                logger.trace("LRG check for conversion to refseq completed", validation)
                # Additional Incorrectly input variant capture training
                """
                Evolving list of common mistakes, see sections below
                """
                # NM_ .g
                if (re.search('^NM_', variant) or re.search('^NR_', variant)) and re.search(':g.', variant):
                    suggestion = input.replace(':g.', ':c.')
                    error = 'Transcript reference sequence input as genomic (g.) reference sequence. Did you mean ' + suggestion + '?'
                    validation['warnings'] = validation['warnings'] + ': ' + error
                    logger.warning(error)
                    continue
                # NR_ c.
                if re.search('^NR_', input) and re.search(':c.', input):
                    suggestion = input.replace(':c.', ':n.')
                    error = 'Non-coding transcript reference sequence input as coding (c.) reference sequence. Did you mean ' + suggestion + '?'
                    validation['warnings'] = validation['warnings'] + ': ' + error
                    logger.warning(error)
                    continue
                # NM_ n.
                if re.search('^NM_', input) and re.search(':n.', input):
                    suggestion = input.replace(':n.', ':c.')
                    error = 'Coding transcript reference sequence input as non-coding transcript (n.) reference sequence. Did you mean ' + suggestion + '?'
                    validation['warnings'] = validation['warnings'] + ': ' + error
                    logger.warning(error)
                    continue

                # NM_ NC_ NG_ NR_ p.
                if (re.search('^NM_', variant) or re.search('^NR_', variant) or re.search('^NC_', variant) or re.search(
                        '^NG_', variant)) and re.search(':p.', variant):
                    issue_link = 'http://varnomen.hgvs.org/recommendations/protein/'
                    error = 'Using a nucleotide reference sequence (NM_ NR_ NG_ NC_) to specify protein-level (p.) variation is not HGVS compliant. Please select an appropriate protein reference sequence (NP_)'
                    validation['warnings'] = validation['warnings'] + ': ' + error
                    logger.warning(error)
                    continue

                # NG_ c or NC_c..
                if (re.search('^NG_', variant) or re.search('^NC_', variant)) and re.search(':c.', variant):
                    suggestion = ': For additional assistance, submit ' + str(variant) + ' to VariantValidator'
                    error = 'NG_:c.PositionVariation descriptions should not be used unless a transcript reference sequence has also been provided e.g. NG_(NM_):c.PositionVariation' + suggestion
                    validation['warnings'] = validation['warnings'] + ': ' + error
                    logger.warning(error)
                    continue

                logger.trace("Passed 'common mistakes' catcher", validation)
                # Primary validation of the input
                """
                An evolving set of variant structure and content searches which identify 
                and warn users about inappropriate use of HGVS
                Primarily, this code filters out variants that cannot realistically be 
                auto corrected and will cause the downstream functions to return errors
                """
                input_parses = hp.parse_hgvs_variant(input)
                if input_parses.type == 'g':
                    if re.match('^NC_', input_parses.ac) or re.match('^NG_', input_parses.ac) or re.match('^NT_',
                                                                                                          input_parses.ac) or re.match(
                        '^NW_', input_parses.ac):
                        pass
                    else:
                        error = 'Invalid reference sequence identifier (' + input_parses.ac + ')'
                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                        logger.warning(error)
                        continue
                    try:
                        vr.validate(input_parses)
                    except hgvs.exceptions.HGVSError as e:
                        error = str(e)
                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                        logger.warning(error)
                        continue
                    except Exception as e:
                        error = str(e)
                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                        logger.warning(error)
                        continue
                    # Additional test
                    try:
                        hn.normalize(input_parses)
                    except hgvs.exceptions.HGVSError as e:
                        error = str(e)
                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                        logger.warning(error)
                        continue
                    else:
                        exceptPass()

                elif input_parses.type == 'c':
                    if re.search('\*', str(input_parses)) or re.search('c.\-', str(input_parses)):
                        # Catch variation in UTRs
                        # These should be in the sequence so can be directly validated. Need to pass to n.
                        try:
                            vr.validate(input_parses)
                        except hgvs.exceptions.HGVSError as e:
                            error = str(e)
                            if re.search('datums is ill-defined', error):
                                called_ref = input_parses.posedit.edit.ref
                                try:
                                    to_n = evm.c_to_n(input_parses)
                                except hgvs.exceptions.HGVSInvalidVariantError as e:
                                    error = str(e)
                                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                    logger.warning(error)
                                    continue
                                actual_ref = to_n.posedit.edit.ref
                                if called_ref != actual_ref:
                                    error = 'Variant reference (' + called_ref + ') does not agree with reference sequence (' + actual_ref + ')'
                                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                    logger.warning(error)
                                    continue
                                else:
                                    input_parses.posedit.edit.ref = ''
                                    variant = str(input_parses)
                            else:
                                if re.search('bounds', error) or re.search('intronic variant', error):
                                    try:
                                        hn.normalize(input_parses)
                                    except hgvs.exceptions.HGVSError as e:
                                        exceptPass()
                                    if re.search('bounds', str(e)):
                                        try:
                                            identity_info = hdp.get_tx_identity_info(input_parses.ac)
                                            ref_start = identity_info[3]
                                            ref_end = identity_info[4]
                                            if re.match('-', str(
                                                    input_parses.posedit.pos.start)) and input_parses.posedit.pos.start.offset == 0:
                                                # upstream positions
                                                boundary = int('-' + str(ref_start))
                                                remainder = int(str(input_parses.posedit.pos.start)) - boundary
                                                input_parses.posedit.pos.start.base = boundary
                                                input_parses.posedit.pos.start.offset = remainder
                                            if re.match('-', str(
                                                    input_parses.posedit.pos.end)) and input_parses.posedit.pos.end.offset == 0:
                                                boundary = int('-' + str(ref_start))
                                                remainder = int(str(input_parses.posedit.pos.end)) - boundary
                                                input_parses.posedit.pos.end.base = boundary
                                                input_parses.posedit.pos.end.offset = remainder
                                            if re.match('\*', str(
                                                    input_parses.posedit.pos.start)) and input_parses.posedit.pos.start.offset == 0:
                                                # downstream positions
                                                tot_end_pos = str(input_parses.posedit.pos.start).replace('*', '')
                                                ts_seq = sf.fetch_seq(input_parses.ac)
                                                boundary = len(ts_seq) - ref_end
                                                input_parses.posedit.pos.start.base = boundary
                                                offset = int(tot_end_pos) - int(boundary)
                                                input_parses.posedit.pos.start.offset = offset
                                            if re.match('\*', str(
                                                    input_parses.posedit.pos.end)) and input_parses.posedit.pos.end.offset == 0:
                                                tot_end_pos = str(input_parses.posedit.pos.end).replace('*', '')
                                                ts_seq = sf.fetch_seq(input_parses.ac)
                                                boundary = len(ts_seq) - ref_end
                                                input_parses.posedit.pos.end.base = boundary
                                                offset = int(tot_end_pos) - int(boundary)
                                                input_parses.posedit.pos.end.offset = offset

                                            # Create a lose vm instance
                                            lose_vm = hgvs.variantmapper.VariantMapper(hdp,
                                                                                       replace_reference=True,
                                                                                       prevalidation_level=None
                                                                                       )


                                            report_gen = va_func.myevm_t_to_g(input_parses, hdp, no_norm_evm,
                                                                              primary_assembly, lose_vm, hp, hn, sf, nr_vm)
                                            error = 'Using a transcript reference sequence to specify a variant position that lies outside of the reference sequence is not HGVS-compliant: Instead use ' + valstr(
                                                report_gen)
                                        except Exception as e:
                                            exceptPass()
                                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                        logger.warning(str(error))
                                        continue
                                    else:
                                        pass
                                else:
                                    pass

                        try:
                            input_parses = evm.c_to_n(input_parses)
                        except hgvs.exceptions.HGVSError as e:
                            error = str(e)
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            logger.warning(str(e))
                            continue

                        if re.search('n.1-', str(input_parses)):
                            input_parses = evm.n_to_c(input_parses)
                            error = 'Using a transcript reference sequence to specify a variant position that lies outside of the reference sequence is not HGVS-compliant. Instead use '
                            genomic_position = va_func.myevm_t_to_g(input_parses, hdp, no_norm_evm, primary_assembly,
                                                                    vm, hp, hn, sf, nr_vm)
                            error = error + valstr(genomic_position)
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            logger.warning(str(error))
                            continue
                        else:
                            pass

                        # Re-map input_parses back to c. variant
                        input_parses = evm.n_to_c(input_parses)

                        # Intronic positions in UTRs
                        if re.search('\d\-\d', str(input_parses)) or re.search('\d\+\d', str(input_parses)):
                            # Can we go c-g-c
                            try:
                                to_genome = va_func.myevm_t_to_g(input_parses, hdp, no_norm_evm, primary_assembly, vm,
                                                                 hp, hn, sf, nr_vm)
                                to_tx = evm.g_to_t(to_genome, input_parses.ac)
                            except hgvs.exceptions.HGVSInvalidIntervalError as e:
                                error = str(e)
                                if re.search('bounds', error):
                                    try:
                                        identity_info = hdp.get_tx_identity_info(input_parses.ac)
                                        ref_start = identity_info[3]
                                        ref_end = identity_info[4]
                                        if re.match('-', str(input_parses.posedit.pos.start)):
                                            # upstream positions
                                            boundary = int('-' + str(ref_start))
                                            remainder = int(str(input_parses.posedit.pos.start)) - boundary
                                            input_parses.posedit.pos.start.base = boundary
                                            input_parses.posedit.pos.start.offset = remainder
                                        if re.match('-', str(input_parses.posedit.pos.end)):
                                            boundary = int('-' + str(ref_start))
                                            remainder = int(str(input_parses.posedit.pos.end)) - boundary
                                            input_parses.posedit.pos.end.base = boundary
                                            input_parses.posedit.pos.end.offset = remainder
                                        if re.match('\*', str(input_parses.posedit.pos.start)):
                                            # downstream positions
                                            tot_end_pos = str(input_parses.posedit.pos.start).replace('*', '')
                                            ts_seq = sf.fetch_seq(input_parses.ac)
                                            boundary = len(ts_seq) - ref_end
                                            input_parses.posedit.pos.start.base = boundary
                                            te1, te2 = tot_end_pos.split('+')
                                            tot_end_pos = int(te1) + int(te2)
                                            offset = int(tot_end_pos) - int(boundary)
                                            input_parses.posedit.pos.start.offset = offset
                                        if re.match('\*', str(input_parses.posedit.pos.end)):
                                            tot_end_pos = str(input_parses.posedit.pos.end).replace('*', '')
                                            ts_seq = sf.fetch_seq(input_parses.ac)
                                            boundary = len(ts_seq) - ref_end
                                            input_parses.posedit.pos.end.base = boundary
                                            te1, te2 = tot_end_pos.split('+')
                                            tot_end_pos = int(te1) + int(te2)
                                            offset = int(tot_end_pos) - int(boundary)
                                            input_parses.posedit.pos.end.offset = offset

                                        report_gen = va_func.myevm_t_to_g(input_parses, hdp, no_norm_evm,
                                                                          primary_assembly, lose_vm, hp, hn, sf, nr_vm)
                                        error = 'Using a transcript reference sequence to specify a variant position that lies outside of the reference sequence is not HGVS-compliant. Instead use ' + valstr(
                                            report_gen)
                                    except Exception as e:
                                        exceptPass()
                                else:
                                    pass
                                validation['warnings'] = validation['warnings'] + ': ' + str(
                                    error)
                                logger.warning(str(error))
                                continue

                            except hgvs.exceptions.HGVSDataNotAvailableError as e:
                                error = str(e)
                                if 'Alignment is incomplete' in error:
                                    e_list = error.split('~')
                                    gens = []
                                    for el in e_list:
                                        el_l = el.split('/')
                                        if el_l[-1] == '':
                                            continue
                                        gens.append(el_l[-1])
                                    acs = '; '.join(gens)
                                    error = 'Cannot map ' + valstr(
                                        input_parses) + ' to a genomic position. ' + input_parses.ac + ' can only be partially aligned to genomic reference sequences ' + acs
                                    validation['warnings'] = validation['warnings'] + ': ' + str(
                                        error)
                                    logger.warning(str(error))
                                    continue

                    elif re.search('\d\-', str(input_parses)) or re.search('\d\+', str(input_parses)):
                        # Quick look at syntax validation
                        try:
                            vr.validate(input_parses)
                        except hgvs.exceptions.HGVSInvalidVariantError as e:
                            error = str(e)
                            if re.search('bounds', error):
                                try:
                                    report_gen = va_func.myevm_t_to_g(input_parses, hdp, no_norm_evm, primary_assembly,
                                                                      lose_vm, hp, hn, sf, nr_vm)
                                except hgvs.exceptions.HGVSError as e:
                                    exceptPass()
                                else:
                                    error = 'Using a transcript reference sequence to specify a variant position that lies outside of the reference sequence is not HGVS-compliant. Instead use ' + valstr(
                                        report_gen)
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                            elif re.search('insertion length must be 1', error):
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                            elif re.search('base start position must be <= end position', error):
                                correction = copy.deepcopy(input_parses)
                                st = input_parses.posedit.pos.start
                                ed = input_parses.posedit.pos.end
                                correction.posedit.pos.start = ed
                                correction.posedit.pos.end = st
                                error = error + ': Did you mean ' + str(correction) + '?'
                                # error = 'Interval start position ' + str(input_parses.posedit.pos.start) + ' > interval end position ' + str(input_parses.posedit.pos.end)
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue

                        # Create a specific minimal evm with no normalizer and no replace_reference
                        # Have to use this method due to potential multi chromosome error, note, normalizes but does not replace sequence
                        try:
                            output = va_func.noreplace_myevm_t_to_g(input_parses, evm, hdp, primary_assembly, vm, hn,
                                                                    hp, sf, no_norm_evm)
                        except hgvs.exceptions.HGVSDataNotAvailableError as e:
                            tx_ac = input_parses.ac
                            try:
                                gene_symbol = va_dbCrl.data.get_gene_symbol_from_transcriptID(tx_ac)
                            except:
                                gene_symbol = 'none'
                            if gene_symbol == 'none':
                                error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
                            else:
                                error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' or ' + gene_symbol + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            logger.warning(str(error))
                            continue
                        except ValueError as e:
                            error = str(e)
                            if re.search('> end', error):
                                error = 'Interval start position ' + str(
                                    input_parses.posedit.pos.start) + ' > interval end position ' + str(
                                    input_parses.posedit.pos.end)
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                        except hgvs.exceptions.HGVSInvalidVariantError as e:
                            error = str(e)
                            if re.search('base start position must be <= end position', error):
                                correction = copy.deepcopy(input_parses)
                                st = input_parses.posedit.pos.start
                                ed = input_parses.posedit.pos.end
                                correction.posedit.pos.start = ed
                                correction.posedit.pos.end = st
                                error = error + ': Did you mean ' + str(correction) + '?'
                                error = 'Interval start position ' + str(
                                    input_parses.posedit.pos.start) + ' > interval end position ' + str(
                                    input_parses.posedit.pos.end)
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                            else:
                                error = str(e)
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue

                        try:
                            evm.g_to_t(output, input_parses.ac)
                        except hgvs.exceptions.HGVSError as e:
                            error = str(e)
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            logger.warning(str(error))
                            continue

                        try:
                            vr.validate(output)
                        except hgvs.exceptions.HGVSError as e:
                            error = str(e)
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            logger.warning(str(error))
                            continue

                    else:
                        # All other variation
                        try:
                            vr.validate(input_parses)
                        except hgvs.exceptions.HGVSUnsupportedOperationError:
                            exceptPass()
                        except hgvs.exceptions.HGVSInvalidVariantError as e:
                            error = str(e)
                            """
                            #Phil: Honestly not sure what the purpose of any of these is, we act the same regardless of what
                            #kind of error it is.
                            if re.search('Length implied by coordinates', error):
                                # Applies to del and inv
                                # NOTE, there has been no normalization at all so this error is valid here
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            # Will apply to > del and inv
                            if re.search('does not agree with reference sequence', error):
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            # ensures x_y for insertions
                            if re.search('insertion length must be 1', error):
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            # Boundary issue
                            if re.search('Variant coordinate is out of the bound of CDS region', error):
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            """
                            # This catches errors in introns
                            if re.search('base start position must be <= end position', error):
                                correction = copy.deepcopy(input_parses)
                                st = input_parses.posedit.pos.start
                                ed = input_parses.posedit.pos.end
                                correction.posedit.pos.start = ed
                                correction.posedit.pos.end = st
                                error = error + ': Did you mean ' + str(correction) + '?'
                                error = 'Interval start position ' + str(
                                    input_parses.posedit.pos.start) + ' > interval end position ' + str(
                                    input_parses.posedit.pos.end)
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            logger.warning(str(error))
                            continue

                        except hgvs.exceptions.HGVSDataNotAvailableError as e:
                            error = e
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            logger.warning(str(error))
                            continue
                        except hgvs.exceptions.HGVSError as e:
                            error = str(e)
                            if re.search('bounds', error):
                                error = error + ' (' + input_parses.ac + ')'
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                            exceptPass()


                elif input_parses.type == 'n':
                    if re.search('\+', str(input_parses)) or re.search('\-', str(input_parses)):
                        # Catch variation in UTRs
                        # These should be in the sequence so can be directly validated. Need to pass to n.
                        try:
                            vr.validate(input_parses)
                        except hgvs.exceptions.HGVSError as e:
                            error = str(e)
                            if re.search('intronic variant', error):
                                pass
                            elif re.search('datums is ill-defined', error):
                                called_ref = input_parses.posedit.edit.ref
                                to_n = evm.c_to_n(input_parses)
                                actual_ref = to_n.posedit.edit.ref
                                if called_ref != actual_ref:
                                    error = 'Variant reference (' + called_ref + ') does not agree with reference sequence (' + actual_ref + ')'
                                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                    logger.warning(str(error))
                                    continue
                                else:
                                    input_parses.posedit.edit.ref = ''
                                    variant = str(input_parses)

                            elif re.search('base must be >=1 for datum = SEQ_START or CDS_END', error):
                                error = 'The given coordinate is outside the bounds of the reference sequence.'

                                try:
                                    if re.match('-', str(input_parses.posedit.pos.start)):
                                        # upstream positions
                                        boundary = 1
                                        remainder = int(str(input_parses.posedit.pos.start)) - boundary
                                        remainder = remainder + 1
                                        input_parses.posedit.pos.start.base = boundary
                                        input_parses.posedit.pos.start.offset = remainder
                                    if re.match('-', str(input_parses.posedit.pos.end)):
                                        boundary = 1
                                        remainder = int(str(input_parses.posedit.pos.end)) - boundary
                                        remainder = remainder + 1
                                        input_parses.posedit.pos.end.base = boundary
                                        input_parses.posedit.pos.end.offset = remainder
                                    report_gen = va_func.myevm_t_to_g(input_parses, hdp, no_norm_evm, primary_assembly,
                                                                      lose_vm, hp, hn, sf, nr_vm)
                                    error = 'Using a transcript reference sequence to specify a variant position that lies outside of the reference sequence is not HGVS-compliant. Instead use ' + valstr(
                                        report_gen)
                                except Exception as e:
                                    exceptPass()
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                            else:
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue

                    if re.search('n.1-', str(input_parses)):
                        error = 'Using a transcript reference sequence to specify a variant position that lies outside of the reference sequence is not HGVS-compliant. Instead use '
                        genomic_position = va_func.myevm_t_to_g(input_parses, hdp, no_norm_evm, primary_assembly, vm,
                                                                hp, hn, sf, nr_vm)
                        error = error + valstr(genomic_position)
                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                        logger.warning(str(error))
                        continue
                    else:
                        pass

                    if re.search('\d\-', str(input_parses)) or re.search('\d\+', str(input_parses)):
                        # Quick look at syntax validation
                        try:
                            vr.validate(input_parses)
                        except hgvs.exceptions.HGVSInvalidVariantError as e:
                            error = str(e)
                            if re.search('bounds', error):
                                try:
                                    report_gen = va_func.myevm_t_to_g(input_parses, hdp, no_norm_evm, primary_assembly,
                                                                      lose_vm, hp, hn, sf, nr_vm)
                                except hgvs.exceptions.HGVSError as e:
                                    exceptPass()
                                else:
                                    error = 'Using a transcript reference sequence to specify a variant position that lies outside of the reference sequence is not HGVS-compliant. Instead use ' + valstr(
                                        report_gen)
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                            elif re.search('insertion length must be 1', error):
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                            elif re.search('base start position must be <= end position', error):
                                correction = copy.deepcopy(input_parses)
                                st = input_parses.posedit.pos.start
                                ed = input_parses.posedit.pos.end
                                correction.posedit.pos.start = ed
                                correction.posedit.pos.end = st
                                error = error + ': Did you mean ' + str(correction) + '?'
                                # error = 'Interval start position ' + str(input_parses.posedit.pos.start) + ' > interval end position ' + str(input_parses.posedit.pos.end)
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                            elif re.search('Cannot validate sequence of an intronic variant', error):
                                try:
                                    test_g = va_func.myevm_t_to_g(input_parses, hdp, no_norm_evm, primary_assembly, vm,
                                                                  hp, hn, sf, nr_vm)
                                    back_to_n = evm.g_to_t(test_g, input_parses.ac)
                                except hgvs.exceptions.HGVSError as e:
                                    error = str(e)
                                    if re.search('bounds', error):
                                        report_gen = va_func.myevm_t_to_g(input_parses, hdp, no_norm_evm,
                                                                          primary_assembly, lose_vm, hp, hn, sf, nr_vm)
                                        error = 'Using a transcript reference sequence to specify a variant position that lies outside of the reference sequence is not HGVS-compliant. Instead use ' + valstr(
                                            report_gen)
                                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                        logger.warning(str(error))
                                        continue
                                else:
                                    exceptPass()

                        # Create a specific minimal evm with no normalizer and no replace_reference
                        # Have to use this method due to potential multi chromosome error, note, normalizes but does not replace sequence
                        try:
                            output = va_func.noreplace_myevm_t_to_g(input_parses, evm, hdp, primary_assembly, vm, hn,
                                                                    hp, sf, no_norm_evm)
                        except hgvs.exceptions.HGVSDataNotAvailableError as e:
                            tx_ac = input_parses.ac
                            try:
                                gene_symbol = va_dbCrl.data.get_gene_symbol_from_transcriptID(tx_ac)
                            except:
                                gene_symbol = 'none'
                            if gene_symbol == 'none':
                                error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
                            else:
                                error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' or ' + gene_symbol + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            logger.warning(str(error))
                            continue
                        except ValueError as e:
                            error = str(e)
                            if re.search('> end', error):
                                error = 'Interval start position ' + str(
                                    input_parses.posedit.pos.start) + ' > interval end position ' + str(
                                    input_parses.posedit.pos.end)
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                        except hgvs.exceptions.HGVSInvalidVariantError as e:
                            error = str(e)
                            if re.search('base start position must be <= end position', error):
                                correction = copy.deepcopy(input_parses)
                                st = input_parses.posedit.pos.start
                                ed = input_parses.posedit.pos.end
                                correction.posedit.pos.start = ed
                                correction.posedit.pos.end = st
                                error = error + ': Did you mean ' + str(correction) + '?'
                                error = 'Interval start position ' + str(
                                    input_parses.posedit.pos.start) + ' > interval end position ' + str(
                                    input_parses.posedit.pos.end)
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                        try:
                            vr.validate(output)
                        except hgvs.exceptions.HGVSError as e:
                            error = str(e)
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            logger.warning(str(error))
                            continue

                    else:
                        # All other variation
                        try:
                            vr.validate(input_parses)
                        except hgvs.exceptions.HGVSUnsupportedOperationError:

                            exceptPass()
                        except hgvs.exceptions.HGVSInvalidVariantError as e:
                            error = str(e)
                            """
                            if re.search('Length implied by coordinates', error):
                                # Applies to del and inv
                                # NOTE, there has been no normalization at all so this error is valid here
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                continue
                            # Will apply to > del and inv
                            if re.search('does not agree with reference sequence', error):
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                continue
                            # ensures x_y for insertions
                            if re.search('insertion length must be 1', error):
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                continue
                            # Boundary issue
                            if re.search('Variant coordinate is out of the bound of CDS region', error):
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                continue
                            """
                            # This catches errors in introns
                            if re.search('base start position must be <= end position', error):
                                correction = copy.deepcopy(input_parses)
                                st = input_parses.posedit.pos.start
                                ed = input_parses.posedit.pos.end
                                correction.posedit.pos.start = ed
                                correction.posedit.pos.end = st
                                error = error + ': Did you mean ' + str(correction) + '?'
                                error = 'Interval start position ' + str(
                                    input_parses.posedit.pos.start) + ' > interval end position ' + str(
                                    input_parses.posedit.pos.end)
                                logger.warning(str(error))
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                continue
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            logger.warning(str(error))
                            continue
                        except hgvs.exceptions.HGVSDataNotAvailableError as e:
                            error = e
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            logger.warning(str(error))
                            continue
                        except hgvs.exceptions.HGVSError as e:
                            error = str(e)
                            if re.search('bounds', error):
                                error = error + ' (' + input_parses.ac + ')'
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                else:
                    pass
                logger.trace("Variant structure and contents searches passed", validation)
                # Mitochondrial variants
                """
                Reformat m. into the new HGVS standard which is now m again!
                """
                if type == ':m.' or re.match('NC_012920.1', str(input_parses.ac)) or re.match('NC_001807.4',
                                                                                              str(input_parses.ac)):
                    hgvs_mito = copy.deepcopy(input_parses)
                    if (re.match('NC_012920.1', str(hgvs_mito.ac)) and hgvs_mito.type == 'g') or (
                            re.match('NC_001807.4', str(hgvs_mito.ac)) and hgvs_mito.type == 'g'):
                        hgvs_mito.type = 'm'
                        caution = ''
                    try:
                        vr.validate(hgvs_mito)
                    except hgvs.exceptions.HGVSError as e:
                        error = caution + ': ' + str(e)
                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                        logger.warning(str(error))
                        continue
                    except KeyError as e:
                        error = caution + ': Currently unable to validate ' + hgvs_mito.ac + ' sequence variation'
                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                        logger.warning(str(error))
                        continue
                    else:
                        # Any transcripts?
                        rel_var = va_func.relevant_transcripts(hgvs_mito, evm, hdp, alt_aln_method, reverse_normalizer, hp)
                        hgvs_genomic = copy.deepcopy(hgvs_mito)
                        if len(rel_var) == 0:
                            validation['genomic_g'] = valstr(hgvs_mito)
                            validation['description'] = 'Homo sapiens mitochondrion, complete genome'
                            logger.info('Homo sapiens mitochondrion, complete genome')
                            continue
                        # Currently we are not expecting this path to be activated because not m. transcripts seem to be NM_
                        # This route may throw up errors in the future
                        else:
                            pass

                # handle :p.
                if type == ':p.':
                    error = 'false'
                    # Try to validate the variant
                    try:
                        hgvs_object = hp.parse_hgvs_variant(variant)
                    except hgvs.exceptions.HGVSError as e:
                        error = str(e)
                    try:
                        vr.validate(hgvs_object)
                    except hgvs.exceptions.HGVSError as e:
                        error = str(e)
                    if error != 'false':
                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                        logger.warning(str(error))
                        continue
                    else:
                        # Get accurate descriptions from the relevant databases
                        # RefSeq databases
                        if alt_aln_method != 'genebuild':
                            # Gene description  - requires GenBank search to get all the required info, i.e. transcript variant ID
                            # accession number
                            hgvs_object = hp.parse_hgvs_variant(variant)
                            accession = hgvs_object.ac
                            # Look for the accession in our database
                            # Connect to database and send request
                            record = va_func.entrez_efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
                            try:
                                description = record.description
                            except:
                                description = 'Unable to recover the description of ' + accession + ' from Entrez'
                            try:
                                vr.validate(hgvs_object)
                            except hgvs.exceptions.HGVSError as e:
                                error = str(e)
                            else:
                                error = str(
                                    hgvs_object) + ' is HGVS compliant and contains a valid reference amino acid description'
                            reason = 'Protein level variant descriptions are not fully supported due to redundancy in the genetic code'
                            validation['warnings'] = validation['warnings'] + ': ' + str(reason) + ': ' + str(error)
                            validation['protein'] = str(hgvs_object)
                            logger.warning(str(reason) + ": " + str(error))
                            continue

                # handle :r.
                """
                convert r, into c.
                """
                trapped_input = input
                if type == ':r.':
                    hgvs_input = hp.parse_hgvs_variant(input)  # Traps the hgvs variant of r. for further use
                    # Change to coding variant
                    type = ':c.'
                    # Change input to reflect!
                    try:
                        hgvs_c = va_func.hgvs_r_to_c(hgvs_input)
                    except hgvs.exceptions.HGVSDataNotAvailableError as e:
                        error = str(e)
                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                        logger.warning(str(error))
                        continue
                    input = str(hgvs_c)
                    variant = str(hgvs_c)

                # COLLECT gene symbol, name and ACCESSION INFORMATION
                # Gene symbol
                logger.trace("Handled mitochondrial variants", validation)
                """
                Identifies the transcript reference sequence name and HGNC gene symbol
                """
                if (type != ':g.'):
                    error = 'false'
                    hgvs_vt = hp.parse_hgvs_variant(variant)
                    try:
                        tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
                    except hgvs.exceptions.HGVSError as e:
                        error = str(e)
                        logger.warning(error)
                    if error != 'false':
                        error = 'Please inform UTA admin of the following error: ' + str(error)
                        issue_link = "https://bitbucket.org/biocommons/uta/issues?status=new&status=open"
                        reason = "VariantValidator cannot recover information for transcript " + str(
                            hgvs_vt.ac) + ' beacuse it is not available in the Universal Transcript Archive'
                        validation['warnings'] = validation['warnings'] + ': ' + str(reason)
                        logger.warning(str(reason) + ": " + str(error))
                        continue
                    else:
                        # Get hgnc Gene name from command
                        hgnc = tx_id_info[6]
                        issue_link = 'false'

                    # ACCESS THE GENE INFORMATION RECORDS ON THE UTA DATABASE
                    # Refseq accession
                    tx_for_gene = va_func.tx_for_gene(hgnc, hdp)
                    refseq_ac = va_func.ng_extract(tx_for_gene)

                    # Additional gene info
                    gene_info = hdp.get_gene_info(hgnc)
                    # Chromosomal location
                    try:
                        maploc = gene_info[1]
                    except:
                        maploc = ''
                    chr_loc = ("Chromosome location: " + maploc)

                    # Get accurate transcript descriptions from the relevant databases
                    # RefSeq databases
                    if alt_aln_method != 'genebuild':
                        # Gene description  - requires GenBank search to get all the required info, i.e. transcript variant ID
                        # accession number
                        hgvs_object = hp.parse_hgvs_variant(variant)
                        accession = hgvs_object.ac
                        # Look for the accession in our database
                        # Connect to database and send request
                        entry = va_dbCrl.data.in_entries(accession, 'transcript_info')

                        # Analyse the returned data and take the necessary actions
                        # If the error key exists
                        if 'error' in entry:
                            # Open a hgvs exception log file in append mode
                            error = entry['description']
                            validation['warnings'] = validation['warnings'] + ': ' + str(
                                error) + ': A Database error occurred, please contact admin'
                            logger.warning(str(error) + ": A Database error occurred, please contact admin")
                            continue

                        # If the accession key is found
                        elif 'accession' in entry:
                            description = entry['description']
                            # If the current entry is too old
                            if entry['expiry'] == 'true':
                                dbaction = 'update'
                                try:
                                    entry = va_btch.data_add(input=input, alt_aln_method=alt_aln_method,
                                                             accession=accession, dbaction=dbaction, hp=hp, evm=evm,
                                                             hdp=hdp)
                                except hgvs.exceptions.HGVSError as e:
                                    error = 'Transcript %s is not currently supported' % (accession)
                                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                    logger.warning(str(error))
                                    continue
                                except Exception as e:
                                    error = 'Unable to assign transcript identity records to ' + accession + ', potentially an obsolete record :'
                                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                    logger.warning(str(error))
                                    continue
                                hgnc_gene_info = entry['description']
                            else:
                                hgnc_gene_info = entry['description']
                        # If the none key is found add the description to the database
                        elif 'none' in entry:
                            dbaction = 'insert'
                            try:
                                entry = va_btch.data_add(input=input, alt_aln_method=alt_aln_method,
                                                         accession=accession, dbaction=dbaction, hp=hp, evm=evm,
                                                         hdp=hdp)
                            except Exception as e:
                                logger.warning(str(e))
                                error = 'Unable to assign transcript identity records to ' + accession + ', potentially an obsolete record :'
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                            hgnc_gene_info = entry['description']

                        # If no correct keys are found
                        else:
                            # Open a hgvs exception log file in append mode
                            error = 'Unknown error type'
                            validation['warnings'] = validation['warnings'] + ': ' + str(
                                error) + ': A Database error occurred, please contact admin'
                            logger.warning(str(error))
                            continue

                    # Ensembl databases
                    else:
                        # accession number
                        hgvs_object = hp.parse_hgvs_variant(variant)
                        accession = hgvs_object.ac
                        # Look for the accession in our database
                        # Connect to database and send request
                        entry = va_dbCrl.data.in_entries(accession, 'transcript_info')

                        # Analyse the returned data and take the necessary actions
                        # If the error key exists
                        if 'error' in entry:
                            # Open a hgvs exception log file in append mode
                            error = entry['description']
                            validation['warnings'] = validation['warnings'] + ': ' + str(
                                error) + ': A Database error occurred, please contact admin'
                            logger.warning(str(error))
                            continue

                        # If the accession key is found
                        elif 'accession' in entry:
                            description = entry['description']
                            # If the current entry is too old
                            if entry['expiry'] == 'true':
                                dbaction = 'update'
                                entry = va_btch.data_add(input=input, alt_aln_method=alt_aln_method,
                                                         accession=accession, dbaction=dbaction, hp=hp, evm=evm,
                                                         hdp=hdp)
                                hgnc_gene_info = entry['description']
                            else:
                                hgnc_gene_info = entry['description']
                        # If the none key is found add the description to the database
                        elif 'none' in entry:
                            dbaction = 'insert'
                            try:
                                entry = va_btch.data_add(input=input, alt_aln_method=alt_aln_method,
                                                         accession=accession, dbaction=dbaction, hp=hp, evm=evm,
                                                         hdp=hdp)
                            except Exception as e:
                                logger.warning(str(e))
                                error = 'Unable to assign transcript identity records to ' + accession + ', potentially an obsolete record :'
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                            hgnc_gene_info = entry['description']

                        # If no correct keys are found
                        else:
                            # Open a hgvs exception log file in append mode
                            error = 'Unknown error type'
                            validation['warnings'] = validation['warnings'] + ': ' + str(
                                error) + ': A Database error occurred, please contact admin'
                            logger.warning(str(error))
                            continue

                # Genomic type variants will need to be mapped to transcripts
                """
                The following section is used to project genomic variants accurately onto 
                all relevant transcripts
                """

                if (type == ':g.'):
                    g_query = hp.parse_hgvs_variant(variant)

                    # Genomic coordinates can be validated immediately
                    error = 'false'
                    try:
                        vr.validate(g_query)
                    except hgvs.exceptions.HGVSError as e:
                        error = str(e)
                    except KeyError:
                        error = 'Reference sequence ' + hgvs_genomic.ac + ' is either not supported or does not exist'
                    if error != 'false':
                        reason = 'Invalid variant description'
                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                        logger.warning(str(error))
                        continue
                    else:
                        pass

                    # Set test to see if Norm alters the coords
                    g_test = hn.normalize(g_query)

                    # Perform test
                    if g_query.posedit.pos != g_test.posedit.pos:
                        # validation['warnings'] = validation['warnings'] + ': ' + 'Input variant description normalized to ' + str(g_test)
                        hgvs_genomic = g_test
                    else:
                        hgvs_genomic = g_query

                    # Collect rel_var
                    # rel_var is a keyworded list of relevant transcripts with associated coding variants
                    """
                    Initial simple projection from the provided g. position all overlapping
                    transcripts
                    """
                    rel_var = va_func.relevant_transcripts(hgvs_genomic, evm, hdp, alt_aln_method, reverse_normalizer, hp)

                    # Double check rel_vars have not been missed when mapping from a RefSeqGene
                    if len(rel_var) != 0 and re.match('NG_', str(hgvs_genomic.ac)):
                        for var in rel_var:
                            hgvs_coding_variant = hp.parse_hgvs_variant(var)
                            try:
                                hgvs_genomic = va_func.myevm_t_to_g(hgvs_coding_variant, hdp, no_norm_evm,
                                                                    primary_assembly, vm, hp, hn, sf, nr_vm)
                            except hgvs.exceptions.HGVSError as e:
                                try_rel_var = []
                            else:
                                try_rel_var = va_func.relevant_transcripts(hgvs_genomic, evm, hdp, alt_aln_method,
                                                                           reverse_normalizer, hp)
                            if len(try_rel_var) > len(rel_var):
                                rel_var = try_rel_var
                                break
                            else:
                                continue

                    #  Tripple check this assumption by querying the gene position database
                    if len(rel_var) == 0:
                        vcf_dict = va_H2V.hgvs2vcf(hgvs_genomic, primary_assembly, reverse_normalizer, sf)
                        not_di = str(hgvs_genomic.ac) + ':g.' + str(vcf_dict['pos']) + '_' + str(
                            int(vcf_dict['pos']) + (len(vcf_dict['ref']) - 1)) + 'del' + vcf_dict['ref'] + 'ins' + \
                                 vcf_dict['alt']
                        hgvs_not_di = hp.parse_hgvs_variant(not_di)
                        rel_var = va_func.relevant_transcripts(hgvs_not_di, evm, hdp, alt_aln_method,
                                                               reverse_normalizer, hp)

                    # list return statements
                    """
                    If mapping to transcripts has been unsuccessful, provide relevant details
                    """
                    if len(rel_var) == 0:

                        # Check for NG_
                        rsg = re.compile('^NG_')
                        if rsg.search(variant):
                            # parse
                            hgvs_refseqgene = hp.parse_hgvs_variant(variant)
                            # Convert to chromosomal position
                            refseqgene_data = va_g2g.rsg_to_chr(hgvs_refseqgene, primary_assembly, hn, vr)
                            # There should only ever be one description returned
                            refseqgene_data = refseqgene_data[0]

                            # Extract data
                            if refseqgene_data['valid'] == 'true':
                                input = refseqgene_data['hgvs_genomic']
                                # re_submit
                                # Tag the line so that it is not written out
                                validation['warnings'] = validation[
                                                             'warnings'] + ': ' + variant + ' automapped to genome position ' + str(
                                    input)
                                query = {'quibble': input, 'id': validation['id'], 'warnings': validation['warnings'],
                                         'description': '', 'coding': '', 'coding_g': '', 'genomic_r': '',
                                         'genomic_g': '', 'protein': '', 'write': 'true',
                                         'primary_assembly': primary_assembly, 'order': ordering}
                                coding = 'intergenic'
                                batch_list.append(query)
                            else:
                                error = 'Mapping unavailable for RefSeqGene ' + variant + ' using alignment method = ' + alt_aln_method
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue

                        # Chromosome build is not supported or intergenic???
                        else:
                            sfm = va_scb.supported_for_mapping(hgvs_genomic.ac, primary_assembly)
                            if sfm == 'true':
                                try:
                                    vr.validate(hgvs_genomic)
                                except hgvs.exceptions.HGVSError as e:
                                    error = str(e)
                                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                    logger.warning(str(error))
                                    continue
                                else:
                                    # Map to RefSeqGene if available
                                    refseqgene_data = va_g2g.chr_to_rsg(hgvs_genomic, hn, vr)
                                    rsg_data = ''
                                    # Example {'gene': 'NTHL1', 'hgvs_refseqgene': 'NG_008412.1:g.3455_3464delCAAACACACA', 'valid': 'true'}
                                    for data in refseqgene_data:
                                        if data['valid'] == 'true':
                                            data['hgvs_refseqgene'] = hp.parse_hgvs_variant(data['hgvs_refseqgene'])
                                            data['hgvs_refseqgene'] = valstr(data['hgvs_refseqgene'])
                                            rsg_data = rsg_data + data['hgvs_refseqgene'] + ' (' + data['gene'] + '), '

                                    error = 'No transcripts found that fully overlap the described variation in the genomic sequence'
                                    # set output type flag
                                    set_output_type_flag = 'intergenic'
                                    # set genomic and where available RefSeqGene outputs
                                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                    validation['genomic_g'] = valstr(hgvs_genomic)
                                    validation['genomic_r'] = str(rsg_data.split('(')[0])
                                    logger.warning(str(error))
                                    continue
                            else:
                                error = 'Please ensure the requested chromosome version relates to a supported genome build. Supported genome builds are: GRCh37, GRCh38, hg19 and hg38'
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue

                    else:
                        # Tag the line so that it is not written out
                        validation['write'] = 'false'

                        """
                        Gap aware projection from g. to c.
                        """

                        # Set variables for problem specific warnings
                        gapped_alignment_warning = ''
                        corrective_action_taken = ''
                        gapped_transcripts = ''
                        auto_info = ''

                        # Create a pseudo VCF so that normalization can be applied and a delins can be generated
                        hgvs_genomic_variant = hgvs_genomic
                        # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
                        reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic_variant)
                        hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

                        # VCF
                        vcf_dict = va_H2V.hgvs2vcf(reverse_normalized_hgvs_genomic, primary_assembly,
                                                   reverse_normalizer, sf)
                        chr = vcf_dict['chr']
                        pos = vcf_dict['pos']
                        ref = vcf_dict['ref']
                        alt = vcf_dict['alt']

                        # Generate an end position
                        end = str(int(pos) + len(ref) - 1)
                        pos = str(pos)

                        # take a look at the input genomic variant for potential base salvage
                        stash_ac = vcf_dict['chr']
                        stash_pos = int(vcf_dict['pos'])
                        stash_ref = vcf_dict['ref']
                        stash_alt = vcf_dict['alt']
                        stash_end = end
                        # Re-Analyse genomic positions
                        if re.match('NG_', str(stash_input)):
                            c = hp.parse_hgvs_variant(rel_var[0])
                            if hasattr(c.posedit.edit, 'ref') and c.posedit.edit.ref is not None:
                                c.posedit.edit.ref = c.posedit.edit.ref.upper()
                            if hasattr(c.posedit.edit, 'alt') and c.posedit.edit.alt is not None:
                                c.posedit.edit.alt = c.posedit.edit.alt.upper()
                            stash_input = va_func.myevm_t_to_g(c, hdp, no_norm_evm, primary_assembly, vm, hp, hn, sf,
                                                               nr_vm)
                        if re.match('NC_', str(stash_input)) or re.match('NT_', str(stash_input)) or re.match('NW_',
                                                                                                              str(
                                                                                                                      stash_input)):
                            try:
                                hgvs_stash = hp.parse_hgvs_variant(stash_input)
                            except:
                                hgvs_stash = stash_input
                            if hasattr(hgvs_stash.posedit.edit, 'ref') and hgvs_stash.posedit.edit.ref is not None:
                                hgvs_stash.posedit.edit.ref = hgvs_stash.posedit.edit.ref.upper()
                            if hasattr(hgvs_stash.posedit.edit, 'alt') and hgvs_stash.posedit.edit.alt is not None:
                                hgvs_stash.posedit.edit.alt = hgvs_stash.posedit.edit.alt.upper()

                            stash_ac = hgvs_stash.ac
                            # MAKE A NO NORM HGVS2VCF
                            stash_dict = va_H2V.pos_lock_hgvs2vcf(hgvs_stash, primary_assembly, reverse_normalizer, sf)
                            stash_ac = hgvs_stash.ac
                            stash_pos = int(stash_dict['pos'])
                            stash_ref = stash_dict['ref']
                            stash_alt = stash_dict['alt']
                            # Generate an end position
                            stash_end = str(stash_pos + len(stash_ref) - 1)

                        # Store a not real deletion insertion
                        stored_hgvs_not_delins = hp.parse_hgvs_variant(str(
                            hgvs_genomic_5pr.ac) + ':' + hgvs_genomic_5pr.type + '.' + pos + '_' + end + 'del' + ref + 'ins' + alt)
                        stash_hgvs_not_delins = hp.parse_hgvs_variant(
                            stash_ac + ':' + hgvs_genomic_5pr.type + '.' + str(
                                stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)

                        # Set non-valid caution to false
                        non_valid_caution = 'false'

                        # make an empty rel_var
                        nw_rel_var = []

                        # loop through rel_var and amend where required
                        for var in rel_var:
                            # Store the current hgvs:c. description
                            saved_hgvs_coding = hp.parse_hgvs_variant(var)

                            # Remove un-selected transcripts
                            if select_transcripts != 'all':
                                tx_ac = saved_hgvs_coding.ac
                                # If it's in the selected tx dict, keep it
                                if tx_ac.split('.')[0] in select_transcripts_dict.keys():
                                    pass
                                # If not get rid of it!
                                else:
                                    continue

                            # Get orientation of the gene wrt genome and a list of exons mapped to the genome
                            ori = va_func.tx_exons(tx_ac=saved_hgvs_coding.ac, alt_ac=hgvs_genomic_5pr.ac,
                                                   alt_aln_method=alt_aln_method, hdp=hdp)
                            orientation = int(ori[0]['alt_strand'])
                            intronic_variant = 'false'

                            if orientation == -1:
                                # position genomic at its most 5 prime position
                                try:
                                    query_genomic = reverse_normalizer.normalize(hgvs_genomic)
                                except:
                                    query_genomic = hgvs_genomic
                                # Map to the transcript ant test for movement
                                try:
                                    hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
                                except hgvs.exceptions.HGVSError as e:
                                    hgvs_seek_var = saved_hgvs_coding
                                else:
                                    seek_var = valstr(hgvs_seek_var)
                                    seek_ac = str(hgvs_seek_var.ac)
                                if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                                        saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                                        hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                                        saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset):
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
                            else:
                                seek_var = valstr(hgvs_seek_var)
                                seek_ac = str(hgvs_seek_var.ac)
                                if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                                        saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                                        hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                                        saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset):
                                    pass
                                else:
                                    hgvs_seek_var = saved_hgvs_coding

                            try:
                                intron_test = hn.normalize(hgvs_seek_var)
                            except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                                error = str(e)
                                if re.match('Normalization of intronic variants is not supported', error) or re.match(
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
                                if re.search('\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\d+\-', str(
                                        hgvs_seek_var.posedit.pos)) or re.search('\*\d+\+', str(
                                    hgvs_seek_var.posedit.pos)) or re.search('\*\d+\-',
                                                                             str(hgvs_seek_var.posedit.pos)):
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

                            if re.search('\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\d+\-', str(
                                    hgvs_seek_var.posedit.pos)) or re.search('\*\d+\+', str(
                                hgvs_seek_var.posedit.pos)) or re.search('\*\d+\-', str(hgvs_seek_var.posedit.pos)):
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

                            # If exonic, process
                            if intronic_variant != 'true':
                                # map form reverse normalized g. to c.
                                hgvs_from_5n_g = no_norm_evm.g_to_t(hgvs_genomic_5pr, saved_hgvs_coding.ac)

                                # Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping transcript lengths
                                disparity_deletion_in = ['false', 'false']
                                if stored_hgvs_not_delins != '':
                                    # Refresh hgvs_not_delins from stored_hgvs_not_delins
                                    hgvs_not_delins = copy.deepcopy(stored_hgvs_not_delins)
                                    # This test will only occur in dup of single base, insertion or substitution
                                    if not re.search('_', str(hgvs_not_delins.posedit.pos)):
                                        if re.search('dup', hgvs_genomic_5pr.posedit.edit.type) or re.search('ins',
                                                                                                             hgvs_genomic_5pr.posedit.edit.type):
                                            # For gap in chr, map to t. - but becaouse we have pushed to 5 prime by norm, add 1 to end pos
                                            plussed_hgvs_not_delins = copy.deepcopy(hgvs_not_delins)
                                            plussed_hgvs_not_delins.posedit.pos.end.base = plussed_hgvs_not_delins.posedit.pos.end.base + 1
                                            plussed_hgvs_not_delins.posedit.edit.ref = ''
                                            transcript_variant = no_norm_evm.g_to_t(plussed_hgvs_not_delins,
                                                                                    str(saved_hgvs_coding.ac))
                                            if ((
                                                    transcript_variant.posedit.pos.end.base - transcript_variant.posedit.pos.start.base) > (
                                                    hgvs_genomic_5pr.posedit.pos.end.base - hgvs_genomic_5pr.posedit.pos.start.base)):
                                                if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                                                    hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                    start = hgvs_not_delins.posedit.pos.start.base - 1
                                                    end = hgvs_not_delins.posedit.pos.end.base
                                                    ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                                    hgvs_not_delins.posedit.edit.ref = ref_bases
                                                    hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                       :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                             1:] + ref_bases[1:]
                                                elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search(
                                                        'del', str(hgvs_genomic_5pr.posedit.edit)):
                                                    hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                elif re.search('ins',
                                                               str(hgvs_genomic_5pr.posedit.edit)) and not re.search(
                                                    'del', str(hgvs_genomic_5pr.posedit.edit)):
                                                    hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                    start = hgvs_not_delins.posedit.pos.start.base - 1
                                                    end = hgvs_not_delins.posedit.pos.end.base
                                                    ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                                    hgvs_not_delins.posedit.edit.ref = ref_bases
                                                    hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                       :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                             1:] + ref_bases[1:]
                                            else:
                                                if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                                                    hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                    start = hgvs_not_delins.posedit.pos.start.base - 1
                                                    end = hgvs_not_delins.posedit.pos.end.base
                                                    ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                                    hgvs_not_delins.posedit.edit.ref = ref_bases
                                                    hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                       :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                             1:] + ref_bases[1:]
                                                elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search(
                                                        'del', str(hgvs_genomic_5pr.posedit.edit)):
                                                    hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                elif re.search('ins',
                                                               str(hgvs_genomic_5pr.posedit.edit)) and not re.search(
                                                    'del', str(hgvs_genomic_5pr.posedit.edit)):
                                                    hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                    start = hgvs_not_delins.posedit.pos.start.base - 1
                                                    end = hgvs_not_delins.posedit.pos.end.base
                                                    ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                                    hgvs_not_delins.posedit.edit.ref = ref_bases
                                                    hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                       :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                             1:] + ref_bases[1:]
                                        else:
                                            pass
                                    else:
                                        pass

                                    try:
                                        tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, saved_hgvs_coding.ac)
                                    except hgvs.exceptions.HGVSInvalidIntervalError:
                                        tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_genomic_5pr, saved_hgvs_coding.ac)
                                    except hgvs.exceptions.HGVSError:
                                        if str(e) == 'start or end or both are beyond the bounds of transcript record':
                                            tx_hgvs_not_delins = saved_hgvs_coding

                                    # Create normalized version of tx_hgvs_not_delins
                                    rn_tx_hgvs_not_delins = copy.deepcopy(tx_hgvs_not_delins)
                                    # Check for +ve base and adjust
                                    if (re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)) or re.search('\-',
                                                                                                                   str(
                                                                                                                       rn_tx_hgvs_not_delins.posedit.pos.start))) and (
                                            re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)) or re.search(
                                        '\-', str(rn_tx_hgvs_not_delins.posedit.pos.end))):
                                        # Remove offsetting to span the gap
                                        rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                        rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                                        rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
                                        rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                        try:
                                            rn_tx_hgvs_not_delins.posedit.edit.alt = ''
                                        except:
                                            exceptPass()
                                    elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                                        # move tx end base to next available non-offset base
                                        rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
                                        rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                                        rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                        if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                            test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                        else:
                                            test_tx_var = rn_tx_hgvs_not_delins
                                        # re-make genomic and tx
                                        hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm,
                                                                               primary_assembly, vm, hp, hn, sf, nr_vm)
                                        rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                   str(saved_hgvs_coding.ac))
                                    elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                                        # move tx start base to previous available non-offset base
                                        rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                        rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                        if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                            test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                        else:
                                            test_tx_var = rn_tx_hgvs_not_delins
                                        # re-make genomic and tx
                                        hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm,
                                                                               primary_assembly, vm, hp, hn, sf, nr_vm)
                                        rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                   str(saved_hgvs_coding.ac))
                                        rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                    #                                     else:
                                    #                                         pass

                                    # Check for -ve base and adjust
                                    elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search('\-',
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
                                            exceptPass()
                                    elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                                        # move tx end base back to next available non-offset base
                                        rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                                        # Delete the ref
                                        rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                        # Add the additional base to the ALT
                                        start = rn_tx_hgvs_not_delins.posedit.pos.end.base - 1
                                        end = rn_tx_hgvs_not_delins.posedit.pos.end.base
                                        ref_bases = sf.fetch_seq(str(tx_hgvs_not_delins.ac), start, end)
                                        rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + ref_bases
                                        if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                            test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                        else:
                                            test_tx_var = rn_tx_hgvs_not_delins
                                        # re-make genomic and tx
                                        hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm,
                                                                               primary_assembly, vm, hp, hn, sf, nr_vm)
                                        rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                   str(saved_hgvs_coding.ac))
                                    elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                                        # move tx start base to previous available non-offset base
                                        rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                        rn_tx_hgvs_not_delins.posedit.pos.start.base = rn_tx_hgvs_not_delins.posedit.pos.start.base - 1
                                        rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                        if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                            test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                        else:
                                            test_tx_var = rn_tx_hgvs_not_delins
                                        # re-make genomic and tx
                                        hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm,
                                                                               primary_assembly, vm, hp, hn, sf, nr_vm)
                                        rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                   str(saved_hgvs_coding.ac))
                                        rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                    else:
                                        pass

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
                                        hgvs_stash_t = vm.g_to_t(stash_hgvs_not_delins, saved_hgvs_coding.ac)
                                        if len(stash_hgvs_not_delins.posedit.edit.ref) > len(
                                                hgvs_stash_t.posedit.edit.ref):
                                            try:
                                                hn.normalize(hgvs_stash_t)
                                            except:
                                                exceptPass()
                                            else:
                                                gap_length = len(stash_hgvs_not_delins.posedit.edit.ref) - len(
                                                    hgvs_stash_t.posedit.edit.ref)
                                                disparity_deletion_in = ['transcript', gap_length]
                                                try:
                                                    tx_hgvs_not_delins = vm.c_to_n(hgvs_stash_t)
                                                except:
                                                    tx_hgvs_not_delins = hgvs_stash_t
                                                hgvs_not_delins = stash_hgvs_not_delins
                                        elif hgvs_stash_t.posedit.pos.start.offset != 0 or hgvs_stash_t.posedit.pos.end.offset != 0:
                                            disparity_deletion_in = ['transcript', 'Requires Analysis']
                                            try:
                                                tx_hgvs_not_delins = vm.c_to_n(hgvs_stash_t)
                                            except:
                                                tx_hgvs_not_delins = hgvs_stash_t
                                            hgvs_not_delins = stash_hgvs_not_delins
                                            hgvs_genomic_5pr = stash_hgvs_not_delins
                                        else:
                                            pass

                                # Final sanity checks
                                try:
                                    vm.g_to_t(hgvs_not_delins, tx_hgvs_not_delins.ac)
                                except Exception as e:
                                    if str(e) == 'start or end or both are beyond the bounds of transcript record':
                                        hgvs_not_delins = saved_hgvs_coding
                                        disparity_deletion_in = ['false', 'false']
                                    logger.warning(str(e))
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
                                            hgvs_not_delins = saved_hgvs_coding
                                            disparity_deletion_in = ['false', 'false']
                                        elif re.match('Normalization of intronic variants is not supported', error):
                                            # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                                            disparity_deletion_in = ['transcript', 'Requires Analysis']
                                    logger.warning(error)
                                # Pre-processing of tx_hgvs_not_delins
                                try:
                                    if tx_hgvs_not_delins.posedit.edit.alt is None:
                                        tx_hgvs_not_delins.posedit.edit.alt = ''
                                except Exception as e:
                                    if str(e) == "'Dup' object has no attribute 'alt'":
                                        tx_hgvs_not_delins_delins_from_dup = tx_hgvs_not_delins.ac + ':' + tx_hgvs_not_delins.type + '.' + str(
                                            tx_hgvs_not_delins.posedit.pos.start) + '_' + str(
                                            tx_hgvs_not_delins.posedit.pos.end) + 'del' + tx_hgvs_not_delins.posedit.edit.ref + 'ins' + tx_hgvs_not_delins.posedit.edit.ref + tx_hgvs_not_delins.posedit.edit.ref
                                        tx_hgvs_not_delins = hp.parse_hgvs_variant(tx_hgvs_not_delins_delins_from_dup)

                                # GAP IN THE TRANSCRIPT DISPARITY DETECTED
                                if disparity_deletion_in[0] == 'transcript':
                                    gap_position = ''
                                    gapped_alignment_warning = str(
                                        hgvs_genomic_5pr) + ' does not represent a true variant because it is an artefact of aligning the transcripts listed below with genome build ' + primary_assembly

                                    # ANY VARIANT WHOLLY WITHIN THE GAP
                                    if (re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) or re.search('\-',
                                                                                                                str(
                                                                                                                    tx_hgvs_not_delins.posedit.pos.start))) and (
                                            re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) or re.search('\-',
                                                                                                                  str(
                                                                                                                      tx_hgvs_not_delins.posedit.pos.end))):
                                        gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
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
                                                tx_gap_fill_variant = hp.parse_hgvs_variant(
                                                    tx_gap_fill_variant_delins_from_dup)

                                        # Identify which half of the NOT-intron the start position of the variant is in
                                        if re.search('\-', str(tx_gap_fill_variant.posedit.pos.start)):
                                            tx_gap_fill_variant.posedit.pos.start.base = tx_gap_fill_variant.posedit.pos.start.base - 1
                                            tx_gap_fill_variant.posedit.pos.start.offset = int('0')  # int('+1')
                                            tx_gap_fill_variant.posedit.pos.end.offset = int('0')  # int('-1')
                                            tx_gap_fill_variant.posedit.edit.alt = ''
                                            tx_gap_fill_variant.posedit.edit.ref = ''
                                        elif re.search('\+', str(tx_gap_fill_variant.posedit.pos.start)):
                                            tx_gap_fill_variant.posedit.pos.start.offset = int('0')  # int('+1')
                                            tx_gap_fill_variant.posedit.pos.end.base = tx_gap_fill_variant.posedit.pos.end.base + 1
                                            tx_gap_fill_variant.posedit.pos.end.offset = int('0')  # int('-1')
                                            tx_gap_fill_variant.posedit.edit.alt = ''
                                            tx_gap_fill_variant.posedit.edit.ref = ''

                                        try:
                                            tx_gap_fill_variant = vm.n_to_c(tx_gap_fill_variant)
                                        except:
                                            exceptPass()
                                        genomic_gap_fill_variant = vm.t_to_g(tx_gap_fill_variant,
                                                                             reverse_normalized_hgvs_genomic.ac)
                                        genomic_gap_fill_variant.posedit.edit.alt = genomic_gap_fill_variant.posedit.edit.ref

                                        try:
                                            c_tx_hgvs_not_delins = vm.n_to_c(tx_hgvs_not_delins)
                                        except Exception:
                                            c_tx_hgvs_not_delins = copy.copy(tx_hgvs_not_delins)
                                        genomic_gap_fill_variant_alt = vm.t_to_g(c_tx_hgvs_not_delins,
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
                                                genomic_gap_fill_variant = hp.parse_hgvs_variant(
                                                    genomic_gap_fill_variant_delins_from_dup)
                                                genomic_gap_fill_variant_alt_delins_from_dup = genomic_gap_fill_variant_alt.ac + ':' + genomic_gap_fill_variant_alt.type + '.' + str(
                                                    genomic_gap_fill_variant_alt.posedit.pos.start.base) + '_' + str(
                                                    genomic_gap_fill_variant_alt.posedit.pos.end.base) + 'del' + genomic_gap_fill_variant_alt.posedit.edit.ref + 'ins' + genomic_gap_fill_variant_alt.posedit.edit.ref + genomic_gap_fill_variant_alt.posedit.edit.ref
                                                genomic_gap_fill_variant_alt = hp.parse_hgvs_variant(
                                                    genomic_gap_fill_variant_alt_delins_from_dup)

                                        # Correct insertion alts
                                        if genomic_gap_fill_variant_alt.posedit.edit.type == 'ins':
                                            append_ref = sf.fetch_seq(genomic_gap_fill_variant_alt.ac,
                                                                      genomic_gap_fill_variant_alt.posedit.pos.start.base - 1,
                                                                      genomic_gap_fill_variant_alt.posedit.pos.end.base)
                                            genomic_gap_fill_variant_alt.posedit.edit.alt = append_ref[
                                                                                                0] + genomic_gap_fill_variant_alt.posedit.edit.alt + \
                                                                                            append_ref[1]

                                        # Split the reference and replacing alt sequence into a dictionary
                                        reference_bases = list(genomic_gap_fill_variant.posedit.edit.ref)
                                        if genomic_gap_fill_variant_alt.posedit.edit.alt is not None:
                                            alternate_bases = list(genomic_gap_fill_variant_alt.posedit.edit.alt)
                                        else:
                                            # Deletions with no ins
                                            pre_alternate_bases = list(genomic_gap_fill_variant_alt.posedit.edit.ref)
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
                                        for integer in range(genomic_gap_fill_variant_alt.posedit.pos.start.base,
                                                             genomic_gap_fill_variant_alt.posedit.pos.end.base + 1, 1):
                                            if integer == alt_start:
                                                alt_base_dict[integer] = str(''.join(alternate_bases))
                                            else:
                                                alt_base_dict[integer] = 'X'

                                        # Generate the alt sequence
                                        alternate_sequence_bases = []
                                        for integer in range(genomic_gap_fill_variant.posedit.pos.start.base,
                                                             genomic_gap_fill_variant.posedit.pos.end.base + 1, 1):
                                            if integer in alt_base_dict.keys():
                                                alternate_sequence_bases.append(alt_base_dict[integer])
                                            else:
                                                alternate_sequence_bases.append(ref_base_dict[integer])
                                        alternate_sequence = ''.join(alternate_sequence_bases)
                                        alternate_sequence = alternate_sequence.replace('X', '')

                                        # Add the new alt to the gap fill variant and generate transcript variant
                                        genomic_gap_fill_variant.posedit.edit.alt = alternate_sequence
                                        hgvs_refreshed_variant = vm.g_to_t(genomic_gap_fill_variant,
                                                                           tx_gap_fill_variant.ac)

                                        # Set warning
                                        gap_size = str(len(genomic_gap_fill_variant.posedit.edit.ref) - 2)
                                        disparity_deletion_in[1] = [gap_size]
                                        auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                                            stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + gap_size + ' genomic base(s) that fail to align to transcript ' + str(
                                            tx_hgvs_not_delins.ac)
                                        non_valid_caution = 'true'

                                        # Alignment position
                                        for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                        if re.match('NM_', str(for_location_c)):
                                            for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                        if re.match('\-', str(for_location_c.posedit.pos.start.offset)):
                                            gps = for_location_c.posedit.pos.start.base - 1
                                            gpe = for_location_c.posedit.pos.start.base
                                        else:
                                            gps = for_location_c.posedit.pos.start.base
                                            gpe = for_location_c.posedit.pos.start.base + 1
                                        gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                                        auto_info = auto_info + '%s' % (gap_position)
                                    else:
                                        if tx_hgvs_not_delins.posedit.pos.start.offset == 0 and tx_hgvs_not_delins.posedit.pos.end.offset == 0:
                                            # In this instance, we have identified a transcript gap but the n. version of
                                            # the transcript variant but do not have a position which actually hits the gap,
                                            # so the variant likely spans the gap, and is not picked up by an offset.
                                            try:
                                                c1 = vm.n_to_c(tx_hgvs_not_delins)
                                            except:
                                                c1 = tx_hgvs_not_delins
                                            g1 = nr_vm.t_to_g(c1, hgvs_genomic.ac)
                                            g3 = nr_vm.t_to_g(c1, hgvs_genomic.ac)
                                            g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                            ng2 = hn.normalize(g2)
                                            g3.posedit.pos.end.base = g3.posedit.pos.start.base + (
                                                    len(g3.posedit.edit.ref) - 1)
                                            try:
                                                c2 = vm.g_to_t(g3, c1.ac)
                                                if c2.posedit.pos.start.offset == 0 and c2.posedit.pos.end.offset == 0:
                                                    pass
                                                else:
                                                    tx_hgvs_not_delins = c2
                                                    try:
                                                        tx_hgvs_not_delins = vm.c_to_n(tx_hgvs_not_delins)
                                                    except hgvs.exceptions.HGVSError:
                                                        exceptPass()
                                            except hgvs.exceptions.HGVSInvalidVariantError:
                                                exceptPass()

                                        if re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                                                '\+', str(tx_hgvs_not_delins.posedit.pos.end)):
                                            auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                                                stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                                                disparity_deletion_in[
                                                    1]) + ' genomic base(s) that fail to align to transcript ' + str(
                                                tx_hgvs_not_delins.ac)
                                            non_valid_caution = 'true'
                                            try:
                                                c2 = vm.n_to_c(tx_hgvs_not_delins)
                                            except:
                                                c2 = tx_hgvs_not_delins
                                            c1 = copy.deepcopy(c2)
                                            c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                                            c1.posedit.pos.start.offset = 0
                                            c1.posedit.pos.end = c2.posedit.pos.start
                                            c1.posedit.edit.ref = ''
                                            c1.posedit.edit.alt = ''
                                            if orientation != -1:
                                                g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                g1.posedit.edit.alt = g1.posedit.edit.ref
                                            else:
                                                g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                g2.posedit.edit.alt = g2.posedit.edit.ref
                                            reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                            alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                            g3 = copy.deepcopy(g1)
                                            g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                            g3.posedit.edit.ref = reference
                                            g3.posedit.edit.alt = alternate
                                            c3 = vm.g_to_t(g3, c1.ac)
                                            hgvs_refreshed_variant = c3
                                            # Alignment position
                                            for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                            if re.match('NM_', str(for_location_c)):
                                                for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                                gps = for_location_c.posedit.pos.start.base
                                                gpe = for_location_c.posedit.pos.start.base + 1
                                            gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                                            # Warn update
                                            auto_info = auto_info + '%s' % (gap_position)
                                        elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search(
                                                '\+', str(tx_hgvs_not_delins.posedit.pos.start)):
                                            auto_info = auto_info + 'Genome position ' + str(
                                                stored_hgvs_not_delins.ac) + ':g.' + str(
                                                stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(
                                                disparity_deletion_in[1]) + '-bp gap in transcript ' + str(
                                                tx_hgvs_not_delins.ac)
                                            gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
                                            non_valid_caution = 'true'
                                            try:
                                                c1 = vm.n_to_c(tx_hgvs_not_delins)
                                            except:
                                                c1 = tx_hgvs_not_delins
                                            c2 = copy.deepcopy(c1)
                                            c2.posedit.pos.start = c1.posedit.pos.end
                                            c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                                            c2.posedit.pos.end.offset = 0
                                            c2.posedit.edit.ref = ''
                                            c2.posedit.edit.alt = ''
                                            if orientation != -1:
                                                g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                g2.posedit.edit.alt = g2.posedit.edit.ref
                                            else:
                                                g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                g1.posedit.edit.alt = g1.posedit.edit.ref
                                            reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                            alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                            g3 = copy.deepcopy(g1)
                                            g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                            g3.posedit.edit.ref = reference
                                            g3.posedit.edit.alt = alternate
                                            c3 = vm.g_to_t(g3, c1.ac)
                                            hgvs_refreshed_variant = c3
                                            # Alignment position
                                            for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                            if re.match('NM_', str(for_location_c)):
                                                for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                            gps = for_location_c.posedit.pos.end.base
                                            gpe = for_location_c.posedit.pos.end.base + 1
                                            gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                                            # Warn update
                                            auto_info = auto_info + '%s' % (gap_position)
                                        elif re.search('\-',
                                                       str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                                            '\-', str(tx_hgvs_not_delins.posedit.pos.end)):
                                            auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                                                stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                                                disparity_deletion_in[
                                                    1]) + ' genomic base(s) that fail to align to transcript ' + str(
                                                tx_hgvs_not_delins.ac)
                                            non_valid_caution = 'true'
                                            try:
                                                c2 = vm.n_to_c(tx_hgvs_not_delins)
                                            except:
                                                c2 = tx_hgvs_not_delins
                                            c1 = copy.deepcopy(c2)
                                            c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                                            c1.posedit.pos.start.offset = 0
                                            c1.posedit.pos.end = c2.posedit.pos.start
                                            c1.posedit.edit.ref = ''
                                            c1.posedit.edit.alt = ''
                                            if orientation != -1:
                                                g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                g1.posedit.edit.alt = g1.posedit.edit.ref
                                            else:
                                                g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                g2.posedit.edit.alt = g2.posedit.edit.ref
                                            reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                            alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                            g3 = copy.deepcopy(g1)
                                            g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                            g3.posedit.edit.ref = reference
                                            g3.posedit.edit.alt = alternate
                                            c3 = vm.g_to_t(g3, c1.ac)
                                            hgvs_refreshed_variant = c3
                                            # Alignment position
                                            for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                            if re.match('NM_', str(for_location_c)):
                                                for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                            gps = for_location_c.posedit.pos.start.base - 1
                                            gpe = for_location_c.posedit.pos.start.base
                                            gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                                            # Warn update
                                            auto_info = auto_info + '%s' % (gap_position)
                                        elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search(
                                                '\-', str(tx_hgvs_not_delins.posedit.pos.start)):
                                            auto_info = auto_info + 'Genome position ' + str(
                                                stored_hgvs_not_delins.ac) + ':g.' + str(
                                                stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(
                                                disparity_deletion_in[1]) + '-bp gap in transcript ' + str(
                                                tx_hgvs_not_delins.ac)
                                            gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
                                            non_valid_caution = 'true'
                                            try:
                                                c1 = vm.n_to_c(tx_hgvs_not_delins)
                                            except:
                                                c1 = tx_hgvs_not_delins
                                            c2 = copy.deepcopy(c1)
                                            c2.posedit.pos.start = c1.posedit.pos.end
                                            c2.posedit.pos.end.base = c1.posedit.pos.end.base
                                            c2.posedit.pos.end.offset = 0
                                            c2.posedit.edit.ref = ''
                                            c2.posedit.edit.alt = ''
                                            g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                            c2 = vm.g_to_t(g2, c2.ac)
                                            reference = c1.posedit.edit.ref + c2.posedit.edit.ref[1:]
                                            alternate = c1.posedit.edit.alt + c2.posedit.edit.ref[1:]
                                            c3 = copy.deepcopy(c1)
                                            c3.posedit.pos.end = c2.posedit.pos.end
                                            c3.posedit.edit.ref = ''  # reference
                                            c3.posedit.edit.alt = alternate
                                            hgvs_refreshed_variant = c3
                                            # Alignment position
                                            for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                            if re.match('NM_', str(for_location_c)):
                                                for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                            gps = for_location_c.posedit.pos.end.base - 1
                                            gpe = for_location_c.posedit.pos.end.base
                                            gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                                            # Warn update
                                            auto_info = auto_info + '%s' % (gap_position)
                                        else:
                                            auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                                                stored_hgvs_not_delins.posedit.pos) + ' contains ' + str(
                                                disparity_deletion_in[
                                                    1]) + ' genomic base(s) that fail to align to transcript ' + str(
                                                tx_hgvs_not_delins.ac) + '\n'
                                            hgvs_refreshed_variant = tx_hgvs_not_delins
                                            gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)

                                # GAP IN THE CHROMOSOME
                                elif disparity_deletion_in[0] == 'chromosome':
                                    # Set warning variables
                                    gap_position = ''
                                    gapped_alignment_warning = str(
                                        hgvs_genomic_5pr) + ' does not represent a true variant because it is an artefact of aligning the transcripts listed below with genome build ' + primary_assembly
                                    hgvs_refreshed_variant = tx_hgvs_not_delins
                                    # Warn
                                    auto_info = auto_info + str(hgvs_refreshed_variant.ac) + ':c.' + str(
                                        hgvs_refreshed_variant.posedit.pos) + ' contains ' + str(disparity_deletion_in[
                                                                                                     1]) + ' transcript base(s) that fail to align to chromosome ' + str(
                                        hgvs_genomic.ac) + '\n'
                                    gapped_transcripts = gapped_transcripts + str(hgvs_refreshed_variant.ac) + ' '
                                else:
                                    # Try the push
                                    hgvs_stash = copy.deepcopy(stash_hgvs_not_delins)
                                    stash_ac = hgvs_stash.ac
                                    # Make a hard left and hard right not delins g.
                                    stash_dict_right = va_H2V.hard_right_hgvs2vcf(hgvs_stash, primary_assembly, hn, sf)
                                    stash_pos_right = int(stash_dict_right['pos'])
                                    stash_ref_right = stash_dict_right['ref']
                                    stash_alt_right = stash_dict_right['alt']
                                    stash_end_right = str(stash_pos_right + len(stash_ref_right) - 1)
                                    stash_hgvs_not_delins_right = hp.parse_hgvs_variant(
                                        stash_ac + ':' + hgvs_stash.type + '.' + str(
                                            stash_pos_right) + '_' + stash_end_right + 'del' + stash_ref_right + 'ins' + stash_alt_right)
                                    stash_dict_left = va_H2V.hard_left_hgvs2vcf(hgvs_stash, primary_assembly,
                                                                                reverse_normalizer, sf)
                                    stash_pos_left = int(stash_dict_left['pos'])
                                    stash_ref_left = stash_dict_left['ref']
                                    stash_alt_left = stash_dict_left['alt']
                                    stash_end_left = str(stash_pos_left + len(stash_ref_left) - 1)
                                    stash_hgvs_not_delins_left = hp.parse_hgvs_variant(
                                        stash_ac + ':' + hgvs_stash.type + '.' + str(
                                            stash_pos_left) + '_' + stash_end_left + 'del' + stash_ref_left + 'ins' + stash_alt_left)
                                    # Map in-situ to the transcript left and right
                                    try:
                                        tx_hard_right = vm.g_to_t(stash_hgvs_not_delins_right, saved_hgvs_coding.ac)
                                    except Exception as e:
                                        tx_hard_right = saved_hgvs_coding
                                    else:
                                        normalize_stash_right = hn.normalize(stash_hgvs_not_delins_right)
                                        if str(normalize_stash_right.posedit) == str(stash_hgvs_not_delins.posedit):
                                            tx_hard_right = saved_hgvs_coding
                                    try:
                                        tx_hard_left = vm.g_to_t(stash_hgvs_not_delins_left, saved_hgvs_coding.ac)
                                    except Exception as e:
                                        tx_hard_left = saved_hgvs_coding
                                    else:
                                        normalize_stash_left = hn.normalize(stash_hgvs_not_delins_left)
                                        if str(normalize_stash_left.posedit) == str(stash_hgvs_not_delins.posedit):
                                            tx_hard_left = saved_hgvs_coding
                                    # The Logic - Currently limited to genome gaps
                                    if len(stash_hgvs_not_delins_right.posedit.edit.ref) < len(
                                            tx_hard_right.posedit.edit.ref):
                                        tx_hard_right = hn.normalize(tx_hard_right)
                                        gap_position = ''
                                        gapped_alignment_warning = str(
                                            hgvs_genomic_5pr) + ' may be an artefact of aligning the transcripts listed below with genome build ' + primary_assembly
                                        hgvs_refreshed_variant = tx_hard_right
                                        gapped_transcripts = gapped_transcripts + str(tx_hard_right.ac) + ' '
                                    elif len(stash_hgvs_not_delins_left.posedit.edit.ref) < len(
                                            tx_hard_left.posedit.edit.ref):
                                        tx_hard_left = hn.normalize(tx_hard_left)
                                        gap_position = ''
                                        gapped_alignment_warning = str(
                                            hgvs_genomic_5pr) + ' may be an artefact of aligning the transcripts listed below with genome build ' + primary_assembly
                                        hgvs_refreshed_variant = tx_hard_left
                                        gapped_transcripts = gapped_transcripts + str(tx_hard_left.ac) + ' '
                                    else:
                                        # Keep the same by re-setting rel_var
                                        hgvs_refreshed_variant = saved_hgvs_coding

                                # Edit the output
                                if re.match('NM_', str(hgvs_refreshed_variant.ac)) and not re.search('c', str(
                                        hgvs_refreshed_variant.type)):
                                    hgvs_refreshed_variant = evm.n_to_c(hgvs_refreshed_variant)
                                else:
                                    pass
                                try:
                                    hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
                                    if hgvs_refreshed_variant.posedit.edit.type == 'delins' and \
                                            hgvs_refreshed_variant.posedit.edit.ref[-1] == \
                                            hgvs_refreshed_variant.posedit.edit.alt[-1]:
                                        hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[
                                                                                  0:-1]
                                        hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[
                                                                                  0:-1]
                                        hgvs_refreshed_variant.posedit.pos.end.base = hgvs_refreshed_variant.posedit.pos.end.base - 1
                                        hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
                                    elif hgvs_refreshed_variant.posedit.edit.type == 'delins' and \
                                            hgvs_refreshed_variant.posedit.edit.ref[0] == \
                                            hgvs_refreshed_variant.posedit.edit.alt[0]:
                                        hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[
                                                                                  1:]
                                        hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[
                                                                                  1:]
                                        hgvs_refreshed_variant.posedit.pos.start.base = hgvs_refreshed_variant.posedit.pos.start.base + 1
                                        hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
                                except Exception as e:
                                    error = str(e)
                                    # Ensure the final variant is not intronic nor does it cross exon boundaries
                                    if re.match('Normalization of intronic variants is not supported',
                                                error) or re.match(
                                        'Unsupported normalization of variants spanning the exon-intron boundary',
                                        error):
                                        hgvs_refreshed_variant = saved_hgvs_coding
                                    else:
                                        pass
                                    exceptPass()
                                # Send to empty nw_rel_var
                                nw_rel_var.append(hgvs_refreshed_variant)

                            # Otherwise these variants need to be set
                            else:
                                corrective_action_taken = ''
                                gapped_alignment_warning = ''
                                # Send to empty nw_rel_var
                                nw_rel_var.append(saved_hgvs_coding)

                        # Warn the user that the g. description is not valid
                        if gapped_alignment_warning != '':
                            if disparity_deletion_in[0] == 'transcript':
                                corrective_action_taken = 'Automap has deleted  ' + str(
                                    disparity_deletion_in[1]) + ' bp from chromosomal reference sequence ' + str(
                                    hgvs_genomic.ac) + ' to ensure perfect alignment with transcript reference sequence(s)' + gapped_transcripts
                            if disparity_deletion_in[0] == 'chromosome':
                                corrective_action_taken = 'Automap has added  ' + str(
                                    disparity_deletion_in[1]) + ' bp to chromosomal reference sequence ' + str(
                                    hgvs_genomic.ac) + ' to ensure perfect alignment with transcript reference sequence(s) ' + gapped_transcripts

                        # Add additional data to the front of automap
                        if auto_info != '':
                            automap = auto_info + '\n' + automap

                        rel_var = copy.deepcopy(nw_rel_var)

                        # Set the values and append to batch_list
                        for c_description in rel_var:
                            query = {'quibble': str(c_description), 'id': validation['id'],
                                     'warnings': validation['warnings'], 'description': '', 'coding': '',
                                     'coding_g': '', 'genomic_r': '', 'genomic_g': '', 'protein': '', 'write': 'true',
                                     'primary_assembly': primary_assembly, 'order': ordering}
                            batch_list.append(query)
                        logger.warning("Continue reached when mapping transcript types to variants")
                        # Call next description
                        continue
                # TYPE = :c.

                if type == ':c.' or type == ':n.':

                    # Flag for validation
                    valid = 'false'
                    # Collect information for genomic level validation
                    obj = hp.parse_hgvs_variant(variant)

                    tx_ac = obj.ac

                    # Do we keep it?
                    if select_transcripts != 'all':
                        if tx_ac in select_transcripts_dict_plus_version.keys():
                            pass
                        # If not get rid of it!
                        else:
                            # By marking it as Do Not Write and continuing through the validation loop
                            validation['write'] = 'false'
                            continue
                    else:
                        pass

                    # Set a cross_variant object
                    cross_variant = 'false'
                    # Se rec_var to '' so it can be updated later
                    rec_var = ''
                    try:
                        to_g = va_func.myevm_t_to_g(obj, hdp, no_norm_evm, primary_assembly, vm, hp, hn, sf, nr_vm)
                        genomic_ac = to_g.ac
                    except hgvs.exceptions.HGVSDataNotAvailableError as e:
                        if (re.search('~', str(e)) and re.search('Alignment is incomplete', str(e))) or re.match(
                                "No relevant genomic mapping options available", str(e)):
                            reason = 'Unable to map the input variant onto a genomic position'
                            if (re.search('~', str(e)) and re.search('Alignment is incomplete', str(e))):
                                error_list = str(e).split('~')[:-1]
                                combos = [
                                    'Full alignment data between the specified transcript reference sequence and all GRCh37 and GRCh38 genomic reference sequences (including alternate chromosome assemblies, patches and RefSeqGenes) are not available: Consequently the input variant description cannot be fully validated and is not supported: Use the Gene to Transcripts function to determine whether an updated transcript reference sequence is available']  # Partial alignment data is available for the following genomic reference sequences: ']
                                error = '; '.join(combos)
                                error = error.replace(': ;', ': ')
                            else:
                                error = str(e)
                                error = error + ': Consequently the input variant description cannot be fully validated and is not supported: Use the Gene to Transcripts function to determine whether an updated transcript reference sequence is available'
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            logger.warning(str(error))
                            continue
                        try:
                            gene_symbol = va_dbCrl.data.get_gene_symbol_from_transcriptID(tx_ac)
                        except:
                            gene_symbol = 'none'
                        if gene_symbol == 'none':
                            error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
                        else:
                            error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' or ' + gene_symbol + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                        logger.warning(str(error))
                        continue
                    except TypeError as e:
                        try:
                            gene_symbol = va_dbCrl.data.get_gene_symbol_from_transcriptID(tx_ac)
                        except:
                            gene_symbol = 'none'
                        if gene_symbol == 'none':
                            error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
                        else:
                            error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' or ' + gene_symbol + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                        logger.warning(str(error))
                        continue

                    # Get orientation of the gene wrt genome and a list of exons mapped to the genome
                    ori = va_func.tx_exons(tx_ac=tx_ac, alt_ac=genomic_ac, alt_aln_method=alt_aln_method, hdp=hdp)
                    orientation = int(ori[0]['alt_strand'])
                    intronic_variant = 'false'

                    # Collect variant sequence information via normalisation (normalizer) or if intronic via mapping
                    # INTRONIC OFFSETS - Required for Exon table
                    # Variable to collect offset to exon boundary
                    ex_offset = 0
                    plus = re.compile("\d\+\d")  # finds digit + digit
                    minus = re.compile("\d\-\d")  # finds digit - digit

                    geno = re.compile(':g.')
                    if plus.search(input) or minus.search(input):
                        es = re.compile('error')
                        if es.search(str(to_g)):
                            if alt_aln_method != 'genebuild':
                                error = "If the following error message does not address the issue and the problem persists please contact admin: " + to_g
                                reason = "An error has occurred"
                                excep = "%s -- %s -- %s\n" % (time.ctime(), reason, variant)
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue

                            else:
                                error = "If the following error message does not address the issue and the problem persists please contact admin: " + to_g
                                reason = "An error has occurred"
                                excep = "%s -- %s -- %s\n" % (time.ctime(), reason, variant)
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue

                        else:
                            # Insertions at exon boundaries are miss-handled by vm.g_to_t
                            if (
                                    obj.posedit.edit.type == 'ins' and obj.posedit.pos.start.offset == 0 and obj.posedit.pos.end.offset != 0) or (
                                    obj.posedit.edit.type == 'ins' and obj.posedit.pos.start.offset != 0 and obj.posedit.pos.end.offset == 0):
                                variant = str(obj)
                            else:
                                # Normalize was I believe to replace ref. Mapping does this anyway
                                # to_g = hn.normalize(to_g)
                                variant = str(va_func.myevm_g_to_t(evm, to_g, tx_ac))
                                tx_ac = ''

                    elif geno.search(input):
                        if plus.search(variant) or minus.search(variant):
                            to_g = va_func.genomic(variant, no_norm_evm, hp, hdp, primary_assembly, vm, hn, sf, nr_vm)
                            es = re.compile('error')
                            if es.search(str(to_g)):
                                if alt_aln_method != 'genebuild':
                                    error = "If the following error message does not address the issue and the problem persists please contact admin: " + to_g
                                    reason = "An error has occurred"
                                    excep = "%s -- %s -- %s\n" % (time.ctime(), reason, variant)
                                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                    logger.warning(str(error))
                                    continue

                                else:
                                    error = "If the following error message does not address the issue and the problem persists please contact admin: " + to_g
                                    reason = "An error has occurred"
                                    excep = "%s -- %s -- %s\n" % (time.ctime(), reason, variant)
                                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                    logger.warning(str(error))
                                    continue
                        else:
                            # Insertions at exon boundaries are miss-handled by vm.g_to_t
                            if (
                                    obj.posedit.edit.type == 'ins' and obj.posedit.pos.start.offset == 0 and obj.posedit.pos.end.offset != 0) or (
                                    obj.posedit.edit.type == 'ins' and obj.posedit.pos.start.offset != 0 and obj.posedit.pos.end.offset == 0):
                                variant = str(obj)
                            else:
                                # Normalize was I believe to replace ref. Mapping does this anyway
                                # to_g = hn.normalize(to_g)
                                variant = str(va_func.myevm_g_to_t(evm, to_g, tx_ac))
                                tx_ac = ''

                    else:
                        # Normalize the variant
                        error = 'false'
                        try:
                            h_variant = hn.normalize(obj)
                        except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                            error = str(e)
                            if re.match('Unsupported normalization of variants spanning the exon-intron boundary',
                                        error):
                                h_variant = obj
                                variant = variant
                                caution = 'This coding sequence variant description spans at least one intron'
                                automap = 'Use of the corresponding genomic sequence variant descriptions may be invalid. Please refer to https://www35.lamp.le.ac.uk/recommendations/'
                                validation['warnings'] = validation['warnings'] + ': ' + str(caution) + ': ' + str(
                                    automap)
                                logger.warning(str(caution) + ": " + str(automap))
                        else:
                            variant = str(h_variant)

                        tx_ac = ''
                        # Create a crosser (exon boundary crossed) variant
                        crossed_variant = str(evm._maybe_normalize(obj))
                        if variant == crossed_variant:
                            cross_variant = 'false'
                        else:
                            hgvs_crossed_variant = evm._maybe_normalize(obj)
                            cross_variant = [
                                "Coding sequence allowing for exon boundary crossing (default = no crossing)",
                                crossed_variant, hgvs_crossed_variant.ac]
                            cr_available = 'true'

                        # control of cross_variant
                        if boundary == 'false':
                            cross_variant = 'false'

                            error = va_func.validate(variant, hp=hp, vr=vr)
                            if error == 'false':
                                valid = 'true'
                            else:
                                excep = "%s -- %s -- %s\n" % (time.ctime(), error, variant)
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue

                    # Tackle the plus intronic offset
                    cck = 'false'
                    if (plus.search(input)):
                        # Regular expression catches the start of the interval only based on .00+00 pattern
                        inv_start = re.compile("\.\d+\+\d")
                        if (inv_start.search(input)):
                            # Find pattern e.g. +0000 and assign to a variable
                            off_value = re.search(r"(\+\d+)", input)
                            off_value = off_value.group(1)
                            # Integerise the value and assign to ex_offset
                            ex_offset = int(off_value)
                            cck = 'true'
                    if (minus.search(input)):
                        # Regular expression catches the start of the interval only based on .00-00 pattern
                        inv_start = re.compile("\.\d+\-\d")
                        if (inv_start.search(input)):
                            # Find pattern e.g. -0000 and assign to a variable
                            off_value = re.search(r"(\-\d+)", input)
                            off_value = off_value.group(1)
                            # Integerise the value and assign to ex_offset
                            ex_offset = int(off_value)
                            cck = 'true'

                    # COORDINATE CHECKER
                    # hgvs will handle incorrect coordinates so need to automap errors
                    # Make sure any input intronic coordinates are correct
                    # Get the desired transcript
                    pat_r = re.compile(':r.')
                    pat_g = re.compile(':g.')
                    if cck == 'true':
                        dl = re.compile('del')
                        # This should only ever hit coding and RNA variants
                        if dl.search(variant):
                            # RNA
                            if pat_r.search(trapped_input):

                                coding = va_func.coding(variant, hp)
                                trans_acc = coding.ac
                                # c to Genome coordinates - Map the variant to the genome
                                pre_var = va_func.genomic(variant, no_norm_evm, hp, hdp, primary_assembly, vm, hn, sf,
                                                          nr_vm)
                                # genome back to C coordinates
                                post_var = va_func.myevm_g_to_t(evm, pre_var, trans_acc)

                                test = hp.parse_hgvs_variant(input)
                                if post_var.posedit.pos.start.base != test.posedit.pos.start.base or post_var.posedit.pos.end.base != test.posedit.pos.end.base:
                                    caution = 'The entered coordinates do not agree with the intron/exon boundaries for the selected transcript:'
                                    automap = 'Automap has corrected the coordinates to match the intron/exon boundaries for the selected transcript'
                                    # automapping of variant completed
                                    # Change to rna variant
                                    posedit = query.posedit
                                    posedit = posedit.lower()
                                    query.posedit = posedit
                                    query.type = 'r'
                                    post_var = str(query)
                                    automap = trapped_input + ' automapped to ' + str(post_var)
                                    validation['warnings'] = validation['warnings'] + ': ' + str(caution) + ': ' + str(
                                        automap)
                                    relevant = "Select the automapped transcript and click Submit to analyse"
                                    rel_var = []
                                    rel_var.append(post_var)
                                    # Add gene symbols to the link
                                    cp_rel = copy.copy(rel_var)
                                    del rel_var[:]
                                    for accessions in cp_rel:
                                        error = 'false'
                                        hgvs_vt = hp.parse_hgvs_variant(str(accessions))
                                        try:
                                            tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
                                        except hgvs.exceptions.HGVSError as e:
                                            error = str(e)
                                        if error != 'false':
                                            accessions = ['', str(hgvs_vt)]
                                            rel_var.append(accessions)

                                        else:
                                            # Get hgnc Gene name from command
                                            data = va_func.hgnc_rest(path="/search/prev_symbol/" + tx_id_info[6])
                                            if data['error'] != 'false':
                                                error = data['error']
                                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                                logger.warning(str(error))
                                                continue
                                            else:
                                                # Set the hgnc name correctly
                                                # If the name is correct no record will be found
                                                if int(data['record']['response']['numFound']) == 0:
                                                    current = tx_id_info[6]
                                                else:
                                                    current = data['record']['response']['docs'][0]['symbol']
                                            accessions = [str(current), str(hgvs_vt)]
                                            rel_var.append(accessions)
                                            # Add gene symbols to the link
                                            cp_rel = copy.copy(rel_var)
                                            del rel_var[:]
                                            for accessions in cp_rel:
                                                error = 'false'
                                                hgvs_vt = hp.parse_hgvs_variant(str(accessions[1]))
                                                try:
                                                    tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
                                                except hgvs.exceptions.HGVSError as e:
                                                    error = str(e)
                                                if error != 'false':
                                                    accessions = ['', str(hgvs_vt)]
                                                    rel_var.append(accessions)
                                                else:
                                                    # Get hgnc Gene name from command
                                                    data = va_func.hgnc_rest(
                                                        path="/search/prev_symbol/" + tx_id_info[6])
                                                    if data['error'] != 'false':
                                                        error = data['error']
                                                        validation['warnings'] = validation['warnings'] + ': ' + str(
                                                            error)
                                                        logger.warning(str(error))
                                                        continue
                                                    else:
                                                        # Set the hgnc name correctly
                                                        # If the name is correct no record will be found
                                                        if int(data['record']['response']['numFound']) == 0:
                                                            current = tx_id_info[6]
                                                        else:
                                                            current = data['record']['response']['docs'][0]['symbol']
                                                    accessions = [str(current), str(hgvs_vt)]
                                                    rel_var.append(accessions)
                                    # Kill current line and append for re-submission
                                    # Tag the line so that it is not written out
                                    validation['write'] = 'false'
                                    # Set the values and append to batch_list
                                    query = {'quibble': valstr(hgvs_vt), 'id': validation['id'], 'warnings': automap,
                                             'description': '', 'coding': '', 'coding_g': '', 'genomic_r': '',
                                             'genomic_g': '', 'protein': '', 'write': 'true',
                                             'primary_assembly': primary_assembly, 'order': ordering}
                                    batch_list.append(query)

                            # Coding
                            else:
                                coding = va_func.coding(variant, hp)
                                trans_acc = coding.ac
                                # c to Genome coordinates - Map the variant to the genome
                                pre_var = hp.parse_hgvs_variant(variant)
                                try:
                                    pre_var = va_func.myevm_t_to_g(pre_var, hdp, no_norm_evm, primary_assembly, vm, hp,
                                                                   hn, sf, nr_vm)
                                except:
                                    e = sys.exc_info()[1]
                                    error = str(e)
                                    reason = 'Input coordinates may be invalid'
                                    if error == 'expected from_start_i <= from_end_i':
                                        error = 'Automap is unable to correct the input exon/intron boundary coordinates, please check your variant description'
                                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                        continue
                                    else:
                                        exceptPass()
                                else:
                                    exceptPass()
                                # genome back to C coordinates
                                try:
                                    post_var = va_func.myevm_g_to_t(evm, pre_var, trans_acc)
                                except hgvs.exceptions.HGVSError as error:
                                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                    logger.warning(str(error))
                                    continue
                                query = post_var
                                test = hp.parse_hgvs_variant(input)
                                if post_var.posedit.pos.start.base != test.posedit.pos.start.base or post_var.posedit.pos.end.base != test.posedit.pos.end.base:
                                    caution = 'The entered coordinates do not agree with the intron/exon boundaries for the selected transcript:'
                                    automap = 'Automap has corrected the coordinates to match the intron/exon boundaries for the selected transcript'
                                    # automapping of variant completed
                                    automap = trapped_input + ' automapped to ' + str(post_var)
                                    validation['warnings'] = str(validation['warnings']) + str(caution) + ': ' + str(
                                        automap)
                                    relevant = "Select the automapped transcript and click Submit to analyse"
                                    rel_var = []
                                    rel_var.append(post_var)
                                    # Add gene symbols to the link
                                    cp_rel = copy.copy(rel_var)
                                    del rel_var[:]
                                    for accessions in cp_rel:
                                        error = 'false'
                                        hgvs_vt = hp.parse_hgvs_variant(str(accessions))
                                        try:
                                            tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
                                        except hgvs.exceptions.HGVSError as e:
                                            error = str(e)
                                        if error != 'false':
                                            accessions = ['', str(hgvs_vt)]
                                            rel_var.append(accessions)
                                        else:
                                            # Get hgnc Gene name from command
                                            data = va_func.hgnc_rest(path="/search/prev_symbol/" + tx_id_info[6])
                                            if data['error'] != 'false':
                                                error = data['error']
                                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                                logger.warning(str(error))
                                                continue

                                            else:
                                                # Set the hgnc name correctly
                                                # If the name is correct no record will be found
                                                if int(data['record']['response']['numFound']) == 0:
                                                    current = tx_id_info[6]
                                                else:
                                                    current = data['record']['response']['docs'][0]['symbol']
                                            accessions = [str(current), str(hgvs_vt)]
                                            rel_var.append(accessions)
                                            # Add gene symbols to the link
                                            cp_rel = copy.copy(rel_var)
                                            del rel_var[:]
                                            for accessions in cp_rel:
                                                error = 'false'
                                                hgvs_vt = hp.parse_hgvs_variant(str(accessions[1]))
                                                try:
                                                    tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
                                                except hgvs.exceptions.HGVSError as e:
                                                    error = str(e)
                                                if error != 'false':
                                                    accessions = ['', str(hgvs_vt)]
                                                    rel_var.append(accessions)
                                                else:
                                                    # Get hgnc Gene name from command
                                                    data = va_func.hgnc_rest(
                                                        path="/search/prev_symbol/" + tx_id_info[6])
                                                    if data['error'] != 'false':
                                                        error = data['error']
                                                        validation['warnings'] = validation['warnings'] + ': ' + str(
                                                            error)
                                                        logger.warning(str(error))
                                                        continue

                                                    else:
                                                        # Set the hgnc name correctly
                                                        # If the name is correct no record will be found
                                                        if int(data['record']['response']['numFound']) == 0:
                                                            current = tx_id_info[6]
                                                        else:
                                                            current = data['record']['response']['docs'][0]['symbol']
                                                    accessions = [str(current), str(hgvs_vt)]
                                                    rel_var.append(accessions)
                                    # Kill current line and append for re-submission
                                    # Tag the line so that it is not written out
                                    validation['write'] = 'false'
                                    # Set the values and append to batch_list
                                    query = {'quibble': valstr(hgvs_vt), 'id': validation['id'], 'warnings': automap,
                                             'description': '', 'coding': '', 'coding_g': '', 'genomic_r': '',
                                             'genomic_g': '', 'protein': '', 'write': 'true',
                                             'primary_assembly': primary_assembly, 'order': ordering}
                                    batch_list.append(query)

                        else:
                            if pat_r.search(trapped_input):
                                coding = va_func.coding(variant, hp)
                                trans_acc = coding.ac
                                # c to Genome coordinates - Map the variant to the genome
                                pre_var = va_func.genomic(variant, no_norm_evm, hp, hdp, primary_assembly, vm, hn, sf,
                                                          nr_vm)
                                # genome back to C coordinates
                                post_var = va_func.myevm_g_to_t(evm, pre_var, trans_acc)

                                test = hp.parse_hgvs_variant(input)
                                if post_var.posedit.pos.start.base != test.posedit.pos.start.base or post_var.posedit.pos.end.base != test.posedit.pos.end.base:
                                    caution = 'The entered coordinates do not agree with the intron/exon boundaries for the selected transcript:'
                                    automap = 'Automap has corrected the coordinates to match the intron/exon boundaries for the selected transcript'
                                    # automapping of variant completed
                                    # Change to rna variant
                                    posedit = query.posedit
                                    posedit = posedit.lower()
                                    query.posedit = posedit
                                    query.type = 'r'
                                    post_var = str(query)
                                    automap = input + ' automapped to ' + post_var
                                    validation['warnings'] = validation['warnings'] + ': ' + str(caution) + ': ' + str(
                                        automap)
                                    relevant = "Select the automapped transcript and click Submit to analyse"
                                    rel_var = []
                                    rel_var.append(post_var)
                                    # Add gene symbols to the link
                                    cp_rel = copy.copy(rel_var)
                                    del rel_var[:]
                                    for accessions in cp_rel:
                                        error = 'false'
                                        hgvs_vt = hp.parse_hgvs_variant(str(accessions))
                                        try:
                                            tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
                                        except hgvs.exceptions.HGVSError as e:
                                            error = str(e)
                                        if error != 'false':
                                            accessions = ['', str(hgvs_vt)]
                                            rel_var.append(accessions)
                                        else:
                                            # Get hgnc Gene name from command
                                            data = va_func.hgnc_rest(path="/search/prev_symbol/" + tx_id_info[6])
                                            if data['error'] != 'false':
                                                error = data['error']
                                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                                logger.warning(str(error))
                                                continue

                                            else:
                                                # Set the hgnc name correctly
                                                # If the name is correct no record will be found
                                                if int(data['record']['response']['numFound']) == 0:
                                                    current = tx_id_info[6]
                                                else:
                                                    current = data['record']['response']['docs'][0]['symbol']
                                            accessions = [str(current), str(hgvs_vt)]
                                            rel_var.append(accessions)
                                    # Kill current line and append for re-submission
                                    # Tag the line so that it is not written out
                                    validation['write'] = 'false'
                                    # Set the values and append to batch_list
                                    query = {'quibble': valstr(hgvs_vt), 'id': validation['id'], 'warnings': automap,
                                             'description': '', 'coding': '', 'coding_g': '', 'genomic_r': '',
                                             'genomic_g': '', 'protein': '', 'write': 'true',
                                             'primary_assembly': primary_assembly, 'order': ordering}
                                    batch_list.append(query)

                            else:
                                coding = va_func.coding(variant, hp)
                                trans_acc = coding.ac
                                # c to Genome coordinates - Map the variant to the genome
                                pre_var = va_func.genomic(variant, no_norm_evm, hp, hdp, primary_assembly, vm, hn, sf,
                                                          nr_vm)
                                # genome back to C coordinates
                                post_var = va_func.myevm_g_to_t(evm, pre_var, trans_acc)

                                test = hp.parse_hgvs_variant(input)
                                if post_var.posedit.pos.start.base != test.posedit.pos.start.base or post_var.posedit.pos.end.base != test.posedit.pos.end.base:
                                    caution = 'The entered coordinates do not agree with the intron/exon boundaries for the selected transcript:'
                                    automap = 'Automap has corrected the coordinates to match the intron/exon boundaries for the selected transcript'
                                    # automapping of variant completed
                                    automap = str(trapped_input) + ' automapped to ' + str(post_var)
                                    validation['warnings'] = validation['warnings'] + ': ' + str(caution) + ': ' + str(
                                        automap)
                                    relevant = "Select the automapped transcript and click Submit to analyse"
                                    rel_var = []
                                    rel_var.append(post_var)
                                    # Add gene symbols to the link
                                    cp_rel = copy.copy(rel_var)
                                    del rel_var[:]
                                    for accessions in cp_rel:
                                        error = 'false'
                                        hgvs_vt = hp.parse_hgvs_variant(str(accessions))
                                        try:
                                            tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
                                        except hgvs.exceptions.HGVSError as e:
                                            error = str(e)
                                        if error != 'false':
                                            accessions = ['', str(hgvs_vt)]
                                            rel_var.append(accessions)
                                        else:
                                            # Get hgnc Gene name from command
                                            data = va_func.hgnc_rest(path="/search/prev_symbol/" + tx_id_info[6])
                                            if data['error'] != 'false':
                                                reason = 'Cannot currently display the required information:'
                                                error = data['error']
                                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                                logger.warning(str(error))
                                                continue

                                            else:
                                                # Set the hgnc name correctly
                                                # If the name is correct no record will be found
                                                if int(data['record']['response']['numFound']) == 0:
                                                    current = tx_id_info[6]
                                                else:
                                                    current = data['record']['response']['docs'][0]['symbol']
                                            accessions = [str(current), str(hgvs_vt)]
                                            rel_var.append(accessions)
                                    # Kill current line and append for re-submission
                                    # Tag the line so that it is not written out
                                    validation['write'] = 'false'
                                    # Set the values and append to batch_list
                                    query = {'quibble': valstr(hgvs_vt), 'id': validation['id'], 'warnings': automap,
                                             'description': '', 'coding': '', 'coding_g': '', 'genomic_r': '',
                                             'genomic_g': '', 'protein': '', 'write': 'true',
                                             'primary_assembly': primary_assembly, 'order': ordering}
                                    batch_list.append(query)


                    # If cck not true
                    elif pat_r.search(trapped_input):
                        # set input hgvs object
                        hgvs_rna_input = hp.parse_hgvs_variant(
                            trapped_input)  # Traps the hgvs variant of r. for further use
                        inp = str(va_func.hgvs_r_to_c(hgvs_rna_input))
                        # Regex
                        plus = re.compile("\d\+\d")  # finds digit + digit
                        minus = re.compile("\d\-\d")  # finds digit - digit
                        if plus.search(input) or minus.search(input):
                            to_g = va_func.genomic(inp, no_norm_evm, hp, hdp, primary_assembly, vm, hn, sf, nr_vm)
                            es = re.compile('error')
                            if es.search(str(to_g)):
                                if alt_aln_method != 'genebuild':
                                    error = "If the following error message does not address the issue and the problem persists please contact admin: " + to_g
                                    reason = "An error has occurred"
                                    excep = "%s -- %s -- %s\n" % (time.ctime(), reason, variant)
                                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                    logger.warning(str(error))
                                    continue

                                else:
                                    error = "If the following error message does not address the issue and the problem persists please contact admin: " + to_g
                                    reason = "An error has occurred"
                                    excep = "%s -- %s -- %s\n" % (time.ctime(), reason, variant)
                                    validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                    logger.warning(str(error))
                                    continue

                            else:
                                # Set variants pre and post genomic norm
                                hgvs_inp = va_func.myevm_g_to_t(evm, to_g, tx_ac=obj.ac)
                                to_g = hn.normalize(to_g)
                                hgvs_otp = va_func.myevm_g_to_t(evm, to_g, tx_ac=obj.ac)
                                tx_ac = ''
                        else:
                            # Set variants pre and post RNA norm
                            hgvs_inp = hp.parse_hgvs_variant(inp)
                            try:
                                hgvs_otp = hn.normalize(hgvs_inp)
                            except hgvs.exceptions.HGVSError as e:
                                hgvs_otp = hgvs_inp
                            tx_ac = ''

                        # Set remaining variables
                        redit = str(hgvs_otp.posedit.edit)
                        redit = redit.lower()
                        hgvs_otp.posedit.edit = redit
                        otp = str(hgvs_otp)
                        query = str(hgvs_otp.posedit.pos)
                        test = str(hgvs_inp.posedit.pos)
                        query = query.replace('T', 'U')
                        query = query.replace('ENSU', 'ENST')
                        test = test.replace('T', 'U')
                        test = test.replace('ENSU', 'ENST')
                        output = otp.replace(':c.', ':r.')
                        # Apply coordinates test
                        if query != test:
                            caution = 'The variant description ' + input + ' requires alteration to comply with HGVS variant nomenclature:'
                            automap = 'Automap has corrected the variant description'
                            # automapping of variant completed
                            automap = trapped_input + ' automapped to ' + output
                            validation['warnings'] = validation['warnings'] + ': ' + str(caution) + ': ' + str(automap)
                            relevant = "Select the automapped transcript and click Submit to analyse"
                            rel_var = []
                            rel_var.append(output)
                            # Add gene symbols to the link
                            cp_rel = copy.copy(rel_var)
                            del rel_var[:]
                            for accessions in cp_rel:
                                error = 'false'
                                hgvs_vt = hp.parse_hgvs_variant(str(accessions))
                                try:
                                    tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
                                except hgvs.exceptions.HGVSError as e:
                                    error = str(e)
                                if error != 'false':
                                    accessions = ['', str(hgvs_vt)]
                                    rel_var.append(accessions)
                                else:
                                    # Get hgnc Gene name from command
                                    data = va_func.hgnc_rest(path="/search/prev_symbol/" + tx_id_info[6])
                                    if data['error'] != 'false':
                                        error = data['error']
                                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                        logger.warning(str(error))
                                        continue

                                    else:
                                        # Set the hgnc name correctly
                                        # If the name is correct no record will be found
                                        if int(data['record']['response']['numFound']) == 0:
                                            current = tx_id_info[6]
                                        else:
                                            current = data['record']['response']['docs'][0]['symbol']
                                    accessions = [str(current), str(hgvs_vt)]
                                    rel_var.append(accessions)
                            # Kill current line and append for re-submission
                            # Tag the line so that it is not written out
                            validation['write'] = 'false'
                            # Set the values and append to batch_list
                            query = {'quibble': valstr(hgvs_vt), 'id': validation['id'], 'warnings': automap,
                                     'description': '', 'coding': '', 'coding_g': '', 'genomic_r': '', 'genomic_g': '',
                                     'protein': '', 'write': 'true', 'primary_assembly': primary_assembly,
                                     'order': ordering}
                            batch_list.append(query)

                    elif pat_g.search(input):
                        pass

                    else:
                        query = hp.parse_hgvs_variant(variant)
                        test = hp.parse_hgvs_variant(input)
                        if query.posedit.pos != test.posedit.pos:
                            caution = 'The variant description ' + input + ' requires alteration to comply with HGVS variant nomenclature:'
                            automap = 'Automap has corrected the variant description'
                            # automapping of variant completed
                            automap = str(test) + ' automapped to ' + str(query)
                            validation['warnings'] = validation['warnings'] + ': ' + str(caution) + ': ' + str(automap)
                            relevant = "Select the automapped transcript and click Submit to analyse"
                            rel_var = []
                            rel_var.append(query)
                            # Add gene symbols to the link
                            cp_rel = copy.copy(rel_var)
                            del rel_var[:]
                            for accessions in cp_rel:
                                error = 'false'
                                hgvs_vt = hp.parse_hgvs_variant(str(accessions))
                                try:
                                    tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
                                except hgvs.exceptions.HGVSError as e:
                                    error = str(e)
                                if error != 'false':
                                    accessions = ['', str(hgvs_vt)]
                                    rel_var.append(accessions)
                                else:
                                    # Get hgnc Gene name from command
                                    data = va_func.hgnc_rest(path="/search/prev_symbol/" + tx_id_info[6])
                                    if data['error'] != 'false':
                                        reason = 'Cannot currently display the required information:'
                                        error = data['error']
                                        validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                        logger.warning(str(error))
                                        continue

                                    else:
                                        # Set the hgnc name correctly
                                        # If the name is correct no record will be found
                                        if int(data['record']['response']['numFound']) == 0:
                                            current = tx_id_info[6]
                                        else:
                                            current = data['record']['response']['docs'][0]['symbol']
                                    accessions = [str(current), str(hgvs_vt)]
                                    rel_var.append(accessions)
                            # Kill current line and append for re-submission
                            # Tag the line so that it is not written out
                            validation['write'] = 'false'
                            # Set the values and append to batch_list
                            query = {'quibble': valstr(hgvs_vt), 'id': validation['id'], 'warnings': automap,
                                     'description': '', 'coding': '', 'coding_g': '', 'genomic_r': '', 'genomic_g': '',
                                     'protein': '', 'write': 'true', 'primary_assembly': primary_assembly,
                                     'order': ordering}
                            batch_list.append(query)

                    # VALIDATION of intronic variants
                    pre_valid = hp.parse_hgvs_variant(input)
                    post_valid = hp.parse_hgvs_variant(variant)
                    if valid == 'false':
                        error = 'false'
                        genomic_validation = str(
                            va_func.genomic(input, no_norm_evm, hp, hdp, primary_assembly, vm, hn, sf, nr_vm))
                        del_end = re.compile('\ddel$')
                        delins = re.compile('delins')
                        inv = re.compile('inv')
                        if valstr(pre_valid) != valstr(post_valid):
                            if type != ':g.':
                                if caution == '':
                                    caution = valstr(pre_valid) + ' automapped to ' + valstr(post_valid)
                                else:
                                    pass
                                validation['warnings'] = validation['warnings'] + ': ' + str(caution)
                                logger.warning(str(caution))
                            else:
                                pass
                        else:
                            pass

                        # Apply validation to intronic variant descriptions (should be valid but make sure)
                        error = va_func.validate(genomic_validation, hp=hp, vr=vr)
                        if error == 'false':
                            valid = 'true'
                        else:

                            excep = "%s -- %s -- %s\n" % (time.ctime(), error, variant)
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            continue

                    if valid == 'true':
                        var_tab = 'true'
                        cores = "HGVS-compliant variant descriptions" + warning

                        # v0.1a1 edit
                        if valstr(pre_valid) != valstr(post_valid):
                            if type == ':g.':
                                if caution == '':
                                    caution = valstr(pre_valid) + ' automapped to ' + valstr(post_valid)
                                else:
                                    pass
                                validation['warnings'] = validation['warnings'] + ': ' + str(caution)
                            else:
                                pass
                        else:
                            pass

                        # COLLECT VARIANT DESCRIPTIONS
                        ##############################

                        # Coding sequence - BASED ON NORMALIZED VARIANT IF EXONIC
                        hgvs_coding = va_func.coding(variant, hp)
                        boundary = re.compile('exon-intron boundary')
                        spanning = re.compile('exon/intron')

                        try:
                            hgvs_coding = hn.normalize(hgvs_coding)
                        except hgvs.exceptions.HGVSError as e:
                            error = str(e)

                        # Gap compensating code status
                        gap_compensation = True

                        # Gap gene black list
                        try:
                            gene_symbol = va_dbCrl.data.get_gene_symbol_from_transcriptID(hgvs_coding.ac)
                        except Exception:
                            exceptPass()
                        else:
                            # If the gene symbol is not in the list, the value False will be returned
                            gap_compensation = gapGenes.gap_black_list(gene_symbol)

                        # Intron spanning variants
                        if re.search('boundary', str(error)) or re.search('spanning', str(error)):
                            try:
                                hgvs_coding = evm._maybe_normalize(hgvs_coding)
                                gap_compensation = False
                            except hgvs.exceptions.HGVSError as error:
                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                logger.warning(str(error))
                                continue
                        else:
                            pass

                        # Warn status
                        logger.warning("gap_compensation_1 = " + str(gap_compensation))
                        coding = valstr(hgvs_coding)

                        # RNA sequence
                        hgvs_rna = copy.deepcopy(hgvs_coding)
                        hgvs_rna = va_func.hgvs_c_to_r(hgvs_rna)
                        rna = str(hgvs_rna)

                        # Genomic sequence
                        hgvs_genomic = va_func.myevm_t_to_g(hgvs_coding, hdp, no_norm_evm, primary_assembly, vm, hp, hn,
                                                            sf, nr_vm)
                        final_hgvs_genomic = hgvs_genomic

                        # genomic_possibilities
                        # 1. take the simple 3 pr normalized hgvs_genomic
                        # 2. Lock in hgvs_genomic at its most 5 prime position wrt genome
                        hgvs_genomic_possibilities = []

                        # Loop out gap finding code under these circumstances!
                        if gap_compensation is True:
                            logger.warning('g_to_t gap code 1 active')
                            rn_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic)
                            hgvs_genomic_possibilities.append(rn_hgvs_genomic)
                            if orientation != -1:
                                try:
                                    chromosome_normalized_hgvs_coding = reverse_normalizer.normalize(hgvs_coding)
                                except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                                    error = str(e)
                                    chromosome_normalized_hgvs_coding = hgvs_coding
                            else:
                                try:
                                    chromosome_normalized_hgvs_coding = hn.normalize(hgvs_coding)
                                except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                                    error = str(e)
                                    chromosome_normalized_hgvs_coding = hgvs_coding

                            most_3pr_hgvs_genomic = va_func.myvm_t_to_g(chromosome_normalized_hgvs_coding,
                                                                        hgvs_genomic.ac, no_norm_evm, vm, hp, hn, sf,
                                                                        nr_vm)
                            hgvs_genomic_possibilities.append(most_3pr_hgvs_genomic)

                            # Push from side to side to try pick up odd placements
                            # MAKE A NO NORM HGVS2VCF
                            # First to the right
                            hgvs_stash = copy.deepcopy(hgvs_coding)
                            try:
                                hgvs_stash = no_norm_evm.c_to_n(hgvs_stash)
                            except:
                                exceptPass()
                            try:
                                stash_ac = hgvs_stash.ac
                                stash_dict = va_H2V.hard_right_hgvs2vcf(hgvs_stash, primary_assembly, hn, sf)
                                stash_pos = int(stash_dict['pos'])
                                stash_ref = stash_dict['ref']
                                stash_alt = stash_dict['alt']
                                # Generate an end position
                                stash_end = str(stash_pos + len(stash_ref) - 1)
                                # make a not real deletion insertion
                                stash_hgvs_not_delins = hp.parse_hgvs_variant(
                                    stash_ac + ':' + hgvs_stash.type + '.' + str(
                                        stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)
                                try:
                                    stash_hgvs_not_delins = no_norm_evm.n_to_c(stash_hgvs_not_delins)
                                except:
                                    exceptPass()
                                    # Store a tx copy for later use
                                test_stash_tx_right = copy.deepcopy(stash_hgvs_not_delins)
                                # stash_genomic = vm.t_to_g(test_stash_tx_right, hgvs_genomic.ac)
                                stash_genomic = va_func.myvm_t_to_g(test_stash_tx_right, hgvs_genomic.ac, no_norm_evm,
                                                                    vm, hp, hn, sf, nr_vm)
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
                                    hgvs_reform_ident = hp.parse_hgvs_variant(reform_ident)
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
                                test_stash_tx_right = copy.deepcopy(hgvs_coding)
                                exceptPass()
                            # Intronic positions not supported. Will cause a Value Error
                            except ValueError:
                                test_stash_tx_right = copy.deepcopy(hgvs_coding)
                                exceptPass()

                            # Then to the left
                            hgvs_stash = copy.deepcopy(hgvs_coding)
                            try:
                                hgvs_stash = no_norm_evm.c_to_n(hgvs_stash)
                            except:
                                exceptPass()
                            try:
                                stash_ac = hgvs_stash.ac
                                stash_dict = va_H2V.hard_left_hgvs2vcf(hgvs_stash, primary_assembly, reverse_normalizer,
                                                                       sf)
                                stash_pos = int(stash_dict['pos'])
                                stash_ref = stash_dict['ref']
                                stash_alt = stash_dict['alt']
                                # Generate an end position
                                stash_end = str(stash_pos + len(stash_ref) - 1)
                                # make a not real deletion insertion
                                stash_hgvs_not_delins = hp.parse_hgvs_variant(
                                    stash_ac + ':' + hgvs_stash.type + '.' + str(
                                        stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)
                                try:
                                    stash_hgvs_not_delins = no_norm_evm.n_to_c(stash_hgvs_not_delins)
                                except:
                                    exceptPass()
                                    # Store a tx copy for later use
                                test_stash_tx_left = copy.deepcopy(stash_hgvs_not_delins)
                                # stash_genomic = vm.t_to_g(test_stash_tx_left, hgvs_genomic.ac)
                                stash_genomic = va_func.myvm_t_to_g(test_stash_tx_left, hgvs_genomic.ac, no_norm_evm,
                                                                    vm, hp, hn, sf, nr_vm)
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
                                    hgvs_reform_ident = hp.parse_hgvs_variant(reform_ident)
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
                                test_stash_tx_left = copy.deepcopy(hgvs_coding)
                                exceptPass()
                            except ValueError:
                                test_stash_tx_left = copy.deepcopy(hgvs_coding)
                                exceptPass()

                            # direct mapping from reverse_normalized transcript insertions in the delins format
                            try:
                                if hgvs_coding.posedit.edit.type == 'ins':
                                    most_5pr_hgvs_transcript_variant = copy.deepcopy(hgvs_coding)
                                    most_3pr_hgvs_transcript_variant = reverse_normalizer.normalize(hgvs_coding)
                                    try:
                                        n_3pr = vm.c_to_n(most_3pr_hgvs_transcript_variant)
                                        n_5pr = vm.c_to_n(most_5pr_hgvs_transcript_variant)
                                    except:
                                        n_3pr = most_3pr_hgvs_transcript_variant
                                        n_5pr = most_5pr_hgvs_transcript_variant
                                    # Make into a delins by adding the ref bases to the variant ref and alt
                                    pr3_ref = sf.fetch_seq(hgvs_coding.ac, n_3pr.posedit.pos.start.base - 1,
                                                           n_3pr.posedit.pos.end.base)
                                    pr5_ref = sf.fetch_seq(hgvs_coding.ac, n_5pr.posedit.pos.start.base - 1,
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
                                    genomic_from_most_3pr_hgvs_transcript_variant = vm.t_to_g(
                                        most_3pr_hgvs_transcript_variant, hgvs_genomic.ac)
                                    genomic_from_most_5pr_hgvs_transcript_variant = vm.t_to_g(
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
                                            genomic_from_most_3pr_hgvs_transcript_variant = hp.parse_hgvs_variant(
                                                genomic_from_most_3pr_hgvs_transcript_variant_delins_from_dup)

                                    try:
                                        if most_3pr_hgvs_transcript_variant.posedit.edit.alt is None:
                                            most_3pr_hgvs_transcript_variant.posedit.edit.alt = ''
                                    except Exception as e:
                                        if str(e) == "'Dup' object has no attribute 'alt'":
                                            most_3pr_hgvs_transcript_variant_delins_from_dup = most_3pr_hgvs_transcript_variant.ac + ':' + most_3pr_hgvs_transcript_variant.type + '.' + str(
                                                most_3pr_hgvs_transcript_variant.posedit.pos.start.base) + '_' + str(
                                                most_3pr_hgvs_transcript_variant.posedit.pos.end.base) + 'del' + most_3pr_hgvs_transcript_variant.posedit.edit.ref + 'ins' + most_3pr_hgvs_transcript_variant.posedit.edit.ref + most_3pr_hgvs_transcript_variant.posedit.edit.ref
                                            most_3pr_hgvs_transcript_variant = hp.parse_hgvs_variant(
                                                most_3pr_hgvs_transcript_variant_delins_from_dup)

                                    try:
                                        if genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.alt is None:
                                            genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.alt = ''
                                    except Exception as e:
                                        if str(e) == "'Dup' object has no attribute 'alt'":
                                            genomic_from_most_5pr_hgvs_transcript_variant_delins_from_dup = genomic_from_most_5pr_hgvs_transcript_variant.ac + ':' + genomic_from_most_5pr_hgvs_transcript_variant.type + '.' + str(
                                                genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.start.base) + '_' + str(
                                                genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.end.base) + 'del' + genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.ref + 'ins' + genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.ref + genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.ref
                                            genomic_from_most_5pr_hgvs_transcript_variant = hp.parse_hgvs_variant(
                                                genomic_from_most_5pr_hgvs_transcript_variant_delins_from_dup)

                                    try:
                                        if most_5pr_hgvs_transcript_variant.posedit.edit.alt is None:
                                            most_5pr_hgvs_transcript_variant.posedit.edit.alt = ''
                                    except Exception as e:
                                        if str(e) == "'Dup' object has no attribute 'alt'":
                                            most_5pr_hgvs_transcript_variant_delins_from_dup = most_5pr_hgvs_transcript_variant.ac + ':' + most_5pr_hgvs_transcript_variant.type + '.' + str(
                                                most_5pr_hgvs_transcript_variant.posedit.pos.start.base) + '_' + str(
                                                most_5pr_hgvs_transcript_variant.posedit.pos.end.base) + 'del' + most_5pr_hgvs_transcript_variant.posedit.edit.ref + 'ins' + most_5pr_hgvs_transcript_variant.posedit.edit.ref + most_5pr_hgvs_transcript_variant.posedit.edit.ref
                                            most_5pr_hgvs_transcript_variant = hp.parse_hgvs_variant(
                                                most_5pr_hgvs_transcript_variant_delins_from_dup)

                                    if len(genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.alt) < len(
                                            most_3pr_hgvs_transcript_variant.posedit.edit.alt):
                                        hgvs_genomic_possibilities.append(genomic_from_most_3pr_hgvs_transcript_variant)
                                    if len(genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.alt) < len(
                                            most_5pr_hgvs_transcript_variant.posedit.edit.alt):
                                        hgvs_genomic_possibilities.append(genomic_from_most_5pr_hgvs_transcript_variant)

                            except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                                error = str(e)
                                if re.match('Normalization of intronic variants is not supported', error) or re.match(
                                        'Unsupported normalization of variants spanning the exon-intron boundary',
                                        error):
                                    pass

                            logger.info('\nGENOMIC POSSIBILITIES')
                            for possibility in hgvs_genomic_possibilities:
                                if possibility == '':
                                    logger.info('X')
                                else:
                                    logger.info(valstr(possibility))

                            logger.info('\n')

                            # Set variables for problem specific warnings
                            gapped_alignment_warning = ''
                            corrective_action_taken = ''
                            gapped_transcripts = ''
                            auto_info = ''

                            # Mark as not disparity detected
                            disparity_deletion_in = ['false', 'false']

                            # Loop through to see if a gap can be located
                            # Set the variables required for corrective normalization
                            possibility_counter = 0
                            suppress_c_normalization = 'false'  # Applies to boundary crossing normalization

                            # Copy a version of hgvs_genomic_possibilities
                            for possibility in hgvs_genomic_possibilities:
                                possibility_counter = possibility_counter + 1

                                # Loop out stash possibilities which will not spot gaps so are empty
                                if possibility == '':
                                    continue

                                # Use VCF generation code to push hgvs_genomic as for 5 prime as possible to uncover gaps
                                hgvs_genomic_variant = copy.deepcopy(possibility)
                                stored_hgvs_genomic_variant = copy.deepcopy(hgvs_genomic_variant)

                                # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
                                try:
                                    reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic_variant)
                                except hgvs.exceptions.HGVSError as e:
                                    # Strange error caused by gap in genomic
                                    error = str(e)
                                    if re.search('base start position must be <= end position', error):
                                        if hgvs_genomic.posedit.edit.type == 'delins':
                                            start = hgvs_genomic.posedit.pos.start.base
                                            end = hgvs_genomic.posedit.pos.end.base
                                            lhb = sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                                            rhb = sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                                            hgvs_genomic.posedit.edit.ref = lhb + rhb
                                            hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                                            hgvs_genomic.posedit.pos.start.base = end
                                            hgvs_genomic.posedit.pos.end.base = start
                                            reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic)
                                        if hgvs_genomic.posedit.edit.type == 'del':
                                            start = hgvs_genomic.posedit.pos.start.base
                                            end = hgvs_genomic.posedit.pos.end.base
                                            lhb = sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                                            rhb = sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                                            hgvs_genomic.posedit.edit.ref = lhb + rhb
                                            hgvs_genomic.posedit.edit.alt = lhb + rhb
                                            hgvs_genomic.posedit.pos.start.base = end
                                            hgvs_genomic.posedit.pos.end.base = start
                                            reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic)
                                    if re.search('insertion length must be 1', error):
                                        if hgvs_genomic.posedit.edit.type == 'ins':
                                            start = hgvs_genomic.posedit.pos.start.base
                                            end = hgvs_genomic.posedit.pos.end.base
                                            ref_bases = sf.fetch_seq(str(hgvs_genomic.ac), start - 1, end)
                                            lhb = sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                                            rhb = sf.fetch_seq(str(hgvs_genomic.ac), start, end)
                                            hgvs_genomic.posedit.edit.ref = lhb + rhb
                                            hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                                            reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic)

                                hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)
                                # Store a copy for later use
                                stored_hgvs_genomic_5pr = copy.deepcopy(hgvs_genomic_5pr)

                                # Create VCF
                                vcf_dict = va_H2V.hgvs2vcf(reverse_normalized_hgvs_genomic, primary_assembly,
                                                           reverse_normalizer, sf)
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
                                stored_hgvs_not_delins = hp.parse_hgvs_variant(str(
                                    hgvs_genomic_5pr.ac) + ':' + hgvs_genomic_5pr.type + '.' + pos + '_' + end + 'del' + ref + 'ins' + alt)
                                v = [chr, pos, ref, alt]

                                # Detect intronic variation using normalization
                                intronic_variant = 'false'

                                # Save a copy of current hgvs_coding
                                try:
                                    saved_hgvs_coding = no_norm_evm.g_to_t(stored_hgvs_not_delins, hgvs_coding.ac)
                                except hgvs.exceptions.HGVSInvalidIntervalError as e:
                                    if str(e) == 'start or end or both are beyond the bounds of transcript record':
                                        saved_hgvs_coding = hgvs_coding
                                        intronic_variant = 'true'
                                        continue
                                    else:
                                        saved_hgvs_coding = no_norm_evm.g_to_n(reverse_normalized_hgvs_genomic,
                                                                               hgvs_coding.ac)

                                # Look for normalized variant options that do not match hgvs_coding
                                if orientation == -1:
                                    # position genomic at its most 5 prime position
                                    try:
                                        query_genomic = reverse_normalizer.normalize(hgvs_genomic)
                                    except:
                                        query_genomic = hgvs_genomic
                                    # Map to the transcript and test for movement
                                    try:
                                        hgvs_seek_var = evm.g_to_t(query_genomic, hgvs_coding.ac)
                                    except hgvs.exceptions.HGVSError as e:
                                        hgvs_seek_var = saved_hgvs_coding
                                    else:
                                        seek_var = valstr(hgvs_seek_var)
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
                                    # Map to the transcript ant test for movement
                                    try:
                                        hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
                                    except hgvs.exceptions.HGVSError as e:
                                        hgvs_seek_var = saved_hgvs_coding
                                    else:
                                        seek_var = valstr(hgvs_seek_var)
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
                                    if re.search('\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\d+\-', str(
                                            hgvs_seek_var.posedit.pos)) or re.search('\*\d+\+', str(
                                        hgvs_seek_var.posedit.pos)) or re.search('\*\d+\-',
                                                                                 str(hgvs_seek_var.posedit.pos)):
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

                                if re.search('\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\d+\-', str(
                                        hgvs_seek_var.posedit.pos)) or re.search('\*\d+\+', str(
                                    hgvs_seek_var.posedit.pos)) or re.search('\*\d+\-', str(hgvs_seek_var.posedit.pos)):
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
                                            if re.search('dup', hgvs_genomic_5pr.posedit.edit.type) or re.search('ins',
                                                                                                                 hgvs_genomic_5pr.posedit.edit.type):
                                                # For gap in chr, map to t. - but becaouse we have pushed to 5 prime by norm, add 1 to end pos
                                                plussed_hgvs_not_delins = copy.deepcopy(hgvs_not_delins)
                                                plussed_hgvs_not_delins.posedit.pos.end.base = plussed_hgvs_not_delins.posedit.pos.end.base + 1
                                                plussed_hgvs_not_delins.posedit.edit.ref = ''
                                                transcript_variant = no_norm_evm.g_to_t(plussed_hgvs_not_delins,
                                                                                        str(saved_hgvs_coding.ac))
                                                if ((
                                                        transcript_variant.posedit.pos.end.base - transcript_variant.posedit.pos.start.base) > (
                                                        hgvs_genomic_5pr.posedit.pos.end.base - hgvs_genomic_5pr.posedit.pos.start.base)):
                                                    if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                                                        hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                        start = hgvs_not_delins.posedit.pos.start.base - 1
                                                        end = hgvs_not_delins.posedit.pos.end.base
                                                        ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                                        hgvs_not_delins.posedit.edit.ref = ref_bases
                                                        hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                           :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                                 1:] + ref_bases[1:]
                                                    elif re.search('ins',
                                                                   str(hgvs_genomic_5pr.posedit.edit)) and re.search(
                                                        'del', str(hgvs_genomic_5pr.posedit.edit)):
                                                        hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                    elif re.search('ins',
                                                                   str(
                                                                       hgvs_genomic_5pr.posedit.edit)) and not re.search(
                                                        'del', str(hgvs_genomic_5pr.posedit.edit)):
                                                        hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                        start = hgvs_not_delins.posedit.pos.start.base - 1
                                                        end = hgvs_not_delins.posedit.pos.end.base
                                                        ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                                        hgvs_not_delins.posedit.edit.ref = ref_bases
                                                        hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                           :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                                 1:] + ref_bases[1:]
                                                else:
                                                    if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                                                        hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                        start = hgvs_not_delins.posedit.pos.start.base - 1
                                                        end = hgvs_not_delins.posedit.pos.end.base
                                                        ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                                        hgvs_not_delins.posedit.edit.ref = ref_bases
                                                        hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                           :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                                 1:] + ref_bases[1:]
                                                    elif re.search('ins',
                                                                   str(hgvs_genomic_5pr.posedit.edit)) and re.search(
                                                        'del', str(hgvs_genomic_5pr.posedit.edit)):
                                                        hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                    elif re.search('ins',
                                                                   str(
                                                                       hgvs_genomic_5pr.posedit.edit)) and not re.search(
                                                        'del', str(hgvs_genomic_5pr.posedit.edit)):
                                                        hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                        start = hgvs_not_delins.posedit.pos.start.base - 1
                                                        end = hgvs_not_delins.posedit.pos.end.base
                                                        ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                                        hgvs_not_delins.posedit.edit.ref = ref_bases
                                                        hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                           :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                                 1:] + ref_bases[1:]
                                            else:
                                                pass
                                        else:
                                            pass
                                        try:
                                            tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                    saved_hgvs_coding.ac)
                                        except hgvs.exceptions.HGVSInvalidIntervalError:
                                            tx_hgvs_not_delins = no_norm_evm.g_to_n(reverse_normalized_hgvs_genomic,
                                                                                    saved_hgvs_coding.ac)
                                        # Create normalized version of tx_hgvs_not_delins
                                        rn_tx_hgvs_not_delins = copy.deepcopy(tx_hgvs_not_delins)

                                        # Check for +1 base and adjust
                                        if re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search(
                                                '\+',
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
                                                pass

                                        elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                                            # move tx end base to next available non-offset base
                                            rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
                                            rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                                            rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                            if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                                test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                            else:
                                                test_tx_var = rn_tx_hgvs_not_delins
                                            # re-make genomic and tx
                                            hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm,
                                                                                   primary_assembly, vm, hp, hn, sf,
                                                                                   nr_vm)

                                            rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                       str(saved_hgvs_coding.ac))

                                        # tx_hgvs_not_delins = rn_tx_hgvs_not_delins
                                        elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                                            # move tx start base to previous available non-offset base
                                            rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                            rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                            if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                                test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                            else:
                                                test_tx_var = rn_tx_hgvs_not_delins
                                            # re-make genomic and tx
                                            hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm,
                                                                                   primary_assembly, vm, hp, hn, sf,
                                                                                   nr_vm)
                                            rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                       str(saved_hgvs_coding.ac))
                                            rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                        #                                         else:
                                        #                                             pass

                                        # Check for -ve base and adjust
                                        elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search(
                                                '\-',
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
                                                pass
                                        elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                                            rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                                            # Delete the ref
                                            rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                            # Add the additional base to the ALT
                                            start = rn_tx_hgvs_not_delins.posedit.pos.end.base - 1
                                            end = rn_tx_hgvs_not_delins.posedit.pos.end.base
                                            ref_bases = sf.fetch_seq(str(tx_hgvs_not_delins.ac), start, end)
                                            rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + ref_bases
                                            if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                                test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                            else:
                                                test_tx_var = rn_tx_hgvs_not_delins
                                            # re-make genomic and tx
                                            hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm,
                                                                                   primary_assembly, vm, hp, hn, sf,
                                                                                   nr_vm)
                                            rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                       str(saved_hgvs_coding.ac))
                                        elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                                            # move tx start base to previous available non-offset base
                                            rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                            rn_tx_hgvs_not_delins.posedit.pos.start.base = rn_tx_hgvs_not_delins.posedit.pos.start.base - 1
                                            rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                            if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                                test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                            else:
                                                test_tx_var = rn_tx_hgvs_not_delins
                                            # re-make genomic and tx
                                            hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm,
                                                                                   primary_assembly, vm, hp, hn, sf,
                                                                                   nr_vm)
                                            rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                       str(saved_hgvs_coding.ac))
                                            rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                        else:
                                            pass

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
                                            for internal_possibility in hgvs_genomic_possibilities:
                                                if internal_possibility == '':
                                                    continue

                                                hgvs_t_possibility = vm.g_to_t(internal_possibility, hgvs_coding.ac)
                                                if hgvs_t_possibility.posedit.edit.type == 'ins':
                                                    try:
                                                        hgvs_t_possibility = vm.c_to_n(hgvs_t_possibility)
                                                    except:
                                                        exceptPass()
                                                    ins_ref = sf.fetch_seq(hgvs_t_possibility.ac,
                                                                           hgvs_t_possibility.posedit.pos.start.base - 1,
                                                                           hgvs_t_possibility.posedit.pos.start.base + 1)
                                                    try:
                                                        hgvs_t_possibility = vm.n_to_c(hgvs_t_possibility)
                                                    except:
                                                        exceptPass()
                                                    hgvs_t_possibility.posedit.edit.ref = ins_ref
                                                    hgvs_t_possibility.posedit.edit.alt = ins_ref[
                                                                                              0] + hgvs_t_possibility.posedit.edit.alt + \
                                                                                          ins_ref[1]
                                                if internal_possibility.posedit.edit.type == 'ins':
                                                    ins_ref = sf.fetch_seq(internal_possibility.ac,
                                                                           internal_possibility.posedit.pos.start.base - 1,
                                                                           internal_possibility.posedit.pos.end.base)
                                                    internal_possibility.posedit.edit.ref = ins_ref
                                                    internal_possibility.posedit.edit.alt = ins_ref[
                                                                                                0] + internal_possibility.posedit.edit.alt + \
                                                                                            ins_ref[1]

                                                if len(hgvs_t_possibility.posedit.edit.ref) < len(
                                                        internal_possibility.posedit.edit.ref):
                                                    gap_length = len(internal_possibility.posedit.edit.ref) - len(
                                                        hgvs_t_possibility.posedit.edit.ref)
                                                    re_capture_tx_variant = ['transcript', gap_length,
                                                                             hgvs_t_possibility]
                                                    hgvs_not_delins = internal_possibility
                                                    hgvs_genomic_5pr = internal_possibility
                                                    break

                                            if re_capture_tx_variant != []:
                                                try:
                                                    tx_hgvs_not_delins = vm.c_to_n(re_capture_tx_variant[2])
                                                except:
                                                    tx_hgvs_not_delins = re_capture_tx_variant[2]
                                                disparity_deletion_in = re_capture_tx_variant[0:-1]
                                            else:
                                                pass

                                    # 'At hgvs_genomic'
                                    # Final sanity checks
                                    try:
                                        vm.g_to_t(hgvs_not_delins, tx_hgvs_not_delins.ac)
                                    except Exception as e:
                                        if str(e) == 'start or end or both are beyond the bounds of transcript record':
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
                                            elif re.match('Normalization of intronic variants is not supported', error):
                                                # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                                                disparity_deletion_in = ['transcript', 'Requires Analysis']

                                    # amend_RefSeqGene = 'false'
                                    # Recreate hgvs_genomic
                                    if disparity_deletion_in[0] == 'transcript':
                                        hgvs_genomic = hgvs_not_delins

                                    # Find oddly placed gaps where the tx variant is encompassed in the gap
                                    if disparity_deletion_in[0] == 'false' and (
                                            possibility_counter == 3 or possibility_counter == 4):
                                        rg = reverse_normalizer.normalize(hgvs_not_delins)
                                        rtx = vm.g_to_t(rg, tx_hgvs_not_delins.ac)
                                        fg = hn.normalize(hgvs_not_delins)
                                        ftx = vm.g_to_t(fg, tx_hgvs_not_delins.ac)
                                        if (rtx.posedit.pos.start.offset == 0 and rtx.posedit.pos.end.offset == 0) and (
                                                ftx.posedit.pos.start.offset != 0 and ftx.posedit.pos.end.offset != 0):
                                            exons = hdp.get_tx_exons(ftx.ac, hgvs_not_delins.ac, alt_aln_method)
                                            exonic = False
                                            for ex_test in exons:
                                                if ftx.posedit.pos.start.base in range(ex_test[6], ex_test[
                                                    7]) and ftx.posedit.pos.end.base in range(ex_test[6], ex_test[7]):
                                                    exonic = True
                                            if exonic is True:
                                                hgvs_not_delins = fg
                                                hgvs_genomic = fg
                                                hgvs_genomic_5pr = fg
                                                try:
                                                    tx_hgvs_not_delins = vm.c_to_n(ftx)
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
                                            tx_hgvs_not_delins = hp.parse_hgvs_variant(
                                                tx_hgvs_not_delins_delins_from_dup)

                                    # GAP IN THE TRANSCRIPT DISPARITY DETECTED
                                    if disparity_deletion_in[0] == 'transcript':
                                        # Suppress intron boundary crossing due to non-intron intron based c. seq annotations
                                        suppress_c_normalization = 'true'
                                        # amend_RefSeqGene = 'true'
                                        # ANY VARIANT WHOLLY WITHIN THE GAP
                                        if (re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) or re.search(
                                                '\-',
                                                str(
                                                    tx_hgvs_not_delins.posedit.pos.start))) and (
                                                re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) or re.search(
                                            '\-',
                                            str(
                                                tx_hgvs_not_delins.posedit.pos.end))):

                                            gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)

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
                                                    tx_gap_fill_variant = hp.parse_hgvs_variant(
                                                        tx_gap_fill_variant_delins_from_dup)

                                                    # Identify which half of the NOT-intron the start position of the variant is in
                                            if re.search('\-', str(tx_gap_fill_variant.posedit.pos.start)):
                                                tx_gap_fill_variant.posedit.pos.start.base = tx_gap_fill_variant.posedit.pos.start.base - 1
                                                tx_gap_fill_variant.posedit.pos.start.offset = int('0')  # int('+1')
                                                tx_gap_fill_variant.posedit.pos.end.offset = int('0')  # int('-1')
                                                tx_gap_fill_variant.posedit.edit.alt = ''
                                                tx_gap_fill_variant.posedit.edit.ref = ''
                                            elif re.search('\+', str(tx_gap_fill_variant.posedit.pos.start)):
                                                tx_gap_fill_variant.posedit.pos.start.offset = int('0')  # int('+1')
                                                tx_gap_fill_variant.posedit.pos.end.base = tx_gap_fill_variant.posedit.pos.end.base + 1
                                                tx_gap_fill_variant.posedit.pos.end.offset = int('0')  # int('-1')
                                                tx_gap_fill_variant.posedit.edit.alt = ''
                                                tx_gap_fill_variant.posedit.edit.ref = ''

                                            try:
                                                tx_gap_fill_variant = vm.n_to_c(tx_gap_fill_variant)
                                            except:
                                                exceptPass()
                                            genomic_gap_fill_variant = vm.t_to_g(tx_gap_fill_variant,
                                                                                 reverse_normalized_hgvs_genomic.ac)
                                            genomic_gap_fill_variant.posedit.edit.alt = genomic_gap_fill_variant.posedit.edit.ref

                                            try:
                                                c_tx_hgvs_not_delins = vm.n_to_c(tx_hgvs_not_delins)
                                            except Exception:
                                                c_tx_hgvs_not_delins = copy.copy(tx_hgvs_not_delins)
                                            genomic_gap_fill_variant_alt = vm.t_to_g(c_tx_hgvs_not_delins,
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
                                                    genomic_gap_fill_variant = hp.parse_hgvs_variant(
                                                        genomic_gap_fill_variant_delins_from_dup)
                                                    genomic_gap_fill_variant_alt_delins_from_dup = genomic_gap_fill_variant_alt.ac + ':' + genomic_gap_fill_variant_alt.type + '.' + str(
                                                        genomic_gap_fill_variant_alt.posedit.pos.start.base) + '_' + str(
                                                        genomic_gap_fill_variant_alt.posedit.pos.end.base) + 'del' + genomic_gap_fill_variant_alt.posedit.edit.ref + 'ins' + genomic_gap_fill_variant_alt.posedit.edit.ref + genomic_gap_fill_variant_alt.posedit.edit.ref
                                                    genomic_gap_fill_variant_alt = hp.parse_hgvs_variant(
                                                        genomic_gap_fill_variant_alt_delins_from_dup)

                                            # Correct insertion alts
                                            if genomic_gap_fill_variant_alt.posedit.edit.type == 'ins':
                                                append_ref = sf.fetch_seq(genomic_gap_fill_variant_alt.ac,
                                                                          genomic_gap_fill_variant_alt.posedit.pos.start.base - 1,
                                                                          genomic_gap_fill_variant_alt.posedit.pos.end.base)
                                                genomic_gap_fill_variant_alt.posedit.edit.alt = append_ref[
                                                                                                    0] + genomic_gap_fill_variant_alt.posedit.edit.alt + \
                                                                                                append_ref[1]

                                            # Split the reference and replacing alt sequence into a dictionary
                                            reference_bases = list(genomic_gap_fill_variant.posedit.edit.ref)
                                            if genomic_gap_fill_variant_alt.posedit.edit.alt is not None:
                                                alternate_bases = list(genomic_gap_fill_variant_alt.posedit.edit.alt)
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
                                            for integer in range(genomic_gap_fill_variant_alt.posedit.pos.start.base,
                                                                 genomic_gap_fill_variant_alt.posedit.pos.end.base + 1,
                                                                 1):
                                                if integer == alt_start:
                                                    alt_base_dict[integer] = str(''.join(alternate_bases))
                                                else:
                                                    alt_base_dict[integer] = 'X'

                                            # Generate the alt sequence
                                            alternate_sequence_bases = []
                                            for integer in range(genomic_gap_fill_variant.posedit.pos.start.base,
                                                                 genomic_gap_fill_variant.posedit.pos.end.base + 1, 1):
                                                if integer in alt_base_dict.keys():
                                                    alternate_sequence_bases.append(alt_base_dict[integer])
                                                else:
                                                    alternate_sequence_bases.append(ref_base_dict[integer])
                                            alternate_sequence = ''.join(alternate_sequence_bases)
                                            alternate_sequence = alternate_sequence.replace('X', '')

                                            # Add the new alt to the gap fill variant and generate transcript variant
                                            genomic_gap_fill_variant.posedit.edit.alt = alternate_sequence
                                            hgvs_refreshed_variant = vm.g_to_t(genomic_gap_fill_variant,
                                                                               tx_gap_fill_variant.ac)

                                            # Set warning
                                            gap_size = str(len(genomic_gap_fill_variant.posedit.edit.ref) - 2)
                                            disparity_deletion_in[1] = [gap_size]
                                            auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                                                stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + gap_size + ' genomic base(s) that fail to align to transcript ' + str(
                                                tx_hgvs_not_delins.ac)
                                            non_valid_caution = 'true'

                                            # Alignment position
                                            for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                            if re.match('NM_', str(for_location_c)):
                                                for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                            if re.match('\-', str(for_location_c.posedit.pos.start.offset)):
                                                gps = for_location_c.posedit.pos.start.base - 1
                                                gpe = for_location_c.posedit.pos.start.base
                                            else:
                                                gps = for_location_c.posedit.pos.start.base
                                                gpe = for_location_c.posedit.pos.start.base + 1
                                            gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                                            auto_info = auto_info + '%s' % (gap_position)

                                        else:
                                            if tx_hgvs_not_delins.posedit.pos.start.offset == 0 and tx_hgvs_not_delins.posedit.pos.end.offset == 0:
                                                # In this instance, we have identified a transcript gap but the n. version of
                                                # the transcript variant but do not have a position which actually hits the gap,
                                                # so the variant likely spans the gap, and is not picked up by an offset.
                                                try:
                                                    c1 = vm.n_to_c(tx_hgvs_not_delins)
                                                except:
                                                    c1 = tx_hgvs_not_delins
                                                g1 = nr_vm.t_to_g(c1, hgvs_genomic.ac)
                                                g3 = nr_vm.t_to_g(c1, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                ng2 = hn.normalize(g2)
                                                g3.posedit.pos.end.base = g3.posedit.pos.start.base + (
                                                        len(g3.posedit.edit.ref) - 1)
                                                try:
                                                    c2 = vm.g_to_t(g3, c1.ac)
                                                    if c2.posedit.pos.start.offset == 0 and c2.posedit.pos.end.offset == 0:
                                                        pass
                                                    else:
                                                        tx_hgvs_not_delins = c2
                                                        try:
                                                            tx_hgvs_not_delins = vm.c_to_n(tx_hgvs_not_delins)
                                                        except hgvs.exceptions.HGVSError:
                                                            exceptPass()
                                                except hgvs.exceptions.HGVSInvalidVariantError:
                                                    exceptPass()

                                            if re.search('\+',
                                                         str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                                                '\+', str(tx_hgvs_not_delins.posedit.pos.end)):
                                                auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                                                    stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                                                    disparity_deletion_in[
                                                        1]) + ' genomic base(s) that fail to align to transcript ' + str(
                                                    tx_hgvs_not_delins.ac)
                                                non_valid_caution = 'true'
                                                try:
                                                    c2 = vm.n_to_c(tx_hgvs_not_delins)
                                                except:
                                                    c2 = tx_hgvs_not_delins
                                                c1 = copy.deepcopy(c2)
                                                c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                                                c1.posedit.pos.start.offset = 0
                                                c1.posedit.pos.end = c2.posedit.pos.start
                                                c1.posedit.edit.ref = ''
                                                c1.posedit.edit.alt = ''
                                                if orientation != -1:
                                                    g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                    g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                    g1.posedit.edit.alt = g1.posedit.edit.ref
                                                else:
                                                    g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                    g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                    g2.posedit.edit.alt = g2.posedit.edit.ref
                                                reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                                alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                                g3 = copy.deepcopy(g1)
                                                g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                                g3.posedit.edit.ref = reference
                                                g3.posedit.edit.alt = alternate
                                                c3 = vm.g_to_t(g3, c1.ac)
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
                                            elif re.search('\+',
                                                           str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search(
                                                '\+', str(tx_hgvs_not_delins.posedit.pos.start)):
                                                auto_info = auto_info + 'Genome position ' + str(
                                                    stored_hgvs_not_delins.ac) + ':g.' + str(
                                                    stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(
                                                    disparity_deletion_in[1]) + '-bp gap in transcript ' + str(
                                                    tx_hgvs_not_delins.ac)
                                                gapped_transcripts = gapped_transcripts + ' ' + str(
                                                    tx_hgvs_not_delins.ac)
                                                non_valid_caution = 'true'
                                                try:
                                                    c1 = vm.n_to_c(tx_hgvs_not_delins)
                                                except:
                                                    c1 = tx_hgvs_not_delins
                                                c2 = copy.deepcopy(c1)
                                                c2.posedit.pos.start = c1.posedit.pos.end
                                                c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                                                c2.posedit.pos.end.offset = 0
                                                c2.posedit.edit.ref = ''
                                                c2.posedit.edit.alt = ''
                                                if orientation != -1:
                                                    g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                    g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                    g2.posedit.edit.alt = g2.posedit.edit.ref
                                                else:
                                                    g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                    g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                    g1.posedit.edit.alt = g1.posedit.edit.ref
                                                reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                                alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                                g3 = copy.deepcopy(g1)
                                                g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                                g3.posedit.edit.ref = reference
                                                g3.posedit.edit.alt = alternate
                                                c3 = vm.g_to_t(g3, c1.ac)
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
                                            elif re.search('\-',
                                                           str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                                                '\-', str(tx_hgvs_not_delins.posedit.pos.end)):
                                                auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                                                    stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                                                    disparity_deletion_in[
                                                        1]) + ' genomic base(s) that fail to align to transcript ' + str(
                                                    tx_hgvs_not_delins.ac)
                                                non_valid_caution = 'true'
                                                try:
                                                    c2 = vm.n_to_c(tx_hgvs_not_delins)
                                                except:
                                                    c2 = tx_hgvs_not_delins
                                                c1 = copy.deepcopy(c2)
                                                c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                                                c1.posedit.pos.start.offset = 0
                                                c1.posedit.pos.end = c2.posedit.pos.start
                                                c1.posedit.edit.ref = ''
                                                c1.posedit.edit.alt = ''
                                                if orientation != -1:
                                                    g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                    g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                    g1.posedit.edit.alt = g1.posedit.edit.ref
                                                else:
                                                    g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                    g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                    g2.posedit.edit.alt = g2.posedit.edit.ref
                                                reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                                alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                                g3 = copy.deepcopy(g1)
                                                g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                                g3.posedit.edit.ref = reference
                                                g3.posedit.edit.alt = alternate
                                                c3 = vm.g_to_t(g3, c1.ac)
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
                                            elif re.search('\-',
                                                           str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search(
                                                '\-', str(tx_hgvs_not_delins.posedit.pos.start)):
                                                auto_info = auto_info + 'Genome position ' + str(
                                                    stored_hgvs_not_delins.ac) + ':g.' + str(
                                                    stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(
                                                    disparity_deletion_in[1]) + '-bp gap in transcript ' + str(
                                                    tx_hgvs_not_delins.ac)
                                                gapped_transcripts = gapped_transcripts + ' ' + str(
                                                    tx_hgvs_not_delins.ac)
                                                non_valid_caution = 'true'
                                                try:
                                                    c1 = vm.n_to_c(tx_hgvs_not_delins)
                                                except:
                                                    c1 = tx_hgvs_not_delins
                                                c2 = copy.deepcopy(c1)
                                                c2.posedit.pos.start = c1.posedit.pos.end
                                                c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                                                c2.posedit.pos.end.offset = 0
                                                c2.posedit.edit.ref = ''
                                                c2.posedit.edit.alt = ''
                                                if orientation != -1:
                                                    g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                    g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                    g2.posedit.edit.alt = g2.posedit.edit.ref
                                                else:
                                                    g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                    g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                    g1.posedit.edit.alt = g1.posedit.edit.ref
                                                reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                                alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                                g3 = copy.deepcopy(g1)
                                                g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                                g3.posedit.edit.ref = reference
                                                g3.posedit.edit.alt = alternate
                                                c3 = vm.g_to_t(g3, c1.ac)
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
                                                auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                                                    stored_hgvs_not_delins.posedit.pos) + ' contains ' + str(
                                                    disparity_deletion_in[
                                                        1]) + ' genomic base(s) that fail to align to transcript ' + str(
                                                    tx_hgvs_not_delins.ac) + '\n'
                                                tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(
                                                    tx_hgvs_not_delins.posedit.edit.ref) - 1
                                                hgvs_refreshed_variant = tx_hgvs_not_delins

                                    # GAP IN THE CHROMOSOME
                                    elif disparity_deletion_in[0] == 'chromosome':
                                        suppress_c_normalization = 'true'
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
                                    if re.match('NM_', str(hgvs_refreshed_variant.ac)) and not re.search('c', str(
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
                                            logger.warning(error)
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
                                    hgvs_genomic = va_func.myvm_t_to_g(hgvs_refreshed_variant, hgvs_genomic.ac,
                                                                       no_norm_evm, vm, hp, hn, sf, nr_vm)
                                    if hgvs_genomic.posedit.edit.type == 'identity':
                                        re_c = vm.g_to_t(hgvs_genomic, hgvs_refreshed_variant.ac)
                                        if (hn.normalize(re_c)) != (hn.normalize(hgvs_refreshed_variant)):
                                            shuffle_left_g = copy.copy(hgvs_genomic)
                                            shuffle_left_g.posedit.edit.ref = ''
                                            shuffle_left_g.posedit.edit.alt = ''
                                            shuffle_left_g.posedit.pos.start.base = shuffle_left_g.posedit.pos.start.base - 1
                                            shuffle_left_g.posedit.pos.end.base = shuffle_left_g.posedit.pos.end.base - 1
                                            shuffle_left_g = reverse_normalizer.normalize(shuffle_left_g)
                                            re_c = vm.g_to_t(shuffle_left_g, hgvs_refreshed_variant.ac)
                                            if (hn.normalize(re_c)) != (hn.normalize(hgvs_refreshed_variant)):
                                                hgvs_genomic = shuffle_left_g

                                # If it is intronic, these vairables will not have been set
                                else:
                                    # amend_RefSeqGene = 'false'
                                    no_normalized_c = 'false'

                                # Break if gap has been detected
                                if disparity_deletion_in[0] != 'false':
                                    break

                            # Warn user about gapping
                            if auto_info != '':
                                info_lines = auto_info.split('\n')
                                info_keys = {}
                                for information in info_lines:
                                    info_keys[information] = ''
                                info_out = []
                                info_out.append(
                                    'The displayed variants may be artefacts of aligning ' + hgvs_coding.ac + ' with genome build ' + primary_assembly)
                                for ky in info_keys.keys():
                                    info_out.append(ky)
                                auto_info = '\n'.join(info_out)
                                auto_info = auto_info + '\nCaution should be used when reporting the displayed variant descriptions: If you are unsure, please contact admin'
                                auto_info = str(auto_info.replace('\n', ': '))
                                validation['warnings'] = validation['warnings'] + ': ' + str(auto_info)
                                logger.warning(str(auto_info))
                            # Normailse hgvs_genomic
                            try:
                                hgvs_genomic = hn.normalize(hgvs_genomic)
                            except hgvs.exceptions.HGVSError as e:
                                # Strange error caused by gap in genomic
                                error = str(e)

                                if re.search('base start position must be <= end position', error) and \
                                        disparity_deletion_in[0] == 'chromosome':
                                    if hgvs_genomic.posedit.edit.type == 'delins':
                                        start = hgvs_genomic.posedit.pos.start.base
                                        end = hgvs_genomic.posedit.pos.end.base
                                        lhb = sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                                        rhb = sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                                        hgvs_genomic.posedit.edit.ref = lhb + rhb
                                        hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                                        hgvs_genomic.posedit.pos.start.base = end
                                        hgvs_genomic.posedit.pos.end.base = start
                                        hgvs_genomic = hn.normalize(hgvs_genomic)
                                    if hgvs_genomic.posedit.edit.type == 'del':
                                        start = hgvs_genomic.posedit.pos.start.base
                                        end = hgvs_genomic.posedit.pos.end.base
                                        lhb = sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                                        rhb = sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                                        hgvs_genomic.posedit.edit.ref = lhb + rhb
                                        hgvs_genomic.posedit.edit.alt = lhb + rhb
                                        hgvs_genomic.posedit.pos.start.base = end
                                        hgvs_genomic.posedit.pos.end.base = start
                                        hgvs_genomic = hn.normalize(hgvs_genomic)
                            genomic = valstr(hgvs_genomic)

                        else:
                            stored_hgvs_genomic_variant = hgvs_genomic
                            suppress_c_normalization = 'false'
                            gapped_alignment_warning = ''
                            auto_info = ''
                            genomic = valstr(hgvs_genomic)

                        # Create pseudo VCF based on amended hgvs_genomic
                        hgvs_genomic_variant = hgvs_genomic
                        # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
                        reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic_variant)

                        hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

                        # Create vcf
                        vcf_dict = va_H2V.hgvs2vcf(reverse_normalized_hgvs_genomic, primary_assembly,
                                                   reverse_normalizer, sf)
                        chr = vcf_dict['chr']
                        pos = vcf_dict['pos']
                        ref = vcf_dict['ref']
                        alt = vcf_dict['alt']

                        # Create a VCF call
                        vcf_component_list = [str(chr), str(pos), str(ref), (alt)]
                        vcf_genomic = '-'.join(vcf_component_list)

                        # DO NOT DELETE
                        # Generate an end position
                        end = str(int(pos) + len(ref) - 1)
                        pos = str(pos)

                        # DO NOT DELETE
                        stored_hgvs_not_delins = hp.parse_hgvs_variant(str(
                            hgvs_genomic_5pr.ac) + ':' + hgvs_genomic_5pr.type + '.' + pos + '_' + end + 'del' + ref + 'ins' + alt)

                        # Apply gap code to re-format hgvs_coding
                        # Store the current hgvs:c. description
                        saved_hgvs_coding = copy.deepcopy(hgvs_coding)

                        # Get orientation of the gene wrt genome and a list of exons mapped to the genome
                        ori = va_func.tx_exons(tx_ac=saved_hgvs_coding.ac, alt_ac=hgvs_genomic_5pr.ac,
                                               alt_aln_method=alt_aln_method, hdp=hdp)
                        orientation = int(ori[0]['alt_strand'])

                        # Look for normalized variant options that do not match hgvs_coding
                        hgvs_genomic = copy.deepcopy(hgvs_genomic_variant)
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
                                seek_var = valstr(hgvs_seek_var)
                                seek_ac = str(hgvs_seek_var.ac)
                            if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
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
                            # Map to the transcript ant test for movement
                            try:
                                hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
                            except hgvs.exceptions.HGVSError as e:
                                hgvs_seek_var = saved_hgvs_coding
                            else:
                                seek_var = valstr(hgvs_seek_var)
                                seek_ac = str(hgvs_seek_var.ac)
                            if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                                    saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                                    hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                                    saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset):
                                pass
                            else:
                                hgvs_seek_var = saved_hgvs_coding

                        # Loop out gap finding code under these circumstances!
                        logger.warning("gap_compensation_2 = " + str(gap_compensation))
                        if gap_compensation is True:
                            logger.warning('g_to_t gap code 2 active')
                            # is it in an exon?
                            is_it_in_an_exon = 'no'
                            for exon in ori:
                                genomic_start = int(exon['alt_start_i'])
                                genomic_end = int(exon['alt_end_i'])
                                # Take from stored copy
                                # hgvs_genomic_5pr = copy.deepcopy(stored_hgvs_genomic_5pr)
                                if (
                                        hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                                        hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                                    is_it_in_an_exon = 'yes'
                            if is_it_in_an_exon == 'yes':
                                # map form reverse normalized g. to c.
                                hgvs_from_5n_g = no_norm_evm.g_to_t(hgvs_genomic_5pr, saved_hgvs_coding.ac)

                                # Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping transcript lengths
                                disparity_deletion_in = ['false', 'false']
                                if stored_hgvs_not_delins != '':
                                    # Refresh hgvs_not_delins from stored_hgvs_not_delins
                                    hgvs_not_delins = copy.deepcopy(stored_hgvs_not_delins)
                                    # This test will only occur in dup of single base, insertion or substitution
                                    if not re.search('_', str(hgvs_not_delins.posedit.pos)):
                                        if re.search('dup', hgvs_genomic_5pr.posedit.edit.type) or re.search('ins',
                                                                                                             hgvs_genomic_5pr.posedit.edit.type):
                                            # For gap in chr, map to t. - but becaouse we have pushed to 5 prime by norm, add 1 to end pos
                                            plussed_hgvs_not_delins = copy.deepcopy(hgvs_not_delins)
                                            plussed_hgvs_not_delins.posedit.pos.end.base = plussed_hgvs_not_delins.posedit.pos.end.base + 1
                                            plussed_hgvs_not_delins.posedit.edit.ref = ''
                                            transcript_variant = no_norm_evm.g_to_t(plussed_hgvs_not_delins,
                                                                                    str(saved_hgvs_coding.ac))
                                            if ((
                                                    transcript_variant.posedit.pos.end.base - transcript_variant.posedit.pos.start.base) > (
                                                    hgvs_genomic_5pr.posedit.pos.end.base - hgvs_genomic_5pr.posedit.pos.start.base)):
                                                if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                                                    hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                    start = hgvs_not_delins.posedit.pos.start.base - 1
                                                    end = hgvs_not_delins.posedit.pos.end.base
                                                    ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                                    hgvs_not_delins.posedit.edit.ref = ref_bases
                                                    hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                       :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                             1:] + ref_bases[1:]
                                                elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search(
                                                        'del', str(hgvs_genomic_5pr.posedit.edit)):
                                                    hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                elif re.search('ins',
                                                               str(hgvs_genomic_5pr.posedit.edit)) and not re.search(
                                                    'del', str(hgvs_genomic_5pr.posedit.edit)):
                                                    hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                    start = hgvs_not_delins.posedit.pos.start.base - 1
                                                    end = hgvs_not_delins.posedit.pos.end.base
                                                    ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                                    hgvs_not_delins.posedit.edit.ref = ref_bases
                                                    hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                       :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                             1:] + ref_bases[1:]
                                            else:
                                                if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                                                    hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                    start = hgvs_not_delins.posedit.pos.start.base - 1
                                                    end = hgvs_not_delins.posedit.pos.end.base
                                                    ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                                    hgvs_not_delins.posedit.edit.ref = ref_bases
                                                    hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                       :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                             1:] + ref_bases[1:]
                                                elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search(
                                                        'del', str(hgvs_genomic_5pr.posedit.edit)):
                                                    hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                elif re.search('ins',
                                                               str(hgvs_genomic_5pr.posedit.edit)) and not re.search(
                                                    'del', str(hgvs_genomic_5pr.posedit.edit)):
                                                    hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                                                    start = hgvs_not_delins.posedit.pos.start.base - 1
                                                    end = hgvs_not_delins.posedit.pos.end.base
                                                    ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                                                    hgvs_not_delins.posedit.edit.ref = ref_bases
                                                    hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                                                       :1] + hgvs_not_delins.posedit.edit.alt[
                                                                                             1:] + ref_bases[1:]
                                        else:
                                            pass
                                    else:
                                        pass

                                    hard_fail = 'false'
                                    try:
                                        tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, saved_hgvs_coding.ac)
                                    except Exception as e:
                                        if str(e) == 'start or end or both are beyond the bounds of transcript record':
                                            tx_hgvs_not_delins = hgvs_coding
                                            hard_fail = 'true'

                                    # Create normalized version of tx_hgvs_not_delins
                                    rn_tx_hgvs_not_delins = copy.deepcopy(tx_hgvs_not_delins)
                                    # Check for +ve base and adjust
                                    if re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search('\+',
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
                                            exceptPass()

                                    elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                                        # move tx end base to next available non-offset base
                                        rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
                                        rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                                        rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                        if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                            test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                        else:
                                            test_tx_var = rn_tx_hgvs_not_delins
                                        # re-make genomic and tx
                                        hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm,
                                                                               primary_assembly, vm, hp, hn, sf, nr_vm)
                                        rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                   str(saved_hgvs_coding.ac))
                                    elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                                        rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                        # move tx start base to previous available non-offset base
                                        rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                        if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                            test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                        else:
                                            test_tx_var = rn_tx_hgvs_not_delins
                                        # re-make genomic and tx
                                        hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm,
                                                                               primary_assembly, vm, hp, hn, sf, nr_vm)
                                        rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                   str(saved_hgvs_coding.ac))
                                        rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                    #                                     else:
                                    #                                         pass

                                    # Check for -ve base and adjust
                                    elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search('\-',
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
                                            exceptPass()
                                    elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                                        # move tx end base back to next available non-offset base
                                        # rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base - 1
                                        rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                                        # Delete the ref
                                        rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                        # Add the additional base to the ALT
                                        start = rn_tx_hgvs_not_delins.posedit.pos.end.base - 1
                                        end = rn_tx_hgvs_not_delins.posedit.pos.end.base
                                        ref_bases = sf.fetch_seq(str(tx_hgvs_not_delins.ac), start, end)
                                        rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + ref_bases
                                        if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                            test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                        else:
                                            test_tx_var = rn_tx_hgvs_not_delins
                                        # re-make genomic and tx
                                        hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm,
                                                                               primary_assembly, vm, hp, hn, sf, nr_vm)
                                        rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                   str(saved_hgvs_coding.ac))
                                    elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                                        # move tx start base to previous available non-offset base
                                        rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                        rn_tx_hgvs_not_delins.posedit.pos.start.base = rn_tx_hgvs_not_delins.posedit.pos.start.base - 1
                                        rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                        if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                            test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                        else:
                                            test_tx_var = rn_tx_hgvs_not_delins
                                        # re-make genomic and tx
                                        hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm,
                                                                               primary_assembly, vm, hp, hn, sf, nr_vm)
                                        rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                   str(saved_hgvs_coding.ac))
                                        rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                    else:
                                        pass

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
                                        for internal_possibility in hgvs_genomic_possibilities:

                                            if internal_possibility == '':
                                                continue

                                            hgvs_t_possibility = vm.g_to_t(internal_possibility, hgvs_coding.ac)
                                            if hgvs_t_possibility.posedit.edit.type == 'ins':
                                                try:
                                                    hgvs_t_possibility = vm.c_to_n(hgvs_t_possibility)
                                                except:
                                                    exceptPass()
                                                ins_ref = sf.fetch_seq(hgvs_t_possibility.ac,
                                                                       hgvs_t_possibility.posedit.pos.start.base - 1,
                                                                       hgvs_t_possibility.posedit.pos.start.base + 1)
                                                try:
                                                    hgvs_t_possibility = vm.n_to_c(hgvs_t_possibility)
                                                except:
                                                    exceptPass()
                                                hgvs_t_possibility.posedit.edit.ref = ins_ref
                                                hgvs_t_possibility.posedit.edit.alt = ins_ref[
                                                                                          0] + hgvs_t_possibility.posedit.edit.alt + \
                                                                                      ins_ref[1]
                                            if internal_possibility.posedit.edit.type == 'ins':
                                                ins_ref = sf.fetch_seq(internal_possibility.ac,
                                                                       internal_possibility.posedit.pos.start.base - 1,
                                                                       internal_possibility.posedit.pos.end.base)
                                                internal_possibility.posedit.edit.ref = ins_ref
                                                internal_possibility.posedit.edit.alt = ins_ref[
                                                                                            0] + internal_possibility.posedit.edit.alt + \
                                                                                        ins_ref[1]

                                            if len(hgvs_t_possibility.posedit.edit.ref) < len(
                                                    internal_possibility.posedit.edit.ref):
                                                gap_length = len(internal_possibility.posedit.edit.ref) - len(
                                                    hgvs_t_possibility.posedit.edit.ref)
                                                re_capture_tx_variant = ['transcript', gap_length, hgvs_t_possibility]
                                                hgvs_not_delins = internal_possibility
                                                hgvs_genomic_5pr = internal_possibility
                                                break

                                        if re_capture_tx_variant != []:
                                            try:
                                                tx_hgvs_not_delins = vm.c_to_n(re_capture_tx_variant[2])
                                            except:
                                                tx_hgvs_not_delins = re_capture_tx_variant[2]
                                            disparity_deletion_in = re_capture_tx_variant[0:-1]
                                        else:
                                            pass

                                # Final sanity checks
                                try:
                                    vm.g_to_t(hgvs_not_delins, tx_hgvs_not_delins.ac)
                                except Exception as e:
                                    if str(e) == 'start or end or both are beyond the bounds of transcript record':
                                        logger.warning(str(e))
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
                                            logger.warning(error)
                                            continue
                                        elif re.match('Normalization of intronic variants is not supported', error):
                                            # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                                            disparity_deletion_in = ['transcript', 'Requires Analysis']

                                if hard_fail == 'true':
                                    disparity_deletion_in = ['false', 'false']

                                # Recreate hgvs_genomic
                                if disparity_deletion_in[0] == 'transcript':
                                    hgvs_genomic = hgvs_not_delins

                                # Pre-processing of tx_hgvs_not_delins
                                try:
                                    if tx_hgvs_not_delins.posedit.edit.alt is None:
                                        tx_hgvs_not_delins.posedit.edit.alt = ''
                                except Exception as e:
                                    if str(e) == "'Dup' object has no attribute 'alt'":
                                        tx_hgvs_not_delins_delins_from_dup = tx_hgvs_not_delins.ac + ':' + tx_hgvs_not_delins.type + '.' + str(
                                            tx_hgvs_not_delins.posedit.pos.start) + '_' + str(
                                            tx_hgvs_not_delins.posedit.pos.end) + 'del' + tx_hgvs_not_delins.posedit.edit.ref + 'ins' + tx_hgvs_not_delins.posedit.edit.ref + tx_hgvs_not_delins.posedit.edit.ref
                                        tx_hgvs_not_delins = hp.parse_hgvs_variant(tx_hgvs_not_delins_delins_from_dup)


                                # GAP IN THE TRANSCRIPT DISPARITY DETECTED
                                if disparity_deletion_in[0] == 'transcript':
                                    gap_position = ''
                                    gapped_alignment_warning = str(
                                        hgvs_genomic_5pr) + ' does not represent a true variant because it is an artefact of aligning the transcripts listed below with genome build ' + primary_assembly

                                    # ANY VARIANT WHOLLY WITHIN THE GAP
                                    if (re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) or re.search('\-',
                                                                                                                str(
                                                                                                                    tx_hgvs_not_delins.posedit.pos.start))) and (
                                            re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) or re.search('\-',
                                                                                                                  str(
                                                                                                                      tx_hgvs_not_delins.posedit.pos.end))):
                                        gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)

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
                                                tx_gap_fill_variant = hp.parse_hgvs_variant(
                                                    tx_gap_fill_variant_delins_from_dup)

                                        # Identify which half of the NOT-intron the start position of the variant is in
                                        if re.search('\-', str(tx_gap_fill_variant.posedit.pos.start)):
                                            tx_gap_fill_variant.posedit.pos.start.base = tx_gap_fill_variant.posedit.pos.start.base - 1
                                            tx_gap_fill_variant.posedit.pos.start.offset = int('0')  # int('+1')
                                            tx_gap_fill_variant.posedit.pos.end.offset = int('0')  # int('-1')
                                            tx_gap_fill_variant.posedit.edit.alt = ''
                                            tx_gap_fill_variant.posedit.edit.ref = ''
                                        elif re.search('\+', str(tx_gap_fill_variant.posedit.pos.start)):
                                            tx_gap_fill_variant.posedit.pos.start.offset = int('0')  # int('+1')
                                            tx_gap_fill_variant.posedit.pos.end.base = tx_gap_fill_variant.posedit.pos.end.base + 1
                                            tx_gap_fill_variant.posedit.pos.end.offset = int('0')  # int('-1')
                                            tx_gap_fill_variant.posedit.edit.alt = ''
                                            tx_gap_fill_variant.posedit.edit.ref = ''

                                        try:
                                            tx_gap_fill_variant = vm.n_to_c(tx_gap_fill_variant)
                                        except:
                                            exceptPass()
                                        genomic_gap_fill_variant = vm.t_to_g(tx_gap_fill_variant,
                                                                             reverse_normalized_hgvs_genomic.ac)
                                        genomic_gap_fill_variant.posedit.edit.alt = genomic_gap_fill_variant.posedit.edit.ref

                                        try:
                                            c_tx_hgvs_not_delins = vm.n_to_c(tx_hgvs_not_delins)
                                        except Exception:
                                            c_tx_hgvs_not_delins = copy.copy(tx_hgvs_not_delins)
                                        genomic_gap_fill_variant_alt = vm.t_to_g(c_tx_hgvs_not_delins,
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
                                                genomic_gap_fill_variant = hp.parse_hgvs_variant(
                                                    genomic_gap_fill_variant_delins_from_dup)
                                                genomic_gap_fill_variant_alt_delins_from_dup = genomic_gap_fill_variant_alt.ac + ':' + genomic_gap_fill_variant_alt.type + '.' + str(
                                                    genomic_gap_fill_variant_alt.posedit.pos.start.base) + '_' + str(
                                                    genomic_gap_fill_variant_alt.posedit.pos.end.base) + 'del' + genomic_gap_fill_variant_alt.posedit.edit.ref + 'ins' + genomic_gap_fill_variant_alt.posedit.edit.ref + genomic_gap_fill_variant_alt.posedit.edit.ref
                                                genomic_gap_fill_variant_alt = hp.parse_hgvs_variant(
                                                    genomic_gap_fill_variant_alt_delins_from_dup)

                                        # Correct insertion alts
                                        if genomic_gap_fill_variant_alt.posedit.edit.type == 'ins':
                                            append_ref = sf.fetch_seq(genomic_gap_fill_variant_alt.ac,
                                                                      genomic_gap_fill_variant_alt.posedit.pos.start.base - 1,
                                                                      genomic_gap_fill_variant_alt.posedit.pos.end.base)
                                            genomic_gap_fill_variant_alt.posedit.edit.alt = append_ref[
                                                                                                0] + genomic_gap_fill_variant_alt.posedit.edit.alt + \
                                                                                            append_ref[1]

                                        # Split the reference and replacing alt sequence into a dictionary
                                        reference_bases = list(genomic_gap_fill_variant.posedit.edit.ref)
                                        if genomic_gap_fill_variant_alt.posedit.edit.alt is not None:
                                            alternate_bases = list(genomic_gap_fill_variant_alt.posedit.edit.alt)
                                        else:
                                            # Deletions with no ins
                                            pre_alternate_bases = list(genomic_gap_fill_variant_alt.posedit.edit.ref)
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
                                        for integer in range(genomic_gap_fill_variant_alt.posedit.pos.start.base,
                                                             genomic_gap_fill_variant_alt.posedit.pos.end.base + 1, 1):
                                            if integer == alt_start:
                                                alt_base_dict[integer] = str(''.join(alternate_bases))
                                            else:
                                                alt_base_dict[integer] = 'X'

                                        # Generate the alt sequence
                                        alternate_sequence_bases = []
                                        for integer in range(genomic_gap_fill_variant.posedit.pos.start.base,
                                                             genomic_gap_fill_variant.posedit.pos.end.base + 1, 1):
                                            if integer in alt_base_dict.keys():
                                                alternate_sequence_bases.append(alt_base_dict[integer])
                                            else:
                                                alternate_sequence_bases.append(ref_base_dict[integer])
                                        alternate_sequence = ''.join(alternate_sequence_bases)
                                        alternate_sequence = alternate_sequence.replace('X', '')

                                        # Add the new alt to the gap fill variant and generate transcript variant
                                        genomic_gap_fill_variant.posedit.edit.alt = alternate_sequence
                                        hgvs_refreshed_variant = vm.g_to_t(genomic_gap_fill_variant,
                                                                           tx_gap_fill_variant.ac)

                                        # Set warning
                                        gap_size = str(len(genomic_gap_fill_variant.posedit.edit.ref) - 2)
                                        disparity_deletion_in[1] = [gap_size]
                                        auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                                            stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + gap_size + ' genomic base(s) that fail to align to transcript ' + str(
                                            tx_hgvs_not_delins.ac)
                                        non_valid_caution = 'true'

                                        # Alignment position
                                        for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                        if re.match('NM_', str(for_location_c)):
                                            for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                        if re.match('\-', str(for_location_c.posedit.pos.start.offset)):
                                            gps = for_location_c.posedit.pos.start.base - 1
                                            gpe = for_location_c.posedit.pos.start.base
                                        else:
                                            gps = for_location_c.posedit.pos.start.base
                                            gpe = for_location_c.posedit.pos.start.base + 1
                                        gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                                        auto_info = auto_info + '%s' % (gap_position)

                                    else:
                                        if tx_hgvs_not_delins.posedit.pos.start.offset == 0 and tx_hgvs_not_delins.posedit.pos.end.offset == 0:
                                            # In this instance, we have identified a transcript gap but the n. version of
                                            # the transcript variant but do not have a position which actually hits the gap,
                                            # so the variant likely spans the gap, and is not picked up by an offset.
                                            try:
                                                c1 = vm.n_to_c(tx_hgvs_not_delins)
                                            except:
                                                c1 = tx_hgvs_not_delins
                                            g1 = nr_vm.t_to_g(c1, hgvs_genomic.ac)
                                            g3 = nr_vm.t_to_g(c1, hgvs_genomic.ac)
                                            g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                            ng2 = hn.normalize(g2)
                                            g3.posedit.pos.end.base = g3.posedit.pos.start.base + (
                                                    len(g3.posedit.edit.ref) - 1)
                                            try:
                                                c2 = vm.g_to_t(g3, c1.ac)
                                                if c2.posedit.pos.start.offset == 0 and c2.posedit.pos.end.offset == 0:
                                                    pass
                                                else:
                                                    tx_hgvs_not_delins = c2
                                                    try:
                                                        tx_hgvs_not_delins = vm.c_to_n(tx_hgvs_not_delins)
                                                    except hgvs.exceptions.HGVSError:
                                                        exceptPass()
                                            except hgvs.exceptions.HGVSInvalidVariantError:
                                                exceptPass()

                                        if re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                                                '\+', str(tx_hgvs_not_delins.posedit.pos.end)):
                                            auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                                                stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                                                disparity_deletion_in[
                                                    1]) + ' genomic base(s) that fail to align to transcript ' + str(
                                                tx_hgvs_not_delins.ac)
                                            non_valid_caution = 'true'
                                            try:
                                                c2 = vm.n_to_c(tx_hgvs_not_delins)
                                            except:
                                                c2 = tx_hgvs_not_delins
                                            c1 = copy.deepcopy(c2)
                                            c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                                            c1.posedit.pos.start.offset = 0
                                            c1.posedit.pos.end = c2.posedit.pos.start
                                            c1.posedit.edit.ref = ''
                                            c1.posedit.edit.alt = ''
                                            if orientation != -1:
                                                g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                g1.posedit.edit.alt = g1.posedit.edit.ref
                                            else:
                                                g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                g2.posedit.edit.alt = g2.posedit.edit.ref
                                            reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                            alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                            g3 = copy.deepcopy(g1)
                                            g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                            g3.posedit.edit.ref = reference
                                            g3.posedit.edit.alt = alternate
                                            c3 = vm.g_to_t(g3, c1.ac)
                                            hgvs_refreshed_variant = c3
                                            # Alignment position
                                            for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                            if re.match('NM_', str(for_location_c)):
                                                for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                                gps = for_location_c.posedit.pos.start.base
                                                gpe = for_location_c.posedit.pos.start.base + 1
                                            gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                                            # Warn update
                                            auto_info = auto_info + '%s' % (gap_position)
                                        elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search(
                                                '\+', str(tx_hgvs_not_delins.posedit.pos.start)):
                                            auto_info = auto_info + 'Genome position ' + str(
                                                stored_hgvs_not_delins.ac) + ':g.' + str(
                                                stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(
                                                disparity_deletion_in[1]) + '-bp gap in transcript ' + str(
                                                tx_hgvs_not_delins.ac)
                                            gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
                                            non_valid_caution = 'true'
                                            try:
                                                c1 = vm.n_to_c(tx_hgvs_not_delins)
                                            except:
                                                c1 = tx_hgvs_not_delins
                                            c2 = copy.deepcopy(c1)
                                            c2.posedit.pos.start = c1.posedit.pos.end
                                            c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                                            c2.posedit.pos.end.offset = 0
                                            c2.posedit.edit.ref = ''
                                            c2.posedit.edit.alt = ''
                                            if orientation != -1:
                                                g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                g2.posedit.edit.alt = g2.posedit.edit.ref
                                            else:
                                                g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                g1.posedit.edit.alt = g1.posedit.edit.ref
                                            reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                            alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                            g3 = copy.deepcopy(g1)
                                            g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                            g3.posedit.edit.ref = reference
                                            g3.posedit.edit.alt = alternate
                                            c3 = vm.g_to_t(g3, c1.ac)
                                            hgvs_refreshed_variant = c3
                                            # Alignment position
                                            for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                            if re.match('NM_', str(for_location_c)):
                                                for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                            gps = for_location_c.posedit.pos.end.base
                                            gpe = for_location_c.posedit.pos.end.base + 1
                                            gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                                            # Warn update
                                            auto_info = auto_info + '%s' % (gap_position)
                                        elif re.search('\-',
                                                       str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                                            '\-', str(tx_hgvs_not_delins.posedit.pos.end)):
                                            auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                                                stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                                                disparity_deletion_in[
                                                    1]) + ' genomic base(s) that fail to align to transcript ' + str(
                                                tx_hgvs_not_delins.ac)
                                            non_valid_caution = 'true'
                                            try:
                                                c2 = vm.n_to_c(tx_hgvs_not_delins)
                                            except:
                                                c2 = tx_hgvs_not_delins
                                            c1 = copy.deepcopy(c2)
                                            c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                                            c1.posedit.pos.start.offset = 0
                                            c1.posedit.pos.end = c2.posedit.pos.start
                                            c1.posedit.edit.ref = ''
                                            c1.posedit.edit.alt = ''
                                            if orientation != -1:
                                                g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                g1.posedit.edit.alt = g1.posedit.edit.ref
                                            else:
                                                g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                g2.posedit.edit.alt = g2.posedit.edit.ref
                                            reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                            alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                            g3 = copy.deepcopy(g1)
                                            g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                            g3.posedit.edit.ref = reference
                                            g3.posedit.edit.alt = alternate
                                            c3 = vm.g_to_t(g3, c1.ac)
                                            hgvs_refreshed_variant = c3
                                            # Alignment position
                                            for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                            if re.match('NM_', str(for_location_c)):
                                                for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                            gps = for_location_c.posedit.pos.start.base - 1
                                            gpe = for_location_c.posedit.pos.start.base
                                            gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                                            # Warn update
                                            auto_info = auto_info + '%s' % (gap_position)
                                        elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search(
                                                '\-', str(tx_hgvs_not_delins.posedit.pos.start)):
                                            auto_info = auto_info + 'Genome position ' + str(
                                                stored_hgvs_not_delins.ac) + ':g.' + str(
                                                stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(
                                                disparity_deletion_in[1]) + '-bp gap in transcript ' + str(
                                                tx_hgvs_not_delins.ac)
                                            gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
                                            non_valid_caution = 'true'
                                            try:
                                                c1 = vm.n_to_c(tx_hgvs_not_delins)
                                            except:
                                                c1 = tx_hgvs_not_delins
                                            c2 = copy.deepcopy(c1)
                                            c2.posedit.pos.start = c1.posedit.pos.end
                                            c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                                            c2.posedit.pos.end.offset = 0
                                            c2.posedit.edit.ref = ''
                                            c2.posedit.edit.alt = ''
                                            if orientation != -1:
                                                g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                g2.posedit.edit.alt = g2.posedit.edit.ref
                                            else:
                                                g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                g1.posedit.edit.alt = g1.posedit.edit.ref
                                            reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                            alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                            g3 = copy.deepcopy(g1)
                                            g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                            g3.posedit.edit.ref = reference
                                            g3.posedit.edit.alt = alternate
                                            c3 = vm.g_to_t(g3, c1.ac)
                                            hgvs_refreshed_variant = c3
                                            # Alignment position
                                            for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                                            if re.match('NM_', str(for_location_c)):
                                                for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                                            gps = for_location_c.posedit.pos.end.base - 1
                                            gpe = for_location_c.posedit.pos.end.base
                                            gap_position = ' between positions c.' + str(gps) + '_' + str(gpe) + '\n'
                                            # Warn update
                                            auto_info = auto_info + '%s' % (gap_position)
                                        else:
                                            auto_info = auto_info + str(stored_hgvs_not_delins.ac) + ':g.' + str(
                                                stored_hgvs_not_delins.posedit.pos) + ' contains ' + str(
                                                disparity_deletion_in[
                                                    1]) + ' genomic base(s) that fail to align to transcript ' + str(
                                                tx_hgvs_not_delins.ac) + '\n'
                                            hgvs_refreshed_variant = tx_hgvs_not_delins
                                            gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)

                                # GAP IN THE CHROMOSOME

                                elif disparity_deletion_in[0] == 'chromosome':
                                    # Set warning variables
                                    gap_position = ''
                                    gapped_alignment_warning = str(
                                        hgvs_genomic_5pr) + ' does not represent a true variant because it is an artefact of aligning the transcripts listed below with genome build ' + primary_assembly
                                    hgvs_refreshed_variant = tx_hgvs_not_delins
                                    # Warn
                                    auto_info = auto_info + str(hgvs_refreshed_variant.ac) + ':c.' + str(
                                        hgvs_refreshed_variant.posedit.pos) + ' contains ' + str(disparity_deletion_in[
                                                                                                     1]) + ' transcript base(s) that fail to align to chromosome ' + str(
                                        hgvs_genomic.ac) + '\n'
                                    gapped_transcripts = gapped_transcripts + str(hgvs_refreshed_variant.ac) + ' '
                                else:
                                    # Keep the same by re-setting rel_var
                                    hgvs_refreshed_variant = saved_hgvs_coding

                                # Edit the output
                                if re.match('NM_', str(hgvs_refreshed_variant.ac)) and not re.search('c', str(
                                        hgvs_refreshed_variant.type)):
                                    hgvs_refreshed_variant = evm.n_to_c(hgvs_refreshed_variant)
                                else:
                                    pass
                                try:
                                    hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
                                    if hgvs_refreshed_variant.posedit.edit.type == 'delins' and \
                                            hgvs_refreshed_variant.posedit.edit.ref[-1] == \
                                            hgvs_refreshed_variant.posedit.edit.alt[-1]:
                                        hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[
                                                                                  0:-1]
                                        hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[
                                                                                  0:-1]
                                        hgvs_refreshed_variant.posedit.pos.end.base = hgvs_refreshed_variant.posedit.pos.end.base - 1
                                        hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
                                    elif hgvs_refreshed_variant.posedit.edit.type == 'delins' and \
                                            hgvs_refreshed_variant.posedit.edit.ref[0] == \
                                            hgvs_refreshed_variant.posedit.edit.alt[0]:
                                        hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[
                                                                                  1:]
                                        hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[
                                                                                  1:]
                                        hgvs_refreshed_variant.posedit.pos.start.base = hgvs_refreshed_variant.posedit.pos.start.base + 1
                                        hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
                                except Exception as e:
                                    error = str(e)
                                    # Ensure the final variant is not intronic nor does it cross exon boundaries
                                    if re.match('Normalization of intronic variants is not supported',
                                                error) or re.match(
                                        'Unsupported normalization of variants spanning the exon-intron boundary',
                                        error):
                                        hgvs_refreshed_variant = saved_hgvs_coding
                                    else:
                                        pass

                                # Sort out equality to equality c. events where the code will add 2 additional bases
                                if hgvs_coding.posedit.edit.type == 'identity' and hgvs_refreshed_variant.posedit.edit.type == 'identity':  # and len(hgvs_refreshed_variant.posedit.edit.ref) ==  (len(hgvs_coding.posedit.edit.ref) + 2):
                                    pass
                                else:
                                    hgvs_coding = copy.deepcopy(hgvs_refreshed_variant)
                                coding = valstr(hgvs_coding)
                                variant = coding

                        # OBTAIN THE RefSeqGene coordinates
                        # Attempt 1 = UTA
                        sequences_for_tx = hdp.get_tx_mapping_options(hgvs_coding.ac)
                        recovered_rsg = []

                        for sequence in sequences_for_tx:
                            if re.match('^NG_', sequence[1]):
                                recovered_rsg.append(sequence[1])
                        recovered_rsg.sort()
                        recovered_rsg.reverse()

                        if len(recovered_rsg) > 0 and 'NG_' in recovered_rsg[0]:
                            refseqgene_ac = recovered_rsg[0]
                        else:
                            refseqgene_ac = ''

                        # Given the difficulties with mapping to and from RefSeqGenes, we now solely rely on UTA
                        if refseqgene_ac != '':
                            hgvs_refseq = vm.t_to_g(hgvs_coding, refseqgene_ac)
                            # Normalize the RefSeqGene Variant to the correct position
                            try:
                                hgvs_refseq = hn.normalize(hgvs_refseq)
                            except Exception as e:
                                # if re.search('insertion length must be 1', error):
                                hgvs_refseq = 'RefSeqGene record not available'
                                refseq = 'RefSeqGene record not available'
                                hgvs_refseq_ac = 'RefSeqGene record not available'
                                pass
                            else:
                                refseq = valstr(hgvs_refseq)
                                hgvs_refseq_ac = hgvs_refseq.ac
                        else:
                            hgvs_refseq = 'RefSeqGene record not available'
                            refseq = 'RefSeqGene record not available'
                            hgvs_refseq_ac = 'RefSeqGene record not available'

                        # Predicted effect on protein
                        protein_dict = va_func.myc_to_p(hgvs_coding, evm, hdp, hp, hn, vm, sf, re_to_p=False)
                        if protein_dict['error'] == '':
                            hgvs_protein = protein_dict['hgvs_protein']
                            protein = str(hgvs_protein)
                        else:
                            error = protein_dict['error']
                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                            if error == 'Cannot identify an in-frame Termination codon in the variant mRNA sequence':
                                hgvs_protein = protein_dict['hgvs_protein']
                                protein = str(hgvs_protein)
                            else:
                                logger.error(error)
                                continue

                        # Gene orientation wrt genome
                        ori = va_func.tx_exons(tx_ac=hgvs_coding.ac, alt_ac=hgvs_genomic.ac,
                                               alt_aln_method=alt_aln_method, hdp=hdp)
                        ori = int(ori[0]['alt_strand'])

                        # Look for normalized variant options that do not match hgvs_coding
                        # boundary crossing normalization
                        # Re-Save the required variants
                        hgvs_seek_var = copy.deepcopy(hgvs_coding)
                        saved_hgvs_coding = copy.deepcopy(hgvs_coding)

                        if ori == -1:
                            # position genomic at its most 5 prime position
                            try:
                                query_genomic = reverse_normalizer.normalize(hgvs_genomic)
                            except:
                                query_genomic = hgvs_genomic
                            # Map to the transcript and test for movement
                            try:
                                hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
                            except hgvs.exceptions.HGVSError as e:
                                hgvs_seek_var = saved_hgvs_coding
                            else:
                                seek_var = valstr(hgvs_seek_var)
                                seek_ac = str(hgvs_seek_var.ac)
                            if saved_hgvs_coding.posedit.edit.type != hgvs_seek_var.posedit.edit.type:
                                rec_var = 'false'
                                hgvs_seek_var = saved_hgvs_coding
                                seek_var = valstr(hgvs_seek_var)
                                seek_ac = str(hgvs_seek_var.ac)
                            elif suppress_c_normalization == 'true':
                                rec_var = 'false'
                                hgvs_seek_var = saved_hgvs_coding
                                seek_var = valstr(hgvs_seek_var)
                                seek_ac = str(hgvs_seek_var.ac)
                            elif (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                                    saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                                    hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                                    saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
                                try:
                                    automap = valstr(saved_hgvs_coding) + ' normalized to ' + valstr(hgvs_seek_var)
                                    hgvs_coding = hgvs_seek_var
                                    coding = valstr(hgvs_coding)
                                    validation['warnings'] = validation['warnings'] + ': ' + automap
                                    rng = hn.normalize(query_genomic)
                                except NotImplementedError:
                                    pass
                                try:
                                    c_for_p = vm.g_to_t(rng, hgvs_coding.ac)
                                except hgvs.exceptions.HGVSInvalidIntervalError as e:
                                    c_for_p = seek_var
                                try:
                                    # Predicted effect on protein
                                    protein_dict = va_func.myc_to_p(c_for_p, evm, hdp, hp, hn, vm, sf, re_to_p=False)
                                    if protein_dict['error'] == '':
                                        hgvs_protein = protein_dict['hgvs_protein']
                                        protein = str(hgvs_protein)
                                    else:
                                        error = protein_dict['error']
                                        if error == 'Cannot identify an in-frame Termination codon in the variant mRNA sequence':
                                            hgvs_protein = protein_dict['hgvs_protein']
                                            validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                    # Replace protein description in vars table
                                    protein = str(hgvs_protein)
                                except NotImplementedError:
                                    exceptPass()
                            else:
                                # Double check protein position by normalize genomic, and normalize back to c. for normalize or not to normalize issue
                                coding = valstr(hgvs_coding)

                        elif ori != -1:
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
                            else:
                                seek_var = valstr(hgvs_seek_var)
                                seek_ac = str(hgvs_seek_var.ac)
                            if saved_hgvs_coding.posedit.edit.type != hgvs_seek_var.posedit.edit.type:
                                rec_var = 'false'
                                hgvs_seek_var = saved_hgvs_coding
                                seek_var = valstr(hgvs_seek_var)
                                seek_ac = str(hgvs_seek_var.ac)
                            elif suppress_c_normalization == 'true':
                                rec_var = 'false'
                                hgvs_seek_var = saved_hgvs_coding
                                seek_var = valstr(hgvs_seek_var)
                                seek_ac = str(hgvs_seek_var.ac)
                            elif (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                                    saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                                    hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                                    saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
                                try:
                                    automap = valstr(saved_hgvs_coding) + ' normalized to ' + valstr(hgvs_seek_var)
                                    hgvs_coding = hgvs_seek_var
                                    coding = valstr(hgvs_coding)
                                    validation['warnings'] = validation['warnings'] + ': ' + automap
                                except NotImplementedError:
                                    exceptPass()
                            else:
                                # Double check protein position by reverse_norm genomic, and normalize back to c. for normalize or not to normalize issue
                                coding = valstr(hgvs_coding)
                                rng = reverse_normalizer.normalize(query_genomic)
                                try:
                                    # Diagram where - = intron and E = Exon

                                    # 3 prime
                                    # ---------EEEEEEEEEEEEEEEEE-----------
                                    #          <
                                    # Result, normalize of new variant will baulk at intronic
                                    # 5 prime
                                    #                          <
                                    # Result, normalize of new variant will be happy
                                    c_for_p = vm.g_to_t(rng, hgvs_coding.ac)
                                    try:
                                        hn.normalize(c_for_p)
                                    except hgvs.exceptions.HGVSError as e:
                                        exceptPass()
                                    else:
                                        # hgvs_protein = va_func.protein(str(c_for_p), evm, hp)
                                        protein_dict = va_func.myc_to_p(c_for_p, evm, hdp, hp, hn, vm, sf,
                                                                        re_to_p=False)
                                        if protein_dict['error'] == '':
                                            hgvs_protein = protein_dict['hgvs_protein']
                                            protein = str(hgvs_protein)
                                        else:
                                            error = protein_dict['error']
                                            if error == 'Cannot identify an in-frame Termination codon in the variant mRNA sequence':
                                                hgvs_protein = protein_dict['hgvs_protein']
                                                validation['warnings'] = validation['warnings'] + ': ' + str(error)
                                        # Replace protein description in vars table
                                        protein = str(hgvs_protein)
                                except Exception:
                                    exceptPass()

                        # Check for up-to-date transcript version
                        updated_transcript_variant = 'None'
                        tx_id_info = hdp.get_tx_identity_info(hgvs_coding.ac)
                        uta_gene_symbol = tx_id_info[6]
                        tx_for_gene = hdp.get_tx_for_gene(uta_gene_symbol)
                        ac_root, ac_version = hgvs_coding.ac.split('.')
                        version_tracking = '0'
                        update = ''
                        for accession in tx_for_gene:
                            try:
                                if re.match(ac_root, accession[3]):
                                    query_version = accession[3].split('.')[1]
                                    if int(query_version) > int(ac_version) and int(query_version) > int(
                                            version_tracking):
                                        version_tracking = query_version
                                        update = accession[3]
                            except ValueError:
                                exceptPass()

                        if update != '':
                            hgvs_updated = copy.deepcopy(hgvs_coding)
                            hgvs_updated.ac = update
                            try:
                                vr.validate(hgvs_updated)
                            # Updated reference sequence
                            except hgvs.exceptions.HGVSError as e:
                                error = str(e)
                                if re.search('does not agree with reference sequence', str(error)):
                                    match = re.findall('\(([GATC]+)\)', error)
                                    new_ref = match[1]
                                    hgvs_updated.posedit.edit.ref = new_ref
                                    vr.validate(hgvs_updated)
                                    updated_transcript_variant = hgvs_updated
                                else:
                                    pass
                            updated_transcript_variant = hgvs_updated
                            validation['warnings'] = validation[
                                                         'warnings'] + ': ' + 'A more recent version of the selected reference sequence ' + hgvs_coding.ac + ' is available (' + updated_transcript_variant.ac + ')' + ': ' + str(
                                updated_transcript_variant) + ' MUST be fully validated prior to use in reports: select_variants=' + valstr(
                                updated_transcript_variant)

                # Set the data
                set_output_type_flag = 'gene'
                validation['description'] = hgnc_gene_info
                validation['coding'] = str(hgvs_coding)
                validation['genomic_r'] = str(hgvs_refseq)
                validation['genomic_g'] = str(hgvs_genomic)
                validation['protein'] = str(hgvs_protein)
                validation['primary_assembly'] = primary_assembly
                if gap_compensation is True:
                    validation['test_stash_tx_left'] = test_stash_tx_left
                    validation['test_stash_tx_right'] = test_stash_tx_right
                # finish timing
                logger.traceEnd(validation)
            # Report errors to User and VV admin
            except KeyboardInterrupt:
                raise
            except:
                set_output_type_flag = 'error'
                error = 'Validation error'
                validation['warnings'] = str(error)
                exc_type, exc_value, last_traceback = sys.exc_info()
                te = traceback.format_exc()
                tbk = [str(exc_type), str(exc_value), str(te)]
                er = str('\n'.join(tbk))
                logger.error(str(exc_type) + " " + str(exc_value))
                logger.debug(er)

                continue

        # Outside the for loop
        ######################
        logger.trace("End of for loop")
        # order the rows
        # from operator import itemgetter
        by_order = sorted(batch_list, key=itemgetter('order'))

        for valid in by_order:
            if 'write' in valid.keys():
                if valid['write'] == 'true':
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
                    warnings = valid['warnings']
                    warnings = re.sub('del[GATC][GATC][GATC][GATC]+', 'del', warnings)
                    warnings = re.sub('^: ', '', warnings)
                    warnings = re.sub('::', ':', warnings)

                    # Submitted variant
                    submitted = valid['id']

                    # Genomic sequence variation
                    genomic_variant = valid['genomic_g']

                    # genomic accession
                    if genomic_variant != '':
                        hgvs_genomic_variant = hp.parse_hgvs_variant(genomic_variant)
                        genomic_variant = valstr(hgvs_genomic_variant)
                        genomic_accession = hgvs_genomic_variant.ac
                    else:
                        genomic_accession = ''

                    # RefSeqGene variation
                    refseqgene_variant = valid['genomic_r']
                    refseqgene_variant = refseqgene_variant.strip()
                    if re.search('RefSeqGene', refseqgene_variant) or refseqgene_variant == '':
                        warnings = warnings + ': ' + refseqgene_variant
                        refseqgene_variant = ''
                        lrg_variant = ''
                        hgvs_refseqgene_variant = 'false'
                    else:
                        hgvs_refseqgene_variant = hp.parse_hgvs_variant(refseqgene_variant)
                        rsg_ac = va_dbCrl.data.get_lrgID_from_RefSeqGeneID(str(hgvs_refseqgene_variant.ac))
                        if rsg_ac[0] == 'none':
                            lrg_variant = ''
                        else:
                            hgvs_lrg = copy.deepcopy(hgvs_refseqgene_variant)
                            hgvs_lrg.ac = rsg_ac[0]
                            lrg_variant = valstr(hgvs_lrg)
                            if rsg_ac[1] == 'public':
                                pass
                            else:
                                warnings = warnings + ': The current status of ' + str(
                                    hgvs_lrg.ac) + ' is pending therefore changes may be made to the LRG reference sequence'

                    # Transcript sequence variation
                    tx_variant = valid['coding']
                    if tx_variant != '':
                        if '(' in tx_variant and ')' in tx_variant:
                            tx_variant = tx_variant.split('(')[1]
                            tx_variant = tx_variant.replace(')', '')

                        # transcript accession
                        hgvs_tx_variant = hp.parse_hgvs_variant(tx_variant)
                        tx_variant = valstr(hgvs_tx_variant)
                        hgvs_transcript_variant = hp.parse_hgvs_variant(tx_variant)
                        transcript_accession = hgvs_transcript_variant.ac

                        # Handle LRG
                        lrg_status = 'public'
                        lrg_transcript = va_dbCrl.data.get_lrgTranscriptID_from_RefSeqTranscriptID(transcript_accession)
                        if lrg_transcript == 'none':
                            lrg_transcript_variant = ''
                        else:
                            # Note - LRG availability is dependant on UTA containing the data. In some
                            # instances we will be able to display the LRG_tx without being able to
                            # display the LRG gene data

                            # if not re.search('RefSeqGene', refseqgene_variant) or refseqgene_variant != '':
                            # if hgvs_refseqgene_variant != 'RefSeqGene record not available' and hgvs_refseqgene_variant != 'false':
                            try:
                                hgvs_lrg_t = vm.g_to_t(hgvs_refseqgene_variant, transcript_accession)
                                hgvs_lrg_t.ac = lrg_transcript
                                lrg_transcript_variant = valstr(hgvs_lrg_t)
                            except:
                                if hgvs_transcript_variant.posedit.pos.start.offset == 0 and hgvs_transcript_variant.posedit.pos.end.offset == 0:
                                    hgvs_lrg_t = copy.copy(hgvs_transcript_variant)
                                    hgvs_lrg_t.ac = lrg_transcript
                                    lrg_transcript_variant = valstr(hgvs_lrg_t)
                                else:
                                    lrg_transcript_variant = ''
                    else:
                        transcript_accession = ''
                        lrg_transcript_variant = ''

                    # Look for intronic variants
                    if transcript_accession != '' and genomic_accession != '':
                        # Remove del bases
                        str_transcript = valstr(hgvs_transcript_variant)
                        hgvs_transcript_variant = hp.parse_hgvs_variant(str_transcript)
                        try:
                            vr.validate(hgvs_transcript_variant)
                        except hgvs.exceptions.HGVSError as e:
                            error = str(e)
                            if re.search('intronic variant', error):
                                genome_context_transcript_variant = genomic_accession + '(' + transcript_accession + '):c.' + str(
                                    hgvs_transcript_variant.posedit)
                                if refseqgene_variant != '':
                                    hgvs_refseqgene_variant = hp.parse_hgvs_variant(refseqgene_variant)
                                    refseqgene_accession = hgvs_refseqgene_variant.ac
                                    hgvs_coding_from_refseqgene = vm.g_to_t(hgvs_refseqgene_variant,
                                                                            hgvs_transcript_variant.ac)
                                    hgvs_coding_from_refseqgene = valstr(hgvs_coding_from_refseqgene)
                                    hgvs_coding_from_refseqgene = hp.parse_hgvs_variant(hgvs_coding_from_refseqgene)
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
                    predicted_protein_variant = valid['protein']
                    if re.match('NP_', predicted_protein_variant):
                        rs_p, pred_prot_posedit = predicted_protein_variant.split(':')
                        lrg_p = va_dbCrl.data.get_lrgProteinID_from_RefSeqProteinID(rs_p)
                        if re.match('LRG', lrg_p):
                            predicted_protein_variant = rs_p + '(' + lrg_p + '):' + pred_prot_posedit

                    # Gene
                    if transcript_accession != '':
                        try:
                            gene_symbol = va_dbCrl.data.get_gene_symbol_from_transcriptID(transcript_accession)
                        except:
                            gene_symbol = 'Unable to verify gene symbol for ' + str(transcript_accession)
                    else:
                        gene_symbol = ''

                    # Transcript description
                    transcript_description = valid['description']

                    # Stashed variants
                    if 'test_stash_tx_left' not in validation:
                        pass
                    else:
                        test_stash_tx_left = validation['test_stash_tx_left']
                    if 'test_stash_tx_right' not in validation:
                        pass
                    else:
                        test_stash_tx_right = validation['test_stash_tx_right']

                    # Multiple genomic variants
                    # multi_gen_vars = []
                    if tx_variant != '':
                        hgvs_coding = hp.parse_hgvs_variant(str(tx_variant))
                        # Gap gene black list
                        try:
                            gene_symbol = va_dbCrl.data.get_gene_symbol_from_transcriptID(hgvs_coding.ac)
                        except Exception:
                            exceptPass()
                        else:
                            # If the gene symbol is not in the list, the value False will be returned
                            gap_compensation = gapGenes.gap_black_list(gene_symbol)

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
                            exceptPass()

                        # Warn gap code status
                        logger.warning("gap_compensation_3 = " + str(gap_compensation))
                        multi_g = []
                        multi_list = []
                        mapping_options = hdp.get_tx_mapping_options(hgvs_coding.ac)
                        for alt_chr in mapping_options:
                            if (re.match('NC_', alt_chr[1]) or re.match('NT_', alt_chr[1]) or re.match('NW_',
                                                                                                       alt_chr[1])) and \
                                    alt_chr[2] == alt_aln_method:
                                multi_list.append(alt_chr[1])

                        for alt_chr in multi_list:
                            try:
                                # Re set ori
                                ori = va_func.tx_exons(tx_ac=hgvs_coding.ac, alt_ac=alt_chr,
                                                       alt_aln_method=alt_aln_method, hdp=hdp)
                                orientation = int(ori[0]['alt_strand'])
                                hgvs_alt_genomic = va_func.myvm_t_to_g(hgvs_coding, alt_chr, no_norm_evm, vm, hp, hn,
                                                                       sf, nr_vm)
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

                                    most_3pr_hgvs_genomic = va_func.myvm_t_to_g(chromosome_normalized_hgvs_coding,
                                                                                alt_chr,
                                                                                no_norm_evm, vm, hp, hn, sf, nr_vm)
                                    hgvs_genomic_possibilities.append(most_3pr_hgvs_genomic)

                                    # First to the right
                                    hgvs_stash = copy.deepcopy(hgvs_coding)
                                    try:
                                        hgvs_stash = no_norm_evm.c_to_n(hgvs_stash)
                                    except:
                                        exceptPass()
                                    try:
                                        stash_ac = hgvs_stash.ac
                                        stash_dict = va_H2V.hard_right_hgvs2vcf(hgvs_stash, primary_assembly, hn, sf)
                                        stash_pos = int(stash_dict['pos'])
                                        stash_ref = stash_dict['ref']
                                        stash_alt = stash_dict['alt']
                                        # Generate an end position
                                        stash_end = str(stash_pos + len(stash_ref) - 1)
                                        # make a not real deletion insertion
                                        stash_hgvs_not_delins = hp.parse_hgvs_variant(
                                            stash_ac + ':' + hgvs_stash.type + '.' + str(
                                                stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)
                                        try:
                                            stash_hgvs_not_delins = no_norm_evm.n_to_c(stash_hgvs_not_delins)
                                        except:
                                            exceptPass()
                                            # Store a tx copy for later use
                                        test_stash_tx_right = copy.deepcopy(stash_hgvs_not_delins)
                                        stash_genomic = va_func.myvm_t_to_g(test_stash_tx_right, hgvs_alt_genomic.ac,
                                                                            no_norm_evm, vm, hp, hn, sf, nr_vm)
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
                                            hgvs_reform_ident = hp.parse_hgvs_variant(reform_ident)
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
                                        exceptPass()
                                    except ValueError:
                                        exceptPass()

                                    # Then to the left
                                    hgvs_stash = copy.deepcopy(hgvs_coding)
                                    try:
                                        hgvs_stash = no_norm_evm.c_to_n(hgvs_stash)
                                    except:
                                        exceptPass()
                                    try:
                                        stash_ac = hgvs_stash.ac
                                        stash_dict = va_H2V.hard_left_hgvs2vcf(hgvs_stash, primary_assembly,
                                                                               reverse_normalizer, sf)
                                        stash_pos = int(stash_dict['pos'])
                                        stash_ref = stash_dict['ref']
                                        stash_alt = stash_dict['alt']
                                        # Generate an end position
                                        stash_end = str(stash_pos + len(stash_ref) - 1)
                                        # make a not real deletion insertion
                                        stash_hgvs_not_delins = hp.parse_hgvs_variant(
                                            stash_ac + ':' + hgvs_stash.type + '.' + str(
                                                stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)
                                        try:
                                            stash_hgvs_not_delins = no_norm_evm.n_to_c(stash_hgvs_not_delins)
                                        except:
                                            exceptPass()
                                            # Store a tx copy for later use
                                        test_stash_tx_left = copy.deepcopy(stash_hgvs_not_delins)
                                        stash_genomic = va_func.myvm_t_to_g(test_stash_tx_left, hgvs_alt_genomic.ac,
                                                                            no_norm_evm, vm, hp, hn, sf, nr_vm)
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
                                            hgvs_reform_ident = hp.parse_hgvs_variant(reform_ident)
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
                                        exceptPass()
                                    except ValueError:
                                        exceptPass()

                                    # direct mapping from reverse_normalized transcript insertions in the delins format
                                    try:
                                        if hgvs_coding.posedit.edit.type == 'ins':
                                            most_5pr_hgvs_transcript_variant = copy.deepcopy(hgvs_coding)
                                            most_3pr_hgvs_transcript_variant = reverse_normalizer.normalize(hgvs_coding)
                                            try:
                                                n_3pr = vm.c_to_n(most_3pr_hgvs_transcript_variant)
                                                n_5pr = vm.c_to_n(most_5pr_hgvs_transcript_variant)
                                            except:
                                                n_3pr = most_3pr_hgvs_transcript_variant
                                                n_5pr = most_5pr_hgvs_transcript_variant
                                            # Make into a delins by adding the ref bases to the variant ref and alt
                                            pr3_ref = sf.fetch_seq(hgvs_coding.ac, n_3pr.posedit.pos.start.base - 1,
                                                                   n_3pr.posedit.pos.end.base)
                                            pr5_ref = sf.fetch_seq(hgvs_coding.ac, n_5pr.posedit.pos.start.base - 1,
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
                                            genomic_from_most_3pr_hgvs_transcript_variant = vm.t_to_g(
                                                most_3pr_hgvs_transcript_variant, hgvs_genomic.ac)
                                            genomic_from_most_5pr_hgvs_transcript_variant = vm.t_to_g(
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
                                                    genomic_from_most_3pr_hgvs_transcript_variant = hp.parse_hgvs_variant(
                                                        genomic_from_most_3pr_hgvs_transcript_variant_delins_from_dup)

                                            try:
                                                if most_3pr_hgvs_transcript_variant.posedit.edit.alt is None:
                                                    most_3pr_hgvs_transcript_variant.posedit.edit.alt = ''
                                            except Exception as e:
                                                if str(e) == "'Dup' object has no attribute 'alt'":
                                                    most_3pr_hgvs_transcript_variant_delins_from_dup = most_3pr_hgvs_transcript_variant.ac + ':' + most_3pr_hgvs_transcript_variant.type + '.' + str(
                                                        most_3pr_hgvs_transcript_variant.posedit.pos.start.base) + '_' + str(
                                                        most_3pr_hgvs_transcript_variant.posedit.pos.end.base) + 'del' + most_3pr_hgvs_transcript_variant.posedit.edit.ref + 'ins' + most_3pr_hgvs_transcript_variant.posedit.edit.ref + most_3pr_hgvs_transcript_variant.posedit.edit.ref
                                                    most_3pr_hgvs_transcript_variant = hp.parse_hgvs_variant(
                                                        most_3pr_hgvs_transcript_variant_delins_from_dup)

                                            try:
                                                if genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.alt is None:
                                                    genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.alt = ''
                                            except Exception as e:
                                                if str(e) == "'Dup' object has no attribute 'alt'":
                                                    genomic_from_most_5pr_hgvs_transcript_variant_delins_from_dup = genomic_from_most_5pr_hgvs_transcript_variant.ac + ':' + genomic_from_most_5pr_hgvs_transcript_variant.type + '.' + str(
                                                        genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.start.base) + '_' + str(
                                                        genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.end.base) + 'del' + genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.ref + 'ins' + genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.ref + genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.ref
                                                    genomic_from_most_5pr_hgvs_transcript_variant = hp.parse_hgvs_variant(
                                                        genomic_from_most_5pr_hgvs_transcript_variant_delins_from_dup)

                                            try:
                                                if most_5pr_hgvs_transcript_variant.posedit.edit.alt is None:
                                                    most_5pr_hgvs_transcript_variant.posedit.edit.alt = ''
                                            except Exception as e:
                                                if str(e) == "'Dup' object has no attribute 'alt'":
                                                    most_5pr_hgvs_transcript_variant_delins_from_dup = most_5pr_hgvs_transcript_variant.ac + ':' + most_5pr_hgvs_transcript_variant.type + '.' + str(
                                                        most_5pr_hgvs_transcript_variant.posedit.pos.start.base) + '_' + str(
                                                        most_5pr_hgvs_transcript_variant.posedit.pos.end.base) + 'del' + most_5pr_hgvs_transcript_variant.posedit.edit.ref + 'ins' + most_5pr_hgvs_transcript_variant.posedit.edit.ref + most_5pr_hgvs_transcript_variant.posedit.edit.ref
                                                    most_5pr_hgvs_transcript_variant = hp.parse_hgvs_variant(
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
                                        exceptPass()

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
                                                    lhb = sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                                                    rhb = sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                                                    hgvs_genomic.posedit.edit.ref = lhb + rhb
                                                    hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                                                    hgvs_genomic.posedit.pos.start.base = end
                                                    hgvs_genomic.posedit.pos.end.base = start
                                                    reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(
                                                        hgvs_genomic)
                                                if hgvs_genomic.posedit.edit.type == 'del':
                                                    start = hgvs_genomic.posedit.pos.start.base
                                                    end = hgvs_genomic.posedit.pos.end.base
                                                    lhb = sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                                                    rhb = sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
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
                                                    ref_bases = sf.fetch_seq(str(hgvs_genomic.ac), start - 1, end)
                                                    lhb = sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                                                    rhb = sf.fetch_seq(str(hgvs_genomic.ac), start, end)
                                                    hgvs_genomic.posedit.edit.ref = lhb + rhb
                                                    hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                                                    reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(
                                                        hgvs_genomic)

                                        hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)
                                        # Store a copy for later use
                                        stored_hgvs_genomic_5pr = copy.deepcopy(hgvs_genomic_5pr)

                                        # Make VCF
                                        vcf_dict = va_H2V.hgvs2vcf(reverse_normalized_hgvs_genomic, primary_assembly,
                                                                   reverse_normalizer, sf)
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
                                        stored_hgvs_not_delins = hp.parse_hgvs_variant(str(
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
                                                seek_var = valstr(hgvs_seek_var)
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
                                            seek_var = valstr(hgvs_seek_var)
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
                                            if re.search('\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\d+\-',
                                                                                                               str(
                                                                                                                   hgvs_seek_var.posedit.pos)) or re.search(
                                                '\*\d+\+', str(
                                                    hgvs_seek_var.posedit.pos)) or re.search('\*\d+\-', str(
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
                                                                ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start,
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
                                                                ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start,
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
                                                                ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start,
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
                                                                ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start,
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
                                                if re.search('\+',
                                                             str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search(
                                                    '\+',
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
                                                        exceptPass()

                                                elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                                                    # move tx end base to next available non-offset base
                                                    rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
                                                    rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                                                    rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                                    if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                                        test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                                    else:
                                                        test_tx_var = rn_tx_hgvs_not_delins
                                                    # re-make genomic and tx
                                                    hgvs_not_delins = va_func.myvm_t_to_g(test_tx_var, alt_chr,
                                                                                          no_norm_evm, vm, hp, hn, sf,
                                                                                          nr_vm)
                                                    rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                               str(
                                                                                                   saved_hgvs_coding.ac))
                                                elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                                                    # move tx start base to previous available non-offset base
                                                    rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                                    rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                                    if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                                        test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                                    else:
                                                        test_tx_var = rn_tx_hgvs_not_delins
                                                    # re-make genomic and tx
                                                    hgvs_not_delins = va_func.myvm_t_to_g(test_tx_var, alt_chr,
                                                                                          no_norm_evm, vm, hp, hn, sf,
                                                                                          nr_vm)
                                                    rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                               str(
                                                                                                   saved_hgvs_coding.ac))
                                                    rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                                #                                                 else:
                                                #                                                     pass

                                                # Check for -ve base and adjust
                                                elif re.search('\-',
                                                               str(
                                                                   rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search(
                                                    '\-',
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
                                                        exceptPass()
                                                elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                                                    # move tx end base back to next available non-offset base
                                                    rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                                                    # Delete the ref
                                                    rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                                    # Add the additional base to the ALT
                                                    start = rn_tx_hgvs_not_delins.posedit.pos.end.base - 1
                                                    end = rn_tx_hgvs_not_delins.posedit.pos.end.base
                                                    ref_bases = sf.fetch_seq(str(tx_hgvs_not_delins.ac), start, end)
                                                    rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + ref_bases
                                                    if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                                        test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                                    else:
                                                        test_tx_var = rn_tx_hgvs_not_delins
                                                    # re-make genomic and tx
                                                    hgvs_not_delins = va_func.myvm_t_to_g(test_tx_var, alt_chr,
                                                                                          no_norm_evm, vm, hp, hn, sf,
                                                                                          nr_vm)
                                                    rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                               str(
                                                                                                   saved_hgvs_coding.ac))
                                                elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                                                    # move tx start base to previous available non-offset base
                                                    rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                                    rn_tx_hgvs_not_delins.posedit.pos.start.base = rn_tx_hgvs_not_delins.posedit.pos.start.base - 1
                                                    rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                                                    if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                                                        test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                                                    else:
                                                        test_tx_var = rn_tx_hgvs_not_delins
                                                    # re-make genomic and tx
                                                    hgvs_not_delins = va_func.myvm_t_to_g(test_tx_var, alt_chr,
                                                                                          no_norm_evm, vm, hp, hn, sf,
                                                                                          nr_vm)
                                                    rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                               str(
                                                                                                   saved_hgvs_coding.ac))
                                                    rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                                                else:
                                                    exceptPass()

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
                                                        hgvs_t_possibility = vm.g_to_t(possibility, hgvs_coding.ac)
                                                        if hgvs_t_possibility.posedit.edit.type == 'ins':
                                                            try:
                                                                hgvs_t_possibility = vm.c_to_n(hgvs_t_possibility)
                                                            except:
                                                                continue
                                                            if hgvs_t_possibility.posedit.pos.start.offset != 0 or hgvs_t_possibility.posedit.pos.end.offset != 0:
                                                                continue
                                                            ins_ref = sf.fetch_seq(hgvs_t_possibility.ac,
                                                                                   hgvs_t_possibility.posedit.pos.start.base - 1,
                                                                                   hgvs_t_possibility.posedit.pos.start.base + 1)
                                                            try:
                                                                hgvs_t_possibility = vm.n_to_c(hgvs_t_possibility)
                                                            except:
                                                                continue
                                                            hgvs_t_possibility.posedit.edit.ref = ins_ref
                                                            hgvs_t_possibility.posedit.edit.alt = ins_ref[
                                                                                                      0] + hgvs_t_possibility.posedit.edit.alt + \
                                                                                                  ins_ref[1]
                                                        if possibility.posedit.edit.type == 'ins':
                                                            ins_ref = sf.fetch_seq(possibility.ac,
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
                                                            tx_hgvs_not_delins = vm.c_to_n(re_capture_tx_variant[2])
                                                        except:
                                                            tx_hgvs_not_delins = re_capture_tx_variant[2]
                                                        disparity_deletion_in = re_capture_tx_variant[0:-1]
                                                    else:
                                                        pass

                                            # Final sanity checks
                                            try:
                                                vm.g_to_t(hgvs_not_delins, tx_hgvs_not_delins.ac)
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
                                                rtx = vm.g_to_t(rg, tx_hgvs_not_delins.ac)
                                                fg = hn.normalize(hgvs_not_delins)
                                                ftx = vm.g_to_t(fg, tx_hgvs_not_delins.ac)
                                                if (
                                                        rtx.posedit.pos.start.offset == 0 and rtx.posedit.pos.end.offset == 0) and (
                                                        ftx.posedit.pos.start.offset != 0 and ftx.posedit.pos.end.offset != 0):
                                                    exons = hdp.get_tx_exons(ftx.ac, hgvs_not_delins.ac, alt_aln_method)
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
                                                            tx_hgvs_not_delins = vm.c_to_n(ftx)
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
                                                    tx_hgvs_not_delins = hp.parse_hgvs_variant(
                                                        tx_hgvs_not_delins_delins_from_dup)

                                            if disparity_deletion_in[0] == 'transcript':
                                                # amend_RefSeqGene = 'true'
                                                # ANY VARIANT WHOLLY WITHIN THE GAP
                                                if (re.search('\+',
                                                              str(tx_hgvs_not_delins.posedit.pos.start)) or re.search(
                                                    '\-', str(tx_hgvs_not_delins.posedit.pos.start))) and (
                                                        re.search('\+',
                                                                  str(tx_hgvs_not_delins.posedit.pos.end)) or re.search(
                                                    '\-', str(tx_hgvs_not_delins.posedit.pos.end))):
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
                                                            tx_gap_fill_variant = hp.parse_hgvs_variant(
                                                                tx_gap_fill_variant_delins_from_dup)

                                                    # Identify which half of the NOT-intron the start position of the variant is in
                                                    if re.search('\-', str(tx_gap_fill_variant.posedit.pos.start)):
                                                        tx_gap_fill_variant.posedit.pos.start.base = tx_gap_fill_variant.posedit.pos.start.base - 1
                                                        tx_gap_fill_variant.posedit.pos.start.offset = int(
                                                            '0')  # int('+1')
                                                        tx_gap_fill_variant.posedit.pos.end.offset = int(
                                                            '0')  # int('-1')
                                                        tx_gap_fill_variant.posedit.edit.alt = ''
                                                        tx_gap_fill_variant.posedit.edit.ref = ''
                                                    elif re.search('\+', str(tx_gap_fill_variant.posedit.pos.start)):
                                                        tx_gap_fill_variant.posedit.pos.start.offset = int(
                                                            '0')  # int('+1')
                                                        tx_gap_fill_variant.posedit.pos.end.base = tx_gap_fill_variant.posedit.pos.end.base + 1
                                                        tx_gap_fill_variant.posedit.pos.end.offset = int(
                                                            '0')  # int('-1')
                                                        tx_gap_fill_variant.posedit.edit.alt = ''
                                                        tx_gap_fill_variant.posedit.edit.ref = ''

                                                    try:
                                                        tx_gap_fill_variant = vm.n_to_c(tx_gap_fill_variant)
                                                    except:
                                                        exceptPass()
                                                    genomic_gap_fill_variant = vm.t_to_g(tx_gap_fill_variant,
                                                                                         reverse_normalized_hgvs_genomic.ac)
                                                    genomic_gap_fill_variant.posedit.edit.alt = genomic_gap_fill_variant.posedit.edit.ref

                                                    try:
                                                        c_tx_hgvs_not_delins = vm.n_to_c(tx_hgvs_not_delins)
                                                    except Exception:
                                                        c_tx_hgvs_not_delins = copy.copy(tx_hgvs_not_delins)
                                                    genomic_gap_fill_variant_alt = vm.t_to_g(c_tx_hgvs_not_delins,
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
                                                            genomic_gap_fill_variant = hp.parse_hgvs_variant(
                                                                genomic_gap_fill_variant_delins_from_dup)
                                                            genomic_gap_fill_variant_alt_delins_from_dup = genomic_gap_fill_variant_alt.ac + ':' + genomic_gap_fill_variant_alt.type + '.' + str(
                                                                genomic_gap_fill_variant_alt.posedit.pos.start.base) + '_' + str(
                                                                genomic_gap_fill_variant_alt.posedit.pos.end.base) + 'del' + genomic_gap_fill_variant_alt.posedit.edit.ref + 'ins' + genomic_gap_fill_variant_alt.posedit.edit.ref + genomic_gap_fill_variant_alt.posedit.edit.ref
                                                            genomic_gap_fill_variant_alt = hp.parse_hgvs_variant(
                                                                genomic_gap_fill_variant_alt_delins_from_dup)

                                                    # Correct insertion alts
                                                    if genomic_gap_fill_variant_alt.posedit.edit.type == 'ins':
                                                        append_ref = sf.fetch_seq(genomic_gap_fill_variant_alt.ac,
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
                                                        if integer in alt_base_dict.keys():
                                                            alternate_sequence_bases.append(alt_base_dict[integer])
                                                        else:
                                                            alternate_sequence_bases.append(ref_base_dict[integer])
                                                    alternate_sequence = ''.join(alternate_sequence_bases)
                                                    alternate_sequence = alternate_sequence.replace('X', '')

                                                    # Add the new alt to the gap fill variant and generate transcript variant
                                                    genomic_gap_fill_variant.posedit.edit.alt = alternate_sequence
                                                    hgvs_refreshed_variant = vm.g_to_t(genomic_gap_fill_variant,
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
                                                    if re.match('\-', str(for_location_c.posedit.pos.start.offset)):
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
                                                            c1 = vm.n_to_c(tx_hgvs_not_delins)
                                                        except:
                                                            c1 = tx_hgvs_not_delins
                                                        g1 = nr_vm.t_to_g(c1, hgvs_genomic.ac)
                                                        g3 = nr_vm.t_to_g(c1, hgvs_genomic.ac)
                                                        g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                        ng2 = hn.normalize(g2)
                                                        g3.posedit.pos.end.base = g3.posedit.pos.start.base + (
                                                                len(g3.posedit.edit.ref) - 1)
                                                        try:
                                                            c2 = vm.g_to_t(g3, c1.ac)
                                                            if c2.posedit.pos.start.offset == 0 and c2.posedit.pos.end.offset == 0:
                                                                pass
                                                            else:
                                                                tx_hgvs_not_delins = c2
                                                                try:
                                                                    tx_hgvs_not_delins = vm.c_to_n(tx_hgvs_not_delins)
                                                                except hgvs.exceptions.HGVSError:
                                                                    exceptPass()
                                                        except hgvs.exceptions.HGVSInvalidVariantError:
                                                            exceptPass()

                                                    if re.search('\+', str(
                                                            tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                                                        '\+',
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
                                                            c2 = vm.n_to_c(tx_hgvs_not_delins)
                                                        except:
                                                            c2 = tx_hgvs_not_delins
                                                        c1 = copy.deepcopy(c2)
                                                        c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                                                        c1.posedit.pos.start.offset = 0
                                                        c1.posedit.pos.end = c2.posedit.pos.start
                                                        c1.posedit.edit.ref = ''
                                                        c1.posedit.edit.alt = ''
                                                        if orientation != -1:
                                                            g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                            g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                            g1.posedit.edit.alt = g1.posedit.edit.ref
                                                        else:
                                                            g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                            g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                            g2.posedit.edit.alt = g2.posedit.edit.ref
                                                        reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                                        alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                                        g3 = copy.deepcopy(g1)
                                                        g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                                        g3.posedit.edit.ref = reference
                                                        g3.posedit.edit.alt = alternate
                                                        c3 = vm.g_to_t(g3, c1.ac)
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
                                                    elif re.search('\+', str(
                                                            tx_hgvs_not_delins.posedit.pos.end)) and not re.search('\+',
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
                                                            c1 = vm.n_to_c(tx_hgvs_not_delins)
                                                        except:
                                                            c1 = tx_hgvs_not_delins
                                                        c2 = copy.deepcopy(c1)
                                                        c2.posedit.pos.start = c1.posedit.pos.end
                                                        c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                                                        c2.posedit.pos.end.offset = 0
                                                        c2.posedit.edit.ref = ''
                                                        c2.posedit.edit.alt = ''
                                                        if orientation != -1:
                                                            g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                            g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                            g2.posedit.edit.alt = g2.posedit.edit.ref
                                                        else:
                                                            g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                            g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                            g1.posedit.edit.alt = g1.posedit.edit.ref
                                                        reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                                        alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                                        g3 = copy.deepcopy(g1)
                                                        g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                                        g3.posedit.edit.ref = reference
                                                        g3.posedit.edit.alt = alternate
                                                        c3 = vm.g_to_t(g3, c1.ac)
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
                                                    elif re.search('\-', str(
                                                            tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                                                        '\-',
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
                                                            c2 = vm.n_to_c(tx_hgvs_not_delins)
                                                        except:
                                                            c2 = tx_hgvs_not_delins
                                                        c1 = copy.deepcopy(c2)
                                                        c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                                                        c1.posedit.pos.start.offset = 0
                                                        c1.posedit.pos.end = c2.posedit.pos.start
                                                        c1.posedit.edit.ref = ''
                                                        c1.posedit.edit.alt = ''
                                                        if orientation != -1:
                                                            g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                            g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                            g1.posedit.edit.alt = g1.posedit.edit.ref
                                                        else:
                                                            g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                            g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                            g2.posedit.edit.alt = g2.posedit.edit.ref
                                                        reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                                        alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                                        g3 = copy.deepcopy(g1)
                                                        g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                                        g3.posedit.edit.ref = reference
                                                        g3.posedit.edit.alt = alternate
                                                        c3 = vm.g_to_t(g3, c1.ac)
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
                                                    elif re.search('\-', str(
                                                            tx_hgvs_not_delins.posedit.pos.end)) and not re.search('\-',
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
                                                            c1 = vm.n_to_c(tx_hgvs_not_delins)
                                                        except:
                                                            c1 = tx_hgvs_not_delins
                                                        c2 = copy.deepcopy(c1)
                                                        c2.posedit.pos.start = c1.posedit.pos.end
                                                        c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                                                        c2.posedit.pos.end.offset = 0
                                                        c2.posedit.edit.ref = ''
                                                        c2.posedit.edit.alt = ''
                                                        if orientation != -1:
                                                            g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                            g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                            g2.posedit.edit.alt = g2.posedit.edit.ref
                                                        else:
                                                            g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                                                            g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                                                            g1.posedit.edit.alt = g1.posedit.edit.ref
                                                        reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                                                        alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                                                        g3 = copy.deepcopy(g1)
                                                        g3.posedit.pos.end.base = g2.posedit.pos.end.base
                                                        g3.posedit.edit.ref = reference
                                                        g3.posedit.edit.alt = alternate
                                                        c3 = vm.g_to_t(g3, c1.ac)
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
                                            hgvs_alt_genomic = va_func.myvm_t_to_g(hgvs_refreshed_variant, alt_chr,
                                                                                   no_norm_evm, vm, hp, hn, sf, nr_vm)
                                            if hgvs_alt_genomic.posedit.edit.type == 'identity':
                                                re_c = vm.g_to_t(hgvs_alt_genomic, hgvs_refreshed_variant.ac)
                                                if (hn.normalize(re_c)) != (hn.normalize(hgvs_refreshed_variant)):
                                                    shuffle_left_g = copy.copy(hgvs_alt_genomic)
                                                    shuffle_left_g.posedit.edit.ref = ''
                                                    shuffle_left_g.posedit.edit.alt = ''
                                                    shuffle_left_g.posedit.pos.start.base = shuffle_left_g.posedit.pos.start.base - 1
                                                    shuffle_left_g.posedit.pos.end.base = shuffle_left_g.posedit.pos.end.base - 1
                                                    shuffle_left_g = reverse_normalizer.normalize(shuffle_left_g)
                                                    re_c = vm.g_to_t(shuffle_left_g, hgvs_refreshed_variant.ac)
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
                                                lhb = sf.fetch_seq(str(hgvs_alt_genomic.ac), end - 1, end)
                                                rhb = sf.fetch_seq(str(hgvs_alt_genomic.ac), start - 1, start)
                                                hgvs_alt_genomic.posedit.edit.ref = lhb + rhb
                                                hgvs_alt_genomic.posedit.edit.alt = lhb + hgvs_alt_genomic.posedit.edit.alt + rhb
                                                hgvs_alt_genomic.posedit.pos.start.base = end
                                                hgvs_alt_genomic.posedit.pos.end.base = start
                                                hgvs_alt_genomic = hn.normalize(hgvs_alt_genomic)
                                            if hgvs_alt_genomic.posedit.edit.type == 'del':
                                                start = hgvs_alt_genomic.posedit.pos.start.base
                                                end = hgvs_alt_genomic.posedit.pos.end.base
                                                lhb = sf.fetch_seq(str(hgvs_alt_genomic.ac), end - 1, end)
                                                rhb = sf.fetch_seq(str(hgvs_alt_genomic.ac), start - 1, start)
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
                            multi_g.sort()
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
                            for build in genome_builds:
                                test = va_scb.supported_for_mapping(alt_gen_var.ac, build)
                                if test == 'true':
                                    try:
                                        vcf_dict = va_H2V.report_hgvs2vcf(alt_gen_var, build, reverse_normalizer, sf)
                                    except hgvs.exceptions.HGVSInvalidVariantError as e:
                                        continue
                                    # Identify primary assembly positions
                                    if re.match('NC_', alt_gen_var.ac):
                                        if re.match('GRC', build):
                                            primary_genomic_dicts[build.lower()] = {
                                                'hgvs_genomic_description': valstr(alt_gen_var),
                                                'vcf': {'chr': vcf_dict['grc_chr'],
                                                        'pos': vcf_dict['pos'],
                                                        'ref': vcf_dict['ref'],
                                                        'alt': vcf_dict['alt']
                                                        }
                                            }

                                        else:
                                            primary_genomic_dicts[build.lower()] = {
                                                'hgvs_genomic_description': valstr(alt_gen_var),
                                                'vcf': {'chr': vcf_dict['ucsc_chr'],
                                                        'pos': vcf_dict['pos'],
                                                        'ref': vcf_dict['ref'],
                                                        'alt': vcf_dict['alt']
                                                        }
                                            }
                                        if build == 'GRCh38':
                                            vcf_dict = va_H2V.report_hgvs2vcf(alt_gen_var, 'hg38', reverse_normalizer,
                                                                              sf)
                                            primary_genomic_dicts['hg38'] = {
                                                'hgvs_genomic_description': valstr(alt_gen_var),
                                                'vcf': {'chr': vcf_dict['ucsc_chr'],
                                                        'pos': vcf_dict['pos'],
                                                        'ref': vcf_dict['ref'],
                                                        'alt': vcf_dict['alt']
                                                        }
                                            }

                                        continue

                                    else:
                                        if re.match('GRC', build):
                                            dict = {build.lower(): {'hgvs_genomic_description': valstr(alt_gen_var),
                                                                    'vcf': {'chr': vcf_dict['grc_chr'],
                                                                            'pos': vcf_dict['pos'],
                                                                            'ref': vcf_dict['ref'],
                                                                            'alt': vcf_dict['alt']
                                                                            }
                                                                    }
                                                    }
                                        else:
                                            dict = {build.lower(): {'hgvs_genomic_description': valstr(alt_gen_var),
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
                                            vcf_dict = va_H2V.report_hgvs2vcf(alt_gen_var, 'hg38', reverse_normalizer,
                                                                              sf)
                                            dict = {'hg38': {'hgvs_genomic_description': valstr(alt_gen_var),
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
                        if primary_assembly.lower() not in primary_genomic_dicts.keys():
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
                            refseqgene_variant = valstr(hgvs_refseqgene_variant)
                        except:
                            exceptPass()

                    # Add single letter AA code to protein descriptions
                    predicted_protein_variant_dict = {"tlr": str(predicted_protein_variant), "slr": ''}
                    if predicted_protein_variant != '':
                        if not 'Non-coding :n.' in predicted_protein_variant:
                            try:
                                format_p = predicted_protein_variant
                                format_p = re.sub('\(LRG_.+?\)', '', format_p)
                                re_parse_protein = hp.parse_hgvs_variant(format_p)
                                re_parse_protein_singleAA = output_formatter.single_letter_protein(re_parse_protein)
                                predicted_protein_variant_dict["slr"] = str(re_parse_protein_singleAA)
                            except hgvs.exceptions.HGVSParseError:
                                exceptPass()
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
                    ref_records = external.get_urls(dict_out)
                    if ref_records != {}:
                        dict_out['reference_sequence_records'] = ref_records

                    # Append to a list for return
                    batch_out.append(dict_out)
                else:
                    continue
            else:
                continue

        """
        Structure the output into dictionaries rather than a list with descriptive keys
        and a validation type flag
        """
        logger.trace("Populating output dictionary")
        # Create output dictionary
        validation_output = {'flag': None}

        # For gene outputs, i.e. those that hit transcripts
        # dotter = ''
        if set_output_type_flag == 'gene':
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
                            validation_obsolete_counter = validation_obsolete_counter +1
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
        if set_output_type_flag == 'warning':
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
        if set_output_type_flag == 'intergenic':
            validation_output['flag'] = 'intergenic'
            for valid_v in batch_out:
                validation_intergenic_counter = validation_intergenic_counter + 1
                identification_key = 'Intergenic_Variant_%s' % (str(validation_intergenic_counter))

                # Attempt to liftover between genome builds
                # Note: pyliftover uses the UCSC liftOver tool.
                # https://pypi.org/project/pyliftover/
                genomic_position_info = valid_v['primary_assembly_loci']
                for g_p_key in genomic_position_info.keys():

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
                    lifted_response = lift_over(genomic_position_info[g_p_key]['hgvs_genomic_description'], build_from, build_to, hn, vm, vr, hdp, hp, reverse_normalizer, sf, evm)

                    # Sort the respomse into primary assembly and ALT
                    primary_assembly_loci = {}
                    alt_genomic_loci = []
                    for build_key, accession_dict in lifted_response.iteritems():
                        try:
                            accession_key = accession_dict.keys()[0]
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
            metadata["logs"] = logString
        #metadata["variant"] = batch_variant
        #metadata["assembly"] = selected_assembly
        #metadata["transcripts"] = select_transcripts
        #metadata['seqrepo_directory'] = HGVS_SEQREPO_DIR
        #metadata['uta_url'] = UTA_DB_URL
        #metadata['py_liftover_directory'] = PYLIFTOVER_DIR
        #metadata['variantvalidator_data_url'] = VALIDATOR_DB_URL
        #metadata['entrez_id'] = ENTREZ_ID
        metadata['variantvalidator_version'] = VERSION
        metadata['variantvalidator_hgvs_version'] = hgvs_version
        metadata['uta_schema'] = str(hdp.data_version())
        metadata['seqrepo_db'] = HGVS_SEQREPO_DIR.split('/')[-1]
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
        # raise variantValidatorError('Validation error')
        # Return
        # return
        logger.critical(str(exc_type) + " " + str(exc_value))
        logger.debug(str(er))

# Generates a list of transcript (UTA supported) and transcript names from a gene symbol or RefSeq transcript ID
def gene2transcripts(query):
    input = query
    input = input.upper()
    if re.search('\d+ORF\d+', input):
        input = input.replace('ORF', 'orf')
    # Quick check for blank form
    if input == '':
        caution = {'error': 'Please enter HGNC gene name or transcript identifier (NM_, NR_, or ENST)'}
        return caution
    else:
        hgnc = input
        if re.match('NM_', hgnc) or re.match('NR_', hgnc):  # or re.match('ENST', hgnc):
            try:
                tx_info = hdp.get_tx_identity_info(hgnc)
                hgnc = tx_info[6]
            except hgvs.exceptions.HGVSError as e:
                caution = {'error': str(e)}
                return caution

        # First perform a search against the input gene symbol or the symbol inferred from UTA
        initial = va_func.hgnc_rest(path="/fetch/symbol/" + hgnc)
        # Check for a record
        if str(initial['record']['response']['numFound']) != '0':
            current_sym = hgnc
            previous = initial
        # No record found, is it a previous symbol?
        else:
            # Look up current name
            current = va_func.hgnc_rest(path="/search/prev_symbol/" + hgnc)
            # Look for historic names
            # If historic names = 0
            if str(current['record']['response']['numFound']) == '0':
                current_sym = hgnc
            else:
                current_sym = current['record']['response']['docs'][0]['symbol']
            # Look up previous symbols and gene name
            # Re-set the previous variable
            previous = va_func.hgnc_rest(path="/fetch/symbol/" + current_sym)

        # Extract the relevant data
        try:
            previous_sym = previous['record']['response']['docs'][0]['prev_symbol'][0]
        except:
            previous_sym = current_sym

        # Get gene name
        try:
            gene_name = previous['record']['response']['docs'][0]['name']
        except:
            # error = current_sym + ' is not a valid HGNC gene symbol'
            gene_name = 'Gene symbol %s not found in the HGNC database of human gene names www.genenames.org' % query
            return {'error': gene_name}

        # Look up previous name
        try:
            previous_name = previous['record']['response']['docs'][0]['prev_name'][0]
        except:
            previous_name = gene_name

        # Get transcripts
        tx_for_gene = hdp.get_tx_for_gene(current_sym)
        if len(tx_for_gene) == 0:
            tx_for_gene = hdp.get_tx_for_gene(previous_sym)
        if len(tx_for_gene) == 0:
            tx_for_gene = {'error': 'Unable to retrieve data from the UTA, please contact admin'}
            return tx_for_gene

        # Loop through each transcript and get the relevant transcript description
        genes_and_tx = []
        recovered_dict = {}
        for line in tx_for_gene:
            if re.match('^NM_', line[3]) or re.match('^NR_', line[3]):
                # Transcript ID
                tx = line[3]
                tx_description = va_dbCrl.data.get_transcript_description(tx)
                if tx_description == 'none':
                    va_dbCrl.data.update_transcript_info_record(tx, hdp)
                    tx_description = va_dbCrl.data.get_transcript_description(tx)
                # Check for duplicates
                if tx in recovered_dict.keys():
                    continue
                else:
                    try:
                        # Add to recovered_dict
                        recovered_dict[tx] = ''
                        genes_and_tx.append([tx, tx_description, line[1] + 1, line[2]])
                    except:
                        # Add to recovered_dict
                        recovered_dict[tx] = ''
                        genes_and_tx.append([tx, tx_description, 'not applicable', 'not applicable'])
                    # LRG information
                    lrg_transcript = va_dbCrl.data.get_lrgTranscriptID_from_RefSeqTranscriptID(tx)
                    if lrg_transcript == 'none':
                        pass
                    else:
                        genes_and_tx.append([lrg_transcript, tx_description, line[1] + 1, line[2]])

        cp_genes_and_tx = copy.deepcopy(genes_and_tx)
        genes_and_tx = []
        for tx in cp_genes_and_tx:
            tx_d = {'reference': tx[0],
                    'description': tx[1],
                    'coding_start': tx[2] + 1,
                    'coding_end': tx[3]
                    }
            genes_and_tx.append(tx_d)

        # Return data table
        g2d_data = {'current_symbol': current_sym,
                    'previous_symbol': previous_sym,
                    'current_name': gene_name,
                    'previous_name': previous_name,
                    'transcripts': genes_and_tx
                    }

        return g2d_data


# Fetch reference sequence from a HGVS variant description
def hgvs2ref(query):
    logger.info('Fetching reference sequence for ' + query)
    # Dictionary to store the data
    reference = {'variant': query,
                 'start_position': '',
                 'end_position': '',
                 'warning': '',
                 'sequence': '',
                 'error': ''}
    # Step 1: parse the query. Dictionary the parse error if parsing fails
    try:
        input_hgvs_query = hp.parse_hgvs_variant(query)
    except Exception as e:
        reference['error'] = str(e)
    # Step 2: If the variant is a c., it needs to transferred to n.
    try:
        hgvs_query = vm.c_to_n(input_hgvs_query)
    except:
        hgvs_query = input_hgvs_query

    # For transcript reference sequences
    if hgvs_query.type == 'c' or hgvs_query.type == 'n':
        # Step 4: Check for intronic sequence
        if hgvs_query.posedit.pos.start.offset != 0 and hgvs_query.posedit.pos.end.offset != 0:
            reference['warning'] = 'Intronic sequence variation: Use genomic reference sequence'
        elif hgvs_query.posedit.pos.start.offset != 0 or hgvs_query.posedit.pos.end.offset != 0:
            reference['warning'] = 'Partial intronic sequence variation: Returning exonic and/or UTR sequence only'

            # Step 3: split the variant description into the parts required for seqfetching
            accession = hgvs_query.ac
            start = hgvs_query.posedit.pos.start.base - 1
            end = hgvs_query.posedit.pos.end.base

            # Step 5: try and fetch the sequence using SeqFetcher. Dictionary an error if this fails
            try:
                sequence = sf.fetch_seq(accession, start, end)
            except Exception as e:
                reference['error'] = str(e)
                exc_type, exc_value, last_traceback = sys.exc_info()
                te = traceback.format_exc()
                # tr = ''.join(traceback.format_stack())
                tbk = [str(exc_type), str(exc_value), str(te)]
                er = '\n'.join(tbk)
                logger.info(str(exc_type) + " " + str(exc_value))
                logger.debug(er)
            else:
                reference['start_position'] = str(input_hgvs_query.posedit.pos.start.base)
                reference['end_position'] = str(input_hgvs_query.posedit.pos.end.base)
                reference['sequence'] = sequence
        else:
            # Step 3: split the variant description into the parts required for seqfetching
            accession = hgvs_query.ac
            start = hgvs_query.posedit.pos.start.base - 1
            end = hgvs_query.posedit.pos.end.base

            # Step 5: try and fetch the sequence using SeqFetcher. Dictionary an error if this fails
            try:
                sequence = sf.fetch_seq(accession, start, end)
            except Exception as e:
                reference['error'] = str(e)
                exc_type, exc_value, last_traceback = sys.exc_info()
                te = traceback.format_exc()
                # tr = ''.join(traceback.format_stack())
                tbk = [str(exc_type), str(exc_value), str(te)]
                er = '\n'.join(tbk)
                logger.info(er)
            else:
                reference['start_position'] = str(input_hgvs_query.posedit.pos.start.base)
                reference['end_position'] = str(input_hgvs_query.posedit.pos.end.base)
                reference['sequence'] = sequence

    # Genomic reference sequence
    elif hgvs_query.type == 'g' or hgvs_query.type == 'p':
        # Step 3: split the variant description into the parts required for seqfetching
        accession = hgvs_query.ac
        start = hgvs_query.posedit.pos.start.base - 1
        end = hgvs_query.posedit.pos.end.base

        # Step 5: try and fetch the sequence using SeqFetcher. Dictionary an error if this fails
        try:
            sequence = sf.fetch_seq(accession, start, end)
        except Exception as e:
            reference['error'] = str(e)
            exc_type, exc_value, last_traceback = sys.exc_info()
            te = traceback.format_exc()
            # tr = ''.join(traceback.format_stack())
            tbk = [str(exc_type), str(exc_value), str(te)]
            er = '\n'.join(tbk)
            logger.info(str(exc_type) + " " + str(exc_value))
            logger.debug(er)
        else:
            reference['start_position'] = str(input_hgvs_query.posedit.pos.start.base)
            reference['end_position'] = str(input_hgvs_query.posedit.pos.end.base)
            reference['sequence'] = sequence

    # Return the resulting reference sequence or error message
    return reference


def update_vv_data():
    import sys
    #logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    # import update modules
    import mysql_refSeqGene_noMissmatch
    import compile_lrg_data
    # Update refSeqGene Primary assembly alignment data
    mysql_refSeqGene_noMissmatch.update()
    # Update LRG records
    compile_lrg_data.update()

# <LICENSE>
# Copyright (C) 2018  Peter Causey-Freeman, University of Leicester
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
