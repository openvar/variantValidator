import re
import hgvs
from .vvLogging import logger
from . import vvHGVS
from .variant import Variant
from . import vvChromosomes
from . import vvFunctions as fn
from . import gapped_mapping


def collect_transcript_info(variant, validator):
    """Collect transcript information for the variant"""

    if variant.reftype == ':g.':
        toskip = from_genomic(variant, validator)
    else:
        toskip = from_non_genomic(variant, validator)

    return toskip


def get_transcript_info(variant, validator):
    """Collect transcript information from a non-genomic variant"""

    hgvs_vt = validator.hp.parse_hgvs_variant(str(variant.hgvs_formatted))
    try:
        tx_id_info = validator.hdp.get_tx_identity_info(str(hgvs_vt.ac))
    except hgvs.exceptions.HGVSError as e:
        error = 'Please inform UTA admin of the following error: ' + str(e)
        reason = "VariantValidator cannot recover information for transcript " + str(
            hgvs_vt.ac) + ' because it is not available in the Universal Transcript Archive'
        variant.warnings += ': ' + str(reason)
        logger.warning(str(reason) + ": " + str(error))
        return True
    else:
        # Get hgnc Gene name from command
        hgnc = tx_id_info[6]

    # ACCESS THE GENE INFORMATION RECORDS ON THE UTA DATABASE
    # Refseq accession
    tx_for_gene = validator.tx_for_gene(hgnc, validator.hdp)
    refseq_ac = validator.ng_extract(tx_for_gene)

    # Get accurate transcript descriptions from the relevant databases
    # RefSeq databases
    if validator.alt_aln_method != 'genebuild':
        # Gene description  - requires GenBank search to get all the required info, i.e. transcript variant ID
        # accession number
        hgvs_object = validator.hp.parse_hgvs_variant(str(variant.hgvs_formatted))
        accession = hgvs_object.ac
        # Look for the accession in our database
        # Connect to database and send request
        entry = validator.db.in_entries(accession, 'transcript_info')

        # Analyse the returned data and take the necessary actions
        # If the error key exists
        if 'error' in entry:
            # Open a hgvs exception log file in append mode
            error = entry['description']
            variant.warnings += ': ' + str(
                error) + ': A Database error occurred, please contact admin'
            logger.warning(str(error) + ": A Database error occurred, please contact admin")
            return True

        # If the accession key is found
        elif 'accession' in entry:
            # If the current entry is too old
            if entry['expiry'] == 'true':
                try:
                    entry = validator.db.data_add(accession=accession, validator=validator)
                except hgvs.exceptions.HGVSError:
                    error = 'Transcript %s is not currently supported' % (accession)
                    variant.warnings += ': ' + error
                    logger.warning(error)
                    return True
                except Exception:
                    error = 'Unable to assign transcript identity records to ' + accession + ', potentially an obsolete record :'
                    variant.warnings += ': ' + error
                    logger.warning(error)
                    return True
                variant.description = entry['description']
            else:
                variant.description = entry['description']
        # If the none key is found add the description to the database
        elif 'none' in entry:
            try:
                entry = validator.db.data_add(accession=accession, validator=validator)
            except Exception as e:
                logger.warning(str(e))
                error = 'Unable to assign transcript identity records to ' + accession + ', potentially an obsolete record :'
                variant.warnings += ': ' + error
                logger.warning(error)
                return True
            variant.description = entry['description']

        # If no correct keys are found
        else:
            # Open a hgvs exception log file in append mode
            error = 'Unknown error type'
            variant.warnings += ': ' + error + ': A Database error occurred, please contact admin'
            logger.warning(error)
            return True

    # Ensembl databases
    else:
        # accession number
        hgvs_object = validator.hp.parse_hgvs_variant(str(variant.hgvs_formatted))
        accession = hgvs_object.ac
        # Look for the accession in our database
        # Connect to database and send request
        entry = validator.db.in_entries(accession, 'transcript_info')

        # Analyse the returned data and take the necessary actions
        # If the error key exists
        if 'error' in entry:
            # Open a hgvs exception log file in append mode
            error = entry['description']
            variant.warnings += ': ' + str(
                error) + ': A Database error occurred, please contact admin'
            logger.warning(str(error))
            return True

        # If the accession key is found
        elif 'accession' in entry:
            # If the current entry is too old
            if entry['expiry'] == 'true':
                entry = validator.db.data_add(accession=accession, validator=validator)
                variant.description = entry['description']
            else:
                variant.description = entry['description']
        # If the none key is found add the description to the database
        elif 'none' in entry:
            try:
                entry = validator.db.data_add(accession=accession, validator=validator)
            except Exception as e:
                logger.warning(str(e))
                error = 'Unable to assign transcript identity records to ' + accession + ', potentially an obsolete record :'
                variant.warnings += ': ' + error
                logger.warning(error)
                return True
            variant.description = entry['description']

        # If no correct keys are found
        else:
            # Open a hgvs exception log file in append mode
            error = 'Unknown error type'
            variant.warnings += ': ' + error + ': A Database error occurred, please contact admin'
            logger.warning(error)
            return True
    return False


