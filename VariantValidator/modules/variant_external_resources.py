"""
variant_external_resources.py contains functions to take the validated
hgvs genomic variant from VariantValidator https://variantvalidator.org/
and to link to external resources using the correct identifier for the
resource’s database
"""

# import modules
import requests
from requests.exceptions import ConnectionError
import sys
import time
import logging
from logging import handlers

"""
Error logging
"""
# set log to the same name as the module
log = logging.getLogger('variant_external_resources')

# screen logging
log.setLevel(logging.INFO)

# Log with a rotating file-handler. This sets the maximum size of the log
# to 0.5Mb and allows two additional logs.
# The logs are then deleted and replaced in rotation
logHandler = handlers.RotatingFileHandler('variant_external_resources.log', maxBytes=500000, backupCount=2)

# We want to minimise the amount of information we log to capturing bugs
logHandler.setLevel(logging.ERROR)
log.addHandler(logHandler)


# function to log exceptions
def log_exception(msg, err):
    """
    Function to log error message and error
    :param msg:
    :param err:
    :return:
    """

    # Create the message and log
    message = 'error occurred at %s with error message: %s and err %s' % \
              (time.ctime(), msg, err)
    log.exception(message, exc_info=True)


"""
Register custom exceptions
"""


class CustomException(Exception):
    """Raised when NCBI has returned an exception with a 200 success code"""


def get_external_resource_links(validated_hgvs_variant):
    """
    function takes incoming validated HGVS variant from VariantValidator
    https://variantvalidator.org/, constructs the parameters for calling
    the NCBI Variation Services APIs and makes 3 calls to function
    get_info_from_variation_services with the parameters:

    api_service_type: NCBI Variant Services groups its APIs into types e.g. SPDI, HGVS, RefSNP
    api_service_name: the name of the API being called

    1st call:
    api service type is HGVS
    api service name is contextuals
    retrieve contextual alleles in SPDI syntax

    2nd call:
    api service type is SPDI
    api service name is rsids
    gets the rsid* for the SPDI

    *Note: rsid is the identifier used to access dbSNP

    3rd call:
    api service type is RefSNP
    api service name is None
    gets the accession number for ClinVar* which is stored in refsnp_snapshot/
    allele_annotations/clinical/accession_version

    *Note: The ClinVar RCV record that links this particular allele to this
     particular phenotype. For more information,
     see http://www.ncbi.nlm.nih.gov/clinvar/intro/

    :param validated_hgvs_variant:
    :return: url_dict: urls for dbSNP and ClinVar where the identifiers
            have been retrieved from NCBI Variation Services
    """

    # check empty string hasn't been passed in
    try:
        if validated_hgvs_variant == '':
            raise ValueError
    except ValueError as err:
        msg = 'Please enter a valid HGVS Genomic Variant'
        # log message and exit gracefully with error code = 1
        log_exception(msg, err)
        raise SystemExit(1)

    # set rsid and accession_version to None.
    # If found, they will be set to the retrieved value
    rsid = None
    accession_version = None

    # 1st call
    api_service_type = 'hgvs'

    # retrieves the contextual alleles equivalent to the HGVS notation input
    api_service_name = '/contextuals?'

    # 1st call to get_info_from_variation_services function
    resp_dict = get_info_from_variation_services(validated_hgvs_variant,
                                                 api_service_type, api_service_name)
    try:

        # NCBI Variation services returns errors in a dictionary
        if 'error' in resp_dict:
            raise CustomException(Exception)
        elif 'data' in resp_dict:

            # construct the SPDI string from the response - we have a dictionary
            # containing a sub-dictionary containing a list
            spdi_id = resp_dict['data']['spdis'][0]['seq_id'] + ':' + \
                      str(resp_dict['data']['spdis'][0]['position']) + ':' + \
                      resp_dict['data']['spdis'][0]['deleted_sequence'] + ':' + \
                      resp_dict['data']['spdis'][0]['inserted_sequence']

    except CustomException:
        msg = "An occurred: " + resp_dict['error']['message']
        err = str(resp_dict['error']['code'])

        # log message and exit gracefully with error code = 1
        log_exception(msg, err)
        raise SystemExit(1)

    # 2nd call to SPDI apis
    api_service_type = 'spdi'

    # Lookup the RSIDs (if any) associated with the input allele
    api_service_name = '/rsids'

    # 2nd call to get_info_from_variation_services function
    resp_dict2 = get_info_from_variation_services(spdi_id, api_service_type, api_service_name)
    try:

        # NCBI Variation services returns errors in a dictionary

        if 'error' in resp_dict2:
            # if error code is 404 not found, don't raise an exception, continue
            if not str(resp_dict2['error']['code']) == '404':
                raise CustomException(Exception)

        elif 'data' in resp_dict2:
            # construct the rsid string - we have a dictionary containing a
            # sub-dictionary containing a list
            rsid = resp_dict2['data']['rsids'][0]

    except CustomException:
        msg = "An error occurred with error code: " + resp_dict2['error']['message']
        err = str(resp_dict2['error']['code'])
        # log message and exit gracefully with error code = 1
        log_exception(msg, err)
        raise SystemExit(1)

    # 3rd call to RefSNP api
    api_service_type = 'refsnp'

    # 3rd call to get_info_from_variation_services function
    # api service_name isn't required for this call so is passed as None
    # only make the call if an rsid was retrieved earlier
    if rsid != None:
        resp_dict3 = get_info_from_variation_services(rsid, api_service_type, None)
        
        try:

            # NCBI Variation services returns errors in a dictionary
            if 'error' in resp_dict3:
                # if error code is 404 not found, we don't want to raise as an exception, want to continue
                if not str(resp_dict3['error']['code']) == '404':
                    raise CustomException

            elif 'primary_snapshot_data' in resp_dict3:

                # primary snap shot data is a dictionary and allele annotations is a list
                allele_annotations = resp_dict3['primary_snapshot_data']['allele_annotations']

                # loop through the allele_annotations list and find clinical entries
                for d in allele_annotations:

                    # check where there is a clinical entry then retrieve the accession number
                    if d['clinical']:
                        accession_version = d['clinical'][0]['accession_version']

        except CustomException:
            msg = "An occurred with error code: " + resp_dict3['error']['message']
            err = str(resp_dict3['error']['code'])
            # log message and exit gracefully with error code = 1
            log_exception(msg, err)
            raise SystemExit(1)

    # construct the URLs to the external resources and add them to a dictionary
    url_dict = {}

    # only add to the dictionary if we have the identifier for the
    # external resource
    if rsid:
        url_dict['dbsnp'] = 'https://www.ncbi.nlm.nih.gov/snp/' + str(rsid)
    if accession_version:
        url_dict['clinvar'] = 'http://www.ncbi.nlm.nih.gov/clinvar/' + accession_version

    # return dictionary containing the 2 URLS
    return url_dict


def get_info_from_variation_services(service_identifier, service_type, service_name):
    """
    called by get_external_resource_links. This function makes the calls
    to NCBI Variation Services

    :param service_identifier: identifier passed into NCBI variation
        services API
    :param service_type: NCBI Variant Services groups its APIs into types
        e.g. SPDI, HGVS, RefSNP
    :param service_name: the name of the API being called
    :return: JSON dictionary with the response from the API call
    """

    # build URL for NCBI Variation Service API - the base URL is always the same
    spdi_base_url = 'https://api.ncbi.nlm.nih.gov/variation/v0/'

    # if service_name is none, it's because we don't need to specify a service
    if service_name == None:
        spdi_url = spdi_base_url + service_type + '/' + str(service_identifier)
    else:
        spdi_url = spdi_base_url + service_type + '/' + service_identifier + service_name
    try:

        # call HCBI Variant Services
        resp = requests.get(spdi_url)
    except ConnectionError:
        msg = 'https://api.ncbi.nlm.nih.gov/variation/v0/ is currently unavailable'

        # add to the log and exit gracefully with error code 1
        log.exception(msg, exc_info=True)
        raise SystemExit(1)

    # return the response from NCBI Variant Services as a JSON dictionary
    resp_dict = resp.json()

    return resp_dict


# remove once integrated into VariantValidator
if __name__ == '__main__':
    get_external_resource_links(sys.argv[1])

# <LICENSE>
# Copyright (C) 2016-2023 VariantValidator Contributors
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
