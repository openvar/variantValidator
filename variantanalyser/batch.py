# -*- coding: utf-8 -*-

# Module containing functions that use hgvs to return batch variant data

# IMPORT REQUIRED PYTHON MODULES
import re
import os
import sys
import requests

# Set up paths
# BATCH_ROOT = os.path.dirname(os.path.abspath(__file__))

# IMPORT HGVS MODULES
import hgvs.exceptions

# Import Biopython modules
from Bio.Seq import Seq
								
# HGNC rest variables
import httplib2 as http
import json
try:
 	from urlparse import urlparse
except ImportError:
 	from urllib.parse import urlparse

# Import validator functions
import dbControls.data
# import functions

# function for adding information to database
def data_add(input, alt_aln_method, accession, dbaction, hp, evm, hdp):
	# Add accurate transcript descriptions to the database
	# RefSeq databases		
	# Get the Entrez (GenBank) file
	dbControls.data.update_transcript_info_record(accession, hdp)
	entry = dbControls.data.in_entries(accession.split('.')[0], 'transcript_info')
	return entry
			