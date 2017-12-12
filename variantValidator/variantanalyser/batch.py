# -*- coding: utf-8 -*-

"""
batch.py

Contains the link code required to update the transcript_info table when VariantValidator
identifies an out-of-date entry

"""

# Import validator functions
import dbControls.data

# function for adding information to database
def data_add(input, alt_aln_method, accession, dbaction, hp, evm, hdp):
	# Add accurate transcript descriptions to the database
	# RefSeq databases		
	# Get the Entrez (GenBank) file
	dbControls.data.update_transcript_info_record(accession, hdp)
	entry = dbControls.data.in_entries(accession.split('.')[0], 'transcript_info')
	return entry
			
# <LICENSE>

# </LICENSE>