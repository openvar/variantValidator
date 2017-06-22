# -*- coding: utf-8 -*-

# Module containing functions that use hgvs to return batch variant data

# IMPORT REQUIRED PYTHON MODULES
import re
import os

# Set up paths
BATCH_ROOT = os.path.dirname(os.path.abspath(__file__))

# IMPORT HGVS MODULES
import hgvs.exceptions

# Import Biopython modules
from Bio.Seq import Seq
								
# Enables Ensembl Rest variables
import requests, sys
# HGNC rest variables
import httplib2 as http
import json
try:
 	from urlparse import urlparse
except ImportError:
 	from urllib.parse import urlparse

# function for adding information to database
def data_add(input, alt_aln_method, accession, dbaction, hp, evm, hdp):
	
	# Import validator functions
	import dbControls.data
	import functions
	
	# Add accurate transcript descriptions to the database
	# RefSeq databases
	if alt_aln_method != 'genebuild':		
		# Get the Entrez (GenBank) file
		dbControls.data.update_transcript_info_record(accession, hdp)
		entry = dbControls.data.in_entries(accession.split('.')[0], 'transcript_id')
		return entry
			
	# Ensembl databases
	else:
		#return render_template('bootstrap/variantError.html', title=title, user=input, error=accession)
		# THIS STAGE IS CONTROLLED BY CASCADING IF STATEMENTS
		# Use EnsemblRest to search the Ensembl database
		decoded = functions.ensembl_rest(ext = "/overlap/id/" + accession + "?feature=transcript", primary_assembly=primary_assembly)
		if decoded['error'] != 'false':
			hgnc_gene_info = 'Cannot currently display gene information: ' + decoded['error']
		else:
			# Locate the correct transcript from the records
			for features in decoded['record']:
				if str(features['id']) == accession:
					# Extract the transcript version
					tx_variant = features['external_name']
					# Extract the ENSG it originated from
					ensg = features['Parent'] 
					ensembl_gene = ensg
			decoded = ''
			# Search the HGNC (HGNCrest) database for the correct hgnc gene symbol (Ensemble is out of date!!!)
			data = functions.hgnc_rest(path = "/search/ensembl_gene_id/" + ensg)
			if data['error'] != 'false':
				hgnc_gene_info = 'Cannot currently display gene information: ' + decoded['error']
				data = ''
			else:
				# Set the hgnc name correctly
				hgnc = data['record']['response']['docs'][0]['symbol']
				# Now re-search the HGNC database with the correct hgnc name to get the correct description
				data = functions.hgnc_rest(path = "/fetch/symbol/" + hgnc)
				if data['error'] != 'false':
					hgnc_gene_info = 'Cannot currently display gene information: ' + decoded['error']
					data = ''
				else:
					description = data['record']['response']['docs'][0]['name']
					# Check for multiple transcripts to edit the description
					error = 'false'
					var_in = hp.parse_hgvs_variant(input)
					try:
						var_g = evm.c_to_g(var_in)
					except hgvs.exceptions.HGVSError as e:
						error = e
					if error == 'false':
						# return render_template('bootstrap/variantError.html', title=error, user='Data entry failed', error=input, reason=var_g)
						multiple = functions.relevant_transcripts(var_g, evm)
						if len(multiple) > 1:
							# Put together a refseq style description
							desc = 'Homo sapiens ' + description + " (" + hgnc + "), transcript variant " + tx_variant + ", mRNA." 
						else:
							desc = 'Homo sapiens ' + description + " (" + hgnc + "), mRNA." 
						# Set the description to hgnc_gene_info
						data = ''
						data_added = 'false'
						data_added = data.add_entry(accession, desc, 'transcript_id')
						if data_added == 'true':
							entry = data.in_entries(accession, 'transcript_id')
							return entry

