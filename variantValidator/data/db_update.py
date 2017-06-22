# This application uses the flask framework which incorporates:
	# jinja2 template engine
	# werkzeug url engine
	
# IMPORT PYTHON MODULES
import re
import time
import textwrap
import sqlite3
from contextlib import closing
import copy
from decimal import *

# IMPORT FLASK MODULES
from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, make_response
# from flaskext.mail import Mail
from flask_mail import Mail
from flask_mail import Message

# IMPORT HGVS MODULES
import hgvs
import hgvs.parser
import hgvs.dataproviders.uta		
import hgvs.variantmapper
import hgvs.validator
import hgvs.exceptions
import hgvs.variant
import hgvs.location
import hgvs.posedit
import hgvs.edit
import hgvs.variant

# IMPORT variantanalyser functions
import variantanalyser.functions
import variantanalyser.links
import variantanalyser.data

# Import Biopython modules
from Bio.Seq import Seq


# PRE COMPILE HGVS VARIABLES
# Connect to the UTA database and create easy variant mapper based on alignment method
hdp = hgvs.dataproviders.uta.connect()
# From the hgvs parser import, create an instance of hgvs.parser.Parser
hp = hgvs.parser.Parser() 			# From the hgvs parser import, create an instance of hgvs.parser.Parser	
# Validator
vr = hgvs.validator.Validator(hdp=hdp)
# Variant mapper
vm = hgvs.variantmapper.VariantMapper(hdp)

# database configuration
DATABASE = '/local/hgvsWeb_db/descriptions.db'
DEBUG = True
SECRET_KEY = '24.11.79.29.11.80.21.08.14'
USERNAME = 'admin'
PASSWORD = 'tr33house'

# keys for localhost.
RECAPTCHA_PUBLIC_KEY = '6Lc0QRETAAAAAPzrFv_n21bleKADYb9VdwnP5tt_'
RECAPTCHA_PRIVATE_KEY = '6Lc0QRETAAAAALO_veijAZ4-sDUGre2hTgm9gLfu'

# CREATE APP 
app = Flask(__name__)
# configure
app.config.from_object(__name__)

# Create Mail instance
mail = Mail(app)


# APP FUNCTIONS
###############

# connect to database
def connect_db():
    return sqlite3.connect(app.config['DATABASE'])
    
# Initialise the database
def init_db():
    with closing(connect_db()) as db:
        with app.open_resource('schema.sql', mode='r') as f:
            db.cursor().executescript(f.read())
        db.commit()


# Request database connections
##############################
@app.before_request
def before_request():
    g.db = connect_db()

@app.teardown_request
def teardown_request(exception):
    db = getattr(g, 'db', None)
    if db is not None:
        db.close()

# Database functions
####################

def in_entries(connected, accession):
    row = []
    cur = connected.execute('select accession, description from entries where accession = ?',
    	[accession])
    row = cur.fetchone()
    if row is not None:
		entry['accession'] = row[0]
		entry['description'] = row[1]
		return entry
    else:
    	return entry
    	
    	
def add_entry(acn, desc, connected, tm):
    #if not session.get('logged_in'):
    #    abort(401)
    connected.execute('insert into entries (accession, description, updated) values (?, ?, ?)',
    								[acn, desc, tm])
    connected.commit()
    # flash('Updating our database, please allow additional time')
    return 'true'

# Other functions
# Function to convert tri coded amino acids to unicoded
def mon_to_num(month):
	datecode = {
   'Jan' : 1, 'Feb' : 2, 'Mar' : 3, 'Apr' : 4,
    'May' : 5, 'Jun' : 6, 'Jul' : 7 , 'Aug' : 8,
    'Sep' : 9 , 'Oct' : 10, 'Nov' : 11, 'Dec' : 12
	}
	
	get_value = datecode.get(month)
	
	return get_value
	


# Route for HOME
@app.route('/')
def index():
   	title = "Home"
	# Input Primary genome assembly
	primary_assembly = 'GRCh37'
	return render_template('bootstrap/home.html', title=title, primary_assembly=primary_assembly)



# Route for database fetch
@app.route('/fetch/', methods=['POST', 'GET'])
def analysis():
	
	
	# Get current time
	tm = time.ctime()
	#>>> tm
	#'Fri Nov 20 15:17:54 2015'
	
	# Split the c.time string into a list
	tms = tm.split()
	#>>> tms[0]
	#'Fri'
	
	# Get the year, month and day-number
	yr = tms[4]
	mo = tms[1]
	dy = tms[2]
	
	# Convert the month to number
	mo = mon_to_num(mo)
	
	
	# Connect to database and send request
	entry = variantanalyser.data.in_entries(connected=g.db, accession=accession)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	title = "Validation"
	# hgvs errors
	error = 'false'
	
	# COLLECT DATA FROM FORM AND FORMAT
	# Input Primary genome assembly
	primary_assembly = 'GRCh37'
	
	# Flask function request.form[arguement associated with keyword "variant"]
	input = request.args.get('variant', '')
	alt_aln_method = request.args.get('alignment', '')
	# return render_template('bootstrap/variantError.html', title=title, user='giraffe_poo', error=alt_aln_method, reason=input)
	# INITIAL USER INPUT FORMATTING
	# Removes whitespace from the ends of the string
	# Removes gene name (if given)
	# Identifies variant type
	# Returns a dictionary containing the formated input string and the variant type
	# Accepts c, g, n, r currently
	formatted = variantanalyser.functions.user_input(input)
	
	# HTML Template variables
	vars = []
	hgnc = ''
	refseq_gene = ''
	exon_table = ''
	relevant = ''
	warning = ''
	automap = 'false'
	caution = ''
	vmapped = 'false'
	coords = 'false'
	ensembl_gene = 'false'
	hgnc_gene_info = 'false'
	issue_link = 'false'
	
	# Check the initial validity of the input
	if formatted == 'invalid':
		return render_template('bootstrap/variantError.html', title=title, user=input, error='Input not in HGVS format', reason = 'Invalid Variant description')
	else:
		variant = formatted['variant']
		input = formatted['variant']
		type = formatted['type']
	
	# Switch off inversions: currently unsupported 
	#inversion = re.compile('inv')
	#if inversion.search(variant):
	#	return render_template('bootstrap/variantError.html', title=title, user=input, error='Inversions currently unsupported', reason = 'Unable to process')
	# Switch off conversions: currently unsupported 
	conversion = re.compile('con')
	if conversion.search(variant):
		return render_template('bootstrap/variantError.html', title=title, user=input, error='Gene conversions currently unsupported', reason = 'Unable to process')

	
	# Primary check that hgvs will accept the variant 
	error = 'false'
	try:
		hp.parse_hgvs_variant(variant)
	except hgvs.exceptions.HGVSError as e:
		error = e
	if error == 'false':
		pass
		#return render_template('bootstrap/variantError.html', title=title, user=input, error=input, reason = variant)
	else:
		# Open a hgvs exception log file in append mode
		fo = open('/local/logs/hgvsExceptions.txt', 'a')
		excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
		fo.write(excep)
		fo.close()
		# Log the RefSeq Gene as unavailable
		hgvs_refseq = 'Data Unavailable'
		refseq = 'Data Unavailable'
		hgvs_refseq_ac = 'Data Unavailable'
		return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Variant description')
	
	# Ensure ENS sequences are always genebuild
	ens = re.compile('^ENS')
	if ens.search(variant):
		if alt_aln_method != 'genebuild':
			caution = '"'+alt_aln_method+'"' + " alignments are not available for Ensembl sequences"
			alt_aln_method = 'genebuild'
			automap = 'Automap has set the alignment method to "genebuild"'
			# No additional information required
			vmapped = 'false'
			#error = 'giraffe'
			#return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Variant description')
	# Allows bridging between alignment types - ENS to refseq
	elif alt_aln_method == 'genebuild' and type != ':g.':
		error = 'genebuild alignments available for Ensembl (ENS) sequences only'
		return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid alignment option selected')
	else:
		pass
	
	# Create easy variant mapper (over variant mapper) and splign locked evm
	evm = hgvs.variantmapper.EasyVariantMapper(hdp, primary_assembly, alt_aln_method)
	
	
	# handle :r.
	if type == ':r.':
		variant = variant.replace(':r.', ':n.')
		variant = variant.replace('U', 'T')
		# Note, map to :c. currently. May need to look into non-coding transcripts
		# Convert to hgvs object
		hgvs_n = hp.parse_hgvs_variant(variant)
		# try to map to :c.
		error = 'false'
		try:
			hgvs_c = evm.n_to_c(hgvs_n)
		except hgvs.exceptions.HGVSError as e:
			error = e
		if error == 'false':
			# Re-set input and variant
			# input = str(hgvs_c)
			variant = str(hgvs_c)
			type = ':c.'
		else:
			return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid alignment option selected')
	
	# Convert and handle :n.
	# Note, map to :c. currently. May need to look into non-coding transcripts
	if type == ':n.':
		# Convert to hgvs object
		hgvs_n = hp.parse_hgvs_variant(variant)
		# try to map to :c.
		error = 'false'
		try:
			hgvs_c = evm.n_to_c(hgvs_n)
		except hgvs.exceptions.HGVSError as e:
			error = e
		if error == 'false':
			# Re-set input and variant
			old_in = input
			old_var = variant
			input = str(hgvs_c)
			variant = str(hgvs_c)
			type = ':c.'
			caution = 'The use of n. to describe a variant with respect to a protein-coding transcript does not comply with the HGVS variant nomenclature' 
			automap = old_in + ' automapped to ' + variant
		else:
			# Keep everything the same
			pass
		
	
	# Deal with inputs of NG with :c. type
	# Setup relevant transcript string 
	relevant = ''
	ng = re.compile("^NG_")
	ensg = re.compile("^ENSG")
	if type == ':c.':
		#return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Giraffes on the line')
		# If NG found with c type variant
		if ng.search(variant) or ensg.search(variant):
			# Big problem, need to get an transcript but transcript info not available from coding coordinates
			# Step 1, accession number
			# return render_template('bootstrap/variantError.html', title=title, user=input, error=variant)
			hgvs_variant = hp.parse_hgvs_variant(variant)
			accession = hgvs_variant.ac
			# return render_template('bootstrap/variantError.html', title=title, user=input, error=accession)					
			if ng.search(variant):
			# Get the Entrez (GenBank) file
				record =  variantanalyser.functions.entrez_efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
				# Pass the features to an array
				features = record.features
				record = ''
				# Loop through the features and pull out the transcript IDs and append to list
				transcript_ids = []
				# element = 0
				for attributes in features:
					if "transcript_id" in attributes.qualifiers:
						transcript_ids.append(str(attributes.qualifiers["transcript_id"][0]))
					else:
						pass
				# Get the co-ordinates from the original variant
				# Assign the beginning to a variable
				var_str = re.search(r"(NG_\d+\.\d+:c.)", variant)
				var_str = var_str.group(1)
				# Strip away the beginning
				coords = variant.replace(var_str, '')
			if ensg.search(variant):
				# Use EnsemblRest to search the Ensembl database
				decoded = variantanalyser.functions.ensembl_rest(ext = "/overlap/id/" + accession + "?feature=transcript", primary_assembly=primary_assembly)
				if decoded['error'] != 'false':
					reason = 'Cannot currently display the required information:'
					error =  decoded['error']
					return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason=reason)
				else:
					# Locate the correct transcript from the records
					transcript_ids = []
					for attributes in decoded['record']:
						if 'transcript_id' in attributes:
							transcript_ids.append(str(attributes['transcript_id']))
						else:
							pass
				# Get the co-ordinates from the original variant
				# Assign the beginning to a variable
				var_str = re.search(r"(ENSG\d+:c.)", variant)
				var_str = var_str.group(1)
				# Strip away the beginning
				coords = variant.replace(var_str, '')
			# If there is only one transcript then the coordinates will be correct
			if len(transcript_ids) == 1:
				variant = transcript_ids[0] + ':c.' + coords
				caution = 'The format Gene:c.description is valid but not recommended:'
				automap = 'Gene ' + accession + ' encodes a single transcript variant:\n\n Automap has submitted the appropriate transcript level variant to VariantAnalyser'
				# automapping of variant completed
				decoded = ''
				vmapped = 'true'
			else:	
				if ng.search(variant):
					# Recovers additional transcript information from GenBank
					descriptive = []
					for attributes in features:
						if "transcript_id" in attributes.qualifiers:
							tx_descriptive = (str(attributes.qualifiers["transcript_id"][0])) + '(' + str(attributes.qualifiers["gene"][0]) + ': ' + str(attributes.qualifiers["product"][0]) + ')'
							descriptive.append(tx_descriptive)
						else:
							pass
					caution = 'The format Gene:c.description is valid but not recommended'
					automap = 'Gene ' + accession + ' encodes multiple transcript variants \n\n Automap has recovered all available transcripts from the RefSeqGene database: \n\nOther transcripts may exist in the Universal Transcript Archive  \n\nSelect the transcript you require from the list below: coordinates will be added automatically'
					relevant = 'Select the transcript you require and click Submit to analyse' 
					# Give the available transcripts according to Genbank
					rel_var = descriptive
					# Set the variant coordinates 
					coords = ':c.' + coords
					# The web page will list the transcripts available (according to Genbank) and on click 
					# Will add the coordinates and submit
					features = ''
					return render_template('bootstrap/variantAnalysis.html', title=title, user=input,
					relevant=relevant, variants=rel_var, exon_table=exon_table, vars=vars, hgnc='', caution=caution, automap=automap, coords=coords, ensembl_gene=ensembl_gene)
				if ensg.search(variant):
					features = ''
					# Recovers additional transcript information from Ensembl
					# return render_template('bootstrap/variantError.html', title=title, user=input, error=decoded['record'][0])
					descriptive = []
					for attributes in decoded['record']:
						if "transcript_id" in attributes:
							tx_descriptive = attributes['transcript_id'] + '('
							# Search the HGNC (HGNCrest) database for the correct hgnc gene symbol (Ensemble is out of date!!!)
							data = variantanalyser.functions.hgnc_rest(path = "/search/ensembl_gene_id/" + accession)
							if data['error'] != 'false':
								reason = 'Cannot currently display the required information:'
								error =  decoded['error']
								return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason=reason)
							else:
								# Set the hgnc name correctly
								hgnc = data['record']['response']['docs'][0]['symbol']
								# Assemble the description
								tx_descriptive = tx_descriptive + hgnc + ': ' + attributes['external_name'] + ')'
								descriptive.append(tx_descriptive)
						else:
							pass
					caution = 'Multiple transcript overlap the input Ensembl Gene:'
					automap = 'Automap has recovered all available transcripts from the GenBank database: \n\nOther transcripts may exist in the Universal Transcript Archive  \n\nSelect the transcript you require from the list below: coordinates will be added automatically'
					relevant = 'Select the transcript you require and click Submit to analyse' 
					# Give the available transcripts according to Genbank
					rel_var = descriptive
					# Set the variant coordinates 
					coords = ':c.' + coords
					# The web page will list the transcripts available (according to Genbank) and on click 
					# Will add the coordinates and submit
					decoded = ''
					return render_template('bootstrap/variantAnalysis.html', title=title, user=input,
					relevant=relevant, variants=rel_var, exon_table=exon_table, vars=vars, hgnc='', caution=caution, automap=automap, coords=coords, ensembl_gene=ensembl_gene)
	
	# COLLECT gene symbol, name and ACCESSION INFORMATION
	# Gene symbol
	if (type != ':g.'):
		error = 'false'
		hgvs_vt= hp.parse_hgvs_variant(variant)
		try:
			tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
		except hgvs.exceptions.HGVSError as e:
			error = e
		if error != 'false':
			# Open a hgvs exception log file in append mode
			fo = open('/local/logs/hgvsExceptions.txt', 'a')
			excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
			fo.write(excep)
			fo.close()
			error = 'Please inform UTA admin of the following error: ' + str(error)
			issue_link = "https://bitbucket.org/biocommons/uta/issues?status=new&status=open"
			reason = "VariantValidator can not recover information for transcript " + str(hgvs_vt.ac) + ' beacuse it is not available in the Universal Transcript Archive'
			return render_template('bootstrap/variantErrorLink.html', title=title, user=input, error=error, reason = reason, issue_link=issue_link)
		else:
			# Get hgnc Gene name from command
			hgnc = tx_id_info[6]
		
		# ACCESS THE GENE INFORMATION RECORDS ON THE UTA DATABASE 
		# Refseq accession
		tx_for_gene = variantanalyser.functions.tx_for_gene(hgnc, hdp)
		refseq_ac = variantanalyser.functions.ng_extract(tx_for_gene)
		# return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = tx_for_gene)
		
		# Additional gene info
		gene_info = hdp.get_gene_info(hgnc)
		# Chromosomal location
		maploc = gene_info[1]
		chr_loc = ("Chromosome location: " + maploc)
		
		# Get accurate transcript descriptions from the relevant databases
		# RefSeq databases
		if alt_aln_method != 'genebuild':		
			# Gene description  - requires GenBank search to get all the required info, i.e. transcript variant ID
			# accession number
			hgvs_object = hp.parse_hgvs_variant(variant)	
			accession = hgvs_object.ac
			#return render_template('bootstrap/variantError.html', title=title, user=input, error=accession)
			# Look for the accession in our database

			description = entry['description']
			# return render_template('bootstrap/variantError.html', title=accession, user=description, error=accession)
			
			# If the description is present the database search found the accession so continue
			if description != 'false':
				hgnc_gene_info = description
				# return render_template('bootstrap/variantError.html', title='giraffe', user=accession, error=description)
			else:
				return redirect(url_for('add_loading', input=input, alt_aln_method=alt_aln_method, accession=accession))
				
		# Ensembl databases
		else:
			# accession number
			hgvs_coding = variantanalyser.functions.coding(variant, hp)	
			accession = hgvs_coding.ac
			#return render_template('bootstrap/variantError.html', title=title, user=input, error=accession)
			# Look for the accession in our database
			entry = variantanalyser.data.in_entries(connected=g.db, accession=accession)
			description = entry['description']
			# return render_template('bootstrap/variantError.html', title=accession, user=description, error=accession)	
			# If the description is present the database search found the accession so continue
			if description != 'false':
				hgnc_gene_info = description
			else:
				return redirect(url_for('add_loading', input=input, alt_aln_method=alt_aln_method, accession=accession))
		
	# MAP ALL VARIANTS TO CODING VARIANTS
	if type == ':r.':
		try:
			c_from_r = variantanalyser.functions.r_to_c(variant, evm, hp)
		except hgvs.exceptions.HGVSError as e:
			error = e
		if error == 'false':
			c_from_r = variantanalyser.functions.r_to_c(variant, evm, hp)
			variant = c_from_r['variant']
			type = c_from_r['type']
		else:
			# Open a hgvs exception log file in append mode
			fo = open('/local/logs/hgvsExceptions.txt', 'a')
			excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
			fo.write(excep)
			fo.close()
			return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Variant description')
		
	# Genomic type variants will need to be mapped to transcripts
	# Requires user input
	if (type == ':g.'):
		hgvs_genomic = variantanalyser.functions.genomic(variant, evm, hp)
		genomic = str(hgvs_genomic)
		
		# Collect relevant transcripts
		if relevant == '':
			caution = 'Multiple transcripts overlap the input genomic interval:'
			automap = "Automap has compiled a list of transcript level variants that overlap the input genomic interval: \n\nPlease select a variant from the list and click Submit to analyse"
			relevant = 'Transcripts that overlap the input genomic interval'
					
		# rel_var is a keyworded list of relevant transcripts with associated coding variants 
		rel_var = variantanalyser.functions.relevant_transcripts(hgvs_genomic, evm)
		
		# list return statements
		if len(rel_var) == 0:
			# Check for NG_
			rsg = re.compile('^NG_')
			if rsg.search(variant):
				error = 'Mapping unavailable for RefSeqGene ' + variant + ' using alignment method = ' + alt_aln_method
				reason = 'Please contact admin'
				return render_template('bootstrap/variantError.html', title=title, user=input, error=reason, reason = error)
		
		elif len(rel_var) == 1:
			# If there is only 1 transcript, autoselect and submit - Change type to c
			variant = rel_var[0]
			caution = 'A single transcript overlaps the input genomic interval:'
			automap = 'Automap has submitted the appropriate transcript level variant to VariantAnalyser'
			# Variant automapping completed
			vmapped = 'true'
			type = ':c.'
			
			# Get the additional required information
			error = 'false'
			hgvs_vt= hp.parse_hgvs_variant(variant)
			try:
				tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
			except hgvs.exceptions.HGVSError as e:
				error = e
			if error != 'false':
				# Open a hgvs exception log file in append mode
				fo = open('/local/logs/hgvsExceptions.txt', 'a')
				excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
				fo.write(excep)
				fo.close()
				error = 'Please inform UTA admin of the following error: ' + str(error)
				issue_link = "https://bitbucket.org/biocommons/uta/issues?status=new&status=open"
				reason = "VariantValidator can not recover information for transcript " + str(hgvs_vt.ac) + ' beacuse it is not available in the Universal Transcript Archive'
				return render_template('bootstrap/variantErrorLink.html', title=title, user=input, error=error, reason = reason, issue_link=issue_link)
			else:
				# Get hgnc Gene name from command
				hgnc = tx_id_info[6]
	
			# Refseq accession
			tx_for_gene = variantanalyser.functions.tx_for_gene(hgnc, hdp)
			refseq_ac = variantanalyser.functions.ng_extract(tx_for_gene)
	
			# gene name
			gene_info = hdp.get_gene_info(hgnc)
			desc = gene_info[2]
			maploc = gene_info[1]
			chr_loc = ("Chromosome location: " + maploc)
	
			hgnc_gene_info = (hgnc + ": " + desc)
		else:
			#return render_template('bootstrap/variantError.html', title=title, user=input, error=str(rel_var), reason = rel_var)
			# Render the required table
			return render_template('bootstrap/variantAnalysis.html', title=title, user=input,
			relevant=relevant, variants=rel_var, exon_table=exon_table, vars=vars, hgnc=hgnc, caution=caution, automap=automap, coords=coords)

	# TYPE = :c.
	
	if type == ':c.':
		# Flag for validation 
		valid = 'false'

		# Auto info collect
		obj = hp.parse_hgvs_variant(variant)
		tx_ac = obj.ac
		to_g = str(variantanalyser.functions.genomic(variant, evm, hp))
		if to_g == 'error':
			if alt_aln_method != 'genebuild':
				error = 'Please try again later. If the problem persists please contact admin'
				reason = "Currently unable to retrieve required data from NCBI"
				# Open a hgvs exception log file in append mode
				fo = open('/local/logs/hgvsExceptions.txt', 'a')
				excep = "%s -- %s -- %s\n" %(time.ctime(), reason, variant)
				fo.write(excep)
				fo.close()
				return render_template('bootstrap/variantError.html', title=title, user=input, reason=reason, error = error)
			else:
				error = 'Please try again later. If the problem persists please contact admin'
				reason = "Currently unable to retrieve required data from NCBI"
				# Open a hgvs exception log file in append mode
				fo = open('/local/logs/hgvsExceptions.txt', 'a')
				excep = "%s -- %s -- %s\n" %(time.ctime(), reason, variant)
				fo.write(excep)
				fo.close()
				return render_template('bootstrap/variantError.html', title=title, user=input, reason=reason, error = error)
		else:
			pass	
		variant = str(variantanalyser.functions.g_to_c(var_g=to_g, tx_ac=tx_ac, hp=hp, evm=evm))
		tx_ac = ''
		# return render_template('bootstrap/variantError.html', title=title, user=input, error=variant)
		
		# Switch off the deletion validator if the deleted bases are not stated
		# This allows validation of valid deletions that the validator will not permit
		deletion = re.compile("del")
		duplication = re.compile("dup")
		inversion = re.compile("inv")
		if deletion.search(variant) or duplication.search(variant) or inversion.search(variant):
			if type == ':c.':
				hgvs_variant = variantanalyser.functions.coding(variant, hp)
			
			if str(hgvs_variant.posedit.edit.type) == 'del':
				valid = 'true'
			if str(hgvs_variant.posedit.edit.type) == 'delins':
				valid = 'true'
			if str(hgvs_variant.posedit.edit.type) == 'dup':
				valid = 'true'
			if str(hgvs_variant.posedit.edit.type) == 'inv':
				valid = 'true'
		
		
		# INTRONIC OFFSETS - Required for Exon table
		# Variable to collect offset to exon boundary
		ex_offset = 0
		plus = re.compile("\d\+\d") 	# finds digit + digit
		minus = re.compile("\d\-\d")	# finds digit - digit		
		
		cck = 'false'
		# Tackle the plus intronic offset
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
		# return render_template('bootstrap/variantError.html', title=title, user=input, error=variant)
		
		# COORDINATE CHECKER
		# hgvs will handle incorrect coordinates so need to automap errors
		# Make sure any input intronic coordinates are correct
		# Get the desired transcript
		pat_r = re.compile(':r.')
		pat_g = re.compile(':g.')
		if cck == 'true':
			dl = re.compile('del')
			if dl.search(variant):
				coding = variantanalyser.functions.coding(variant, hp)
				trans_acc = coding.ac
				# c to Genome coordinates - Map the variant to the genome
				pre_var = variantanalyser.functions.genomic(variant=variant, evm=evm, hp=hp)
				# genome back to C coordinates
				post_var = variantanalyser.functions.g_to_c(var_g=str(pre_var), tx_ac=trans_acc, hp=hp, evm=evm)
				query = hp.parse_hgvs_variant(post_var)
				test = hp.parse_hgvs_variant(input)
				if query.posedit.pos.start.base != test.posedit.pos.start.base or query.posedit.pos.end.base != test.posedit.pos.end.base:
					caution = 'The entered coordinates do not agree with the intron/exon boundaries for the selected transcript:'
					automap = 'Automap has corrected the coordinates to match the intron/exon boundaries for the selected transcript\n\n(Genomic coordinates have not been changed)'
					# automapping of variant completed
					automap = automap + ':\n\n  ' + input + ' auto-mapped to ' + post_var
					relevant = "Select the automapped transcript and click Submit to analyse"
					rel_var = []
					rel_var.append(post_var)
					# Return an error warning via variant analysis
					return render_template('bootstrap/variantAnalysis.html', title=title, user=input,
					relevant=relevant, variants=rel_var, exon_table=exon_table, vars=vars,
					caution=caution, automap=automap, coords=coords, ensembl_gene=ensembl_gene)
			else:
				coding = variantanalyser.functions.coding(variant, hp)
				trans_acc = coding.ac
				# c to Genome coordinates - Map the variant to the genome
				pre_var = variantanalyser.functions.genomic(variant=variant, evm=evm, hp=hp)
				# genome back to C coordinates
				post_var = variantanalyser.functions.g_to_c(var_g=str(pre_var), tx_ac=trans_acc, hp=hp, evm=evm)
				query = hp.parse_hgvs_variant(post_var)
				test = hp.parse_hgvs_variant(input)
				if query.posedit.pos.start.base != test.posedit.pos.start.base or query.posedit.pos.end.base != test.posedit.pos.end.base:
					caution = 'The entered coordinates do not agree with the intron/exon boundaries for the selected transcript:'
					automap = 'Automap has corrected the coordinates to match the intron/exon boundaries for the selected transcript\n\n(Genomic coordinates have not been changed)'
					# automapping of variant completed
					automap = automap + ':\n\n  ' + input + ' auto-mapped to ' + post_var
					relevant = "Select the automapped transcript and click Submit to analyse"
					rel_var = []
					rel_var.append(post_var)
					# Return an error warning via variant analysis
					return render_template('bootstrap/variantAnalysis.html', title=title, user=input,
					relevant=relevant, variants=rel_var, exon_table=exon_table, vars=vars,
					caution=caution, automap=automap, coords=coords, ensembl_gene=ensembl_gene)
		
		elif pat_r.search(input):
				# change to n.
				inp = input.replace(':r.', ':n.')
				inp = inp.replace('U', 'T') 
				hgvs_inp = hp.parse_hgvs_variant(inp)
				inp_ac = hgvs_inp.ac
				gn = evm.n_to_g(hgvs_inp)
				hgvs_otp = evm.g_to_n(gn, inp_ac)
				otp = str(hgvs_otp)
				query = str(hgvs_otp.posedit.pos)
				test = str(hgvs_inp.posedit.pos)
				query = query.replace('T', 'U')
				test = test.replace('T', 'U')
				output = otp.replace(':n.', ':r.')
				if query != test:
					caution = 'The coordinates for variant ' + input + ' require alteration to comply with HGVS variant nomenclature:'
					automap = 'Automap has corrected the coordinates'
					# automapping of variant completed
					automap = automap + ':\n\n  ' + input + ' auto-mapped to ' + output
					relevant = "Select the automapped transcript and click Submit to analyse"
					rel_var = []
					rel_var.append(output)
					# Return an error warning via variant analysis
					return render_template('bootstrap/variantAnalysis.html', title=title, user=input,
					relevant=relevant, variants=rel_var, exon_table=exon_table, vars=vars,
					caution=caution, automap=automap, coords=coords, ensembl_gene=ensembl_gene)
		
		elif pat_g.search(input):
			pass
			
		else:
			# return render_template('bootstrap/variantError.html', title=title, user=input, error=variant, reason = input)
			post_var = variantanalyser.functions.coding(variant, hp)
			trans_acc = post_var.ac
			query = post_var
			test = hp.parse_hgvs_variant(input)
			# return render_template('bootstrap/variantError.html', title=title, user=input, error=test.posedit.pos, reason = query.posedit.pos)
			if query.posedit.pos != test.posedit.pos:
				caution = 'The coordinates for variant ' + input + ' require alteration to comply with HGVS variant nomenclature:'
				automap = 'Automap has corrected the coordinates'
				# automapping of variant completed
				automap = automap + ':\n\n  ' + str(input) + ' auto-mapped to ' + str(post_var)
				relevant = "Select the automapped transcript and click Submit to analyse"
				rel_var = []
				rel_var.append(post_var)
				# Return an error warning via variant analysis
				return render_template('bootstrap/variantAnalysis.html', title=title, user=input,
				relevant=relevant, variants=rel_var, exon_table=exon_table, vars=vars,
				caution=caution, automap=automap, coords=coords, ensembl_gene=ensembl_gene)
		
		
		# VALIDATION
		# Validate non - del/indel variants	
		if valid == 'false':
			# var_n = str(variantanalyser.functions.rna(variant, evm, hp))
			error = variantanalyser.functions.validate(input=to_g, hp=hp, vr=vr)
			if error == 'false':
				valid = 'true'
			else:
				# Open a hgvs exception log file in append mode
				fo = open('/local/logs/hgvsExceptions.txt', 'a')
				excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
				fo.write(excep)
				fo.close()
				return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Variant description')			
		
		if valid == 'true':
			var_tab = 'true'
			cores = "HGVS-compliant variant descriptions"  + warning
			
			# COLLECT VARIANT DESCRIPTIONS
			# Coding sequence
			hgvs_coding = variantanalyser.functions.coding(variant, hp)
			coding = str(hgvs_coding)
			
			# Transcript start
			hgvs_rna = variantanalyser.functions.rna(variant, evm, hp)
			rna = str(hgvs_rna)
			rna = rna.replace(':n.', ':r.')
			rna = rna.replace('T', 'U')

			# Genomic coordinates
			try:
				hgvs_genomic = variantanalyser.functions.genomic(variant, evm, hp)
			except hgvs.exceptions.HGVSError as e:
				error = e
			if error == 'false':
				hgvs_genomic = variantanalyser.functions.genomic(variant, evm, hp)
				genomic = str(hgvs_genomic)
			else:
				if alt_aln_method == 'genebuild':
					error = 'genebuild alignments available for Ensembl (ENS) sequences only'
					return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Alignment Selection')
				if (alt_aln_method == 'splign') or (alt_aln_method == 'blat'):
					error = 'splign and blat alignments unavailable for Ensembl (ENS) sequences'	
					return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Alignment Selection')
			
			# OBTAIN THE RefSeq gene coordinates
			# Attempt 1 = UTA
			ref_g_dict = variantanalyser.functions.refseq(variant, vm, refseq_ac, hp, evm, hdp)
			
			if ref_g_dict['error'] == 'false':
				hgvs_refseq = ref_g_dict['ref_g']	
				refseq = str(hgvs_refseq)
				hgvs_refseq_ac = hgvs_refseq.ac
			else:
				# Log the unavailable data prior to compensation
				#return render_template('bootstrap/variantError.html', title=title, user=input, error=error)
				# Open a hgvs exception log file in append mode
				error = ref_g_dict['error']
				fo = open('/local/logs/hgvsExceptions.txt', 'a')
				excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
				fo.write(excep)
				fo.close()
				
				# RETRIEVE THE REQUIRED INFORMATION FROM GENBANK RATHER THAN THE UTA
				# Get the accession of the input (should be linked to coding)
				nm_acc = hgvs_coding.ac
				# Accession for the genomic region
				gen_acc = hgvs_genomic.ac
				
				if alt_aln_method != 'genebuild':
					# Get the Entrez (GenBank) file
					record =  variantanalyser.functions.entrez_efetch(db="nuccore", id=nm_acc, rettype="gb", retmode="text")
					# Get the gene ID
					features = record.features
					record = ''
					# Loop through the features and pull out the transcript IDs and append to list
					for attributes in features:
						if 'gene' in attributes.qualifiers:
							gene_id  = attributes.qualifiers['gene'][0]
							break
						else:
							pass
					features = ''
					rsg_hit = 'false'
					# Search the file for the RefSeqGene ID
					with open('/local/logs/gene_RefSeqGene') as inf:
						for line in inf:
							line = line.rstrip()
							cells = line.split('\t')
							if cells[2] == gene_id:
								rsg_id = cells[3]
								rsg_hit = 'true'
								inf.close							
								break
							else:
								rgs_hit = 'false'
								pass
					# Close the file
					inf.close()
				else:
					rsg_hit = 'false'
				if rsg_hit == 'false':
					#return render_template('bootstrap/variantError.html', title=title, user=input, error=error)
					# Open a hgvs exception log file in append mode
					fo = open('/local/logs/hgvsExceptions.txt', 'a')
					error = 'data could not be located in the Entrez database'
					fo = open('/local/logs/hgvsExceptions.txt', 'a')
					excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
					fo.write(excep)
					fo.close()
					# Log the RefSeq Gene as unavailable
					hgvs_refseq = 'Data Unavailable'
					refseq = 'Data Unavailable'
					hgvs_refseq_ac = 'Data Unavailable'
				else:
					# Map the NG to NC and get the coordinates
					ngnc = variantanalyser.functions.ng_to_nc(rsg_id, primary_assembly, gen_acc)
 					# Make a fresh hgvs genomic variant as the base for the new refseq variant
 					hgvs_refseq = variantanalyser.functions.genomic(variant, evm, hp)
 					# Get the orientation
 					ori = variantanalyser.functions.tx_exons(tx_ac=nm_acc, alt_ac=gen_acc, alt_aln_method=alt_aln_method, hdp=hdp)
 					ori = int(ori[0]['alt_strand'])
 					# return render_template('bootstrap/variantError.html', title=title, user=input, error=ori, reason = ori)
 					# Alter the start base and end base
 					if ori != -1:
 						#return render_template('bootstrap/variantError.html', title=title, user=input, error=ori, reason = ori)
 						hgvs_refseq.posedit.pos.start.base = hgvs_refseq.posedit.pos.start.base - int(ngnc['rsg_start'])
						hgvs_refseq.posedit.pos.end.base = hgvs_refseq.posedit.pos.end.base - int(ngnc['rsg_start'])
					else:
						# Need to reverse complement any del/ins sequences
						edit = hgvs_refseq.posedit.edit
						del_rc = 'false'
						ins_rc = 'false'
						deletion = 'false'
						insertion = 'false'
						if re.search(r"((del[GATCUgatcu]+))", str(edit)):
							deletion = 'true'
							dlbs = re.search(r"((del[GATCUgatcu]+))", str(edit))
							del_bases = dlbs.group(1)
							del_bases = del_bases.replace('del', '')
							del_rc = 'true'
						if re.search(r"((ins[GATCUgatcu]+))", str(edit)):
							insertion = 'true'
							inbs = re.search(r"((ins[GATCUgatcu]+))", str(edit))
							ins_bases = inbs.group(1)
							ins_bases = ins_bases.replace('ins', '')
							ins_rc = 'true'
						if del_rc == 'true':
							del_revcomp = variantanalyser.functions.revcomp(bases = del_bases)	
						if ins_rc == 'true':
							ins_revcomp = variantanalyser.functions.revcomp(bases = ins_bases)							
						
						if deletion == 'true' and insertion == 'true':
							hgvs_refseq.posedit.edit = 'del' + del_revcomp + 'ins' + ins_revcomp
						elif deletion == 'true' and insertion == 'false':
							hgvs_refseq.posedit.edit = 'del' + del_revcomp
						elif deletion == 'false' and insertion == 'true':
							hgvs_refseq.posedit.edit = 'ins' + ins_revcomp
						else:
							pass
						# return render_template('bootstrap/variantError.html', title=title, user=input, error=del_bases, reason = str(hgvs_refseq))
						start_base = (int(ngnc['rsg_end'])-int(ngnc['rsg_start'])) - (hgvs_refseq.posedit.pos.end.base-int(ngnc['rsg_start']) -1)
						end_base = (int(ngnc['rsg_end'])-int(ngnc['rsg_start'])) - (hgvs_refseq.posedit.pos.start.base-int(ngnc['rsg_start']) -1)
						hgvs_refseq.posedit.pos.start.base = start_base
						hgvs_refseq.posedit.pos.end.base = end_base
					# Change the accession
					hgvs_refseq.ac = rsg_id
					refseq = str(hgvs_refseq)
					hgvs_refseq_ac = hgvs_refseq.ac
			
			# Predicted effect on protein
			# Translation of inversions currently supported by hgvs - So let's tackle it manually
			inversion = re.compile('inv')
			if inversion.search(variant):
				# Collect the deleted sequence
				del_seq = str(hgvs_coding.posedit.edit)
				del_seq = del_seq.replace('inv', '')
				# Make the inverted sequence
				my_seq = Seq(del_seq)
				inv_seq = my_seq.reverse_complement()
				# Collect the associated protein
				ass_prot = hdp.get_pro_ac_for_tx_ac(hgvs_coding.ac)
				# Intronic inversions go down as uncertain
				int_pl = re.compile('\+')
				int_mi = re.compile('\-')
				if int_pl.search(variant) or int_mi.search(variant):
					# Make the variant
					# return render_template('bootstrap/variantError.html', title=title, user=input, error='poo', reason = 'giraffe')
					#edit = hgvs.edit.AARefAlt(ref='None', alt='None', uncertain='True')
					#posedit = hgvs.posedit.PosEdit(pos='None', edit=edit)
					hgvs_protein = hgvs.variant.SequenceVariant(ac=ass_prot, type='p', posedit='(?)')
					protein = str(hgvs_protein) 
					#return render_template('bootstrap/variantError.html', title=title, user=input, error='poo', reason = protein)
				else:
					# Need to obtain the cds_start
					inf = variantanalyser.functions.tx_identity_info(variant, hdp)
					cds_start = inf[3]
					# return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = str(cds_start))

					# Extract the reference coding sequence from the UTA database
					try:
						ref_seq = variantanalyser.functions.sequence_extractor(ac=hgvs_coding.ac, hdp=hdp)
					except hgvs.exceptions.HGVSError as e:
						error = e
					if error != 'false':
						# Open a hgvs exception log file in append mode
						fo = open('/local/logs/hgvsExceptions.txt', 'a')
						excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
						fo.write(excep)
						fo.close()
						return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Variant description')
					else: 
						pass
					# Create the variant coding sequence
					var_seq = variantanalyser.links.coding_inversion(ref_seq, del_seq, inv_seq, interval_start=hgvs_coding.posedit.pos.start.base+cds_start, interval_end=hgvs_coding.posedit.pos.end.base+cds_start)
					#return render_template('bootstrap/variantError.html', title=title, user=del_seq, error=var_seq, reason = ref_seq)
				
					# Translate the sequences
					prot_ref_seq = variantanalyser.links.translate(ref_seq, cds_start)
					if prot_ref_seq == 'error':
						return render_template('bootstrap/variantError.html', title=title, user=input, error='Translation start codon (ATG) not found at CDS start', reason = 'Request Unavailable')
					prot_var_seq = variantanalyser.links.translate(var_seq, cds_start)
					if prot_var_seq == 'error':
						# Does the edit affect the start codon?
						if hgvs_coding.posedit.pos.start.base >= 1 and hgvs_coding.posedit.pos.start.base <= 3:
							hgvs_protein = hgvs.variant.SequenceVariant(ac=ass_prot, type='p', posedit='(Met1?)')
							protein = str(hgvs_protein)
						else:
							return render_template('bootstrap/variantError.html', title=title, user=input, error='Translation start codon (ATG) not found at CDS start', reason = 'Request Unavailable')
					else:				
				
						# Gather the required information regarding variant interval and sequences
						pro_inv_info = variantanalyser.links.pro_inv_info(prot_ref_seq, prot_var_seq)
						#return render_template('bootstrap/variantError.html', title=title, user=del_seq, error=pro_inv_info, reason = pro_inv_info)

					# Apply tests
						#info = {
						#'variant' : 'true',
						#'prot_del_seq' : '',
						#'prot_ins_seq' : '',
						#'edit_start' : 0,
						#'edit_end' : 0,
						#'terminate' : 'false',
						#'ter_pos' : 0,
						#'error' : 'false'
						#}
					
						if pro_inv_info['error'] == 'true':
							return render_template('bootstrap/variantError.html', title=title, user=input, error='Please contact admin', reason = 'Translation error occurred')
						elif pro_inv_info['variant'] != 'true':
							# Make the variant
							hgvs_protein = hgvs.variant.SequenceVariant(ac=ass_prot, type='p', posedit='(=)')
							protein = str(hgvs_protein)
						else:
							if pro_inv_info['terminate'] == 'true':
								end = 'Ter' + str(pro_inv_info['ter_pos'])
								pro_inv_info['prot_ins_seq'].replace('*', end) 
							iv = hgvs.location.Interval(start=pro_inv_info['edit_start'], end=pro_inv_info['edit_end'])
							# Note for hgvs to continue working, we beed to take the format delXXXinsyyy
							# Need to recode the single letter del and ins sequences
							del_thr = variantanalyser.links.one_to_three(pro_inv_info['prot_del_seq'])
							ins_thr = variantanalyser.links.one_to_three(pro_inv_info['prot_ins_seq'])
							# Make the edit
							del_len = len(del_thr)
							from_aa = del_thr[0:3]
							to_aa = del_thr[del_len-3:]
							posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + str(pro_inv_info['edit_end']) + 'delins' + ins_thr + ')'
							# edit = hgvs.edit.AARefAlt(ref=pro_inv_info['prot_del_seq'], alt=pro_inv_info['prot_ins_seq'])
							#posedit = hgvs.posedit.PosEdit(pos=iv,edit=edit)
							no = 'false'
							try:
								hgvs_p = hgvs.variant.SequenceVariant(ac=ass_prot, type='p', posedit=posedit)
							except hgvs.exceptions.HGVSError as e:
								no = e
							if no == 'false':
								prot = str(hgvs_p)
								hgvs_protein = variantanalyser.functions.hgvs_protein(prot, hp)
								protein = str(hgvs_protein)
							else:
								fo = open('/local/logs/hgvsExceptions.txt', 'a')
								excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
								fo.write(excep)
								fo.close()
								return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Inversion validation requires additional testing. Please contact admin stating the input variant description')


						#return render_template('bootstrap/variantError.html', title=title, user=str(hgvs_protein), error=protein, reason = protein)
				
			else: 
				hgvs_protein = variantanalyser.functions.protein(variant, evm, hp)
				protein = str(hgvs_protein)
			# Assess whether the protein edit is known
			unknown = re.compile('\?')
			# Obtain the edit
			# Search the edit for the ?
			if unknown.search(protein):
				puk = 'true'
			else:
				puk = 'false'
			# Is the unknown protein within coding region?
			pls = re.compile("\+")
			mns = re.compile("\-")
			cdn = ''
			if pls.search(variant) and puk == 'true':
				cdn = 'true'
			elif mns.search(variant) and puk == 'true':
				cdn = 'true'
			else:
				cdn='false'
			
			# External links
			genbank = "http://www.ncbi.nlm.nih.gov/nuccore/"
			ensembl = "http://www.ensembl.org/id/"
			

			# PROCESS VARIANT INFORMATION INTO A TABLE
			# List containing title, variant and Accessions for each variant
			if (alt_aln_method == 'splign') or (alt_aln_method == 'blat'):
				vars = [
							['Genome (GRCh37) (:g.)', genomic, hgvs_genomic.ac, 'genbank', 'genome', cdn],
							['Coding sequence (:c.)', coding, hgvs_coding.ac, 'genbank', 'coding', cdn],
							['Transcript start (:r.)', rna, hgvs_rna.ac, 'genbank', 'transcript', cdn],
							['RefSeqGene (:g.)', refseq, hgvs_refseq_ac, 'genbank', 'refseq', cdn],
							['Predicted effect on protein (:p.)',	protein, hgvs_protein.ac, 'genbank', 'protein', puk]
			 			]
			
			if (alt_aln_method == 'genebuild'):
				# Need to obtain the refseq gene to allowfor full scope of alignments
				error = 'false'
				hgvs_vt= hp.parse_hgvs_variant(variant)
				try:
					tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
				except hgvs.exceptions.HGVSError as e:
					error = e
				if error != 'false':
					# Open a hgvs exception log file in append mode
					fo = open('/local/logs/hgvsExceptions.txt', 'a')
					excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
					fo.write(excep)
					fo.close()
					error = 'Please inform UTA admin of the following error: ' + str(error)
					issue_link = "https://bitbucket.org/biocommons/uta/issues?status=new&status=open"
					reason = "VariantValidator can not recover information for transcript " + str(hgvs_vt.ac) + ' beacuse it is not available in the Universal Transcript Archive'
					return render_template('bootstrap/variantErrorLink.html', title=title, user=input, error=error, reason = reason, issue_link=issue_link)
				else:
					# Get hgnc Gene name from command
					cds_start = tx_id_info[3]
					cds_end = tx_id_info[4]
	
				# Refseq accession
				tx_for_gene = variantanalyser.functions.tx_for_gene(hgnc, hdp)
				refseq_ac = variantanalyser.functions.ng_extract(tx_for_gene)
				
				# RefSeq gene coordinates
				ref_g_dict = variantanalyser.functions.refseq(variant, vm, refseq_ac, hp, evm, hdp)
				
				if ref_g_dict['error'] == 'false':
					hgvs_refseq = ref_g_dict['ref_g']	
				else:
					#return render_template('bootstrap/variantError.html', title=title, user=input, error=error)
					# Open a hgvs exception log file in append mode
					error = ref_g_dict['error']
					fo = open('/local/logs/hgvsExceptions.txt', 'a')
					excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
					fo.write(excep)
					fo.close()
					# Log the RefSeq Gene as unavailable
					hgvs_refseq = 'Data Unavailable'
					refseq = 'Data Unavailable'
					hgvs_refseq_ac = 'Data Unavailable'

				
				vars = [
							['Genome (GRCh37) (:g.)', genomic, hgvs_genomic.ac, 'genbank', 'genome', cdn],
							['Coding sequence (:c.)', coding, hgvs_coding.ac, 'ensembl', 'coding', cdn],
							['Transcript start (:r.)', rna, hgvs_rna.ac, 'ensembl', 'transcript', cdn],
							['RefSeqGene (:g.)', refseq, hgvs_refseq_ac, 'genbank', 'refseq', cdn],
							['Predicted effect on protein (:p.)',	protein, hgvs_protein.ac, 'ensembl', 'protein', puk]
			 			]
			
			# return render_template('bootstrap/variantError.html', title=title, user=input, error=vars[2][5], reason = 'Mr Giraffe')
			
			
			# COMPILE THE EXON TABLE
			# Align the exons: coding to genomic
			tx_ac = hgvs_coding.ac
			alt_ac = hgvs_genomic.ac
		
			tx_exons = variantanalyser.functions.tx_exons(tx_ac, alt_ac, alt_aln_method, hdp)
		
			# Generate exon table
			# Collect the interval positions required and integer the values
			tx_int_start = str(hgvs_coding.posedit.pos.start.base + hgvs_coding.posedit.pos.start.offset)
			tx_int_start = int(tx_int_start)
			cd_pos = str(hgvs_rna.posedit.pos.start.base)
			cd_pos = int(cd_pos)
		
			# Handle 0 position and -ve offsets
			
			if tx_int_start <= 0:
				cd_offset = tx_int_start - cd_pos - ex_offset + 1
			else:
				cd_offset = tx_int_start - cd_pos -ex_offset 

			gen_int_start = str(hgvs_genomic.posedit.pos.start.base)
			gen_int_start = int(gen_int_start)
		
			gen_int_end = str(hgvs_genomic.posedit.pos.end.base)
			gen_int_end = int(gen_int_end)
		
			exon_table = variantanalyser.functions.exon_table(tx_exons, gen_int_start, gen_int_end, cd_offset)
			
			# Is the second row of the exon_table required?
			end_row = 'false'
			if exon_table['start_info']['id'] != exon_table['end_info']['id']:
				end_row = 'true' 
			
			# COLLECT OVERLAPPING RELEVANT TRANSCRIPTS
			# Collect relevant transcripts
			pos = str(hgvs_genomic.posedit.pos)	
			relevant = "Transcript level variants overlapping genomic coordinates " + str(hgvs_genomic.ac) + ":" + pos
			# rel_var is a keyworded list of relevant transcripts with associated coding variants 
			rel_var = variantanalyser.functions.relevant_transcripts(hgvs_genomic, evm)

	# TYPE = :n.
	
	if type == ':n.':
		# Flag for validation 
		valid = 'false'

		# Auto info collect
		obj = hp.parse_hgvs_variant(variant)
		tx_ac = obj.ac
		to_g = str(variantanalyser.functions.genomic(variant, evm, hp))
		if to_g == 'error':
			if alt_aln_method != 'genebuild':
				error = 'Please try again later. If the problem persists please contact admin'
				reason = "Currently unable to retrieve required data from NCBI"
				# Open a hgvs exception log file in append mode
				fo = open('/local/logs/hgvsExceptions.txt', 'a')
				excep = "%s -- %s -- %s\n" %(time.ctime(), reason, variant)
				fo.write(excep)
				fo.close()
				return render_template('bootstrap/variantError.html', title=title, user=input, reason=reason, error = error)
			else:
				error = 'Please try again later. If the problem persists please contact admin'
				reason = "Currently unable to retrieve required data from NCBI"
				# Open a hgvs exception log file in append mode
				fo = open('/local/logs/hgvsExceptions.txt', 'a')
				excep = "%s -- %s -- %s\n" %(time.ctime(), reason, variant)
				fo.write(excep)
				fo.close()
				return render_template('bootstrap/variantError.html', title=title, user=input, reason=reason, error = error)
		else:
			pass	
		variant = str(variantanalyser.functions.g_to_n(var_g=to_g, tx_ac=tx_ac, hp=hp, evm=evm))
		tx_ac = ''
		# return render_template('bootstrap/variantError.html', title=title, user=input, error=variant)
		
		# VALIDATION
		# Validate non - del/indel variants	
		if valid == 'false':
			# var_n = str(variantanalyser.functions.rna(variant, evm, hp))
			error = variantanalyser.functions.validate(input=to_g, hp=hp, vr=vr)
			if error == 'false':
				valid = 'true'
			else:
				# Open a hgvs exception log file in append mode
				fo = open('/local/logs/hgvsExceptions.txt', 'a')
				excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
				fo.write(excep)
				fo.close()
				return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Variant description')			
		
		if valid == 'true':
			var_tab = 'true'
			cores =  "HGVS-compliant variant descriptions"  + warning
			
			# COLLECT VARIANT DESCRIPTIONS
			# Coding sequence
			
			# Transcript start
			hgvs_rna = hp.parse_hgvs_variant(variant)
			rna = str(hgvs_rna)

			# Genomic coordinates
			try:
				hgvs_genomic = variantanalyser.functions.genomic(variant, evm, hp)
			except hgvs.exceptions.HGVSError as e:
				error = e
			if error == 'false':
				hgvs_genomic = variantanalyser.functions.genomic(variant, evm, hp)
				genomic = str(hgvs_genomic)
			else:
				if alt_aln_method == 'genebuild':
					error = 'genebuild alignments available for Ensembl (ENS) sequences only'
					return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Alignment Selection')
				if (alt_aln_method == 'splign') or (alt_aln_method == 'blat'):
					error = 'splign and blat alignments unavailable for Ensembl (ENS) sequences'	
					return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Alignment Selection')
			
			# Format what is left into the web page!			
			# External links
			genbank = "http://www.ncbi.nlm.nih.gov/nuccore/"
			ensembl = "http://www.ensembl.org/id/"
			cdn = 'false'

			# PROCESS VARIANT INFORMATION INTO A TABLE
			# List containing title, variant and Accessions for each variant
			if (alt_aln_method == 'splign') or (alt_aln_method == 'blat'):
				vars = [
							['Genome (GRCh37) (:g.)', genomic, hgvs_genomic.ac, 'genbank', 'genome', cdn],
							['Transcript start (:n.)', rna, hgvs_rna.ac, 'genbank', 'transcript', cdn]
			 			]
			
			if (alt_aln_method == 'genebuild'):
				vars = [
							['Genome (GRCh37) (:g.)', genomic, hgvs_genomic.ac, 'genbank', 'genome', cdn],
							['Transcript start (:n.)', rna, hgvs_rna.ac, 'ensembl', 'transcript', cdn]
			 			]

			# COLLECT OVERLAPPING RELEVANT TRANSCRIPTS
			# Collect relevant transcripts	
			pos = str(hgvs_genomic.posedit.pos)	
			relevant = "Transcript level variants overlapping genomic coordinates " + str(hgvs_genomic.ac) + ":" + pos
			# rel_var is a keyworded list of relevant transcripts with associated coding variants 
			rel_var = variantanalyser.functions.relevant_transcripts(hgvs_genomic, evm)
			
			# Set the unrequired variables so that the form fills
			end_row = ''
			
	# Render the required template
	if (alt_aln_method == 'splign') or (alt_aln_method == 'blat') or (alt_aln_method == 'genebuild'):
		if automap != 'false' and vmapped != 'false':
			automap = automap + ':\n\n  ' + input + ' auto-mapped to ' + variant
		else:
			pass 
		return render_template('bootstrap/variantAnalysis.html', title=title, hgnc=hgnc_gene_info,
		refseq_gene=chr_loc, cores=cores, user=input, vars=vars, relevant=relevant,
		variants=rel_var, exon_table=exon_table, alt_aln_method=alt_aln_method, end_row=end_row, 
		caution=caution, genbank=genbank, ensembl=ensembl, automap=automap, coords=coords, ensembl_gene=ensembl_gene)
		

# Route for variantlinks
@app.route('/variantlinks/', methods=['POST', 'GET'])
def links():
	
	title = "Sequences"
	# hgvs errors
	error = 'false'
	
	# Check for requirement once complete
	valid = 'true'
	message = ''
	auto_correct = ''
	hgnc = ''
	refseq_gene = ''
	hgnc_gene_info = ''
	chr_loc = ''
	cds_start = 0
	cds_end = 0
	alt_start_i = 0
	alt_end_i = 0
	ins_len = 0
	ref_seq = ''
	frame = ''
	edit_type = ''
	newmethod = ''
	
	# Set the current offset for the alignments
	plus_minus = 100
	
	# COLLECT DATA FROM FORM AND FORMAT
	# Input Primary genome assembly
	primary_assembly = 'GRCh37'
	
	# Handle form data
	# Example input URL
	# https://www22.lamp.le.ac.uk/hgvs/variantlinks/?variant=NM_000088.3%3Ac.589G%3ET&alignment=splign&alt_acc=NG_007400.1:g.8638G>T&request=reference
	input = request.args.get('variant', '')
	alt_aln_method = request.args.get('alignment', '')
	output = request.args.get('output', '')
	origin = request.args.get('origin', '')
	bridge = request.args.get('bridge', '')
	
	# INITIAL USER INPUT FORMATTING
	# Removes whitespace from the ends of the string
	# Removes gene name (if given)
	# Identifies variant type
	# Returns a dictionary containing the formated input string and the variant type
	# Accepts c, g, n, r currently: Protein added 15.07.15
	formatted = variantanalyser.functions.user_input(input)

	# Check the initial validity of the input
	if formatted == 'invalid':
		return render_template('bootstrap/variantError.html', title=title, user=input, error='Input not in HGVS format', reason = 'Invalid Variant description')
	else:
		variant = formatted['variant']
		type = formatted['type']
		if type == ':r.':
			type = ':n.'
			variant = variant.replace(':r.', ':n.')
			variant = variant.replace('U', 'T')
	
	# Deal with ENS inputs where Genebuild is not selected
	ens = re.compile("^ENS")
	if (ens.search(variant)):
		if (alt_aln_method != 'genebuild'):
			alt_aln_method = 'genebuild'
			message = ' (autoselect \'genebuild\' alignment for Ensembl sequences)'

	# Deal with inputs of NG with :c. type
	# Setup relevant transcript string 
	relevant = ''
	ng = re.compile("^NG_")
	if type == ':c.':
		# If NG found with c type variant
		if ng.search(variant):
			#set type to g
			type = ':g.'
			# replace variant c to g
			variant = variant.replace(":c.", ":g.")
			message = "RefSeqGene accessions are genomic - select a corresponding transcript below"
	
	# Create easy variant mapper (over variant mapper)
	evm = hgvs.variantmapper.EasyVariantMapper(hdp, primary_assembly, alt_aln_method)
	
	# COLLECT gene symbol, name and ACCESSION INFORMATION
	# Gene symbol
	if (type != ':g.') and (type != ':p.'):
	
		tx_id_info = variantanalyser.functions.tx_identity_info(variant, hdp)
		cds_start = tx_id_info[3]
		cds_end = tx_id_info[4]
		
		if tx_id_info != None:
			# Get hgnc Gene name from command
			hgnc = tx_id_info[6]
		
		# Get accurate transcript descriptions from the relevant databases
		# RefSeq databases
		if alt_aln_method != 'genebuild':		
			# Gene description  - requires GenBank search to get all the required info, i.e. transcript variant ID
			# accession number
			if type == ':c.':
				hgvs_coding = variantanalyser.functions.coding(variant, hp)	
				accession = hgvs_coding.ac
			if type == ':n.':
				hgvs_coding = variantanalyser.functions.hgvs_rna(variant, hp)	
				accession = hgvs_coding.ac
			#return render_template('bootstrap/variantError.html', title=title, user=variant, error=type)
			# Get the Entrez (GenBank) file
			#record = variantanalyser.functions.entrez_efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
			# Extract the information required
			#desc = record.description
			#record = ''
			#return render_template('bootstrap/variantError.html', title=title, user=input, error=desc)
			# hgnc_gene_info = (hgnc + ": " + desc)
			entry = variantanalyser.data.in_entries(connected=g.db, accession=accession)
			description = entry['description']
			
			hgnc_gene_info = description	
		# Ensembl databases
		else:
			# accession number
			if type == ':c.':
				hgvs_coding = variantanalyser.functions.coding(variant, hp)	
				accession = hgvs_coding.ac
			if type == ':r.':
				hgvs_coding = variantanalyser.functions.hgvs_rna(variant, hp)	
				accession = hgvs_coding.ac
			#return render_template('bootstrap/variantError.html', title=title, user=input, error=accession)
			# THIS STAGE IS CONTROLLED BY CASCADING IF STATEMENTS
			# Use EnsemblRest to search the Ensembl database
			decoded = variantanalyser.functions.ensembl_rest(ext = "/overlap/id/" + accession + "?feature=transcript", primary_assembly=primary_assembly)
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
				data = variantanalyser.functions.hgnc_rest(path = "/search/ensembl_gene_id/" + ensg)
				if data['error'] != 'false':
					hgnc_gene_info = 'Cannot currently display gene information: ' + decoded['error']
					data = ''
				else:
					# Set the hgnc name correctly
					hgnc = data['record']['response']['docs'][0]['symbol']
					# Now re-search the HGNC database with the correct hgnc name to get the correct description
					data = variantanalyser.functions.hgnc_rest(path = "/fetch/symbol/" + hgnc)
					if data['error'] != 'false':
						hgnc_gene_info = 'Cannot currently display gene information: ' + decoded['error']
						data = ''
					else:
						description = data['record']['response']['docs'][0]['name']
						# Check for multiple transcripts to edit the description
						if type == ':c.':
							try:
								var_g = variantanalyser.functions.genomic(variant, evm, hp)
							except hgvs.exceptions.HGVSError as e:
								error = e
						if type == ':n.':
							c_from_r = variantanalyser.functions.r_to_c(variant, evm, hp)
							c_type = c_from_r['variant']
							try:
								var_g = variantanalyser.functions.genomic(variant=str(c_type), evm=evm, hp=hp)
							except hgvs.exceptions.HGVSError as e:
								error = e
						if error == 'false':
							if type == ':c.':
								var_g = variantanalyser.functions.genomic(variant, evm, hp)
							if type == ':n.':
								var_g = variantanalyser.functions.genomic(variant=str(c_type), evm=evm, hp=hp)
							multiple = variantanalyser.functions.relevant_transcripts(var_g, evm)
							if len(multiple) > 1:
								# Put together a refseq style description
								desc = 'Homo sapiens ' + description + " (" + hgnc + "), transcript variant " + tx_variant + ", mRNA." 
							else:
								desc = 'Homo sapiens ' + description + " (" + hgnc + "), mRNA." 
							# Set the description to hgnc_gene_info
							hgnc_gene_info = desc
							data = ''
						else:
							# Open a hgvs exception log file in append mode
							fo = open('/local/logs/hgvsExceptions.txt', 'a')
							excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
							fo.write(excep)
							fo.close()
							data = ''
							return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Variant description')


		# gene name
		gene_info = hdp.get_gene_info(hgnc)
		maploc = gene_info[1]
		chr_loc = ("Chromosome location: " + maploc)
	
		
	# Coding sequence requests
	##########################
	if type == ':c.':
		# Return hgvs object of the variant and RNA
		if valid == 'true':
			hgvs_query = variantanalyser.functions.coding(variant, hp)
		
			# Map to the RNA variant to keep the co-ordinates true to the Accession
			hgvs_reference = variantanalyser.functions.rna(variant, evm, hp)
			
			# Assign the accession string to on object
			tx_ac = hgvs_query.ac
			
			# Extract the sequence from the UTA database
			try:
				ac_seq = variantanalyser.functions.sequence_extractor(ac=tx_ac, hdp=hdp)
			except hgvs.exceptions.HGVSError as e:
				error = e
			
			if error != 'false':
				# Open a hgvs exception log file in append mode
				fo = open('/local/logs/hgvsExceptions.txt', 'a')
				excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
				fo.write(excep)
				fo.close()
				return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Variant description')
			
			else: 
				#ac_seq = variantanalyser.functions.sequence_extractor(ac=tx_ac, hdp=hdp)
				pass
				
				# Mark coding in uppercase and UTRs in lowercase
				five_utr = ac_seq[0:cds_start]
				five_utr = five_utr.lower()
				cds = ac_seq[cds_start: cds_end]
				three_utr = ac_seq[cds_end: ]
				three_utr = three_utr.lower()
				# Re join the sequence
				ac_seq = five_utr + cds + three_utr
			
			# Obtain the correct variant descriptions
			
	# Transcript sequence requests
	##############################
	if type == ':n.':
		# Return hgvs object of the variant and RNA
		if valid == 'true':
			# Use reference as coordinates are correct
			hgvs_reference = variantanalyser.functions.hgvs_rna(variant, hp)
			
			# Assign the accession string to on object
			tx_ac = hgvs_reference.ac
			
			# Extract the sequence from the UTA database
			try:
				ac_seq = variantanalyser.functions.sequence_extractor(ac=tx_ac, hdp=hdp)
			except hgvs.exceptions.HGVSError as e:
				error = e
			
			if error != 'false':
				# Open a hgvs exception log file in append mode
				fo = open('/local/logs/hgvsExceptions.txt', 'a')
				excep = "%s -- %s\n" %(time.ctime(), error)
				fo.write(excep)
				fo.close()
				return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Variant description')
			
			else: 
				ac_seq = variantanalyser.functions.sequence_extractor(ac=tx_ac, hdp=hdp)
				
				# Mark coding in uppercase and UTRs in lowercase
				five_utr = ac_seq[0:cds_start]
				five_utr = five_utr.lower()
				cds = ac_seq[cds_start: cds_end]
				three_utr = ac_seq[cds_end: ]
				three_utr = three_utr.lower()
				# Re join the sequence
				ac_seq = five_utr + cds + three_utr


	# RefSeqGene requests
	#####################
	if type == ':g.':
		# Return hgvs object of the variant and RNA
		if valid == 'true':
			# Use reference as coordinates are correct
			hgvs_reference = variantanalyser.functions.hgvs_refseq(variant, hp)
			
			# Assign the accession string to on object
			tx_ac = hgvs_reference.ac
			
			# Extract the sequence from the UTA database
			try:
				ac_seq = variantanalyser.functions.sequence_extractor(ac=tx_ac, hdp=hdp)
			except hgvs.exceptions.HGVSError as e:
				error = e
			
			if error != 'false':
				# Open a hgvs exception log file in append mode
				fo = open('/local/logs/hgvsExceptions.txt', 'a')
				excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
				fo.write(excep)
				fo.close()
				return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Variant description')
			
			else: 
				ac_seq = variantanalyser.functions.sequence_extractor(ac=tx_ac, hdp=hdp)
			
			# Annotate the introns/exons wrt origin variant
			if origin != '':
				
				# First, lower case the sequence string
				ac_seq = ac_seq.lower()
				#return render_template('bootstrap/variantError.html', title=title, user=input, error=ac_seq, reason = 'Request Unavailable')
				# hgvs variant the origin variant
				origin = hp.parse_hgvs_variant(origin)
				
				# Assign the origin variant accession string to on object
				alt_ac = origin.ac
				
				# Genomic type variants will need to be mapped to transcripts to extract HGNC info
				# hgvs variant the genomic variant
				hgvs_genomic = variantanalyser.functions.genomic(variant, evm, hp)
		
				tx_id_info = variantanalyser.functions.tx_id_info(alt_ac, hdp)
		
				hgnc = tx_id_info[6]

				# Refseq accession
				tx_for_gene = variantanalyser.functions.tx_for_gene(hgnc, hdp)
				refseq_ac = variantanalyser.functions.ng_extract(tx_for_gene)
	
				# Get accurate descriptions from the relevant databases
				# RefSeq databases	
				# Gene description  - requires GenBank search to get all the required info, i.e. transcript variant ID
				# accession number
				hgvs_refseq = variantanalyser.functions.hgvs_refseq(variant, hp)	
				accession = hgvs_refseq.ac
				# return render_template('bootstrap/variantError.html', title=title, user=input, error=accession)
				# Get the Entrez (GenBank) file
				#record = variantanalyser.functions.entrez_efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
				# Extract the information required
				#desc = record.description
				#record = ''
				#return render_template('bootstrap/variantError.html', title=title, user=input, error=desc)
				# hgnc_gene_info = (hgnc + ": " + desc)
				entry = variantanalyser.data.in_entries(connected=g.db, accession=accession)
				description = entry['description']
				hgnc_gene_info = description
				
				# gene location
				gene_info = hdp.get_gene_info(hgnc)
				maploc = gene_info[1]
				maploc = maploc + ": (exon/intron boundary annotations = " + str(origin.ac) + ", alignment = " + alt_aln_method +")"

				# Splign alignments already available for NM yo NG accessions (refseq)
				h_ex = re.compile("^hgvs Exception:")
				if alt_aln_method == 'splign':
					# For this function switch alt_ac and tx_ac to keep the formatting in line
					# The try function in the function will catch missing refseq to NM alignments
					tx_exons = variantanalyser.functions.tx_exons(tx_ac=alt_ac, alt_ac=tx_ac, alt_aln_method=alt_aln_method, hdp=hdp)
					# return render_template('bootstrap/variantErrorLink.html', title=title, user=input, error=str(tx_exons), reason = 'giraffe', issue_link='none')
					if tx_exons == 'error':
						newmethod = 'splign'
					elif h_ex.search(str(tx_exons)):
						# hgvs_vt = hp.parse_hgvs_variant(variant)
						tx_exons = tx_exons.replace('hgvs Exception:', '')
						# Open a hgvs exception log file in append mode
						fo = open('/local/logs/hgvsExceptions.txt', 'a')
						excep = "%s -- %s -- %s\n" %(time.ctime(), tx_exons, variant)
						fo.write(excep)
						fo.close()
						error = 'Please inform UTA admin of the following error: ' + tx_exons
						issue_link = "https://bitbucket.org/biocommons/uta/issues?status=new&status=open"
						reason = "VariantValidator can not recover information for transcript " + str(alt_ac) + ' beacuse it is not available in the Universal Transcript Archive'
						return render_template('bootstrap/variantErrorLink.html', title=title, user=input, error=error, reason = reason, issue_link=issue_link)
					else:
						# Set up a list to contain the features
						features = []
						#crds = []
						#crde = []	
				
						# Pick out the five prime UTR and add to the features list
						cds_start = tx_exons[0]['alt_start_i']
						#crds.append(cds_start)
						five_utr = ac_seq[0:cds_start]
						features.append(five_utr)
				
						# Set exon counter to 0 and total exons
						exon_counter = 0
						total_exons = len(tx_exons)
				
						# Loop through tx_exons, note, because the reference will be a refseq gene, orientation always = sense
						for exons in tx_exons:
							# Pick out start and end coordinates for the current exon
							alt_start_i = int(tx_exons[exon_counter]['alt_start_i'])
							alt_end_i = int(tx_exons[exon_counter]['alt_end_i'])
							#crds.append(alt_start_i)
							#crde.append(alt_end_i)
							# Extract the exon sequence and append to features
							exon = ac_seq[alt_start_i:alt_end_i]
							exon = exon.upper()
							features.append(exon)
							# Extract the intron sequence and append to features
							if exon_counter + 1 < total_exons:
								# Pick out the start of the next exon
								intron_end = tx_exons[exon_counter + 1]['alt_start_i']
								#crde.append(intron_end)
								intron = ac_seq[alt_end_i:intron_end]
								features.append(intron)
							else:
								three_utr = ac_seq[alt_end_i:]
								features.append(intron)
							# Add 1 to the exon_counter
							exon_counter = exon_counter + 1
						#return render_template('bootstrap/variantError.html', title=title, user=input, error=str(crde))
				
				# Alignments of non-splign coordinates via genomic coordinates (antisense then sense)
				if alt_aln_method != 'splign' or newmethod == 'splign':
					# return render_template('bootstrap/variantError.html', title=title, user=input, error=ac_seq)
					# hgvs obnject the genomic coordinates and the refseqgene coordinates
					gencor = variantanalyser.functions.hgvs_genomic(variant=bridge, hp=hp)
					refcor = variantanalyser.functions.hgvs_refseq(variant, hp)
					ori = str(origin)
					oricor = variantanalyser.functions.coding(variant=ori, hp=hp)
					
					# Obtain the genomic refseq starting coordinates
					genstart = int(gencor.posedit.pos.start.base)
					refstart = int(refcor.posedit.pos.start.base)
					
					# return render_template('bootstrap/variantError.html', title=title, user=input, error=str(ofs))
					# Obtain the origin accession and tx exons to the genomic coordinates
					tx_ac = oricor.ac
					alt_ac = gencor.ac
					tx_exons = variantanalyser.functions.tx_exons(tx_ac, alt_ac, alt_aln_method, hdp)
					
					# Gene orientation - if antisense
					if tx_exons[0]['alt_strand'] == -1:
						ofs = genstart + refstart - 1
						# return render_template('bootstrap/variantError.html', title=title, user=input, error=ac_seq)
						# Set up a list to contain the features
						features = []
						#crds = []
						#crde = []	
						# Pick out the five prime UTR and add to the features list
						cds_start = ofs - tx_exons[0]['alt_end_i']
						#crds.append(cds_start)
						five_utr = ac_seq[0:cds_start]
						features.append(five_utr)
				
						# Set exon counter to 0 and total exons
						exon_counter = 0
						total_exons = len(tx_exons)
						# return render_template('bootstrap/variantError.html', title=title, user=input, error=ac_seq)
						# Loop through tx_exons, note, because the reference will be a refseq gene, orientation always = sense
						for exons in tx_exons:
				
							# Pick out start and end coordinates for the current exon
							alt_start_i = ofs - int(tx_exons[exon_counter]['alt_end_i'])
							alt_end_i = ofs - int(tx_exons[exon_counter]['alt_start_i'])
							#crds.append(alt_start_i)
							#crde.append(alt_end_i)
						
							# Extract the exon sequence and append to features
							exon = ac_seq[alt_start_i:alt_end_i]
							exon = exon.upper()
							features.append(exon)
					
							# Extract the intron sequence and append to features
							if exon_counter + 1 < total_exons:
								# Pick out the start of the next exon
								intron_end = ofs - tx_exons[exon_counter + 1]['alt_end_i']
								#crde.append(intron_end)
								intron = ac_seq[alt_end_i:intron_end]
								features.append(intron)
							else:
								three_utr = ac_seq[alt_end_i:]
								features.append(intron)
					
							# Add 1 to the exon_counter
							exon_counter = exon_counter + 1
						# return render_template('bootstrap/variantError.html', title=title, user=input, error=ac_seq)

					else: # If sense
						ofs = genstart - refstart

						# Set up a list to contain the features
						features = []
						#crds = []
						#crde = []	
						# Pick out the five prime UTR and add to the features list
						cds_start = tx_exons[0]['alt_start_i'] - ofs
						#crds.append(cds_start)
						five_utr = ac_seq[0:cds_start]
						features.append(five_utr)
				
						# Set exon counter to 0 and total exons
						exon_counter = 0
						total_exons = len(tx_exons)
				
						# Loop through tx_exons, note, because the reference will be a refseq gene, orientation always = sense
						for exons in tx_exons:
				
							# Pick out start and end coordinates for the current exon
							alt_start_i = int(tx_exons[exon_counter]['alt_start_i']) - ofs
							alt_end_i = int(tx_exons[exon_counter]['alt_end_i']) - ofs
							#crds.append(alt_start_i)
							#crde.append(alt_end_i)
						
							# Extract the exon sequence and append to features
							exon = ac_seq[alt_start_i:alt_end_i]
							exon = exon.upper()
							features.append(exon)
					
							# Extract the intron sequence and append to features
							if exon_counter + 1 < total_exons:
								# Pick out the start of the next exon
								intron_end = tx_exons[exon_counter + 1]['alt_start_i'] - ofs
								#crde.append(intron_end)
								intron = ac_seq[alt_end_i:intron_end]
								features.append(intron)
							else:
								three_utr = ac_seq[alt_end_i:]
								features.append(intron)
					
							# Add 1 to the exon_counter
							exon_counter = exon_counter + 1
						#return render_template('bootstrap/variantError.html', title=title, user=input, error=str(crde))

				# Join the features back into a sequence string
				ac_seq = ''.join(features)
				# return render_template('bootstrap/variantError.html', title=title, user=input, error=ac_seq, reason='the giraffe pood!')
				
	# Protein requests
	##################
	if type == ':p.':
		
		# Return hgvs object of the protein variant as reference
		if valid == 'true':
			# Use reference as coordinates are correct
			hgvs_reference = variantanalyser.functions.hgvs_protein(variant, hp)
			
			# Assign the accession string to on object
			tx_ac = hgvs_reference.ac
			
			# Validator doesn't work so extract the sequence only
			ref_seq = variantanalyser.functions.sequence_extractor(ac=tx_ac, hdp=hdp)
			
			# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=ref_seq)
			
			# Use the origin.ac to get the gene and map location coordinates
			tx_id_info = variantanalyser.functions.tx_identity_info(origin, hdp)

			# Get the gene name and CDS coordinates
			hgnc = tx_id_info[6]
			cds_start = tx_id_info[3]
			cds_end = tx_id_info[4]
	
			# Refseq accession
			tx_for_gene = variantanalyser.functions.tx_for_gene(hgnc, hdp)
			refseq_ac = variantanalyser.functions.ng_extract(tx_for_gene)
			
			# Get accurate descriptions from the relevant databases
			# RefSeq databases	
			# Gene description  - requires GenBank search to get all the required info, i.e. transcript variant ID
			# accession number
			hgvs_protein = variantanalyser.functions.hgvs_protein(variant, hp)	
			accession = hgvs_protein.ac
			# return render_template('bootstrap/variantError.html', title=title, user=input, error=accession)
			# Get the Entrez (GenBank) file
			#record = variantanalyser.functions.entrez_efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
			# Extract the information required
			#desc = record.description
			#record = ''
			#return render_template('bootstrap/variantError.html', title=title, user=input, error=desc)
			# hgnc_gene_info = (hgnc + ": " + desc)
			entry = variantanalyser.data.in_entries(connected=g.db, accession=accession)
			description = entry['description']
			hgnc_gene_info = description
				
			# gene location
			gene_info = hdp.get_gene_info(hgnc)
			maploc = gene_info[1]
			
			# Set the query as the DNA sequence for the genomic DNA using the origin
			#hgvs_query = variantanalyser.functions.rna(origin, evm, hp)
			c_typ = re.compile(':c.')
			if c_typ.search(origin):
				hgvs_query = variantanalyser.functions.rna(origin,evm, hp)
			else: 
				hgvs_query = hp.parse_hgvs_variant(origin)
			
			# Assign the accession string to on object
			alt_ac = hgvs_query.ac
			ac_seq = variantanalyser.functions.sequence_extractor(ac=alt_ac, hdp=hdp)
			# Set the query sequence to the reference so that it is carried through to the editing stages
			# Editing is based on the RNA level DNA sequence from the variant
			# This will be followed by translation
			# Only if the mutation is unambiguous
			ambiguity = re.compile("\?")
			if ambiguity.search(variant):
				return render_template('bootstrap/variantError.html', title=title, user=input, error='Variant sequences and alignments can only be returned for unambiguous variant descriptions', reason = 'Request Unavailable')
			else:	
				hgvs_reference = hgvs_query
			

	##############################
	# Reference Accession requests
	##############################
		
	# Request = reference, output the reference only fasta	
	if (output == 'reference') or (output == 'alignment'):
			
		# UTA protein sequences do not end * where as translations will
		if type == ':p.':
				# Add the termination codon *
				ref_seq = ref_seq + '*'
		
		# For alignments, need to maintain the reference sequence
		if output == 'alignment':
			if type != ':p.':
				ref_seq = ac_seq
		
		if output == 'reference':
			# To obtain the variant sequence if protein sequence
			# set variant ac_seq  to the coding sequence so it can be modified and translated
			if type == ':p.':
				ac_seq = ref_seq
		
			# Format a valid fasta title and add to the list
			# return render_template('bootstrap/variantError.html', title=title, user=input, error=ac_seq, reason='the giraffe weed')
			tx_ac_fasta_title = "> " + tx_ac +  "\n"
			# return render_template('bootstrap/variantError.html', title=title, user=input, error=ac_seq, reason='The giraffe is watching you' )
			#tx_ac_fasta = variantanalyser.links.reference(tx_ac_fasta_title, ac_seq)
			tx_ac_fasta = []
			tx_ac_fasta.append(tx_ac_fasta_title)
			
			# Split the sequence strings into lists of 60 base chunks
			seq_list = variantanalyser.links.nsplit(ac_seq, n=60)
			# loop and append as required
			for lines in seq_list:
				tx_ac_fasta.append(lines + '\n')
			
			# Add web page headers
			seq_type = "Fasta sequence"
			tx_ac = tx_ac + " " + message
			
			return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=tx_ac_fasta, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method='none') 


	############################
	# Variant accession Requests
	############################
		
	# Request = reference, output the reference only fasta	
	if (output == 'variant') or (output == 'alignment'):

		if output == 'variant':
			# Format a valid fasta title and add to the list
			tx_ac_fasta_title = "> " + input + " : Variant alignment method = " + alt_aln_method + "\n"

		if output == 'alignment':
			# Format a valid fasta title and add to the list
			tx_ac_fasta_title = "> " + input + " aligned to " + tx_ac + " : Variant alignment method = " + alt_aln_method + "\n"
		
		# Collect the interval starting and ending positions (based on RNA) also the edit (based on coding)
		interval_start = int(hgvs_reference.posedit.pos.start.base)
		interval_end = int(hgvs_reference.posedit.pos.end.base)
		# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=str(hgvs_reference))
		
		# Handle coding sequence offsetting
		if (type == ':c.'):
			edit = str(hgvs_query.posedit.edit)
		else:
			edit = str(hgvs_reference.posedit.edit)
		
		# SIMPLE DELETIONS - assumes validation has been completed to check sequence where required
		###########################################################################################
		# Look for deletion at the end of the edit string
		deletion = re.compile("del$")
		deletion_b = re.compile("del[GATCgatc]+$")
			
		if deletion.search(edit) or deletion_b.search(edit):
			#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta='giraffe')
			tx_ac_fasta = variantanalyser.links.sim_del(tx_ac_fasta_title, ac_seq, interval_start, interval_end, variant)
			
			# Translate the returned sequence if protein was requested
			if type == ':p.':
				ed_seq = tx_ac_fasta[1]
				tx_ac_fasta[1] = ''
				tx_ac_fasta[1] = variantanalyser.links.translate(ed_seq, cds_start)
				if tx_ac_fasta == 'error':
					return render_template('bootstrap/variantError.html', title=title, user=input, error='Translation start codon (ATG) not found at CDS start', reason = 'Request Unavailable')
			
			if output == 'variant':
				seq_type = "Fasta sequence"
				alt_aln_method = alt_aln_method + " " + message
				tx_ac = variant + " " + message
				
				# Split the sequence strings into lists of 60 base chunks
				ac_seq = tx_ac_fasta.pop()
				seq_list = variantanalyser.links.nsplit(ac_seq, n=60)
				# loop and append as required
				for lines in seq_list:
					tx_ac_fasta.append(lines + '\n')
		
				return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=tx_ac_fasta, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)

			if output == 'alignment':
				#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta='Beware of the giraffe')
				# Collect variant sequence
				var_seq = tx_ac_fasta[1]
				var_seq = var_seq.strip()
				var_seq = var_seq.replace('\n', '')
				
				seq_type = "Sequence alignment"
				alt_aln_method = alt_aln_method + " " + message
				tx_ac = input + " to " + tx_ac + " " + message

				# Alignment module applied here
				aln = []
				aln.append(tx_ac_fasta_title)
				
				# Re position the interval start and end for the protein variant
				# flag for edit_type
				if type == ':p.':
					#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta='giraffe')
					# hgvs object of protein variant
					var_p = variantanalyser.functions.hgvs_protein(variant, hp)
					interval_start = int(var_p.posedit.pos.start.base)
					interval_end = int(var_p.posedit.pos.end.base)
					# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=str(interval_end), seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
					hgvs_protein = variantanalyser.functions.hgvs_protein(variant, hp) 
					pr_ed = str(hgvs_protein.posedit.edit)
					# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=pr_ed)
					dels = re.compile('^del$')
					dels_b = re.compile('^del[GATCgatc]+$')
					delins = re.compile('delins')
					delins_b = re.compile("del[GATCgatc]+ins")
					if dels.search(pr_ed):
						frame = 'true'
						edit_type = ''
					elif dels_b.search(pr_ed):
						frame = 'true'
						edit_type = ''
					elif delins.search(pr_ed):
						frame = 'true'
						edit_type = 'delins'
						ins_len = re.search(r"(delins[A-Za-z]+)", variant)
						ins_len = ins_len.group(1)
						ins_len = ins_len.replace('delins', '')
						ins_len = len(ins_len) / 3
					elif delins_b.search(pr_ed):
						ins_len = re.search(r"(delins[A-Za-z]+)", variant)
						ins_len = ins_len.group(1)
						ins_len = ins_len.replace('delins', '')
						frame = 'true'
						edit_type = 'delins'
					else:
						frame = 'false'
						edit_type = 'sim_prot'				
				
				# Deal with intronic positions
				# Create sequences to be aligned.
				# Intronic positions need to be dealt with
				pl = re.compile('\+')
				mi = re.compile('\-')
				if (pl.search(variant)) or (mi.search(variant)):
					if interval_start == interval_end:	
						#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta='Beware of the giraffe')
						# Create sequences to be aligned.
						begin = interval_start-plus_minus
						if begin < 0:
							begin = 0
						end = interval_end + plus_minus
						if end > interval_end + len(ac_seq[interval_end:]):
							end = interval_end + len(ac_seq[interval_end:])
							# poo = len(ac_seq[interval_end:])
						align1 = var_seq[begin:end]
						align2 = ref_seq[begin:end]
						edit_type = 'sim_align'
						#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=edit_type)
					else:
						# Create sequences to be aligned.
						#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta='Beware of the giraffe')
						begin = interval_start-plus_minus
						if begin < 0:
							begin = 0
						end = interval_end + plus_minus
						if end > interval_end + len(ac_seq[interval_end:]):
							end = interval_end + len(ac_seq[interval_end:])
							# poo = len(ac_seq[interval_end:])
						align1 = var_seq[begin:end -(interval_end - interval_start) -1]
						align2 = ref_seq[begin:end]
						edit_type = 'sim_del'
						# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=str(interval_end), seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
				else:
					# Create sequences to be aligned.
					begin = interval_start-plus_minus
					if begin < 0:
						begin = 0
					end = interval_end + plus_minus
					#giraffe = end
					if end > interval_end + len(ac_seq[interval_end:]):
						end = interval_end + len(ac_seq[interval_end:])
						# poo = len(ac_seq[interval_end:])				
					align1 = var_seq[begin:end -(interval_end - interval_start) -1]
					align2 = ref_seq[begin:end]
					#aligns = align1 + '\n\n' + align2
					#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=aligns, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
					# Assign edit type if not already done
					if edit_type == '':
						edit_type = 'sim_del'
					# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=align2, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)

				# Send data to function that creates the alignment
				alignment = variantanalyser.links.format_alignment(align1, align2, edit_type, interval_start, interval_end, begin, end, ins_len, type, frame)
				# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=alignment, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)

				# Capture the length of the sequences
				seq_len = len(alignment[0])
				lin_len = len(alignment[1])
				ref_len = len(alignment[2])
				
				# All strings should be the same length
				if (seq_len == lin_len) and (seq_len == ref_len) and (lin_len == ref_len):
					# Split the strings into lists of 60 base chunks
					seq_list = variantanalyser.links.nsplit(alignment[0], n=60)
					lin_list = variantanalyser.links.nsplit(alignment[1], n=60)
					ref_list = variantanalyser.links.nsplit(alignment[2], n=60)
					
					# element counter
					chunk = 0
					
					# loop and append as required
					for lines in seq_list:
						aln.append(seq_list[chunk] + "\n")
						aln.append(lin_list[chunk] + "\n")
						aln.append(ref_list[chunk] + "\n")
						chunk = chunk + 1
								
					return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=aln, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
			
		# SIMPLE INSERTIONS: assumes validation has been completed to check sequence where required
		###########################################################################################
		# Look for ins within the string but omit delins
		insertion = re.compile('ins[GATCUgatcu]+$')	
		deletion = re.compile("del")
		if insertion.search(edit):
			if deletion.search(edit):
				pass
			else:
				tx_ac_fasta_dict = variantanalyser.links.sim_ins(tx_ac_fasta_title, ac_seq, interval_start, interval_end, edit, variant)
				# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=tx_ac_fasta_dict)
				tx_ac_fasta = tx_ac_fasta_dict['tx_ac_fasta']
				ins_len = tx_ac_fasta_dict['edit_len']
				
				# Translate the returned sequence if protein was requested
				if type == ':p.':
					ed_seq = tx_ac_fasta[1]
					#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=ed_seq)
					tx_ac_fasta[1] = ''
					tx_ac_fasta[1] = variantanalyser.links.translate(ed_seq, cds_start)
					if tx_ac_fasta == 'error':
						return render_template('bootstrap/variantError.html', title=title, user=input, error='Translation start codon (ATG) not found at CDS start', reason = 'Request Unavailable')
			
				if output == 'variant':
					seq_type = "Fasta sequence"
					alt_aln_method = alt_aln_method + " " + message
					tx_ac = variant + " " + message
					
					# Split the sequence strings into lists of 60 base chunks
					ac_seq = tx_ac_fasta.pop()
					seq_list = variantanalyser.links.nsplit(ac_seq, n=60)
					# loop and append as required
					for lines in seq_list:
						tx_ac_fasta.append(lines + '\n')
		
					return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=tx_ac_fasta, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)

				if output == 'alignment':
				
					# Collect variant sequence
					var_seq = tx_ac_fasta[1]
					var_seq = var_seq.strip()
					var_seq = var_seq.replace('\n', '')
				
					seq_type = "Sequence alignment"
					alt_aln_method = alt_aln_method + " " + message
					tx_ac = input + " to " + tx_ac + " " + message

					# Alignment module applied here
					aln = []
					aln.append(tx_ac_fasta_title)
			
					# Re position the interval start and end for the protein variant
					if type == ':p.':
						# hgvs object of protein variant
						var_p = variantanalyser.functions.hgvs_protein(variant, hp)
						interval_start = int(var_p.posedit.pos.start.base)
						interval_end = int(var_p.posedit.pos.end.base)
						# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=str(interval_end), seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
						ins = re.compile('\dins[A-Za-z]+')
						if ins.search(variant):
							ins_len = re.search(r"(ins[A-Za-z]+)", variant)
							ins_len = ins_len.group(1)
							ins_len = ins_len.replace('ins', '')
							ins_len = len(ins_len) / 3
							frame = 'true'
						else:
							ins_len = 0
							edit_type = 'sim_prot'
							frame = 'false'			

					# Deal with intronic positions
					# Create sequences to be aligned.
					# Intronic positions need to be dealt with
					pl = re.compile('\+')
					mi = re.compile('\-')
					if (pl.search(variant)) or (mi.search(variant)):
						if interval_start == interval_end:	
							begin = interval_start-plus_minus
							if begin < 0:
								begin = 0
							end = interval_end + plus_minus
							if end > interval_end + len(ac_seq[interval_end:]):
								end = interval_end + len(ac_seq[interval_end:])
								# poo = len(ac_seq[interval_end:])
							align1 = var_seq[begin:end]
							align2 = ref_seq[begin:end]
							edit_type = 'sim_align'
							#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=edit_type)
						else:
							# Create sequences to be aligned.
							begin = interval_start-plus_minus
							if begin < 0:
								begin = 0
							end = interval_end + plus_minus
							if end > interval_end + len(ac_seq[interval_end:]):
								end = interval_end + len(ac_seq[interval_end:])
								# poo = len(ac_seq[interval_end:])
							align1 = var_seq[begin:end]
							align2 = ref_seq[begin:end]
							edit_type = 'sim_del'
							# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=str(interval_end), seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
					else:
						# Create sequences to be aligned.
						begin = interval_start-plus_minus
						if begin < 0:
							begin = 0
						end = interval_end + plus_minus
						if end > interval_end + len(ac_seq[interval_end:]):
							end = interval_end + len(ac_seq[interval_end:])
							# poo = len(ac_seq[interval_end:])
						align1 = var_seq[begin:end]
						align2 = ref_seq[begin:end - ins_len]
						if edit_type == '':
							edit_type = 'sim_ins'
						# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=str(interval_end), seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
				
					# Send data to function that creates the alignment
					alignment = variantanalyser.links.format_alignment(align1, align2, edit_type, interval_start, interval_end, begin, end, ins_len, type, frame)
					#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=alignment, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
				
					# Capture the length of the sequences
					seq_len = len(alignment[0])
					lin_len = len(alignment[1])
					ref_len = len(alignment[2])
				
					# All strings should be the same length
					if (seq_len == lin_len) and (seq_len == ref_len) and (lin_len == ref_len):
						# Split the strings into lists of 60 base chunks
						seq_list = variantanalyser.links.nsplit(alignment[0], n=60)
						lin_list = variantanalyser.links.nsplit(alignment[1], n=60)
						ref_list = variantanalyser.links.nsplit(alignment[2], n=60)
					
						# element counter
						chunk = 0
					
						# loop and append as required
						for lines in seq_list:
							aln.append(seq_list[chunk] + "\n")
							aln.append(lin_list[chunk] + "\n")
							aln.append(ref_list[chunk] + "\n")
							chunk = chunk + 1
								
						return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=aln, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)

		
		# DELINS - assumes validation has been completed to check sequence where required
		#################################################################################
		# Look for deletion at the end of the edit string
		delins = re.compile("del[GATCgatc]+ins")
		delins_b = re.compile("delins")
		#return render_template('bootstrap/variantError.html', title=title, user=input, error=edit, reason = 'Giraffe')
			
		if delins.search(edit) or delins_b.search(edit):
			# return render_template('bootstrap/variantError.html', title=title, user=input, error=edit, reason = 'Giraffe')
			# Compile the fasta and recover the edit length
			tx_ac_fasta_dict = variantanalyser.links.delins(tx_ac_fasta_title, ac_seq, interval_start, interval_end, edit, variant, type)
			#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=tx_ac_fasta_dict)
			tx_ac_fasta = tx_ac_fasta_dict['tx_ac_fasta']
			ins_len = tx_ac_fasta_dict['edit_len']

			# Translate the returned sequence if protein was requested
			if type == ':p.':
				ed_seq = tx_ac_fasta[1]
				tx_ac_fasta[1] = ''
				tx_ac_fasta[1] = variantanalyser.links.translate(ed_seq, cds_start)
				if tx_ac_fasta == 'error':
					return render_template('bootstrap/variantError.html', title=title, user=input, error='Translation start codon (ATG) not found at CDS start', reason = 'Request Unavailable')
			
			if output == 'variant':
				seq_type = "Fasta sequence"
				alt_aln_method = alt_aln_method + " " + message
				tx_ac = variant + " " + message
			
				# Split the sequence strings into lists of 60 base chunks
				ac_seq = tx_ac_fasta.pop()
				seq_list = variantanalyser.links.nsplit(ac_seq, n=60)
				# loop and append as required
				for lines in seq_list:
					tx_ac_fasta.append(lines + '\n')
				
				
				return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=tx_ac_fasta, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
			
			if output == 'alignment':
				
				# Collect variant sequence
				var_seq = tx_ac_fasta[1]
				var_seq = var_seq.strip()
				var_seq = var_seq.replace('\n', '')
				
				seq_type = "Sequence alignment"
				alt_aln_method = alt_aln_method + " " + message
				tx_ac = input + " to " + tx_ac + " " + message

				# Alignment module applied here
				aln = []
				aln.append(tx_ac_fasta_title)
			
				# Re position the interval start and end for the protein variant
				if type == ':p.':
					# hgvs object of protein variant
					# return render_template('bootstrap/variantError.html', title=title, user=input, error='Bum', reason = 'Giraffe')
					var_p = variantanalyser.functions.hgvs_protein(variant, hp)
					interval_start = int(var_p.posedit.pos.start.base)
					interval_end = int(var_p.posedit.pos.end.base)
					delin = re.compile('\ddelins')
					ins = re.compile('\dins')
					# Simple fs search
					fs = re.compile('fs')
					if delin.search(variant):
						ins_len = re.search(r"(delins[A-Za-z]+)", variant)
						ins_len = ins_len.group(1)
						ins_len = ins_len.replace('delins', '')
						ins_len = len(ins_len) / 3
						frame = 'true'
						edit_type = ''
						# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=str(ins_len), seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
					elif ins.search(variant):
						ins_len = re.search(r"(ins[A-Za-z]+)", variant)
						ins_len = ins_len.group(1)
						ins_len = ins_len.replace('ins', '')
						ins_len = len(ins_len) / 3
						frame = 'true'
						edit_type = 'sim_ins'
					elif fs.search(str(var_p)):
						# return render_template('bootstrap/variantError.html', title=title, user=input, error=var_p, reason = 'Giraffe')
						edit_type = 'sim_prot'
						frame = 'false'
					else:
						ins_len = 0
						edit_type = 'sim_prot'
						frame = 'false'
					
						# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=frame, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)

				# Create sequences to be aligned.
				begin = interval_start-plus_minus
				if begin < 0:
					begin = 0
				end = interval_end + plus_minus
				if end > interval_end + len(ac_seq[interval_end:]):
					end = interval_end + len(ac_seq[interval_end:])
					# poo = len(ac_seq[interval_end:])
				
				# pl mi intronic only intervals will have a zero insert length
				# Set to simple alignment 
				if ins_len == 0 and edit_type != 'sim_prot':
					align1 = var_seq[begin:end]		
					align2 = ref_seq[begin:end]
					edit_type = 'sim_align'
				
				else:
					align1 = var_seq[begin:end + ins_len]		
					align2 = ref_seq[begin:end + (interval_end - interval_start +1)]
				
				# Set the edit_type if not already set
				if edit_type == '':
					edit_type = 'delins'

				alignment = variantanalyser.links.format_alignment(align1, align2, edit_type, interval_start, interval_end, begin, end, ins_len, type, frame)
				#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=alignment, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
				
				# Capture the length of the sequences
				seq_len = len(alignment[0])
				lin_len = len(alignment[1])
				ref_len = len(alignment[2])
				
				# All strings should be the same length
				if (seq_len == lin_len) and (seq_len == ref_len) and (lin_len == ref_len):
					# Split the strings into lists of 60 base chunks
					seq_list = variantanalyser.links.nsplit(alignment[0], n=60)
					lin_list = variantanalyser.links.nsplit(alignment[1], n=60)
					ref_list = variantanalyser.links.nsplit(alignment[2], n=60)
					
					# element counter
					chunk = 0
					
					# loop and append as required
					for lines in seq_list:
						aln.append(seq_list[chunk] + "\n")
						aln.append(lin_list[chunk] + "\n")
						aln.append(ref_list[chunk] + "\n")
						chunk = chunk + 1
								
					return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=aln, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)


		# DUPLICATIONS - Note, alignment duplications might be merged with insertions
		##############
			
		dup = re.compile ("^dup")
			
		if dup.search(edit):
			tx_ac_fasta_dict = variantanalyser.links.dupn(tx_ac_fasta_title, ac_seq, interval_start, interval_end, edit, variant)
			#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=tx_ac_fasta_dict)
			tx_ac_fasta = tx_ac_fasta_dict['tx_ac_fasta']
			ins_len = tx_ac_fasta_dict['edit_len']	
			#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=str(ins_len))
			
			# Translate the returned sequence if protein was requested
			if type == ':p.':
				ed_seq = tx_ac_fasta[1]
				#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=ed_seq)
				tx_ac_fasta[1] = ''
				tx_ac_fasta[1] = variantanalyser.links.translate(ed_seq, cds_start)
				if tx_ac_fasta == 'error':
					return render_template('bootstrap/variantError.html', title=title, user=input, error='Translation start codon (ATG) not found at CDS start', reason = 'Request Unavailable')
			
			if output == 'variant':
				seq_type = "Fasta sequence"
				alt_aln_method = alt_aln_method + " " + message
				tx_ac = variant + " " + message
				
				# Split the sequence strings into lists of 60 base chunks
				ac_seq = tx_ac_fasta.pop()
				seq_list = variantanalyser.links.nsplit(ac_seq, n=60)
				# loop and append as required
				for lines in seq_list:
					tx_ac_fasta.append(lines + '\n')
		
				return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=tx_ac_fasta, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)

			if output == 'alignment':
				
				# Collect variant sequence
				var_seq = tx_ac_fasta[1]
				var_seq = var_seq.strip()
				var_seq = var_seq.replace('\n', '')
				
				seq_type = "Sequence alignment"
				alt_aln_method = alt_aln_method + " " + message
				tx_ac = input + " to " + tx_ac + " " + message

				# Alignment module applied here
				aln = []
				aln.append(tx_ac_fasta_title)
			
				# Re position the interval start and end for the protein variant
				if type == ':p.':
					# hgvs object of protein variant
					var_p = variantanalyser.functions.hgvs_protein(variant, hp)
					interval_start = int(var_p.posedit.pos.start.base)
					interval_end = int(var_p.posedit.pos.end.base)
					# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=str(interval_end), seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
					ins = re.compile('\dins[A-Za-z]+')
					duplication = re.compile('dup')
					if ins.search(variant):
						ins_len = re.search(r"(ins[A-Za-z]+)", variant)
						ins_len = ins_len.group(1)
						ins_len = ins_len.replace('ins', '')
						ins_len = len(ins_len) / 3
						frame = 'true'
					if duplication.search(variant):
						ins_len = ins_len/3
						frame = 'true'
					else:
						ins_len = 0
						edit_type = 'sim_prot'
						frame = 'false'			

				# Deal with intronic positions
				# Create sequences to be aligned.
				# Intronic positions need to be dealt with
				pl = re.compile('\+')
				mi = re.compile('\-')
				if (pl.search(variant)) or (mi.search(variant)):
					if interval_start == interval_end:	
						begin = interval_start-plus_minus
						if begin < 0:
							begin = 0
						end = interval_end + plus_minus
						if end > interval_end + len(ac_seq[interval_end:]):
							end = interval_end + len(ac_seq[interval_end:])
							# poo = len(ac_seq[interval_end:])
						align1 = var_seq[begin:end]
						align2 = ref_seq[begin:end]
						edit_type = 'sim_align'
						#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=edit_type)
					else:
						# Create sequences to be aligned.
						begin = interval_start-plus_minus
						if begin < 0:
							begin = 0
						end = interval_end + plus_minus
						if end > interval_end + len(ac_seq[interval_end:]):
							end = interval_end + len(ac_seq[interval_end:])
							# poo = len(ac_seq[interval_end:])
						align1 = var_seq[begin:end]
						align2 = ref_seq[begin:end]
						edit_type = 'sim_del'
						# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=str(interval_end), seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
				else:
					# Create sequences to be aligned.
					begin = interval_start-plus_minus
					if begin < 0:
						begin = 0
					end = interval_end + plus_minus
					if end > interval_end + len(ac_seq[interval_end:]):
						end = interval_end + len(ac_seq[interval_end:])
						# poo = len(ac_seq[interval_end:])
					align1 = var_seq[begin:end]
					align2 = ref_seq[begin:end - ins_len]
					if edit_type == '':
						edit_type = 'dupn'
					#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=alignment, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
				
				# Send data to function that creates the alignment
				alignment = variantanalyser.links.format_alignment(align1, align2, edit_type, interval_start, interval_end, begin, end, ins_len, type, frame)
				#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=alignment, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
				
				# Capture the length of the sequences
				seq_len = len(alignment[0])
				lin_len = len(alignment[1])
				ref_len = len(alignment[2])
				
				# All strings should be the same length
				if (seq_len == lin_len) and (seq_len == ref_len) and (lin_len == ref_len):
					# Split the strings into lists of 60 base chunks
					seq_list = variantanalyser.links.nsplit(alignment[0], n=60)
					lin_list = variantanalyser.links.nsplit(alignment[1], n=60)
					ref_list = variantanalyser.links.nsplit(alignment[2], n=60)
					
					# element counter
					chunk = 0
					
					# loop and append as required
					for lines in seq_list:
						aln.append(seq_list[chunk] + "\n")
						aln.append(lin_list[chunk] + "\n")
						aln.append(ref_list[chunk] + "\n")
						chunk = chunk + 1
								
					return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=aln, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
	
		
		# Inversions - assumes validation has been completed to check sequence where required
		#####################################################################################
		# Look for deletion at the end of the edit string
		inversion = re.compile("inv")
			
		if inversion.search(edit):

			#return render_template('bootstrap/variantError.html', title=title, user=input, error=edit, reason = 'Giraffe')

			
			if type == ':n.' or type == ':c.':
				hgvs_variant = hp.parse_hgvs_variant(origin) # Origin is the coding variant used to populate the variant table
				# Need to obtain the cds_start
				inf = variantanalyser.functions.tx_identity_info(variant, hdp)
				cds_start = inf[3]
			elif type == ':p.':
				hgvs_variant = hp.parse_hgvs_variant(origin) # Origin is the coding variant used to populate the variant table
				inf = variantanalyser.functions.tx_identity_info(origin, hdp)
				cds_start = inf[3]
			else:
				hgvs_variant = hp.parse_hgvs_variant(variant)
				cds_start = 0
			
			# INFORMATION REQUIRED FOR THE EDIT
			# Collect the deleted sequence
			del_seq = str(hgvs_variant.posedit.edit)
			del_seq = del_seq.replace('inv', '')
			# Make the inverted sequence
			my_seq = Seq(del_seq)
			inv_seq = my_seq.reverse_complement()

			# Make the variant sequence
			if type == ':p.':
				# Re Extract the reference coding sequence from the UTA database
				try:
					tp_seq = variantanalyser.functions.sequence_extractor(ac=hgvs_variant.ac, hdp=hdp)
				except hgvs.exceptions.HGVSError as e:
					error = e
				if error != 'false':
					# Open a hgvs exception log file in append mode
					fo = open('/local/logs/hgvsExceptions.txt', 'a')
					excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
					fo.write(excep)
					fo.close()
					return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Variant description')
				else: 
					pass
				fasta_seq = variantanalyser.links.coding_inversion(tp_seq, del_seq, inv_seq, interval_start=hgvs_variant.posedit.pos.start.base+cds_start, interval_end=hgvs_variant.posedit.pos.end.base+cds_start)
			else:
				# Pass the sequence
				ref_seq = ac_seq
				fasta_seq = variantanalyser.links.coding_inversion(ac_seq, del_seq, inv_seq, interval_start=hgvs_variant.posedit.pos.start.base+cds_start, interval_end=hgvs_variant.posedit.pos.end.base+cds_start)
				tx_ac_fasta = [tx_ac_fasta_title, fasta_seq]
				#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=tx_ac_fasta, seq_type=fasta_seq, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
			
			# tx_ac_fasta = tx_ac_fasta_dict['tx_ac_fasta']
			# ins_len = tx_ac_fasta_dict['edit_len']

			# Translate the returned sequence if protein was requested
			if type == ':p.':
				fasta_seq = variantanalyser.links.translate(fasta_seq, cds_start)
				#return render_template('bootstrap/variantError.html', title=title, user=input, error='Translation start codon (ATG) not found at CDS start', reason = )
				if fasta_seq == 'error':
					return render_template('bootstrap/variantError.html', title=title, user=input, error='Translation start codon (ATG) not found at CDS start', reason = 'Request Unavailable')
			
			if output == 'variant':
				seq_type = "Fasta sequence"
				alt_aln_method = alt_aln_method + " " + message
				tx_ac = variant + " " + message
			
				# Split the sequence strings into lists of 60 base chunks
				ac_seq = tx_ac_fasta.pop()
				seq_list = variantanalyser.links.nsplit(ac_seq, n=60)
				# loop and append as required
				for lines in seq_list:
					tx_ac_fasta.append(lines + '\n')
				
				
				return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=tx_ac_fasta, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
			
			if output == 'alignment':
				
				# Collect variant sequence
				var_seq = fasta_seq
				var_seq = var_seq.strip()
				# var_seq = var_seq.replace('\n', '')
				
				seq_type = "Sequence alignment"
				alt_aln_method = alt_aln_method + " " + message
				tx_ac = input + " to " + tx_ac + " " + message

				# Alignment module applied here
				aln = []
				aln.append(tx_ac_fasta_title)
			
				if type == ':p.':
					# Get the start and end positions
					#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=crds, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
					hgvs_protein = variantanalyser.functions.hgvs_protein(variant, hp)
					interval_start = hgvs_protein.posedit.pos.start.base
					interval_end = hgvs_protein.posedit.pos.end.base
					# return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=str(interval_end), seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
				else:
					interval_start = hgvs_variant.posedit.pos.start.base + cds_start
					interval_end = hgvs_variant.posedit.pos.end.base + cds_start
				#crds = str(interval_start) + '  ' + str(interval_end)
				#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=crds, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
				# Create sequences to be aligned.
				begin = interval_start-plus_minus
				if begin < 0:
					begin = 0
				end = interval_end + plus_minus
				if end > interval_end + len(fasta_seq[interval_end:]):
					end = interval_end + len(fasta_seq[interval_end:])
					# poo = len(ac_seq[interval_end:])
				
				# pl mi intronic only intervals will have a zero insert length
				# Set to simple alignment 
				#if ins_len == 0 and edit_type != 'sim_prot':
				#	align1 = var_seq[begin:end]		
				#	align2 = ref_seq[begin:end]
				#	edit_type = 'sim_align'
				
				else:
					align1 = fasta_seq[begin:end]		
					align2 = ref_seq[begin:end]
					aligns = align1 + '\n\n' + align2
				
				#crds = str(begin) + '  ' + str(end)
				#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=crds, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)

				# Set the edit_type if not already set - Note, Inversions will always be equal therefore point align is fine
				edit_type = 'point'
				#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=aligns, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
				alignment = variantanalyser.links.format_alignment(align1, align2, edit_type, interval_start, interval_end, begin, end, ins_len, type, frame)
				#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=alignment, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)
				
				# Capture the length of the sequences
				seq_len = len(alignment[0])
				lin_len = len(alignment[1])
				ref_len = len(alignment[2])
				
				# All strings should be the same length
				if (seq_len == lin_len) and (seq_len == ref_len) and (lin_len == ref_len):
					# Split the strings into lists of 60 base chunks
					seq_list = variantanalyser.links.nsplit(alignment[0], n=60)
					lin_list = variantanalyser.links.nsplit(alignment[1], n=60)
					ref_list = variantanalyser.links.nsplit(alignment[2], n=60)
					
					# element counter
					chunk = 0
					
					# loop and append as required
					for lines in seq_list:
						aln.append(seq_list[chunk] + "\n")
						aln.append(lin_list[chunk] + "\n")
						aln.append(ref_list[chunk] + "\n")
						chunk = chunk + 1
								
					return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=aln, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method)


		# SIMPLE POINT MUTATIONS
		########################
		pm = re.compile("[GATCUgatcu]\>[GATCUgatcu]")
			
		if pm.search(edit):
			# Find pattern and assign to a variable
			tx_ac_fasta_dict = variantanalyser.links.point(tx_ac_fasta_title, ac_seq, interval_start, edit, variant)
			tx_ac_fasta = tx_ac_fasta_dict['tx_ac_fasta']
			flag = tx_ac_fasta_dict['flag']
			
			# Translate the returned sequence if protein was requested
			if type == ':p.':
				ed_seq = tx_ac_fasta[1]
				tx_ac_fasta[1] = ''
				tx_ac_fasta[1] = variantanalyser.links.translate(ed_seq, cds_start)
				if tx_ac_fasta == 'error':
					return render_template('bootstrap/variantError.html', title=title, user=input, error='Translation start codon (ATG) not found at CDS start', reason = 'Request Unavailable')
			
			if output == 'variant':
				seq_type = "Fasta sequence"
				alt_aln_method = alt_aln_method + " " + message
				tx_ac = variant + " " + message
				
				# Split the sequence strings into lists of 60 base chunks
				ac_seq = tx_ac_fasta.pop()
				seq_list = variantanalyser.links.nsplit(ac_seq, n=60)
				# loop and append as required
				for lines in seq_list:
					tx_ac_fasta.append(lines + '\n')
			
				return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=tx_ac_fasta, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method, flag=flag)
			
			if output == 'alignment':
				
				# Collect variant sequence
				var_seq = tx_ac_fasta[1]
				var_seq = var_seq.strip()
				var_seq = var_seq.replace('\n', '')
				
				seq_type = "Sequence alignment"
				alt_aln_method = alt_aln_method + " " + message
				tx_ac = input + " to " + tx_ac + " " + message

				# Alignment module applied here
				aln = []
				aln.append(tx_ac_fasta_title)
				
				# Re position the interval start and end for the protein variant
				if type == ':p.':
					eq = re.compile('\=')
					# hgvs object of protein variant
					var_p = variantanalyser.functions.hgvs_protein(variant, hp)
					interval_start = int(var_p.posedit.pos.start.base)
					interval_end = int(var_p.posedit.pos.end.base)
					ter = re.compile('Ter')
					if (ter.search(variant)):
						frame = 'false'
						edit_type = 'sim_prot'
				
				# Create sequences to be aligned.
				begin = interval_start-plus_minus
				if begin < 0:
					begin = 0
				end = interval_end + plus_minus
				if end > interval_end + len(ac_seq[interval_end:]):
					end = interval_end + len(ac_seq[interval_end:])
					# poo = len(ac_seq[interval_end:])

				align1 = var_seq[begin:end]
				align2 = ref_seq[begin:end]
				if edit_type == '':
					edit_type = 'point'

				# Deal with intronic positions
				# Create sequences to be aligned.
				# Intronic positions need to be dealt with
				pl = re.compile('\+')
				mi = re.compile('\-')
				if (pl.search(variant)) or (mi.search(variant)):
					edit_type = 'sim_align'

				# Send data to function that creates the alignment
				alignment = variantanalyser.links.format_alignment(align1, align2, edit_type, interval_start, interval_end, begin, end, ins_len, type, frame)
				#return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=alignment, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method, flag=flag)
				
				# Capture the length of the sequences
				seq_len = len(alignment[0])
				lin_len = len(alignment[1])
				ref_len = len(alignment[2])
				
				# All strings should be the same length
				if (seq_len == lin_len) and (seq_len == ref_len) and (lin_len == ref_len):
					# Split the strings into lists of 60 base chunks
					seq_list = variantanalyser.links.nsplit(alignment[0], n=60)
					lin_list = variantanalyser.links.nsplit(alignment[1], n=60)
					ref_list = variantanalyser.links.nsplit(alignment[2], n=60)
					
					# element counter
					chunk = 0
					
					# loop and append as required
					for lines in seq_list:
						aln.append(seq_list[chunk] + "\n")
						aln.append(lin_list[chunk] + "\n")
						aln.append(ref_list[chunk] + "\n")
						chunk = chunk + 1
								
				return render_template('bootstrap/variantLinks.html', title=title, tx_ac_fasta=aln, seq_type=seq_type, tx_ac=tx_ac, alt_aln_method=alt_aln_method, flag=flag)
	
	# If a variant has not been captured, return an error response
	# Open a hgvs exception log file in append mode
	fo = open('/local/logs/hgvsExceptions.txt', 'a')
	error = 'Variant not recognised'
	excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
	fo.write(excep)
	fo.close()
	return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Something went wrong. Please contact admin stating the input variant and time of input')			
			

# Route for links loader page
@app.route('/links_loading/', methods=['POST', 'GET'])
def links_loading():
	
	title = "loading"
	user = "loading the requested data"
	# Handle form data
	# Example input URL
	# https://www22.lamp.le.ac.uk/hgvs/variantlinks/?variant=NM_000088.3%3Ac.589G%3ET&alignment=splign&alt_acc=NG_007400.1:g.8638G>T&request=reference
	input = request.args.get('variant', '')
	alt_aln_method = request.args.get('alignment', '')
	output = request.args.get('output', '')
	origin = request.args.get('origin', '')
	bridge = request.args.get('bridge', '')
	
	# Re assemble the URL
	# url = '/variantlinks/?variant='+variant+'&alignment='+alt_aln_method+'&output='+output+'&origin='+origin+'&bridge='+bridge
	
	# Tell the loading page which loader is being called
	loader = 'links_loading'
	# Sent the data to an auto submit form loading 
	return render_template('bootstrap/loading.html', title=title, user=user, loader=loader, input=input, alt_aln_method=alt_aln_method, output=output, origin=origin, bridge=bridge)


# Route for adding loader page
@app.route('/add_loading/', methods=['POST', 'GET'])
def add_loading():
	
	title = "db_adding"
	user = "Adding the transcript description to our database"
	# Handle form data
	# Example input URL
	# https://www22.lamp.le.ac.uk/hgvs/variantlinks/?variant=NM_000088.3%3Ac.589G%3ET&alignment=splign&alt_acc=NG_007400.1:g.8638G>T&request=reference
	input = request.args.get('input', '')
	alt_aln_method = request.args.get('alt_aln_method', '')
	accession = request.args.get('accession', '')
	
	# Tell the loading page which loader is being called
	loader = 'add_loading'
	# Sent the data to an auto submit form loading 
	return render_template('bootstrap/loading.html', title=title, user=user, loader=loader, input=input, alt_aln_method=alt_aln_method, accession=accession)


# Route for adding to database
@app.route('/data_add/', methods=['POST', 'GET'])
def data_add():
	input = request.args.get('variant', '')
	alt_aln_method = request.args.get('alignment', '')
	accession = request.args.get('accession', '')
	title = "db_adding"

	# return render_template('bootstrap/variantError.html', title=title, user=input, error=alt_aln_method, reason=accession)
	
	# Add accurate transcript descriptions to the database
	# RefSeq databases
	if alt_aln_method != 'genebuild':		
		# Get the Entrez (GenBank) file
		# return render_template('bootstrap/variantError.html', title=title, user=input, error=alt_aln_method, reason=accession)
		record = variantanalyser.functions.entrez_efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
		# Extract the information required
		desc = ''
		desc = record.description
		record = ''
		#return render_template('bootstrap/variantError.html', title=input, user=desc, error=accession)
		# Add the record to the database
		# Connect to database and send request
		tm = time.ctime()
		# return render_template('bootstrap/variantError.html', title=title, user=tm, error=alt_aln_method, reason=accession)
		data_added = 'false'
		data_added = variantanalyser.data.add_entry(acn=accession, desc=desc, connected=g.db, tm=tm)
		if data_added == 'true':
			return redirect(url_for('analysis', variant=input, alignment=alt_aln_method))
		else:
			return render_template('bootstrap/variantError.html', title=error, user='Data entry failed', error='Please contact admin if the problem persists', reason='Please re-submit')
			
	# Ensembl databases
	else:
		#return render_template('bootstrap/variantError.html', title=title, user=input, error=accession)
		# THIS STAGE IS CONTROLLED BY CASCADING IF STATEMENTS
		# Use EnsemblRest to search the Ensembl database
		decoded = variantanalyser.functions.ensembl_rest(ext = "/overlap/id/" + accession + "?feature=transcript", primary_assembly=primary_assembly)
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
			data = variantanalyser.functions.hgnc_rest(path = "/search/ensembl_gene_id/" + ensg)
			if data['error'] != 'false':
				hgnc_gene_info = 'Cannot currently display gene information: ' + decoded['error']
				data = ''
			else:
				# Set the hgnc name correctly
				hgnc = data['record']['response']['docs'][0]['symbol']
				# Now re-search the HGNC database with the correct hgnc name to get the correct description
				data = variantanalyser.functions.hgnc_rest(path = "/fetch/symbol/" + hgnc)
				if data['error'] != 'false':
					hgnc_gene_info = 'Cannot currently display gene information: ' + decoded['error']
					data = ''
				else:
					description = data['record']['response']['docs'][0]['name']
					# Check for multiple transcripts to edit the description
					try:
						var_g = variantanalyser.functions.genomic(variant, evm, hp)
					except hgvs.exceptions.HGVSError as e:
						error = e
					if error == 'false':
						var_g = variantanalyser.functions.genomic(variant, evm, hp)
						multiple = variantanalyser.functions.relevant_transcripts(var_g, evm)
						if len(multiple) > 1:
							# Put together a refseq style description
							desc = 'Homo sapiens ' + description + " (" + hgnc + "), transcript variant " + tx_variant + ", mRNA." 
						else:
							desc = 'Homo sapiens ' + description + " (" + hgnc + "), mRNA." 
						# Set the description to hgnc_gene_info
						data = ''
						tm = time.ctime()
						data_added = 'false'
						data_added = variantanalyser.data.add_entry(acn=accession, desc=desc, connected=g.db, tm=tm)
						if data_added == 'true':
							return redirect(url_for('analysis', variant=input, alignment=alt_aln_method))
						else:
							return render_template('bootstrap/variantError.html', title=error, user='Data entry failed', error='Please contact admin if the problem persists', reason='Please re-submit')
					else:
						# Open a hgvs exception log file in append mode
						fo = open('/local/logs/hgvsExceptions.txt', 'a')
						excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
						fo.write(excep)
						fo.close()
						data = ''
						return render_template('bootstrap/variantError.html', title=title, user=input, error=error, reason = 'Invalid Variant description')


# Route for contact admin
@app.route('/contact_admin/')
def contact_admin(form=None):
	if form is None:
		form = ContactForm()
	#name = session.get("name", [])
	#email = session.get("email", [])
	#variant = session.get("variant", [])
	#query = session.get("query", [])

	title = 'ContactAdmin'
	sub = 'Please fill in all mandatory fields'
	valid = 'true'
	# return render_template('bootstrap/contact_admin.html', title=title, sub=sub, valid=valid, name=name, email=email, variant=variant, query=query, form=form)
	return render_template('bootstrap/contact_admin.html', title=title, sub=sub, valid=valid, form=form)


# Route email admin response
@app.route('/email_admin/', methods=['POST'])
def email_admin():

	if request.method == 'POST':
		# Collect POST data from form
		# Collect POST data from form
		form = ContactForm()
		if form.validate_on_submit():
			name = form.name.data
			email = form.email.data
			variant = form.variant.data
			query = form.query.data
			honey = 'false'
			honey = form.honey.data
			# return render_template('bootstrap/variantError.html', title='title', user=email, error=name, reason = query)
		
		# If the honeypot is filled in, set up a spam query string
		if honey != 'false':
			query = ''
			query = 'penis viagra bingo'
			
		# if request.form["name"] != '' and request.form["email"] != '' and request.form["query"] != '':
			#name = request.form["name"]
			#email = request.form["email"]
			#variant = request.form["variant"]
			#query = request.form["query"]
	
			# Log queries
			fo = open('/local/logs/email_log.txt', 'a')
			excep = "%s \t %s \t %s \t %s \t %s\n" %(time.ctime(), name, email, query, variant)
			fo.write(excep)
			fo.close()

			# Some basic spam filtering
			at = re.compile('@')
			dot = re.compile('\.')
			# Is the at present?
			if at.search(email):
				elements = email.split('@')
				if len(elements) > 2:
					title = 'bad_request'
					return render_template('bootstrap/bad_request.html', title=title, user=email, error='Blocked email address', reason = 'Email address is not valid')
				else:
					domain = ''.join(email.split('@')[1])
					# alias = ''.join(email.split('@')[0])
				
				# Is the dot present?
				if dot.search(domain):
					# List of bad domain names stored in text file
					# open the file as a list
					bad_domains = []
					with open("/local/logs/bad_domain.txt", 'r') as fo:
						bad_domains = fo.read().splitlines()
 					fo.close()

					# Check for bad domain and remove
					for d in bad_domains:
						if d != domain:
							pass
						else:
							title = 'bad_request'
							return render_template('bootstrap/bad_request.html', title=title, user=email, error='Blocked email address', reason = 'Email address is not valid')

					# Import spam filter
					# Import spam filter and configure
					from reverend.thomas import Bayes
					import spam_filter.training
					import spam_filter.spam_filter

					spam_db='/local/logs/filter.bayes'
					guesser = Bayes()

					# Load the spam filter database
					spam_filter.training.load(spam_db, guesser)

					
					
					
					
					
					
					
					# Train spam filter
					####################
					spammy = 'congratulations sexy viagra penis won ppi compensation bingo spam vanquis'
					hammy = 'variant description genomic coding'
					spam_filter.training.train_spam(spammy, guesser, spam_db)
					spam_filter.training.train_ham(hammy, guesser, spam_db)	
					
					
					
					
					
					
					
					
								
						
					# Is the query content bad?
					spam_tot = 0
					ham_tot = 0
					total = 0
					cln = ''
					text = query.split()
					for word in text:
						result = spam_filter.spam_filter.guess(word, guesser)
						if result[0] > 0:
							ham_tot = ham_tot + 1
						if result[1] > 0:
							spam_tot = spam_tot + 1
						total = total + 1
					
					# Set the sum to use decimals
					getcontext().prec = 25
					hams = Decimal(ham_tot) / Decimal(total)
					spams = Decimal(spam_tot) / Decimal(total)
					# return render_template('bootstrap/bad_request.html', title='giraffe', user='f', error=spams, reason = hams)
						
					# figure out if it's spam or not!	
					cln = 'good'
					if (spams) > (hams):
						cln = 'bad'
						
					if cln == 'good': 
						# Compile the message
						message = variant + '\n' + query + '\n' + name
						subject = 'hgvsWeb query' + ' ' + '(' + variant + ')'
	
						# Format into msg
						msg = Message(recipients=["pjf9@le.ac.uk"],	
		            					sender=email,
		            					body=message,
        	        					subject=subject)
                	
						# Send the email
						mail.send(msg)
	
						# Load variantanalyser with a message
						title = 'Email sent to admin'
						primary_assembly = 'GRCh37'
						caution = 'Message sent to admin'
						return render_template('bootstrap/variantAnalyser.html', title=title, caution=caution, primary_assembly=primary_assembly, message=message)	
					else:
						# Check against a list of good domains
						# Open good domains
						with open("/local/logs/good_domain.txt", 'r') as fo:
							good_domains = fo.read().splitlines()
 						fo.close()
						for d in good_domains:
							if d == domain:
								mark_bad = 'false'
								break
							else:
								mark_bad = 'true'
								
						if mark_bad == 'true':
							# Add the spam filtered domain name to the bad_domain list
							fo = open('/local/logs/bad_domain.txt', 'a')
							excep = "%s \n" %(domain)
							fo.write(excep)
							fo.close()
							# Alert admin that new bad domain is being added
							message = 'A new domain name ' + domain + ' has beed added to bad_domain.txt.\n\nThe domain name originated from ' + email + ' \n\n' + query + '. \n\n'
							subject = 'VariantValidator.org: bad domain name added' + ' ' + '(' + domain + ')'
	
							# Format into msg
							msg = Message(recipients=["pjf9@le.ac.uk"],	
		            						sender=email,
		            						body=message,
        	    							subject=subject)
                	
							# Send the email
							mail.send(msg)
							
							# Report that it has been filtered
							title = 'bad_request'
							return render_template('bootstrap/bad_request.html', title=title, user=email, error='Email submission denied', reason = 'Email is considered to be Spam')
						else:
							# Send the spam mail to check for errors
							# Alert admin that new bad domain is being added
							message = 'A suspected spam email has been sent from domain ' + '\n\nThe domain name originated from ' + email + '\n\n' + query + '. \n\n'
							subject = 'VariantValidator.org: Suslicious email recieved' + ' ' + '(' + email + ')'
	
							# Format into msg
							msg = Message(recipients=["pjf9@le.ac.uk"],	
		            						sender=email,
		            						body=message,
        	    							subject=subject)
                	
							# Send the email
							mail.send(msg)
							
							# Report that it has been filtered
							title = 'bad_request'
							return render_template('bootstrap/bad_request.html', title=title, user=email, error='Email submission denied', reason = 'Email is considered to be Spam')
						
				# bad email address
				else:	
					title = 'bad_request'
					return render_template('bootstrap/bad_request.html', title=title, user=email, error='Email submission denied', reason = 'Email address is not valid')
			# bad email address
			else:	
				title = 'bad_request'
				return render_template('bootstrap/bad_request.html', title=title, user=email, error='Email submission denied', reason = 'Email address is not valid')
		# incomplete mandatory fields
		else:
			form = ContactForm()
			title = 'ContactAdmin'
			sub = 'Please fill in ALL mandatory fields'
			valid = 'false'
			return render_template('bootstrap/contact_admin.html', title=title, sub=sub, valid=valid, form=form)
			# return render_template('bootstrap/variantError.html', title=title, user='poo', error='giraffe', reason = 'Invalid Variant description')


# Route for fasta download page
@app.route('/fasta_download/', methods=['POST'])
def fasta_download():
	# Collect the form data
	header = request.form["header"]
	sequence = request.form["sequence"]
	
	sequence = sequence.replace('[', '')
	sequence = sequence.replace(']', '')
	sequence = sequence.replace("'", "")
	sequence = sequence.replace('\\n', '')
	sequence = sequence.replace(',', '')
	sequence = sequence.replace(' ', '\n') 
	
	# String file together
	file = header.strip() + '\n' + sequence
	
	# Make response
	response = make_response(file)
	
	# Set the download header rather than just print on browser
	response.headers["Content-Disposition"] = "attachment; filename=hgvsWeb_download.txt"
	
	return response
	
# Route for alignment download page
@app.route('/alnignment_download/', methods=['POST'])
def alignment_download():
	# Collect the form data
	alignment=request.form["alignment"]
	head=request.form["head"]
	
	# Format the head and remove unwanted characters
	hd = str(head)
	hd = hd.replace("[u'", "")
	
	# Format the alignment and remove unwanted characters
	aln=str(alignment)

	aln = aln.replace(']', '')
	aln = aln.replace('[', '')
	aln = aln.replace("'", "")
	aln = aln.replace(', ', '')
	aln = aln.replace("\\n",'\n')
	
	# Assemble the file
	file = hd + aln
	
	
	# Make response
	response = make_response(file)
	
	# Set the download header rather than just print on browser
	response.headers["Content-Disposition"] = "attachment; filename=hgvsWeb_download.txt"
	
	return response



# Run the server (if file is called directly by python, server is internal dev server)
if __name__ == '__main__':
    app.debug=True
    app.run()
	
	
	
	
	
	
	
	
	
	