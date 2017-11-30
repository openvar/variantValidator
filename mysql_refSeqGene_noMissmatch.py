import re
import os
import urllib2
import copy
import variantanalyser
import variantanalyser.dbControls
import variantanalyser.dbControls.data as db_data

def update():
	print 'Updating RefSeqGene no Missmatch MySQL data'
	# Set os path
	# Set up os paths data and log folders 
	ROOT = os.path.dirname(os.path.abspath(__file__))

	# Download data from RefSeqGene
	# Download data
	rsg = urllib2.Request('http://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene')
	response = urllib2.urlopen(rsg)
	rsg_file = response.read()
	rsg_data_line = rsg_file.split('\n')
	rsg_data = []
	for data in rsg_data_line:
		rsg_data.append(data)

	# Download data
	grch37 = urllib2.Request('http://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/GCF_000001405.25_refseqgene_alignments.gff3')
	response = urllib2.urlopen(grch37)
	grch37_file = response.read()
	grch37_data_line = grch37_file.split('\n')
	grch37_align_data = []
	for data in grch37_data_line:
		grch37_align_data.append(data)

	# Download data
	grch38 = urllib2.Request('http://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/GCF_000001405.28_refseqgene_alignments.gff3')
	response = urllib2.urlopen(grch38)
	grch38_file = response.read()
	grch38_data_line = grch38_file.split('\n')
	grch38_align_data = []
	for data in grch38_data_line:
		grch38_align_data.append(data)

	# Open Lists
	#rsg_data = open(os.path.join(ROOT, 'gene_RefSeqGene'), 'r')
	rsg_id_info = []
	# grch37_align_data = open(os.path.join(ROOT, 'GCF_000001405.25_refseqgene_alignments.gff3'), 'r')
	grch37_align = []
	# grch38_align_data = open(os.path.join(ROOT, 'GCF_000001405.28_refseqgene_alignments.gff3'), 'r')
	grch38_align = []

	# Place the required data from each file into a dictionary
	hash = re.compile('#')
	for line in rsg_data:
		if hash.search(line):
			pass
		else:
			line = line.strip()
			info = line.split()
			if len(info) == 0:
				pass
			else:
				dict = {'symbol' : info[2], 'rsg_id' : info[3], 'gene_id' : info[1]} 
				rsg_id_info.append(dict)

	# Create dictionary to store RefSeqGene and gene symbol data NOTE RefSeqGene ID stored without version number!
	rsg_to_symbol = {}
	# Collect the data
	for ent in rsg_id_info:
		rsg_id = copy.deepcopy(ent['rsg_id'])
		rsg_id = rsg_id.split('.')[0]
		rsg_to_symbol[rsg_id] = {'symbol' : ent['symbol'], 'gene_id' : ent['gene_id']}

	# Count total number of NG to NC mappings
	total_rsg_to_nc = 0
	total_rsg_to_nc_rejected = 0
	for line in grch37_align_data:
		# Count NG_ to NC_ and remove the entries we don't care about!
		if re.search('NC_', line) and re.search('NG_', line):
			total_rsg_to_nc = total_rsg_to_nc + 1
		else:
			continue
		if hash.search(line):
			pass
		elif not re.search('gap_count=0', line):
			if re.search('NC_', line) and re.search('NG_', line):
				total_rsg_to_nc_rejected = total_rsg_to_nc_rejected + 1
			# print line
			pass
		else:
			line = line.strip()
			info = line.split('\t')
			if len(info) != 9:
				pass
			else:
				metrics = info[8].split(';')
				id_ori = metrics[1].replace('Target=', '')
				id_ori_list = id_ori.split()
				dict = {'rsg_id' : id_ori_list[0], 'chr_id' : info[0], 'rsg_start' : info[3], 'rsg_end' : info[4], 'ori' : id_ori_list[3]} 
				grch37_align.append(dict)

	for line in grch38_align_data:
		if re.search('NC_', line) and re.search('NG_', line):
			total_rsg_to_nc = total_rsg_to_nc + 1
		else:
			continue
		if hash.search(line):
			pass
		elif not re.search('gap_count=0', line):
			if re.search('NC_', line) and re.search('NG_', line):
				total_rsg_to_nc_rejected = total_rsg_to_nc_rejected + 1
			# print line
			pass
		else:
			line = line.strip()
			info = line.split('\t')
			if len(info) != 9:
				pass
			else:		
				metrics = info[8].split(';')
				id_ori = metrics[1].replace('Target=', '')
				id_ori_list = id_ori.split()
				dict = {'rsg_id' : id_ori_list[0], 'chr_id' : info[0], 'rsg_start' : info[3], 'rsg_end' : info[4], 'ori' : id_ori_list[3]} 
				grch38_align.append(dict)

	# Create a data array containing the database
	db = []
	# map line
	for line in grch37_align:
		ml = []
		link = line['rsg_id']
		ml.append(link)
		ml.append(line['chr_id'])
		ml.append('GRCh37')
		ml.append(line['rsg_start'])
		ml.append(line['rsg_end'])
		ml.append(line['ori'])
		# Add the additional data from rsg_id_info
		for data in rsg_id_info:
			if link == data['rsg_id']:
				ml.append(data['symbol'])
				ml.append(data['gene_id'])
			else:
				continue
		# Create the entry and append to db
		db.append(ml)

	for line in grch38_align:
		ml = []
		link = line['rsg_id']
		ml.append(link)
		ml.append(line['chr_id'])
		ml.append('GRCh38')
		ml.append(line['rsg_start'])
		ml.append(line['rsg_end'])
		ml.append(line['ori'])
		# Add the additional data from rsg_id_info
		for data in rsg_id_info:
			if link == data['rsg_id']:
				ml.append(data['symbol'])
				ml.append(data['gene_id'])
			else:
				continue
		# Create the entry and append to db
		db.append(ml)

	# Known missing identifiers
	known = {
			'NG_021289.1' : {'symbol' : 'CFAP47', 'gene_id' : '286464'},
			'NG_027707.1' : {'symbol' : 'DUX4L1', 'gene_id' : '22947'},
			'NG_033266.1' : {'symbol' : 'DSE', 'gene_id': '29940'}
			}
	
	# Known Obsolete identifiers
	obsolete = {
				'NG_016553.1' : 'OBSOLETE'
				}
		
	# Identify lines with missing data e.g. gene symbols
	for line in db:
		try:
			line[6]
		except IndexError:
			try:
				identifier = copy.deepcopy(line[0])
				identifier = identifier.split('.')[0]
				line.append(rsg_to_symbol[identifier]['symbol'])
				line.append(rsg_to_symbol[identifier]['gene_id'])
			except KeyError:
				try:
					line.append(known[line[0]]['symbol'])
					line.append(known[line[0]]['gene_id'])
				except KeyError:
					check = obsolete[line[0]]
					print str(line[0]) + ' : ' + check	


	# Open a text file to be used as a simple database and write the database
	# rsg_db = open(os.path.join(ROOT, 'rsg_chr_db.txt'), 'w')

	to_mysql = []
	for line in db:
		if line[0] in obsolete.keys():
			continue
		# Only gap-less RefSeqGenes will have passed. The rest will be alternatively curated
		write = []
		# Take the mapping data
		write = copy.deepcopy(line[0:6])
		# add RSG ranges
		write.append('1')
		end_rsg = int(line[4]) - int(line[3]) + 1
		end_rsg = str(end_rsg)
		write.append(end_rsg)
		# Create block data chr then rsg
		chr_block = str(line[3]) + '-' + str(line[4])
		write.append(chr_block)
		rsg_block = str(write[6]) + '-' + str(write[7])
		write.append(rsg_block)
		# Add gene ID and Gene symbol(s)
		write.append(line[7])
		write.append(line[6])
		# write_me = '\t'.join(write)
		# rsg_db.write(write_me + '\n')
		del write[6]
		to_mysql.append(write)

	# Set up code to write to database
	for line in to_mysql:
		current_symbol = db_data.get_gene_symbol_from_refSeqGeneID(line[0])
		if line[10] == current_symbol:
			pass
		else:
			if current_symbol != 'none':
				line[10] = current_symbol
			else:
				pass		
		db_data.update_refSeqGene_loci(line)		


	# Close database
	# rsg_db.close()

	print 'Total NG_ to NC_ alignments = ' + str(total_rsg_to_nc)
	print 'Gapps within NG_ to NC_ alignments = ' + str(total_rsg_to_nc_rejected)

	print 'complete'
	return







