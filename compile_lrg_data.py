import urllib2
import re
# import untangle
import variantanalyser
import variantanalyser.dbControls
import variantanalyser.dbControls.data as data
import time

def update():
	print 'Updating LRG lookup tables'
	lr2rs_download = urllib2.Request('http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_transcripts_xrefs.txt')
	# Open and read
	lr2rs_data = urllib2.urlopen(lr2rs_download)
	lr2rs = lr2rs_data.read()
	# List the data
	lr2rs = lr2rs.strip()
	lr2rs = lr2rs.split('\n')

	# Download 
	lrg_status_download = urllib2.Request('http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_GRCh38.txt')
	# Open and read
	lrg_status_data = urllib2.urlopen(lrg_status_download)
	lrg_status = lrg_status_data.read()
	# List the data
	lrg_status = lrg_status.strip()
	lrg_status = lrg_status.split('\n')

	# Download 
	rs2lr_download = urllib2.Request('http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene')
	# Open and read
	rs2lr_data = urllib2.urlopen(rs2lr_download)
	rs2lr = rs2lr_data.read()
	# List the data
	rs2lr = rs2lr.strip()
	rs2lr = rs2lr.split('\n')

	# Download LRG transcript (_t) to LRG Protein (__p) data file
	lr_t2p_downloaded = urllib2.Request('http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_proteins_RefSeq.txt')
	# Open and read
	lr_t2p_data = urllib2.urlopen(lr_t2p_downloaded)
	lr_t2p = lr_t2p_data.read()
	# List the data
	lr_t2p = lr_t2p.strip()
	lr_t2p = lr_t2p.split('\n')

	# Dictionary the status by LRG_ID
	lrg_status_dict = {}
	# Compile dictionary
	for line in lrg_status:
		if re.search('^#', line):
			continue
		else:
			list = line.split()
			lrgID = list[0]
			stat = list[2]
			lrg_status_dict[lrgID] = stat

	# Required lookup tables
	# LRG_ID	GeneSymbol	RefSeqGeneID	status	
	# LRG_ID	RefSeqTranscriptID
	# LRG_T2LRG_P

	print 'Update LRG and LRG_transcript lookup tables'
	# Populate lists lrg_rs_lookup (LRG to RefSeqGene) and lrg_t2nm_ (LRG Transcript to RefSeq Transcript)
	for line in lr2rs:
		if re.search('^#', line):
			continue
		else:
			list = line.split()
			# Assign objects
			lrg_id = list[0]
			symbol = list[1]
			rsgid = list[2]
			lrg_tx = str(list[0]) + str(list[3])
			rstid = list[4]
			status = lrg_status_dict[lrg_id]
			# pass data to relevant lists
			# lrg_rs_lookup
			lrg_rs_lookup = [lrg_id, symbol, rsgid, status]

			# update LRG to RefSeqGene database
			data.update_lrg_rs_lookup(lrg_rs_lookup)
		
			# lrg_t2nm_
			lrgtx_to_rstID = [lrg_tx, rstid]
			# update database
			data.update_lrgt_rst(lrgtx_to_rstID)

	print 'Update LRG protein lookup table'
	# Populate LRG protein RefSeqProtein lokup table
	for line in lr_t2p:
		if re.search('^#', line):
			continue
		else:
			list = line.split()
			# Assign objects
			lrg_p = list[0]
			rs_p = list[1]
			# update LRG to RefSeqGene database
			data.update_lrg_p_rs_p_lookup(lrg_p, rs_p)
	
	print 'LRG lookup tables updated'
	return			
		

	




		
		
		
		
		
		
		
		
		
			