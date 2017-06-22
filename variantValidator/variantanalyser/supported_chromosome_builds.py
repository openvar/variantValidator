def supported_for_mapping(ac, primary_assembly):	
	GRCh37 = {
			"NC_000001.10" : 'true',
			"NC_000002.11" : 'true',
			"NC_000003.11" : 'true',
			"NC_000004.11" : 'true',
			"NC_000005.9" : 'true',
			"NC_000006.11" : 'true',
			"NC_000007.13" : 'true',
			"NC_000008.10" : 'true',
			"NC_000009.11" : 'true',
			"NC_000010.10" : 'true',
			"NC_000011.9" : 'true',
			"NC_000012.11" : 'true',
			"NC_000013.10" : 'true',
			"NC_000014.8" : 'true',
			"NC_000015.9" : 'true',
			"NC_000016.9" : 'true',
			"NC_000017.10" : 'true',
			"NC_000018.9" : 'true',
			"NC_000019.9" : 'true',
			"NC_000020.10" : 'true',
			"NC_000021.8" : 'true',
			"NC_000022.10" : 'true',
			"NC_000023.10" : 'true',
			"NC_000024.9" : 'true'	
				}
	GRCh38 = {
			"NC_000001.11" : 'true',
			"NC_000002.12" : 'true',
			"NC_000003.12" : 'true',
			"NC_000004.12" : 'true',
			"NC_000005.10" : 'true',
			"NC_000006.12" : 'true',
			"NC_000007.14" : 'true',
			"NC_000008.11" : 'true',
			"NC_000009.12" : 'true',
			"NC_000010.11" : 'true',
			"NC_000011.10" : 'true',
			"NC_000012.12" : 'true',
			"NC_000013.11" : 'true',
			"NC_000014.9" : 'true',
			"NC_000015.10" : 'true',
			"NC_000016.10" : 'true',
			"NC_000017.11" : 'true',
			"NC_000018.10" : 'true',
			"NC_000019.10" : 'true',
			"NC_000020.11" : 'true',
			"NC_000021.9" : 'true',
			"NC_000022.11" : 'true',
			"NC_000023.11" : 'true',
			"NC_000024.10" : 'true'	
				}
	
	if primary_assembly == 'GRCh37':
		sfm = GRCh37.get(ac)
	elif primary_assembly == 'GRCh38':
		sfm = GRCh38.get(ac)
	else:
		sfm = 'false'
	
	
	if sfm == 'true':
		pass
	else:
		sfm = 'false'
	
	return sfm
	
	
def to_accession(chr_num, primary_assembly):
	# Available genome builds
	GRCh37 = {
			"1" : "NC_000001.10",
			"2" : "NC_000002.11",
			"3" : "NC_000003.11",
			"4" : "NC_000004.11",
			"5" : "NC_000005.9",
			"6" : "NC_000006.11",
			"7" : "NC_000007.13",
			"8" : "NC_000008.10",
			"9" : "NC_000009.11",
			"10" : "NC_000010.10",
			"11" : "NC_000011.9",
			"12" : "NC_000012.11",
			"13" : "NC_000013.10",
			"14" : "NC_000014.8",
			"15" : "NC_000015.9",
			"16" : "NC_000016.9",
			"17" : "NC_000017.10",
			"18" : "NC_000018.9",
			"19" : "NC_000019.9",
			"20" : "NC_000020.10",
			"21" : "NC_000021.8",
			"22" : "NC_000022.10",
			"23" : "NC_000023.10",
			"24" : "NC_000024.9",	
			"x" : "NC_000023.10",
			"y" : "NC_000024.9",
			"X" : "NC_000023.10",
			"Y" : "NC_000024.9"
			}

	GRCh38 = {
			"1" : "NC_000001.11",
			"2" : "NC_000002.12",
			"3" : "NC_000003.12",
			"4" : "NC_000004.12",
			"5" : "NC_000005.10",
			"6" : "NC_000006.12",
			"7" : "NC_000007.14",
			"8" : "NC_000008.11",
			"9" : "NC_000009.12",
			"10" : "NC_000010.11",
			"11" : "NC_000011.10",
			"12" : "NC_000012.12",
			"13" : "NC_000013.11",
			"14" : "NC_000014.9",
			"15" : "NC_000015.10",
			"16" : "NC_000016.10",
			"17" : "NC_000017.11",
			"18" : "NC_000018.10",
			"19" : "NC_000019.10",
			"20" : "NC_000020.11",
			"21" : "NC_000021.9",
			"22" : "NC_000022.11",
			"23" : "NC_000023.11",
			"24" : "NC_000024.10",
			"x" : "NC_000023.11",
			"y" : "NC_000024.10",
			"X" : "NC_000023.11",
			"Y" : "NC_000024.10"				
			}
	
	# Convert call line to rs line
	if primary_assembly == 'GRCh37':
		chr_accession = GRCh37.get(chr_num)
	if primary_assembly == 'GRCh38':
		chr_accession = GRCh38.get(chr_num)
	#if primary_assembly == 'NCBI36':
	#	chr_accession] = NCBI36.get(chr_num)
	# Return rs 
	return chr_accession
	
def to_chr_num(accession):
	# Available genome builds
	chr_num_convert = {
			"NC_000001" : "1",
		    "NC_000002" : "2",
			"NC_000003" : "3",
			"NC_000004" : "4",
			"NC_000005" : "5",
			"NC_000006" : "6",
			"NC_000007" : "7",
			"NC_000008" : "8",
			"NC_000009" : "9",
			"NC_000010" : "10",
			"NC_000011" : "11",
			"NC_000012" : "12",
			"NC_000013" : "13",
			"NC_000014" : "14",
			"NC_000015" : "15",
			"NC_000016" : "16",
			"NC_000017" : "17",
			"NC_000018" : "18",
			"NC_000019" : "19",
			"NC_000020" : "20",
			"NC_000021" : "21",
			"NC_000022" : "22",
			"NC_000023" : "X",
			"NC_000024" : "Y"
			}
	accession = accession.split('.')[0]
	chr_num = chr_num_convert.get(accession)
	return chr_num