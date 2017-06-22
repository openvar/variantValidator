# Defining reference sequence type from accession
def ref_type_assign(accession):
	import re
	if re.match('NC_', accession) or re.match('NG_', accession) or re.match('NT_', accession) or re.match('NW_', accession):
		ref_type = ':g.'		
	elif re.match('NM_', accession):
		ref_type = ':c.'
	elif re.match('NR_', accession):
		ref_type = ':n.'
	elif re.match('NP_', accession):
		ref_type = ':p.'	
	return ref_type