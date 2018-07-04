# -*- coding: utf-8 -*-
"""
ref_seq_type.py

Simple function that assignes the correct reference sequence type (c., n., p., g.) to
reference sequences

# Example
ref_type_assign(accession)
"""

# Defining reference sequence type from accession
import re
from variantanalyser import dbControls


def ref_type_assign(accession):
    if re.match('NC_', accession) or re.match('NG_', accession) or re.match('NT_', accession) or re.match('NW_',
                                                                                                          accession):
        ref_type = ':g.'
    elif re.match('NM_', accession):
        ref_type = ':c.'
    elif re.match('NR_', accession):
        ref_type = ':n.'
    elif re.match('NP_', accession):
        ref_type = ':p.'
    elif re.match('LRG_', accession):
        if re.search('t', accession):
            refseqtranscript_reference = dbControls.data.get_RefSeqTranscriptID_from_lrgTranscriptID(accession)
            if re.match('NM_', refseqtranscript_reference):
                ref_type = ':c.'
            else:
                ref_type = ':n.'
        elif re.search('_p', accession):
            ref_type = ':p.'
        else:
            ref_type = ':g.'
    return ref_type

# <LICENSE>

# </LICENSE>
