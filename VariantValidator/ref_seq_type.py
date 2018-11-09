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
# Copyright (C) 2018  Peter Causey-Freeman, University of Leicester
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
