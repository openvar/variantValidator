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
    entry = dbControls.data.in_entries(accession, 'transcript_info')
    return entry

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
