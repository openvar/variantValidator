# -*- coding: utf-8 -*-

"""
data.py

Contains all database functions

Takes requests via the functions and makes the appropriate MySQL queries via the relevant
query type
"""

# import database modules
import dbquery
import dbinsert
import dbupdate
import dbfetchone
import dbfetchall

# Import python modules
import os
import sys
import re

# Needs functions from variantanalyser - directory above, unless in a single directory
try:
    import variantanalyser.functions as functions
except ImportError:	
    parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.sys.path.insert(0,parentdir)
    import functions
except AttributeError:	
    parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.sys.path.insert(0,parentdir)
    import functions

# Retrieve transcript information
def in_entries(entry, table): 

    # Use dbquery.py to connect to mysql and return the necessary data
    if table == 'transcript_info':
        row = dbquery.query_with_fetchone(entry, table)

        if row[0] == 'error':
            data = {
                'error' : 'false',
                'description': 'false'
                }

            data['error'] = row[0]
            data['description'] = row[1]

        elif row[0] == 'none':
            data = {
                'none' : 'false',
                'description': 'false'
                }

            data['none'] = row[0]
            data['description'] = row[1]

        else:
            data = {
                'accession' : 'false',
                'description': 'false',
                'updated' : 'false',
                'expiry' : 'false'
                }

            data['accession'] = row[0]
            data['description'] = row[1]
            data['variant'] = row[2]
            data['version'] = row[3]
            data['hgnc_symbol'] = row[4]
            data['uta_symbol'] = row[5]
            data['updated'] = row[6]
            data['expiry'] = row[7]

    return data

# Add new entry  	
def add_entry(entry, data, table):
    success = dbinsert.insert(entry, data, table)
    return success
    
def insert_transcript_loci(add_data, primary_assembly):
    success = dbinsert.insert_transcript_loci(add_data, primary_assembly)
    return success


# Update entries
def update_entry(entry, data, table):
    success = dbupdate.update(entry, data, table)
    return success

def update_transcript_info_record(accession, hdp):

    # Search Entrez for corresponding record for the RefSeq ID
    try:
        record = functions.entrez_efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        version = record.id
        description = record.description
        variant = '0'

        if re.search('transcript variant', description):
            tv = re.search('transcript variant \w+', description)
            tv = str(tv.group(0))
            tv = tv.replace('transcript variant', '')
            variant = tv.strip()
            variant = variant.upper() # Some tv descriptions are a or A
        else:
            variant = '0'

        # Get information from UTA
        try:
            uta_info = hdp.get_tx_identity_info(version)
        except:
            version_ac_ver = version.split('.')
            version = version_ac_ver[0] + '.' + str(int(version_ac_ver[1]) - 1)
            uta_info = hdp.get_tx_identity_info(version)

        uta_symbol = str(uta_info[6])

        # First perform a search against the input gene symbol or the symbol inferred from UTA
        initial = functions.hgnc_rest(path = "/fetch/symbol/" + uta_symbol)
        # Check for a record
        if str(initial['record']['response']['numFound']) != '0':
            hgnc_symbol = uta_symbol
        # No record found, is it a previous symbol?
        else:
            # Search hgnc rest to see if symbol is out of date
            rest_data = functions.hgnc_rest(path = "/search/prev_symbol/" + uta_symbol)
            # If the name is correct no record will be found
            if rest_data['error'] == 'false':
                if int(rest_data['record']['response']['numFound']) == 0:
                    hgnc_symbol = uta_info[6]
                else:
                    hgnc_symbol = rest_data['record']['response']['docs'][0]['symbol']
            else:
                hgnc_symbol = 'unassigned'

    # List of connection error types. May need to be expanded.
    # Outcome - Put off update for 3 months!
    except Exception as e:
        if str(e) == '<urlopen error [Errno -2] Name or service not known>':
            # Issues with DNSSEC for the nih.gov
            previous_entry = in_entries(accession, 'transcript_info')
            accession = accession
            description = previous_entry['description']
            variant = previous_entry['variant']
            version = previous_entry['version']
            hgnc_symbol = previous_entry['hgnc_symbol']
            uta_symbol = previous_entry['uta_symbol']

    # Query information
    # query_info = [accession, description, variant, version, hgnc_symbol, uta_symbol]
    query_info = [version, description, variant, version, hgnc_symbol, uta_symbol]
    table='transcript_info'

    # Update the transcript_info table (needs plugging in)
    returned_data = in_entries(version, table)
    # If the entry is not in the database add it
    if 'none' in returned_data:
        add_entry(version, query_info, table)
    # If the data in the entry has changed, update it
    else:
        update_entry(version, query_info, table)
    return

def update_refSeqGene_loci(rsg_data):
    # First query the database
    # import dbfetchone
    entry_exists = dbfetchone.get_refSeqGene_data_by_refSeqGeneID(rsg_data[0], rsg_data[2])
    if entry_exists[0] == 'none':
        # import dbinsert
        dbinsert.insert_refSeqGene_data(rsg_data)
    else:
        # import dbupdate
        dbupdate.update_refSeqGene_data(rsg_data)
    return

def update_lrg_rs_lookup(lrg_rs_lookup):
    # First query the database
    rsgID = dbfetchone.get_RefSeqGeneID_from_lrgID(lrg_rs_lookup[0])
    if rsgID == 'none':
        # import dbinsert
        dbinsert.insert_RefSeqGeneID_from_lrgID(lrg_rs_lookup)
        return
    else:
        return

def update_lrgt_rst(lrgtx_to_rstID):
    # First query the database
    rstID = dbfetchone.get_RefSeqTranscriptID_from_lrgTranscriptID(lrgtx_to_rstID[0])
    if rstID == 'none':
        # import dbinsert
        dbinsert.insert_LRG_transcript_data(lrgtx_to_rstID)
        return
    else:
        return

def update_lrg_p_rs_p_lookup(lrg_p, rs_p):
    # First query the database
    rspID = dbfetchone.get_RefSeqProteinID_from_lrgProteinID(lrg_p)
    if rspID == 'none':
        # import dbinsert
        dbinsert.insert_LRG_protein_data(lrg_p, rs_p)
        return
    else:
        return

# Direct methods (GET)
def get_transcript_info_for_gene(gene_symbol):
    rows = dbfetchall.get_transcript_info_for_gene(gene_symbol)
    return rows

def get_uta_symbol(gene_symbol):
    # returns the UTA gene symbol when HGNC gene symbol is input
    utaSymbol = str(dbfetchone.get_utaSymbol(gene_symbol)[0])
    return utaSymbol

def get_hgnc_symbol(gene_symbol):
    # returns the HGNC gene symbol when UTA gene symbol is input
    hgncSymbol = str(dbfetchone.get_hgncSymbol(gene_symbol)[0])
    return hgncSymbol

def get_transcript_description(transcript_id):
    # returns the transcript description for a given transcript
    tx_description = dbfetchone.get_transcript_description(transcript_id)
    return tx_description

def get_gene_symbol_from_transcriptID(transcript_id):
    # returns gene symbol for a given transcript ID
    gene_symbol = dbfetchone.get_gene_symbol_from_transcriptID(transcript_id)
    return gene_symbol

def get_gene_symbol_from_refSeqGeneID(refSeqGeneID):
    # Returns the databases most up-to-date gene symbol for a given NG_ ID
    gene_symbol = dbfetchone.get_gene_symbol_from_refSeqGeneID(refSeqGeneID)
    return gene_symbol

def get_g_to_g_info():
    # Recovers the g_to_g data table
    table = dbfetchall.get_g_to_g_info()
    return table

def get_all_transcriptID():
    # Returns a list of transcript IDs in our database
    table = dbfetchall.get_all_transcriptID()
    return table

def get_RefSeqGeneID_from_lrgID(lrgID):
    # Get the relevant RefSeqGeneID for a given LRG ID
    rsgID = dbfetchone.get_RefSeqGeneID_from_lrgID(lrgID)
    return rsgID

def get_RefSeqTranscriptID_from_lrgTranscriptID(lrg_txID):
    rstID = dbfetchone.get_RefSeqTranscriptID_from_lrgTranscriptID(lrg_txID)
    return rstID

def get_lrgTranscriptID_from_RefSeqTranscriptID(rstID):
    lrg_tx = dbfetchone.get_lrgTranscriptID_from_RefSeqTranscriptID(rstID)
    return lrg_tx

def get_lrgProteinID_from_RefSeqProteinID(rs_p):
    lrg_p = dbfetchone.get_lrgProteinID_from_RefSeqProteinID(rs_p)
    return lrg_p

def get_lrgID_from_RefSeqGeneID(rsgID):	
    lrgID = dbfetchone.get_lrgID_from_RefSeqGeneID(rsgID)
    return lrgID

def get_refseqgene_info(refseqgene_id, primary_assembly):
    refseqgene_info = dbfetchone.get_refseqgene_info(refseqgene_id, primary_assembly)
    return refseqgene_info

def get_LRG_data_from_LRGid(lrg_id):
	LRG_data = dbfetchone.get_LRG_data_from_LRGid(lrg_id)
	return LRG_data

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
