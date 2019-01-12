import mysql.connector
from mysql.connector.pooling import MySQLConnectionPool
from vvLogging import logger
from vvFunctions import entrez_efetch,hgnc_rest,handleCursor
from vvDBInsert import vvDBInsert
from vvDBGet import vvDBGet
import re
import os

class vvDatabase:
    # This class contains and handles the mysql connections for the variant validator database.
    def __init__(self,val,dbConfig):
        self.conn = mysql.connector.pooling.MySQLConnectionPool(pool_size=10, **dbConfig)
        # self.cursor will be none UNLESS you're wrapping a function in @handlecursor, which automatically opens and
        # closes connections for you.
        self.cursor=None
        self.dbConfig=dbConfig
        # Construct database URL
        #'mysqlx://vvadmin:var1ant@127.0.0.1/validator'
        self.path="mysqlx://"+dbConfig["user"]+":"+dbConfig["password"]+"@"+dbConfig["host"]+"/"+dbConfig["database"]
        os.environ["VALIDATOR_DB_URL"]=self.path
        self.val=val
        self.insert = vvDBInsert(self.conn,self.cursor) # contains dbinsert, dbupdate
        self.get = vvDBGet(self.conn,self.cursor)       # contains dbfetchone, dbfetchall

    # from dbquery
    @handleCursor
    def query_with_fetchone(self,entry, table):
        #if table == 'transcript_info':
        query = "SELECT refSeqID, description, transcriptVariant, currentVersion, hgncSymbol, utaSymbol, updated, IF(updated < NOW() - INTERVAL 3 MONTH , 'true', 'false') FROM transcript_info WHERE refSeqID = '%s'" %(entry)
        self.cursor.execute(query)
        row = self.cursor.fetchone()
        if row is None:
            row = ['none', 'No data']
            logger.debug("No data returned from query "+str(query))
        return row
    # From data
    # Retrieve transcript information
    def in_entries(self,entry, table):
        # Use dbquery.py to connect to mysql and return the necessary data
        data={}
        if table == 'transcript_info':
            row = self.query_with_fetchone(entry, table)
            if row[0] == 'error':
                data['error'] = row[0]
                data['description'] = row[1]
            elif row[0] == 'none':
                data['none'] = row[0]
                data['description'] = row[1]
            else:
                data['accession'] = row[0]
                data['description'] = row[1]
                data['variant'] = row[2]
                data['version'] = row[3]
                data['hgnc_symbol'] = row[4]
                data['uta_symbol'] = row[5]
                data['updated'] = row[6]
                data['expiry'] = row[7]
        return data
    def update_transcript_info_record(self,accession, hdp):
        # Search Entrez for corresponding record for the RefSeq ID
        # Prime these entries, just in case.
        previous_entry = self.in_entries(accession, 'transcript_info')
        accession = accession
        description = previous_entry['description']
        variant = previous_entry['variant']
        version = previous_entry['version']
        hgnc_symbol = previous_entry['hgnc_symbol']
        uta_symbol = previous_entry['uta_symbol']
        try:
            record = entrez_efetch(self.val,db="nucleotide", id=accession, rettype="gb", retmode="text")
            version = record.id
            description = record.description
            variant = '0'

            if 'transcript variant' in description:
                tv = re.search('transcript variant \w+', description)
                tv = str(tv.group(0))
                tv = tv.replace('transcript variant', '')
                variant = tv.strip()
                variant = variant.upper()  # Some tv descriptions are a or A
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
            initial = hgnc_rest(path = "/fetch/symbol/" + uta_symbol)
            # Check for a record
            if str(initial['record']['response']['numFound']) != '0':
                hgnc_symbol = uta_symbol
            # No record found, is it a previous symbol?
            else:
                # Search hgnc rest to see if symbol is out of date
                rest_data = hgnc_rest(path = "/search/prev_symbol/" + uta_symbol)
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
            if not str(e) == '<urlopen error [Errno -2] Name or service not known>':
                # Issues with DNSSEC for the nih.gov
                raise

        # Query information
        # query_info = [accession, description, variant, version, hgnc_symbol, uta_symbol]
        query_info = [version, description, variant, version, hgnc_symbol, uta_symbol]
        table='transcript_info'

        # Update the transcript_info table (needs plugging in)
        returned_data = self.in_entries(version, table)
        # If the entry is not in the database add it
        if 'none' in returned_data:
            self.insert.add_entry(version, query_info, table)
        # If the data in the entry has changed, update it
        else:
            self.insert.update_entry(version, query_info, table)
        return

    def update_refSeqGene_loci(self,rsg_data):
        # First query the database
        entry_exists = self.get.get_refSeqGene_data_by_refSeqGeneID(rsg_data[0], rsg_data[2])
        if entry_exists[0] == 'none':
            self.insert.insert_refSeqGene_data(rsg_data)
        else:
            self.insert.update_refSeqGene_data(rsg_data)
    def update_lrg_rs_lookup(self,lrg_rs_lookup):
        # First query the database
        rsgID = self.get.get_RefSeqGeneID_from_lrgID(lrg_rs_lookup[0])
        if rsgID == 'none':
            self.insert.insert_RefSeqGeneID_from_lrgID(lrg_rs_lookup)
    def update_lrgt_rst(self,lrgtx_to_rstID):
        # First query the database
        rstID = self.get.get_RefSeqTranscriptID_from_lrgTranscriptID(lrgtx_to_rstID[0])
        if rstID == 'none':
            self.insert.insert_LRG_transcript_data(lrgtx_to_rstID)
    def update_lrg_p_rs_p_lookup(self,lrg_p, rs_p):
        # First query the database
        rspID = self.get.get_RefSeqProteinID_from_lrgProteinID(lrg_p)
        if rspID == 'none':
            self.insert.insert_LRG_protein_data(lrg_p, rs_p)
