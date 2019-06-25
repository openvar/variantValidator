from .logger import Logger
from . import utils
from .utils import handleCursor
from . import vvDBInsert
import re


class Database(vvDBInsert.Mixin):
    """
    This class contains and handles the mysql connections for the variant validator database.

    It now uses mixins, and the order of inheritance is
    vvDBInit.Mixin
       v
    vvDBGet.Mixin
       v
    vvDBInsert.Mixin
       v
    vvDatabase
    """
    # from dbquery
    @handleCursor
    def query_with_fetchone(self, entry):
        query = "SELECT refSeqID, description, transcriptVariant, currentVersion, hgncSymbol, utaSymbol, updated, " \
                "IF(updated < NOW() - INTERVAL 3 MONTH , 'true', 'false') FROM transcript_info WHERE " \
                "refSeqID = '%s'" % entry
        self.cursor.execute(query)
        row = self.cursor.fetchone()
        if row is None:
            row = ['none', 'No data']
            Logger.debug("No data returned from query " + str(query))
        return row

    # From data
    def data_add(self, accession, validator):
        """
        # Add accurate transcript descriptions to the database
        :param accession:
        :param validator:
        :return:
        """
        self.update_transcript_info_record(accession, validator)
        entry = self.in_entries(accession, 'transcript_info')
        return entry

    def in_entries(self, entry, table):
        """
        Retrieve transcript information
        :param entry:
        :param table:
        :return:
        """
        data = {}
        if table == 'transcript_info':
            row = self.query_with_fetchone(entry)
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

    def update_transcript_info_record(self, accession, validator):
        """
        Search Entrez for corresponding record for the RefSeq ID
        """
        # Prime these entries, just in case.
        previous_entry = self.in_entries(accession, 'transcript_info')
        accession = accession
        if 'none' not in previous_entry.keys():
            description = previous_entry['description']
            variant = previous_entry['variant']
            version = previous_entry['version']
            hgnc_symbol = previous_entry['hgnc_symbol']
            uta_symbol = previous_entry['uta_symbol']
        try:
            record = validator.entrez_efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            version = record.id
            description = record.description
            variant = '0'

            if 'transcript variant' in description:
                tv = re.search(r'transcript variant \w+', description)
                tv = str(tv.group(0))
                tv = tv.replace('transcript variant', '')
                variant = tv.strip()
                variant = variant.upper()  # Some tv descriptions are a or A
            else:
                variant = '0'

            # Get information from UTA
            try:
                uta_info = validator.hdp.get_tx_identity_info(version)
            except:
                version_ac_ver = version.split('.')
                version = version_ac_ver[0] + '.' + str(int(version_ac_ver[1]) - 1)
                uta_info = validator.hdp.get_tx_identity_info(version)

            uta_symbol = str(uta_info[6])

            # First perform a search against the input gene symbol or the symbol inferred from UTA
            initial = utils.hgnc_rest(path="/fetch/symbol/" + uta_symbol)
            # Check for a record
            if str(initial['record']['response']['numFound']) != '0':
                hgnc_symbol = uta_symbol
            # No record found, is it a previous symbol?
            else:
                # Search hgnc rest to see if symbol is out of date
                rest_data = utils.hgnc_rest(path="/search/prev_symbol/" + uta_symbol)
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
        table = 'transcript_info'

        # Update the transcript_info table (needs plugging in)
        returned_data = self.in_entries(version, table)
        # If the entry is not in the database add it
        if 'none' in returned_data:
            self.insert(version, query_info, table)
        # If the data in the entry has changed, update it
        else:
            self.update(version, query_info)
        return

    def update_refseqgene_loci(self, rsg_data):
        # First query the database
        entry_exists = self.get_refseq_data_by_refseq_id(rsg_data[0], rsg_data[2])
        if entry_exists[0] == 'none':
            self.insert_refseq_gene_data(rsg_data)
        else:
            self.update_refseq_gene_data(rsg_data)

    def update_lrg_rs_lookup(self, lrg_rs_lookup):
        # First query the database
        rsg_id = self.get_refseq_id_from_lrg_id(lrg_rs_lookup[0])
        if rsg_id == 'none':
            self.insert_refseq_gene_id_from_lrg_id(lrg_rs_lookup)

    def update_lrgt_rst(self, lrgtx_to_rst_id):
        # First query the database
        rst_id = self.get_refseq_transcript_id_from_lrg_transcript_id(lrgtx_to_rst_id[0])
        if rst_id == 'none':
            self.insert_lrg_transcript_data(lrgtx_to_rst_id)

    def update_lrg_p_rs_p_lookup(self, lrg_p, rs_p):
        # First query the database
        rsp_id = self.get_refseq_protein_id_from_lrg_protein_id(lrg_p)
        if rsp_id == 'none':
            self.insert_lrg_protein_data(lrg_p, rs_p)

    def ref_type_assign(self, accession):
        if 'NC_' in accession or 'NG_' in accession or 'NT_' in accession or 'NW_' in accession:
            ref_type = ':g.'
        elif accession.startswith('NM_'):
            ref_type = ':c.'
        elif accession.startswith('NR_'):
            ref_type = ':n.'
        elif accession.startswith('NP_'):
            ref_type = ':p.'
        elif accession.startswith('LRG_'):
            if 't' in accession:
                refseqtranscript_reference = self.get_refseq_transcript_id_from_lrg_transcript_id(accession)
                if refseqtranscript_reference.startswith('NM_'):
                    ref_type = ':c.'
                else:
                    ref_type = ':n.'
            elif '_p' in accession:
                ref_type = ':p.'
            else:
                ref_type = ':g.'
        else:
            # shouldn't reach this point
            raise Exception('Unable to recognise accession')
        return ref_type
