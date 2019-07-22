from . import utils
from .utils import handleCursor
from . import vvDBInsert
import re
import hgvs.exceptions
import logging
import json

logger = logging.getLogger(__name__)


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
            logger.debug("No data returned from query " + str(query))
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

        try:
            record = validator.entrez_efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        except IOError:
            raise utils.DatabaseConnectionError("Cannot retrieve data from NCBI Entrez")

        version = record.id
        description = record.description

        genbank_symbol = str(record.features[1].qualifiers['gene'][0])

        # Although it is obsolete, might still be in UTA database so would work in our case
        # if 'comment' in record.annotations:
        #     comment = record.annotations['comment']
        #     if 'WARNING' in comment and 'this sequence was replaced by' in comment:
        #         raise utils.ObsoleteSeqError("Sequence is obsolete in NCBI Entrez record")

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
        except hgvs.exceptions.HGVSDataNotAvailableError:
            version_ac_ver = version.split('.')
            version = version_ac_ver[0] + '.' + str(int(version_ac_ver[1]) - 1)
            try:
                uta_info = validator.hdp.get_tx_identity_info(version)
            except hgvs.exceptions.HGVSDataNotAvailableError:
                raise utils.DatabaseConnectionError("Cannot retrieve data from UTA database")

        uta_symbol = str(uta_info[6])
        symbol = uta_symbol
        if uta_symbol == '':
            # raise utils.ObsoleteSeqError("Cannot find UTA symbol, accession is likely obsolete")
            uta_symbol = 'unassigned'
            symbol = genbank_symbol

        hgnc_symbol = symbol

        try:
            # First perform a search against the input gene symbol or the symbol inferred from UTA
            initial = utils.hgnc_rest(path="/fetch/symbol/" + symbol)

            # Check for a record
            if str(initial['record']['response']['numFound']) == '0':
                # Search hgnc rest to see if symbol is out of date
                rest_data = utils.hgnc_rest(path="/search/prev_symbol/" + symbol)
                # If the name is correct no record will be found
                if rest_data['error'] == 'false' and int(rest_data['record']['response']['numFound']) != 0:
                    hgnc_symbol = rest_data['record']['response']['docs'][0]['symbol']
                    initial = utils.hgnc_rest(path="/fetch/symbol/" + hgnc_symbol)

            if hgnc_symbol != 'unassigned' and int(initial['record']['response']['numFound']) != 0:
                docs = initial['record']['response']['docs'][0]
                hgnc_id = ''
                entrez_id = ''
                ensembl_gene_id = ''
                omim_id = json.dumps([])
                ucsc_id = ''
                vega_id = ''
                ccds_id = json.dumps([])

                if 'hgnc_id' in docs:
                    hgnc_id = docs['hgnc_id']
                if 'entrez_id' in docs:
                    entrez_id = docs['entrez_id']
                if 'ensembl_gene_id' in docs:
                    ensembl_gene_id = docs['ensembl_gene_id']
                if 'omim_id' in docs:
                    omim_id = json.dumps(docs['omim_id'])
                if 'ucsc_id' in docs:
                    ucsc_id = docs['ucsc_id']
                if 'vega_id' in docs:
                    vega_id = docs['vega_id']
                if 'ccds_id' in docs:
                    ccds_id = json.dumps(docs['ccds_id'])

                gene_stable_ids = {
                    "hgnc_id": hgnc_id,
                    "entrez_id": entrez_id,
                    "ensembl_gene_id": ensembl_gene_id,
                    "omim_id": omim_id,
                    "ucsc_id": ucsc_id,
                    "vega_id": vega_id,
                    "ccds_id": ccds_id,
                    "hgnc_symbol": hgnc_symbol

                }
                gene_id_info = self.get_stable_gene_id_from_hgnc_id(gene_stable_ids["hgnc_id"])
                if gene_id_info[1] != 'No data':
                    self.update_gene_stable_ids(gene_stable_ids)
                else:
                    self.insert_gene_stable_ids(gene_stable_ids)

        except Exception as e:
            logger.debug("Except pass, %s", e)
            logger.info("Unable to connect to HGNC with symbol %s", symbol)

        # Query information
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

# <LICENSE>
# Copyright (C) 2019 VariantValidator Contributors
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
