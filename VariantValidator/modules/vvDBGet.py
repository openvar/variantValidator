import logging
from .utils import handleCursor
from . import vvDBInit

logger = logging.getLogger(__name__)


class Mixin(vvDBInit.Mixin):
    """
    Most of the functions in DBGet generate queries for retrieving data from the databases.
    """

    @handleCursor
    def execute(self, query):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        cursor.execute(query)
        row = cursor.fetchone()
        if row is None:
            logger.debug("No data returned from query " + str(query))
            row = ['none', 'No data']

        # Close conn
        cursor.close()
        conn.close()
        return row

    @handleCursor
    def execute_all(self, query):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        cursor.execute(query)
        rows = cursor.fetchall()
        if not rows:
            logger.debug("No data returned from query " + str(query))
            rows = ['none', 'No data']

        # Close conn
        cursor.close()
        conn.close()
        return rows

    # from dbfetchone
    def get_uta(self, gene_symbol):
        query = "SELECT utaSymbol FROM transcript_info WHERE hgncSymbol = '%s'" % gene_symbol
        return self.execute(query)

    def get_hgnc(self, gene_symbol):
        query = "SELECT hgncSymbol FROM transcript_info WHERE utaSymbol = '%s'" % gene_symbol
        return self.execute(query)

    def get_transcript_description(self, transcript_id):
        query = "SELECT description FROM transcript_info WHERE refSeqID = '%s'" % transcript_id
        return str(self.execute(query)[0])

    def get_transcript_annotation(self, transcript_id):
        query = "SELECT transcriptVariant FROM transcript_info WHERE refSeqID = '%s'" % transcript_id
        return str(self.execute(query)[0])

    def get_gene_symbol_from_transcript_id(self, transcript_id):
        query = "SELECT hgncSymbol FROM transcript_info WHERE refSeqID = '%s'" % transcript_id
        return str(self.execute(query)[0])

    def get_refseq_data_by_refseq_id(self, refseq_id, genome_build):
        query = "SELECT refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos, orientation, totalLength, " \
                "chrPos, rsgPos, entrezID, hgncSymbol FROM refSeqGene_loci WHERE refSeqGeneID = '%s' " \
                "AND genomeBuild = '%s'" % (refseq_id, genome_build)
        return self.execute(query)

    def get_gene_symbol_from_refseq_id(self, refseq_id):
        query = "SELECT hgncSymbol FROM refSeqGene_loci WHERE refSeqGeneID = '%s'" % refseq_id
        return self.execute(query)[0]

    def get_refseq_id_from_lrg_id(self, lrg_id):
        query = "SELECT RefSeqGeneID FROM LRG_RSG_lookup WHERE lrgID = '%s'" % lrg_id
        return self.execute(query)[0]

    def get_refseq_transcript_id_from_lrg_transcript_id(self, lrg_tx_id):
        query = "SELECT RefSeqTranscriptID FROM LRG_transcripts WHERE LRGtranscriptID = '%s'" % lrg_tx_id
        return self.execute(query)[0]

    def get_lrg_transcript_id_from_refseq_transcript_id(self, rst_id):
        query = "SELECT LRGtranscriptID FROM LRG_transcripts WHERE RefSeqTranscriptID = '%s'" % rst_id
        return self.execute(query)[0]

    def get_lrg_id_from_refseq_gene_id(self, rsg_id):
        query = "SELECT lrgID, status FROM LRG_RSG_lookup WHERE RefSeqGeneID = '%s'" % rsg_id
        return self.execute(query)

    def get_refseqgene_info(self, refseqgene_id, primary_assembly):
        query = "SELECT refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos FROM refSeqGene_loci " \
                "WHERE refSeqGeneID = '%s' AND genomeBuild = '%s'" % (refseqgene_id, primary_assembly)
        return self.execute(query)

    def get_refseq_protein_id_from_lrg_protein_id(self, lrg_p):
        query = "SELECT RefSeqProteinID FROM LRG_proteins WHERE LRGproteinID = '%s'" % lrg_p
        return self.execute(query)[0]

    def get_lrg_protein_id_from_ref_seq_protein_id(self, rs_p):
        query = "SELECT LRGproteinID FROM LRG_proteins WHERE  RefSeqProteinID = '%s'" % rs_p
        return self.execute(query)[0]

    def get_lrg_data_from_lrg_id(self, lrg_id):
        query = "SELECT * FROM LRG_RSG_lookup WHERE lrgID = '%s'" % lrg_id
        return self.execute(query)

    def get_transcript_info_for_gene(self, gene_symbol):
        query = "SELECT refSeqID, description, transcriptVariant, currentVersion, hgncSymbol, utaSymbol, " \
                "updated, IF(updated < NOW() - INTERVAL 3 MONTH , 'true', 'false') FROM transcript_info " \
                "WHERE hgncSymbol = '%s'" % gene_symbol
        return self.execute_all(query)

    def get_g_to_g_info(self):
        query = "SELECT refSeqGeneID, refSeqChromosomeID, startPos, endPos, orientation, hgncSymbol, " \
                "genomeBuild FROM refSeqGene_loci"
        return self.execute_all(query)

    def get_all_transcript_id(self):
        query = "SELECT refSeqID FROM transcript_info"
        return self.execute_all(query)

    def get_stable_gene_id_info(self, hgnc_symbol):
        query = "SELECT * FROM stableGeneIds WHERE hgnc_symbol = '%s'" % hgnc_symbol
        return self.execute(query)

    def get_stable_gene_id_from_hgnc_id(self, hgnc_id):
        query = "SELECT * FROM stableGeneIds WHERE hgnc_id = '%s'" % hgnc_id
        return self.execute(query)

    def get_db_version(self):
        """
        :return: current version of the Validator database
        """
        query = "SELECT current_version FROM version"
        return self.execute(query)

    # Direct methods (GET)
    def get_uta_symbol(self, gene_symbol):
        # returns the UTA gene symbol when HGNC gene symbol is input
        return str(self.get_uta(gene_symbol)[0])

    def get_hgnc_symbol(self, gene_symbol):
        # returns the HGNC gene symbol when UTA gene symbol is input
        return str(self.get_hgnc(gene_symbol)[0])

    # from external.py
    def get_urls(self, dict_out):
        # Provide direct links to reference sequence records
        # Add urls
        report_urls = {}
        if 'NM_' in dict_out['hgvs_transcript_variant'] or 'NR_' in dict_out['hgvs_transcript_variant']:
            report_urls['transcript'] = 'https://www.ncbi.nlm.nih.gov' \
                                        '/nuccore/%s' % dict_out['hgvs_transcript_variant'].split(':')[0]
        if 'NP_' in str(dict_out['hgvs_predicted_protein_consequence']['slr']):
            report_urls['protein'] = 'https://www.ncbi.nlm.nih.gov' \
                                     '/nuccore/%s' % str(
                dict_out['hgvs_predicted_protein_consequence']['slr']).split(':')[0]
        if 'NG_' in dict_out['hgvs_refseqgene_variant']:
            report_urls['refseqgene'] = 'https://www.ncbi.nlm.nih.gov' \
                                        '/nuccore/%s' % dict_out['hgvs_refseqgene_variant'].split(':')[0]
        if 'LRG' in dict_out['hgvs_lrg_variant']:
            lrg_id = dict_out['hgvs_lrg_variant'].split(':')[0]
            lrg_data = self.get_lrg_data_from_lrg_id(lrg_id)
            lrg_status = str(lrg_data[4])
            if lrg_status == 'public':
                report_urls['lrg'] = 'http://ftp.ebi.ac.uk/pub' \
                                     '/databases/lrgex/%s.xml' % dict_out['hgvs_lrg_variant'].split(':')[0]
            else:
                report_urls['lrg'] = 'http://ftp.ebi.ac.uk' \
                                     '/pub/databases/lrgex' \
                                     '/pending/%s.xml' % dict_out['hgvs_lrg_variant'].split(':')[0]
        # Ensembl needs to be added at a later date
        # "http://www.ensembl.org/id/" ? What about historic versions?????

        return report_urls

# <LICENSE>
# Copyright (C) 2016-2021 VariantValidator Contributors
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
