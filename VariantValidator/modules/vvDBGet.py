import logging
from .utils import handleCursor
from . import vvDBInit

logger = logging.getLogger(__name__)
LRG_TX_LINK = {}

class Mixin(vvDBInit.Mixin):
    """
    Most of the functions in DBGet generate queries for retrieving data from the databases.
    """

    @handleCursor
    def execute(self, *query_args):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        cursor.execute(*query_args)
        row = cursor.fetchone()
        if row is None:
            logger.debug("No data returned from query " + str(query_args))
            row = ['none', 'No data']

        # Close conn
        cursor.close()
        conn.close()
        return row

    @handleCursor
    def execute_all(self, *query_args):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        cursor.execute(*query_args)
        rows = cursor.fetchall()
        if not rows:
            logger.debug("No data returned from query " + str(query_args))
            rows = ['none', 'No data']

        # Close conn
        cursor.close()
        conn.close()
        return rows

    # from dbfetchone
    def get_uta(self, gene_symbol):
        query = "SELECT utaSymbol FROM transcript_info WHERE hgncSymbol = %s"
        return self.execute(query,(gene_symbol,))

    def get_hgnc(self, gene_symbol):
        query = "SELECT hgncSymbol FROM transcript_info WHERE utaSymbol = %s"
        return self.execute(query,(gene_symbol,))

    def get_transcript_description(self, transcript_id):
        query = "SELECT description FROM transcript_info WHERE refSeqID = %s"
        return str(self.execute(query,(transcript_id,))[0])

    def get_transcript_annotation(self, transcript_id):
        query = "SELECT transcriptVariant FROM transcript_info WHERE refSeqID = %s"
        return str(self.execute(query,(transcript_id,))[0])

    def get_gene_symbol_from_transcript_id(self, transcript_id):
        query = "SELECT hgncSymbol FROM transcript_info WHERE refSeqID = %s"
        return str(self.execute(query,(transcript_id,))[0])

    def get_refseq_data_by_refseq_id(self, refseq_id, genome_build):
        query = "SELECT refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos, orientation, totalLength, " \
                "chrPos, rsgPos, entrezID, hgncSymbol FROM refSeqGene_loci WHERE refSeqGeneID = %s " \
                "AND genomeBuild = %s"
        return self.execute(query,(refseq_id, genome_build))

    def get_gene_symbol_from_refseq_id(self, refseq_id):
        query = "SELECT hgncSymbol FROM refSeqGene_loci WHERE refSeqGeneID = %s"
        return self.execute(query,(refseq_id,))[0]

    def get_refseq_id_from_lrg_id(self, lrg_id):
        query = "SELECT RefSeqGeneID FROM LRG_RSG_lookup WHERE lrgID = %s"
        return self.execute(query,(lrg_id,))[0]

    def get_refseq_transcript_id_from_lrg_transcript_id(self, lrg_tx_id):
        query = "SELECT RefSeqTranscriptID FROM LRG_transcripts WHERE LRGtranscriptID = %s"
        return self.execute(query,(lrg_tx_id,))[0]

    def get_lrg_transcript_id_from_refseq_transcript_id(self, rst_id):
        if not LRG_TX_LINK:
            query = "SELECT RefSeqTranscriptID,LRGtranscriptID FROM LRG_transcripts"
            lrg_dat = self.execute_all(query)
            for dat in lrg_dat:
                LRG_TX_LINK[dat[0]] = dat[1]
        return LRG_TX_LINK.get(rst_id,'none')

    def get_lrg_id_from_refseq_gene_id(self, rsg_id):
        query = "SELECT lrgID, status FROM LRG_RSG_lookup WHERE RefSeqGeneID = %s"
        return self.execute(query,(rsg_id,))

    def get_refseqgene_info(self, refseqgene_id, primary_assembly):
        query = "SELECT refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos FROM refSeqGene_loci " \
                "WHERE refSeqGeneID = %s AND genomeBuild = %s"
        return self.execute(query,(refseqgene_id, primary_assembly))

    def get_refseq_protein_id_from_lrg_protein_id(self, lrg_p):
        query = "SELECT RefSeqProteinID FROM LRG_proteins WHERE LRGproteinID = %s"
        return self.execute(query,(lrg_p,))[0]

    def get_lrg_protein_id_from_ref_seq_protein_id(self, rs_p):
        query = "SELECT LRGproteinID FROM LRG_proteins WHERE  RefSeqProteinID = %s"
        return self.execute(query,(rs_p,))[0]

    def get_lrg_data_from_lrg_id(self, lrg_id):
        query = "SELECT * FROM LRG_RSG_lookup WHERE lrgID = %s"
        return self.execute(query,(lrg_id,))

    def get_transcript_info_for_gene(self, gene_symbol):
        query = "SELECT refSeqID, description, transcriptVariant, currentVersion, hgncSymbol, utaSymbol, " \
                "updated, IF(updated < NOW() - INTERVAL 3 MONTH , 'true', 'false') FROM transcript_info " \
                "WHERE hgncSymbol = %s"
        return self.execute_all(query,(gene_symbol,))

    def get_g_to_g_info(self, rsg_id=None, gen_id=None, start=None, end=None):
        """
        Return a set of RSG to genome mapping data.
        Can be all such mappings or can be limited to either data for a specific
        RSG or for a specific genomic ref id, and optional location.
        """
        query = "SELECT refSeqGeneID, refSeqChromosomeID, startPos, endPos, orientation, hgncSymbol, " \
                "genomeBuild FROM refSeqGene_loci"
        if rsg_id:
            query = query + f" WHERE refSeqGeneID= %s "
            query_vals = (rsg_id,)
        elif gen_id:
            query = query + f" WHERE refSeqChromosomeID= %s "
            query_vals = (gen_id,)
            if start:
                query = query + f"AND startPos <= %s "
                query_vals = query_vals + (str(start),)
            if end:
                query = query + f"AND endPos >= %s "
                query_vals = query_vals + (str(end),)
        return self.execute_all(query,query_vals)

    def get_all_transcript_id(self):
        query = "SELECT refSeqID FROM transcript_info"
        return self.execute_all(query)

    def get_stable_gene_id_info(self, hgnc_symbol):
        query = "SELECT * FROM stableGeneIds WHERE hgnc_symbol = %s"
        return self.execute(query,(hgnc_symbol,))

    def get_stable_gene_id_from_hgnc_id(self, hgnc_id):
        query = "SELECT * FROM stableGeneIds WHERE hgnc_id = %s"
        return self.execute(query,(hgnc_id,))

    def get_transcripts_from_annotations(self, statement):
        testval = "%" + statement + "%"
        query = "SELECT * FROM transcript_info WHERE transcriptVariant LIKE %s"

        # print("My query is: ", query)
        return self.execute_all(query,(testval,))

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

        # Refseq
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
        
        # Ensembl
        # When selectec_assembly is GRCh37
        if 'ENST' in dict_out['hgvs_transcript_variant'] and str(dict_out['selected_assembly']).lower() == 'grch37':
            report_urls['transcript'] = 'https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?' \
                                        'db=core;t=%s' % dict_out['hgvs_transcript_variant'].split(':')[0]
        if 'ENSP' in str(dict_out['hgvs_predicted_protein_consequence']['slr']) and str(dict_out['selected_assembly']).lower() == 'grch37':
            report_urls['protein'] = 'https://grch37.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?' \
                                     'db=core;p=%s' % str(
                                        dict_out['hgvs_predicted_protein_consequence']['slr']).split(':')[0]
        
        # When selected_assembly is GRCh38
        if 'ENST' in dict_out['hgvs_transcript_variant'] and str(dict_out['selected_assembly']).lower() == 'grch38':
            report_urls['transcript'] = 'https://www.ensembl.org/Homo_sapiens/Transcript/Summary?' \
                                        'db=core;t=%s' % dict_out['hgvs_transcript_variant'].split(':')[0]
        if 'ENSP' in str(dict_out['hgvs_predicted_protein_consequence']['slr']) and str(dict_out['selected_assembly']).lower() == 'grch38':
            report_urls['protein'] = 'https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?' \
                                     'db=core;p=%s' % str(
                                        dict_out['hgvs_predicted_protein_consequence']['slr']).split(':')[0]                        
        # "http://www.ensembl.org/id/" ? What about historic versions?????

        return report_urls

# <LICENSE>
# Copyright (C) 2016-2025 VariantValidator Contributors
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
