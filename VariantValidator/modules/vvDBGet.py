from vvFunctions import handleCursor
from vvLogging import logger

class vvDBGet:
    def __init__(self,conn,cursor):
        # These are inherited by reference from the vvDatabase object.
        self.conn=conn
        self.cursor=cursor

    @handleCursor
    def execute(self,query):
        self.cursor.execute(query)
        row = self.cursor.fetchone()
        if row is None:
            logger.debug("No data returned from query "+str(query))
            row = ['none', 'No data']
        return row
    @handleCursor
    def executeAll(self,query):
        self.cursor.execute(query)
        rows = self.cursor.fetchone()
        if rows==[]:
            logger.debug("No data returned from query "+str(query))
            row = ['none', 'No data']
        return rows
    # from dbfetchone
    def get_utaSymbol(self,gene_symbol):
        query= "SELECT utaSymbol FROM transcript_info WHERE hgncSymbol = '%s'" %(gene_symbol)
        return self.execute(query)
    def get_hgncSymbol(self,gene_symbol):
        query= "SELECT hgncSymbol FROM transcript_info WHERE utaSymbol = '%s'" %(gene_symbol)
        return self.execute(query)
    def get_transcript_description(self,transcript_id):
        query= "SELECT description FROM transcript_info WHERE refSeqID = '%s'" %(transcript_id)
        return str(self.execute(query)[0])
    def get_gene_symbol_from_transcriptID(self,transcript_id):
        query = "SELECT hgncSymbol FROM transcript_info WHERE refSeqID = '%s'" %(transcript_id)
        return str(self.execute(query)[0])
    def get_refSeqGene_data_by_refSeqGeneID(self,refSeqGeneID, genomeBuild):
        query = "SELECT refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos, orientation, totalLength, chrPos, rsgPos, entrezID, hgncSymbol FROM refSeqGene_loci WHERE refSeqGeneID = '%s' AND genomeBuild = '%s'" %(refSeqGeneID, genomeBuild)
        return self.execute(query)
    def get_gene_symbol_from_refSeqGeneID(self,refSeqGeneID):
        query = "SELECT hgncSymbol FROM refSeqGene_loci WHERE refSeqGeneID = '%s'" %(refSeqGeneID)
        return self.execute(query)[0]
    def get_RefSeqGeneID_from_lrgID(self,lrgID):
        query = "SELECT RefSeqGeneID FROM LRG_RSG_lookup WHERE lrgID = '%s'" %(lrgID)
        return self.execute(query)[0]
    def get_RefSeqTranscriptID_from_lrgTranscriptID(self,lrgtxID):
        query = "SELECT RefSeqTranscriptID FROM LRG_transcripts WHERE LRGtranscriptID = '%s'" %(lrgtxID)
        return self.execute(query)[0]
    def	get_lrgTranscriptID_from_RefSeqTranscriptID(self,rstID):
        query = "SELECT LRGtranscriptID FROM LRG_transcripts WHERE RefSeqTranscriptID = '%s'" %(rstID)
        return self.execute(query)[0]
    def get_lrgID_from_RefSeqGeneID(self,rsgID):
        query = "SELECT lrgID, status FROM LRG_RSG_lookup WHERE RefSeqGeneID = '%s'" %(rsgID)
        return self.execute(query)
    def get_refseqgene_info(self,refseqgene_id, primary_assembly):
        query = "SELECT refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos FROM refSeqGene_loci WHERE refSeqGeneID = '%s' AND genomeBuild = '%s'" %(refseqgene_id, primary_assembly)
        return self.execute(query)
    def get_RefSeqProteinID_from_lrgProteinID(self,lrg_p):
        query = "SELECT RefSeqProteinID FROM LRG_proteins WHERE LRGproteinID = '%s'" %(lrg_p)
        return self.execute(query)[0]
    def get_lrgProteinID_from_RefSeqProteinID(self,rs_p):
        query = "SELECT LRGproteinID FROM LRG_proteins WHERE  RefSeqProteinID = '%s'" %(rs_p)
        return self.execute(query)[0]
    def get_LRG_data_from_LRGid(self,lrg_id):
        query = "SELECT * FROM LRG_RSG_lookup WHERE lrgID = '%s'" %(lrg_id)
        return self.execute(query)
    #from dbfetchall
    def get_transcript_info_for_gene(self,gene_symbol):
        query = "SELECT refSeqID, description, transcriptVariant, currentVersion, hgncSymbol, utaSymbol, updated, IF(updated < NOW() - INTERVAL 3 MONTH , 'true', 'false') FROM transcript_info WHERE hgncSymbol = '%s'" %(gene_symbol)
        return self.executeAll(query)
    def get_g_to_g_info(self):
        query = "SELECT refSeqGeneID, refSeqChromosomeID, startPos, endPos, orientation, hgncSymbol, genomeBuild FROM refSeqGene_loci"
        return self.executeAll(query)
    def get_all_transcriptID(self):
        query = "SELECT refSeqID FROM transcript_info"
        return self.executeAll(query)
    # Direct methods (GET)
    def get_uta_symbol(self,gene_symbol):
        # returns the UTA gene symbol when HGNC gene symbol is input
        return str(self.get_utaSymbol(gene_symbol)[0])
    def get_hgnc_symbol(self,gene_symbol):
        # returns the HGNC gene symbol when UTA gene symbol is input
        return str(self.get_hgncSymbol(gene_symbol)[0])
    # from external.py
    def get_urls(self,dict_out):
        # Provide direct links to reference sequence records
        # Add urls
        report_urls = {}
        if 'NM_' in dict_out['hgvs_transcript_variant'] or 'NR_' in dict_out['hgvs_transcript_variant']:
            report_urls['transcript'] = 'https://www.ncbi.nlm.nih.gov' \
                                        '/nuccore/%s' % dict_out['hgvs_transcript_variant'].split(':')[0]
        if 'NP_' in dict_out['hgvs_predicted_protein_consequence']['slr']:
            report_urls['protein'] = 'https://www.ncbi.nlm.nih.gov' \
                                     '/nuccore/%s' % str(dict_out['hgvs_predicted_protein_consequence']['slr']).split(':')[0]
        if 'NG_' in dict_out['hgvs_refseqgene_variant']:
            report_urls['refseqgene'] = 'https://www.ncbi.nlm.nih.gov' \
                                        '/nuccore/%s' % dict_out['hgvs_refseqgene_variant'].split(':')[0]
        if 'LRG' in dict_out['hgvs_lrg_variant']:
            lrg_id = dict_out['hgvs_lrg_variant'].split(':')[0]
            lrg_data = self.get_LRG_data_from_LRGid(lrg_id)
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

