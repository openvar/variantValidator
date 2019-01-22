from vvFunctions import handleCursor

class vvDBInsert:
    '''
    This object is a function container for inserting objects into the database.
    '''
    def __init__(self,db):
        # These are inherited by reference from the vvDatabase object.
        self.db=db
    # Add new entry
    def add_entry(self,entry, data, table):
        return self.insert(entry, data, table)
    def insert_transcript_loci(self,add_data, primary_assembly):
        return self.insert_transcript_loci(add_data, primary_assembly)

    #from dbinsert
    @handleCursor
    def insert(self,entry, data, table):
        # MySQL queries
        if table == 'transcript_info':
            accession = entry
            description = data[1]
            variant = data[2]
            version = data[3]
            hgnc_symbol = data[4]
            uta_symbol	= data[5]
            query = "INSERT INTO transcript_info(refSeqID, description, transcriptVariant, currentVersion, hgncSymbol, utaSymbol, updated) VALUES(%s,%s, %s, %s, %s, %s, NOW())"
            self.db.cursor.execute(query, (accession, description, variant, version, hgnc_symbol, uta_symbol))
        # Query report
        if self.db.cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'

        # Commit and close connection
        self.db.conn.commit()
        return success
    @handleCursor
    def insert_refSeqGene_data(self,rsg_data):
        query = "INSERT INTO refSeqGene_loci(refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos, orientation, totalLength, chrPos, rsgPos, entrezID, hgncSymbol, updated) VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, NOW())"
        self.db.cursor.execute(query, (rsg_data[0], rsg_data[1], rsg_data[2], rsg_data[3], rsg_data[4], rsg_data[5], rsg_data[6], rsg_data[7], rsg_data[8], rsg_data[9], rsg_data[10]))
        # Query report
        if self.db.cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'
        # Commit and close connection
        self.db.conn.commit()
        return success
    @handleCursor
    def insert_RefSeqGeneID_from_lrgID(self,lrg_rs_lookup):
        query = "INSERT INTO LRG_RSG_lookup(lrgID, hgncSymbol, RefSeqGeneID, status) VALUES(%s,%s,%s,%s)"
        self.db.cursor.execute(query, (lrg_rs_lookup[0], lrg_rs_lookup[1], lrg_rs_lookup[2], lrg_rs_lookup[3]))
        # Query report
        if self.db.cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'
        # Commit and close connection
        self.db.conn.commit()
        return success
    @handleCursor
    def insert_LRG_transcript_data(self,lrgtx_to_rstID):
        query = "INSERT INTO LRG_transcripts(LRGtranscriptID, RefSeqTranscriptID) VALUES(%s,%s)"
        self.db.cursor.execute(query, (lrgtx_to_rstID[0], lrgtx_to_rstID[1]))
        # Query report
        if self.db.cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'

        # Commit and close connection
        self.db.conn.commit()
        return success
    @handleCursor
    def insert_LRG_protein_data(self,lrg_p, rs_p):
        query = "INSERT INTO LRG_proteins(LRGproteinID, RefSeqProteinID) VALUES(%s,%s)"
        self.db.cursor.execute(query, (lrg_p, rs_p))
        # Query report
        if self.db.cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'

        # Commit and close connection
        self.db.conn.commit()
        return success
    # from dbupdate
    @handleCursor
    def update(self,entry, data, table):
        # MySQL queries
        #if table == 'transcript_info':
        accession = entry
        description = data[1]
        variant = data[2]
        version = data[3]
        hgnc_symbol = data[4]
        uta_symbol	= data[5]
        query = "UPDATE transcript_info SET description=%s, transcriptVariant=%s, currentVersion=%s, hgncSymbol=%s, utaSymbol=%s, updated=NOW() WHERE refSeqID = %s"
        self.db.cursor.execute(query, (description, variant, version, hgnc_symbol, uta_symbol, accession))
        success = 'true'
        self.db.conn.commit()
        return success
        # 'true'??? check this.
    @handleCursor
    def update_refSeqGene_data(self,rsg_data):
        query = "UPDATE refSeqGene_loci SET hgncSymbol=%s, updated=NOW() WHERE refSeqGeneID=%s"
        self.db.cursor.execute(query, (rsg_data[10], rsg_data[0]))
        success = 'true'
        self.db.conn.commit()
        return success
    # Update entries
    def update_entry(self,entry, data, table):
        return self.update(entry, data, table)
