from .utils import handleCursor
from . import vvDBGet


class Mixin(vvDBGet.Mixin):
    """
    This object is a function container for inserting objects into the database.
    """

    @handleCursor
    def insert(self, entry, data, table):
        # MySQL queries
        if table == 'transcript_info':
            accession = entry
            description = data[1]
            variant = data[2]
            version = data[3]
            hgnc_symbol = data[4]
            uta_symbol = data[5]
            query = "INSERT INTO transcript_info(refSeqID, description, transcriptVariant, currentVersion, " \
                    "hgncSymbol, utaSymbol, updated) VALUES(%s,%s, %s, %s, %s, %s, NOW())"
            self.cursor.execute(query, (accession, description, variant, version, hgnc_symbol, uta_symbol))
        # Query report
        if self.cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'

        # Commit and close connection (?close?)
        self.conn.commit()
        return success

    @handleCursor
    def insert_refseq_gene_data(self, rsg_data):
        query = "INSERT INTO refSeqGene_loci(refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos, " \
                "orientation, totalLength, chrPos, rsgPos, entrezID, hgncSymbol, updated) " \
                "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, NOW())"
        self.cursor.execute(query, (rsg_data[0], rsg_data[1], rsg_data[2], rsg_data[3], rsg_data[4], rsg_data[5],
                                    rsg_data[6], rsg_data[7], rsg_data[8], rsg_data[9], rsg_data[10]))
        # Query report
        if self.cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'
        # Commit and close connection
        self.conn.commit()
        return success

    @handleCursor
    def insert_refseq_gene_id_from_lrg_id(self, lrg_rs_lookup):
        query = "INSERT INTO LRG_RSG_lookup(lrgID, hgncSymbol, RefSeqGeneID, status) VALUES (%s,%s,%s,%s)"
        self.cursor.execute(query, (lrg_rs_lookup[0], lrg_rs_lookup[1], lrg_rs_lookup[2], lrg_rs_lookup[3]))
        # Query report
        if self.cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'
        # Commit and close connection
        self.conn.commit()
        return success

    @handleCursor
    def insert_lrg_transcript_data(self, lrgtx_to_rst_id):
        query = "INSERT INTO LRG_transcripts(LRGtranscriptID, RefSeqTranscriptID) VALUES (%s,%s)"
        self.cursor.execute(query, (lrgtx_to_rst_id[0], lrgtx_to_rst_id[1]))
        # Query report
        if self.cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'

        # Commit and close connection
        self.conn.commit()
        return success

    @handleCursor
    def insert_lrg_protein_data(self, lrg_p, rs_p):
        query = "INSERT INTO LRG_proteins(LRGproteinID, RefSeqProteinID) VALUES (%s,%s)"
        self.cursor.execute(query, (lrg_p, rs_p))
        # Query report
        if self.cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'

        # Commit and close connection
        self.conn.commit()
        return success

    @handleCursor
    def insert_gene_stable_ids(self, data):
        query = "INSERT INTO stableGeneIds(hgnc_id, hgnc_symbol, entrez_id, ensembl_gene_id, omim_id, ucsc_id, " \
                "vega_id, ccds_ids) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)"
        self.cursor.execute(query, (
            data['hgnc_id'],
            data['hgnc_symbol'],
            data['entrez_id'],
            data['ensembl_gene_id'],
            data['omim_id'],
            data['ucsc_id'],
            data['vega_id'],
            data['ccds_id']
        ))

        if self.cursor.lastrowid:
            success = 'true'
        else:
            success = 'unknown error'

        self.conn.commit()
        return success

    @handleCursor
    def update(self, entry, data):
        accession = entry
        description = data[1]
        variant = data[2]
        version = data[3]
        hgnc_symbol = data[4]
        uta_symbol = data[5]
        query = "UPDATE transcript_info SET description=%s, transcriptVariant=%s, currentVersion=%s, hgncSymbol=%s, " \
                "utaSymbol=%s, updated=NOW() WHERE refSeqID = %s"
        self.cursor.execute(query, (description, variant, version, hgnc_symbol, uta_symbol, accession))
        success = 'true'
        self.conn.commit()
        return success

    @handleCursor
    def update_refseq_gene_data(self, rsg_data):
        query = "UPDATE refSeqGene_loci SET hgncSymbol=%s, updated=NOW() WHERE refSeqGeneID=%s"
        self.cursor.execute(query, (rsg_data[10], rsg_data[0]))
        success = 'true'
        self.conn.commit()
        return success

    @handleCursor
    def update_gene_stable_ids(self, gene_stable_ids):

        # Insert or update combined statement
        query = "UPDATE stableGeneIds SET hgnc_symbol=%s, entrez_id=%s, ensembl_gene_id=%s, omim_id=%s, ucsc_id=%s, " \
                "vega_id=%s, ccds_ids=%s WHERE hgnc_id=%s"

        self.cursor.execute(query, (
            gene_stable_ids["hgnc_symbol"],
            gene_stable_ids["entrez_id"],
            gene_stable_ids["ensembl_gene_id"],
            gene_stable_ids["omim_id"],
            gene_stable_ids["ucsc_id"],
            gene_stable_ids["vega_id"],
            gene_stable_ids["ccds_id"],
            gene_stable_ids["hgnc_id"]
        ))
        success = 'true'
        self.conn.commit()
        return success

