from .utils import handleCursor
from . import vvDBGet


class Mixin(vvDBGet.Mixin):
    """
    This object is a function container for inserting objects into the database.
    """

    @handleCursor
    def insert(self, entry, data, table):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

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
            cursor.execute(query, (accession, description, variant, version, hgnc_symbol, uta_symbol))
        # Query report
        if cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'

        # Commit and close connection
        conn.commit()
        cursor.close()
        conn.close()
        return success

    @handleCursor
    def insert_refseq_gene_data(self, rsg_data):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = "INSERT INTO refSeqGene_loci(refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos, " \
                "orientation, totalLength, chrPos, rsgPos, entrezID, hgncSymbol, updated) " \
                "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, NOW())"
        cursor.execute(query, (rsg_data[0], rsg_data[1], rsg_data[2], int(rsg_data[3]), int(rsg_data[4]), rsg_data[5],
                                    int(rsg_data[6]), rsg_data[7], rsg_data[8], int(rsg_data[9]), rsg_data[10]))
        # Query report
        if cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'

        # Commit and close connection
        conn.commit()
        cursor.close()
        conn.close()
        return success

    @handleCursor
    def insert_refseq_gene_id_from_lrg_id(self, lrg_rs_lookup):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = "INSERT INTO LRG_RSG_lookup(lrgID, hgncSymbol, RefSeqGeneID, status) VALUES (%s,%s,%s,%s)"
        cursor.execute(query, (lrg_rs_lookup[0], lrg_rs_lookup[1], lrg_rs_lookup[2], lrg_rs_lookup[3]))
        # Query report
        if cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'

        # Commit and close connection
        conn.commit()
        cursor.close()
        conn.close()
        return success

    @handleCursor
    def insert_lrg_transcript_data(self, lrgtx_to_rst_id):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = "INSERT INTO LRG_transcripts(LRGtranscriptID, RefSeqTranscriptID) VALUES (%s,%s)"
        cursor.execute(query, (lrgtx_to_rst_id[0], lrgtx_to_rst_id[1]))
        # Query report
        if cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'

        # Commit and close connection
        conn.commit()
        cursor.close()
        conn.close()
        return success

    @handleCursor
    def insert_lrg_protein_data(self, lrg_p, rs_p):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = "INSERT INTO LRG_proteins(LRGproteinID, RefSeqProteinID) VALUES (%s,%s)"
        cursor.execute(query, (lrg_p, rs_p))
        # Query report
        if cursor.lastrowid:
            success = 'true'
        else:
            success = 'Unknown error'

        # Commit and close connection
        conn.commit()
        cursor.close()
        conn.close()
        return success

    @handleCursor
    def insert_gene_stable_ids(self, data):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = "INSERT INTO stableGeneIds(hgnc_id, hgnc_symbol, entrez_id, ensembl_gene_id, omim_id, ucsc_id, " \
                "vega_id, ccds_ids) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)"
        cursor.execute(query, (
            data['hgnc_id'],
            data['hgnc_symbol'],
            data['entrez_id'],
            data['ensembl_gene_id'],
            data['omim_id'],
            data['ucsc_id'],
            data['vega_id'],
            data['ccds_id']
        ))

        if cursor.lastrowid:
            success = 'true'
        else:
            success = 'unknown error'

        # Commit and close connection
        conn.commit()
        cursor.close()
        conn.close()
        return success

    @handleCursor
    def update(self, entry, data):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        accession = entry
        description = data[1]
        variant = data[2]
        version = data[3]
        hgnc_symbol = data[4]
        uta_symbol = data[5]

        query = "UPDATE transcript_info SET description=%s, transcriptVariant=%s, currentVersion=%s, hgncSymbol=%s, " \
                "utaSymbol=%s, updated=NOW() WHERE refSeqID = %s"
        cursor.execute(query, (description, variant, version, hgnc_symbol, uta_symbol, accession))
        success = 'true'

        # Commit and close connection
        conn.commit()
        cursor.close()
        conn.close()
        return success

    @handleCursor
    def update_refseq_gene_data(self, rsg_data):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = "UPDATE refSeqGene_loci SET hgncSymbol=%s, updated=NOW() WHERE refSeqGeneID=%s"
        cursor.execute(query, (rsg_data[10], rsg_data[0]))
        success = 'true'

        # Commit and close connection
        conn.commit()
        cursor.close()
        conn.close()
        return success

    @handleCursor
    def update_gene_stable_ids(self, gene_stable_ids):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        # Insert or update combined statement
        query = "UPDATE stableGeneIds SET hgnc_symbol=%s, entrez_id=%s, ensembl_gene_id=%s, omim_id=%s, ucsc_id=%s, " \
                "vega_id=%s, ccds_ids=%s WHERE hgnc_id=%s"

        cursor.execute(query, (
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

        # Commit and close connection
        conn.commit()
        cursor.close()
        conn.close()
        return success

    @handleCursor
    def update_db_version(self, db_version):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        db_version = [db_version]
        query = "UPDATE version SET current_version=%s"
        cursor.execute(query, db_version)
        success = 'true'

        # Commit and close connection
        conn.commit()
        cursor.close()
        conn.close()
        return success

# ---------------------------------------------------------------------------
# SQLite INSERT backend
# ---------------------------------------------------------------------------

from VariantValidator.modules.vvDBGet import SQLiteDBGet as _SQLiteDBGet


class SQLiteDBInsert(_SQLiteDBGet):
    """SQLite equivalent of the MySQL Mixin INSERT/UPDATE queries.

    SQL translation rules applied:
    * ``%s``                       → ``?``
    * ``NOW()``                    → ``datetime('now')``
    * ``INSERT INTO ... ON DUPLICATE KEY UPDATE ...``
                                   → ``INSERT OR REPLACE INTO ...``
    * For UPDATE statements: same substitutions.
    """

    # ------------------------------------------------------------------
    # transcript_info
    # ------------------------------------------------------------------

    def insert(self, entry, data, table):
        """Insert a row into *table* (currently only 'transcript_info')."""
        if table == 'transcript_info':
            accession = entry
            description = data[1]
            variant = data[2]
            version = data[3]
            hgnc_symbol = data[4]
            uta_symbol = data[5]
            query = (
                "INSERT OR REPLACE INTO transcript_info"
                "(refSeqID, description, transcriptVariant, currentVersion,"
                " hgncSymbol, utaSymbol, updated)"
                " VALUES(?, ?, ?, ?, ?, ?, datetime('now'))"
            )
            conn = self.get_conn()
            cursor = conn.cursor()
            cursor.execute(query, (accession, description, variant, version, hgnc_symbol, uta_symbol))
            success = 'true' if cursor.lastrowid else 'Unknown error'
            conn.commit()
            cursor.close()
            conn.close()
            return success
        raise ValueError(f"Unsupported table for insert(): {table!r}")

    def insert_transcript_info(
        self,
        refseq_id,
        description,
        transcript_variant,
        current_version,
        hgnc_symbol,
        uta_symbol,
    ):
        """INSERT OR REPLACE a transcript_info row."""
        query = (
            "INSERT OR REPLACE INTO transcript_info"
            "(refSeqID, description, transcriptVariant, currentVersion,"
            " hgncSymbol, utaSymbol, updated)"
            " VALUES(?, ?, ?, ?, ?, ?, datetime('now'))"
        )
        conn = self.get_conn()
        cursor = conn.cursor()
        cursor.execute(query, (refseq_id, description, transcript_variant, current_version, hgnc_symbol, uta_symbol))
        success = 'true' if cursor.lastrowid else 'Unknown error'
        conn.commit()
        cursor.close()
        conn.close()
        return success

    def update(self, entry, data):
        """UPDATE a transcript_info row."""
        accession = entry
        description = data[1]
        variant = data[2]
        version = data[3]
        hgnc_symbol = data[4]
        uta_symbol = data[5]
        query = (
            "UPDATE transcript_info SET description=?, transcriptVariant=?,"
            " currentVersion=?, hgncSymbol=?, utaSymbol=?, updated=datetime('now')"
            " WHERE refSeqID = ?"
        )
        conn = self.get_conn()
        cursor = conn.cursor()
        cursor.execute(query, (description, variant, version, hgnc_symbol, uta_symbol, accession))
        conn.commit()
        cursor.close()
        conn.close()
        return 'true'

    # ------------------------------------------------------------------
    # refSeqGene_loci
    # ------------------------------------------------------------------

    def insert_refseq_gene_data(self, rsg_data):
        """INSERT OR REPLACE into refSeqGene_loci."""
        query = (
            "INSERT OR REPLACE INTO refSeqGene_loci"
            "(refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos,"
            " orientation, totalLength, chrPos, rsgPos, entrezID, hgncSymbol, updated)"
            " VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, datetime('now'))"
        )
        conn = self.get_conn()
        cursor = conn.cursor()
        cursor.execute(query, (
            rsg_data[0], rsg_data[1], rsg_data[2],
            int(rsg_data[3]), int(rsg_data[4]),
            rsg_data[5], int(rsg_data[6]),
            rsg_data[7], rsg_data[8],
            int(rsg_data[9]), rsg_data[10],
        ))
        success = 'true' if cursor.lastrowid else 'Unknown error'
        conn.commit()
        cursor.close()
        conn.close()
        return success

    def update_refseq_gene_data(self, rsg_data):
        """UPDATE refSeqGene_loci hgncSymbol."""
        query = (
            "UPDATE refSeqGene_loci SET hgncSymbol=?, updated=datetime('now')"
            " WHERE refSeqGeneID=?"
        )
        conn = self.get_conn()
        cursor = conn.cursor()
        cursor.execute(query, (rsg_data[10], rsg_data[0]))
        conn.commit()
        cursor.close()
        conn.close()
        return 'true'

    # ------------------------------------------------------------------
    # LRG tables
    # ------------------------------------------------------------------

    def insert_refseq_gene_id_from_lrg_id(self, lrg_rs_lookup):
        """INSERT OR REPLACE into LRG_RSG_lookup."""
        query = (
            "INSERT OR REPLACE INTO LRG_RSG_lookup"
            "(lrgID, hgncSymbol, RefSeqGeneID, status)"
            " VALUES(?, ?, ?, ?)"
        )
        conn = self.get_conn()
        cursor = conn.cursor()
        cursor.execute(query, (lrg_rs_lookup[0], lrg_rs_lookup[1], lrg_rs_lookup[2], lrg_rs_lookup[3]))
        success = 'true' if cursor.lastrowid else 'Unknown error'
        conn.commit()
        cursor.close()
        conn.close()
        return success

    def insert_lrg_transcript_data(self, lrgtx_to_rst_id):
        """INSERT OR REPLACE into LRG_transcripts."""
        query = (
            "INSERT OR REPLACE INTO LRG_transcripts"
            "(LRGtranscriptID, RefSeqTranscriptID)"
            " VALUES(?, ?)"
        )
        conn = self.get_conn()
        cursor = conn.cursor()
        cursor.execute(query, (lrgtx_to_rst_id[0], lrgtx_to_rst_id[1]))
        success = 'true' if cursor.lastrowid else 'Unknown error'
        conn.commit()
        cursor.close()
        conn.close()
        return success

    def insert_lrg_protein_data(self, lrg_p, rs_p):
        """INSERT OR REPLACE into LRG_proteins."""
        query = (
            "INSERT OR REPLACE INTO LRG_proteins"
            "(LRGproteinID, RefSeqProteinID)"
            " VALUES(?, ?)"
        )
        conn = self.get_conn()
        cursor = conn.cursor()
        cursor.execute(query, (lrg_p, rs_p))
        success = 'true' if cursor.lastrowid else 'Unknown error'
        conn.commit()
        cursor.close()
        conn.close()
        return success

    # ------------------------------------------------------------------
    # stableGeneIds
    # ------------------------------------------------------------------

    def insert_gene_stable_ids(self, data):
        """INSERT OR REPLACE into stableGeneIds."""
        query = (
            "INSERT OR REPLACE INTO stableGeneIds"
            "(hgnc_id, hgnc_symbol, entrez_id, ensembl_gene_id, omim_id, ucsc_id, vega_id, ccds_ids)"
            " VALUES(?, ?, ?, ?, ?, ?, ?, ?)"
        )
        conn = self.get_conn()
        cursor = conn.cursor()
        cursor.execute(query, (
            data['hgnc_id'],
            data['hgnc_symbol'],
            data['entrez_id'],
            data['ensembl_gene_id'],
            data['omim_id'],
            data['ucsc_id'],
            data['vega_id'],
            data['ccds_id'],
        ))
        success = 'true' if cursor.lastrowid else 'unknown error'
        conn.commit()
        cursor.close()
        conn.close()
        return success

    def update_gene_stable_ids(self, gene_stable_ids):
        """UPDATE stableGeneIds."""
        query = (
            "UPDATE stableGeneIds SET hgnc_symbol=?, entrez_id=?, ensembl_gene_id=?,"
            " omim_id=?, ucsc_id=?, vega_id=?, ccds_ids=?"
            " WHERE hgnc_id=?"
        )
        conn = self.get_conn()
        cursor = conn.cursor()
        cursor.execute(query, (
            gene_stable_ids["hgnc_symbol"],
            gene_stable_ids["entrez_id"],
            gene_stable_ids["ensembl_gene_id"],
            gene_stable_ids["omim_id"],
            gene_stable_ids["ucsc_id"],
            gene_stable_ids["vega_id"],
            gene_stable_ids["ccds_id"],
            gene_stable_ids["hgnc_id"],
        ))
        conn.commit()
        cursor.close()
        conn.close()
        return 'true'

    # ------------------------------------------------------------------
    # version table
    # ------------------------------------------------------------------

    def update_db_version(self, db_version):
        """UPDATE version table."""
        query = "UPDATE version SET current_version=?"
        conn = self.get_conn()
        cursor = conn.cursor()
        cursor.execute(query, [db_version])
        conn.commit()
        cursor.close()
        conn.close()
        return 'true'


# <LICENSE>
# Copyright (C) 2016-2026 VariantValidator Contributors
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
