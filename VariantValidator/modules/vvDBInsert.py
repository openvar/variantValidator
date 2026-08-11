from .utils import handleCursor
from . import vvDBGet


class Mixin(
    vvDBGet.Mixin):
    """
    Object is a function container for inserting objects into the database.
    """

    @handleCursor  # Decorated function
    def insert(self,
               entry,
               data,
               table):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        if table == "transcript_info":
            query = (
                "INSERT INTO transcript_info("
                "refSeqID, description, transcriptVariant, currentVersion, "
                "hgncSymbol, utaSymbol, updated"
                ") VALUES (%s, %s, %s, %s, %s, %s, NOW())"
            )
            cursor.execute(
                query,
                (
                    entry,
                    data[1],
                    data[2],
                    data[3],
                    data[4],
                    data[5],
                ),
            )

        if cursor.lastrowid:
            success = "true"
        else: # Error
            success = "Unknown error"

        conn.commit()
        cursor.close()
        conn.close()

        return success # return

    @handleCursor  # Decorated function
    def insert_refseq_gene_data(self, rsg_data):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = (
            "INSERT INTO refSeqGene_loci("
            "refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos, "
            "orientation, totalLength, chrPos, rsgPos, entrezID, hgncSymbol, "
            "updated"
            ") VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, NOW())"
        )
        cursor.execute(
            query,
            (
                rsg_data[0],
                rsg_data[1],
                rsg_data[2],
                int(rsg_data[3]),
                int(rsg_data[4]),
                rsg_data[5],
                int(rsg_data[6]),
                rsg_data[7],
                rsg_data[8],
                int(rsg_data[9]),
                rsg_data[10],
            ),
        )

        if cursor.lastrowid:
            success = "true"
        else: # Error
            success = "Unknown error"

        conn.commit()
        cursor.close()
        conn.close()

        return success  # return

    @handleCursor  # Decorated function
    def insert_refseq_gene_id_from_lrg_id(self, lrg_rs_lookup):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = (
            "INSERT INTO LRG_RSG_lookup("
            "lrgID, hgncSymbol, RefSeqGeneID, status"
            ") VALUES (%s, %s, %s, %s)"
        )
        cursor.execute(
            query,
            (
                lrg_rs_lookup[0],
                lrg_rs_lookup[1],
                lrg_rs_lookup[2],
                lrg_rs_lookup[3],
            ),
        )

        if cursor.lastrowid:
            success = "true"
        else: # Error
            success = "Unknown error"

        conn.commit()
        cursor.close()
        conn.close()

        return success # return

    @handleCursor  # Decorated function
    def insert_lrg_transcript_data(self, lrgtx_to_rst_id):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = (
            "INSERT INTO LRG_transcripts("
            "LRGtranscriptID, RefSeqTranscriptID"
            ") VALUES (%s, %s)"
        )
        cursor.execute(
            query,
            (
                lrgtx_to_rst_id[0],
                lrgtx_to_rst_id[1],
            ),
        )

        if cursor.lastrowid:
            success = "true"
        else: # Errro
            success = "Unknown error"

        conn.commit()
        cursor.close()
        conn.close()

        return success # Return

    @handleCursor  # Decorated function
    def insert_lrg_protein_data(self, lrg_p, rs_p):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = (
            "INSERT INTO LRG_proteins("
            "LRGproteinID, RefSeqProteinID"
            ") VALUES (%s, %s)"
        )
        cursor.execute(
            query,
            (lrg_p, rs_p),
        )

        if cursor.lastrowid:
            success = "true"
        else: # Error
            success = "Unknown error"

        conn.commit()
        cursor.close()
        conn.close()

        return success # Return

    @handleCursor  # Decorated function
    def insert_gene_stable_ids(self, data):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = (
            "INSERT INTO stableGeneIds("
            "hgnc_id, hgnc_symbol, entrez_id, ensembl_gene_id, omim_id, "
            "ucsc_id, vega_id, ccds_ids"
            ") VALUES (%s, %s, %s, %s, %s, %s, %s, %s)"
        )
        cursor.execute(
            query,
            (
                data["hgnc_id"],
                data["hgnc_symbol"],
                data["entrez_id"],
                data["ensembl_gene_id"],
                data["omim_id"],
                data["ucsc_id"],
                data["vega_id"],
                data["ccds_id"],
            ),
        )

        if cursor.lastrowid:
            success = "true"
        else:
            success = "Unknown error"

        conn.commit()
        cursor.close()
        conn.close()

        return success

    @handleCursor  # Decorated function
    def update(self, entry, data):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = (
            "UPDATE transcript_info SET "
            "description=%s, transcriptVariant=%s, currentVersion=%s, "
            "hgncSymbol=%s, utaSymbol=%s, updated=NOW() "
            "WHERE refSeqID=%s"
        )
        cursor.execute(
            query,
            (
                data[1],
                data[2],
                data[3],
                data[4],
                data[5],
                entry,
            ),
        )

        conn.commit()
        cursor.close()
        conn.close()

        return "true"

    @handleCursor  # Decorated function
    def update_refseq_gene_data(self, rsg_data):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = (
            "UPDATE refSeqGene_loci SET "
            "hgncSymbol=%s, updated=NOW() "
            "WHERE refSeqGeneID=%s"
        )
        cursor.execute(
            query,
            (
                rsg_data[10],
                rsg_data[0],
            ),
        )

        conn.commit()
        cursor.close()
        conn.close()

        return "true"

    @handleCursor  # Decorated function
    def update_gene_stable_ids(self, gene_stable_ids):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = (
            "UPDATE stableGeneIds SET "
            "hgnc_symbol=%s, entrez_id=%s, ensembl_gene_id=%s, "
            "omim_id=%s, ucsc_id=%s, vega_id=%s, ccds_ids=%s "
            "WHERE hgnc_id=%s"
        )
        cursor.execute(
            query,
            (
                gene_stable_ids["hgnc_symbol"],
                gene_stable_ids["entrez_id"],
                gene_stable_ids["ensembl_gene_id"],
                gene_stable_ids["omim_id"],
                gene_stable_ids["ucsc_id"],
                gene_stable_ids["vega_id"],
                gene_stable_ids["ccds_id"],
                gene_stable_ids["hgnc_id"],
            ),
        )

        conn.commit()
        cursor.close()
        conn.close()

        return "true"

    @handleCursor  # Decorated function
    def update_db_version(self, db_version):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        query = "UPDATE version SET current_version=%s"
        cursor.execute(
            query,
            (db_version,),
        )

        conn.commit()
        cursor.close()
        conn.close()

        return "true"


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
