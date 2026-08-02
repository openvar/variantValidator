import logging
from collections import OrderedDict

from VariantValidator import settings

from .utils import handleCursor
from . import vvDBInit


logger = logging.getLogger(__name__)

LRG_TX_LINK = {}

DB_GET_CACHE = OrderedDict()
_CACHE_MISS = object()


def _get_cached(key):
    if not settings.vvDB_GET_CACHE:
        return _CACHE_MISS

    try:
        value = DB_GET_CACHE.pop(key)
    except KeyError:
        return _CACHE_MISS

    DB_GET_CACHE[key] = value
    return value


def _set_cached(key, value):
    if not settings.vvDB_GET_CACHE:
        return value

    DB_GET_CACHE[key] = value
    DB_GET_CACHE.move_to_end(key)

    while len(DB_GET_CACHE) > settings.vvDB_GET_CACHE_SIZE:
        DB_GET_CACHE.popitem(last=False)

    return value


class Mixin(vvDBInit.Mixin):
    """
    Most of the functions in DBGet generate queries for retrieving data
    from the databases.
    """

    @handleCursor
    def execute(self, *query_args):
        attempts = 3

        for attempt in range(attempts):
            conn = self.get_conn()
            cursor = self.get_cursor(conn)

            try:
                cursor.execute(*query_args)
                row = cursor.fetchone()

                if row is None:
                    logger.debug(
                        "No data returned from query %s",
                        query_args,
                    )
                    row = ["none", "No data"]

                return row

            except Exception as e:
                logger.info(
                    "MySQL error (attempt %s/%s): %s",
                    attempt + 1,
                    attempts,
                    e,
                )

                if attempt < attempts - 1:
                    try:
                        conn.reconnect(
                            attempts=1,
                            delay=0,
                        )
                    except Exception:
                        pass
                else:
                    raise

            finally:
                try:
                    cursor.close()
                    conn.close()
                except Exception:
                    pass

    @handleCursor
    def execute_write(self, *query_args):
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        try:
            cursor.execute(*query_args)
            conn.commit()
        finally:
            try:
                cursor.close()
                conn.close()
            except Exception:
                pass

    @handleCursor
    def execute_all(self, *query_args):
        attempts = 3

        for attempt in range(attempts):
            conn = self.get_conn()
            cursor = self.get_cursor(conn)

            try:
                cursor.execute(*query_args)
                rows = cursor.fetchall()

                if not rows:
                    logger.debug(
                        "No data returned from query %s",
                        query_args,
                    )
                    rows = [["none", "No data"]]

                return rows

            except Exception as e:
                logger.info(
                    "MySQL error (attempt %s/%s): %s",
                    attempt + 1,
                    attempts,
                    e,
                )

                if attempt < attempts - 1:
                    try:
                        conn.reconnect(
                            attempts=1,
                            delay=0,
                        )
                    except Exception:
                        pass
                else:
                    raise

            finally:
                try:
                    cursor.close()
                    conn.close()
                except Exception:
                    pass

    # From dbfetchone.
    def get_uta(self, gene_symbol):
        query = (
            "SELECT utaSymbol FROM transcript_info "
            "WHERE hgncSymbol = %s"
        )
        return self.execute(
            query,
            (gene_symbol,),
        )

    def get_hgnc(self, gene_symbol):
        query = (
            "SELECT hgncSymbol FROM transcript_info "
            "WHERE utaSymbol = %s"
        )
        return self.execute(
            query,
            (gene_symbol,),
        )

    def get_transcript_description(self, transcript_id):
        key = ("transcript_description", transcript_id)
        cached = _get_cached(key)

        if cached is not _CACHE_MISS:
            return cached

        query = (
            "SELECT description FROM transcript_info "
            "WHERE refSeqID = %s"
        )
        result = str(
            self.execute(
                query,
                (transcript_id,),
            )[0]
        )

        return _set_cached(key, result)

    def get_transcript_annotation(self, transcript_id):
        key = ("transcript_annotation", transcript_id)
        cached = _get_cached(key)

        if cached is not _CACHE_MISS:
            return cached

        query = (
            "SELECT transcriptVariant FROM transcript_info "
            "WHERE refSeqID = %s"
        )
        result = str(
            self.execute(
                query,
                (transcript_id,),
            )[0]
        )

        return _set_cached(key, result)

    def get_gene_symbol_from_transcript_id(self, transcript_id):
        key = ("transcript_gene_symbol", transcript_id)
        cached = _get_cached(key)

        if cached is not _CACHE_MISS:
            return cached

        query = (
            "SELECT hgncSymbol FROM transcript_info "
            "WHERE refSeqID = %s"
        )
        result = str(
            self.execute(
                query,
                (transcript_id,),
            )[0]
        )

        return _set_cached(key, result)

    def get_refseq_data_by_refseq_id(
            self,
            refseq_id,
            genome_build,
    ):
        query = (
            "SELECT refSeqGeneID, refSeqChromosomeID, genomeBuild, "
            "startPos, endPos, orientation, totalLength, chrPos, "
            "rsgPos, entrezID, hgncSymbol "
            "FROM refSeqGene_loci "
            "WHERE refSeqGeneID = %s "
            "AND genomeBuild = %s"
        )
        return self.execute(
            query,
            (refseq_id, genome_build),
        )

    def get_gene_symbol_from_refseq_id(self, refseq_id):
        key = ("gene_symbol_from_refseq_id", refseq_id)
        cached = _get_cached(key)

        if cached is not _CACHE_MISS:
            return cached

        query = (
            "SELECT hgncSymbol FROM refSeqGene_loci "
            "WHERE refSeqGeneID = %s"
        )
        result = self.execute(
            query,
            (refseq_id,),
        )[0]

        return _set_cached(key, result)

    def get_refseq_id_from_lrg_id(self, lrg_id):
        query = (
            "SELECT RefSeqGeneID FROM LRG_RSG_lookup "
            "WHERE lrgID = %s"
        )
        return self.execute(
            query,
            (lrg_id,),
        )[0]

    def get_refseq_transcript_id_from_lrg_transcript_id(
            self,
            lrg_tx_id,
    ):
        query = (
            "SELECT RefSeqTranscriptID FROM LRG_transcripts "
            "WHERE LRGtranscriptID = %s"
        )
        return self.execute(
            query,
            (lrg_tx_id,),
        )[0]

    def get_lrg_transcript_id_from_refseq_transcript_id(
            self,
            rst_id,
    ):
        if not LRG_TX_LINK:
            query = (
                "SELECT RefSeqTranscriptID, LRGtranscriptID "
                "FROM LRG_transcripts"
            )
            lrg_data = self.execute_all(query)

            for refseq_id, lrg_id in lrg_data:
                LRG_TX_LINK[refseq_id] = lrg_id

        return LRG_TX_LINK.get(
            rst_id,
            "none",
        )

    def get_lrg_id_from_refseq_gene_id(self, rsg_id):
        key = ("lrg_id_from_refseq_gene_id", rsg_id)
        cached = _get_cached(key)

        if cached is not _CACHE_MISS:
            return cached

        query = (
            "SELECT lrgID, status FROM LRG_RSG_lookup "
            "WHERE RefSeqGeneID = %s"
        )
        result = self.execute(
            query,
            (rsg_id,),
        )

        return _set_cached(key, result)

    def get_refseqgene_info(
            self,
            refseqgene_id,
            primary_assembly,
    ):
        query = (
            "SELECT refSeqGeneID, refSeqChromosomeID, genomeBuild, "
            "startPos, endPos "
            "FROM refSeqGene_loci "
            "WHERE refSeqGeneID = %s "
            "AND genomeBuild = %s"
        )
        return self.execute(
            query,
            (refseqgene_id, primary_assembly),
        )

    def get_refseq_protein_id_from_lrg_protein_id(self, lrg_p):
        query = (
            "SELECT RefSeqProteinID FROM LRG_proteins "
            "WHERE LRGproteinID = %s"
        )
        return self.execute(
            query,
            (lrg_p,),
        )[0]

    def get_lrg_protein_id_from_ref_seq_protein_id(self, rs_p):
        key = ("lrg_protein_id_from_refseq_protein_id", rs_p)
        cached = _get_cached(key)

        if cached is not _CACHE_MISS:
            return cached

        query = (
            "SELECT LRGproteinID FROM LRG_proteins "
            "WHERE RefSeqProteinID = %s"
        )
        result = self.execute(
            query,
            (rs_p,),
        )[0]

        return _set_cached(key, result)

    def get_lrg_data_from_lrg_id(self, lrg_id):
        query = (
            "SELECT * FROM LRG_RSG_lookup "
            "WHERE lrgID = %s"
        )
        return self.execute(
            query,
            (lrg_id,),
        )

    def get_transcript_info_for_gene(self, gene_symbol):
        query = (
            "SELECT refSeqID, description, transcriptVariant, "
            "currentVersion, hgncSymbol, utaSymbol, updated, "
            "IF(updated < NOW() - INTERVAL 3 MONTH, 'true', 'false') "
            "FROM transcript_info "
            "WHERE hgncSymbol = %s"
        )
        return self.execute_all(
            query,
            (gene_symbol,),
        )

    def get_g_to_g_info(
            self,
            rsg_id=None,
            gen_id=None,
            start=None,
            end=None,
    ):
        """
        Return a set of RSG to genome mapping data.

        Can be all such mappings or can be limited to either data for a
        specific RSG or for a specific genomic reference ID, with an
        optional location.
        """
        query = (
            "SELECT refSeqGeneID, refSeqChromosomeID, startPos, "
            "endPos, orientation, hgncSymbol, genomeBuild "
            "FROM refSeqGene_loci"
        )
        query_vals = ()

        if rsg_id:
            query += " WHERE refSeqGeneID = %s"
            query_vals = (rsg_id,)

        elif gen_id:
            query += " WHERE refSeqChromosomeID = %s"
            query_vals = (gen_id,)

            if start is not None:
                query += " AND startPos <= %s"
                query_vals += (start,)

            if end is not None:
                query += " AND endPos >= %s"
                query_vals += (end,)

        return self.execute_all(
            query,
            query_vals,
        )

    def get_all_transcript_id(self):
        query = "SELECT refSeqID FROM transcript_info"
        return self.execute_all(query)

    def get_stable_gene_id_info(self, hgnc_symbol):
        key = ("stable_gene_id_info", hgnc_symbol)
        cached = _get_cached(key)

        if cached is not _CACHE_MISS:
            return cached

        query = (
            "SELECT * FROM stableGeneIds "
            "WHERE hgnc_symbol = %s"
        )
        result = self.execute(
            query,
            (hgnc_symbol,),
        )

        return _set_cached(key, result)

    def get_stable_gene_id_from_hgnc_id(self, hgnc_id):
        key = ("stable_gene_id_from_hgnc_id", hgnc_id)
        cached = _get_cached(key)

        if cached is not _CACHE_MISS:
            return cached

        query = (
            "SELECT * FROM stableGeneIds "
            "WHERE hgnc_id = %s"
        )
        result = self.execute(
            query,
            (hgnc_id,),
        )

        return _set_cached(key, result)

    def get_transcripts_from_annotations(self, statement):
        query = (
            "SELECT * FROM transcript_info "
            "WHERE transcriptVariant LIKE %s"
        )
        return self.execute_all(
            query,
            (f"%{statement}%",),
        )

    def get_db_version(self):
        """
        Return the current version of the Validator database.
        """
        query = "SELECT current_version FROM version"
        return self.execute(query)

    # Direct methods (GET).
    def get_uta_symbol(self, gene_symbol):
        # Return the UTA gene symbol when an HGNC gene symbol is input.
        return str(
            self.get_uta(gene_symbol)[0]
        )

    def get_hgnc_symbol(self, gene_symbol):
        # Return the HGNC gene symbol when a UTA gene symbol is input.
        return str(
            self.get_hgnc(gene_symbol)[0]
        )

    # From external.py.
    def get_urls(self, dict_out):
        """
        Provide direct links to reference sequence records.
        """
        report_urls = {}

        transcript_variant = dict_out["hgvs_transcript_variant"]
        transcript_accession = transcript_variant.split(":", 1)[0]

        protein_consequence = str(
            dict_out["hgvs_predicted_protein_consequence"]["slr"]
        )
        protein_accession = protein_consequence.split(":", 1)[0]

        refseqgene_variant = dict_out["hgvs_refseqgene_variant"]
        refseqgene_accession = refseqgene_variant.split(":", 1)[0]

        lrg_variant = dict_out["hgvs_lrg_variant"]
        lrg_id = lrg_variant.split(":", 1)[0]

        selected_assembly = str(
            dict_out["selected_assembly"]
        ).lower()

        # RefSeq.
        if transcript_accession.startswith(("NM_", "NR_")):
            report_urls["transcript"] = (
                "https://www.ncbi.nlm.nih.gov/"
                f"nuccore/{transcript_accession}"
            )

        if protein_accession.startswith("NP_"):
            report_urls["protein"] = (
                "https://www.ncbi.nlm.nih.gov/"
                f"nuccore/{protein_accession}"
            )

        if refseqgene_accession.startswith("NG_"):
            report_urls["refseqgene"] = (
                "https://www.ncbi.nlm.nih.gov/"
                f"nuccore/{refseqgene_accession}"
            )

        if lrg_id.startswith("LRG"):
            lrg_data = self.get_lrg_data_from_lrg_id(
                lrg_id
            )

            # LRG identifiers may not have a corresponding lookup
            # record in the database. In this case execute() returns
            # ["none", "No data"], so only attempt to determine the
            # publication status when a valid row has been returned.
            if (
                lrg_data
                and lrg_data[0] != "none"
                and len(lrg_data) > 4
            ):
                lrg_status = str(lrg_data[4])

                if lrg_status == "public":
                    report_urls["lrg"] = (
                        "http://ftp.ebi.ac.uk/pub/"
                        f"databases/lrgex/{lrg_id}.xml"
                    )
                else:
                    report_urls["lrg"] = (
                        "http://ftp.ebi.ac.uk/pub/"
                        "databases/lrgex/pending/"
                        f"{lrg_id}.xml"
                    )

        # Ensembl.
        if selected_assembly == "grch37":
            if transcript_accession.startswith("ENST"):
                report_urls["transcript"] = (
                    "https://grch37.ensembl.org/"
                    "Homo_sapiens/Transcript/Summary?"
                    f"db=core;t={transcript_accession}"
                )

            if protein_accession.startswith("ENSP"):
                report_urls["protein"] = (
                    "https://grch37.ensembl.org/"
                    "Homo_sapiens/Transcript/ProteinSummary?"
                    f"db=core;p={protein_accession}"
                )

        elif selected_assembly == "grch38":
            if transcript_accession.startswith("ENST"):
                report_urls["transcript"] = (
                    "https://www.ensembl.org/"
                    "Homo_sapiens/Transcript/Summary?"
                    f"db=core;t={transcript_accession}"
                )

            if protein_accession.startswith("ENSP"):
                report_urls["protein"] = (
                    "https://www.ensembl.org/"
                    "Homo_sapiens/Transcript/ProteinSummary?"
                    f"db=core;p={protein_accession}"
                )

        return report_urls


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
