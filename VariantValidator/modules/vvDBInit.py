import random
try:
    import mariadb
except ModuleNotFoundError:
    try:
        import mysql.connector
        from mysql.connector.pooling import MySQLConnectionPool
    except ModuleNotFoundError:
        mariadb = None
        mysql = None


class Mixin:
    """
    A mixin containing the database initialisation routines.
    """
    def __init__(self, db_config):
        self.pool = None
        self.dbConfig = db_config
        self.init_db()

    def __del__(self):
        if self.pool:
            self.pool = None

    def init_db(self):
        try:
            self.pool = mysql.connector.pooling.MySQLConnectionPool(pool_size=5, **self.dbConfig)
        except NameError:
            ran_num = random.random()
            name_my_pool = 'pool%s' % str(ran_num)
            try:
                self.pool = mariadb.ConnectionPool(pool_size=5,
                                                   pool_name=name_my_pool,
                                                   pool_reset_connection=False,
                                                   host=self.dbConfig['host'],
                                                   user=self.dbConfig['user'],
                                                   port=int(self.dbConfig['port']),
                                                   password=self.dbConfig['password'],
                                                   database=self.dbConfig['database']
                                                   )
            except mariadb.ProgrammingError:
                ran_num = random.random()
                name_my_pool = 'pool%s' % str(ran_num)
                self.pool = mariadb.ConnectionPool(pool_size=5,
                                                   pool_name=name_my_pool,
                                                   pool_reset_connection=False,
                                                   host=self.dbConfig['host'],
                                                   user=self.dbConfig['user'],
                                                   port=int(self.dbConfig['port']),
                                                   password=self.dbConfig['password'],
                                                   database=self.dbConfig['database']
                                                   )

    def get_conn(self):
        try:
            conn = self.pool.get_connection()
        except Exception:
            self.init_db()
            conn = self.pool.get_connection()
        return conn

    def get_cursor(self, conn):
        try:
            cursor = conn.cursor(buffered=True)
        except Exception:
            self.init_db()
            self.get_conn()
            cursor = conn.cursor(buffered=True)
        return cursor

# ---------------------------------------------------------------------------
# SQLite backend (new)
# ---------------------------------------------------------------------------

import sqlite3 as _sqlite3

_SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS transcript_info (
    refSeqID        TEXT PRIMARY KEY,
    description     TEXT,
    transcriptVariant TEXT,
    currentVersion  TEXT,
    hgncSymbol      TEXT,
    utaSymbol       TEXT,
    updated         TEXT
);

CREATE TABLE IF NOT EXISTS refSeqGene_loci (
    refSeqGeneID        TEXT,
    refSeqChromosomeID  TEXT,
    genomeBuild         TEXT,
    startPos            INTEGER,
    endPos              INTEGER,
    orientation         TEXT,
    totalLength         INTEGER,
    chrPos              TEXT,
    rsgPos              TEXT,
    entrezID            INTEGER,
    hgncSymbol          TEXT,
    updated             TEXT,
    PRIMARY KEY (refSeqGeneID, genomeBuild)
);

CREATE TABLE IF NOT EXISTS LRG_RSG_lookup (
    lrgID           TEXT PRIMARY KEY,
    hgncSymbol      TEXT,
    RefSeqGeneID    TEXT,
    status          TEXT
);

CREATE TABLE IF NOT EXISTS LRG_transcripts (
    LRGtranscriptID     TEXT PRIMARY KEY,
    RefSeqTranscriptID  TEXT
);

CREATE TABLE IF NOT EXISTS LRG_proteins (
    LRGproteinID    TEXT PRIMARY KEY,
    RefSeqProteinID TEXT
);

CREATE TABLE IF NOT EXISTS stableGeneIds (
    hgnc_id         TEXT PRIMARY KEY,
    hgnc_symbol     TEXT,
    entrez_id       TEXT,
    ensembl_gene_id TEXT,
    omim_id         TEXT,
    ucsc_id         TEXT,
    vega_id         TEXT,
    ccds_ids        TEXT
);

CREATE TABLE IF NOT EXISTS version (
    current_version TEXT
);
"""


class SQLiteDBInit:
    """SQLite connection — file-based, no server required."""

    def __init__(self, db_path: str):
        self.db_path = db_path
        self._init_schema()

    def _init_schema(self):
        conn = _sqlite3.connect(self.db_path, check_same_thread=False)
        conn.executescript(_SCHEMA_SQL)
        conn.commit()
        conn.close()

    def get_conn(self):
        return _sqlite3.connect(self.db_path, check_same_thread=False)

    def get_cursor(self, conn):
        return conn.cursor()


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
