import pytest
from unittest import TestCase
from VariantValidator.modules.vvDatabase import Database
from VariantValidator import update_vv_db


class TestUpdate(TestCase):
    """
    I want these tests to run first, and in order so that I can make this testcase part of the preparation for
    running VV on travis.
    """

    def count_rows(self, db, table):
        query = "SELECT COUNT(*) FROM %s" % table
        row = db.execute(query)
        print(row)
        return row[0]

    def test_connection(self):
        db_conn = update_vv_db.connect()

        self.assertIsInstance(db_conn, Database)
        conn = db_conn.get_conn()
        try:
            self.assertTrue(conn.is_connected())
        except AttributeError:
            pass
            
    def test_deletion(self):
        db = update_vv_db.connect()
        initial_count = self.count_rows(db, 'LRG_transcripts')
        if initial_count > 0:
            pytest.skip("Already have data so not going to run this test here.")
        db.insert_refseq_gene_data(['id', 'chr', 'genome', '0', '1', '1', '10', '2', '2', '1', 'hgnc', 'False'])

        count = self.count_rows(db, 'refSeqGene_loci')
        self.assertGreaterEqual(count, 1)

        update_vv_db.delete()
        db = update_vv_db.connect()
        for table in ['transcript_info', 'refSeqGene_loci', 'LRG_transcripts', 'LRG_proteins', 'LRG_RSG_lookup']:
            print(table)
            count = self.count_rows(db, table)
            self.assertEqual(count, 0)

    def test_update(self):
        db = update_vv_db.connect()
        initial_count = self.count_rows(db, 'LRG_transcripts')
        if initial_count > 0:
            pytest.skip("Already have data so not going to run this test here.")

        update_vv_db.update()
        db = update_vv_db.connect()
        for table in ['refSeqGene_loci', 'LRG_transcripts', 'LRG_proteins', 'LRG_RSG_lookup']:
            print(table)
            count = self.count_rows(db, table)
            self.assertGreater(count, 0)

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
