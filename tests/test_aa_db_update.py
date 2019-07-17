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
        self.assertTrue(db_conn.conn.is_connected())

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
        for table in ['refSeqGene_loci', 'LRG_transcripts', 'LRG_proteins', 'LRG_RSG_lookup']:
            print(table)
            count = self.count_rows(db, table)
            self.assertGreater(count, 0)


