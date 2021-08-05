import unittest
import shutil
import os
import subprocess
import sys
import pytest
from configparser import ConfigParser


class TestConfigSetUp(unittest.TestCase):
    """
    Will test the configuration set up works as it should
    """

    @classmethod
    def setUpClass(cls):
        cls.filename = os.path.join(os.path.expanduser('~'), '.variantvalidator')

    def setUp(self):
        new_filename = self.filename + '_ori'
        shutil.move(self.filename, new_filename)
        print('Moved file')

    def insert_blank(self):
        subprocess.check_output(['python', '-c', 'import VariantValidator'])

    def open_config(self):
        self.config = ConfigParser()
        self.config.read(self.filename)

    def write_config(self):
        with open(self.filename, 'w') as fh:
            self.config.write(fh)

    def test_no_config_file(self):
        if 'VariantValidator' in list(sys.modules.keys()):
            pytest.skip("VariantValidator already imported")
        self.assertFalse(os.path.exists(self.filename))
        with self.assertRaises(SystemExit):
            import VariantValidator

    def test_no_config_file_msg(self):
        self.assertFalse(os.path.exists(self.filename))
        output = subprocess.check_output(['python', '-c', 'import VariantValidator'])
        print(output)
        self.assertTrue('Welcome to VariantValidator' in output.decode())
        self.assertTrue('Please edit this file' in output.decode())

    def test_unchanged_file(self):
        if 'VariantValidator' in list(sys.modules.keys()):
            pytest.skip("VariantValidator already imported")
        self.assertFalse(os.path.exists(self.filename))
        self.insert_blank()
        self.assertTrue(os.path.exists(self.filename))
        with self.assertRaises(SystemExit):
            import VariantValidator

    def test_unchanged_file_msg(self):
        self.insert_blank()
        self.assertTrue(os.path.exists(self.filename))
        output = subprocess.check_output(['python', '-c', 'import VariantValidator'])
        print(output)
        self.assertTrue('MySQL' in output.decode())
        self.assertTrue('Please edit your configuration' in output.decode())

    def test_changed_mysql(self):
        if 'VariantValidator' in list(sys.modules.keys()):
            pytest.skip("VariantValidator already imported")
        self.insert_blank()
        self.open_config()
        self.assertEqual(self.config['mysql']['user'], 'USERNAME')
        self.config['mysql']['user'] = 'myusername'
        self.assertEqual(self.config['mysql']['password'], 'PASSWORD')
        self.config['mysql']['password'] = 'mypass'
        self.write_config()

        self.assertTrue(os.path.exists(self.filename))
        with self.assertRaises(SystemExit):
            import VariantValidator

    def test_changed_mysql_msg(self):
        self.insert_blank()
        self.open_config()
        self.assertEqual(self.config['mysql']['user'], 'USERNAME')
        self.config['mysql']['user'] = 'myusername'
        self.assertEqual(self.config['mysql']['password'], 'PASSWORD')
        self.config['mysql']['password'] = 'mypass'
        self.write_config()

        output = subprocess.check_output(['python', '-c', 'import VariantValidator'])
        print(output)
        self.assertTrue('PostgreSQL' in output.decode())
        self.assertTrue('Please edit your configuration' in output.decode())

    def test_changed_postgres(self):
        if 'VariantValidator' in list(sys.modules.keys()):
            pytest.skip("VariantValidator already imported")
        self.insert_blank()
        self.open_config()
        self.config['mysql']['user'] = 'myusername'
        self.config['mysql']['password'] = 'mypass'

        self.assertEqual(self.config['postgres']['user'], 'USERNAME')
        self.assertEqual(self.config['postgres']['password'], 'PASSWORD')
        self.config['postgres']['user'] = 'me'
        self.config['postgres']['password'] = 'pass'
        self.write_config()

        self.assertTrue(os.path.exists(self.filename))
        with self.assertRaises(SystemExit):
            import VariantValidator

    def test_changed_postgres_msg(self):
        self.insert_blank()
        self.open_config()
        self.config['mysql']['user'] = 'myusername'
        self.config['mysql']['password'] = 'mypass'

        self.assertEqual(self.config['postgres']['user'], 'USERNAME')
        self.assertEqual(self.config['postgres']['password'], 'PASSWORD')
        self.config['postgres']['user'] = 'me'
        self.config['postgres']['password'] = 'pass'
        self.write_config()

        output = subprocess.check_output(['python', '-c', 'import VariantValidator'])
        self.assertTrue('Seqrepo' in output.decode())
        self.assertTrue('Please edit your configuration' in output.decode())

    def test_zz_changed_seqrepo(self):
        """
        Test is named as such so it runs last - as it will successfully import VariantValidator
        :return:
        """
        self.insert_blank()
        self.open_config()
        self.config['mysql']['user'] = 'myusername'
        self.config['mysql']['password'] = 'mypass'
        self.config['postgres']['user'] = 'me'
        self.config['postgres']['password'] = 'pass'

        self.assertEqual(self.config['seqrepo']['location'], '/PATH/TO/SEQREPO')
        self.config['seqrepo']['location'] = 'here'
        self.write_config()

        self.assertTrue(os.path.exists(self.filename))
        try:
            import VariantValidator
        except SystemExit:
            self.fail('SystemExit raised on Import')

    def tearDown(self):
        original = os.path.join(os.path.expanduser('~'), '.variantvalidator')
        new_filename = original + '_ori'
        shutil.move(new_filename, original)
        print('Moved file back')


class TestConfigValues(unittest.TestCase):
    """
    This class will test the config values that we're using, and that they are being read into VV correctly.
    """

    @classmethod
    def setUpClass(cls):
        cls.filename = os.path.join(os.path.expanduser('~'), '.variantvalidator')

    def setUp(self):
        self.original = self.filename + '_ori'
        shutil.copy(self.filename, self.original)
        config = ConfigParser()
        config.read(self.filename)
        self.config = config

    def write_config(self):
        with open(self.filename, 'w') as fh:
            self.config.write(fh)

    def test_file_structure(self):
        self.assertEqual(self.config.sections(), ['mysql', 'seqrepo', 'postgres',  'logging', 'Entrez'])
        self.assertEqual(list(self.config['mysql']), ['host', 'port', 'database', 'user', 'password', 'version'])
        self.assertEqual(list(self.config['seqrepo']), ['version', 'location'])
        self.assertEqual(list(self.config['postgres']), ['host', 'port', 'database', 'version', 'user', 'password'])
        self.assertEqual(list(self.config['logging']), ['log', 'console', 'file'])
        self.assertEqual(list(self.config['Entrez']), ['email', 'api_key'])

    def test_file_contents(self):
        self.assertNotEqual(self.config['mysql']['user'], 'USERNAME')
        self.assertNotEqual(self.config['mysql']['password'], 'PASSWORD')

        #self.assertEqual(self.config['seqrepo']['version'], '2018-08-21')
        path = os.path.join(self.config['seqrepo']['location'], self.config['seqrepo']['version'])
        self.assertTrue(os.path.exists(path))

        self.assertEqual(self.config['postgres']['version'], 'vvta_2021_2')
        self.assertNotEqual(self.config['postgres']['user'], 'USERNAME')
        self.assertNotEqual(self.config['postgres']['password'], 'PASSWORD')

        self.assertIsInstance(self.config['logging'].getboolean('log'), bool)
        self.assertIn(self.config['logging']['console'].upper(), ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG'])
        self.assertIn(self.config['logging']['file'].upper(), ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG'])

        self.assertRegex(self.config['Entrez']['email'], r'\w+@\w+.\w+')

    def test_file_parsing(self):
        import VariantValidator

        vv = VariantValidator.Validator()

        self.assertEqual(self.config['mysql']['user'], vv.dbConfig['user'])
        self.assertEqual(self.config['mysql']['password'], vv.dbConfig['password'])
        self.assertEqual(self.config['mysql']['host'], vv.dbConfig['host'])
        self.assertEqual(self.config['mysql']['database'], vv.dbConfig['database'])

        self.assertEqual(vv.seqrepoPath,
                         os.path.join(self.config['seqrepo']['location'], self.config['seqrepo']['version']))

        self.assertEqual(vv.utaPath, "postgresql://%s:%s@%s:%s/%s/%s" % (
            self.config["postgres"]["user"],
            self.config["postgres"]["password"],
            self.config['postgres']['host'],
            self.config['postgres']['port'],
            self.config['postgres']['database'],
            self.config['postgres']['version']
        ))

        self.assertEqual(vv.entrez_email, self.config['Entrez']['email'])
        if self.config['Entrez']['api_key'] == 'YOUR_API_KEY':
            self.assertEqual(vv.entrez_api_key, None)
        else:
            self.assertEqual(vv.entrez_api_key, self.config['Entrez']['api_key'])

    def tearDown(self):
        shutil.move(self.original, self.filename)

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
