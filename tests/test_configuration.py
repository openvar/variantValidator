import unittest
import shutil
import os
import subprocess
import sys
import pytest
from configparser import ConfigParser
import tempfile

class TestConfigSetUp(unittest.TestCase):
    """
    Will test the configuration set up works as it should
    """

    @classmethod
    def setUpClass(cls):
        cls.filename = os.path.join(
            tempfile.gettempdir(),
            f".variantvalidator_{cls.__name__}"
        )

        os.environ["VARIANTVALIDATOR_TEST_CONFIG"] = cls.filename

        if os.path.exists(cls.filename):
            os.remove(cls.filename)

    @classmethod
    def tearDownClass(cls):
        os.environ.pop("VARIANTVALIDATOR_TEST_CONFIG", None)

        for f in (cls.filename, cls.filename + "_ori"):
            if os.path.exists(f):
                os.remove(f)

    def setUp(self):
        if os.path.exists(self.filename):
            shutil.move(self.filename, self.filename + "_ori")
        print("Moved file")

    def tearDown(self):
        backup = self.filename + "_ori"
        if os.path.exists(backup):
            shutil.move(backup, self.filename)

    def insert_blank(self):
        subprocess.check_output(
            ['python', '-c', 'import VariantValidator'],
            env=os.environ.copy()
        )

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
        output = subprocess.check_output(
            ['python', '-c', 'import VariantValidator'],
            env=os.environ.copy()
        )
        print(output)
        self.assertIn('Welcome to VariantValidator', output.decode())
        self.assertIn('Please edit this file', output.decode())

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
        output = subprocess.check_output(
            ['python', '-c', 'import VariantValidator'],
            env=os.environ.copy()
        )
        print(output)
        self.assertIn('MySQL', output.decode())
        self.assertIn('Please edit your configuration', output.decode())

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

        output = subprocess.check_output(
            ['python', '-c', 'import VariantValidator'],
            env=os.environ.copy()
        )

        print(output)

        self.assertIn('PostgreSQL', output.decode())
        self.assertIn('Please edit your configuration', output.decode())

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

        output = subprocess.check_output(
            ['python', '-c', 'import VariantValidator'],
            env=os.environ.copy()
        )

        self.assertIn('Seqrepo', output.decode())
        self.assertIn('Please edit your configuration', output.decode())

    def test_zz_changed_seqrepo(self):
        """
        Test is named as such so it runs last - as it will successfully import VariantValidator
        """
        self.insert_blank()
        self.open_config()

        self.config['mysql']['user'] = 'myusername'
        self.config['mysql']['password'] = 'mypass'

        self.config['postgres']['user'] = 'me'
        self.config['postgres']['password'] = 'pass'

        self.assertEqual(
            self.config['seqrepo']['location'],
            '/PATH/TO/SEQREPO'
        )

        self.config['seqrepo']['location'] = 'here'

        self.write_config()

        self.assertTrue(os.path.exists(self.filename))

        try:
            import VariantValidator
        except SystemExit:
            self.fail("SystemExit raised on Import")


class TestConfigValues(unittest.TestCase):
    """
    This class will test the config values that we're using, and that they are being read into VV correctly.
    """

    @classmethod
    def setUpClass(cls):
        cls.real_filename = os.path.join(
            os.path.expanduser('~'),
            '.variantvalidator'
        )

        cls.filename = os.path.join(
            tempfile.gettempdir(),
            f".variantvalidator_{cls.__name__}"
        )

        # Copy the user's working config to an isolated temporary location
        shutil.copy(cls.real_filename, cls.filename)

        os.environ["VARIANTVALIDATOR_TEST_CONFIG"] = cls.filename

    @classmethod
    def tearDownClass(cls):
        os.environ.pop("VARIANTVALIDATOR_TEST_CONFIG", None)

        for f in (cls.filename, cls.filename + "_ori"):
            if os.path.exists(f):
                os.remove(f)

    def setUp(self):
        self.original = self.filename + "_ori"
        shutil.copy(self.filename, self.original)

        self.config = ConfigParser()
        self.config.read(self.filename)

    def write_config(self):
        with open(self.filename, "w") as fh:
            self.config.write(fh)

    def test_file_structure(self):

        try:
            self.assertCountEqual(
                self.config.sections(),
                ['mysql', 'seqrepo', 'postgres', 'logging', 'Entrez']
            )
        except AssertionError:
            self.assertCountEqual(
                self.config.sections(),
                ['mysql', 'seqrepo', 'postgres', 'logging', 'Entrez', 'auth']
            )

        self.assertCountEqual(
            list(self.config['mysql']),
            ['host', 'port', 'database', 'user', 'password', 'version']
        )

        self.assertCountEqual(
            list(self.config['seqrepo']),
            ['version', 'location', 'require_threading']
        )

        self.assertCountEqual(
            list(self.config['postgres']),
            ['host', 'port', 'database', 'version', 'user', 'password']
        )

        self.assertCountEqual(
            list(self.config['logging']),
            ['log', 'console', 'file']
        )

        self.assertCountEqual(
            list(self.config['Entrez']),
            ['email', 'api_key']
        )

    def test_file_contents(self):
        self.assertNotEqual(self.config['mysql']['user'], 'USERNAME')
        self.assertNotEqual(self.config['mysql']['password'], 'PASSWORD')

        path = os.path.join(
            self.config['seqrepo']['location'],
            self.config['seqrepo']['version']
        )

        self.assertTrue(os.path.exists(path))

        self.assertNotEqual(self.config['postgres']['user'], 'USERNAME')
        self.assertNotEqual(self.config['postgres']['password'], 'PASSWORD')

        self.assertIsInstance(
            self.config['logging'].getboolean('log'),
            bool
        )

        self.assertIn(
            self.config['logging']['console'].upper(),
            ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG']
        )

        self.assertIn(
            self.config['logging']['file'].upper(),
            ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG']
        )

        if self.config['Entrez']['email'] != "OPTIONAL":
            self.assertRegex(
                self.config['Entrez']['email'],
                r'\w+@\w+\.\w+'
            )

    def test_file_parsing(self):
        import VariantValidator

        vv = VariantValidator.Validator()

        self.assertEqual(self.config['mysql']['user'], vv.dbConfig['user'])
        self.assertEqual(self.config['mysql']['password'], vv.dbConfig['password'])
        self.assertEqual(self.config['mysql']['host'], vv.dbConfig['host'])
        self.assertEqual(self.config['mysql']['database'], vv.dbConfig['database'])

        if 'unix_socket' in vv.dbConfig:
            self.assertEqual(
                self.config['mysql']['unix_socket'],
                vv.dbConfig['unix_socket']
            )

        self.assertEqual(
            vv.seqrepoPath,
            os.path.join(
                self.config['seqrepo']['location'],
                self.config['seqrepo']['version']
            )
        )

        host_or_socketfile = self.config['postgres']['host'].replace('/', '%2F')

        uta_path = "postgresql://%s:%s@%s:%s/%s/%s" % (
            self.config["postgres"]["user"],
            self.config["postgres"]["password"],
            host_or_socketfile,
            self.config['postgres']['port'],
            self.config['postgres']['database'],
            self.config['postgres']['version']
        )

        self.assertEqual(vv.utaPath, uta_path)

        self.assertEqual(
            vv.entrez_email,
            self.config['Entrez']['email']
        )

        if self.config['Entrez']['api_key'] == 'YOUR_API_KEY':
            self.assertEqual(vv.entrez_api_key, None)
        else:
            self.assertEqual(
                vv.entrez_api_key,
                self.config['Entrez']['api_key']
            )

    def tearDown(self):
        shutil.move(self.original, self.filename)

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
