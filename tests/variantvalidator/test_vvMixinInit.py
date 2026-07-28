import os
from unittest import TestCase
from unittest.mock import patch
from unittest.mock import MagicMock

from VariantValidator.modules.vvMixinInit import (
    Mixin,
    InitialisationError,
)
from VariantValidator import version


class TestVVMixinInit(TestCase):

    @patch("VariantValidator.modules.vvMixinInit.os.path.exists")
    @patch("VariantValidator.modules.vvMixinInit.settings.get_config_dir")
    def test_missing_configuration_file(
        self,
        mock_get_config_dir,
        mock_exists,
    ):
        mock_get_config_dir.return_value = "/tmp/config.ini"
        mock_exists.return_value = False

        with self.assertRaises(InitialisationError) as err:
            Mixin()

        self.assertIn(
            "Configuration file not found",
            str(err.exception),
        )
        self.assertIn(
            "/tmp/config.ini",
            str(err.exception),
        )

    def test_del_clears_pool(self):

        mixin = Mixin.__new__(Mixin)
        mixin.pool = object()

        mixin.__del__()

        self.assertIsNone(mixin.pool)

    def test_del_without_pool(self):

        mixin = Mixin.__new__(Mixin)

        mixin.__del__()

        self.assertIsNone(getattr(mixin, "pool", None))

    @patch("VariantValidator.modules.vvMixinInit.os.path.exists")
    @patch("VariantValidator.modules.vvMixinInit.settings.get_config_dir")
    @patch("VariantValidator.modules.vvMixinInit.ConfigParser")
    @patch("VariantValidator.modules.vvMixinInit.Database")
    def test_database_version_mismatch(
        self,
        mock_database,
        mock_configparser,
        mock_get_config_dir,
        mock_exists,
    ):
        mock_get_config_dir.return_value = "/tmp/config.ini"
        mock_exists.return_value = True

        config = mock_configparser.return_value
        config.read.return_value = None

        config.__getitem__.side_effect = lambda key: {
            "Entrez": {
                "email": "test@test.com",
                "api_key": "YOUR_API_KEY",
            },
            "seqrepo": {
                "version": "VV_SR",
                "require_threading": "True",
                "location": "/tmp",
            },
            "mysql": {
                "version": "vvdb_1",
                "user": "user",
                "password": "pass",
                "host": "localhost",
                "port": "3306",
                "database": "vvdb",
            },
            "postgres": {
                "user": "user",
                "password": "pass",
                "host": "localhost",
                "port": "5432",
                "database": "uta",
                "version": "uta_1",
            },
        }[key]

        config.get.return_value = False

        db = mock_database.return_value
        db.get_db_version.return_value = ("wrong_version",)

        with self.assertRaises(InitialisationError):
            Mixin()

    @patch("VariantValidator.modules.vvMixinInit.os.path.exists")
    @patch("VariantValidator.modules.vvMixinInit.settings.get_config_dir")
    @patch("VariantValidator.modules.vvMixinInit.ConfigParser")
    @patch("VariantValidator.modules.vvMixinInit.Database")
    @patch("VariantValidator.modules.vvMixinInit.vvhgvs")
    def test_check_same_thread_false_and_unix_socket(
        self,
        mock_vvhgvs,
        mock_database,
        mock_configparser,
        mock_get_config_dir,
        mock_exists,
    ):

        mock_get_config_dir.return_value = "/tmp/config.ini"
        mock_exists.return_value = True

        config = mock_configparser.return_value
        config.read.return_value = None

        config.__getitem__.side_effect = lambda key: {
            "Entrez": {
                "email": "test@test.com",
                "api_key": "YOUR_API_KEY",
            },
            "seqrepo": {
                "version": "VV_SR",
                "require_threading": "False",
                "location": "/tmp",
            },
            "mysql": {
                "version": "vvdb_1",
                "user": "user",
                "password": "pass",
                "host": "localhost",
                "port": "3306",
                "database": "vvdb",
            },
            "postgres": {
                "user": "user",
                "password": "pass",
                "host": "localhost",
                "port": "5432",
                "database": "uta",
                "version": "uta_1",
            },
        }[key]

        config.get.return_value = "/tmp/mysql.sock"

        db = mock_database.return_value
        db.get_db_version.return_value = ("vvdb_1",)

        mock_hdp = MagicMock()
        mock_hdp.data_version.return_value = "uta_1"

        mock_vvhgvs.__version__ = "4.0.0"
        mock_vvhgvs.dataproviders.uta.connect.return_value = mock_hdp
        mock_vvhgvs.parser.Parser.return_value = MagicMock()
        mock_vvhgvs.validator.Validator.return_value = MagicMock()
        mock_vvhgvs.variantmapper.VariantMapper.return_value = MagicMock()
        mock_vvhgvs.dataproviders.seqfetcher.SeqFetcher.return_value = MagicMock()
        mock_vvhgvs.normalizer.Normalizer.return_value = MagicMock()

        mixin = Mixin()

        self.assertTrue(mixin.check_same_thread)
        self.assertEqual(
            mixin.dbConfig["unix_socket"],
            "/tmp/mysql.sock",
        )

    @patch("VariantValidator.modules.vvMixinInit.os.path.exists")
    @patch("VariantValidator.modules.vvMixinInit.settings.get_config_dir")
    @patch("VariantValidator.modules.vvMixinInit.ConfigParser")
    @patch("VariantValidator.modules.vvMixinInit.Database")
    @patch("VariantValidator.modules.vvMixinInit.vvhgvs")
    def test_entrez_api_key_and_defaults(
        self,
        mock_vvhgvs,
        mock_database,
        mock_configparser,
        mock_get_config_dir,
        mock_exists,
    ):

        mock_get_config_dir.return_value = "/tmp/config.ini"
        mock_exists.return_value = True

        config = mock_configparser.return_value
        config.read.return_value = None

        config.__getitem__.side_effect = lambda key: {
            "Entrez": {
                "email": "test@test.com",
                "api_key": "abcdef123",
            },
            "seqrepo": {
                "version": "VV_SR",
                "require_threading": "True",
                "location": "/tmp",
            },
            "mysql": {
                "version": "vvdb_1",
                "user": "user",
                "password": "pass",
                "host": "localhost",
                "port": "3306",
                "database": "vvdb",
            },
            "postgres": {
                "user": "user",
                "password": "pass",
                "host": "localhost",
                "port": "5432",
                "database": "uta",
                "version": "uta_1",
            },
        }[key]

        config.get.return_value = False

        db = mock_database.return_value
        db.get_db_version.return_value = ("vvdb_1",)

        mock_hdp = MagicMock()
        mock_hdp.data_version.return_value = "uta_1"

        mock_vvhgvs.__version__ = "4.0.0"
        mock_vvhgvs.dataproviders.uta.connect.return_value = mock_hdp
        mock_vvhgvs.parser.Parser.return_value = MagicMock()
        mock_vvhgvs.validator.Validator.return_value = MagicMock()
        mock_vvhgvs.variantmapper.VariantMapper.return_value = MagicMock()
        mock_vvhgvs.dataproviders.seqfetcher.SeqFetcher.return_value = MagicMock()
        mock_vvhgvs.normalizer.Normalizer.return_value = MagicMock()

        mixin = Mixin()

        self.assertEqual(mixin.entrez_api_key, "abcdef123")
        self.assertFalse(mixin.check_same_thread)
        self.assertEqual(mixin.version, version.__version__)
        self.assertEqual(mixin.releasedVersion, version._is_released_version)
        self.assertEqual(mixin.primary_assembly, "GRCh38")
        self.assertFalse(mixin.testing)
        self.assertIsNone(mixin.selected_assembly)
        self.assertIsNone(mixin.select_transcripts)
        self.assertIsNone(mixin.alt_aln_method)

    @patch("VariantValidator.modules.vvMixinInit.os.path.exists")
    @patch("VariantValidator.modules.vvMixinInit.settings.get_config_dir")
    @patch("VariantValidator.modules.vvMixinInit.ConfigParser")
    @patch("VariantValidator.modules.vvMixinInit.Database")
    @patch("VariantValidator.modules.vvMixinInit.vvhgvs")
    def test_environment_variables(
        self,
        mock_vvhgvs,
        mock_database,
        mock_configparser,
        mock_get_config_dir,
        mock_exists,
    ):

        mock_get_config_dir.return_value = "/tmp/config.ini"
        mock_exists.return_value = True

        config = mock_configparser.return_value
        config.read.return_value = None

        config.__getitem__.side_effect = lambda key: {
            "Entrez": {
                "email": "test@test.com",
                "api_key": "YOUR_API_KEY",
            },
            "seqrepo": {
                "version": "VV_SR",
                "require_threading": "True",
                "location": "/tmp",
            },
            "mysql": {
                "version": "vvdb_1",
                "user": "user",
                "password": "pass",
                "host": "localhost",
                "port": "3306",
                "database": "vvdb",
            },
            "postgres": {
                "user": "user",
                "password": "pass",
                "host": "localhost",
                "port": "5432",
                "database": "uta",
                "version": "uta_1",
            },
        }[key]

        config.get.return_value = False

        db = mock_database.return_value
        db.get_db_version.return_value = ("vvdb_1",)

        mock_hdp = MagicMock()
        mock_hdp.data_version.return_value = "uta_1"

        mock_vvhgvs.__version__ = "4.0.0"
        mock_vvhgvs.dataproviders.uta.connect.return_value = mock_hdp
        mock_vvhgvs.parser.Parser.return_value = MagicMock()
        mock_vvhgvs.validator.Validator.return_value = MagicMock()
        mock_vvhgvs.variantmapper.VariantMapper.return_value = MagicMock()
        mock_vvhgvs.dataproviders.seqfetcher.SeqFetcher.return_value = MagicMock()
        mock_vvhgvs.normalizer.Normalizer.return_value = MagicMock()

        Mixin()

        self.assertEqual(
            os.environ["HGVS_SEQREPO_DIR"],
            "/tmp/VV_SR",
        )

        self.assertIn(
            "postgresql://user:pass@localhost:5432/uta/uta_1",
            os.environ["UTA_DB_URL"],
        )

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
