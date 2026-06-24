import pytest
from unittest.mock import MagicMock, patch

from VariantValidator.modules.vvDBInit import Mixin


@pytest.fixture
def db_config():
    return {
        "host": "localhost",
        "user": "user",
        "port": "3306",
        "password": "password",
        "database": "vv",
    }


def test_del_clears_pool():
    obj = Mixin.__new__(Mixin)
    obj.pool = object()

    obj.__del__()

    assert obj.pool is None


def test_del_with_no_pool():
    obj = Mixin.__new__(Mixin)
    obj.pool = None

    obj.__del__()

    assert obj.pool is None


@patch("VariantValidator.modules.vvDBInit.mysql.connector.pooling.MySQLConnectionPool")
def test_init_db_mysql_pool(mock_pool, db_config):
    obj = Mixin.__new__(Mixin)
    obj.dbConfig = db_config

    obj.init_db()

    mock_pool.assert_called_once_with(
        pool_size=5,
        **db_config
    )


def test_get_conn_success():
    obj = Mixin.__new__(Mixin)

    conn = MagicMock()
    pool = MagicMock()
    pool.get_connection.return_value = conn

    obj.pool = pool

    result = obj.get_conn()

    assert result is conn


def test_get_conn_reconnect():
    obj = Mixin.__new__(Mixin)

    conn = MagicMock()

    pool = MagicMock()
    pool.get_connection.side_effect = [
        Exception("boom"),
        conn,
    ]

    obj.pool = pool
    obj.init_db = MagicMock()

    result = obj.get_conn()

    obj.init_db.assert_called_once()
    assert result is conn


def test_get_cursor_success():
    obj = Mixin.__new__(Mixin)

    cursor = MagicMock()
    conn = MagicMock()
    conn.cursor.return_value = cursor

    result = obj.get_cursor(conn)

    assert result is cursor


def test_get_cursor_reconnect():
    obj = Mixin.__new__(Mixin)

    cursor = MagicMock()

    bad_conn = MagicMock()
    bad_conn.cursor.side_effect = Exception("boom")

    good_conn = MagicMock()
    good_conn.cursor.return_value = cursor

    obj.init_db = MagicMock()
    obj.get_conn = MagicMock(return_value=good_conn)

    result = obj.get_cursor(bad_conn)

    obj.init_db.assert_called_once()
    obj.get_conn.assert_called_once()
    assert result is cursor


@patch("VariantValidator.modules.vvDBInit.Mixin.init_db")
def test_init_called_from_constructor(mock_init, db_config):
    obj = Mixin(db_config)

    assert obj.dbConfig == db_config
    mock_init.assert_called_once()

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