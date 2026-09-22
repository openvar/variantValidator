"""
Tests for database connection configuration in VariantValidator MixinInit.

These tests exercise the real VariantValidator database initialisation code
and the real vvhgvs UTA PostgreSQL dataprovider.

Actual database servers are replaced only at the database-driver boundary:

    MixinInit
        -> vvDBInit.Database
        -> MySQLConnectionPool [fake driver]
        -> UTA_postgresql / uta.connect()
        -> psycopg2 ThreadedConnectionPool [fake driver]

The tests deliberately do not patch the connection-building logic in
MixinInit, vvhgvs, or elastic_pool.py.
"""

import os
from configparser import ConfigParser
from pathlib import Path
from unittest.mock import MagicMock

import pytest

import vvhgvs
from VariantValidator import settings
from VariantValidator.modules.vvDBInit import Mixin
from VariantValidator.modules.vvMixinInit import Mixin as MixinInit


POSTGRES_HOST = "10.73.203.5"
POSTGRES_INSTANCE = "tcshaip:europe-west1:vv-vvta-replica"
POSTGRES_SOCKET = (
    "/cloudsql/tcshaip:europe-west1:vv-vvta-replica"
)
POSTGRES_VERSION = "vvta_2025_02"

MYSQL_HOST = "10.73.203.7"
MYSQL_INSTANCE = "tcshaip:europe-west1:vv-vdb-replica"
MYSQL_SOCKET = (
    "/cloudsql/tcshaip:europe-west1:vv-vdb-replica"
)
MYSQL_VERSION = "vvdb_2026_07"


# ---------------------------------------------------------------------------
# Fake MySQL driver
# ---------------------------------------------------------------------------


class FakeMySQLCursor:
    """Minimal cursor implementation required by vvDBInit.Database."""

    def __init__(self):
        self.executed = []

    def execute(self, query, params=None):
        self.executed.append((query, params))

    def fetchall(self):
        return []

    def fetchone(self):
        # vvDatabase.Database.get_db_version() obtains the VVDb version
        # through the database cursor.
        return (MYSQL_VERSION,)

    def close(self):
        pass


class FakeMySQLConnection:
    """Minimal MySQL connection used by the real vvDBInit code."""

    def __init__(self):
        self.closed = False

    def ping(self, *args, **kwargs):
        return True

    def cursor(self, *args, **kwargs):
        return FakeMySQLCursor()

    def close(self):
        self.closed = True


class FakeMySQLConnectionPool:
    """
    Replacement for mysql.connector.pooling.MySQLConnectionPool.

    The real vvDBInit.Database class is still used. Only the actual database
    server connection is replaced.
    """

    calls = []

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

        FakeMySQLConnectionPool.calls.append(
            {
                "args": args,
                "kwargs": kwargs,
            }
        )

        self.connection = FakeMySQLConnection()

    def get_connection(self):
        return self.connection

    @classmethod
    def reset(cls):
        cls.calls = []


# ---------------------------------------------------------------------------
# Fake PostgreSQL driver
# ---------------------------------------------------------------------------


class FakePostgresCursor:
    """
    Minimal psycopg2 cursor implementation required by UTA.

    UTA requests DictCursor explicitly, so cursor() on the fake connection
    must accept cursor_factory.
    """

    def __init__(self):
        self.executed = []

    def execute(self, query, params=None):
        self.executed.append((query, params))

    def fetchone(self):
        query = self.executed[-1][0] if self.executed else ""

        # UTA checks whether the configured schema exists.
        if "pg_namespace" in query:
            return (True,)

        # UTA schema_version() uses a DictCursor and accesses ['value'].
        # This is the numeric matview version, not the VVTA data/schema name.
        return {"value": "0.9"}

    def fetchall(self):
        return [(POSTGRES_VERSION,)]

    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False


class FakePostgresConnection:
    """Minimal psycopg2 connection used by the real UTA dataprovider."""

    def __init__(self):
        self.closed = False
        self.autocommit = False

    def cursor(self, *args, **kwargs):
        """
        UTA calls:

            conn.cursor(
                cursor_factory=psycopg2.extras.DictCursor
            )
        """
        return FakePostgresCursor()

    def close(self):
        self.closed = True


class FakeRawPsycopg2BoundaryPool:
    """
    Exercise the real UTA ThreadedConnectionPool call down to the raw
    psycopg2.connect boundary.

    The pool itself is fake, but its constructor deliberately calls the
    patched test-only psycopg2.connect function with exactly the arguments
    supplied by UTA. This lets the test prove that no Dave-style global
    psycopg2 argument interceptor is required.
    """

    calls = []

    def __init__(self, minconn, maxconn, **kwargs):
        from vvhgvs.dataproviders import uta

        self.minconn = minconn
        self.maxconn = maxconn
        self.kwargs = kwargs

        FakeRawPsycopg2BoundaryPool.calls.append(
            {
                "minconn": minconn,
                "maxconn": maxconn,
                "kwargs": kwargs,
            }
        )

        # This is the terminal driver boundary for this test.
        # Crucially, no production interceptor sits between UTA and this
        # function. The fake connect implementation below receives exactly
        # what UTA passes.
        self.connection = uta.psycopg2.connect(**kwargs)

    def getconn(self):
        return self.connection

    def putconn(self, connection, close=False):
        if close:
            connection.close()

    def closeall(self):
        pass

    @classmethod
    def reset(cls):
        cls.calls = []


class FakePostgresConnectionPool:
    """
    Replacement for psycopg2.pool.ThreadedConnectionPool.

    The real vvhgvs UTA PostgreSQL implementation is retained. Only the
    actual psycopg2 connection creation is replaced.
    """

    calls = []

    def __init__(self, minconn, maxconn, **kwargs):
        self.minconn = minconn
        self.maxconn = maxconn
        self.kwargs = kwargs

        FakePostgresConnectionPool.calls.append(
            {
                "minconn": minconn,
                "maxconn": maxconn,
                "kwargs": kwargs,
            }
        )

        self.connection = FakePostgresConnection()

    def getconn(self):
        return self.connection

    def putconn(self, connection, close=False):
        if close:
            connection.close()

    def closeall(self):
        pass

    @classmethod
    def reset(cls):
        cls.calls = []


# ---------------------------------------------------------------------------
# Test helpers
# ---------------------------------------------------------------------------


def reset_fake_databases():
    """Reset captured fake-driver state between tests."""

    FakeMySQLConnectionPool.reset()
    FakePostgresConnectionPool.reset()
    FakeRawPsycopg2BoundaryPool.reset()


def make_config(
    *,
    postgres_host=POSTGRES_HOST,
    postgres_port="5432",
    postgres_socket=None,
    postgres_version=POSTGRES_VERSION,
    mysql_host=MYSQL_HOST,
    mysql_port="3306",
    mysql_socket=None,
    mysql_version=MYSQL_VERSION,
):
    """
    Construct a minimal configuration containing the database sections
    required by MixinInit.
    """

    config = ConfigParser()

    config["Entrez"] = {
        "email": "pytest@example.com",
        "api_key": "",
    }

    config["mysql"] = {
        "host": mysql_host,
        "port": mysql_port,
        "database": "validator",
        "user": "vvadmin",
        "password": "var1ant",
        "version": mysql_version,
    }

    if mysql_socket is not None:
        config["mysql"]["unix_socket"] = mysql_socket

    config["postgres"] = {
        "host": postgres_host,
        "port": postgres_port,
        "database": "vvta",
        "user": "uta_admin",
        "password": "uta_admin",
        "version": postgres_version,
    }

    if postgres_socket is not None:
        config["postgres"]["unix_socket"] = postgres_socket

    config["seqrepo"] = {
        "version": "VV_SR_2026_07/master",
        "location": "/data/variantvalidator/",
        "require_threading": "True",
    }

    return config


def write_config(directory: Path, config: ConfigParser):
    """Write a ConfigParser configuration to docker.ini."""

    config_path = directory / "docker.ini"

    with config_path.open("w") as handle:
        config.write(handle)

    return config_path


@pytest.fixture
def database_test_environment(tmp_path, monkeypatch):
    """
    Configure MixinInit to use a temporary docker.ini.

    MixinInit checks the existence of the configured file before loading it.
    """

    config_path = tmp_path / "docker.ini"

    monkeypatch.setattr(
        settings,
        "get_config_dir",
        lambda: str(config_path),
    )

    monkeypatch.setattr(
        "VariantValidator.modules.vvMixinInit.os.path.exists",
        lambda path: True,
    )

    return config_path


@pytest.fixture(autouse=True)
def reset_fake_database_state():
    """Reset captured driver calls around every test."""

    reset_fake_databases()

    yield

    reset_fake_databases()


@pytest.fixture
def fake_mysql_driver(monkeypatch):
    """
    Replace the actual MySQL connection pool used by vvDBInit.Database.
    """

    monkeypatch.setattr(
        "VariantValidator.modules.vvDBInit.MySQLConnectionPool",
        FakeMySQLConnectionPool,
    )

    return FakeMySQLConnectionPool


@pytest.fixture
def raw_psycopg2_driver(monkeypatch):
    """
    Test the raw psycopg2 boundary without Dave's global interceptor.

    The fake driver intentionally rejects the two malformed states that
    Dave's interceptor used to repair:

    - non-integer ports
    - bare colon-delimited Cloud SQL instance names in host

    If MixinInit normalises the configuration correctly, psycopg2 receives
    a connection request that needs no further transformation.
    """

    calls = []

    def fake_connect(*args, **kwargs):
        host = kwargs.get("host")
        port = kwargs.get("port")

        # A real psycopg2 call must already receive a usable integer port.
        if not isinstance(port, int):
            raise AssertionError(
                "MixinInit passed a non-integer PostgreSQL port to "
                "psycopg2; a global interceptor would still be required"
            )

        # A bare Cloud SQL instance name must already have been converted
        # by the Mixin before reaching psycopg2.
        if isinstance(host, str) and ":" in host:
            if not host.startswith("/cloudsql/") and not host.startswith("/"):
                raise AssertionError(
                    "MixinInit passed a bare colon-delimited Cloud SQL "
                    "host to psycopg2; a global interceptor would still "
                    "be required"
                )

        calls.append(
            {
                "args": args,
                "kwargs": dict(kwargs),
            }
        )

        return FakePostgresConnection()

    monkeypatch.setattr(
        "vvhgvs.dataproviders.uta.psycopg2.connect",
        fake_connect,
    )

    monkeypatch.setattr(
        "vvhgvs.dataproviders.uta.psycopg2.pool.ThreadedConnectionPool",
        FakeRawPsycopg2BoundaryPool,
    )

    FakeRawPsycopg2BoundaryPool.reset()

    return calls


@pytest.fixture
def fake_uta_driver(monkeypatch):
    """
    Replace only psycopg2's actual connection-pool boundary.

    The real vvhgvs UTA implementation remains active.
    """

    monkeypatch.setattr(
        "vvhgvs.dataproviders.uta.psycopg2.pool.ThreadedConnectionPool",
        FakePostgresConnectionPool,
    )

    return FakePostgresConnectionPool


@pytest.fixture
def fake_seqrepo(monkeypatch):
    """
    Prevent UTA from requiring a real local SeqRepo installation.

    This replaces the SeqFetcher dependency used by the actual UTA
    dataprovider rather than replacing UTA itself.
    """

    fake_seqfetcher = MagicMock(name="SeqFetcher")

    monkeypatch.setattr(
        "vvhgvs.dataproviders.uta.SeqFetcher",
        fake_seqfetcher,
    )

    monkeypatch.setattr(
        "vvhgvs.dataproviders.seqfetcher.SeqFetcher",
        fake_seqfetcher,
    )

    return fake_seqfetcher


def initialise_mixin(monkeypatch, config_path):
    """
    Initialise the real VariantValidator MixinInit.

    The configuration path is supplied through the actual settings mechanism.
    """

    monkeypatch.setattr(
        settings,
        "get_config_dir",
        lambda: str(config_path),
    )

    return MixinInit()


# ---------------------------------------------------------------------------
# PostgreSQL / UTA tests
# ---------------------------------------------------------------------------


def test_postgres_unix_socket_is_suppressed_for_tcp_host(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    fake_uta_driver,
    fake_seqrepo,
):
    """
    Reproduce Dave's ConfigParser behaviour.

    When postgres.host is a normal TCP host and postgres.unix_socket is also
    configured, unix_socket is suppressed and the configured host is used.
    """

    config = make_config(
        postgres_host=POSTGRES_HOST,
        postgres_port="5432",
        postgres_socket=POSTGRES_SOCKET,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    mixin = initialise_mixin(
        monkeypatch,
        config_path,
    )

    expected_url = (
        "postgresql://uta_admin:uta_admin@"
        "10.73.203.5:5432/vvta/vvta_2025_02"
    )

    assert mixin.utaPath == expected_url
    assert vvhgvs.global_config.uta.pool_max == settings.VVTA_POSTGRES_POOL_SIZE

    assert len(fake_uta_driver.calls) == 1

    call = fake_uta_driver.calls[0]

    assert call["minconn"] == vvhgvs.global_config.uta.pool_min
    assert call["maxconn"] == settings.VVTA_POSTGRES_POOL_SIZE

    assert call["kwargs"]["host"] == POSTGRES_HOST
    assert call["kwargs"]["host"] != POSTGRES_SOCKET
    assert call["kwargs"]["port"] == 5432
    assert call["kwargs"]["database"] == "vvta"
    assert call["kwargs"]["user"] == "uta_admin"
    assert call["kwargs"]["password"] == "uta_admin"

def test_postgres_host_contract_when_no_socket_is_configured(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    fake_uta_driver,
    fake_seqrepo,
):
    """
    When postgres.unix_socket is absent, the configured PostgreSQL host
    must be passed through unchanged.
    """

    # A pre-existing UTA_DB_URL must be replaced by the configured URL.
    # Dave's Mixin explicitly sets UTA_DB_URL during initialisation.
    pre_existing_uta_url = (
        "postgresql://bad:bad@bad-host:9999/bad/bad"
    )
    monkeypatch.setenv(
        "UTA_DB_URL",
        pre_existing_uta_url,
    )

    config = make_config(
        postgres_host=POSTGRES_HOST,
        postgres_port="5432",
        postgres_socket=None,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    mixin = initialise_mixin(
        monkeypatch,
        config_path,
    )

    expected_url = (
        "postgresql://uta_admin:uta_admin@"
        "10.73.203.5:5432/vvta/vvta_2025_02"
    )

    assert mixin.utaPath == expected_url
    assert os.environ["UTA_DB_URL"] == expected_url

    call = fake_uta_driver.calls[0]

    assert call["kwargs"]["host"] == POSTGRES_HOST
    assert call["kwargs"]["port"] == 5432
    assert call["kwargs"]["database"] == "vvta"
    assert call["kwargs"]["user"] == "uta_admin"
    assert call["kwargs"]["password"] == "uta_admin"


def test_postgres_cloud_sql_instance_host_contract(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    fake_uta_driver,
    fake_seqrepo,
):
    """
    A bare Cloud SQL instance string in postgres.host must be converted
    to the /cloudsql Unix socket path.
    """

    config = make_config(
        postgres_host=POSTGRES_INSTANCE,
        postgres_port="5432",
        postgres_socket=None,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    mixin = initialise_mixin(
        monkeypatch,
        config_path,
    )

    expected_url = (
        "postgresql://uta_admin:uta_admin@"
        "%2Fcloudsql%2Ftcshaip%3Aeurope-west1%3Avv-vvta-replica"
        ":5432/vvta/vvta_2025_02"
    )

    assert mixin.utaPath == expected_url

    call = fake_uta_driver.calls[0]

    assert call["kwargs"]["host"] == POSTGRES_SOCKET
    assert call["kwargs"]["port"] == 5432


def test_postgres_cloud_sql_instance_host_is_rewritten_when_socket_is_also_configured(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    fake_uta_driver,
    fake_seqrepo,
):
    """
    Reproduce Dave's combined ConfigParser and psycopg2 behaviour.

    A bare Cloud SQL instance name in postgres.host suppresses unix_socket
    because the host is not already a /cloudsql path. The host is then
    converted to /cloudsql/<instance> before psycopg2 is called.
    """

    different_socket = "/cloudsql/tcshaip:europe-west1:some-other-instance"

    config = make_config(
        postgres_host=POSTGRES_INSTANCE,
        postgres_port="5432",
        postgres_socket=different_socket,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    mixin = initialise_mixin(
        monkeypatch,
        config_path,
    )

    expected_url = (
        "postgresql://uta_admin:uta_admin@"
        "%2Fcloudsql%2Ftcshaip%3Aeurope-west1%3Avv-vvta-replica"
        ":5432/vvta/vvta_2025_02"
    )

    assert mixin.utaPath == expected_url

    call = fake_uta_driver.calls[0]
    assert call["kwargs"]["host"] == POSTGRES_SOCKET
    assert call["kwargs"]["host"] != different_socket

def test_postgres_colon_delimited_port_is_normalised(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    fake_uta_driver,
    fake_seqrepo,
):
    """
    A colon-delimited port value must arrive at psycopg2 as an integer.

    This covers the behaviour that Dave's temporary ConfigParser /
    psycopg2 patch was compensating for.
    """

    config = make_config(
        postgres_host=POSTGRES_INSTANCE,
        postgres_port="europe-west1:vv-vvta-replica:5432",
        postgres_socket=None,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    initialise_mixin(
        monkeypatch,
        config_path,
    )

    call = fake_uta_driver.calls[0]

    assert call["kwargs"]["port"] == 5432
    assert isinstance(call["kwargs"]["port"], int)


def test_postgres_invalid_colon_delimited_port_raises_value_error(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    fake_uta_driver,
    fake_seqrepo,
):
    """
    A colon-delimited PostgreSQL port whose final component is not numeric
    must be rejected by VariantValidator.

    The configured port is authoritative. VariantValidator must not invent
    or substitute a default port.
    """

    config = make_config(
        postgres_host=POSTGRES_INSTANCE,
        postgres_port="europe-west1:vv-vvta-replica:not-a-port",
        postgres_socket=None,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    with pytest.raises(ValueError):
        initialise_mixin(
            monkeypatch,
            config_path,
        )


def test_postgres_cloud_sql_socket_reaches_raw_psycopg2_without_interceptor(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    raw_psycopg2_driver,
    fake_seqrepo,
):
    """
    Prove that the new Mixin handles the complete Cloud SQL PostgreSQL
    connection configuration before psycopg2 is called.

    This deliberately does NOT install Dave's global psycopg2.connect
    interceptor. The test-only psycopg2 boundary rejects malformed host/port
    values, so success proves that the Mixin no longer depends on the global
    patch.
    """

    config = make_config(
        postgres_host=POSTGRES_INSTANCE,
        postgres_port="europe-west1:vv-vvta-replica:5432",
        postgres_socket=None,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    mixin = initialise_mixin(
        monkeypatch,
        config_path,
    )

    assert len(raw_psycopg2_driver) == 1

    call = raw_psycopg2_driver[0]

    assert call["kwargs"]["host"] == POSTGRES_SOCKET
    assert call["kwargs"]["port"] == 5432
    assert isinstance(call["kwargs"]["port"], int)
    assert call["kwargs"]["database"] == "vvta"
    assert call["kwargs"]["user"] == "uta_admin"
    assert call["kwargs"]["password"] == "uta_admin"

    assert mixin.utaPath == (
        "postgresql://uta_admin:uta_admin@"
        "%2Fcloudsql%2Ftcshaip%3Aeurope-west1%3Avv-vvta-replica"
        ":5432/vvta/vvta_2025_02"
    )


def test_postgres_explicit_socket_reaches_raw_psycopg2_without_interceptor(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    raw_psycopg2_driver,
    fake_seqrepo,
):
    """
    An explicit /cloudsql Unix socket is suppressed when postgres.host is a
    normal TCP host, matching Dave's ConfigParser behaviour.
    """

    config = make_config(
        postgres_host=POSTGRES_HOST,
        postgres_port="5432",
        postgres_socket=POSTGRES_SOCKET,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    initialise_mixin(
        monkeypatch,
        config_path,
    )

    assert len(raw_psycopg2_driver) == 1

    call = raw_psycopg2_driver[0]

    assert call["kwargs"]["host"] == POSTGRES_HOST
    assert call["kwargs"]["host"] != POSTGRES_SOCKET
    assert call["kwargs"]["port"] == 5432
    assert isinstance(call["kwargs"]["port"], int)

def test_postgres_tcp_host_reaches_raw_psycopg2_without_interceptor(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    raw_psycopg2_driver,
    fake_seqrepo,
):
    """
    A conventional TCP PostgreSQL host must also reach the raw driver in a
    directly usable form, without relying on a global psycopg2 interceptor.
    """

    config = make_config(
        postgres_host=POSTGRES_HOST,
        postgres_port="5432",
        postgres_socket=None,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    initialise_mixin(
        monkeypatch,
        config_path,
    )

    assert len(raw_psycopg2_driver) == 1

    call = raw_psycopg2_driver[0]

    assert call["kwargs"]["host"] == POSTGRES_HOST
    assert call["kwargs"]["port"] == 5432
    assert isinstance(call["kwargs"]["port"], int)


def test_postgres_pool_size_is_honoured(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    fake_uta_driver,
    fake_seqrepo,
):
    """
    Verify that VVTA_POSTGRES_POOL_SIZE reaches the actual UTA
    ThreadedConnectionPool.
    """

    monkeypatch.setattr(
        settings,
        "VVTA_POSTGRES_POOL_SIZE",
        7,
    )

    config = make_config(
        postgres_host=POSTGRES_HOST,
        postgres_port="5432",
        postgres_socket=None,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    initialise_mixin(
        monkeypatch,
        config_path,
    )

    call = fake_uta_driver.calls[0]

    assert call["maxconn"] == 7


# ---------------------------------------------------------------------------
# MySQL / vvDBInit tests
# ---------------------------------------------------------------------------


def test_mysql_explicit_unix_socket_contract(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    fake_uta_driver,
    fake_seqrepo,
):
    """
    Exercise the real vvDBInit.Database path and verify the exact
    configuration supplied to MySQLConnectionPool.
    """

    config = make_config(
        mysql_host=MYSQL_HOST,
        mysql_port="3306",
        mysql_socket=MYSQL_SOCKET,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    mixin = initialise_mixin(
        monkeypatch,
        config_path,
    )

    assert mixin.dbConfig["host"] == MYSQL_HOST
    assert mixin.dbConfig["port"] == 3306
    assert mixin.dbConfig["unix_socket"] == MYSQL_SOCKET
    assert mixin.dbConfig["use_pure"] is True
    assert mixin.dbConfig["connection_timeout"] == 15

    assert len(fake_mysql_driver.calls) == 1

    call = fake_mysql_driver.calls[0]

    assert call["kwargs"]["host"] == MYSQL_HOST
    assert call["kwargs"]["port"] == 3306
    assert call["kwargs"]["unix_socket"] == MYSQL_SOCKET
    assert call["kwargs"]["use_pure"] is True
    assert call["kwargs"]["connection_timeout"] == 15


def test_mysql_host_contract_without_socket(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    fake_uta_driver,
    fake_seqrepo,
):
    """
    When mysql.unix_socket is absent, the configured host must be
    passed directly to MySQL Connector/Python.
    """

    config = make_config(
        mysql_host=MYSQL_HOST,
        mysql_port="3306",
        mysql_socket=None,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    mixin = initialise_mixin(
        monkeypatch,
        config_path,
    )

    assert mixin.dbConfig["host"] == MYSQL_HOST
    assert mixin.dbConfig["port"] == 3306
    assert "unix_socket" not in mixin.dbConfig

    call = fake_mysql_driver.calls[0]

    assert call["kwargs"]["host"] == MYSQL_HOST
    assert call["kwargs"]["port"] == 3306
    assert "unix_socket" not in call["kwargs"]


def test_mysql_cloud_sql_instance_host_is_not_rewritten(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    fake_uta_driver,
    fake_seqrepo,
):
    """
    MySQL uses its explicit unix_socket parameter.

    The host itself must therefore remain the configured host rather
    than being rewritten to /cloudsql/... by MixinInit.
    """

    config = make_config(
        mysql_host=MYSQL_INSTANCE,
        mysql_port="3306",
        mysql_socket=MYSQL_SOCKET,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    mixin = initialise_mixin(
        monkeypatch,
        config_path,
    )

    assert mixin.dbConfig["host"] == MYSQL_INSTANCE
    assert mixin.dbConfig["unix_socket"] == MYSQL_SOCKET

    call = fake_mysql_driver.calls[0]

    assert call["kwargs"]["host"] == MYSQL_INSTANCE
    assert call["kwargs"]["unix_socket"] == MYSQL_SOCKET


def test_mysql_pool_size_is_honoured(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    fake_uta_driver,
    fake_seqrepo,
):
    """
    Verify that VALIDATOR_MYSQL_POOL_SIZE reaches the actual
    MySQLConnectionPool constructor.
    """

    monkeypatch.setattr(
        settings,
        "VALIDATOR_MYSQL_POOL_SIZE",
        6,
    )

    config = make_config(
        mysql_host=MYSQL_HOST,
        mysql_port="3306",
        mysql_socket=MYSQL_SOCKET,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    initialise_mixin(
        monkeypatch,
        config_path,
    )

    call = fake_mysql_driver.calls[0]

    assert call["kwargs"]["pool_size"] == 6


def test_mysql_port_is_integer(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    fake_uta_driver,
    fake_seqrepo,
):
    """
    Verify that the final MySQL connector configuration contains a
    real integer port.
    """

    config = make_config(
        mysql_host=MYSQL_HOST,
        mysql_port="3306",
        mysql_socket=None,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    mixin = initialise_mixin(
        monkeypatch,
        config_path,
    )

    assert mixin.dbConfig["port"] == 3306
    assert isinstance(mixin.dbConfig["port"], int)

    call = fake_mysql_driver.calls[0]

    assert call["kwargs"]["port"] == 3306
    assert isinstance(call["kwargs"]["port"], int)


# ---------------------------------------------------------------------------
# Combined Cloud Run configuration contract
# ---------------------------------------------------------------------------


def test_mysql_colon_delimited_port_is_normalised(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    fake_uta_driver,
    fake_seqrepo,
):
    """
    A colon-delimited MySQL port must be normalised by VariantValidator
    before the real vvDBInit.Database layer receives the configuration.
    """

    config = make_config(
        mysql_host=MYSQL_INSTANCE,
        mysql_port="europe-west1:vv-vdb-replica:3306",
        mysql_socket=MYSQL_SOCKET,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    initialise_mixin(
        monkeypatch,
        config_path,
    )

    call = fake_mysql_driver.calls[0]

    assert call["kwargs"]["port"] == 3306
    assert isinstance(call["kwargs"]["port"], int)


def test_cloud_run_database_configuration_contract(
    monkeypatch,
    database_test_environment,
    fake_mysql_driver,
    fake_uta_driver,
    fake_seqrepo,
):
    """
    Exercise both production database initialisation paths using the
    configuration shape used by Cloud Run.

    PostgreSQL:
        normal TCP host plus unix_socket configured, following Dave's
        suppression behaviour

    MySQL:
        explicit /cloudsql socket

    No production connection-building code is mocked.
    """

    config = make_config(
        postgres_host=POSTGRES_HOST,
        postgres_port="5432",
        postgres_socket=POSTGRES_SOCKET,
        mysql_host=MYSQL_HOST,
        mysql_port="3306",
        mysql_socket=MYSQL_SOCKET,
    )

    config_path = write_config(
        database_test_environment.parent,
        config,
    )

    mixin = initialise_mixin(
        monkeypatch,
        config_path,
    )

    # MySQL contract.
    assert mixin.dbConfig["host"] == MYSQL_HOST
    assert mixin.dbConfig["port"] == 3306
    assert mixin.dbConfig["unix_socket"] == MYSQL_SOCKET
    assert mixin.dbConfig["use_pure"] is True
    assert mixin.dbConfig["connection_timeout"] == 15

    mysql_call = fake_mysql_driver.calls[0]

    assert mysql_call["kwargs"]["host"] == MYSQL_HOST
    assert mysql_call["kwargs"]["port"] == 3306
    assert mysql_call["kwargs"]["unix_socket"] == MYSQL_SOCKET

    # PostgreSQL contract.
    pg_call = fake_uta_driver.calls[0]

    assert pg_call["kwargs"]["host"] == POSTGRES_HOST
    assert pg_call["kwargs"]["host"] != POSTGRES_SOCKET
    assert pg_call["kwargs"]["port"] == 5432
    assert pg_call["kwargs"]["database"] == "vvta"
    assert pg_call["kwargs"]["user"] == "uta_admin"
    assert pg_call["kwargs"]["password"] == "uta_admin"

    assert mixin.utaPath == (
        "postgresql://uta_admin:uta_admin@"
        "10.73.203.5:5432/vvta/vvta_2025_02"
    )

# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
