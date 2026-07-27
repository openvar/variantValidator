# Installation

VariantValidator is developed and tested primarily on Linux. macOS systems operate similarly.

For most users and developers, we recommend using the **Quick Start** installation, which uses Docker to provide the required VariantValidator databases while installing VariantValidator itself in a local Python environment.

Users who prefer to install and manage all database components themselves can follow the **Native installation** instructions below.

Windows users should see the [Windows Installation Guide](windows_installation.md).

---

# Quick Start (Recommended)

The recommended installation method is to use Docker to provide the Validator, SeqRepo and VVTA databases while installing VariantValidator itself in a local Python environment.

This approach avoids having to install, populate and configure the databases individually while retaining full access to the VariantValidator Python API, command-line tools and development environment.

Follow the **Quick Start** section of the [Docker Installation Guide](docker.md).

If you prefer to run the entire VariantValidator software stack in Docker, see the **Full Docker Installation** section of the [Docker Installation Guide](docker.md).

---

# Native installation

The following instructions describe how to install VariantValidator and its accompanying databases natively on Linux or macOS.

## Pre-requisites

Required:

* MySQL
* Python 3.6 or above
* SQLite version 3.8.0 or above

Optional:

* PostgreSQL version 10.5 or above

---

## Download the source code

To download the VariantValidator source code, clone the repository:

```bash
git clone https://github.com/openvar/variantValidator.git
cd variantValidator/
```

---

## Python environment

When installing VariantValidator we recommend using a virtual environment, as it requires specific versions of several libraries including Python and SQLite.

This can be done using either conda or pip.

### Via conda (Recommended)

After installing conda, create a new virtual environment using:

```bash
conda env create -f environment.yml
conda activate vvenv
```

The packages required for VariantValidator to function are now installed in the `vvenv` environment.

### Via pip

If you already have suitable versions of Python and SQLite installed, you can create a Python virtual environment using:

```bash
python -m venv vvenv
source vvenv/bin/activate
```

The VariantValidator dependencies will be installed when VariantValidator itself is installed below.

---

## Additional steps for MariaDB

If you intend to use MariaDB instead of MySQL, install the MariaDB Python library:

```bash
pip install mariadb
```

You may also need to install MariaDB Connector/C.

See the MariaDB documentation for installation instructions.

!!! note

    These additional steps are only required if you intend to use MariaDB instead of MySQL.

---

## Installing VariantValidator

Ensure that the `vvenv` environment is activated and that you are in the `variantValidator` repository directory.

Install VariantValidator using:

```bash
pip install .
python -m VariantValidator.bin.setup_lovd_syntax_checker
```

---

# Setting up the Validator database

VariantValidator requires a MySQL database named `validator`.

We recommend creating a dedicated database user and password:

```mysql
CREATE USER 'USER'@'HOST' IDENTIFIED WITH mysql_native_password BY 'PASSWORD';
CREATE DATABASE validator;
GRANT SELECT,INSERT,UPDATE,DELETE ON validator.* TO 'USER'@'HOST';
```

Where:

* `USER` is the database username, for example `vvadmin`.
* `HOST` is the MySQL host, usually `127.0.0.1`.
* `PASSWORD` is a unique password for the database user.

!!! note

    We have received reports that some systems require `ALL PRIVILEGES` rather than `SELECT,INSERT,UPDATE,DELETE`.

Download and install the latest pre-populated Validator database.

!!! important

    Check the VariantValidator data repository and ensure that you download the most recent available database.

```bash
wget https://data.variantvalidator.org/vvdata/validator/validator_202x-xx-xx.sql.gz
gunzip validator_202x-xx-xx.sql.gz
mysql validator < validator_202x-xx-xx.sql -u USER -p
```

See the database update documentation for instructions on keeping this database up to date.

If you intend to test the installation using pytest, we recommend doing so before updating the database.

---

# Setting up SeqRepo

VariantValidator requires a local SeqRepo database.

The SeqRepo Python package is installed as a VariantValidator dependency, but the sequence database itself must also be downloaded.

!!! important

    Check the VariantValidator data repository and download the most recent SeqRepo release.

For example:

```bash
mkdir /path/to/seqrepo
cd /path/to/seqrepo

wget https://data.variantvalidator.org/vvdata/vv_seqrepo/VV_SR_20xx_xx.tar
tar -xvf VV_SR_20xx_xx.tar
rm VV_SR_20xx_xx.tar
```

`/path/to/seqrepo` should be replaced with the location where you wish to store the sequence database, for example:

```text
/Users/Shared/seqrepo_dumps/
```

or:

```text
/local/seqrepo
```

---

# Setting up VVTA

VariantValidator requires a local VVTA PostgreSQL database.

First create the database and user account:

```sql
psql

CREATE ROLE <USER> WITH CREATEDB;
ALTER ROLE <USER> WITH LOGIN;
ALTER ROLE <USER> WITH PASSWORD '<PASSWORD>';
CREATE DATABASE vvta WITH OWNER=<USER> TEMPLATE=template0;
```

Where:

* `<USER>` is the PostgreSQL username, for example `uta_admin`.
* `<PASSWORD>` is a unique password for that user.

Download and import the latest VVTA database.

!!! important

    Check the VariantValidator data repository and ensure that you download the most recent available VVTA database.

For PostgreSQL versions earlier than 14:

```bash
wget --output-document=VVTA_202x_xx.noseq.psql.gz https://data.variantvalidator.org/vvdata/vvta/VVTA_202x_xx.noseq.psql.gz
gzip -cdq VVTA_202x_xx.noseq.psql.gz | psql -U <USER> -v ON_ERROR_STOP=1 -d vvta -Eae
```

For PostgreSQL 14 and above:

```bash
wget --output-document=VVTA_202x_xx.noseq.psql.gz https://data.variantvalidator.org/vvdata/vvta/VVTA_202x_xx.noseq.psql.gz
gzip -cdq -k VVTA_202x_xx.noseq.psql.gz | sed 's/anyarray/anycompatiblearray/g' | psql -U <USER> -v ON_ERROR_STOP=1 -d vvta -Eae
```

---

# Configure MySQL for testing

The complete VariantValidator test suite creates a large number of simultaneous database connections. Before running the tests in parallel against a native MySQL installation, we recommend increasing the MySQL `max_connections` setting.

Edit your MySQL configuration file (typically `my.cnf` or `mysqld.cnf`) and set:

```ini
[mysqld]
max_connections = 1000
```

Restart the MySQL server after making this change.

If you do not increase the connection limit, the test suite may fail with MySQL connection errors when tests are run in parallel.

---

# Configuration

Before VariantValidator can be used, the installation must be configured.

Configuration is described in the [Configuration Guide](configuration_cli.md).

---

# Verify the installation

After completing the installation and configuration, we recommend running the complete VariantValidator test suite to verify that the software and databases have been installed correctly.

If you are using the Docker databases from the Quick Start installation, the test suite can be run in parallel:

```bash
pytest -n 4
```

If you performed a fully native installation and have not increased the MySQL connection limit, run the test suite using a single process:

```bash
pytest
```

If `max_connections` has been increased appropriately, the native installation can also run the tests in parallel:

```bash
pytest -n 4
```

The test suite performs functional tests covering VariantValidator, VariantFormatter, the Validator database, SeqRepo, VVTA and external dependencies.

A successful installation should complete the test suite without failures.

---

# Troubleshooting

If you encounter installation or configuration problems, see the [Configuration Troubleshooting Guide](configuration_troubleshooting.md).
