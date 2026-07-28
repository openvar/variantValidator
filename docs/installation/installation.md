<img src="../static/img/logos/VV_logo.png" width="20%" />

# Installation

VariantValidator is developed and tested primarily on Linux. macOS systems operate similarly.

For most users and developers, we recommend using the **Quick Start** installation, which uses Docker to provide the required VariantValidator databases while installing VariantValidator and VariantFormatter in the VariantValidator conda environment.

Users who prefer to install and manage all database components themselves can follow the **Native installation** instructions below.

Windows users should see the [Windows Installation Guide](windows_installation.md).

## See also

- [Docker Installation Guide](docker.md) — Recommended installation using Docker-hosted databases.
- [Windows Installation Guide](windows_installation.md) — Installation guidance for Windows users.
- [Configuration Guide](configuration_cli.md) — Configure VariantValidator after installation.
- [Configuration Troubleshooting Guide](configuration_troubleshooting.md) — Resolve common installation and configuration problems.

---

# Quick Start (Recommended)

The recommended installation method is to use Docker to provide the Validator, SeqRepo and VVTA databases while installing VariantValidator and VariantFormatter in the VariantValidator conda environment.

This approach avoids having to install, populate and configure the databases individually while retaining full access to the VariantValidator and VariantFormatter Python APIs, command-line tools and development environment.

Follow the **Quick Start** section of the [Docker Installation Guide](docker.md).

If you prefer to run the entire VariantValidator software stack in Docker, see the **Full Docker Installation** section of the [Docker Installation Guide](docker.md).

---

# Native installation

The following instructions describe how to install VariantValidator, VariantFormatter and their accompanying databases natively on Linux or macOS.

## Pre-requisites

Required:

* MySQL
* PostgreSQL
* Conda

---

## Download the source code

To download the VariantValidator source code, clone the repository:

```bash
git clone https://github.com/openvar/variantValidator.git
cd variantValidator/
```

The repository contains both VariantValidator and VariantFormatter. They are maintained together but are packaged as independently installable Python distributions.

---

## Create the VariantValidator conda environment

VariantValidator should be installed in the conda environment defined by the repository `environment.yml` file.

Create and activate the environment using:

```bash
conda env create -f environment.yml
conda activate vvenv
```

The required Python, SQLite and other conda-managed dependencies are installed into the `vvenv` environment.

---

## Additional steps for MariaDB

If you intend to use MariaDB instead of MySQL, install the MariaDB Python library into the activated `vvenv` environment:

```bash
pip install mariadb
```

You may also need to install MariaDB Connector/C.

See the MariaDB documentation for installation instructions.

!!! note

    These additional steps are only required if you intend to use MariaDB instead of MySQL.

---

## Installing VariantValidator and VariantFormatter

Ensure that the `vvenv` conda environment is activated and that you are in the `variantValidator` repository directory.

VariantValidator and VariantFormatter are packaged as separate Python distributions within the same repository.

For development, install both distributions in editable mode:

```bash
pip install -e ./packaging/variantvalidator
pip install -e ./packaging/variantformatter
python -m VariantValidator.bin.setup_lovd_syntax_checker
```

VariantFormatter declares VariantValidator as a dependency, while VariantValidator can be installed independently of VariantFormatter.

To install only VariantValidator:

```bash
pip install -e ./packaging/variantvalidator
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

---

# Getting help

VariantValidator has been developed to support a wide range of users, from those new to HGVS nomenclature to experienced clinical scientists and bioinformaticians. If you encounter difficulties installing, configuring or testing VariantValidator, we encourage you to seek assistance.

Before contacting the development team, you may find the following documentation helpful:

- [Docker Installation Guide](docker.md)
- [Configuration Guide](configuration_cli.md)
- [Configuration Troubleshooting Guide](configuration_troubleshooting.md)

If you still require assistance, you can contact the VariantValidator team using our [contact form](https://variantvalidator.org/help/contact/).

Software bugs and feature requests can be reported through the [VariantValidator GitHub issue tracker](https://github.com/openvar/VariantValidator/issues).

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="../static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="../static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>

