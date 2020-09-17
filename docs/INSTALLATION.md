# Installation

These instructions will allow you to install the package and accompanying databases on Linux. Mac OS X computers operate similarly.
For any other systems, or if you cannot install the databases, we recommend installing via [docker](DOCKER.md).

## Pre-requisites

Required:
* MySQL
* Python 3.6 or above
* SQLite version 3.8.0 or above

Optional:
* PostgreSQL version 9.5 or above. 

## Download the source code

To download the VariantValidator source code simply clone the master branch.

```
$ git clone https://github.com/openvar/variantValidator.git
$ cd variantValidator/
```

## Python 3.6 environment

When installing VariantValidator we recommend using a virtual environment, as it requires specific versions of several libraries including python and sqlite. This can be done either via conda **or** pip.

#### Via conda  
After [installing conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) you can create a new virtual environment with the correct python and sqlite versions by running:
```
$ conda env create -f environment.yml
$ conda activate vvenv
```
The packages required for VariantValidator to function are now set up in the environment "vvenv".

#### Via pip

If you already have the right versions of python (>=3.6) and sqlite (>=3.8), then you can use pip to BT create a new virtual envionment and install the packages required for VariantValidator to function.

```
$ python -m venv vvenv
$ source activate vvenv
$ pip install -r requirements.txt
```

## Installing VariantValidator

Hint: your new environment vvenv should still be activated from the previous and you should still be in the /variantValidator directory where setup.py is located.

To install VariantValidator within your virtual environment run:
```
$ python setup.py install
```

## Setting up validator database (MySQL)

A MySQL database is required to run VariantValidator. We recommend creating a user and password specific to the
VariantValidator database, for example:

```
$ mysql CREATE USER '<USER>'@'<HOST>' IDENTIFIED BY '<PASSWORD>';
CREATE DATABASE validator;
GRANT SELECT,INSERT,UPDATE,DELETE ON validator.* TO '<USER>'@'<HOST>';
```
Where:
- <USER> should be a user-name e.g. vvadmin
- <HOST> is the MySQL host ID, usually 127.0.0.1
- <PASSWORD> is a unique password for your database

*Note: We have had reports that on some systems ALL PRIVILEGES may be required rather than SELECT,INSERT,UPDATE,DELETE*

In the `VariantValidator/configuration` folder is a copy of the empty mysql database needed by VariantValidator to run. You need to upload it to the running MySQL database with:
```
$ mysql validator < VariantValidator/configuration/empty_vv_db.sql 
```
However, we highly recommend that you download and and upload our pre-populated database to MySQL. The current version can be accessed as follows

```bash
$ wget --output-document=validator_2019-09-18.sql.gz https://leicester.figshare.com/ndownloader/files/17714429
$ gunzip validator_2019-09-18.sql.gz
$ mysql validator < validator_2019-09-18.sql
```

See the [Manual](MANUAL.md) for instructions on updating this database, which should be done regularly.

If you wish to test your installation using pytest (see below) we recommend that you do this before updating the database. 

## Setting up Seqrepo (SQLite >=3.8)

VariantValidator requires a local SeqRepo database. The seqrepo package has already been installed into the virtual environment, but you'll need to download an actual seqrepo database. This can go anywhere on your system drive.

```
$ mkdir /path/to/seqrepo
$ seqrepo --root-directory /path/to/seqrepo pull -i 2018-08-21
```
where /path/to/seqrepo should be where you install the database e.g. /Users/Shared/seqrepo_dumps/

To check it has downloaded:
```
$ seqrepo --root-directory /path/to/seqrepo list-local-instances
```
where the output should be a list of available seqrepo databases e.g. 2018-08-21

## Setting up UTA database (Optional, PostGreSQL >=9.5)

It's recommended for performance reasons to use a local version of the UTA databases. Alternatively, see below if you wish to use the remote version.
We again recommend creating a specific user account, for example:

```
psql
CREATE ROLE <USER> WITH CREATEDB;
ALTER ROLE <USER> WITH LOGIN;
\password
CREATE DATABASE uta WITH OWNER=<USER> TEMPLATE=template0;
```
Where:
- <USER> should be a user-name e.g. uta_admin
- password is a unique password for your local version of the database

To fill this database, download the gzipped uta genetics database, and upload it into psql.
```
$ wget https://leicester.figshare.com/ndownloader/files/17797259 -O /docker-entrypoint-initdb.d/uta_20180821.sql.gz
$ gzip -cdq uta_20180821.pgd.gz | psql -U uta_admin -v ON_ERROR_STOP=0 -d uta -Eae
```

If you wish to use the remote, public UTA database, see the instructions [here](https://github.com/biocommons/uta#accessing-the-public-uta-instance).

## Configuration

Before using VariantValidator some configuration is required, as described in the [Manual](MANUAL.md).

## Developers

To work on the VariantValidator code, you'll need to install additional dependencies and install VariantValidator in an editable manner. Tests can be run using PyTest.

```bash
cd variantValidator/
pip install requirements_dev.txt
pip install -e .
pytest
```
Note: -e and the full-stop '.' tells pip to install in editable or development mode in the current directory.
  
Please make all Pull Requests to the develop_v3 branch.
