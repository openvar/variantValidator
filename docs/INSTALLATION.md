# Installation

These instructions will allow you to install the package and accompanying databases on Linux. Mac OS X computers operate similarly.
For any other systems, or if you cannot install the databases, we recommend installing via [docker](DOCKER.md).

## Pre-requisites

Required:
* MySQL
* Python 3.6 or above
* SQLite version 3.8.0 or above

Optional:
* PostgreSQL version 10.5 or above. 


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

## Additional steps for running MariaDB
Install the mariadb python library

***Note: Only do this if you intend to run MariaDB instead of MySQL***
```bash
$ pip install mariadb
```
***Additional steps may be required***
Install [MariaDB Connector/C](https://downloads.mariadb.com/Connectors/c/)
Installation instructions can be found [here](https://mariadb.com/kb/en/about-mariadb-connector-c/) 

## Installing VariantValidator

Hint: your new environment vvenv should still be activated from the previous steps and you should still be in the /variantValidator directory where setup.py is located.

To install VariantValidator within your virtual environment run:
```
$ pip install -e .
```

## Setting up validator database (MySQL)

A MySQL database called validator is required to run VariantValidator. We recommend creating a user and password specific to the
validator database, for example:

```mysql
CREATE USER 'USER'@'HOST' IDENTIFIED BY 'PASSWORD';
CREATE DATABASE validator;
GRANT SELECT,INSERT,UPDATE,DELETE ON validator.* TO 'USER'@'HOST';
```
Where:
- USER should be a user-name e.g. vvadmin
- HOST is the MySQL host ID, usually 127.0.0.1
- PASSWORD is a unique password for your database

*Note: We have had reports that on some systems ALL PRIVILEGES may be required rather than SELECT,INSERT,UPDATE,DELETE*

Download and our pre-populated database to MySQL as follows. ***Note: check [here](https://www528.lamp.le.ac.uk/vvdata/validator/) for the most up-to-date version***

```bash
$ wget https://www528.lamp.le.ac.uk/vvdata/validator/validator_2021-07-21.sql.gz
$ gunzip validator_2021-07-21.sql.gz
$ mysql validator < validator_2021-07-21.sql -u HOST -p
```

See the [Manual](MANUAL.md) for instructions on updating this database, which should be done regularly.

If you wish to test your installation using pytest (see below) we recommend that you do this before updating the database. 

## Setting up Seqrepo (SQLite >=3.8)

VariantValidator requires a local SeqRepo database. The seqrepo package has already been installed into the virtual environment, but you'll need to download an actual seqrepo database. This can go anywhere on your system drive.
***Note: check [here](https://www528.lamp.le.ac.uk/vvdata/vv_seqrepo/) for the most up-to-date version where the required file is
 e.g. VV_SR_2021_2.tar and the numbers indicate the creation date i.e. 2021_02 = February 2021***

```
$ mkdir /path/to/seqrepo
$ cd /path/to/seqrepo
$ wget https://www528.lamp.le.ac.uk/vvdata/vv_seqrepo/VV_SR_2021_2.tar
$ tar -xvf VV_SR_2021_2.tar
$ rm VV_SR_2021_2.tar
```
where /path/to/seqrepo should be where you install the database e.g. /Users/Shared/seqrepo_dumps/ or /local/seqrepo


## Setting up VVTA database (PostGreSQL >=10.5)

You will need to install a local version of the VVTA database. 

First create the database and a user account:

```
psql
CREATE ROLE <USER> WITH CREATEDB;
ALTER ROLE <USER> WITH LOGIN;
ALTER ROLE <USER> WITH PASSWORD '<password>';
CREATE DATABASE vvta WITH OWNER=<USER> TEMPLATE=template0;
```
Where:
- \<USER\> should be a user-name e.g. uta_admin
- password is a unique password for user

To fill this database, download the gzipped uta genetics database, and upload it into psql.
```
$ wget --output-document=VVTA_2021_2_noseq.psql.gz https://www528.lamp.le.ac.uk/vvdata/vvta/VVTA_2021_2_noseq.psql.gz
$ gzip -cdq VVTA_2021_2_noseq.psql.gz | psql -U <USER> -v ON_ERROR_STOP=0 -d vvta -Eae
```

## Configuration

Before using VariantValidator some configuration is required, as described in the [Manual](MANUAL.md).

## Developers

To work on the VariantValidator code, you'll need to install additional dependencies and install VariantValidator in an editable manner. Tests can be run using PyTest.

```bash
cd variantValidator/
pip install -r requirements.txt
pip install -r requirements_dev.txt
pip install -e .
pytest
```
  
Please make all Pull Requests to the develop branch. Id you are unsure, contact admin via [issues](https://github.com/openvar/variantValidator/issues)
