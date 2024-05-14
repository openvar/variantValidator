# Installation

These instructions will allow you to install the package and accompanying databases on Windows.

## Pre-requisites

Required:
* MySQL
* Python 3.6 or above
* SQLite version 3.8.0 or above

Optional:
* PostgreSQL version 10.5 or above. 

## Install Linux as Windows subsystem
Windows Subsystem for Linux (WSL) can be installed by following the instructions [here](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-11-with-gui-support#1-overview). This installs Ubuntu as a subsystem. 

The rest of the installation is done via the Ubuntu application. 

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
$ conda activate vvenvcond
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

You need to check the mySQL service is working correctly and connect to the server. It can be done so in the following way:

1. Start the mySQL service

```
$ sudo service mysql start
```

2. Check the status of the service

```
$ sudo service mysql status
```

The output should look something like this:

```
* /usr/bin/mysqladmin  Ver 8.0.29-0ubuntu0.20.04.3 for Linux on x86_64 ((Ubuntu))
Copyright (c) 2000, 2022, Oracle and/or its affiliates.

Oracle is a registered trademark of Oracle Corporation and/or its
affiliates. Other names may be trademarks of their respective
owners.

Server version          8.0.29-0ubuntu0.20.04.3
Protocol version        10
Connection              Localhost via UNIX socket
UNIX socket             /var/run/mysqld/mysqld.sock
Uptime:                 28 min 28 sec

Threads: 2  Questions: 10  Slow queries: 0  Opens: 117  Flush tables: 3  Open tables: 36  Queries per second avg: 0.005
```

3. Start the security script

```
sudo mysql_secure_installation
```

This command prompts you to answer a series of security questions to configure and secure the mySQL server. 

4. Connect to the mySQL server

```
sudo mysql -u root -p
```

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

Download and our pre-populated database to MySQL as follows. 

***Essential Step: check [here](https://www528.lamp.le.ac.uk/vvdata/validator/) and make sure you download and install the most up-to-date version***

```bash
$ wget https://www528.lamp.le.ac.uk/vvdata/validator/validator_2022_04.sql.gz
$ gunzip validator_2022_04.sql.gz
$ mysql validator < validator_2022_04.sql -u USER -p
```

See the [Manual](MANUAL.md) for instructions on updating this database, which should be done regularly.

If you wish to test your installation using pytest (see below) we recommend that you do this before updating the database. 

## Setting up Seqrepo (SQLite >=3.8)

VariantValidator requires a local SeqRepo database. The seqrepo package has already been installed into the virtual environment, but you'll need to download an actual seqrepo database. This can go anywhere on your system drive.

***Essential Step: check [here](https://www528.lamp.le.ac.uk/vvdata/vv_seqrepo/) and make sure you download and install the most up-to-date version where the required file is
 e.g. VV_SR_2021_2.tar and the numbers indicate the creation date i.e. 2021_02 = February 2021***

```
$ mkdir /path/to/seqrepo
$ cd /path/to/seqrepo
$ wget https://www528.lamp.le.ac.uk/vvdata/vv_seqrepo/VV_SR_2022_02.tar
$ tar -xvf VV_SR_2022_02.tar
$ rm VV_SR_2022_02.tar
```
where /path/to/seqrepo should be where you install the database e.g. /Users/Shared/seqrepo_dumps/ or /local/seqrepo


## Setting up VVTA database (PostGreSQL >=10.5)

You need to start the PostGreSQL service before creating the database. It can be done so in the following way:

```
sudo service postgresql start
```

```
su - postgres
psql
```

You will need to install a local version of the VVTA database. 

First create the database and a user account:

```psql
CREATE USER <USER> WITH CREATEDB PASSWORD '<password>';
CREATE DATABASE vvta WITH OWNER=<USER> TEMPLATE=template0;
```

Where:
- \<USER\> should be a user-name e.g. uta_admin
- password is a unique password for user

To fill this database, download the gzipped uta genetics database, and upload it into psql.

***Essential Step: check [here](https://www528.lamp.le.ac.uk/vvdata/vvta/) and make sure you download and install the most up-to-date version***

```
$ wget --output-document=VVTA_2022_02.noseq.sql.gz https://www528.lamp.le.ac.uk/vvdata/vvta/VVTA_2022_02_noseq.sql.gz
$ gzip -cdq VVTA_2022_02.noseq.psql.gz | psql -U <USER> -v ON_ERROR_STOP=0 -d vvta -Eae
```

***Possible error***
If you get a 'Peer authnetication' error whilst doing the gzip command, you may need to update the [pg_hba.conf file](https://www.postgresql.org/docs/8.3/auth-pg-hba-conf.html). Change the method of authentication from peer to md5 and restart the postgresql service. 

More information on fixing it can be found [here](https://stackoverflow.com/a/21889759).

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
