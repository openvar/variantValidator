# Variant Validator Installation

These instructions will allow you to configure the software on Linux and Mac OS X computers.

There are several steps involved in setting up variant validator:
* The python environment must be set up with the correct packages
* The variantValidator files themselves must be downloaded and installed.
* The databases must be downloaded and set up
* The configuration files must be changed to point the validator at those databases.

## Virtual environment

Variant validator currently requires python 2.7.

When installing Variant Validator it is wise to use a virtual environment, as it requires specific versions of several libraries.
We recommend using conda.
```
$ conda create -n VVenv
$ conda activate VVenv
$ conda install -c conda-forge sqlite python=2.7 pytest protobuf=3.5.1 docutils python-daemon httplib2 mysql-connector-python mysql-python 
$ conda install -c auto biotools
$ conda install -c bioconda pyliftover pysam
$ conda install -c conda-forge setuptools numpy
$ pip install hgvs==1.1.3
```
The packages required for variant validator to function are now set up in the environment "VVenv".

## Installing validator code

To clone this software from GIT, use:
```
$ git clone https://github.com/UniOfLeicester/variantValidator.git
$ cd variantValidator/
```
Run the installation script to integrate variant validator with python's site packages.
```
$ python setup.py install
```
For development purposes, you can use
```
$ pip install -e .
```
to ensure any changes you make in the local variant validator folder is reflected in your python site-packages.

## Setting up MySQL

A MySQL database is required to run variantValidator. We recommend creating a user and password specific to the
variant Validator database.

```mysql
CREATE USER 'vvadmin'@'localhost' IDENTIFIED BY 'var1ant';
CREATE DATABASE validator;
GRANT SELECT,INSERT,UPDATE,DELETE ON validator.* TO vvadmin;
```

In the `VariantValidator/data` folder is a copy of the empty mysql database needed by Variant Validator to run. The software will populate it as variants are run. You need to upload it to the running MySQL database with:
```
$ mysql validator < VariantValidator/data/emptyValidatorDump.sql 
```


## Setting up PostGreSQL

It's recommended for performance reasons to use a local varsion of the UTA database. We again recommend creating a specific user account.
```
CREATE ROLE uta_admin WITH CREATEDB;
ALTER ROLE uta_admin WITH LOGIN;
\password
CREATE DATABASE uta WITH OWNER=uta_admin TEMPLATE=template0;
```

To fill this database, download the gzipped uta genetics database, and upload it into psql.
```
$ wget http://dl.biocommons.org/uta/uta_20180821.pgd.gz
$ gzip -cdq uta_20180821.pgd.gz | psql -U uta_admin -v ON_ERROR_STOP=0 -d uta -Eae
```


## Setting up Seqrepo

Similarly, things run much faster with a local SeqRepo database. The seqrepo library is already installed, but you'll need to download an actual sequence repository. These instructions assume you are using your home directory; you can put it anywhere so long as you modify the config.ini file, and environment variables accordingly.
```
$ mkdir seqrepo
$ seqrepo --root-directory ~/seqrepo pull -i 2018-08-21
```
To check it has downloaded:
```
$ seqrepo --root-directory ~/seqrepo list-local-instances
```

## Configuration

See the [manual](MANUAL.md) for configuration instructions.
