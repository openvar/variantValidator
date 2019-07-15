# VariantValidator Installation

These instructions will allow you to configure the software on Linux and Mac OS X computers.

There are several steps involved in setting up VariantValidator:
* The python environment must be set up with the correct packages
* The variantValidator files themselves must installed.
* The databases must be downloaded and set up
* The configuration files must be changed to point VariantValidator at those databases.

## Download the source code

To download the VariantValidator source code simply clone the master branch.

```
$ git clone https://github.com/openvar/variantValidator.git
$ cd variantValidator/
```


## Virtual environment (Python 2.7)

When installing VariantValidator we recommend using a virtual environment, as it requires specific versions of several libraries including python and sqlite. This can be done either via conda or pip.

#### Via conda  
After [installing conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) you can create a new virtual environment with the correct python and sqlite versions by running:
```
$ conda env create -f environment.yml
$ conda activate vvenv
```
The packages required for variant validator to function are now set up in the environment "vvenv".

#### Via pip

If you already have the right versions of python (2.7) and sqlite (>=3.8), then you can use pip to install the remaining packages.

```
$ python -m venv vvenv
$ source activate vvenv
$ pip install -r requirements.txt
```

## Installing VariantValidator

To install VariantValidator within your virtual environment run:
```
$ python setup.py install
```
For development purposes, you can use 
```
$ pip install -e .
```
to ensure any changes you make in the local variant validator folder is reflected in your python site-packages.

## Setting up validator database (MySQL)

A MySQL database is required to run VariantValidator. We recommend creating a user and password specific to the
VariantValidator database.

```mysql
CREATE USER 'vvadmin'@'localhost' IDENTIFIED BY 'var1ant';
CREATE DATABASE validator;
GRANT SELECT,INSERT,UPDATE,DELETE ON validator.* TO 'vvadmin'@'localhost';
```

You can then use either our pre-populated database, or create an empty database that will fill as VariantValidator runs. Note this latter option.
will make the library slower and may return empty values if there is a network connection error.

***We strongly recommend that you download and extract the pre-populated database***

#### Downloading the pre-populated database 

The database is available for download from [figshare](https://doi.org/10.25392/leicester.data.8859317.v1). You can also download the file via the command line:

```bash
wget https://leicester.figshare.com/ndownloader/files/16237784 -O validator_2019-07-10.sql.gz
```

Once downloaded the file needs to be extracted and imported into mysql.

```bash
gunzip validator_2019-07-10.sql.gz
mysql validator < validator_2019-07-10.sql
```

#### Creating an empty database

If you don't wish to use the pre-populated database, in the `VariantValidator/data` folder is a copy of the empty mysql 
database needed by Variant Validator to run. The software will populate it as variants are run. You need to import it into 
MySQL:
```
$ mysql validator < VariantValidator/data/emptyValidatorDump.sql 
```

#### Updating the database contents

The RefSeq and LRG lookup tables may need updating, to do this you'll need to run `bin/update_vdb.py` 
which will download the latest RefSeq data and populate the validator database. 
Note, if you created an empty database you'll need to do this before running Variant Validator.

## Setting up UTA database (PostGreSQL >=9.5)

It's recommended for performance reasons to use a local version of the UTA database. We again recommend creating a specific user account.
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

If you wish to use the remote, public UTA database, see the instructions [here](https://github.com/biocommons/uta#accessing-the-public-uta-instance).

## Setting up Seqrepo (SQLite >=3.8)

VariantValidator requires a local SeqRepo database. The seqrepo library is already installed, but you'll need to download an actual seqrepo database. These instructions assume you are using your home directory; you can put it anywhere so long as you modify the config.ini file, and environment variables accordingly.
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
