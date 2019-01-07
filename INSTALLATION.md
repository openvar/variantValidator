# Variant Validator installation instructions

In these instructions, lines that must be entered at the command prompt are preceded with >, like so:
> ls

These instructions will allow you to configure the software on Linux. Mac OS X computers operate similarly.

There are several steps involved in setting up variant validator:
* The application files themselves must be installed from SVN.
* The python environment must be set up. On a LAMP, only a custom version of Python will do.
* Protobuf must be compiled and installed
* Required python packages need to be installed, too.
* The databases must be downloaded and set up
* The configuration files must be changed to point the validator at those databases.

## Virtual environment

Variant validator currently requires python 2.7.

When installing Variant Validator it is wise to use a virtual environment, as it requires specific versions of several libraries.
First, download and set up conda (in this case miniconda as we don't need all packages)
 > wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
 > bash Miniconda2-latest-Linux-x86_64.sh 
 > echo ". /local/miniconda2/etc/profile.d/conda.sh" >> ~/.bashrc
 > source ~/.bashrc 
Then create the conda environment and install the necessary programs (this should be done in an environment.yml file eventually). Note, installing biotools downgrades the version of setuptools so that needs to be reinstalled before the pip command to install hgvs=1.1.3
 > conda create -n VVenv
 > conda activate VVenv
 > conda install -c conda-forge sqlite python=2.7 protobuf=3.5.1 docutils python-daemon httplib2 mysql-connector-python mysql-python 
 > conda install -c auto biotools
 > conda install -c bioconda pyliftover pysam
 > conda install setuptools numpy
 > conda install -c anaconda pytest
 > pip install hgvs==1.1.3
The packages required for variant validator to function are now set up in the environment "VVenv".

## Installing validator code

To clone this software from GIT, use:
 > git clone https://github.com/pjdp2/variantValidator.git
This'll create a variantValidator folder in the directory you run it in.
 > cd variantValidator
Run the installation script to integrate variant validator with python's site packages.
 > python setup.py install
For development purposes, you can use
 > pip install -e .
to ensure any changes you make in the local variant validator folder is reflected in your python site-packages.

## Setting up MySQL

This step is not optional for getting variant validator to work. Install packages with:
 > sudo apt-get install mysql-server

This will install everything you need and start the database server. Make sure you note down the root account password that you're prompted for during installation!
Check it runs with:
 > sudo service mysql status
If it's not running, use
 > sudo service mysql start
to boot it up.
Enter mysql from any user's shell prompt with
 > mysql -u root -p
This will prompt you for the root password you made earlier. Within MySQL, create the variant validator user:
 > CREATE USER 'vvadmin'@'localhost' IDENTIFIED BY 'var1ant';
You should create the database too
 > CREATE DATABASE validator;
 > USE validator;
Grant access rights to the vvadmin user:
 > GRANT SELECT,INSERT,UPDATE,DELETE ON validator.* TO vvadmin;
Quit mysql with
 > \q
Bye indeed.

You must source a copy of the validator database from somewhere. That'll have to be fixed for release...
Copy it over to a temporary folder (say, temp, in your home directory).
 > scp someone@somewhere~/databases/validator_2018-11-08.sql ~/temp/validator_2018-11-08.sql
Then, upload it to the running MySQL with:
 > mysql -u root -p validator < ~/databases/validator_2018-11-08.sql 
You should log into MySQL and check to see if the database uploaded correctly. Login with vvadmin, password "var1ant".
Then:
 > USE validator;
 > SHOW TABLES;
which should give some good lines.

## Setting up PostGreSQL

It's recommended for performance reasons to use a local varsion of the UTA database. To do this, first install the required packages with:
 > sudo apt-get install postgresql postgresql-contrib
You need to switch to the "postgres" user to make anything work initially.
 > sudo -i -u postgres
Create a new user with a name matching your user account. In my case - pjdp2. When prompted, make yourself a superuser.
 > createuser --interactive
The postgres user doesn't have a unix password, so you'll need to use exit to get your account back.
 > exit
Enter the database with psql. You'll be signed by default into the "postgres" database, which serves as a kind of master database for controlling user accounts.
 > psql postgres
Inside psql, create the uta_admin role, and set the password when prompted to "uta_admin".
 > CREATE ROLE uta_admin WITH CREATEDB;
 > ALTER ROLE uta_admin WITH LOGIN;
 > \password uta_admin
Create an empty uta database
 > CREATE DATABASE uta WITH OWNER=uta_admin TEMPLATE=template0;
That's enough setting up. Quit psql with:
 > \q
Now you're back to your own prompt, download the gzipped uta genetics database, and upload it into psql. You'll be prompted for your password.
 > wget http://dl.biocommons.org/uta/uta_20180821.pgd.gz
 > gzip -cdq uta_20180821.pgd.gz | psql -U uta_admin -v ON_ERROR_STOP=0 -d uta -Eae
The database should now be uploaded. Don't worry, you can access the database uta with uta_admin if it's uploaded by someone else.
If the database returns errors when the validator runs, you will need to change the postgresql authentication methods, by editing
 > pg_hba.conf 
This file lives, on linux, in /etc/postgresql/9.3/main/pg_hba.conf but on other systems you may need to search for it.
Inside the file, you should change all instances of "peer" to "md5".

## Setting up Seqrepo

Similarly, things run much faster with a local SeqRepo database. You've installed the seqrepo package with pip, but you'll need to download an actual sequence repository. These instructions assume you are using your home directory; you can put it anywhere so long as you modify the config.ini file accordingly.
 > mkdir seqrepo
Then make a cup of tea while this command runs:
 > seqrepo --root-directory ~/seqrepo pull -i 2018-08-21
After it finishes downloading, check it installed correctly:
 > seqrepo --root-directory ~/seqrepo list-local-instances

## Configuration

See the file MANUAL.md for configuration instructions.
