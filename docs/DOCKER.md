# Docker

To install VariantValidator via Docker, first ensure you have both docker and docker-compose installed. 
See their [documentation](https://docs.docker.com/compose/install/) for information.

Create a directory collate your cloned repositories. Move into the directory then, clone the repository. 

```bash
$ git clone https://github.com/openvar/variantValidator
```

Once the repository has been cloned, cd into the variantValidator directory that the clone creates.
```bash
$ cd variantValidator/
``` 

If you have cloned the repository previously, update it

```bash
$ git pull
```

## Configure

Edit the file configuration/docker.ini as required

From version 2.0.0 adding an API key is optional.
[Entrez API key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)
As is adding an email address

*Note: configuration can be updated (see below for details)*

## Build the container

*Note: some of these steps take >>1hr to complete depending on the speed of your internet connection, particularly 
compiling SeqRepo*

*Note: Depending on your system setup you may need to use `sudo` or be root to run docker. If you are using `sudo`
you will need to prefix the `docker-compose` commands below with `sudo --preserve-env=HOME`, or else if just
using `sudo` edit the `docker-comose.yml` file to replace `${HOME}` with a location of your choice, making sure to
create the `variantvalidator_data` and `share` folders as needed.*

- Pull images

```bash
$ docker-compose pull
```

- Create a directory for sharing resources between your computer and the container
```bash
$ mkdir ~/variantvalidator_data
$ mkdir ~/variantvalidator_data/share
```
i.e. a directory called share in your home directory

- Edit the `vdb_docker.df` file 

You need to select your chip set e.g. Arm or Intel and remove the relevant hash. Default is intel

```
# For Arm chips e.g. Apple M1
# FROM biarms/mysql:5.7

# For Intel chips
FROM mysql:5.7-debian
```

- Build

```bash
$ docker-compose build --no-cache
```

- Complete build
    - The first time you do this, it will complete the build process, for example, populating the required the databases
    - When this is completed you will need to shutdown the services and re-start (see below)
    - The build takes a while because the  vv databases are large. However, this is a significant improvement on previou
    s versions. Install time is approximately 30 minutes (depending on the speed of you computer and internet connection)
    - The build has completed when you see the message ***"Successfully built <container number string>"***
    - example: "Successfully built fc9b83c8d21fa8bdebd52e0e87b9fde967933a043dace1a31916f8106110c8d8
"
    - Then complete the following steps
```bash
# Create the containers (This only takes a coule of minutes)
$ docker-compose up

# When you see the following message the containers have been created. 
# "vvta_1     | 2021-07-23 16:29:17.590 UTC [1] LOG:  database system is ready to accept connections"
# Initial shut down prior to re-launch and working with VarinatValidator in Docker
ctrl + c
```

### Build errors you may encounter

***If other services on your system conflict with the network resources used to expose the docker versions of MySQL
and or Postgres databases you may encounter an error***

> "ERROR: for vdb Cannot start service vdb: Ports are not available: listen tcp 0.0.0.0:53306: bind: address already in use"

or

> "ERROR: for vvta Cannot start service vvta: Ports are not available: listen tcp 0.0.0.0:54321: bind: address already in use"

If you encounter either of these issues, stop the build by pressing `ctrl+c`

- Reconfigure the ports used in the `docker-comose.yml` file
    - If you do not intend to use the databases except to run VariantValidator you can just remove the 'ports:' and
      'expose:' sections containing the offending port number
    - Otherwise you can fix the problem by moving to an unused port, via advancing the port number mentioned in the
      error message by 1
    - you are recommended to hash (`#`) the conflicting port first, then add the new ports if wanted by advancing the
      external (first) port numbers, as shown below
```yml
services:
  vdb:
    build:
      context: .
      dockerfile: vdb_docker.df
    ports:
      # - "53306:3306"
      - "53307:3306"
    expose:
      # - "53306"
      - "53307"
```
or
```yml
  vvta:
    build:
      context: .
      dockerfile: uta_docker.df
    ports:
      # - "54321:5432"
      - "54322:5432"
    expose:
      # - "54321"
      - "54322"
``` 

- save the file once you have finished editing

- force-recreate the container to apply the new settings, by running the following docker commands

```bash
$ docker-compose down
$ docker-compose up --force-recreate
```

## Checking the installation
go into the container via bash

```bash
$ docker-compose run vv bash
```

Use Pytest to check the integrity of the installation (Recommended)

- First check that the SeqRepo has installed correctly

```
$ ls /usr/local/share/seqrepo/
# returns
VV_SR_2021_2
``` 

```bash
$ cd /app
$ pytest
```

## Launch
You can then launch the docker containers and run them using

```bash
$ docker-compose up
```

Once installed and running it is possible to run just the container containing VariantValidator, either to 
run the validator script

```bash
$ docker-compose run vv variant_validator.py
```

**Example**
```bash
# Note: The variant description must be contained in '' or "". See MANUAL.md for more examples
$ docker-compose run vv variant_validator.py -v 'NC_000017.11:g.50198002C>A' -g GRCh38 -t mane -s individual -f json -m -o stdout
```

**Example 2 - use Python to collect output**
```python
import subprocess
validation = subprocess.run(["docker-compose run vv variant_validator.py -v 'NC_000017.11:g.50198002C>A' -g GRCh38 -t mane -s individual -f json -m -o stdout"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
print(validation.stdout.decode("utf-8"))
```

run python

```bash
$ docker-compose run vv python
```

or go into the container via bash

```bash
$ docker-compose run vv bash
```

Note, that each time one of these commands is run a new container is created. 
For more information on how to use docker-compose see their [documentation](https://docs.docker.com/compose/).

## Accessing the VariantValidator databases externally
It is possible to access both the UTA and Validator databases outside of docker as they expose the
 default PostgreSQL and MySQL ports (5432 and 3306 respectively). It is also possible to 
 access the seqrepo database outside of docker by editing your config file to point at the shared directory.
 
```
~/variantvalidator_data/share/seqrepo
```
 

## Accessing VariantValidator directly through bash and reconfiguring a container post build
The container hosts a full install of VariantValidator. 

To start this version you use the command

```bash
$ docker-compose run vv bash
```

When you are finished exit the container

```bash
$ exit
```

#### What you can do in bash mode

1. Run pytest on VariantValidator, to test the function of your VariantValidator install

2. Run VariantValidator, which can be run on the commandline from within the container
    - Instructions can be found in the VariantValidator [manual](https://github.com/openvar/variantValidator/blob/master/docs/MANUAL.md) under sections **Database updates** and **Operation**


## Developing VariantValidator in Docker
The container has been configured with git installed. This means that you can clone Repos directly into the container

To develop VariantValidator in the container

Start the container 

```bash
$ docker-compose run vv bash
```

ON YOUR COMPUTER change into the share directory

```bash
$ cd ~/share
```

Then create a directory for development

```bash
$ mkdir DevelopmentRepos
$ cd ~/share/DevelopmentRepos
```

Clone the VariantValidator Repo

```bash
$ git clone https://github.com/openvar/variantValidator.git
```

Checkout the develop branch

```bash
$ git checkout develop
$ git pull
```

Create an new branch for your developments

```bash
$ git branch name_of_branch
$ git checkout name_of_branch
```

IN THE CONTAINER, pip install the code so it can be run by the container

```bash
$ cd /usr/local/share/DevelopmentRepos/variantValidator
$ pip install -e . 
```

You can then use the containers Python interpreter to run queries, e.g.

```python
import json
import VariantValidator
vval = VariantValidator.Validator()
variant = 'NM_000088.3:c.589G>T'
genome_build = 'GRCh38'
select_transcripts = 'all'
validate = vval.validate(variant, genome_build, select_transcripts)
validation = validate.format_as_dict(with_meta=True)
print(json.dumps(validation, sort_keys=True, indent=4, separators=(',', ': ')))
```

## Updating variantValidator using docker-compose
Update requires that the vv container is deleted from your system. This is not achieved by removing the container

If you are only running variantValidator in docker, we recommend deleting and re-building all containers

```bash
# Delete all containers
$ docker-compose down
$ docker system prune -a --volumes
```

***Once you have deleted the containers, got to Install and Build***

Alternatively, you may wish to try and force the containers to re-build without deleting

```bash
# Force re-build
$ docker-compose down
$ docker-compose up --force-recreate
```

***If you choose this option, make sure you see the container vv being re-created and all Python packages being
reinstalled in the printed logs, otherwise the container may not actually be rebuilt and the contained modules may not
 update***
