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

***Essential step***

Edit the file configuration/docker.ini
You will need to provide an email address and an 
[Entrez API key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)

*Note: configuration can be updated (see below for details)*

## Install and build
*Note: If you have MySQl and or Postgres databases already running, you will need to alter the ports used in the 
docker-comose.yml file. The relevant section is shown here*
```yml
services:
  vdb:
    build:
      context: .
      dockerfile: vdb_docker.df
    ports:
      - "3306:3306"
    expose:
      - "3306"
  uta:
    build:
      context: .
      dockerfile: uta_docker.df
    ports:
      - "5432:5432"
    expose:
      - "5432"
``` 

*Note: some of these steps take >>1hr to complete depending on the speed of your internet connection, particularly 
compiling SeqRepo*



```bash
# Pull images
$ docker-compose pull

# Build
$ docker-compose build --no-cache

# Build and load vv and databases
# This step can take >>1hour and is complete when you see the message
# - "variantvalidator_seqrepo_1 exited with code 0"
$ docker-compose up

# Shutdown
ctrl + c
```

## Launch
You can then launch the docker containers and run them using

```bash
$ docker-compose up
```

Note, the first time this is run it will download each of the databases including the pre-populated
validator database and could take up >1hr depending on your connection. We do not recommend
running this in the background as you need to see the logs and therefore when the databases are
ready to be used.

Once installed and running it is possible to run just the container containing VariantValidator, either to 
run the validator script

```bash
docker-compose run vv variant_validator.py
```

run python

```bash
docker-compose run vv python
```

or go into the container via bash

```bash
docker-compose run vv bash
```

Note, that each time one of these commands is run a new container is created. 
For more information on how to use docker-compose see their [documentation](https://docs.docker.com/compose/).

It is possible to access both the UTA and Validator databases outside of docker as they expose the
 default PostgreSQL and MySQL ports (5432 and 3306 respectively). In the current set-up it is not possible to 
 access the seqrepo database outside of docker.
 
Finally, it should be noted that the current UTA docker container is not up-to-date and only contains the 
2017-10-26 release. Therefore use caution when interpreting these results, and be advised the
 VariantValidator tests will fail. 