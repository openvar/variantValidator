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

## Build the container

*Note: some of these steps take >>1hr to complete depending on the speed of your internet connection, particularly 
compiling SeqRepo*

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

- Build

```bash
$ docker-compose build --no-cache
```

- Start the container
    - The first time you do this, it will complete the build process, for example, populating the required the databases
    - This step can take >>1hour and is complete when you see the message `rest_variantvalidator_seqrepo_1 exited with code 0"`
    - When this is completed you will need to shutdown the services and re-start (see below)

```bash
$ docker-compose up
```    

- Shutdown services when you want to stop the container

```bash
ctrl + c
```

- Re-start services

```bash
$ docker-compose up
```


### Build errors you may encounter

***If you have MySQL and or Postgres databases already running, you may encounter an error***  

> "ERROR: for vdb  Cannot start service vdb: Ports are not available: listen tcp 0.0.0.0:3306: bind: address already in use" 

If you encounter these issues, stop the build by pressing `ctrl+c`

- Reconfigure the ports used in the `docker-comose.yml` file as shown here
```yml
services:
  vdb:
    build:
      context: .
      dockerfile: vdb_docker.df
    ports:
      # - "33060:3306"
      - "3306:3306"
    expose:
      # - "33060"
      - "3306"
  uta:
    build:
      context: .
      dockerfile: uta_docker.df
    ports:
      - "54320:5432"
    expose:
      - "54320"

``` 
- hash (`#`) the conflicting port and add the new ports as shown above
- force-recreate the container

```bash
$ docker-compose down
$ docker-compose up --force-recreate
```

***You may encounter a build error relating to other unavailable ports***  

> "Cannot start service restvv: Ports are not available: listen tcp 0.0.0.0:8000: bind: address already in use" 

If you encounter these issues, stop the build by pressing `ctrl+c`

- Reconfigure the ports used in the `docker-comose.yml` file as shown here

```yml
  restvv:
    build: .
    depends_on:
      - vdb
      - uta
    volumes:
      - seqdata:/usr/local/share/seqrepo
    ports:
      - "5000:5000"
      # - "8000:8000"
      - "8080:8080"
    expose:
      - "5000"
      # - "8000"
      - 8080
```

- hash (`#`) the conflicting port and add the new ports as shown above
- Change the command in Dockerfile to reflect the changes e.g. `CMD gunicorn  -b 0.0.0.0:8080 app --threads=5 --chdir ./rest_VariantValidator/`
- force-recreate the container

```bash
$ docker-compose down
$ docker-compose up --force-recreate
```

## Accessing the VariantValidator databases externally
It is possible to access both the UTA and Validator databases outside of docker as they expose the
 default PostgreSQL and MySQL ports (5432 and 3306 respectively). In the current set-up it is not possible to 
 access the seqrepo database outside of docker.
 

## Accessing VariantValidator directly through bash and reconfiguring a container post build
The container hosts a full install of VariantValidator. 

To start this version you use the command

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
