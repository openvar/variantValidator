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

1. Run VariantValidator can be run on the commandline from within the container
    - Instructions can be found in the VariantValidator [manual](https://github.com/openvar/variantValidator/blob/master/docs/MANUAL.md) under sections **Database updates** and **Operation**
    
2. Start the REST services in development mode, bound to port 5000 
    - For example, this is useful if you want to develop new methods and test them
    - Note: Under the terms and conditions of our [license](https://github.com/openvar/rest_variantValidator/blob/master/LICENSE.txt) changes to the code and improvements must be made available to the community so that we can integrate them for the good of all our users 
    - See instructions on VariantValidator development in Docker 


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

## Updating rest_variantValidator using docker-compose
Update requires that the restvv container is deleted from your system. This is not achieved by removing the container

If you are only running rest_variantValidator in docker, we recommend deleting and re-building all containers

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

***If you choose this option, make sure you see the container restvv being re-created and all Python packages being 
reinstalled in the printed logs, otherwise the container may not actually be rebuilt and the contained modules may not
 update***