# Docker Installation

VariantValidator supports two Docker-based installation methods.

| Method                       | VariantValidator | Supporting Services |
| ---------------------------- | ---------------- | ------------------- |
| **Docker Quick Start**       | Runs locally     | Runs in Docker      |
| **Full Docker Installation** | Runs in Docker   | Runs in Docker      |

The **Docker Quick Start** runs the required backend services (Validator database, VVTA, and SeqRepo) in Docker while VariantValidator itself is installed and executed locally within a Conda environment. This is the recommended installation method for users who wish to run VariantValidator directly on their local machine while Docker manages the supporting infrastructure.

The **Full Docker Installation** runs both VariantValidator and all supporting services inside Docker, providing a completely containerised environment suitable for isolated deployments, reproducible testing, and continuous integration.

> **Note:** Docker images are only built during the initial installation. Once created, the containers can be started and stopped without rebuilding or reinstalling VariantValidator.

---

# Docker Quick Start

The Quick Start installation runs the required backend services in Docker:

* Validator MySQL database
* VVTA PostgreSQL database
* SeqRepo

VariantValidator itself is installed and executed locally within a Conda environment.

## Requirements

* Docker
* Conda (Miniconda or Anaconda)
* Git

## Clone the repository

```bash
git clone https://github.com/openvar/VariantValidator.git
cd VariantValidator
```

## Create a Docker network

```bash
docker network create variantvalidator-network
```

## Build and start the VVTA PostgreSQL database

```bash
docker build \
    --no-cache \
    -f db_dockerfiles/vvta/Dockerfile \
    -t postgres-vvta \
    db_dockerfiles/vvta

docker run -d \
    --name vv-vvta \
    --network variantvalidator-network \
    --shm-size=2g \
    -p 54320:5432 \
    postgres-vvta
```

## Build and start the Validator MySQL database

```bash
docker build \
    --no-cache \
    -f db_dockerfiles/vdb/Dockerfile \
    -t mysql-validator \
    db_dockerfiles/vdb

docker run -d \
    --name vv-vdb \
    --network variantvalidator-network \
    -p 33060:3306 \
    mysql-validator
```

## Build the SeqRepo image

```bash
docker build \
    --no-cache \
    -f db_dockerfiles/vvsr/Dockerfile \
    -t seqrepo-validator \
    db_dockerfiles/vvsr
```

## Extract the SeqRepo data

Create a temporary container and copy the SeqRepo data to the host.

```bash
sudo mkdir -p /usr/local/share/seqdata

docker create \
    --name vv-seqrepo-temp \
    seqrepo-validator

sudo docker cp \
    vv-seqrepo-temp:/usr/local/share/seqdata/. \
    /usr/local/share/seqdata

docker rm vv-seqrepo-temp
```

## Verify the SeqRepo installation

The directory should contain the extracted SeqRepo data.

```bash
ls /usr/local/share/seqdata
```

## Wait for the databases to initialise

Wait until the database containers have completed their initialisation before continuing.

### Wait for MySQL

Run this command until you see the message "MySQL ready."

```bash
until docker exec vv-vdb \
    mysqladmin ping \
    -u vvadmin \
    -pvar1ant \
    --silent
do
    sleep 5
done

echo "MySQL ready."
```

### Wait for VVTA

Run this command until you see the message "VVTA ready."

```bash
until docker logs vv-vvta 2>&1 | \
    grep -q "PostgreSQL init process complete; ready for start up."
do
    sleep 10
done

echo "VVTA ready."
```

## Create the Conda environment

```bash
conda env create -f environment.yml
```

## Activate the environment

```bash
conda activate vvenv
```

## Configure VariantValidator

```bash
cp configuration/docker-local.ini ~/.variantvalidator
```

## Install VariantValidator

```bash
pip install .
```

## Install the LOVD syntax checker

```bash
python -m VariantValidator.bin.setup_lovd_syntax_checker
```

VariantValidator is now installed and ready to use.

## Test the installation

Run the test suite locally against the Docker-hosted backend services:

```bash
pytest \
    -n 4 \
    tests
```

---

# Full Docker Installation

The Full Docker Installation runs both VariantValidator and its supporting services entirely within Docker.

## Build and start the supporting services

Complete the following sections from the **Docker Quick Start**:

* Create a Docker network
* Build and start the VVTA PostgreSQL database
* Build and start the Validator MySQL database
* Build and start SeqRepo
* Wait for the databases to initialise

## Build the VariantValidator image

```bash
docker build \
    --no-cache \
    -t variantvalidator \
    .
```

## Start VariantValidator

```bash
docker run -d \
    --name variantvalidator \
    --network variantvalidator-network \
    --volumes-from vv-seqrepo:ro \
    variantvalidator
```

## Verify the installation

Run the test suite directly:

```bash
docker exec variantvalidator \
    pytest \
    -n 4 \
    tests
```

Alternatively, open an interactive shell within the container:

```bash
docker exec -it variantvalidator bash
```

Then, from within the container, run:

```bash
pytest \
    -n 4 \
    tests
```

---

# Restarting an Existing Installation

Once the Docker images and containers have been created, they can be restarted without repeating the installation steps.

## Docker Quick Start

Start the supporting services:

```bash
docker start vv-vvta vv-vdb vv-seqrepo
```

Wait a few moments for the database services to initialise, then activate the Conda environment:

```bash
conda activate vvenv
```

VariantValidator is now ready to use.

## Full Docker Installation

Start all containers:

```bash
docker start vv-vvta vv-vdb vv-seqrepo variantvalidator
```

Wait a few moments for the services to initialise.

VariantValidator is now ready to use.

---

# Stopping the Services

## Docker Quick Start

```bash
docker stop vv-vvta vv-vdb vv-seqrepo
```

## Full Docker Installation

```bash
docker stop variantvalidator vv-seqrepo vv-vdb vv-vvta
```

---

# Cleaning Up

To completely remove the installation, stop and remove the containers, then remove the Docker network.

```bash
docker stop variantvalidator vv-seqrepo vv-vdb vv-vvta || true

docker rm variantvalidator vv-seqrepo vv-vdb vv-vvta || true

docker network rm variantvalidator-network || true
```

The Docker images are retained and can be reused the next time you start VariantValidator. To completely remove the installation, including the Docker images, run:

```bash
docker rmi variantvalidator seqrepo-validator mysql-validator postgres-vvta
```
