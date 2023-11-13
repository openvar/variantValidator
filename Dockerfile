# Declare the base container
FROM python:3.11

# Create the WorkDir
WORKDIR /app

# Use an environment variable to control whether to mount the volume
ARG USE_VOLUME=false

# Only copy the local code if USE_VOLUME is not set to true
# This line will be skipped if USE_VOLUME=true
COPY . /app/

# Update apt-get
RUN apt update

# Install apt managed software
RUN apt -y install git \
    postgresql-client \
    sqlite3

# Upgrade pip
RUN pip install --upgrade pip

# Install the tool
RUN pip install .

# Copy the config file into the container home directory
COPY configuration/docker.ini /root/.variantvalidator

# Define the entrypoint as an empty command
ENTRYPOINT []

# Start the container with CMD and set Gunicorn timeout to 600 seconds
CMD ["gunicorn", "-b", "0.0.0.0:8000", "--timeout", "600", "app", "--threads=5", "--chdir", "./rest_VariantValidator/"]
