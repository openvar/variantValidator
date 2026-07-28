# Declare the base image
FROM python:3.12.11

# Create the WorkDir
WORKDIR /app

# Copy the repository into the container
COPY . /app

# Create logging directory
RUN mkdir -p /usr/local/share/logs

# Update apt-get and install system dependencies
RUN apt update && apt -y install \
    git \
    postgresql-client \
    sqlite3 \
    php

# Upgrade pip
RUN pip install --upgrade pip

# Install VariantValidator
RUN pip install ./packaging/variantvalidator

# Install VariantFormatter
RUN pip install ./packaging/variantformatter

# Copy the config file into the container home directory
COPY configuration/docker.ini /root/.variantvalidator

# Set up the LOVD HGVS Syntax Checker
RUN python -m VariantValidator.bin.setup_lovd_syntax_checker

# Set entrypoint
ENTRYPOINT []

# Keep the container running
CMD ["tail", "-f", "/dev/null"]