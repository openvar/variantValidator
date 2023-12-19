# Declare the base image
FROM python:3.11

# Create the WorkDir
WORKDIR /app

# Copy the current directory contents into the container's /app directory
COPY . /app

# Create logging directory
RUN mkdir /usr/local/share/logs

# Update apt-get
RUN apt update

# Install apt managed sofware
RUN apt -y install git \
    postgresql-client \
    sqlite3

# Upgrade pip
RUN pip install --upgrade pip

# Install the app
RUN pip install -e .

# Copy the config file into the container home directory
COPY configuration/docker.ini /root/.variantvalidator

# Set entrypoint
ENTRYPOINT []

# Set command
CMD ["tail", "-f", "/dev/null"]
