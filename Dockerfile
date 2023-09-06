FROM python:3.10

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container's /app directory
COPY . /app

# Create logging directory
RUN mkdir /usr/local/share/logs

# Update apt-get
RUN apt-get update

# Install git
RUN apt-get -y install git

# Updrade pip
RUN pip install --upgrade pip

# Install the tool
RUN pip install -e .

# Copy the config file into the container home directory
COPY configuration/docker.ini /root/.variantvalidator

# Define the entrypoint as an empty command
ENTRYPOINT []

# Start the container with CMD and set Gunicorn timeout to 600 seconds
CMD ["gunicorn", "-b", "0.0.0.0:8000", "--timeout", "600", "app", "--threads=5", "--chdir", "./rest_VariantValidator/"]
