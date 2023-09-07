FROM python:3.10

WORKDIR /app

COPY . /app

# Update apt-get
RUN apt-get update

# Install git
RUN apt-get -y install git

# Upgrade pip
RUN pip install --upgrade pip

RUN pip install .

COPY configuration/docker.ini /root/.variantvalidator

ENTRYPOINT []

CMD ["/bin/bash"]
