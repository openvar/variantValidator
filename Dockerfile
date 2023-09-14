FROM python:3.10

WORKDIR /app

COPY . /app

# Update apt-get
RUN apt update

# Install apt managed sofware
RUN apt -y install git \
    postgresql-client \
    sqlite3

# Upgrade pip
RUN pip install --upgrade pip

RUN pip install .

COPY configuration/docker.ini /root/.variantvalidator

ENTRYPOINT []

CMD ["tail", "-f", "/dev/null"]
