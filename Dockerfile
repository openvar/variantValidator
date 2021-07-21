FROM python:3.6

WORKDIR /app

COPY . /app

# Update apt-get
RUN apt-get update

# Install git
RUN apt-get -y install git

# Updrade pip
RUN pip install --upgrade pip

RUN pip install -r requirements_dev.txt

RUN pip install -e .

COPY configuration/docker.ini /root/.variantvalidator

CMD python3 bin/variant_validator.py
