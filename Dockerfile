FROM python:3.6

WORKDIR /app

COPY . /app

RUN pip install -r requirements_dev.txt

RUN pip install -e .

COPY configuration/docker.ini /root/.variantvalidator

CMD python3 bin/variant_validator.py
