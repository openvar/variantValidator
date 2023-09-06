# Should run on all processors
FROM ubuntu/mysql:8.0-22.04_beta

ENV MYSQL_RANDOM_ROOT_PASSWORD yes
ENV MYSQL_DATABASE validator
ENV MYSQL_USER vvadmin
ENV MYSQL_PASSWORD var1ant

RUN apt-get update && apt-get install -y \
    wget

RUN rm -rf /var/lib/apt/lists/*

RUN wget https://www528.lamp.le.ac.uk/vvdata/validator/validator_2023_08.sql.gz -O /docker-entrypoint-initdb.d/validator_2023_08.sql.gz

