# Should run on all processors
FROM ubuntu/mysql:8.0-22.04_beta

ENV MYSQL_RANDOM_ROOT_PASSWORD yes
ENV MYSQL_DATABASE validator
ENV MYSQL_USER vvadmin
ENV MYSQL_PASSWORD var1ant

RUN apt-get update && apt-get install -y \
    wget

RUN rm -rf /var/lib/apt/lists/*

# Set the max_connections directly in the my.cnf
RUN echo '[mysqld]' >> /etc/mysql/my.cnf && \
    echo 'max_connections=250' >> /etc/mysql/my.cnf

RUN wget https://www528.lamp.le.ac.uk/vvdata/validator/validator_2023_08.sql.gz -O /docker-entrypoint-initdb.d/validator_2023_08.sql.gz

CMD ["mysqld"]
