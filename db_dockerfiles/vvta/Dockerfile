FROM postgres:14.7

ENV POSTGRES_DB=vvta
ENV POSTGRES_USER=uta_admin
ENV POSTGRES_PASSWORD=uta_admin

# Install necessary dependencies
RUN apt-get update && \
    apt-get install -y \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Auto-create and set PostgreSQL configuration
RUN echo "shared_buffers = 2GB" > /docker-entrypoint-initdb.d/postgresql.conf

# Step 1: Download the file
RUN wget https://www528.lamp.le.ac.uk/vvdata/vvta/vvta_2023_05_no_seq.sql.gz -O input_file.sql.gz

# Step 2: Extract the gzipped file
RUN gzip -dq input_file.sql.gz

# Step 3: Use sed to replace text
RUN sed 's/anyarray/anycompatiblearray/g' input_file.sql > modified_file.sql

# Step 4: Compress the modified file
RUN rm input_file.sql
RUN gzip modified_file.sql

# Step 5: Move the modified and compressed file to the desired location
RUN mv modified_file.sql.gz /docker-entrypoint-initdb.d/vvta_2023_05_noseq.sql.gz

