pipeline {
    agent {
        docker { image 'python:3.10' }
    }

    environment {
        CODECOV_TOKEN = '50dd5c2e-4259-4cfa-97a7-b4429e0d179e'
    }

    stages {
        stage("Build VVTA") {
            steps {
                script {
                    sh 'docker network create vvta_network'

                    // Start the PostgreSQL container
                    sh 'docker run --name postgres-vvta -d --network vvta_network -e POSTGRES_DB=vvta -e POSTGRES_USER=uta_admin -e POSTGRES_PASSWORD=uta_admin -p 5432:5432 postgres:14.7'
                    sh 'sleep 10'
                    sh 'echo "shared_buffers = 2GB" > postgres_config.conf'
                    sh 'wget https://www528.lamp.le.ac.uk/vvdata/vvta/vvta_2023_05_no_seq.sql.gz -O input_file.sql.gz'
                    sh 'gzip -dq input_file.sql.gz'
                    sh 'sed "s/anyarray/anycompatiblearray/g" input_file.sql > modified_file.sql'
                    sh 'gzip modified_file.sql'
                    sh 'docker cp postgres_config.conf postgres-vvta:/docker-entrypoint-initdb.d/postgresql.conf'
                    sh 'docker cp modified_file.sql.gz postgres-vvta:/docker-entrypoint-initdb.d/vvta_2023_05_noseq.sql.gz'
                    sh 'rm postgres_config.conf input_file.sql.gz modified_file.sql.gz'
                    sh 'docker restart postgres-vvta'

                    // Expose PostgreSQL container port
                    sh 'docker network connect bridge postgres-vvta'
                }
            }
        }

        stage("Build Validator") {
            steps {
                script {
                    // Start the MySQL container
                    sh 'docker run --name mysql-validator -d --network vvta_network -e MYSQL_RANDOM_ROOT_PASSWORD=yes -e MYSQL_DATABASE=validator -e MYSQL_USER=vvadmin -e MYSQL_PASSWORD=var1ant -p 3306:3306 ubuntu/mysql:8.0-22.04_beta'
                    sh 'sleep 10'
                    sh 'docker exec mysql-validator wget https://www528.lamp.le.ac.uk/vvdata/validator/validator_2023_08.sql.gz -O /docker-entrypoint-initdb.d/validator_2023_08.sql.gz'

                    // Expose MySQL container port
                    sh 'docker network connect bridge mysql-validator'
                }
            }
        }

        stage("Before Install") {
            steps {
                sh 'apt-get update'
                sh 'apt-get install -y wget'
            }
        }

        stage("Install") {
            steps {
                script {
                    sh 'mkdir -p /usr/local/share/seqrepo'
                    sh 'mkdir -p /usr/local/share/logs'
                    sh 'wget --output-document="/usr/local/share/seqrepo/VV_SR_2023_05.tar" https://www528.lamp.le.ac.uk/vvdata/vv_seqrepo/VV_SR_2023_05.tar'
                    sh 'tar -xvf /usr/local/share/seqrepo/VV_SR_2023_05.tar -C /usr/local/share/seqrepo/'
                    sh 'pip install --upgrade pip'
                    sh 'pip install .'
                    sh 'cp configuration/continuous_integration.ini "$HOME"/.variantvalidator'
                }
            }
        }

        stage("Test") {
            steps {
                script {
                    sh 'pytest --cov-report=term --cov=VariantValidator/'
                }
            }
        }

        stage("After Script") {
            steps {
                script {
                    sh 'codecov'
                }
            }
        }
    }
}
