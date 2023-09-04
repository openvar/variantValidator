pipeline {
    agent {
        dockerContainer {
            image 'python:3.10'
        }
    }

    environment {
        CODECOV_TOKEN = '50dd5c2e-4259-4cfa-97a7-b4429e0d179e'
    }

    stages {
        stage("Build VVTA") {
            steps {
                script {
                    sh '-x docker network create vvta_network'

                    // Start the PostgreSQL container
                    sh '-x docker run --name postgres-vvta -d --network vvta_network -e POSTGRES_DB=vvta -e POSTGRES_USER=uta_admin -e POSTGRES_PASSWORD=uta_admin -p 5432:5432 postgres:14.7'
                    sh '-x sleep 10'
                    sh '-x echo "shared_buffers = 2GB" > postgres_config.conf'
                    sh '-x wget https://www528.lamp.le.ac.uk/vvdata/vvta/vvta_2023_05_no_seq.sql.gz -O input_file.sql.gz'
                    sh '-x gzip -dq input_file.sql.gz'
                    sh '-x sed "s/anyarray/anycompatiblearray/g" input_file.sql > modified_file.sql'
                    sh '-x gzip modified_file.sql'
                    sh '-x docker cp postgres_config.conf postgres-vvta:/docker-entrypoint-initdb.d/postgresql.conf'
                    sh '-x docker cp modified_file.sql.gz postgres-vvta:/docker-entrypoint-initdb.d/vvta_2023_05_noseq.sql.gz'
                    sh '-x rm postgres_config.conf input_file.sql.gz modified_file.sql.gz'
                    sh '-x docker restart postgres-vvta'

                    // Expose PostgreSQL container port
                    sh '-x docker network connect bridge postgres-vvta'
                }
            }
        }

        stage("Build Validator") {
            steps {
                script {
                    // Start the MySQL container
                    sh '-x docker run --name mysql-validator -d --network vvta_network -e MYSQL_RANDOM_ROOT_PASSWORD=yes -e MYSQL_DATABASE=validator -e MYSQL_USER=vvadmin -e MYSQL_PASSWORD=var1ant -p 3306:3306 ubuntu/mysql:8.0-22.04_beta'
                    sh '-x sleep 10'
                    sh '-x docker exec mysql-validator wget https://www528.lamp.le.ac.uk/vvdata/validator/validator_2023_08.sql.gz -O /docker-entrypoint-initdb.d/validator_2023_08.sql.gz'

                    // Expose MySQL container port
                    sh '-x docker network connect bridge mysql-validator'
                }
            }
        }

        stage("Before Install") {
            steps {
                sh '-x apt-get update'
                sh '-x apt-get install -y wget'
            }
        }

        stage("Install") {
            steps {
                script {
                    sh '-x mkdir -p /usr/local/share/seqrepo'
                    sh '-x mkdir -p /usr/local/share/logs'
                    sh '-x wget --output-document="/usr/local/share/seqrepo/VV_SR_2023_05.tar" https://www528.lamp.le.ac.uk/vvdata/vv_seqrepo/VV_SR_2023_05.tar'
                    sh '-x tar -xvf /usr/local/share/seqrepo/VV_SR_2023_05.tar -C /usr/local/share/seqrepo/'
                    sh '-x pip install --upgrade pip'
                    sh '-x pip install .'
                    sh '-x cp configuration/continuous_integration.ini "$HOME"/.variantvalidator'
                }
            }
        }

        stage("Test") {
            steps {
                script {
                    sh '-x pytest --cov-report=term --cov=VariantValidator/'
                }
            }
        }

        stage("After Script") {
            steps {
                script {
                    sh '-x codecov'
                }
            }
        }
    }
}
