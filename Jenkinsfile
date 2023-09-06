pipeline {
    agent {
        docker {
            image "docker"
        }
    }
    environment {
        CODECOV_TOKEN = "50dd5c2e-4259-4cfa-97a7-b4429e0d179e"
        CONTAINER_SUFFIX = "${BUILD_NUMBER}"
    }
    stages {
        stage("Clone Repository") {
            steps {
                checkout scm
            }
        }
        stage("Build VVTA PostgreSQL") {
            steps {
                script {
                    docker.image("postgres:14.7").withRun("-p 5432:5432 -d --name postgres-vvta-${CONTAINER_SUFFIX}") { c ->
                        // Container-specific commands for PostgreSQL setup
                        sh 'apt-get update'
                        sh 'apt-get install -y wget'
                        sh 'echo "shared_buffers = 2GB" > /docker-entrypoint-initdb.d/postgresql.conf'
                        sh 'wget https://www528.lamp.le.ac.uk/vvdata/vvta/vvta_2023_05_no_seq.sql.gz -O input_file.sql.gz'
                        sh 'gzip -dq input_file.sql.gz'
                        sh 'sed "s/anyarray/anycompatiblearray/g" input_file.sql > modified_file.sql'
                        sh 'rm input_file.sql'
                        sh 'gzip modified_file.sql'
                        sh 'mv modified_file.sql.gz /docker-entrypoint-initdb.d/vvta_2023_05_noseq.sql.gz'
                    }
                }
            }
        }
        stage("Build Validator MySQL") {
            steps {
                script {
                    docker.image("ubuntu/mysql:8.0-22.04_beta").withRun("-p 3306:33306 -d --name mysql-validator-${CONTAINER_SUFFIX}") { c ->
                        // Container-specific commands for MySQL setup
                        sh 'apt-get update && apt-get install -y wget'
                        sh 'wget https://www528.lamp.le.ac.uk/vvdata/validator/validator_2023_08.sql.gz -O /docker-entrypoint-initdb.d/validator_2023_08.sql.gz'
                    }
                }
            }
        }
        stage("Start Containers") {
            steps {
                sh 'docker start postgres-vvta-${CONTAINER_SUFFIX}'
                sh 'docker start mysql-validator-${CONTAINER_SUFFIX}'
            }
        }
        stage("Build SeqRepo") {
            agent {
                docker {
                    image "ubuntu:22.04"
                    args '-v $WORKSPACE:/workspace'
                }
            }
            steps {
                sh 'apt-get update'
                sh 'apt-get install -y wget'
                sh 'mkdir -p /workspace/seqrepo'
                sh 'wget --output-document=/workspace/seqrepo/VV_SR_2023_05.tar https://www528.lamp.le.ac.uk/vvdata/vv_seqrepo/VV_SR_2023_05.tar'
                sh 'tar -xvf /workspace/seqrepo/VV_SR_2023_05.tar --directory /workspace/seqrepo'
                sh 'rm /workspace/seqrepo/VV_SR_2023_05.tar'
            }
        }
        stage("Build and Test VariantValidator") {
            agent {
                docker {
                    image "python:3.10"
                    args "-p 3306:3306 -p 5432:5432"
                }
            }
            steps {
                sh 'pip install --upgrade pip'
                sh 'pip install .'
                sh 'cp configuration/continuous_integration.ini "$HOME"/.variantvalidator'
                sh 'pytest --cov-report=term --cov=VariantValidator/'
                sh 'codecov'
            }
        }
        stage("Cleanup Docker") {
            steps {
                sh 'docker stop postgres-vvta-${CONTAINER_SUFFIX}'
                sh 'docker rm postgres-vvta-${CONTAINER_SUFFIX}'
                sh 'docker stop mysql-validator-${CONTAINER_SUFFIX}'
                sh 'docker rm mysql-validator-${CONTAINER_SUFFIX}'
                sh 'docker rmi postgres:14.7'
                sh 'docker rmi ubuntu/mysql:8.0-22.04_beta'
            }
        }
    }
}
