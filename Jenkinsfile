pipeline {
    agent {
        docker {
            image "docker"
        }
    }
    environment {
        CODECOV_TOKEN = "50dd5c2e-4259-4cfa-97a7-b4429e0d179e"
        CONTAINER_SUFFIX = "${BUILD_NUMBER}"
        DOCKER_NETWORK = "variantvalidator_docker_network-$CONTAINER_SUFFIX"
        DATA_VOLUME = "docker-shared-space"
    }
    stages {
        stage("Clone Repository Romove dangling docker components and Create Docker Network") {
            steps {
                checkout scm
                sh 'docker system prune -f'
                sh 'docker network create $DOCKER_NETWORK'
            }
        }
        stage("Build and Run VVTA PostgreSQL") {
            steps {
                script {
                    def dockerfile = './db_dockerfiles/vvta/Dockerfile'
                    def vvtaContainer = docker.build("postgres-vvta-${CONTAINER_SUFFIX}", "--no-cache -f ${dockerfile} ./db_dockerfiles/vvta")
                    vvtaContainer.run("-p 5432:5432 -d --name vv-vvta --network $DOCKER_NETWORK --shm-size=2g")
                    sh 'echo Building and running VVTA PostgreSQL'
                }
            }
        }
        stage("Build and Run Validator MySQL") {
            steps {
                script {
                    def dockerfile = './db_dockerfiles/vdb/Dockerfile'
                    def validatorContainer = docker.build("mysql-validator-${CONTAINER_SUFFIX}", "--no-cache -f ${dockerfile} ./db_dockerfiles/vdb")
                    validatorContainer.run("-p 3306:3306 -d --name vv-vdb --network $DOCKER_NETWORK")
                    sh 'echo Building and running Validator MySQL'
                }
            }
        }
        stage("Build and Run SeqRepo") {
            steps {
                script {
                    def dockerfile = './db_dockerfiles/vvsr/Dockerfile'
                    def seqRepoContainer = docker.build("sqlite-seqrepo-${CONTAINER_SUFFIX}", "--no-cache -f ${dockerfile} ./db_dockerfiles/vvsr")
                    seqRepoContainer.run("--network $DOCKER_NETWORK --name vv-seqrepo -v $DATA_VOLUME:/usr/local/share:rw")
                    sh 'echo Building and running SeqRepo'
                }
            }
        }
        stage("Build and Run VariantValidator") {
            steps {
                script {
                    def dockerfile = './Dockerfile'
                    def variantValidatorContainer = docker.build("variantvalidator-${CONTAINER_SUFFIX}", "--no-cache -f ${dockerfile} .")
                    variantValidatorContainer.run("-v $DATA_VOLUME:/usr/local/share:rw -d --name variantvalidator-${CONTAINER_SUFFIX} --network $DOCKER_NETWORK")
                    sh 'echo Building and running VariantValidator'
                }
            }
        }
        stage("Run Pytest and Codecov") {
            steps {
                sh 'docker ps'
                sh 'docker exec variantvalidator-${CONTAINER_SUFFIX} pytest --cov-report=term --cov=VariantValidator/'
                sh 'docker exec variantvalidator-${CONTAINER_SUFFIX} codecov -t $CODECOV_TOKEN -b ${BRANCH_NAME}'
            }
        }
        stage("Cleanup Docker") {
            steps {
                sh 'docker stop vv-vvta'
                sh 'docker rm vv-vvta'
                sh 'docker stop vv-vdb'
                sh 'docker rm vv-vdb'
                sh 'docker stop vv-seqrepo'
                sh 'docker rm vv-seqrepo'
                sh 'docker stop variantvalidator-${CONTAINER_SUFFIX}'
                sh 'docker rm variantvalidator-${CONTAINER_SUFFIX}'
                sh 'docker network rm $DOCKER_NETWORK'
            }
        }
    }
    post {
        success {
            echo 'Pipeline succeeded! Your project is built and tested.'
            emailext(
                subject: 'Pipeline Success',
                body: 'Your Jenkins pipeline has succeeded. Your project is built and tested successfully.',
                recipientProviders: [[$class: 'CulpritsRecipientProvider']],
                to: 'admin@variantvalidator.org'
            )
        }
        failure {
            echo 'Pipeline failed. Please check the logs for details.'
            emailext(
                subject: 'Pipeline Failure',
                body: 'Your Jenkins pipeline has failed. Please check the logs for details.',
                recipientProviders: [[$class: 'CulpritsRecipientProvider']],
                to: 'admin@variantvalidator.org'
            )
        }
    }
}
