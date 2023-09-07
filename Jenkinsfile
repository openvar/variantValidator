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
    }
    stages {
        stage("Clone Repository and Create Docker Network") {
            steps {
                checkout scm
                sh 'docker network create $DOCKER_NETWORK'
            }
        }
        stage("Create Directories on Host") {
            steps {
                sh 'mkdir /var/jenkins_home/workspace/VariantValidator_ci/variantvalidator_data'
                sh 'mkdir /var/jenkins_home/workspace/VariantValidator_ci/variantvalidator_data/share'
                sh 'mkdir /var/jenkins_home/workspace/VariantValidator_ci/variantvalidator_data/share/seqrepo/'
                sh 'mkdir /var/jenkins_home/workspace/VariantValidator_ci/variantvalidator_data/share/logs'
            }
        }
        stage("Where am I") {
            steps {
                sh 'pwd'
                sh 'ls -l'
            }
        }
        stage("Build and Run VVTA PostgreSQL") {
            steps {
                script {
                    def dockerfile = './db_dockerfiles/vvta/Dockerfile'
                    def vvtaContainer = docker.build("postgres-vvta-${CONTAINER_SUFFIX}", "-f ${dockerfile} ./db_dockerfiles/vvta")
                    vvtaContainer.run("-p 5432:5432 -d --name vvta --network $DOCKER_NETWORK")
                    sh 'echo Building and running VVTA PostgreSQL'
                }
            }
        }
        stage("Build and Run Validator MySQL") {
            steps {
                script {
                    def dockerfile = './db_dockerfiles/vdb/Dockerfile'
                    def validatorContainer = docker.build("mysql-validator-${CONTAINER_SUFFIX}", "-f ${dockerfile} ./db_dockerfiles/vdb")
                    validatorContainer.run("-p 3306:3306 -d --name vdb --network $DOCKER_NETWORK")
                    sh 'echo Building and running Validator MySQL'
                }
            }
        }
        stage("Build and Run SeqRepo") {
            steps {
                script {
                    def dockerfile = './db_dockerfiles/vvsr/Dockerfile'
                    def seqRepoContainer = docker.build("sqlite-seqrepo-${CONTAINER_SUFFIX}", "-f ${dockerfile} ./db_dockerfiles/vvsr")
                    seqRepoContainer.run("--network $DOCKER_NETWORK -v /var/jenkins_home/workspace/VariantValidator_ci/variantvalidator_data/share:/usr/local/share -v /var/jenkins_home/workspace/VariantValidator_ci/variantvalidator_data/share/seqrepo:/usr/local/share/seqrepo -v /var/jenkins_home/workspace/VariantValidator_ci/variantvalidator_data/share/logs:/usr/local/share/logs")
                    sh 'echo Building and running SeqRepo'
                }
            }
        }
        stage("Find Seqrepo Mount") {
            steps {
                sh 'pwd'
                sh 'ls -l'
                sh 'ls -l /var/jenkins_home/workspace/VariantValidator_ci/variantvalidator_data/share/seqrepo/'
            }
        }
        stage("Build and Run VariantValidator") {
            steps {
                script {
                    def dockerfile = './Dockerfile'
                    def variantValidatorContainer = docker.build("variantvalidator-${CONTAINER_SUFFIX}", "-f ${dockerfile} .")
                    variantValidatorContainer.run("-v /var/jenkins_home/workspace/VariantValidator_ci/variantvalidator_data/share:/usr/local/share -v /var/jenkins_home/workspace/VariantValidator_ci/variantvalidator_data/share/logs:/usr/local/share/logs -v /var/jenkins_home/workspace/VariantValidator_ci/variantvalidator_data/share/seqrepo:/usr/local/share/seqrepo -d --name variantvalidator-${CONTAINER_SUFFIX} --network $DOCKER_NETWORK")
                    sh 'echo Building and running VariantValidator'
                }
            }
        }
        stage("Run Pytest and Codecov") {
            steps {
                sh 'docker ps'
                sh 'docker exec variantvalidator-${CONTAINER_SUFFIX} pytest --cov-report=term --cov=VariantValidator/'
                sh 'docker exec variantvalidator-${CONTAINER_SUFFIX} codecov'
            }
        }
        stage("Cleanup Docker") {
            steps {
                sh 'docker stop vvta'
                sh 'docker rm vvta'
                sh 'docker stop vdb'
                sh 'docker rm vdb'
                sh 'docker stop sqlite-seqrepo-${CONTAINER_SUFFIX}'
                sh 'docker rm sqlite-seqrepo-${CONTAINER_SUFFIX}'
                sh 'docker stop variantvalidator-${CONTAINER_SUFFIX}'
                sh 'docker rm variantvalidator-${CONTAINER_SUFFIX}'
                sh 'docker network rm $DOCKER_NETWORK'
            }
        }
    }
}
