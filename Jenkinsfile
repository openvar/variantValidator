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
                sh 'chmod a+r ./*'
            }
        }
        stage("Where am I") {
            steps {
                sh 'echo $HOME'
                sh 'pwd'
                sh 'ls -l'
            }
        }
        stage("Build and Run VVTA PostgreSQL") {
            steps {
                script {
                    def vvtaContainer = docker.build("postgres-vvta-${CONTAINER_SUFFIX}", "./db_dockerfiles/vvta/Dockerfile")
                    vvtaContainer.run("-p 5432:5432 -d")
                    sh 'echo Building and running VVTA PostgreSQL'
                }
            }
        }
        stage("Build and Run Validator MySQL") {
            steps {
                script {
                    def validatorContainer = docker.build("mysql-validator-${CONTAINER_SUFFIX}", "./db_dockerfiles/vdb/Dockerfile")
                    validatorContainer.run("-p 3306:3306 -d")
                    sh 'echo Building and running Validator MySQL'
                }
            }
        }
        stage("Build and Run SeqRepo") {
            steps {
                script {
                    def seqRepoContainer = docker.build("sqlite-seqrepo-${CONTAINER_SUFFIX}", "./db_dockerfiles/vvsr/Dockerfile")
                    seqRepoContainer.run("-p 3306:3306 -p 5432:5432 -d")
                    sh 'echo Building and running SeqRepo'
                }
            }
        }
        stage("Build and Run VariantValidator") {
            steps {
                script {
                    def variantValidatorContainer = docker.build("variantvalidator-${CONTAINER_SUFFIX}", "./Dockerfile")
                    variantValidatorContainer.run("-p 3306:3306 -p 5432:5432 -d")
                    sh 'echo Building and running VariantValidator'
                }
            }
        }
        stage("Run Pytest and Codecov") {
            steps {
                // Run pytest and codecov in the variantvalidator container
                sh 'docker ps'
                // sh 'docker exec variantvalidator-${CONTAINER_SUFFIX} pytest --cov-report=term --cov=VariantValidator/'
                // sh 'docker exec variantvalidator-${CONTAINER_SUFFIX} codecov'
            }
        }
        stage("Cleanup Docker") {
            steps {
                sh 'docker stop postgres-vvta-${CONTAINER_SUFFIX}'
                sh 'docker rm postgres-vvta-${CONTAINER_SUFFIX}'
                sh 'docker stop mysql-validator-${CONTAINER_SUFFIX}'
                sh 'docker rm mysql-validator-${CONTAINER_SUFFIX}'
                sh 'docker rmi postgres-vvta-${CONTAINER_SUFFIX}'
                sh 'docker rmi mysql-validator-${CONTAINER_SUFFIX}'
                sh 'docker stop variantvalidator-${CONTAINER_SUFFIX}'
                sh 'docker rm variantvalidator-${CONTAINER_SUFFIX}'
                sh 'docker rmi variantvalidator-${CONTAINER_SUFFIX}'
            }
        }
    }
}
