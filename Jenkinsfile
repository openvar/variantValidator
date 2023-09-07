pipeline {
    agent any

    stages {
        stage("Clone Repository") {
            steps {
                checkout scm
            }
        }
        stage("Build and Run VVTA PostgreSQL") {
            steps {
                script {
                    def vvtaContainer = docker.build("postgres-vvta-${CONTAINER_SUFFIX}", "-f ./db_dockerfiles/vvta/Dockerfile ./db_dockerfiles/vvta")
                    vvtaContainer.inside("-p 5432:5432 -d") {
                        // Steps to run inside the VVTA PostgreSQL container
                        sh 'echo Building and running VVTA PostgreSQL'
                        sh 'docker ps' // Check running containers
                    }
                }
            }
        }
        stage("Build and Run Validator MySQL") {
            steps {
                script {
                    def validatorContainer = docker.build("mysql-validator-${CONTAINER_SUFFIX}", "-f ./db_dockerfiles/vdb/Dockerfile ./db_dockerfiles/vdb")
                    validatorContainer.inside("-p 3306:3306 -d") {
                        // Steps to run inside the Validator MySQL container
                        sh 'echo Building and running Validator MySQL'
                        sh 'docker ps' // Check running containers
                    }
                }
            }
        }
        stage("Build and Run SeqRepo") {
            steps {
                script {
                    def seqRepoContainer = docker.build("sqlite-seqrepo-${CONTAINER_SUFFIX}", "-f ./db_dockerfiles/vvsr/Dockerfile ./db_dockerfiles/vvsr")
                    seqRepoContainer.inside("-p 3306:3306 -p 5432:5432 -d") {
                        // Steps to run inside the SeqRepo container
                        sh 'echo Building and running SeqRepo'
                        sh 'docker ps' // Check running containers
                    }
                }
            }
        }
        stage("Build and Run VariantValidator") {
            steps {
                script {
                    def variantValidatorContainer = docker.build("variantvalidator-${CONTAINER_SUFFIX}", "-f ./Dockerfile .")
                    variantValidatorContainer.inside("-p 3306:3306 -p 5432:5432 -d") {
                        // Steps to run inside the VariantValidator container
                        sh 'echo Building and running VariantValidator'
                        sh 'docker ps' // Check running containers
                    }
                }
            }
        }
        stage("Run Pytest and Codecov") {
            steps {
                // Run pytest and codecov in the VariantValidator container or other necessary steps
                sh 'docker ps' // Check running containers
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
                sh 'docker stop sqlite-seqrepo-${CONTAINER_SUFFIX}'
                sh 'docker rm sqlite-seqrepo-${CONTAINER_SUFFIX}'
                sh 'docker stop variantvalidator-${CONTAINER_SUFFIX}'
                sh 'docker rm variantvalidator-${CONTAINER_SUFFIX}'
            }
        }
    }
}
