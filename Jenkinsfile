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
                node {
                    def postgresImage = docker.build("postgres-vvta-${CONTAINER_SUFFIX}", "./vvta_docker.df")
                }
            }
        }
        stage("Build Validator MySQL") {
            steps {
                node {
                    def mysqlImage = docker.build("mysql-validator-${CONTAINER_SUFFIX}", "./vdb_docker.df")
                }
            }
        }
        stage("Build SeqRepo") {
            steps {
                node {
                    def seqrepoImage = docker.build("sqlite-seqrepo-${CONTAINER_SUFFIX}", "./vvsr_docker.df")
                }
            }
        }
        stage("Build VariantValidator") {
            steps {
                node {
                    def variantvalidatorImage = docker.build("variantvalidator-${CONTAINER_SUFFIX}", "./Dockerfile")
                }
            }
        }
        stage("Run Pytest and Codecov") {
            steps {
                // Add your run and cleanup steps here
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
