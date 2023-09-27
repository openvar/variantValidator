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
        stage("Clone Repository Remove dangling docker components and Create Docker Network") {
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
                    variantValidatorContainer.run("-v $DATA_VOLUME:/usr/local/share:rw -d --name variantvalidator --network $DOCKER_NETWORK")
                    sh 'echo Building and running VariantValidator'
                }
            }
        }
        stage("Run Pytest and Codecov") {
            steps {
                script {
                    sh 'docker ps'

                    def connectionSuccessful = false

                    for (int attempt = 1; attempt <= 5; attempt++) {
                        echo "Attempt $attempt to connect to the database..."
                        def exitCode = sh(script: '''
                            docker exec -e PGPASSWORD=uta_admin variantvalidator psql -U uta_admin -d vvta -h vv-vvta -p 5432
                        ''', returnStatus: true)

                        if (exitCode == 0) {
                            connectionSuccessful = true
                            echo "Connected successfully! Running pytest..."

                            // Run pytest and capture the output
                            def pytestOutput = sh(script: '''
                                docker exec variantvalidator pytest --cov-report=term --cov=VariantValidator/
                            ''', returnStdout: true)

                            // Check for test failures in the pytest output
                            if (pytestOutput.contains("collected") && pytestOutput.contains("failed")) {
                                error "Pytest completed with test failures:\n$pytestOutput"
                            }

                            // Check the exit code
                            def pytestExitCode = sh(script: 'echo \$?', returnStatus: true)

                            if (pytestExitCode != 0) {
                                error "Pytest failed with exit code $pytestExitCode"
                            }

                            docker exec variantvalidator codecov -t $CODECOV_TOKEN -b ${BRANCH_NAME}
                            break
                        }

                        echo "Connection failed. Waiting for 60 seconds before the next attempt..."
                        sleep 60
                    }

                    if (!connectionSuccessful) {
                        error "All connection attempts failed. Exiting..."
                    }
                }
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
                sh 'docker stop variantvalidator'
                sh 'docker rm variantvalidator'
                sh 'docker network rm $DOCKER_NETWORK'
            }
        }
    }
    post {
        failure {
            script {
                currentBuild.result = 'FAILURE' // Mark the build as FAILURE
                echo 'Pipeline failed. Please check the logs for details.'
                def errorMessage = currentBuild.rawBuild.getLog(1000).join('\n')
                echo "Error Message:\n${errorMessage}"

                // Update README badges on failure
                sh 'sed -i "s|\\[![codecov](.*\\)|[![codecov](https://codecov.io/gh/openvar/variantValidator/branch/${BRANCH_NAME}/graph/badge.svg)](https://codecov.io/gh/openvar/variantValidator)|" README.md'
                sh 'sed -i "s|\\[![Build Status](.*\\)|[![Build Status](https://d174-130-88-226-17.ngrok-free.app/buildStatus/icon?job=VariantValidator+CI%2Fci&branch=${BRANCH_NAME})](https://d174-130-88-226-17.ngrok-free.app/job/VariantValidator%20CI/job/ci/)|" README.md'

                // Commit and push to GitHub
                sh 'git commit -am "Update README badges to failure by Jenkins"'
                sh 'git push origin ${BRANCH_NAME}'
            }
        }
        success {
            script {
                currentBuild.result = 'SUCCESS' // Mark the build as SUCCESS
                echo 'Pipeline succeeded! Your project is built and tested.'

                // Update README badges on success
                sh 'sed -i "s|\\[![codecov](.*\\)|[![codecov](https://codecov.io/gh/openvar/variantValidator/branch/${BRANCH_NAME}/graph/badge.svg)](https://codecov.io/gh/openvar/variantValidator)|" README.md'
                sh 'sed -i "s|\\[![Build Status](.*\\)|[![Build Status](https://d174-130-88-226-17.ngrok-free.app/buildStatus/icon?job=VariantValidator+CI%2Fci&branch=${BRANCH_NAME})](https://d174-130-88-226-17.ngrok-free.app/job/VariantValidator%20CI%2Fci/job/ci/)|" README.md'

                // Commit and push to GitHub
                sh 'git commit -am "Update README badges to success by Jenkins"'
                sh 'git push origin ${BRANCH_NAME}'
            }
        }
    }
}
