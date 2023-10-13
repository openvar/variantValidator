pipeline {
    agent {
        docker {
            image "docker:24.0.6-git" // Set the Docker image for the Jenkins agent
        }
    }
    environment {
        CODECOV_TOKEN = "50dd5c2e-4259-4cfa-97a7-b4429e0d179e" // Define an environment variable for the Codecov token
        CONTAINER_SUFFIX = "${BUILD_NUMBER}" // Use the build number as a container suffix for uniqueness
        DOCKER_NETWORK = "variantvalidator_docker_network-$CONTAINER_SUFFIX" // Create a unique Docker network for this build
        DATA_VOLUME = "docker-shared-space" // Define a data volume for shared data
    }
    stages {
        stage("Clone Repository Remove dangling docker components and Create Docker Network") {
            steps {
                checkout scm // Checkout the source code from the configured source code management system
                sh 'docker system prune -f' // Remove unused Docker resources
                sh 'docker network create $DOCKER_NETWORK' // Create a Docker network for containers
            }
        }
        stage("Build and Run VVTA PostgreSQL") {
            steps {
                script {
                    def dockerfile = './db_dockerfiles/vvta/Dockerfile' // Define the Dockerfile path
                    def vvtaContainer = docker.build("postgres-vvta-${CONTAINER_SUFFIX}", "--no-cache -f ${dockerfile} ./db_dockerfiles/vvta")
                    // Build and run a PostgreSQL container for VVTA
                    vvtaContainer.run("-p 5432:5432 -d --name vv-vvta --network $DOCKER_NETWORK --shm-size=2g")
                    sh 'echo Building and running VVTA PostgreSQL' // Display a message
                }
            }
        }
        stage("Build and Run Validator MySQL") {
            steps {
                script {
                    def dockerfile = './db_dockerfiles/vdb/Dockerfile' // Define the Dockerfile path
                    def validatorContainer = docker.build("mysql-validator-${CONTAINER_SUFFIX}", "--no-cache -f ${dockerfile} ./db_dockerfiles/vdb")
                    // Build and run a MySQL container for the Validator
                    validatorContainer.run("-p 3306:3306 -d --name vv-vdb --network $DOCKER_NETWORK")
                    sh 'echo Building and running Validator MySQL' // Display a message
                }
            }
        }
        stage("Build and Run SeqRepo") {
            steps {
                script {
                    def dockerfile = './db_dockerfiles/vvsr/Dockerfile' // Define the Dockerfile path
                    def seqRepoContainer = docker.build("sqlite-seqrepo-${CONTAINER_SUFFIX}", "--no-cache -f ${dockerfile} ./db_dockerfiles/vvsr")
                    // Build and run a SQLite SeqRepo container
                    seqRepoContainer.run("--network $DOCKER_NETWORK --name vv-seqrepo -v $DATA_VOLUME:/usr/local/share:rw")
                    sh 'echo Building and running SeqRepo' // Display a message
                }
            }
        }
        stage("Build and Run VariantValidator") {
            steps {
                script {
                    def dockerfile = './Dockerfile' // Define the Dockerfile path
                    def variantValidatorContainer = docker.build("variantvalidator-${CONTAINER_SUFFIX}", "--no-cache -f ${dockerfile} .")
                    // Build and run the VariantValidator container
                    variantValidatorContainer.run("-v $DATA_VOLUME:/usr/local/share:rw -d --name variantvalidator --network $DOCKER_NETWORK")
                    sh 'echo Building and running VariantValidator' // Display a message
                }
            }
        }
        stage("Run Pytest and Codecov") {
            steps {
                script {
                    sh 'docker ps' // List running Docker containers
                    def connectionSuccessful = false

                    for (int attempt = 1; attempt <= 5; attempt++) {
                        echo "Attempt $attempt to connect to the database..."
                        def exitCode = sh(script: '''
                            docker exec -e PGPASSWORD=uta_admin variantvalidator psql -U uta_admin -d vvta -h vv-vvta -p 5432
                        ''', returnStatus: true)

                        if (exitCode == 0) {
                            connectionSuccessful = true
                            echo "Connected successfully! Running pytest..."

                            // Run pytest
                            sh 'docker exec variantvalidator pytest --cov-report=term --cov=VariantValidator/'

                            // Check for test failures in the captured output
                            if (currentBuild.rawBuild.getLog(2000).join('\n').contains("test summary info") && currentBuild.rawBuild.getLog(2000).join('\n').contains("FAILED")) {
                                error "Pytest completed with test failures"
                            }

                            // Check the Jenkins console log for pytest exit code
                            def pytestExitCode = currentBuild.rawBuild.getLog(2000).find { line -> line =~ /.*Pytest exit code: (\d+).*/ }
                            if (pytestExitCode) {
                                pytestExitCode = Integer.parseInt(pytestExitCode.replaceAll(/.*Pytest exit code: (\d+).*/, '$1'))
                                if (pytestExitCode != 0) {
                                    error "Pytest failed with exit code $pytestExitCode"
                                }
                            }

                            // Run Codecov with the provided token and branch name
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
    }
    post {
        always { // This ensures cleanup is executed regardless of build outcome
            script {
                // Cleanup Docker
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
        failure {
            script {
                currentBuild.result = 'FAILURE' // Mark the build as FAILURE
                echo 'Pipeline failed. Please check the logs for details.'
                // Update README badges on failure
                sh 'sed -i "s|\\[\\![codecov\\](.*)\\]|\\[![codecov](https://codecov.io/gh/openvar/variantValidator/branch/${BRANCH_NAME}/graph/badge.svg)](https://codecov.io/gh/openvar/variantValidator)|" README.md'
                sh 'sed -i "s|\\[\\!\\[Build Status\\](.*)\\]|\\[![Build Status](https://d174-130-88-226-17.ngrok-free.app/buildStatus/icon?job=VariantValidator+CI%2Fci&branch=${BRANCH_NAME})](https://d174-130-88-226-17.ngrok-free.app/job/VariantValidator%20CI%2Fci/job/ci/)|" README.md'
                // Commit and push to GitHub
                sh 'git branch'
                echo "${BRANCH_NAME}"
                sh 'git commit -am "Update README badges to failure by Jenkins"'
                sh 'git push -u origin ${BRANCH_NAME}'
            }
        }
        success {
            script {
                currentBuild.result = 'SUCCESS' // Mark the build as SUCCESS
                echo 'Pipeline succeeded! Your project is built and tested.'
                // Update README badges on success
                sh 'sed -i "s|\\[![codecov](.*)\\]|[![codecov](https://codecov.io/gh/openvar/variantValidator/branch/${BRANCH_NAME}/graph/badge.svg)](https://codecov.io/gh/openvar/variantValidator)|" README.md'
                sh 'sed -i "s|\\[![Build Status](.*\\)|[![Build Status](https://d174-130-88-226-17.ngrok-free.app/buildStatus/icon?job=VariantValidator+CI%2Fci&branch=${BRANCH_NAME})](https://d174-130-88-226-17.ngrok-free.app/job/VariantValidator%20CI%2Fci/job/ci/)|" README.md'
                // Commit and push to GitHub
                sh 'git branch'
                echo "${BRANCH_NAME}"
                sh 'git commit -um "Update README badges to success by Jenkins"'
                sh 'git push -u origin ${BRANCH_NAME}'
            }
        }
    }
}
