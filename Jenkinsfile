pipeline {
    agent {
        docker {
            image "docker:24.0.6-git" // Set the Docker image for the Jenkins agent
        }
    }
    environment {
        CODECOV_TOKEN = credentials('CODECOV_TOKEN') // Use the Codecov token from Jenkins secret
        CONTAINER_SUFFIX = "${BUILD_NUMBER}" // Use the build number as a container suffix for uniqueness
        DOCKER_NETWORK = "variantvalidator_docker_network-$CONTAINER_SUFFIX" // Create a unique Docker network for this build
        DATA_VOLUME = "docker-shared-space" // Define a data volume for shared data
    }
    stages {
        stage("Clone Repository Remove dangling docker components and Create Docker Network") {
            steps {
                checkout scm // Checkout the source code from the configured source code management system
                sh 'docker system prune --all --volumes --force' // Remove unused Docker resources
                sh 'docker network create $DOCKER_NETWORK' // Create a Docker network for containers
            }
        }
        stage("Switch to Git Branch") {
            steps {
                sh "git checkout ${BRANCH_NAME}"
                sh "git pull"
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

                    // Run variantValidatorContainer and Mount the DATA_VOLUME
                    variantValidatorContainer.run("-v $DATA_VOLUME:/usr/local/share:rw -d --name variantvalidator --network $DOCKER_NETWORK")

                    // Display a message indicating that VariantValidator is being built and run
                    sh 'echo Building and running VariantValidator'
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

                            // Run pytest && Run Codecov with the provided token and branch name
                            sh 'docker exec variantvalidator pytest -n 3 --cov=VariantValidator --cov=VariantFormatter --cov-report=term tests/'

                            // Send coverage report to Codecov
                            sh 'docker exec variantvalidator codecov -t $CODECOV_TOKEN -b ${BRANCH_NAME}'

                            // Check for test failures in the captured output
                            if (currentBuild.rawBuild.getLog(2000).join('\n').contains("test summary info") && currentBuild.rawBuild.getLog(2000).join('\n').contains("FAILED")) {
                                failure(message:"Pytest completed with test failures")
                            }

                            // Check the Jenkins console log for pytest exit code
                            def pytestExitCode = currentBuild.rawBuild.getLog(2000).find { line -> line =~ /.*Pytest exit code: (\d+).*/ }
                            if (pytestExitCode) {
                                pytestExitCode = Integer.parseInt(pytestExitCode.replaceAll(/.*Pytest exit code: (\d+).*/, '$1'))
                                if (pytestExitCode != 0) {
                                    failure(message:"Pytest failed with exit code $pytestExitCode")
                                }
                            }
                            break
                        }

                        echo "Connection failed. Waiting for 60 seconds before the next attempt..."
                        sleep 60
                    }

                    if (!connectionSuccessful) {
                        failure(message:"All connection attempts failed. Exiting...")
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
                sh 'docker system prune --all --volumes --force'
            }
        }
    }
}
