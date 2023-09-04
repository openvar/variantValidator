pipeline {
    agent any

    environment {
        PGPORT = '5433'
        MYSQL_PORT = '3306'
        CODECOV_TOKEN = '50dd5c2e-4259-4cfa-97a7-b4429e0d179e'
    }

    stages {
        stage("Before Install") {
            steps {
                script {
                    sh 'sudo mount -o remount,size=50% /var/ramfs'
                    sh 'mysql -e "CREATE DATABASE validator;"'
                    sh 'df -h'
                    sh 'sudo apt-get -y install tabix'
                    sh 'sudo sed -i -e "/local.*peer/s/postgres/all/" -e "s/peer\\|md5/trust/g" /etc/postgresql/12/main/pg_hba.conf'
                    sh 'sudo service postgresql@12-main restart'
                    sleep(3)
                    sh 'createuser -e --createdb uta_admin'
                    sh 'createdb -e vvta -O uta_admin'
                    sh 'psql -d vvta -U postgres -c "CREATE USER ta_user WITH PASSWORD \'read_only\'"'
                    sh 'wget --output-document=vvta_2023_05_noseq.sql.gz https://www528.lamp.le.ac.uk/vvdata/vvta/vvta_2023_05_noseq.sql.gz'
                    sh 'gunzip -c vvta_2023_05_noseq.psql.gz | psql --quiet vvta'
                    sh 'psql -d vvta -U postgres -c "GRANT SELECT ON vvta_2023_05.gene TO public;"'
                    sh 'psql -d vvta -U postgres -c "GRANT SELECT ON ALL TABLES IN SCHEMA public TO ta_user;"'
                    sh 'psql -d vvta -U postgres -c "GRANT SELECT ON vvta_2023_05.tx_def_summary_v TO ta_user;"'
                    sh 'psql -d vvta -U postgres -c "GRANT SELECT ON vvta_2023_05.tx_exon_aln_v TO ta_user;"'
                    sh 'psql -d vvta -U postgres -c "GRANT SELECT ON vvta_2023_05.transcript_lengths_v TO ta_user;"'
                    sh 'psql -d vvta -U postgres -c "GRANT SELECT ON vvta_2023_05.exon_set TO ta_user;"'
                    sh 'cp configuration/continuous_integration.ini "$HOME"/.variantvalidator'
                    sh 'wget --output-document=validator_2023_08.sql.gz https://www528.lamp.le.ac.uk/vvdata/validator/validator_2023_08.sql.gz'
                    sh 'gunzip validator_2023_08.sql.gz'
                }
            }
        }

        stage("Install") {
            steps {
                script {
                    sh 'pip install .'
                    sh 'mkdir "$HOME"/vvta_seqrepo'
                    sh 'wget --output-document="$HOME"/vvta_seqrepo/VV_SR_2023_05.tar https://www528.lamp.le.ac.uk/vvdata/vv_seqrepo/VV_SR_2023_05.tar'
                    sh 'cd "$HOME"/vvta_seqrepo/'
                    sh 'tar -xvf VV_SR_2023_05.tar'
                    sh 'cd -'
                }
            }
        }

        stage("Test") {
            steps {
                script {
                    sh 'pytest --cov-report=term --cov=VariantValidator/'
                }
            }
        }

        stage("After Script") {
            steps {
                script {
                    sh 'codecov'
                }
            }
        }
    }
}
