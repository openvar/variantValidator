ABOUT VARIANTVALIDATOR

VariantValidator is a user-friendly software tool designed to validate the syntax and 
parameters of DNA variant descriptions according to the HGVS Sequence Variant 
Nomenclature. 

VariantValidator ensures that users are guided through the intricacies of the HGVS 
nomenclature, e.g. if the user makes a mistake, VariantValidator automatically corrects 
the mistake if it can, or provides helpful guidance if it cannot. In addition, 
VariantValidator accurately interconverts between transcript variant descriptions and 
genomic variant descriptions in HGVS and Variant Call Format (VCF)

VariantValidator interfaces with the hgvs package to parse, format, and manipulate 
biological sequence variants.  See https://github.com/biocommons/hgvs/ for details of the
hgvs package

VariantValidator is a highly functional platform enabling high-throughput and embeddable
utilisation of functionality of https://variantvalidator.org/


FEATURES

The basic functionality of https://variantvalidator.org/ and VarinantValidator is documented here https://www.ncbi.nlm.nih.gov/pubmed/28967166


OTHER FEATURES

VariantValidator simultaneously and accurately projects genomic sequence variations onto all overlapping transcript reference sequences, and vice-versa

Alternatively, genomic sequence variation can be projected onto a specified single, or specified subset of transcript reference sequences for any given gene

Projection of sequence variations between reference sequences takes account of discrepancies between genomic and transcript reference sequences, thus ensuring an accurate prediction of the effect on encoded proteins for every gene

For sequence variations falling within the open reading frames of genes, VariantValidator automatically projects sequence variants via the transcript reference sequence onto genome builds GRCh38, GRCh37, hg38 and hg19 (HGVS format and VCF components), including projection onto relevant Alternative genomic reference sequences, the composition of which varies between patched GRC genome builds and static hg genome builds

REQUIREMENTS
MySQL
validator MySQL database
User account for the database. 
Default:
user = vvadmin  
password = var1ant

OPTIMAL PERFORMANCE REQUIREMENTS
PostgreSQL version 9.5 or higher (version 10 not supported)
sqlite3 >= 3.8.0

For optimal performance, we recommend local installations of the Universal Transcript Archive (UTA) and SeqRepo. Database information http://dl.biocommons.org/ 

Install UTA 
Instructions at https://github.com/biocommons/uta/
Restore the database http://dl.biocommons.org/uta
gzip -cdq <CURRENT UTA SCHEMA>.pgd.gz | psql -U uta_admin -v ON_ERROR_STOP=0 -d uta –Eae

Install SeqRepo
pip install biocommons.seqrepo
Pull the current SeqRepo database http://dl.biocommons.org/seqrepo
seqrepo --root-directory <PATH TO>/seqrepo pull -i <CURRENT SEQREPO DATABASE>

NOTE FOR MULTITHREADED APPLICATIONS

Connections to SQLite3 need to be modified within SeqRepo. You must have a thread-safe installation of SQLite3. 

Example
sqlite3.connect(self._db_path, check_same_thread=False)

INSTALLATION
pip install --upgrade setuptools, pip
pip install variantValidator_0.1.0.tar.gz

Variant Validator also requires mysql_connector >= 2.1.4. mysql_connector is  OS specific and can be downloaded from https://dev.mysql.com/downloads/connector/python/

REQUIRED ALTERATIONS TO DEPENDENCIES
The following alteration must be made to the hgvs Python package https://github.com/biocommons/hgvs/ 

CONFIGURATION
A config.ini file is located 
<PATH TO>/variantValidator/configuration 
within which URLs for the required databases can be set

Example
[mysql]
host = 127.0.0.1
database = validator
user = vvadmin  
password = var1ant

[EntrezID]
entrezid = <YOUR VariantValidator ACCOUNT EMAIL ADDRESS>

[SeqRepo]
seqrepo_dir = <PATH TO>/seqrepo/<CURRENT SEQREPO DATABASE>

[UTA]
uta_url = postgresql://uta_admin:uta_admin@127.0.0.1/uta/<CURRENT UTA SCHEMA>

Alternatively, these configurations can be set on via environment variables

Example
import os
os.environ['HGVS_SEQREPO_DIR'] = <PATH TO>/seqrepo/<CURRENT SEQREPO DATABASE>
os.environ['UTA_DB_URL'] = 'postgresql://uta_admin:uta_admin@127.0.0.1/uta/<CURRENT UTA SCHEMA>'
VALIDATOR_DB_URL = 'mysqlx://vvadmin:var1ant@127.0.0.1/validator'
os.environ['VALIDATOR_DB_URL'] = VALIDATOR_DB_URL
os.environ ['ENTREZ_ID'] = '<YOUR VariantValidator ACCOUNT EMAIL ADDRESS>'

# optional
os.environ['VALIDATOR_DEBUG'] = 'TRUE'

To file <PATH TO>/hgvs/_data/defaults.ini
Edit the following line as below
[formatting]
max_ref_length = 1000000

TOP LEVEL FUNCTIONS

# Import
import variantValidator
from variantValidator import variantValidator

# Check database configurations example
variantValidator.my_config()

# Validate a sequence variant example
# Accepted formats
# NM_000088.3:c.589G>T
# NC_000017.10:g.48275363C>A
# NG_007400.1:g.8638G>T
# LRG_1:g.8638G>T
# LRG_1t1:c.589G>T
# 17-50198002-C-A (GRCh38)
# chr17:50198002C>A (GRCh38)

variant = ' NM_000088.3:c.589G>T'
selected_assembly = 'GRCh37' # or GRCh37, hg19, hg38
select_transcripts = 'all'
validation = variantValidator.validator(variant, selected_assembly, select_transcripts)

# Specify transcripts example
variant = '1-43212925-C-T'
selected_assembly = 'GRCh37' 

# no transcripts specified
select_transcripts = 'all' 

# Single transcript specified
select_transcripts = 'NM_022356.3' 

# 3 transcripts specified
select_transcripts = 'NM_022356.3| NM_001146289.1| NM_001243246.1' 

validation = variantValidator.validator(variant, selected_assembly, select_transcripts)

# View supported transcripts for a gene example
# HGNC gene symbol https://www.genenames.org/
variantValidator.validator.gene2transcripts ('HTT')
# RefSeq Transcript
variantValidator.validator.gene2transcripts (' NM_002111.8')

# Get reference sequence for HGVS variant description
variantValidator.validator.hgvs2ref('NM_000088.3:c.589_594del')

LISCENSE
See LISCENSE.txt

CITE US
Hum Mutat. 2017 Oct 1. doi: 10.1002/humu.23348

VariantValidator: Accurate validation, mapping and formatting of sequence variation descriptions.

Freeman PJ, Hart RK, Gretton LJ, Brookes AJ, Dalgleish R.

# Copyright (C) 2018  Peter Causey-Freeman, University of Leicester
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>