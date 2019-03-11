# Variant Validator Operation Manual

## Configuration

Presently VariantValidator uses a combination of environment variables to configure itself. The configuration file is in /VariantValidator/configuration/config.ini and should be edited with the current user's details. Specifically, the section:
```
[mysql]
host = 127.0.0.1
database = validator
user = vvadmin  
password = var1ant
```
needs to be changed if the VariantValidator database login details are different.

The section
```
[logging]
string = error file console trace
```
which can be changed to alter the level of verbosity of the validator output. Alternatively you can set the environment variable VALIDATOR_DEBUG to a string of the same format.
The string should contain any of the following words:
* file - Writes the logging output to the "vvLog.txt" file. Without the word "file" in the environment variable, the logs will be posted instead to the console.
* debug - Logs all events, including debugging.
* trace - Used for diagnosis during development.
* info - Information events on the decisions the validator is making are logged.
* warning - Warnings indicate malformed variants. This is the default logging level.
* error - Variants that produce errors are nonsensical to the point where they cannot be validated.
* critical - Fatal errors that crash the validator are logged at this level.
During a test, this is set to maximum verbosity.

The validator itself will set environment variables to allow for the correct operation of HGVS software.

## Operation

Python scripts importing VariantValidator will have to set up a last few configuration variables before they can proceed. These variables must be set in such a way that they don't go out of scope - otherwise the validator won't work.

This example script will validate the variant NM_000088.3:c.589G>T and then print the output as a json file. You might need to change it to point to the correct seqrepo directory.

```
import json
import os
seqrepo_current_version = '2018-08-21'
HGVS_SEQREPO_DIR = '~/seqrepo/' + seqrepo_current_version
os.environ['HGVS_SEQREPO_DIR'] = HGVS_SEQREPO_DIR
uta_current_version = 'uta_20180821'
UTA_DB_URL = 'postgresql://uta_admin:uta_admin@127.0.0.1/uta/' + uta_current_version
os.environ['UTA_DB_URL'] = UTA_DB_URL
from VariantValidator import variantValidator
variantValidator.my_config()
```
From this point onward, 
```
variant = 'NM_000088.3:c.589G>T'
select_transcripts = 'all'
selected_assembly = 'GRCh37'
validation = variantValidator.validator(variant, selected_assembly, select_transcripts)
print json.dumps(validation, sort_keys=True, indent=4, separators=(',', ': '))
```
Much of the script is currently related to setting up environment variables. In future versions, this information will be stored in a local configuration file.

The accepted formats for variants include:
```
NM_000088.3:c.589G>T
NC_000017.10:g.48275363C>A
NG_007400.1:g.8638G>T
LRG_1:g.8638G>T
LRG_1t1:c.589G>T
17-50198002-C-A (GRCh38)
chr17:50198002C>A (GRCh38)
```
Possible assemblies are:
```
GRCh37
GRCh38
hg19
hg38
```
You can select all transcripts by passing 'all', or use multiple transcripts with: `select_transcripts = 'NM_022356.3| NM_001146289.1| NM_001243246.1' `

VariantValidator produces a dictionary output that contain all possible interpretations of the input variant.

View supported transcripts for a gene example: HGNC gene symbol https://www.genenames.org/
```
variantValidator.validator.gene2transcripts('HTT') 
# RefSeq Transcript
variantValidator.validator.gene2transcripts(' NM_002111.8') 

# Get reference sequence for HGVS variant description
variantValidator.validator.hgvs2ref('NM_000088.3:c.589_594del')
```

## Unit testing

VariantValidator is written to be pytest-compatible. Run
`pytest`
in the variant validator root folder, the same as that in which this file resides. The test will take several minutes to complete, but runs through over three hundred common and malformed variants.

Note that you will need to set the environment variables first. 

```
export UTA_DB_URL="postgresql://uta_admin:uta_admin@127.0.0.1/uta/uta_20180821"
export HGVS_SEQREPO_DIR="path/to/seqreo"
pytest
```

