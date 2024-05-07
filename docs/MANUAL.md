# Manual

## Configuration

After first installing VariantValidator, a configuration file will need to be created and edited to contain the database credentials and locations. 
By default the edited configuration will be placed in the users home directory (`~/.variantvalidator`), this location can be changed for all users by editing the `VariantValidator/settings.py` file.
To create this file automatically, run the configuration script installed alongside the package.

```bash
$ python bin/vv_configure.py
```

This will ask you to enter a value for each item in the configuration file. 
The default/existing value is shown in square brackets and will continue to be used 
if you don't enter anything else. The items in the configuration file are:

```text
[mysql]
host = 127.0.0.1
port = 3306
database = validator
port = 3306
user = USERNAME
password = PASSWORD
version = VERSION

[seqrepo]
version = VERSION/master
location = /PATH/TO/SEQREPO
require_threading = False  # Set to True if you want to deploy VV in a web application

[postgres]
host = 127.0.0.1
database = vvta
port = 5432
version = VERSION
user = USERNAME
password = PASSWORD

[logging]
# Levels control verbosity and can be set to "CRITICAL" "ERROR" "WARNING" "INFO" or "DEBUG".
log = True
console = INFO
file = WARNING

[Entrez]
email = YOUR@EMAIL.COM
api_key = YOUR_API_KEY
```

The values in capitals must be replaced for VariantValidator to run, see below for more details.

**Standard formats for configuration**

*Update the dates of each database build to match the version you downloaded*

- mysql database version e.g. vvdb_2022_04 (update the data to the latest build to match your download)
- vvSeqRepo [VERSION/master] e.g. VV_SR_2022_02/master
- VVTA database version e.g. vvta_2022_02

**This script can also be used to uodate your configuration at a later date**

#### Logging

By default VariantValidator will log to both the console and to a file, the output level for each can be set in the configuration file.
The levels control verbosity and can be set to "CRITICAL", "ERROR", "WARNING", "INFO" or "DEBUG". To turn off logging, set the log configuration to "False". The log file name and
log options can be changed for all users by editing the `VariantValidator/settings.py` file. By default the file log is 
set to output in the users home directory (`~/.vv_errorlog`).

#### Entrez

For access to the NCBI Entrez database  you must provide a valid email address in 
the respective configuration setting. Optionally, you can also provide an NCBI API key that will increase the number of requests
made per second. See [this article](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) on how to generate an API key.


## Operation

To run VariantValidator, we have provided the installed script `bin/variant_validator.py`, running this with the flag `-h` shows the running options:

```text
usage: variant_validator.py [-h] -v VARIANT [VARIANT ...]
                            [-g [{GRCh37,GRCh38,hg19,hg38}]]
                            [-t [TRANSCRIPTS]] [-s {individual,batch}]
                            [-f {dict,table,json}] [-o OUTPUT] [-m]

optional arguments:
  -h, --help            show this help message and exit
  -v VARIANT [VARIANT ...], --variant VARIANT [VARIANT ...]
                        Variant(s) to validate
  -g [{GRCh37,GRCh38,hg19,hg38}], --genome [{GRCh37,GRCh38,hg19,hg38}]
                        Genome assembly (default: GRCh37)
  -t [TRANSCRIPTS], --transcripts [TRANSCRIPTS]
                        Transcripts to output results for (default: all)
  -s {individual,batch}, --submission {individual,batch}
                        Submit variants individually or as a single batch
                        validation (default: individual)
  -f {dict,table,json}, --output_format {dict,table,json}
                        Output validations as a list or as a dictionary
                        (default: dict)
  -o OUTPUT, --output OUTPUT
                        Specifies the output file (default: stdout)
  -m, --meta            Also output metadata (default: False)
```

From this script you can run the validator with a number of different input and output options.

**Simple Example**
```bash
# Note: Variant descriptions must always be contained in '' or ""
$ bin/variant_validator.py -v 'NC_000017.11:g.50198002C>A' -g GRCh38 -t mane -s individual -f json -m -o stdout
```

**Capture the output using Python**
```python
import subprocess
validation = subprocess.run("./bin/variant_validator.py -v 'NC_000017.11:g.50198002C>A' -g GRCh38 -t mane -s individual -f json -m -o stdout", stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
print(validation.stdout.decode("utf-8"))
```

You can also import and use the package directly within python. For example:

```python
import VariantValidator
validator = VariantValidator.Validator()

# To validate a variant
output = validator.validate('NM_000088.3:c.589G>T', 'GRCh37', 'all')
# This returns an ValOutput object that can be used to output the results in a number of different ways (dictionary, json or table)
output.format_as_dict(with_meta=True)

# The Validator object also contains other useful methods, such as finding all transcripts from a gene ID/symbol
validator.gene2transcripts('COL1A1')
```

The accepted format for variants include:
```text
NM_000088.3:c.589G>T
NC_000017.10:g.48275363C>A
NG_007400.1:g.8638G>T
LRG_1:g.8638G>T
LRG_1t1:c.589G>T
17-50198002-C-A  # Note this variant is in the context of GRCh38
chr17:50198002C>A  # Note this variant is in the context of GRCh38
```
