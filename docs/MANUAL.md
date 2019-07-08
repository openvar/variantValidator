# Manual

## Configuration

After first installing Variant Validator, a configuration file will need to be created and edited to contain the database credentials and locations. To do this run the configuration script installed alongside the package.

```bash
vv_configure.py
```

This will ask you to enter a value for each item in the configuration file. 
The default/existing value is shown in square brackets and will continue to be used 
if you don't enter anything else. The items in the configuration file are:

```text
[mysql]
host = localhost
database = validator
user = USERNAME
password = PASSWORD

[seqrepo]
version = 2018-08-21
location = /PATH/TO/SEQREPO

[postgres]
host = localhost
database = uta
version = uta_20180821
user = USERNAME
password = PASSWORD

[logging]
level = info
console = true
file = false
trace = false

[EntrezID]
entrezid = admin@variantvalidator.org

[liftover]
location = /path/to/liftover
```

The values in capitals must be replaced for Variant Validator to run.

By default the edited configuration will be placed in the users home directory (`~/.variantvalidator`), this location can be changed for all users by editing the `VariantValidator/settings.py` file.

#####Liftover

If the UCSC Liftover [files](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/) have been previously downloaded their location can be set within the configuration file. By default the necessary files will be downloaded automatically when first requested. 

## Database updates

To import the initial data into the Validator MySQL database, run the following script:

```bash
update_vdb.py
```

This will download the required data to convert between LRG and RefSeq IDs. We recommend re-running this command on a regular basis as changes are continually made to the RefSeq and LRG collections.  

## Operation

To run Variant Validator, we have provided the installed script `variant_validator.py`, running this with the flag `-h` shows the running options:

```text
usage: variant_validator.py [-h] -v VARIANT [VARIANT ...]
                            [-g [{GRCh37,GRCh38,hg19,hg38}]]
                            [-t [TRANSCRIPTS]] [-s {individual,batch}]
                            [-f {dict,list,json}] [-o OUTPUT]

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
  -f {dict,list,json}, --output_format {dict,list,json}
                        Output validations as a list or as a dictionary
                        (default: dict)
  -o OUTPUT, --output OUTPUT
                        Specifies the output file (default: stdout)
```

From this script you can run the validator with a number of different input and output options.

You can also import and use the package directly within python. For example:

```python
import VariantValidator
validator = VariantValidator.Validator()

# To validate a variant
output = validator.validate('NM_000088.3:c.589G>T', 'GRCh37', 'all')
# This returns an ValOutput object that can be used to output the results in a number of different ways
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
17-50198002-C-A (GRCh38)
chr17:50198002C>A (GRCh38)
```

