# Variant Validator Operation Manual

## Configuration
Variant Validator will create a configuration file for each user if it does not detect one, located in the folder
 > ~/.config/VariantValidator/config.ini
This file, freshly created, will be missing the path to the SeqRepo directory which you should fill in after installation accordingly. If the configuration file hasn't been filled in correctly, the validator will exit immediately with an error.

It's possible to use a remote seqrepo directory, at a cost of greatly reduced performane.

The mysql database is configured in this section:
 > [mysql]
 > host = 127.0.0.1
 > database = validator
 > user = vvadmin  
 > password = var1ant
Information here also needs to be changed if the variant validator database login details are different.

The [uta] section also contains path information to the UTA archive, if it's installed.

The section
 > [logging]
contains several headings which can be changed to alter the level of verbosity of the validator output.

* file - If "True", writes the logging output to the "vvLog.txt" file in the current working directory. While useful for diagnostics, logging in this way has permissions issues and will fill up the hard disk of an automated installation quickly.
* level - Can be one of several values. All errors below the selected level of severity will not be logged. By default, and to help with setting things up, the info level statements will be logged, but you should change this to make the validator less talkative in normal use.
** debug - Logs all events, including debugging.
** info - Information events on the decisions the validator is making are logged.
** warning - Warnings indicate malformed variants. This is the default logging level.
** error - Variants that produce errors are nonsensical to the point where they cannot be validated.
** critical - Fatal errors that crash the validator are logged at this level.
* trace - Used for diagnosis during development. Can be set to 'True' if you need to profile the validator code.

The validator itself will set environment variables to allow for the correct operation of HGVS software.

## Operation

Validating variants, provided the software is installed correctly, is as simple as:

> from VariantValidator import Validator
> 
> validator = Validator()
> variant = 'NC_000012.11:g.122064776delG'
> select_transcripts = 'all'
> selected_assembly = 'GRCh37'
> 
> out=Validator().validate(variant, selected_assembly, select_transcripts)

The 'out' object is a simple dictionary containing the genetic information of the validated variant. The simpleTestScript.py will validate this variant and then print the output nicely as a json.

The accepted formats for variants include:
> NM_000088.3:c.589G>T
> NC_000017.10:g.48275363C>A
> NG_007400.1:g.8638G>T
> LRG_1:g.8638G>T
> LRG_1t1:c.589G>T
> 17-50198002-C-A (GRCh38)
> chr17:50198002C>A (GRCh38)

Possible assemblies are:
> GRCh37
> hg19
> hg38

You can select all transcripts by passing 'all', or use multiple transcripts with:
> select_transcripts = 'NM_022356.3| NM_001146289.1| NM_001243246.1' 

## Unit testing

Variant Validator is written to be pytest-compatible. Run
> pytest
in the variant validator testing folder, the same as that in which this file resides. The test will take several minutes to complete, but runs through over three hundred common and malformed variants.

 	
