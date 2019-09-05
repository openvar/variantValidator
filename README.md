# VariantValidator
[![codecov](https://codecov.io/gh/openvar/variantValidator/branch/restructuring_py3/graph/badge.svg)](https://codecov.io/gh/openvar/variantValidator) [![Build Status](https://travis-ci.org/openvar/variantValidator.png)](https://travis-ci.org/openvar/variantValidator)

## About

VariantValidator is a user-friendly software tool designed to validate the syntax and 
parameters of DNA variant descriptions according to the HGVS Sequence Variant 
Nomenclature. 

VariantValidator ensures that users are guided through the intricacies of the HGVS 
nomenclature, e.g. if the user makes a mistake, VariantValidator automatically corrects 
the mistake if it can, or provides helpful guidance if it cannot. In addition, 
VariantValidator accurately inter-converts between transcript variant descriptions and 
genomic variant descriptions in HGVS and Variant Call Format (VCF)

VariantValidator interfaces with the hgvs package to parse, format, and manipulate 
biological sequence variants.  See https://github.com/biocommons/hgvs/ for details of the
hgvs package

VariantValidator is a highly functional platform enabling high-throughput and embeddable
utilisation of functionality of https://variantvalidator.org/

## Features

The basic functionality of https://variantvalidator.org/ and VarinantValidator is documented here https://www.ncbi.nlm.nih.gov/pubmed/28967166

VariantValidator simultaneously and accurately projects genomic sequence variations onto all overlapping transcript reference sequences, and vice-versa

Alternatively, genomic sequence variation can be projected onto a specified single, or specified subset of transcript reference sequences for any given gene

Projection of sequence variations between reference sequences takes account of discrepancies between genomic and transcript reference sequences, thus ensuring an accurate prediction of the effect on encoded proteins for every gene

For sequence variations falling within the open reading frames of genes, VariantValidator automatically projects sequence variants via the transcript reference sequence onto genome builds GRCh38, GRCh37, hg38 and hg19 (HGVS format and VCF components), including projection onto relevant Alternative genomic reference sequences, the composition of which varies between patched GRC genome builds and static hg genome builds

## Pre-requisites

VariantValidator will work locally on Mac OS X or Linux-compatible computers. It can also work within a [docker container](docs/DOCKER.md).

Required software:
* MySQL
* Python 3.6 or above
* SQLite version 3.8.0 or above

Optional software:
* Postgres version 9.5 or above, Postgres 10 is not supported.

For installation instructions please see [INSTALLATION.md](docs/INSTALLATION.md). For the older python 2 version, see previous releases.

# Operation and configuration

Please see [MANUAL.md](docs/MANUAL.md). Note that the latest version is not compatible with previous releases.

## License

Please see [LICENSE.txt](LICENSE.txt)

## Terms and conditions of use
### By continuing to use VariantValidator, you accept the following terms:
All contents of VariantValidator are protected by local and international copyright laws. User inputs are submitted for the purpose of interpreting genetic and clinical information (including clinical genomic and genetic data). Outputs returned may or may not have a causal association with disease phenotypes, irrespective of stated classifications or other information presented by VariantValidator. All information returned by VariantValidator, including variant classifications, is subject to change and there is no warranty, express or implied, as to its accuracy, completeness, or fitness for a particular purpose. Use of VariantValidator and information is subject to User responsibility and discretion. Clinical decisions regarding individual patient care should be carried out in conjunction with a healthcare professional with expertise in the relevant genes and diseases. We do not accept any liability for any injury, loss or damage incurred by use of or reliance on the information provided by VariantValidator.

## Web-services additional Terms and Conditions
### By continuing to use VariantValidator web-services, you accept the following terms:
1. Secure/Remote Access. All access and use of VariantValidator must be made via a secure network.

2. Variations in Content. The University of Leicester reserve the right, in its reasonable and good faith discretion, to remove or modify materials accessed by VariantValidator or outputs provided by VariantValidator because such materials contain errors or could be subject to an infringement or other adverse claim by a third party.

3. Remedial Action. Without limiting the above, University may suspend delivery of the Service if it reasonably determines a User’s failure to comply with this Agreement may cause irreparable harm to it.

4. Service Level. University will use reasonable efforts to provide access to the Service on a continuous 24/7 basis (except for regularly scheduled maintenance when Service may be suspended) and free from viruses or other harmful software. University shall not be liable for any failure or delay or interruption in the Service or failure of any equipment or telecommunications resulting from any cause beyond University’s reasonable control. User is responsible for providing all required information.

4. No Warranty. University makes no warranty that Service is error free or that the use thereof will be uninterrupted and User acknowledges and agrees that the existence of such errors shall not constitute a breach of this Agreement. University disclaims all other warranties with respect to Service, either express or implied, including but not limited to any implied warranties relating to quality, fitness for any particular purpose or ability to achieve a particular result. Limitation of Liability. University’s total liability for any claims, losses, damages or expenses whatsoever and howsoever caused (even if caused by University’s negligence and/or breach of contract) shall be limited to the price paid to University for the products or services that are the subject of the User’s claim.

## Cite us

Hum Mutat. 2017 Oct 1. doi: 10.1002/humu.23348

VariantValidator: Accurate validation, mapping and formatting of sequence variation descriptions.

Freeman PJ, Hart RK, Gretton LJ, Brookes AJ, Dalgleish R.

> <LICENSE>
> Copyright (C) 2019 VariantValidator Contributors
>
> This program is free software: you can redistribute it and/or modify
> it under the terms of the GNU Affero General Public License as
> published by the Free Software Foundation, either version 3 of the
> License, or (at your option) any later version.
>
> This program is distributed in the hope that it will be useful,
> but WITHOUT ANY WARRANTY; without even the implied warranty of
> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
> GNU Affero General Public License for more details.
>
> You should have received a copy of the GNU Affero General Public License
> along with this program.  If not, see <https://www.gnu.org/licenses/>.
> </LICENSE>


