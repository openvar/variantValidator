# Markdown for exon_numbering.py
Authors: Katie Williams (@kwi11iams) and Katherine Winfield (@kjwinfield)

This file explains the exon numbering we were trying to solve, and how are module `exon_numbering.py` goes about fixing this. 

Relevant Sources and Links
--------------------------

- The relevant Git Issue for this project can be found at [https://github.com/openvar/variantValidator/issues/258](https://github.com/openvar/variantValidator/issues/258)
- User requirements were gathered by the Git Issue found at [https://github.com/openvar/variantValidator/issues/251](https://github.com/openvar/variantValidator/issues/251)
- Aspects of this project will be used to enhance code created related to issue [https://github.com/openvar/variantValidator/issues/257](https://github.com/openvar/variantValidator/issues/257), for example, SO terms for splice sites and splice regions require knowledge of the intron exon boundaries 
- Trello Board used by Katherine and Katie to stay organised, and work together in an agile manner by encouraging an iterative approach to our work: [https://trello.com/b/t0ZkUkCK/agile-board](https://trello.com/b/t0ZkUkCK/agile-board)

Outline Exon Numbering Problem
------------------------------
The aim of this project was to incorporate an exon numbering feature into VariantValidator, such that when a variant is searched it will provide an output of the exon position of that variant. This was determined to be a useful feature for users of VariantValidator as it could automate finding the exon/intron location of a patient's variant, and avoid spending extra time having to find this out from VEP, or manually from a genome browser. More information on the User Story can be found at this [Git issue](https://github.com/openvar/variantValidator/issues/258). 

The first thing to determine the context of how the exon numbering should be done. Three suggestions were made:
1. Exon in the context of individual transcripts
2. Exon in the context of the whole gene
3. Classical numbering

The [Git issue](https://github.com/openvar/variantValidator/issues/251) for gathering User requirements determined that, the first proposed method, of exon numbering in the context of individual transcripts was the most useful too for users. This was also confirmed as the preffered method by asking this same question to clinical geneticists in our office. Moreover, it was also highlighted that it would be useful to provide the exons/introns the variant starts and ends in. Hence, these user requests formed the basis of our module we created to solve this problem, `exon_numbering.py`. 

`exon_numbering.py` Module
--------------------------
Due to difficulties trying to install the VariantValidator package and the accompanying databases onto our laptops, our project could not be integrated directly into VariantValidator and its databases. Instead, we created a python module, `exon_numbering.py`, that calls VariantValidator's APIs. 

For documentation on the specifics of the functions used within this module, and the parameters for them, please see the python module itself. Here, I shall outline the functionality and thought processes that govern this module. 

The main function in `exon_numbering.py`, that calls all the other functions defined in the module within itself is called `finds_exon_number(variant, genome_build)`. This function accepts two parameters as an input: a variant defined by a transcript, written in HGVS nomenclature, and the genome build of the variant, defined by either "GRCh38", "GRCh37", "hg19" or "hg18" only. The function will output the exon/intron postions of the start and end of the given variant, for each of the aligned transcripts to the input variant's transcript. 

The functionality occuring within this function is outlined here:
1. The input variant is checked by calling the VV endpoint of the API. 
2. The transcript from the input variant is used to call the VV gene2transcripts endpoints and saves all the transcript information associated with that gene. 
3. The output from 2. is then filtered so it only contains information about the input variant's transcript. This contains the exon structure for the aligned transcripts to the input variant's transcript.
4. The start and end positions are parsed from the input variant nomenclature. For SNVs, start position = end position. 
5. A loop determines the exon/intron location of the start and end points of the variant, for each aligned transcript. This is done using the exon structure found from the gene2transcripts API. 


Test File `exon_numbering_tests.py`
----------------------------------
We have designed some tests that indicate whether an input varaint outputs the expected exon/intron positions. 

Future Integration with VariantValidator
----------------------------------------
The future application of this project will be to use the code from the `exon_numbering.py` module and integrate this functionality directly into VariantValidator. This will be done by using the generated outputs of our function as an additional dictionary in the VariantValidator object, to provide the exons/introns the variant starts and ends in. This will then be displayed on the interactive website so users can make use of this information. 

