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
The aim of this project is to incorporate an exon numbering feature into VariantValidator, such that when a variant is searched VariantValidator will provide the exon or intron position of that variant. As we (both Katherine and Katie) work in clinical genetics labs, we could speak to clinical geneticists, and they confirmed that this would be a useful feature for users of VariantValidator. Automating finding the exon/intron location of a patient's variant would save time, as the clinical geneticist would no longer have to find this out from VEP, or manually from a genome browser. More information on the User Story can be found at this [Git issue](https://github.com/openvar/variantValidator/issues/258). 

The first thing to do was to determine how the exon numbering should be done. Three suggestions were made:
1. Exon number in the context of individual transcripts
2. Exon number in the context of the whole gene
3. Classical numbering

The [Git issue](https://github.com/openvar/variantValidator/issues/251) for gathering user requirements determined that the first proposed method (exon numbering in the context of individual transcripts) would be the most useful for users. When we asked this same question to our clinical geneticist colleagues, they confirmed that this method would be the most helpful way of representing exon numbering information. Additionally, it was highlighted that it would be necessary to provide the exon/intron location for start and end of the variant, as some variants span across an exon boundary. Hence, these user requirements formed the basis of the module we created to solve this problem, `exon_numbering.py`. 

Module `exon_numbering.py` 
--------------------------
Due to difficulties trying to install the VariantValidator package and the accompanying databases onto our laptops, our project could not be integrated directly into VariantValidator and its databases. Instead, we created a python module, `exon_numbering.py`, that calls VariantValidator's APIs. 

For documentation on the specifics of the functions used within this module, and the parameters for them, please see the python module itself. This .md file outlines the functionality and thought processes that govern this module. 

The main function in `exon_numbering.py` is called `finds_exon_number(variant, genome_build)`, and it calls all the other functions defined in the module. This function accepts two parameters as an input:
- A variant written in HGVS nomenclature
- The genome build of the variant, defined by either "GRCh38", "GRCh37", "hg19" or "hg18" only.

This function outputs the exon/intron postions of the start and end of the given variant. It creates a dictionary of these positions for each chromosomal or gene reference sequence aligned to the input variant's transcript. 

The functionality of  `finds_exon_number` is outlined below:
1. The input variant is checked by calling the VV endpoint of the API. 
2. The transcript from the input variant is used to call the VV gene2transcripts endpoints and saves all the transcript information associated with that gene. 
3. The output from 2. is then filtered so it only contains information about the input variant's transcript. This contains the exon structure for the aligned transcripts to the input variant's transcript.
4. The start and end positions are parsed from the input variant nomenclature. For SNVs, start position = end position. 
5. A loop determines the exon/intron location of the start and end points of the variant, for each aligned transcript. This is done using the exon structure found from the gene2transcripts API. 


Module `exon_numbering_tests.py`
----------------------------------
The module `exon_numbering_tests.py` runs automated tests. These tests indicate whether an input variant outputs the expected exon/intron positions. This allows us to test `finds_exon_number` and confirm that it returns the correct output.

Future Integration with VariantValidator
----------------------------------------
The future application of this project will be to integrate the functionality of the `exon_numbering.py` module directly into VariantValidator. This will be done by using the generated outputs of our function as an additional dictionary in the VariantValidator object, to provide the exon/intron location that the variant starts and ends in. This will then be displayed on the interactive website so users can make use of this information. 

