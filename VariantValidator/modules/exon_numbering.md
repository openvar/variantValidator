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
The aim of this project was to incorporate an exon numbering feature into VariantValidator, such that when a variant is searched it will provide an output of the exon position of that variant. The first thing to determine the context of how the exon numbering should be done. Three suggestions were made:
1. Exon in the context of individual transcripts
2. Exon in the context of the whole gene
3. Classical numbering

The [Git issue](https://github.com/openvar/variantValidator/issues/251) for gathering User requirements determined that, the first proposed method, of exon numbering in the context of individual transcripts was the most useful too for users. This was also confirmed as the preffered method by asking this same question to clinical geneticists in our office. 

`exon_numbering.py` Module
--------------------------


Test File `exon_numbering_tests.py`
----------------------------------


Future Integration with VariantValidator
----------------------------------------
The future application of this project will be to use the code from the `exon_numbering.py` module and integrate this functionality directly into VariantValidator. This will be done by using the generated outputs of our function as an additional dictionary in the VariantValidator object, to provide the exons/introns the variant starts and ends in. This will then be displayed on the interactive website so users can make use of this information. 

