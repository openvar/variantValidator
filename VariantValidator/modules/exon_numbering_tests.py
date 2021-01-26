"""
Exon_numbering_tests

Authors: Katie Williams (@kwi11iams) and Katherine Winfield (@kjwinfield)

This code runs tests on exon_numbering.py to check the outputs are as expected

"""
from exon_numbering import finds_exon_number, request_sequence

#Testing
#define some variants to test with 
test_variant_2 = "NM_007294.3:c.1067A>G"
test_variant_1  = 'NM_000088.3:c.642+1GG>G'
test_variant_3 = "NM_000088.3:c.642+1G>A"
test_variant_4 = "NM_000094.3:c.6751-3_6751-2del"
test_variant_5 = "NM_007294.3:c.5426-2del"
test_variant_6 = "NM_000088.3:c.589G>T" #this is an exon boundary and maps to 715
test_variant_7 = "NM_000088.3:c.642del" #this should map to 768
test_variant_8 = "NM_000094.3:c.6751-3_6753del" #starts in intron, ends in exon

#test for our variant
#finds_exon_number(test_variant_2)
#print(finds_exon_number(test_variant_1))
#print(finds_exon_number(test_variant_3))
#print(finds_exon_number(test_variant_4))
#print(finds_exon_number(test_variant_5))
#print(finds_exon_number(test_variant_5))
print(finds_exon_number(test_variant_8))