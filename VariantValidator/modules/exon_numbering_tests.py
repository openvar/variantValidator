"""
Exon_numbering_tests

Authors: Katie Williams (@kwi11iams) and Katherine Winfield (@kjwinfield)

This code runs tests on exon_numbering.py to check the outputs are as expected

"""
from exon_numbering import finds_exon_number, request_sequence

#define some variants to test with 
test_variant_1 = "NM_007294.3:c.1067A>G"
test_variant_2 = "NM_000088.3:c.642+1G>A"
test_variant_3 = "NM_000094.3:c.6751-3_6751-2del"
test_variant_4 = "NM_000088.3:c.589G>T" #this is an exon boundary and maps to 715
test_variant_5 = "NM_000088.3:c.642del" #this should map to 768
test_variant_6 = "NM_000094.3:c.6751-3_6753del" #starts in intron, ends in exon

#run the tests, reporting failure if the output is not what is expected
test_1 = finds_exon_number(test_variant_1)
if test_1['NC_000017.10']['start_exon'] != '14' and test_1['NG_005905.2']['end_exon'] != '10':
    print("Failed: Test variant 1, " + test_variant_1)
else:
    print("Passed: Test Variant 1, " + test_variant_1)
    print("Output dictionary: ")
    print(test_1)

test_2 = finds_exon_number(test_variant_2)
if test_2['NC_000017.10']['start_exon'] != '44i' and test_2['NG_007400.1']['end_exon'] != '8i':
    print("Failed: test variant 2")
else:
    print("Passed: Test Variant 2, " + test_variant_2)
    print("Output dictionary: ")
    print(test_2)

test_3 = finds_exon_number(test_variant_3)
if test_3['NC_000003.11']['start_exon'] != '32i' and test_3['NG_007065.1']['end_exon'] != '85i':
    print("Failed: test variant 3")

test_4 = finds_exon_number(test_variant_4)
if test_4['NC_000017.10']['start_exon'] != '44' and test_4['NG_007400.1']['end_exon'] != '8':
    print("Failed: test variant 4")

test_5 = finds_exon_number(test_variant_5)
if test_5['NC_000017.10']['start_exon'] != '44' and test_5['NG_007400.1']['end_exon'] != '8':
    print("Failed: test variant 5")

test_6 = finds_exon_number(test_variant_6)
if test_6['NC_000003.11']['start_exon'] != '32i' and test_6['NC_000003.11']['end_exon'] != '33':
    print("Failed: test variant 6")
else:
    print("Passed: Test Variant 6, " + test_variant_6)
    print("Output dictionary: ")
    print(test_6)
