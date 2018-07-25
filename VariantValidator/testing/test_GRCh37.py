# -*- coding: utf-8 -*-
"""
This script contains a list of variants required to test the VariantValidator API for consistent validation of variants
prior to committing to a repo.

There are two lists, curated and new
    Curated are variants that have been checked and form the stable testing list
    New are variants which will be written to a file for manual curation.

    Once curated, new variants are removed from the New list and added to the currated list
"""

# Import modules
import json
import io
import os
import sys
try:
    import variantanalyser.functions as functions
except ImportError:
    parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.sys.path.insert(0, parentdir)
    from VariantValidator import variantValidator
except AttributeError:
    parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.sys.path.insert(0, parentdir)
    from VariantValidator import variantValidator


# Error types
class variantValidatorError(Exception):
    pass


def get_curated_variants():
    curated_variant_list = []
    return curated_variant_list


def get_new_variants():
    new_variant_list = ["NC_000016.9:g.2099572TC>T",
                        "NM_000088.3:c.589GG>CT",
                        "NM_000094.3:c.6751-2_6751-3del",
                        "COL5A1:c.5071A>T",
                        "NG_007400.1:c.5071A>T",
                        "chr16:15832508_15832509delinsAC",
                        "NM_000088.3:c.589-1GG>G",
                        "NM_000088.3:c.642+1GT>G",
                        "NM_000088.3:c.589-2AG>G",
                        "NC_000017.10:g.48279242G>T",
                        "NM_000500.7:c.-107-19C>T",
                        "NM_000518.4:c.-130C>T",
                        "NR_138595.1:n.-810C>T",
                        "NR_138595.1:n.1-810C>T",
                        "NC_000017.10:g.48261457_48261463TTATGTT=",
                        "NC_000017.10:g.48275363C>A",
                        "NM_000088.3:c.589-1G>T",
                        "NM_000088.3:c.591_593inv",
                        "11-5248232-T-A",
                        "NG_007400.1(NM_000088.3):c.589-1G>T",
                        "1:150550916G>A",
                        "1-150550916-G-A",
                        "NG_008123.1(LEPRE1_v003):c.2055+18G>A",
                        "NG_008123.1:c.2055+18G>A",
                        "NG_008123.1(NM_022356.3):c.2055+18G>A",
                        "NM_021983.4:c.490G>C",
                        "NM_032470.3:c.4del",
                        "NM_001194958.2:c.20C>A",
                        "NM_000022.2:c.534A>G",
                        "HSCHR6_MHC_SSTO_CTG1-3852542-C-G",
                        "NM_000368.4:c.363+1dupG",
                        "NM_000368.4:c.363dupG",
                        "NM_000089.3:c.1033_1035delGTT",
                        "NM_000089.3:c.1035_1035+2delTGT",
                        "NM_000088.3:c.2023_2028delGCAAGA",
                        "NM_000089.3:c.938-1delG",
                        "NM_000088.3:c.589G=",
                        "NM_000088.3:c.642A=",
                        "NM_000088.3:c.642+1GG>G",
                        "NM_000088.3:c.589-2GG>G",
                        "NM_000088.3:c.589-6_589-5insTTTT",
                        "NM_000088.3:c.642+3_642+4insAAAA",
                        "NM_000088.3:c.589-4_589-3insTT",
                        "NM_000088.3:c.589-8del",
                        "NM_000527.4:c.-187_-185delCTC",
                        "NM_206933.2:c.6317C>G",
                        "NC_000013.10:g.32929387T>C",
                        "NM_015102.3:c.2818-2T>A",
                        "19-41123094-G-GG",
                        "15-72105928-AC-A",
                        "12-122064773-CCCGCCA-C",
                        "12-122064774-CCGCCA-CCGCCA",
                        "12-122064773-CCCGCCACCGCCACCGC-CCCGCCACCGCCGCCGTC",
                        "NC_000012.11:g.122064777C>A",
                        "NC_000012.11:g.122064776delG",
                        "NC_000012.11:g.122064776dupG",
                        "NC_000012.11:g.122064776_122064777insTTT",
                        "NC_000012.11:g.122064772_122064775del",
                        "NC_000012.11:g.122064772_122064775dup",
                        "NC_000012.11:g.122064773_122064774insTTTT",
                        "NC_000012.11:g.122064772_122064777del",
                        "NC_000012.11:g.122064772_122064777dup",
                        "NC_000012.11:g.122064779_122064782dup",
                        "NC_000012.11:g.122064772_122064782del",
                        "NC_000002.11:g.95847041_95847043GCG=",
                        "NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=",
                        "NC_000003.11:g.14561629_14561630GC=",
                        "NC_000003.11:g.14561629_14561630insG",
                        "NC_000004.11:g.140811111_140811122del",
                        "NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=",
                        "NC_000004.11:g.140811117_140811122del",
                        "NC_000004.11:g.140811111_140811117del",
                        "NC_000004.11:g.140811117C>A",
                        "NC_000002.11:g.73675227_73675228insCTC",
                        "9-136132908-T-TC",
                        "9-136132908-TAC-TCA",
                        "9-136132908-TA-TA",
                        "NM_020469.2:c.258delG",
                        "NM_020469.2:c.260_262TGA=",
                        "NM_020469.2:c.261delG",
                        "NM_020469.2:c.261dupG",
                        "NM_020469.2:c.261_262insTT",
                        "NC_000019.10:g.50378563_50378564insTAC",
                        "NC_000019.10:g.50378563_50378564insC",
                        "NC_000019.10:g.50378564_50378565insTACA",
                        "NC_000019.10:g.50378565_50378567dup",
                        "NC_000019.10:g.50378563_50378564=",
                        "NC_000019.10:g.50378563_50378564insTCGG",
                        "NC_000019.10:g.50378563delinsTTAC",
                        "NC_000019.10:g.50378563_50378564insTAAC",
                        "NC_000019.10:g.50378562_50378565del",
                        "NC_000019.10:g.50378562_50378565delinsTC",
                        "NC_000007.14:g.149779575_149779577delinsT",
                        "NC_000007.14:g.149779575_149779577=",
                        "NC_000007.14:g.149779576_149779578del",
                        "NC_000007.14:g.149779577del",
                        "NC_000007.14:g.149779573_149779579del",
                        "NC_000007.14:g.149779573_149779579delinsCA",
                        "NM_000088.3:c.590_591inv",
                        "NM_024989.3:c.1778_1779inv",
                        "NM_032815.3:c.555_556inv",
                        "NM_006138.4:c.3_4inv",
                        "NM_000038.5:c.3927_3928delAAinsTT",
                        "NM_001034853.1:c.2847_2848delAGinsCT",
                        "NM_000088.3:c.4392_*2inv",
                        "NM_000088.3:c.4392_*5inv",
                        "NM_000088.3:c.4390_*7inv",
                        "NM_005732.3:c.2923-5insT",
                        "NM_198283.1(EYS):c.*743120C>T",
                        "NM_133379.4(TTN):c.*265+26591C>T",
                        "NM_000088.3:c.589-2_589-1AG>G",
                        "NM_000088.3:c.642+1_642+2delGTinsG",
                        "NM_004415.3:c.1-1insA",
                        "NM_000273.2:c.1-5028_253del",
                        "NM_002929.2:c.1006C>T",
                        "NR_125367.1:n.167+18165G>A",
                        "NM_006005.3:c.3071_3073delinsTTA",
                        "NM_000089.3:n.1504_1506del",
                        "NC_012920.1:m.1011C>T",
                        "NC_000006.11:g.90403795G=",
                        "1-169519049-T-.",
                        "NC_000005.9:g.35058667_35058668AG=",
                        "NM_000251.1:c.1296_1348del",
                        "NM_000088.3:c.2023_2028del",
                        "NM_000088.3:c.2024_2028+1del",
                        "ENST00000450616.1:n.31+1G>C",
                        "ENST00000491747:c.5071A>T",
                        "NM_000088.3:c.589G>T",
                        "NG_007400.1:g.8638G>T",
                        "LRG_1:g.8638G>T",
                        "LRG_1t1:c.589G>T",
                        "chr16:g.15832508_15832509delinsAC",
                        "NG_012386.1:g.24048dupG",
                        "NM_033517.1:c.1307_1309delCGA",
                        "HG1311_PATCH-33720-CCGA-C",
                        "2-73675227-TCTC-TCTCCTC",
                        "2-73675227-TC-TC",
                        "3-14561627-AG-AGG",
                        "3-14561630-CC-CC",
                        "6-90403795-G-G",
                        "6-90403795-G-A",
                        "6-32012992-CG-C",
                        "17-48275363-C-A",
                        "17-48275364-C-A",
                        "17-48275359-GGA-TCC",
                        "7-94039128-CTTG-C",
                        "9-135800972-AC-ACC",
                        "1-43212925-C-T",
                        "HG987_PATCH-355171-C-A",
                        "20-43252915-T-C",
                        "1-216219781-A-C",
                        "2-209113113-G-A,C,T",
                        "NC_000005.9:g.35058665_35058666CA=",
                        "NC_000002.11:g.73675227_73675229delTCTinsTCTCTC",
                        "NM_000828.4:c.-2dupG",
                        "X-122318386-A-AGG",
                        "NM_000828.4:c.-2G>T",
                        "NM_000828.4:c.-2G=",
                        "X-122318386-A-AT",
                        "NM_000828.4:c.-2_-1insT",
                        "NM_000828.4:c.-3_-2insT",
                        "NM_000828.4:c.-2delGinsTT",
                        "NM_000828.4:c.-2_-1delGCinsTT",
                        "NM_000828.4:c.-3_-2delAGinsTT",
                        "15-72105929-C-C",
                        "15-72105928-AC-ATT",
                        "15-72105928-ACC-ATT",
                        "15-72105927-GACC-GTT",
                        "19-41123093-A-AG",
                        "19-41123093-A-AT",
                        "19-41123093-AG-A",
                        "19-41123093-AG-AG",
                        "NM_012309.4:c.913-5058G>A",
                        "LRG_199t1:c.2376[G>C];[G>C]",
                        "LRG_199t1:c.[2376G>C];[3103del]",
                        "LRG_199t1:c.[4358_4359del;4361_4372del]",
                        "LRG_199t1:c.2376G>C(;)3103del",
                        "LRG_199t1:c.2376[G>C];[(G>C)]",
                        "LRG_199t1:c.[2376G>C];[?]",
                        "LRG_199t1:c.[296T>G;476T=];[476T=](;)1083A>C",
                        "LRG_199t1:c.[296T>G];[476T>C](;)1083A>C(;)1406del",
                        "LRG_199t1:c.[976-20T>A;976-17_976-1dup]",
                        "1-5935162-A-T",
                        "1-12065948-C-T",
                        "1-46655125-CTCAC-C",
                        "1-68912523-TGAGCCAGAG-T",
                        "1-68912526-GCCAGAG-G",
                        "1-109817590-G-T",
                        "1-145597475-GAAGT-G",
                        "1-153791300-CTG-C",
                        "1-156104666-TTGAGAGCCGGCTGGCGGATGCGCT-TCCCC",
                        "1-156108541-G-GG",
                        "1-161279695-T-A",
                        "1-169519049-T-T",
                        "1-226125468-G-A",
                        "10-89623035-CGCA-C",
                        "11-62457852-C-A",
                        "11-108178710-A-AT",
                        "11-111735981-G-A",
                        "12-11023080-C-A",
                        "12-22018712-TC-T",
                        "12-52912946-T-C",
                        "12-103234292-TC-T",
                        "12-103311124-T-C",
                        "12-111064166-G-A",
                        "12-123738430-CA-C",
                        "13-31789169-CT-C",
                        "14-62187287-G-A",
                        "14-62188231-TT-GA",
                        "14-63174827-C-A",
                        "15-42680000-CA-C",
                        "15-42680000-CA-CAA",
                        "15-42703179-T-TTCA",
                        "15-42703179-TAG-TTCATCT",
                        "15-48782203-C-T",
                        "15-72105929-CC-C",
                        "15-89873415-G-A",
                        "16-2103394-C-T",
                        "16-3779300-C-G",
                        "16-5128843-C-G",
                        "16-74808559-C-T",
                        "16-89574804-C-A",
                        "16-89574826-A-C",
                        "16-89574914-G-GT",
                        "16-89574916-C-CGTC",
                        "16-89575009-G-A",
                        "16-89575040-C-A,CA",
                        "16-89576896-A-C",
                        "16-89576930-T-TA,TT",
                        "16-89576931-G-GTG",
                        "16-89598368-CGGCCCCCCCGGCTGTGGGAAGACGCT-C",
                        "16-89613064-AGGAGAGGCG-AT",
                        "16-89613069-AGGCGGGAGA-AT",
                        "16-89613145-C-T",
                        "17-7578194-GCAC-G",
                        "17-7578523-T-TG",
                        "17-17119692-A-C",
                        "17-41197588-GGACA-G",
                        "17-41256884-C-G",
                        "17-42991428-C-A",
                        "17-48252809-A-T",
                        "17-62022709-G-GTC",
                        "17-62022711-C-CT",
                        "17-62023005-G-GGC",
                        "17-62023006-C-A",
                        "17-62034787-G-A",
                        "18-24128261-GTCCTCC-G",
                        "19-15291774-G-A",
                        "19-15311794-A-G",
                        "19-39076592-G-A",
                        "2-50149352-T-C",
                        "2-50847195-G-A",
                        "2-71825797-C-G",
                        "2-166179712-G-C",
                        "2-166183371-A-G",
                        "2-166929889-GTCCAGGTCCT-GAC",
                        "2-166929891-CCAGGTCCT-C",
                        "2-179393504-G-T",
                        "2-185803444-TGCAGCTGCTGCAGCTGCAGCTGCA-T",
                        "2-201950249-G-T",
                        "2-238268730-C-A",
                        "21-43897396-C-T",
                        "22-30064360-G-GCGACGC",
                        "3-10188187-TGTCCCGATAG-T",
                        "3-50402127-T-G",
                        "3-50402890-G-A",
                        "3-57851007-AG-A",
                        "3-122003832-G-C",
                        "4-153332910-C-CAGG",
                        "5-1295183-G-A",
                        "5-77396835-TTTC-T",
                        "5-118811422-GGTGA-G",
                        "5-118811422-GGTGAG-G",
                        "5-131705587-CG-C",
                        "5-148406482-T-C",
                        "6-110036337-T-TCAG",
                        "6-110036337-TGAT-T",
                        "6-152651802-C-A",
                        "6-152737643-C-G",
                        "7-6026775-T-C",
                        "7-55242465-GGAATTAAGAGAAGCA-G",
                        "7-55248992-T-TTCCAGGAAGCCT",
                        "7-75932111-C-A",
                        "7-91652178-A-AAAC",
                        "7-117199644-ATCT-A",
                        "7-140453136-AC-CT",
                        "7-140453136-A-T",
                        "7-140453137-C-T",
                        "7-143013488-A-T",
                        "7-143018934-G-A",
                        "7-143048771-C-T",
                        "8-1871951-C-T",
                        "9-13112056-T-TG",
                        "9-21971208-C-A",
                        "9-35683240-T-TG",
                        "9-135796754-G-A",
                        "HG536_PATCH-10391-AC-A",
                        "HG865_PATCH-33547-G-A",
                        "HG865_PATCH-569441-G-T",
                        "HG865_PATCH-574546-C-T",
                        "HSCHR1_1_CTG31-133178-TAG-T",
                        "HSCHR6_MHC_MANN_CTG1-3848158-T-G",
                        "HSCHR6_MHC_MANN_CTG1-3851043-C-A",
                        "X-70443101-C-T",
                        "X-107845202-GACCACC-GACC,G",
                        "X-153296777-G-A"                      
                        ]
    return new_variant_list


def write_variant_files():
    ROOT = os.path.dirname(os.path.abspath(__file__))
    new_file = os.path.join(ROOT, 'new_variants.txt')
    cur_file = os.path.join(ROOT, 'curated_variants.txt')
    new_variants = get_new_variants()
    curated_variants = get_curated_variants()
    all_new = []
    for nv in new_variants:
        print 'Variant ' + ': ' + nv
        try:
            result = variantValidator.validator(nv, 'GRCh37', 'all')
            all_new.append({nv: result})
            print json.dumps(result, sort_keys=True, indent=4, separators=(',', ': '))
        except variantValidatorError as e:
            sys.exit(str(e))
    all_cur = []
    for cv in curated_variants:
        print 'Variant ' + ': ' + cv
        try:
            result = variantValidator.validator(cv, 'GRCh37', 'all')
            all_cur.append({cv: result})
            print json.dumps(result, sort_keys=True, indent=4, separators=(',', ': '))
        except variantValidatorError as e:
            sys.exit(str(e))

    if all_new != []:
        with io.open(new_file, 'w', encoding='utf-8') as nf:
            nf.write(json.dumps(all_new, sort_keys=True, indent=4, separators=(',', ': '), ensure_ascii=False))
    if all_cur != []:
        with io.open(cur_file, 'w', encoding='utf-8') as cf:
            cf.write(json.dumps(all_cur, sort_keys=True, indent=4, separators=(',', ': '), ensure_ascii=False))
