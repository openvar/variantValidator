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
    from variantValidator import variantValidator
except AttributeError:
    parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.sys.path.insert(0, parentdir)
    from variantValidator import variantValidator


# Error types
class variantValidatorError(Exception):
    pass


def get_curated_variants():
    curated_variant_list = []
    return curated_variant_list


def get_new_variants():
    new_variant_list = ['NC_000016.9:g.2099572TC>T',
                        'NM_000088.3:c.589GG>CT',
                        'NM_000094.3:c.6751-2_6751-3del',
                        'COL5A1:c.5071A>T',
                        'NG_007400.1:c.5071A>T',
                        'chr16:15832508_15832509delinsAC',
                        'NM_000088.3:c.589-1GG>G',
                        'NM_000088.3:c.642+1GT>G',
                        'NM_000088.3:c.589-2AG>G',
                        'NC_000017.10:g.48279242G>T',
                        'NM_000500.7:c.-107-19C>T',
                        'NM_000518.4:c.-130C>T',
                        'NR_138595.1:n.-810C>T',
                        'NR_138595.1:n.1-810C>T',
                        'NC_000017.10:g.48261457_48261463TTATGTT=',
                        'NC_000017.10:g.48275363C>A',
                        'NM_000088.3:c.589-1G>T',
                        'NM_000088.3:c.591_593inv',
                        '11-5248232-T-A',
                        'NG_007400.1(NM_000088.3):c.589-1G>T',
                        '1:150550916G>A',
                        '1-150550916-G-A',
                        'NG_008123.1(LEPRE1_v003):c.2055+18G>A',
                        'NG_008123.1:c.2055+18G>A',
                        'NG_008123.1(NM_022356.3):c.2055+18G>A',
                        'NM_021983.4:c.490G>C',
                        'NM_032470.3:c.4del',
                        'NM_001194958.2:c.20C>A',
                        'NM_000022.2:c.534A>G',
                        'HSCHR6_MHC_SSTO_CTG1-3852542-C-G',
                        'NM_000368.4:c.363+1dupG',
                        'NM_000368.4:c.363dupG',
                        'NM_000089.3:c.1033_1035delGTT',
                        'NM_000089.3:c.1035_1035+2delTGT',
                        'NM_000088.3:c.2023_2028delGCAAGA',
                        'NM_000089.3:c.938-1delG',
                        'NM_000088.3:c.589G=',
                        'NM_000088.3:c.642A=',
                        'NM_000088.3:c.642+1GG>G',
                        'NM_000088.3:c.589-2GG>G',
                        'NM_000088.3:c.589-6_589-5insTTTT',
                        'NM_000088.3:c.642+3_642+4insAAAA',
                        'NM_000088.3:c.589-4_589-3insTT',
                        'NM_000088.3:c.589-8del',
                        'NM_000527.4:c.-187_-185delCTC',
                        'NM_206933.2:c.6317C>G',
                        'NC_000013.10:g.32929387T>C',
                        'NM_015102.3:c.2818-2T>A',
                        '19-41123094-G-GG',
                        '15-72105928-AC-A',
                        '12-122064773-CCCGCCA-C',
                        '12-122064774-CCGCCA-CCGCCA',
                        '12-122064773-CCCGCCACCGCCACCGC-CCCGCCACCGCCGCCGTC',
                        'NC_000012.11:g.122064777C>A',
                        'NC_000012.11:g.122064776delG',
                        'NC_000012.11:g.122064776dupG',
                        'NC_000012.11:g.122064776_122064777insTTT',
                        'NC_000012.11:g.122064772_122064775del',
                        'NC_000012.11:g.122064772_122064775dup',
                        'NC_000012.11:g.122064773_122064774insTTTT',
                        'NC_000012.11:g.122064772_122064777del',
                        'NC_000012.11:g.122064772_122064777dup',
                        'NC_000012.11:g.122064779_122064782dup',
                        'NC_000012.11:g.122064772_122064782del',
                        'NC_000002.11:g.95847041_95847043GCG=',
                        'NC_000017.10:g.5286863_5286889AGTGTTTGGAATTTTCTGTTCATATAG=',
                        'NC_000003.11:g.14561629_14561630GC=',
                        'NC_000003.11:g.14561629_14561630insG',
                        'NC_000004.11:g.140811111_140811122del',
                        'NC_000004.11:g.140811111_140811122CTGCTGCTGCTG=',
                        'NC_000004.11:g.140811117_140811122del',
                        'NC_000004.11:g.140811111_140811117del',
                        'NC_000004.11:g.140811117C>A',
                        'NC_000002.11:g.73675227_73675228insCTC',
                        '9-136132908-T-TC',
                        '9-136132908-TAC-TCA',
                        '9-136132908-TA-TA',
                        'NM_020469.2:c.258delG',
                        'NM_020469.2:c.260_262TGA=',
                        'NM_020469.2:c.261delG',
                        'NM_020469.2:c.261dupG',
                        'NM_020469.2:c.261_262insTT',
                        'NM_000088.3:c.590_591inv',
                        'NM_024989.3:c.1778_1779inv',
                        'NM_032815.3:c.555_556inv',
                        'NM_006138.4:c.3_4inv',
                        'NM_000038.5:c.3927_3928delAAinsTT',
                        'NM_001034853.1:c.2847_2848delAGinsCT',
                        'NM_000088.3:c.4392_*2inv',
                        'NM_000088.3:c.4392_*5inv',
                        'NM_000088.3:c.4390_*7inv',
                        'NM_005732.3:c.2923-5insT',
                        'NM_198283.1(EYS):c.*743120C>T',
                        'NM_133379.4(TTN):c.*265+26591C>T',
                        'NM_000088.3:c.589-2_589-1AG>G',
                        'NM_000088.3:c.642+1_642+2delGTinsG',
                        'NM_004415.3:c.1-1insA',
                        'NM_000273.2:c.1-5028_253del',
                        'NM_002929.2:c.1006C>T',
                        'NR_125367.1:n.167+18165G>A',
                        'NM_006005.3:c.3071_3073delinsTTA',
                        'NM_000089.3:n.1504_1506del',
                        'NC_012920.1:m.1011C>T',
                        'NC_000006.11:g.90403795G=',
                        '1-169519049-T-.',
                        'NC_000005.9:g.35058667_35058668AG=',
                        'NM_000251.1:c.1296_1348del',
                        'NM_000088.3:c.2023_2028del',
                        'NM_000088.3:c.2024_2028+1del',
                        'ENST00000450616.1:n.31+1G>C',
                        'ENST00000491747:c.5071A>T',
                        'NM_000088.3:c.589G>T',
                        'NG_007400.1:g.8638G>T',
                        'LRG_1:g.8638G>T',
                        'LRG_1t1:c.589G>T',
                        'chr16:g.15832508_15832509delinsAC',
                        'NG_012386.1:g.24048dupG',
                        'NM_033517.1:c.1307_1309delCGA',
                        'HG1311_PATCH-33720-CCGA-C',
                        '2-73675227-TCTC-TCTCCTC',
                        '2-73675227-TC-TC',
                        '3-14561627-AG-AGG',
                        '3-14561630-CC-CC',
                        '6-90403795-G-G',
                        '6-90403795-G-A',
                        '6-32012992-CG-C',
                        '17-48275363-C-A',
                        '17-48275364-C-A',
                        '17-48275359-GGA-TCC',
                        '7-94039128-CTTG-C',
                        '9-135800972-AC-ACC',
                        '1-43212925-C-T',
                        'HG987_PATCH-355171-C-A',
                        '20-43252915-T-C',
                        '1-216219781-A-C',
                        '2-209113113-G-A,C,T',
                        'NC_000005.9:g.35058665_35058666CA=',
                        'NC_000002.11:g.73675227_73675229delTCTinsTCTCTC',
                        'NM_000828.4:c.-2dupG',
                        'X-122318386-A-AGG',
                        'NM_000828.4:c.-2G>T',
                        'NM_000828.4:c.-2G=',
                        'X-122318386-A-AT',
                        'NM_000828.4:c.-2_-1insT',
                        'NM_000828.4:c.-3_-2insT',
                        'NM_000828.4:c.-2delGinsTT',
                        'NM_000828.4:c.-2_-1delGCinsTT',
                        'NM_000828.4:c.-3_-2delAGinsTT',
                        '15-72105929-C-C',
                        '15-72105928-AC-ATT',
                        '15-72105928-ACC-ATT',
                        '15-72105927-GACC-GTT',
                        '19-41123093-A-AG',
                        '19-41123093-A-AT',
                        '19-41123093-AG-A',
                        '19-41123093-AG-AG',
                        'NM_012309.4:c.913-5058G>A',
                        'LRG_199t1:c.2376[G>C];[G>C]',
                        'LRG_199t1:c.[2376G>C];[3103del]',
                        'LRG_199t1:c.[4358_4359del;4361_4372del]',
                        'LRG_199t1:c.2376G>C(;)3103del',
                        'LRG_199t1:c.2376[G>C];[(G>C)]',
                        'LRG_199t1:c.[2376G>C];[?]',
                        'LRG_199t1:c.[296T>G;476T=];[476T=](;)1083A>C',
                        'LRG_199t1:c.[296T>G];[476T>C](;)1083A>C(;)1406del',
                        'LRG_199t1:c.[976-20T>A;976-17_976-1dup]',
                        'chr2:g.[29443695G>T];[29443695G>C;29443697A>G]',
                        'chr7:g.87053221C>T',
                        'chr9:g.133738306G>A',
                        'chr2:g.[29443695G>T];[29443695G>C;29443697A>G]'
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
