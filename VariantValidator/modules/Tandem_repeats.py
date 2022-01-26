#### New code thoughts #####

from operator import contains
import re
import json
from turtle import position
from xmlrpc.client import boolean
import VariantValidator
from numpy import var
import vvhgvs

#Match variants that are correctly formatted expanded to deal with variants

def tandem_repeats_handler(var_input: str):
    """[Take string input and return require output for variant validator]
    param: input: str [String of HGVS formatted variant]

    Return: Information about variant (i.e. coordinates for each genome build,
    other info.)
    """
    #old output
    vval = VariantValidator.Validator()
    variant = str(var_input)
    genome_build = 'GRCh37'  # global variable to assign value
    select_transcripts = 'all'
    validate = vval.validate(variant, genome_build, select_transcripts)
    validation = validate.format_as_dict(with_meta=True)
    print(json.dumps(validation, sort_keys=True, indent=4, separators=(',', ': ')))
    #HGVS parser
    # parse the genomic variant into a Python structure
    hdp = vvhgvs.dataproviders.uta.connect()
    hp = vvhgvs.parser.Parser()
    var_hgvs = hp.parse_hgvs_variant(variant)
    print(var_hgvs.posedit.pos)
    print(var_hgvs.posedit)
    hn = vvhgvs.normalizer.Normalizer(hdp)
    hn.normalize(hp.parse_hgvs_variant(variant))


def check_variants(input: str, filemode: boolean):
    REGEX1 = r":(g|c)\.[0-9]+_[0-9]+([ACTG]|[actg])+\[[0-9]+\]"
    REGEX2 = r":(g|c)\.[0-9]+([ACTG]|[actg])+\[[0-9]+\]"  # selects more variants
    #  importing test file
    with open(file=input) as file:
        for line in file:
            s = str(line)
            raw_s = r'{}'.format(s)
            if re.search(pattern=REGEX1, string=raw_s):
                print(f"{raw_s}-True 1")
                tandem_repeats_handler(var_input=raw_s)
            elif re.search(pattern=REGEX2, string=raw_s):
                print(f"{raw_s}-True 2")
                tandem_repeats_handler(var_input=raw_s)
            else:
                print("False")  # change to raise error


def main():
    check_variants(test_file = "/home/rswilson1/Documents/Programming_2021/variantValidator/test_variants.txt")


import re

variant1="LRG_199:g.1ACT[20]"
variant2 = "LRG_199:g.1ACT[20]A"
variant3 = "LRG_199:g.1AC[20]"
variant4 = "LRG_199t1:c.1ACT[20]"
variant5 = "LRG_199t1:c.1AC[20]"
variant6 = "LRG_199t1:c.1act[20]"
variant7 = "LRG199c.1A[1_2]"
variant8 = "LRG_199:g.13ACT[20]"
variant9 = "LRG_199:g.13_25ACTG[1]"
variant10 = "LRG199:g.13_125ACTG[1]"
# Other types not LRG
variant11 = "ENSG00000198947.15:g.1ACT[10]"
variant12 = "ENST00000357033.8:c.13AC[22]"

def check_expanded_repeat(my_variant):
    if "[" or "]" in my_variant:
        if ":" not in my_variant:
            print("Unable to identify a colon (:) in the variant description. A colon is required in HGVS variant")
        else:
            prefix = my_variant.split(":")[0]
            print(f'Variant prefix: {prefix}')
            if re.search('^LRG', my_variant):
                # Check if underscore after LRG is included
                if "_" not in prefix:
                    # Add in underscore between LRG and number
                    prefix = re.sub(r"(?i)(?<=[a-z])(?=\d)",'_', prefix)
                    print(f'Updated prefix: {prefix}')
            # Find whether genomic or coding
            variant_type = re.search(':(.*?)\.', my_variant)
            print(f'Variant type: {variant_type.group(1)}')
            # Get g or c position
            var_position = re.search('\.(.*?)[ACTG]', my_variant)
            print(f'Variant position: {var_position.group(1)}')
            if "_" in var_position.group(1):
                start_range, end_range = var_position.group(1).split("_")
                rep_seq = re.search('\.[0-9]*_[0-9]*(.*?)\[', my_variant)
                print(f'Repeated sequence: {rep_seq.group(1)}')
            else: 
                rep_seq = re.search('\.[0-9]*(.*?)\[', my_variant)
                print(f'Repeat seq without range: {rep_seq.group(1)}')
            # Get number of unit repeats
            dup_no = re.search('\[(.*?)\]', my_variant)
            print(f'Number of unit repeats: {dup_no.group(1)}')
            # Get anything after ] to check
            after_brac = re.search('\](.*)',my_variant)
            print(f'Anything after bracket: {after_brac.group(1)}')

print(check_expanded_repeat(variant4))


#if __name__ == "__main__":
#    main()
