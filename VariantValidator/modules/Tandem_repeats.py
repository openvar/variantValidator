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


my_variant="LRG_199:g.1_3ACT[20]"
variant2 = "LRG_199:g.1ACT[20]A"
variant3 = "LRG_199:g.1AC[20]"
variant4 = "LRG_199t1:c.1ACT[20]"
variant5 = "LRG_199t1:c.1AC[20]"

v6 = "lrg199c.ACT(20)"
v7 = "LRG199c.1A[1_2]"

if "[" or "]" in my_variant:
    if ":" not in my_variant:
        print(": character not found")
    else:
        prefix = my_variant.split(":")[0]
        print(prefix)
        # Find whether genomic or coding
        # Capture everything between : and ., backslash as . special char in re
        variant_type = re.search(':(.*?)\.', my_variant)
        print(variant_type.group(1))
        # Get g or c position
        before = my_variant.split(":")[1]
        #Using the split sequence simplifies the regex and allows for better extraction
        position = re.search(r'\.(.*?)[ACTG]+', before)
        position_selection = position.group(1)
        print(position_selection)
        if "_" in position_selection:
            repeat_seq = re.search(r'[0-9]+_[0-9]+(.*?)\[', before)
        else:
            repeat_seq = re.search(r'[0-9]+(.*?)\[', before)
        print("repseq")
        print(repeat_seq.group(1))
        #after = my_variant.split(".")[1].split("[")
        #print(after)


#if __name__ == "__main__":
#    main()
