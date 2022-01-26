#### New code thoughts #####

from operator import contains
import re
import json
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
    if filemode="file":
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
    else:
        s = str(input)
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


if __name__ == "__main__":
    main()
