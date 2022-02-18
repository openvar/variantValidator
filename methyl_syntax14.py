



"""

create python dictionary of variable test cases

"""

test_cases =  {"test1":"NC_000011.10::g.1999904_1999946|gom",

               "test2":"NC_000011.10::g.1999904_1999946|lom",

               "test3":"NC_000011.10::g.1999904_1999946|met=" ,

               "test4":"NM_000719.7:c.5550G>A",

               "test5":"NM_00256.3:c.2373dupG",

               "test6":"NM_00256.3:c.2373dupG | NM_000719.7:c.5550G>A",

               "test7":"NC_000011.10|gom:g.1999904_1999946",

               "test8":"NM_00256.3:c.2373dupG | NM_000719.7:c.5550G>A | NC_000011.10|gom:g.1999904_1999946 | NC_000011.10::g.1999904_1999946|gom",

               "test9":"NM_00256.3:c.2373dupG | NM_000719.7:c.5550G>A | NC_000011.10::g.1999904_1999946|gom",

               "test10":"NM_00256.3:c.2373dupG | NM_000719.7:c.5550G>A | NC_000011.10::g.1999904_1999946|GOM",

               "test11":"NC_000011.10::g.1999904_1999946|GOM"}



"""

Custom error exceptions

MethylVariantSyntaxError - will raise an error that syntax is not correct for methylation variant

MethylVariantSyntaxCorrect -  check methyl variant syntax are correct and triggers the downstream processes to continue

CaseSensitiveError - outputs an error that methylation variant is in uppercase and must be in lowercase

"""


class MethylVariantSyntaxError(Exception):


    def __init__(self, my_variant):

        self.message = f"{my_variant} is wrong syntax"

        super().__init__(self.message)



class MethylVariantSyntaxCorrect(Exception):

    pass


class CaseSensitiveError(Exception):


    def __init__(self, my_variant):

        self.message = f"the suffix of {my_variant} must be in lowercase"

        super().__init__(self.message)


"""

methyl_syntax_check: Function that loops through methylation variant and checks they are formatted correctly

"""

def methyl_syntax_check(my_variant):


    if my_variant.endswith("gom") or my_variant.endswith("lom") or my_variant.endswith("met="):

        raise MethylVariantSyntaxCorrect(my_variant)

    if my_variant.endswith("GOM") or my_variant.endswith("LOM") or my_variant.endswith("MET="):

        raise CaseSensitiveError(my_variant)

    elif "gom" in my_variant and not my_variant.endswith("gom"):

        raise MethylVariantSyntaxError(my_variant)

    elif "lom" in my_variant and not my_variant.endswith("lom"):

        raise MethylVariantSyntaxError(my_variant)

    elif "met=" in my_variant and not my_variant.endswith("met="):

        raise MethylVariantSyntaxError(my_variant)


"""

get_unique_variant : Function which takes the variant entry stored by variant validator and identifies

the variant that may potentially be methylation variant, this then triggers the methyl_syntax function for sanity check.

"""


def get_unique_variant(variant_entry):

    variant_input = variant_entry.split(" | ")
#    print(variant_input)

    for variant in variant_input:

        try:

            if "|" in variant:
                print(variant)

                methylation_variant = variant

                methyl_syntax_check(methylation_variant)

        except MethylVariantSyntaxCorrect:

            pass

    return True


#get_unique_variant(test_cases["test10"])