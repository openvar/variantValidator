"""
Handling Methylation Variants
Task was to identify methylation variants and stop them from being split at the
"|" character
while raising exceptions on syntax errors
"""


"""
Custom error exceptions

MethylSyntaxError - Will raise an error that syntax is incorrect for the
modified variant
MethylSyntaxCorrect - If modification syntax is correct this triggers
continuation of the downstream processes
CaseSensitiveError - Throws an error if the modified variant is in uppercase
"""


class MethylSyntaxError(Exception):
    def __init__(self, my_variant):
        self.message = f"Modification at incorrect position in {my_variant}"
        super().__init__(self.message)


class MethylSyntaxCorrect(Exception):
    pass


class CaseSensitiveError(Exception):
    def __init__(self, my_variant):
        self.message = f"Capitalised modification in {my_variant} not recognised"
        super().__init__(self.message)


"""
methyl_syntax_check: Function checks that modified variant description is
in correct format
"""


def methyl_syntax_check(my_variant):

    if my_variant.endswith("gom") or my_variant.endswith("lom") or my_variant.endswith("met="):
        raise MethylSyntaxCorrect(my_variant)
    if my_variant.endswith("GOM") or my_variant.endswith("LOM") or my_variant.endswith("MET="):
        raise CaseSensitiveError(my_variant)
    elif "gom" in my_variant and not my_variant.endswith("gom"):
        raise MethylSyntaxError(my_variant)
    elif "lom" in my_variant and not my_variant.endswith("lom"):
        raise MethylSyntaxError(my_variant)
    elif "met=" in my_variant and not my_variant.endswith("met="):
        raise MethylSyntaxError(my_variant)


"""
Create python dictionary of variable test cases.
"""

test_cases =  {"test1":"NC_000011.10::g.1999904_1999946|gom",

               "test2":"NC_000011.10::g.1999904_1999946|lom",

               "test3":"NC_000011.10::g.1999904_1999946|met=",

               "test4":"NM_000719.7:c.5550G>A",

               "test5":"NM_00256.3:c.2373dupG",

               "test6":"NM_00256.3:c.2373dupG | NM_000719.7:c.5550G>A",

               "test7":"NC_000011.10|gom:g.1999904_1999946",

               "test8":"NM_00256.3:c.2373dupG | NM_000719.7:c.5550G>A | NC_000011.10|gom:g.1999904_1999946 | NC_000011.10::g.1999904_1999946|gom",

               "test9":"NM_00256.3:c.2373dupG | NM_000719.7:c.5550G>A | NC_000011.10::g.1999904_1999946|gom",

               "test10":"NM_00256.3:c.2373dupG | NM_000719.7:c.5550G>A | NC_000011.10::g.1999904_1999946|GOM",

               "test11":"NC_000011.10::g.1999904_1999946|GOM"}

"""
get_unique_variant : Function identifies a modified variant and triggers
the methyl_syntax function for sanity check.
"""


def get_unique_variant(variant_entry):

    # Variants split where "|" character occurs, mimicking VariantValidator

    variant_input = variant_entry.split(" | ")
#    print(variant_input)

    for variant in variant_input:
        try:
            if "|" in variant:
                print(variant)  # prints output
                methylation_variant = variant
                methyl_syntax_check(methylation_variant)
        except MethylSyntaxCorrect:
            pass
    return True


""" Runs methyl_syntax_check on test of interest.
We expect to be thrown a CaseSensitiveError when inputting test11: NC_000011.10::g.1999904_1999946|GOM"""

get_unique_variant(test_cases["test11"])

"""The test below will run each of the tests in the dictionary above through
the methyl_syntax_check function"""
"""
for test, variant_entry in test_cases.items():
    try:
        get_unique_variant(variant_entry)
    except MethylSyntaxError:
        print( "Modification at incorrect position in variant description",
        test)
    except CaseSensitiveError:
        print( "Capitalised modification in variant description not
        recognised", test)"""
