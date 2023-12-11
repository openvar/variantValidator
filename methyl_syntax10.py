# batch_variant = the string containing the variant to be validated

def methyl_syntax(my_variant):

    if "|" in my_variant:
        print(my_variant)
        if my_variant.endswith("gom"):
            raise Exception() # Prevents "|"" from splitting variant
        elif my_variant.endswith("GOM"):
            raise print("Capitalised modifications not recognised")
        elif "gom" in my_variant and not my_variant.endswith("gom"):
            raise print("Modification at incorrect position in submitted variant description")

        elif my_variant.endswith("lom"):
            raise Exception() # Prevents "|"" from splitting variant
        elif my_variant.endswith("LOM"):
            raise print("Capitalised modifications not recognised")
        elif "lom" in my_variant and not my_variant.endswith("lom"):
            raise print("Modification at incorrect position in submitted variant description")

        elif my_variant.endswith("met="):
            raise Exception() # Prevents "|"" from splitting variant
        elif my_variant.endswith("MET="):
            raise print("Capitalised modifications not recognised")
        elif "met=" in my_variant and not my_variant.endswith("met="):
            raise print("Modification at incorrect position in submitted variant description")
        else:
            return True

Test1 = "NC_000011.10::g.1999904_1999946|gom"
Test2 = "NC_000011.10::g.1999904_1999946|lom"
Test3 = "NC_000011.10::g.1999904_1999946|met="
Test4 = "NM_000719.7:c.5550G>A"
Test5 = "NM_00256.3:c.2373dupG"
Test6 = "NM_00256.3:c.2373dupG | NM_000719.7:c.5550G>A"
Test7 = "NC_000011.10|gom:g.1999904_1999946"
Test8 = "NC_000011.10::g.1999904_1999946|GOM"
Test9 = "NC_000011.10::g.1999904_1999946|MET="

list_of_tests = [Test1,Test2,Test3,Test4,Test5,Test6,Test7,Test8,Test9]

for i in list_of_tests:
    try:
        methyl_syntax(i)
        batch_queries = i.split('|')
    except:
        batch_queries = [i]

    print(batch_queries)
