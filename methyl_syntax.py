# Pipe syntax error

# batch_variant = the string containing the variant to be validated


def methyl_syntax(my_variant):
    """
    :param my_variant:
    :return: False if no | is detected otherwise raises Exception and provides information as to where the
    | is located
    """

    if "|" in my_variant:
        print(my_variant)
        if "gom" in my_variant:
            raise Exception("message")
        elif "lom" in my_variant:
            raise Exception("message")
        elif "met=" in my_variant:
            raise Exception("message")
        #elif "gom" or "lom" or "met=" in my_variant and "gom" or "lom" or "met=" not in my_variant.end:
        #    raise Error("Variant description is not in accepted format")
        else:
            return True


Test1 = "NC_000011.10::g.1999904_1999946|gom"
Test2 = "NC_000011.10::g.1999904_1999946|lom"
Test3 = "NC_000011.10::g.1999904_1999946|met="
Test4 = "NM_000719.7:c.5550G>A"
Test5 = "NM_00256.3:c.2373dupG"
Test6 = "NM_00256.3:c.2373dupG | NM_000719.7:c.5550G>A"

list_of_tests = [Test1,Test2,Test3,Test4,Test5,Test6]

for i in list_of_tests:
    try:
        methyl_syntax(i)
        batch_queries = i.split('|')
    except:
        batch_queries = [i]

    print(batch_queries)

#test_split = Test6.split('|')
#print(test_split)