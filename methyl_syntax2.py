# Pipe syntax error

# batch_variant = the string containing the variant to be validated


def methyl_syntax(my_variant):
    """
    :param my_variant:
    :return: False if no | is detected otherwise raises Exception and provides information as to where the
    | is located
    @type my_variant: object
    """

    if "|" in my_variant:
        print(my_variant)
        if my_variant.endswith("gom"):
            raise Exception()
        elif "gom" in my_variant and not my_variant.endswith("gom"):
            print("Variant is in incorrect format")
        elif my_variant.endswith("lom"):
            raise Exception()
        elif "lom" in my_variant and not my_variant.endswith("lom"):
            print("Variant is in incorrect format")
        elif my_variant.endswith("met="):
            raise Exception()
        elif "met=" in my_variant and not my_variant.endswith("met="):
            print("Variant is in incorrect format")
        else:
            return True


Test1 = "NC_000011.10::g.1999904_1999946|gom"
Test2 = "NC_000011.10::g.1999904_1999946|lom"
Test3 = "NC_000011.10::g.1999904_1999946|met="
Test4 = "NM_000719.7:c.5550G>A"
Test5 = "NM_00256.3:c.2373dupG"
Test6 = "NM_00256.3:c.2373dupG | NM_000719.7:c.5550G>A"
Test7 = "NC_000011.10|gom:g.1999904_1999946"

list_of_tests = [Test1,Test2,Test3,Test4,Test5,Test6,Test7]

for i in list_of_tests:
    try:
        methyl_syntax(i)
        batch_queries = i.split('|')
    except:
        batch_queries = [i]

    print(batch_queries)
