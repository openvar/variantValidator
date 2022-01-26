import re


#Match variants that are correctly formatted expanded to deal with variants
def check_variants(input: str):
    """[Match variant to defined format to decide if to process or reformat,
    this is now deprecated. check_expanded_repeat() replaces this.]

    Args:
        input (str): [description]
    """
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


#  List of variants to check format and split into constituents.
variant1="LRG_199:g.1ACT[20]"
variant2 = "LRG_199:g.1ACT[20]A"
variant3 = "LRG_199:g.1AC[20]"
variant4 = "LRG_199t1:c.1_3ACT[20]"
variant5 = "LRG_199t1:c.1AC[20]"
variant6 = "LRG_199t1:c.1act[20]"
variant7 = "LRG199c.1A[1_2]"
variant8 = "LRG_199:g.13ACT[20]"
variant9 = "LRG_199:g.13_25ACTG[1]"
variant10 = "LRG199:g.13_125ACTG[1]"
# Other types not LRG
variant11 = "ENSG00000198947.15:g.1ACT[10]"
variant12 = "ENST00000357033.8:c.13AC[22]"
#Missing information accepted
variant13 = "LRG_199t1:c.1_ACT[20]"
# change * to + to only allow variants with range or single location

def check_expanded_repeat(my_variant):
    """
    Summary:
    This takes a variant string and breaks it into its constituents.
    This isolates the constituents with regex.

    Args:
        my_variant ([type]): (Variant string i.e. LRG_199:g.1ACT[20])
    Returns:
    Prints out constituents and assigns them to variables for further processing.
    """
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
                print(start_range)
                print(end_range)
                rep_seq = re.search('\.[0-9]+_[0-9]+(.*?)\[', my_variant)
                print(f'Repeated sequence: {rep_seq.group(1)}')
            else:
                rep_seq = re.search('\.[0-9]+(.*?)\[', my_variant)
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
