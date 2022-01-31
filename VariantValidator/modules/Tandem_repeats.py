"""

Script to check syntax of expanded repeat variants
By Rebecca Locke + Rob Wilson

"""
# Import modules
import re
import logging

#  List of variants to check format and split into constituents
variant1 = "LRG_199t1:c.1ACT[20]"
variant2 = "LRG_199:g.1ACT[20]A"
variant3 = "LRG_199:g.1AC[20]"
variant4 = "LRG_199t1:c.1_3ACT[20]"
variant5 = "LRG_199t1:c.1AC[10]"
variant6 = "LRG_199t1:c.1act[20]"
variant7 = "LRG_199t1:c.1A[12]"
variant8 = "LRG_199:g.13ACT[20]"
variant9 = "LRG_199:g.13_25ACTG[12]"
variant10 = "LRG199t3:c.13_125ACTG[5]"
variant11 = "ENSG00000198947.15:g.1ACT[10]"
variant12 = "ENST00000357033.8:c.13AC[2]"
variant13 = "LRG_199t1:c.1_2ACT[20]"
variant14 = "LRG_199t1:c.20A[A]"
variant15 = "NM_004006.2:c.13AC[22]"

variant17 = "NG_004006.2:g.1_2act[22]"
variant18 = "NM_024312.4:c.2686A[10]"
variant19 = "NM_024312.4:c.1738TA[6]"
variant20 = "LRG199t2:c.1_5C[10]"
variant21 = "LRG199t2:c.1-5C[10]"
variant22 = "NM_024312.1:c.2686A[10]"
# Alleles not supported
variant23 = "LRG_199:g.[123456A>G];[345678G>C]"

variant24 = "LRG_199t1:c.15_20AG[10]"
variant25 = "LRG_199:g.1AG[10]"

# Parse the variant to get relevant parts
def parse_repeat_variant(my_variant):
    """
    Summary:
    This takes a variant string and breaks it into its constituents with regex.

    Args:
        my_variant (string): Variant string e.g. "LRG_199:g.1ACT[20]A"
    Returns:
        prefix (string): Transcript or gene; everything before the first colon, e.g. "LRG_199"
        var_type (string): The variant genomic or coding type e.g. "g"
        var_pos (string): Position of the variant, e.g. "1" or "1_12"
        repeated_seq (string): The repeated sequence e.g. "ACT"
        no_of_repeats (string): The number of repeat units e.g. "20"
        after_the_bracket (string): Captures anything after the number of repeats bracket e.g. "A"
    """
    if "[" or "]" in my_variant:
        assert ":" in my_variant, f"Unable to identify a colon (:) in the variant description {my_variant}. A colon is required in HGVS variant descriptions to separate the reference accession from the reference type i.e. <accession>:<type>. e.g. :c"
        assert ";" not in my_variant, "Alleles not yet supported"
        assert "," not in my_variant, "Alleles not yet supported"
        prefix, suffix = my_variant.split(":")
        # Find reference sequence used (g or c)
        variant_type = re.search('^.*?(.*?)\.', suffix)
        var_type = variant_type.group(1)
        # Get g or c position(s)
        # Extract bit between . and [ e.g. 1ACT
        pos_and_seq = suffix.split(".")[1].split("[")[0]
        assert re.search(
            "[a-z]+", pos_and_seq, re.IGNORECASE), "Please ensure that the repeated sequence is included between the position and number of repeat units, e.g. g.1ACT[20]"
        rep_seq = re.search("[ACTG]+", pos_and_seq, re.IGNORECASE)
        repeated_seq = rep_seq.group()
        # Ensure sign used to indicate range is “_” (underscore), not “-“ (minus)
        if "-" in pos_and_seq:
            pos_and_seq = pos_and_seq.replace('-', '_')
        # Check both ends of range are given
        if "_" in pos_and_seq:
            assert re.search(
                "[0-9]+_[0-9]+", pos_and_seq), "Please ensure the start and the end of the full repeat range is provided, separated by an underscore"
            variant_positions = re.search("[0-9]+_[0-9]+", pos_and_seq)
            var_pos = variant_positions.group()
        else:
            # If just start pos, get digits
            variant_position = re.search("\d+", pos_and_seq)
            var_pos = variant_position.group()
        # Get number of unit repeats
        repeat_no = re.search('\[(.*?)\]', my_variant)
        no_of_repeats = repeat_no.group(1)
        # Get anything after ] to check
        if re.search('\](.*)', my_variant):
            after_brac = re.search('\](.*)', my_variant)
            after_the_bracket = after_brac.group(1)
        else:
            after_the_bracket = ""
        return prefix, var_type, var_pos, repeated_seq, no_of_repeats, after_the_bracket


variant_check = parse_repeat_variant(variant25)

the_prefix = variant_check[0]
variant_type = variant_check[1]
variant_position = variant_check[2]
repeated_sequence = variant_check[3]
number_of_repeats = variant_check[4]
after_bracket = variant_check[5]

def check_transcript_type(prefix):
    """
    Summary:
        Find transcript type. N.B. Future development could instead store the transcript and replace it with RefSeq.
    Args:
        prefix (string): The prefix from parse_variant_repeat
    Returns: 
        None, prints variant type
    Raises:
        NameError: [Error for unknown transcript type.]
    """
    if bool(re.match(r"^LRG", prefix)):
        print("LRG variant")
        # reformat_prefix_LRG(prefix)
    elif bool(re.match(r"^E", prefix)):
        print("Ensembl variant")
    elif bool(re.match(r"^N", prefix)):
        print("RefSeq variant")
    else:
        raise NameError('Unknown transcript type present. \
                        Try RefSeq transcript ID')


check_transcript_type(the_prefix)


def reformat_prefix(prefix):
    if re.match(r'^LRG', prefix):
        if re.match(r'^LRG\d+', prefix):
            prefix = prefix.replace('LRG', 'LRG_')
            print("LRG variant updated to include underscore")
        # Get transcript number
        if "t" in prefix:
            transcript_num = re.search("t(.*?)$", prefix)
            transcript_version = f"t{transcript_num.group(1)}"
    elif re.match(r'^ENST', prefix) or re.match(r'^NM_', prefix):
        assert "." in prefix, "Please ensure the transcript version is included following a '.' after the transcript name e.g. ENST00000357033.8"
    return prefix


the_prefix = reformat_prefix(the_prefix)


def check_genomic_or_coding(prefix, var_type):
    """Takes prefix and works out if variant type should be c. or g. and raises error if incorrect type supplied
    Args:
        prefix (string): The prefix e.g. "LRG_199"
        var_type (string): Variant type genomic or coding e.g. "g"
    """
    if re.match(r'^LRG', prefix):
        if "t" in prefix:
            assert var_type == "c", "Please ensure variant type is coding if an LRG transcript is provided"
        else:
            assert var_type == "g", "Please ensure variant type is genomic if LRG gene is used"
    elif re.match(r'^ENST', prefix):
        assert var_type == "c", "Please ensure variant type is coding if an Ensembl transcript is provided"
    elif re.match(r'^ENSG', prefix):
        assert var_type == "g", "Please ensure variant type is genomic if Ensembl gene is used"
    elif re.match(r'^NM', prefix):
        assert var_type == "c", "Please ensure variant type is coding if a RefSeq transcript is provided"
    elif re.match(r'^NC', prefix):
        assert var_type == "g", "Please ensure variant type is genomic if RefSeq chromosome is used"
    elif re.match(r'^NG', prefix):
        assert var_type == "g", "Please ensure variant type is genomic if RefSeq gene is used"


check_genomic_or_coding(the_prefix, variant_type)

# For variants with the full range of the position given (not only start pos)
def check_positions_given(repeated_sequence, variant_pos, no_of_rep_units):
    """Checks the position range given and updates it if it doesn't match the length of the repeated sequence and number of repeat units when full range is needed
        Args:
        repeated_sequence (string): The repeated sequence e.g. "ACT"
        variant_pos (string): The position of the variant e.g. "1" or "1_5"
        no_of_rep_units (string): The number of repeat units e.g. "20"
    Returns: 
        full_range (string): The full range supplied if correct or the full range updated if inputted range was incorrect, e.g. "1_20"
    """
    start_range, end_range = variant_pos.split("_")
    rep_seq_length = len(repeated_sequence)
    the_range = int(end_range) - int(start_range) + 1
    repeat_length = (rep_seq_length * int(no_of_rep_units))
    if the_range == repeat_length:
        print("Range given matches repeat sequence length and number of repeat units")
        full_range = f"{start_range}_{end_range}"
    else:
        print("Warning: sequence range (X_X) given must match repeat unit sequence length and number of repeat units. Updating the range based on repeat sequence length and number of repeat units")
        new_end_range = int(start_range) + repeat_length - 1
        full_range = f"{start_range}_{new_end_range}"
    return full_range

# # Note Community Consultation is prepared which will suggest to allow only one format where the entire range of the repeated sequence must be indicated, e.g. g.123_191CAG[23]. This small function will give you the range from getting only a start position


def get_range_from_single_pos(repeated_sequence, start_range, no_of_rep_units):
    rep_seq_length = len(repeated_sequence)
    repeat_length = (rep_seq_length * int(no_of_rep_units))
    the_end_range = int(start_range) + repeat_length - 1
    full_range = f"{start_range}_{the_end_range}"
    return full_range

"""exception: using a coding DNA reference sequence (“c.” description) a Repeated sequence variant description can be used only for repeat units with a length which is a multiple of 3, i.e. which can not affect the reading frame. Consequently, use NM_024312.4:c.2692_2693dup and not NM_024312.4:c.2686A[10], use NM_024312.4:c.1741_1742insTATATATA and not NM_024312.4:c.1738TA[6]."""

# This will reformat tandem repeat variants in c. which should be noted as dup or ins as they are not multiples of 3
def reformat_not_multiple_of_three(pref, vartype, position, rep_seq, no_of_repeats):
    reformatted = ""
    rep_seq_length = len(rep_seq)
    # Repeat of 1 base should be a dup with full range given
    if rep_seq_length == 1:
        if "_" in position:
            position = check_positions_given(rep_seq, position, no_of_repeats)
        else:
            position = get_range_from_single_pos(rep_seq, position, no_of_repeats)
        print("Warning: Repeated sequence is not a multiple of three! Updating variant description to a duplication")
        reformatted = f'{pref}:{vartype}.{position}dup'
    # Repeat of 2 bases should be an ins with only first two nts given as range
    elif rep_seq_length >= 2:
        expanded_rep_seq = rep_seq*int(no_of_repeats)
        if not "_" in position:
            second_range = int(position)+1
            position = f"{position}_{second_range}"
        else:
            start, end = position.split("_")
            end = int(start)+1
            position = f"{start}_{end}"
        print("Warning: Repeated sequence is not a multiple of three! Updating variant description to insertion")
        reformatted = f'{pref}:{vartype}.{position}ins{expanded_rep_seq}'
    return reformatted

# Reformat the variant for HGVS consistency

def reformat(var_prefix, the_variant_type, the_var_pos, the_repeated_sequence, the_number_of_repeats, all_after_bracket):
    # Check number of repeat units is integers and that the sequence is A,C,T or G
    assert the_number_of_repeats.isdecimal(
    ), "The number of repeat units included between square brackets must be numeric"
    assert re.search("[actg]+", the_repeated_sequence,
                     re.IGNORECASE), "Please ensure the repeated sequence includes only A, C, T or G"
    # Update the repeated sequence to be upper case
    the_repeated_sequence = the_repeated_sequence.upper()
    if all_after_bracket != "":
        print("No information should be included after the number of repeat units. Mixed repeats are not currently supported.")
    # Reformat c. variants
    if the_variant_type == "c": 
        rep_seq_length = len(the_repeated_sequence)
        if rep_seq_length % 3 != 0:
            final_format = reformat_not_multiple_of_three(var_prefix, the_variant_type, the_var_pos, the_repeated_sequence, the_number_of_repeats)
        else:
            print("Repeat length is consistent with c. type")
            if "_" in the_var_pos:
                the_var_pos = check_positions_given(the_repeated_sequence, the_var_pos, the_number_of_repeats)
            final_format = f"{var_prefix}:{the_variant_type}.{the_var_pos}{the_repeated_sequence}[{the_number_of_repeats}]"
    # Reformat g. variants
    else:
        if "_" in the_var_pos:
            the_var_pos = check_positions_given(the_repeated_sequence, the_var_pos, the_number_of_repeats)
        final_format = f"{var_prefix}:{the_variant_type}.{the_var_pos}{the_repeated_sequence}[{the_number_of_repeats}]"
    return final_format

print(reformat(the_prefix, variant_type, variant_position,
      repeated_sequence, number_of_repeats, after_bracket))


# def check_no_repeats(start_range, end_range, repeated_sequence):
#     rep_seq_length = len(repeated_sequence)
#     the_range = int(end_range) - int(start_range) + 1
#     no_of_units = int(the_range / rep_seq_length)
#     return no_of_units


# if __name__ == "__main__":
#    main()
