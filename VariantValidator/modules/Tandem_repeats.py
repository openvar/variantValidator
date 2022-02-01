"""

Script to check syntax of expanded repeat variants
By Rebecca Locke + Rob Wilson

"""
# Import modules
import re
import logging

class Tandem_repeats:
    def __init__(self, prefix, variant_type, variant_position, repeat_sequence,copy_number, after_the_bracket):
        self.prefix = prefix
        self.variant_type = variant_type
        self.variant_position = variant_position
        self.repeat_sequence = repeat_sequence
        self.copy_number = copy_number
        self.after_the_bracket = after_the_bracket

    @classmethod
    def parse_repeat_variant(cls, variant_str):
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
        if "[" or "]" in variant_str:
            assert ":" in variant_str, f"Unable to identify a colon (:) in the variant description {variant_str}. A colon is required in HGVS variant descriptions to separate the reference accession from the reference type i.e. <accession>:<type>. e.g. :c"
            assert ";" not in variant_str, "Alleles not yet supported"
            assert "," not in variant_str, "Alleles not yet supported"
            prefix, suffix = variant_str.split(":")
            # Find reference sequence used (g or c)
            var_type = re.search('^.*?(.*?)\.', suffix)
            variant_type = var_type.group(1)
            # Get g or c position(s)
            # Extract bit between . and [ e.g. 1ACT
            pos_and_seq = suffix.split(".")[1].split("[")[0]
            assert re.search(
                "[a-z]+", pos_and_seq, re.IGNORECASE), "Please ensure that the repeated sequence is included between the position and number of repeat units, e.g. g.1ACT[20]"
            rep_seq = re.search("[ACTG]+", pos_and_seq, re.IGNORECASE)
            repeat_sequence = rep_seq.group()
            # Ensure sign used to indicate range is “_” (underscore), not “-“ (minus)
            if "-" in pos_and_seq:
                pos_and_seq = pos_and_seq.replace('-', '_')
            # Check both ends of range are given
            if "_" in pos_and_seq:
                assert re.search(
                    "[0-9]+_[0-9]+", pos_and_seq), "Please ensure the start and the end of the full repeat range is provided, separated by an underscore"
                variant_positions = re.search("[0-9]+_[0-9]+", pos_and_seq)
                variant_position = variant_positions.group()
            else:
                # If just start pos, get digits
                variant_position = re.search("\d+", pos_and_seq)
                variant_position = variant_position.group()
            # Get number of unit repeats
            repeat_no = re.search('\[(.*?)\]', variant_str)
            copy_number = repeat_no.group(1)
            # Get anything after ] to check
            if re.search('\](.*)', variant_str):
                after_brac = re.search('\](.*)', variant_str)
                after_the_bracket = after_brac.group(1)
            else:
                after_the_bracket = ""
        return cls(prefix, variant_type, variant_position, repeat_sequence, copy_number, after_the_bracket)

    def check_transcript_type(self):
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
        if bool(re.match(r"^LRG", self.prefix)):
            print("LRG variant")
            # reformat_prefix_LRG(prefix)
        elif bool(re.match(r"^E", self.prefix)):
            print("Ensembl variant")
        elif bool(re.match(r"^N", self.prefix)):
            print("RefSeq variant")
        else:
            raise NameError('Unknown transcript type present. \
                            Try RefSeq transcript ID')

    def reformat_prefix(self):
        if re.match(r'^LRG', self.prefix):
            if re.match(r'^LRG\d+', self.prefix):
                self.prefix = self.prefix.replace('LRG', 'LRG_')
                print("LRG variant updated to include underscore")
            # Get transcript number
            if "t" in self.prefix:
                transcript_num = re.search("t(.*?)$", self.prefix)
                transcript_version = f"t{transcript_num.group(1)}"
        elif re.match(r'^ENS', self.prefix) or re.match(r'^N', self.prefix):
            assert "." in self.prefix, "Please ensure the transcript or gene version is included following a '.' after the transcript or gene name e.g. ENST00000357033.8"
        return self.prefix

    def check_genomic_or_coding(self):
        """Takes prefix and works out if variant type should be c. or g. and raises error if incorrect type supplied
        Args:
            prefix (string): The prefix e.g. "LRG_199"
            var_type (string): Variant type genomic or coding e.g. "g"
        Returns:
            None, gives error if wrong variant type is used
        """
        if re.match(r'^LRG', self.prefix):
            if "t" in self.prefix:
                assert self.variant_type == "c", "Please ensure variant type is coding if an LRG transcript is provided"
            else:
                assert self.variant_type == "g", "Please ensure variant type is genomic if LRG gene is used"
        elif re.match(r'^ENST', self.prefix):
            assert self.variant_type == "c", "Please ensure variant type is coding if an Ensembl transcript is provided"
        elif re.match(r'^ENSG', self.prefix):
            assert self.variant_type == "g", "Please ensure variant type is genomic if Ensembl gene is used"
        elif re.match(r'^NM', self.prefix):
            assert self.variant_type == "c", "Please ensure variant type is coding if a RefSeq transcript is provided"
        elif re.match(r'^NC', self.prefix):
            assert self.variant_type == "g", "Please ensure variant type is genomic if RefSeq chromosome is used"
        elif re.match(r'^NG', self.prefix):
            assert self.variant_type == "g", "Please ensure variant type is genomic if RefSeq gene is used"

    def check_positions_given(self):
        """Checks the position range given and updates it if it doesn't match the length of the repeated sequence and number of repeat units when full range is needed
            Args:
            repeated_sequence (string): The repeated sequence e.g. "ACT"
            variant_pos (string): The position of the variant e.g. "1" or "1_5"
            no_of_rep_units (string): The number of repeat units e.g. "20"
        Returns: 
            full_range (string): The full range supplied if correct or the full range updated if inputted range was incorrect, e.g. "1_20"
        """
        start_range, end_range = self.variant_position.split("_")
        rep_seq_length = len(self.repeat_sequence)
        the_range = int(end_range) - int(start_range) + 1
        repeat_length = (rep_seq_length * int(self.copy_number))
        if the_range == repeat_length:
            print("Range given matches repeat sequence length and number of repeat units")
            full_range = f"{start_range}_{end_range}"
        else:
            print("Warning: sequence range (X_X) given must match repeat unit sequence length and number of repeat units. Updating the range based on repeat sequence length and number of repeat units")
            new_end_range = int(start_range) + repeat_length - 1
            full_range = f"{start_range}_{new_end_range}"
        return full_range

    def get_range_from_single_pos(self):
        rep_seq_length = len(self.repeat_sequence)
        repeat_range = (rep_seq_length * int(self.copy_number))
        the_end_range = int(self.variant_position) + repeat_range - 1
        full_range = f"{self.variant_position}_{the_end_range}"
        return full_range

    """exception: using a coding DNA reference sequence (“c.” description) a Repeated sequence variant description can be used only for repeat units with a length which is a multiple of 3, i.e. which can not affect the reading frame. Consequently, use NM_024312.4:c.2692_2693dup and not NM_024312.4:c.2686A[10], use NM_024312.4:c.1741_1742insTATATATA and not NM_024312.4:c.1738TA[6]."""
    # This will reformat tandem repeat variants in c. which should be noted as dup or ins as they are not multiples of 3
    def reformat_not_multiple_of_three(self):
        reformatted = ""
        rep_seq_length = len(self.repeat_sequence)
        # Repeat of 1 base should be a dup with full range given
        if rep_seq_length == 1:
            if "_" in self.variant_position:
                self.variant_position = self.check_positions_given()
            else:
                self.variant_position = self.get_range_from_single_pos()
            print("Warning: Repeated sequence is coding and not a multiple of three! Updating variant description to a duplication")
            reformatted = f'{self.prefix}:{self.variant_type}.{self.variant_position}dup'
        # Repeat of 2 bases should be an ins with only first two nts given as range
        elif rep_seq_length >= 2:
            expanded_rep_seq = self.repeat_sequence*int(self.copy_number)
            if not "_" in self.variant_position:
                second_range = int(self.variant_position)+1
                position = f"{self.variant_position}_{second_range}"
            else:
                start, end = self.variant_position.split("_")
                end = int(start)+1
                position = f"{start}_{end}"
            print("Warning: Repeated sequence is coding and not a multiple of three! Updating variant description to insertion")
            reformatted = f'{self.prefix}:{self.variant_type}.{position}ins{expanded_rep_seq}'
        return reformatted

    def reformat(self):
        # Check number of repeat units is integers and that the sequence is A,C,T or G
        assert self.copy_number.isdecimal(
        ), "The number of repeat units included between square brackets must be numeric"
        assert re.search("[actg]+", self.repeat_sequence,
                        re.IGNORECASE), "Please ensure the repeated sequence includes only A, C, T or G"
        # Update the repeated sequence to be upper case
        self.repeat_sequence = self.repeat_sequence.upper()
        if self.after_the_bracket != "":
            print("No information should be included after the number of repeat units. Mixed repeats are not currently supported.")
        # Reformat c. variants
        if self.variant_type == "c": 
            rep_seq_length = len(self.repeat_sequence)
            if rep_seq_length % 3 != 0:
                final_format = self.reformat_not_multiple_of_three()
            else:
                print("Repeat length is consistent with c. type")
                if "_" in self.variant_position:
                    self.variant_position = self.check_positions_given()
                final_format = f"{self.prefix}:{self.variant_type}.{self.variant_position}{self.repeat_sequence}[{self.copy_number}]"
        # Reformat g. variants
        else:
            if "_" in self.variant_position:
                self.variant_position = self.check_positions_given()
            final_format = f"{self.prefix}:{self.variant_type}.{self.variant_position}{self.repeat_sequence}[{self.copy_number}]"
        return final_format

# Gives LRG_199t1:c.1_2insACACACACACACACACACACACACACAC
variant1 = "LRG_199t1:c.1_5AC[14]"
# Gives LRG_199:g.1ACT[20]
variant2 = "LRG_199:g.1ACT[20]A"
# Gives LRG_199:g.1AC[20]
variant3 = "LRG_199:g.1AC[20]"
# Gives LRG_199t1:c.1_60ACT[20]
variant4 = "LRG_199t1:c.1_3ACT[20]"
# Gives LRG_199t1:c.1_2insACACACACACACACACACAC
variant5 = "LRG_199t1:c.1AC[10]"
# Gives LRG_199t1:c.1ACT[20]
variant6 = "LRG_199t1:c.1act[20]"
# Gives LRG_199t1:c.1_12dup
variant7 = "LRG_199t1:c.1A[12]"
# Gives LRG_199:g.13ACT[20]
variant8 = "LRG_199:g.13ACT[20]"
# Gives LRG_199:g.13_60ACTG[12]
variant9 = "LRG_199:g.13_25ACTG[12]"
# Gives LRG_199t3:c.13_14insACTGACTGACTGACTGACTG
variant10 = "LRG199t3:c.13_125ACTG[5]"
# Gives ENSG00000198947.15:g.1ACT[10]
variant11 = "ENSG00000198947.15:g.1ACT[10]"
# Gives ENST00000357033.8:c.13_14insACAC
variant12 = "ENST00000357033.8:c.13AC[2]"
# Gives LRG_199t1:c.1_60ACT[20]
variant13 = "LRG_199t1:c.1_2ACT[20]"
# Gives AssertionError: The number of repeat units included between square brackets must be numeric
variant14 = "LRG_199t1:c.20A[A]"
# Gives NM_004006.2:c.13_14insACACACACACACAC
variant15 = "NM_004006.2:c.13AC[7]"
# Gives AssertionError: Unable to identify a colon (:) in the variant description NG_004006.2g.1_2act[22]. A colon is required in HGVS variant descriptions to separate the reference accession from the reference type i.e. <accession>:<type>. e.g. :c
variant16 = "NG_004006.2g.1_2act[22]"
# Gives NG_004006.2:g.1_66ACT[22]
variant17 = "NG_004006.2:g.1_2act[22]"
# Gives NM_024312.4:c.2686_2695dup
variant18 = "NM_024312.4:c.2686A[10]"
# Gives NM_024312.4:c.1738_1739insTATATATATATA
variant19 = "NM_024312.4:c.1738TA[6]"
# Gives LRG_199t2:c.1_10dup
variant20 = "LRG199t2:c.1_5C[10]"
# Gives LRG_199t2:c.1_10dup
variant21 = "LRG199t2:c.1-5C[10]"
# Gives NM_024312.1:c.2686_2695dup
variant22 = "NM_024312.1:c.2686A[10]"
# Gives Error: Alleles not supported
variant23 = "LRG_199:g.[123456A>G];[345678G>C]"
# Gives LRG_199t1:c.15_16insAGAGAGAGAGAGAGAGAGAG
variant24 = "LRG_199t1:c.15_20AG[10]"
# Gives AssertionError: Please ensure the transcript or gene version is included following a '.' after the transcript or gene name e.g. ENST00000357033.8
variant25 = "ENST00000198947:c.1_2AG[10]"


my_variant = Tandem_repeats.parse_repeat_variant(variant1)
Tandem_repeats.check_transcript_type(my_variant)
Tandem_repeats.reformat_prefix(my_variant)
Tandem_repeats.check_genomic_or_coding(my_variant)
print(Tandem_repeats.reformat(my_variant))