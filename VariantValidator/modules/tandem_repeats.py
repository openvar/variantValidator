"""
NAME:          modules/tandem_repeats.py
AUTHORS:       Rebecca Locke (@rklocke) & Robert Wilson (@RSWilson1)
DATE:          18/02/22
INSTITUTION:   University of Manchester/Cambridge University Hospitals
DESCRIPTION:   This script contains the Tandem Repeats class and methods, aiming to ultimately check the syntax and reformat tandem repeat variants for VariantValidator
"""

# Import modules
import json
import re
import logging
import os
import VariantValidator

# Initialise vv class
vval = VariantValidator.Validator()

# Get path of directory the script is run in
CURRENT_DIR = os.path.abspath(os.getcwd())

# Configure logging message format
LOG_FORMAT = "%(asctime)s — %(name)s — %(levelname)s — %(lineno)d — %(message)s"

# Set logging level to debug, define logging format, re-write file each time
logging.basicConfig(
    filename=f"{CURRENT_DIR}/expanded_repeats.log",
    filemode="w",
    level=logging.DEBUG,
    format=LOG_FORMAT,
)

# Set up logger
logger = logging.getLogger()


class TandemRepeats:

    """Class to represent a tandem repeat variant

    Attributes
    ----------
    reference : str
        reference sequence name
    prefix : str
        reference sequence type used
    variant_position : str
        nucleotide position(s)
    repeat_sequence : str
        sequence repeat unit
    copy_number : str
        number of repeat units
    after_the_bracket: str
        anything after the last square bracket
    """

    def __init__(
        self,
        reference,
        prefix,
        variant_position,
        repeat_sequence,
        copy_number,
        after_the_bracket,
        build,
        select_transcripts,
    ):
        """Constructs necessary parts of the TandemRepeats object"""
        self.reference = reference
        self.prefix = prefix
        self.variant_position = variant_position
        self.repeat_sequence = repeat_sequence
        self.copy_number = copy_number
        self.after_the_bracket = after_the_bracket
        self.build = build
        self.select_transcripts = select_transcripts

    @classmethod
    def parse_repeat_variant(cls, variant_str, build, select_transcripts):
        """
        Summary:
            Takes a variant string and breaks it into its constituent parts with regex to be processed in downstream functions, assigns them to class variables.
        Args:
            variant_str : str
                variant string e.g. "LRG_199:g.1ACT[20]A"
        Returns:
            Updates each class attribute:

            reference : str
                Transcript or gene; everything before the first colon, e.g. "LRG_199"
            prefix : str
                The variant genomic or coding type e.g. "g"
            variant_position : str
                Position of the variant, e.g. "1" or "1_12"
            repeat_sequence : str
                The repeated sequence e.g. "ACT"
            copy_number : str
                The number of repeat units e.g. "20"
            after_the_bracket : str
                Captures anything after the number of repeats bracket e.g. "A"
        Example:
            >>>parse_repeat_variant("LRG_199:g.1ACT[20]A")
                "LRG_199", "g", "1", "ACT", "20", "A"
        """
        print(f"Variant entered: {variant_str}")
        logger.info(f"Parsing variant: parse_repeat_variant({variant_str})")
        # Strip any whitespace
        variant_str = variant_str.strip()

        # Check if square brackets included which indicate tandem repeat
        # variant
        if "[" or "]" in variant_str:
            try:
                assert ":" in variant_str, f"Unable to identify a colon (:) in the variant description {variant_str}. A colon is required in HGVS variant descriptions to separate the reference accession from the reference type i.e. <accession>:<type>. e.g. :c"
            except AssertionError:
                logger.critical("A colon is required in the variant description. Ending program")
                raise
            else:
                reference, suffix = variant_str.split(":")
            assert ";" not in variant_str, "Alleles not yet supported"
            assert "," not in variant_str, "Alleles not yet supported"
            # Find reference sequence used (g or c)
            var_type = re.search("^.*?(.*?)\\.", suffix)
            prefix = var_type.group(1).lower()
            # Get g or c position(s)
            # Extract bit between . and [ e.g. 1ACT
            pos_and_seq = suffix.split(".")[1].split("[")[0]
            try:
                assert re.search(
                "[a-z]+", pos_and_seq, re.IGNORECASE), \
                "Please ensure that the repeated sequence is included between the position and number of repeat units, e.g. g.1ACT[20]"
            except AssertionError:
                logger.critical("Unable to identify a repeated sequence between position and number of repeats. Ending program")
                raise
            else:
                rep_seq = re.search("[ACTG]+", pos_and_seq, re.IGNORECASE)
                repeat_sequence = rep_seq.group()
            # Ensure sign used to indicate range is “_” (underscore), not “-“
            if "-" in pos_and_seq:
                pos_and_seq = pos_and_seq.replace("-", "_")
            # Check both ends of range are given
            if "_" in pos_and_seq:
                try:
                    assert re.search(
                            "[0-9]+_[0-9]+", pos_and_seq
                        ), "Please ensure the start and the end of the full repeat range is provided, separated by an underscore"
                except AssertionError:
                    logger.critical("Only one value in the range is provided. Ending program")
                    raise
                else:
                    variant_positions = re.search("[0-9]+_[0-9]+", pos_and_seq)
                    variant_position = variant_positions.group()
            else:
                # If just start pos, get digits
                variant_position = re.search("\\d+", pos_and_seq)
                variant_position = variant_position.group()
            # Get number of unit repeats
            repeat_no = re.search("\\[(.*?)\\]", variant_str)
            copy_number = repeat_no.group(1)
            # Get anything after ] to check
            if re.search("\\](.*)", variant_str):
                after_brac = re.search("\\](.*)", variant_str)
                after_the_bracket = after_brac.group(1)
            else:
                after_the_bracket = ""
        return cls(
            reference,
            prefix,
            variant_position,
            repeat_sequence,
            copy_number,
            after_the_bracket,
            build,
            select_transcripts,
        )

    def check_transcript_type(self):
        """
        Summary:
            Find transcript type. N.B. Future development could instead store the transcript and replace it with RefSeq.
        Args:
            reference (string): The reference from parse_variant_repeat
        Returns:
            None, prints variant type
        Raises:
            NameError: [Error for unknown transcript type.]
        """
        logger.info(
            f"Checking transcript type: check_transcript_type({self.reference})"
        )
        if bool(re.match(r"^LRG", self.reference)):
            logger.info("Variant type: LRG variant")
        elif bool(re.match(r"^E", self.reference)):
            logger.info("Variant type: Ensembl variant")
        elif bool(re.match(r"^N", self.reference)):
            logger.info("Variant type: RefSeq variant")
        else:
            raise NameError(
                "Unknown transcript type present. \
                            Try RefSeq transcript ID"
            )

    def reformat_reference(self):
        """Reformats the reference sequence name"""
        logger.info(f"Reformatting reference: reformat_reference({self.reference})")
        if re.match(r"^LRG", self.reference):
            if re.match(r"^LRG\d+", self.reference):
                self.reference = self.reference.replace("LRG", "LRG_")
                logger.warning("LRG variant updated to include underscore")
            # Get transcript number
            if "t" in self.reference:
                transcript_num = re.search("t(.*?)$", self.reference)
                transcript_version = f"t{transcript_num.group(1)}"
        elif re.match(r"^ENS", self.reference) or re.match(r"^N", self.reference):
            assert (
                "." in self.reference
            ), "Please ensure the transcript or gene version is included following a '.' after the transcript or gene name e.g. ENST00000357033.8"
        return self.reference

    def check_genomic_or_coding(self):
        """Takes reference and works out what prefix type should be used
        Raises:
            AssertionError if wrong prefix type is used
        """
        logger.info(
            f"Checking prefix is consistent with reference: check_genomic_or_coding({self.reference},{self.prefix})"
        )
        if re.match(r"^LRG", self.reference):
            if "t" in self.reference:
                assert (
                    self.prefix == "c"
                ), "Please ensure variant type is coding if an LRG transcript is provided"
            else:
                assert (
                    self.prefix == "g"
                ), "Please ensure variant type is genomic if LRG gene is used"
        elif re.match(r"^ENST", self.reference):
            assert (
                self.prefix == "c"
            ), "Please ensure variant type is coding if an Ensembl transcript is provided"
        elif re.match(r"^ENSG", self.reference):
            assert (
                self.prefix == "g"
            ), "Please ensure variant type is genomic if Ensembl gene is used"
        elif re.match(r"^NM", self.reference):
            assert (
                self.prefix == "c"
            ), "Please ensure variant type is coding if a RefSeq transcript is provided"
        elif re.match(r"^NC", self.reference):
            assert (
                self.prefix == "g"
            ), "Please ensure variant type is genomic if RefSeq chromosome is used"
        elif re.match(r"^NG", self.reference):
            assert (
                self.prefix == "g"
            ), "Please ensure variant type is genomic if RefSeq gene is used"
        elif re.match(r"^NR", self.reference):
            assert (
                self.prefix == "n"
            ), "Please ensure variant type is non-coding if NR transcript is used"

    def check_positions_given(self):
        """Checks the position range given and updates it if it doesn't match the length of the repeated sequence and number of repeat units when full range is needed
        Args:
            repeat_sequence (string): The repeated sequence e.g. "ACT"
            variant_position (string): The position of the variant e.g. "1" or "1_5"
            copy_number (string): The number of repeat units e.g. "20"
        Returns:
            full_range (string): The full range supplied if correct or the full range updated if inputted range was incorrect, e.g. "1_20"
        """
        logger.info(
            f"Checking range given: check_positions_given({self.repeat_sequence},{self.variant_position},{self.copy_number})"
        )
        start_range, end_range = self.variant_position.split("_")
        rep_seq_length = len(self.repeat_sequence)
        the_range = int(end_range) - int(start_range) + 1
        repeat_length = rep_seq_length * int(self.copy_number)
        if the_range == repeat_length:
            logger.info(
                f"Range given ({self.variant_position}) matches repeat sequence length and number of repeat units"
            )
            full_range = f"{start_range}_{end_range}"
        else:
            logger.warning(
                f"Warning: sequence range {self.variant_position} given does not match repeat unit sequence length and number of repeat units. \
                Updating the range based on repeat sequence length and number of repeat units"
            )
            new_end_range = int(start_range) + repeat_length - 1
            full_range = f"{start_range}_{new_end_range}"
        return full_range

    def get_range_from_single_pos(self):
        """
        Gets full range of the variant if this is needed when a single start position is supplied
        """
        logger.info(
            f"Fetching the range from a given single position: get_range_from_single_pos({self.repeat_sequence},{self.copy_number},{self.variant_position}"
        )
        rep_seq_length = len(self.repeat_sequence)
        repeat_range = rep_seq_length * int(self.copy_number)
        the_end_range = int(self.variant_position) + repeat_range - 1
        full_range = f"{self.variant_position}_{the_end_range}"
        return full_range

    def reformat_not_multiple_of_three(self):
        """
        Reformats coding variants (c.) to a dup or ins if they are not a multiple of three
        """
        logger.info(
            f"Reformatting variant as not a multiple of three: reformat_not_multiple_of_three({self.repeat_sequence},{self.variant_position},{self.reference},{self.prefix})"
        )
        reformatted = ""
        rep_seq_length = len(self.repeat_sequence)
        # Repeat of 1 base should be a dup with full range given
        if rep_seq_length == 1:
            if "_" in self.variant_position:
                self.variant_position = self.check_positions_given()
            else:
                self.variant_position = self.get_range_from_single_pos()
            logger.warning(
                f"Warning: Repeated sequence is coding and is of length {rep_seq_length}, not a multiple of three! Updating variant description to a duplication"
            )
            reformatted = f"{self.reference}:{self.prefix}.{self.variant_position}dup"
        # Repeat of 2 bases should be an ins with only first two nts given as
        # range
        elif rep_seq_length >= 2:
            expanded_rep_seq = self.repeat_sequence * int(self.copy_number)
            if "_" not in self.variant_position:
                second_range = int(self.variant_position) + 1
                position = f"{self.variant_position}_{second_range}"
            else:
                start, end = self.variant_position.split("_")
                end = int(start) + 1
                position = f"{start}_{end}"
            logger.warning(
                f"Warning: Repeated sequence is coding and is of length {rep_seq_length}, not a multiple of three! Updating variant description to insertion"
            )
            reformatted = (
                f"{self.reference}:{self.prefix}.{position}ins{expanded_rep_seq}"
            )
        return reformatted

    def reformat(self):
        """Reformats and returns final formatted variant as a string"""
        logger.info(
            f"Reformatting variant: reformat({self.repeat_sequence},{self.after_the_bracket},{self.prefix},{self.variant_position},{self.copy_number})"
        )
        assert (
            self.copy_number.isdecimal()
        ), "The number of repeat units included between square brackets must be numeric"
        assert re.search(
            "[actg]+", self.repeat_sequence, re.IGNORECASE
        ), "Please ensure the repeated sequence includes only A, C, T or G"
        # Update the repeated sequence to be upper case
        self.repeat_sequence = self.repeat_sequence.upper()
        if self.after_the_bracket != "":
            logger.warning(
                f"No information should be included after the number of repeat units. Currently '{self.after_the_bracket}'' is included. This will be removed as mixed repeats are not currently supported."
            )
        # Reformat c. variants
        if self.prefix == "c":
            rep_seq_length = len(self.repeat_sequence)
            if rep_seq_length % 3 != 0:
                final_format = self.reformat_not_multiple_of_three()
            else:
                logger.info("Checked repeat length is consistent with c. type")
                if "_" in self.variant_position:
                    self.variant_position = self.check_positions_given()
                final_format = f"{self.reference}:{self.prefix}.{self.variant_position}{self.repeat_sequence}[{self.copy_number}]"
        # Reformat g. variants
        else:
            if "_" in self.variant_position:
                self.variant_position = self.check_positions_given()
            final_format = f"{self.reference}:{self.prefix}.{self.variant_position}{self.repeat_sequence}[{self.copy_number}]"
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
# Gives AssertionError: The number of repeat units included between square
# brackets must be numeric
variant14 = "LRG_199t1:c.20A[A]"
# Gives NM_004006.2:c.13_14insACACACACACACAC
variant15 = "NM_004006.2:c.13AC[7]"
# Gives AssertionError: Unable to identify a colon (:) in the variant
# description NG_004006.2g.1_2act[22]. A colon is required in HGVS variant
# descriptions to separate the reference accession from the reference type
# i.e. <accession>:<type>. e.g. :c
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
# Gives AssertionError: Please ensure the transcript or gene version is
# included following a '.' after the transcript or gene name e.g.
# ENST00000357033.8
variant25 = "ENST00000198947:c.1_2AG[10]"
# Gives AssertionError: Please ensure that the repeated sequence is
# included between the position and number of repeat units, e.g.
# g.1ACT[20]
variant26 = "ENST00000198947.1:c.1_2[10]"
# Returns ENST00000198947.1:c.1_10dup
variant27 = "ENST00000198947.1:C.1_2A[10]"
variant28 = "LRG_199t1:c.1_2ACT[8]"


def main():
    my_variant = TandemRepeats.parse_repeat_variant(variant28, "GRCh37", "all")
    my_variant.check_transcript_type()
    my_variant.reformat_reference()
    my_variant.check_genomic_or_coding()
    formatted = my_variant.reformat()
    print(f"Variant formatted: {formatted}")

    types_to_put_into_vv = ["ins", "dup"]
    if any(x in formatted for x in types_to_put_into_vv):
        validate = vval.validate(
            formatted, my_variant.build, my_variant.select_transcripts
        ).format_as_dict(with_meta=True)
        print(json.dumps(validate, sort_keys=True, indent=4, separators=(",", ": ")))


if __name__ == "__main__":
    main()
