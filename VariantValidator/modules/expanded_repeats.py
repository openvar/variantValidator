"""
NAME:          modules/expanded_repeats.py
AUTHORS:       Rebecca Locke (@rklocke) & Robert Wilson (@RSWilson1)
DATE:          18/02/22
INSTITUTION:   University of Manchester/Cambridge University Hospitals

DESCRIPTION:   This script contains the TandemRepeats class and methods,
               aiming to check the syntax and
               reformat tandem repeat variants for VariantValidator.
"""

# Importing Modules
import json
import re
import logging
import os
import VariantValidator

# Set up logger
logger = logging.getLogger(__name__)


# Established class for Tandem repeats
class TandemRepeats:
    """
    Class used for create instances of
    expanded repeat variants.
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
        variant_str,
        ref_type
    ):
        """
        This initialised an instance of the class with set class vars.

        Parameters
        ----------
        variant_str : str
            (Variant string i.e. LRG_199:g.1ACT[20])
        build : str
            Which genome reference the variant_string refers to e.g. GRCh37
        select_transcripts : str
            Return all possible transcripts or only select ones e.g. "all"
        Returns
        -------
        None. But a class instance of variant is created.
        """
        self.reference = reference
        self.prefix = prefix
        self.variant_position = variant_position
        self.repeat_sequence = repeat_sequence
        self.copy_number = copy_number
        self.after_the_bracket = after_the_bracket
        self.build = build
        self.select_transcripts = select_transcripts
        self.variant_str = variant_str
        self.ref_type = ref_type

    @classmethod
    def parse_repeat_variant(cls, variant_str, build, select_transcripts):
        """
        Summary
        -------
        Takes a variant string and breaks it into its constituent parts
        with regex to be processed in downstream functions,
        assigns them to class variables.
        Parameters
        ----------
            variant_str : str
                variant string e.g. "LRG_199:g.1ACT[20]A"
        Returns
        -------
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

        Example 1:
            >>>parse_repeat_variant("LRG_199:g.1ACT[20]A")
                "LRG_199", "g", "1", "ACT", "20", "A"
        Example 2:
            >>>parse_repeat_variant("NM_024312.4:c.1_10A[10]")
                "NM_024312.4", "c", "1_10" "A", "10", ""
        """
        logger.info(f"Parsing variant: parse_repeat_variant({variant_str})")
        # Strip any whitespace
        variant_str = variant_str.strip()

        # Check if square brackets included which indicate tandem repeat
        # variant

        if '[' in variant_str or ']' in variant_str:
            try:
                assert ":" in variant_str,\
                    (
                        f"Unable to identify a colon (:)"
                        f"in the variant description {variant_str}. "
                        f"A colon is required in HGVS variant descriptions to separate"
                        f"the reference accession from the reference type"
                        f" i.e. <accession>:<type>. e.g. :c"
                    )
            except AssertionError:
                logger.critical("A colon is required in the variant description. Ending program")
                raise
            else:
                reference, suffix = variant_str.split(":")
            try:
                assert ";" not in variant_str,\
                "A semi-colon is included in variant but alleles are not yet supported"
            except AssertionError:
                logger.critical(
                    "A semi-colon is included but alleles are not yet supported. Ending program"
                    )
                raise
            try:
                assert "," not in variant_str,\
                "A comma is included in variant but alleles are not yet supported"
            except AssertionError:
                logger.critical("A comma is included in variant but alleles "\
                                "are not yet supported. Ending program")
                raise
            # Find reference sequence used (g or c)
            var_type = re.search("^.*?(.*?)\\.", suffix)
            prefix = var_type.group(1).lower()
            # Get g or c position(s)
            # Extract bit between . and [ e.g. 1ACT
            pos_and_seq = suffix.split(".")[1].split("[")[0]
            try:
                assert re.search(
                "[a-z]+", pos_and_seq, re.IGNORECASE), \
                "Please ensure that the repeated sequence is included between "\
                "the position and number of repeat units, e.g. g.1ACT[20]"
            except AssertionError:
                logger.critical("Unable to identify a repeated sequence "\
                                "between position and number of repeats. Ending program")
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
                        ), f"Please ensure the start and the end of the "\
                           f"full repeat range is provided, "\
                           f"separated by an underscore"
                except AssertionError:
                    logger.critical("""Only one value in the range is provided.
                                     Ending program""")
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
            # Save anything after bracket so that mixed repeats are supported in future
            if re.search("\\](.*)", variant_str):
                after_brac = re.search("\\](.*)", variant_str)
                after_the_bracket = after_brac.group(1)
            else:
                after_the_bracket = ""

            ref_type = ""
        else:
            logger.info(
                "Unable to identify a tandem repeat. Ending program."
                "Check Format matches HGVS: "
                "(https://varnomen.hgvs.org/recommendations/DNA/variant/repeated/)"
            )
            return False
            #  This returns False to VV to indicate no tandem repeats present.

        return cls(
            reference,
            prefix,
            variant_position,
            repeat_sequence,
            copy_number,
            after_the_bracket,
            build,
            select_transcripts,
            variant_str,
            ref_type
        )

    def check_transcript_type(self):
        """
        Find transcript type and run relevant function for processing,
        that transcript type.
        N.B. Future development could instead store
        the transcript and replace it with refseq.

        Parameters
        ----------
        self.reference:str
            The reference from parse_variant_repeat
            i.e. Everything before the first colon.
        Returns
        -------
        None, prints variant type.

        Raises:
            Exception: (Error for unknown transcript type.)
        """
        logger.info(
            f"Checking transcript type: check_transcript_type({self.reference})"
        )
        if bool(re.match(r"^LRG", self.reference)):
            logger.info("Variant type: LRG variant")
            self.ref_type = "LRG"
        elif bool(re.match(r"^ENS", self.reference)):
            logger.info("Variant type: Ensembl variant")
            self.ref_type = "Ensembl"
        elif bool(re.match(r"NM", self.reference)):
            logger.info("Variant type: RefSeq variant")
            self.ref_type = "RefSeq"
        else:
            raise Exception(
                "Unknown transcript type present. " \
                "Try using the RefSeq transcript ID")

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
            ), """Please ensure the transcript or gene version is included
                  following a '.' after the transcript 
                  or gene name e.g. ENST00000357033.8"""
        return self.reference

    def check_genomic_or_coding(self):
        """Takes reference and works out what prefix type should be used
        Raises:
            AssertionError if wrong prefix type is used
        """
        logger.info(
            f"Checking prefix is consistent with reference: "\
            f"check_genomic_or_coding({self.reference},{self.prefix})"
        )
        if re.match(r"^LRG", self.reference):
            if "t" in self.reference:
                assert (
                    self.prefix == "c"
                ), """Please ensure variant type is coding 
                      if an LRG transcript is provided"""
            else:
                assert (
                    self.prefix == "g"
                ), "Please ensure variant type is genomic if LRG gene is used"
        elif re.match(r"^ENST", self.reference):
            assert (
                self.prefix == "c"
            ), """Please ensure variant type is coding 
                if an Ensembl transcript is provided"""
        elif re.match(r"^ENSG", self.reference):
            assert (
                self.prefix == "g"
            ), "Please ensure variant type is genomic if Ensembl gene is used"
        elif re.match(r"^NM", self.reference):
            assert (
                self.prefix == "c"
            ), """Please ensure variant type is coding 
                  if a RefSeq transcript is provided"""
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
        """
        Checks the position range given and
        updates it if it doesn't match the length of the repeated sequence
        and number of repeat units when full range is needed.
        Parameters
        ----------
            repeat_sequence (string): The repeated sequence e.g. "ACT"
            variant_position (string): The position of the variant e.g. "1" or "1_5"
            copy_number (string): The number of repeat units e.g. "20"
        Returns
        -------
            full_range (string): The full range supplied if correct or
            the full range updated if inputted range was incorrect, e.g. "1_20"
        """
        logger.info(
            f"Checking range given: "\
            f"check_positions_given({self.repeat_sequence}, "\
            f"{self.variant_position},{self.copy_number})"
        )
        start_range, end_range = self.variant_position.split("_")
        rep_seq_length = len(self.repeat_sequence)
        the_range = int(end_range) - int(start_range) + 1
        repeat_length = rep_seq_length * int(self.copy_number)
        if the_range == repeat_length:
            logger.info(
                f"Range given ({self.variant_position}) matches "\
                f"repeat sequence length and number of repeat units"
            )
            full_range = f"{start_range}_{end_range}"
        else:
            logger.warning(
                f"Warning: sequence range {self.variant_position}"\
                f"given does not match repeat unit sequence length "\
                f"and number of repeat units. Updating the range "\
                f"based on repeat sequence length and number of repeat units"
            )
            new_end_range = int(start_range) + repeat_length - 1
            full_range = f"{start_range}_{new_end_range}"
        return full_range

    def get_range_from_single_pos(self):
        """
        Gets full range of the variant if this is needed
        when a single start position is supplied
        """
        logger.info(
            f"Fetching the range from a given single position: "\
            f"get_range_from_single_pos({self.repeat_sequence},"\
            f"{self.copy_number},{self.variant_position}"
        )
        rep_seq_length = len(self.repeat_sequence)
        repeat_range = rep_seq_length * int(self.copy_number)
        the_end_range = int(self.variant_position) + repeat_range - 1
        full_range = f"{self.variant_position}_{the_end_range}"
        return full_range

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
            ), "Please ensure the transcript or gene version is \
                included following a '.' \
                after the transcript or gene name e.g. ENST00000357033.8"
        return self.reference

    def reformat_not_multiple_of_three(self):
        """
        Reformats coding variants (c.) to a dup or ins
        if they are not a multiple of three
        """
        logger.info(
            f"eformatting variant as not a multiple of three: "\
            f"reformat_not_multiple_of_three({self.repeat_sequence}, "\
            f"{self.variant_position},{self.reference},{self.prefix})"
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
                f"Warning: Repeated sequence is coding "\
                f"and is of length {rep_seq_length}, "\
                f"not a multiple of three! "\
                f"Updating variant description to a duplication"
            )
            reformatted = f"{self.reference}:{self.prefix}.{self.variant_position}dup"
        # Repeat of 2 bases should be an ins
        # with only first two nts given as range
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
                f"Warning: Repeated sequence is coding "\
                f"and is of length {rep_seq_length}, "\
                f"not a multiple of three! "\
                f"Updating variant description to insertion"\
            )
            reformatted = (
                f"{self.reference}:{self.prefix}.{position}ins{expanded_rep_seq}"
            )
        return reformatted

    def reformat(self):
        """Reformats and returns final formatted variant as a string"""
        logger.info(
            f"Reformatting variant: reformat({self.repeat_sequence}, "\
            f"{self.after_the_bracket}, "\
            f"{self.prefix}, {self.variant_position}, {self.copy_number})"
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
                f"No information should be included after "\
                f"the number of repeat units. "\
                f"Currently '{self.after_the_bracket}'' is included. "\
                f"This will be removed as mixed repeats are not currently supported."
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
                #Uncomment if you want to always have range in final format
                # else:
                #     self.variant_position = self.get_range_from_single_pos()
                final_format = f"{self.reference}:{self.prefix}."\
                f"{self.variant_position}{self.repeat_sequence}"\
                f"[{self.copy_number}]"
        # Reformat g. variants
        else:
            if "_" in self.variant_position:
                self.variant_position = self.check_positions_given()
            #Uncomment if you want to always have range in final format
            # else:
            #     self.variant_position = self.get_range_from_single_pos()
            final_format = f"{self.reference}:"\
            f"{self.prefix}.{self.variant_position}"\
            f"{self.repeat_sequence}[{self.copy_number}]"
        return final_format

    def simple_split_string(self):
        """
        Splits the string into two parts divided by the colon (:).
        Useful for segregating the two distict elements of the variant string
        And for troubleshooting.
        Parameters
        ----------
        self.variant_str:str
                (Variant string i.e. LRG_199:g.1ACT[20])

        Returns
        -------
        self.prefix:str
            String for the transcript of the variant. I.e. NM_40091.5
        self.suffix:str
            String for the remaining variant string for further processing.
        """
        self.begining = self.variant_str.split(":")[0]
        self.end = ":" + self.variant_str.split(":")[1]
        return self.begining, self.end


def convert_tandem(variant_str, build, my_all):
    expanded_variant = TandemRepeats.parse_repeat_variant(variant_str, build, my_all)
    if expanded_variant is False:
        return False
    else:
        expanded_variant_string = expanded_variant.reformat()
        return True



# <LICENSE>
# Copyright (C) 2016-2023 VariantValidator Contributors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>
