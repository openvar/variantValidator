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
import re
import logging
import copy
from vvhgvs.assemblymapper import AssemblyMapper
# AlignmentMapper allows us to skip redundant validation & ref filling on c<->n
# the no normalisation flags for AssemblyMapper only skip the validation & post
# map replacement but ref will always be subject to initial pre map filling
import vvhgvs
from vvhgvs.alignmentmapper import AlignmentMapper
from vvhgvs.location import BaseOffsetInterval, Interval

from .hgvs_utils import hgvs_delins_parts_to_hgvs_obj, \
        _hgvs_offset_pos_from_str_in
# Set up logger
logger = logging.getLogger(__name__)


class RepeatSyntaxError(Exception):
    """Raised when the syntax of the expanded repeat is incorrect"""

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
        self.intronic_g_reference = False
        self.prefix = prefix
        self.variant_position = variant_position
        self.repeat_sequence = repeat_sequence
        self.copy_number = copy_number
        self.after_the_bracket = after_the_bracket
        self.build = build
        self.select_transcripts = select_transcripts
        self.variant_str = variant_str
        self._c_to_n_tx_maper = None # only valid for c<->n type mappings
        self.g_strand = 1 # only valid for intronic +/-1
        self.evm = False
        self.no_norm_evm = False
        self.genomic_conversion = None
        self.original_position = None
        self.reference_sequence_bases = None

        # Define the wobble bases map with proper regex
        self._wobble_bases_map = {
            'A': r'A',       # Adenine
            'C': r'C',       # Cytosine
            'G': r'G',       # Guanine
            'T': r'T',       # Thymine
            'U': r'U',       # Uracil (in RNA)
            'R': r'[AG]',    # A or G (puRine)
            'Y': r'[CT]',    # C or T (pYrimidine)
            'S': r'[GC]',    # G or C (Strong interaction)
            'W': r'[AT]',    # A or T (Weak interaction)
            'K': r'[GT]',    # G or T (Keto)
            'M': r'[AC]',    # A or C (aMino)
            'B': r'[CGT]',   # C or G or T (not A)
            'D': r'[AGT]',   # A or G or T (not C)
            'H': r'[ACT]',   # A or C or T (not G)
            'V': r'[ACG]',   # A or C or G (not T)
            'N': r'[ACGT]',  # Any base (A or C or G or T)
        }

    @classmethod
    def parse_repeat_variant(cls, variant_str, build, select_transcripts, validator):
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
            if not ( '[' in variant_str and ']' in variant_str ):
                raise RepeatSyntaxError(
                    f"RepeatSyntaxError: This variant {variant_str} contains a square bracket "
                    "and as such has been detected as a possible tandem repeat, despite this the "
                    "variant in question is missing a matched bracket pair. Please either remove "
                    "the stray square bracket, or add a paired end, which should produce an "
                    "output pair that looks something like [20] or [5].")
            # This along with some of the later stages will fail for miss-formatted variants,
            # we rely on the previous VV code to check for variants missing : or a '.' type.
            # Extra checking will need to be done before using this code if is is being used
            # stand alone from the rest of VV.
            reference, _sep, suffix = variant_str.partition(":")

            # Find reference sequence used (g, n, c etc.)
            prefix, _sep, pos_edit = suffix.partition(".")
            prefix = prefix.lower()
            # Get position/span by extracting the bit between '.' and [ e.g. 1ACT
            pos_and_seq, _sep, post_bracket  = pos_edit.partition("[")
            rep_seq = re.search("[AaCcTtGgUuMmNnRrYyKkSsWwHhBbVvDd]+", pos_and_seq)
            if not rep_seq:
                raise RepeatSyntaxError(
                    "RepeatSyntaxError: Ensure that the repeated sequence is included between "
                    "the variant position and the number of repeat units, e.g. g.1_3ACT[20]")
            for char in pos_and_seq:
                if char.isalpha():
                    if char not in "AaCcTtGgUuMmNnRrYyKkSsWwHhBbVvDd":
                        raise RepeatSyntaxError(
                            "RepeatSyntaxError: Please ensure the repeated sequence includes"
                            " only Aa, Cc, Tt, Gg, Uu or a valid IUPAC nucleotide code from "
                            "https://genome.ucsc.edu/goldenPath/help/iupac.html")
            repeat_sequence = rep_seq.group()
            variant_position = pos_and_seq[:rep_seq.start()]
            if prefix == 'g':
                if '_' in variant_position:
                    start, _sep, end = variant_position.partition('_')
                    start = int(start)
                    end = int(end)
                else:
                    start = int(variant_position)
                    end = start
                variant_position =  Interval(
                    start=vvhgvs.location.SimplePosition(base=start),
                    end=vvhgvs.location.SimplePosition(base=end))
            else:
                if '_' in variant_position:
                    start, _sep, end = variant_position.partition('_')
                    variant_position = _hgvs_offset_pos_from_str_in(
                            start,None,end=end,ref_type=prefix)
                elif variant_position:
                    variant_position = _hgvs_offset_pos_from_str_in(
                            variant_position,1,ref_type=prefix)
                variant_position = BaseOffsetInterval(
                        start=variant_position[0],
                        end=variant_position[1])

            # Get number of unit repeats
            copy_number, _sep, after_the_bracket = post_bracket.partition("]")
            # Save anything after bracket so that mixed repeats are supported in future
            if not after_the_bracket:
                after_the_bracket = ""
        else:
            logger.info(
                "Unable to identify a tandem repeat, if a tandem repeat is "
                "expected then please check that the format matches HGVS: "
                "(https://varnomen.hgvs.org/recommendations/DNA/variant/repeated/)"
            )
            return False
            #  This returns False to VV to indicate no tandem repeats present.

        if "LRG" in reference:
            if "t" in reference:
                reference = validator.db.get_refseq_transcript_id_from_lrg_transcript_id(reference)
            else:
                reference = validator.db.get_refseq_id_from_lrg_id(reference)

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
        )

    def reformat_reference(self):
        """Reformats the reference sequence name"""
        logger.info(f"Reformatting reference: reformat_reference({self.reference})")
        if self.reference.startswith("ENS") or self.reference.startswith("N"):
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
        if self.reference.startswith("ENST"):
            assert (
                self.prefix == "c"
            ), """Please ensure variant type is coding
            if an Ensembl transcript is provided"""
        elif self.reference.startswith("NM"):
            assert (
                self.prefix == "c"
            ), """Please ensure variant type is coding
                  if a RefSeq transcript is provided"""
        elif self.reference.startswith("NC"):
            assert (
                self.prefix == "g"
            ), "Please ensure variant type is genomic if RefSeq chromosome is used"
        elif self.reference.startswith("NG"):
            assert (
                self.prefix == "g"
            ), "Please ensure variant type is genomic if RefSeq gene is used"
        elif self.reference.startswith("NR"):
            assert (
                self.prefix == "n"
            ), "Please ensure variant type is non-coding if NR transcript is used"

    def check_positions_given(self, validator):
        """
        Checks the position range matches the genomic reference
        Assumes that range is n type coordinates, 1 based WRT reference start
        """
        logger.info(
            f"Checking range given: "\
            f"check_positions_given({self.repeat_sequence}, "\
            f"{str(self.variant_position)}, {self.copy_number})"
        )
        ref = self.reference
        if self.intronic_g_reference:
            ref = self.intronic_g_reference
        reference_repeat_sequence = validator.sf.fetch_seq(
                ref,
                self.variant_position.start.base-1,
                self.variant_position.end.base)
        logger.info(f"Reference repeat sequence: {reference_repeat_sequence}")

        # Check if the length of reference_repeat_sequence is a multiple of the length of query_str
        if len(reference_repeat_sequence) % len(self.repeat_sequence) != 0:
            raise RepeatSyntaxError(
                f"RepeatSyntaxError: The stated range {str(self.variant_position)} is not a "
                f"multiple of the length of the provided repeat sequence {self.repeat_sequence}")

        # Create a new string by repeating query_str enough times
        repeated_str = self.repeat_sequence * (len(reference_repeat_sequence)
                                               // len(self.repeat_sequence))

        # Do not replace this regex for a stored version. Critical
        regex = self.build_regex(repeated_str)
        match = regex.search(reference_repeat_sequence)
        try:
            match.group()
            logger.info(f"Regex matched {match.group()}")
            self.reference_sequence_bases = match.group()
            return
        except AttributeError:
            # We previously had a case for an empty repeated_str here, but that should always make
            # a match.group() on a valid input and regex.search() should fail otherwise
            raise RepeatSyntaxError(
                f"RepeatSyntaxError: The repeat sequence does not match the reference sequence at "
                f"the given position {str(self.variant_position)}, expected {repeated_str} but the "
                f"reference is {reference_repeat_sequence} at the specified position")

    def get_valid_n_or_g_range_from_input_pos(self, validator):
        """
        Substitute for get_range_from_single_or_start_pos without full re-build.

        Works by converting c->n or c/n->g without rebuilding, leaves other
        input untouched. Used to keep the original input range, mainly for
        testing. Mutually exclusive with get_range_from_single_or_start_pos,
        using both will break c type inputs. This also avoids ALL checks for
        c/n intronic to g mapping, for the same reasons (testing or deliberate
        maintenance of the original) so should be combined with extra tests for
        the sequence state if it is being used for user visible output.
        """
        pos = self.variant_position
        if isinstance(self.variant_position, BaseOffsetInterval) and (
                pos.start.offset or pos.end.offset):
            if not self.no_norm_evm:
                # if we normalise at the wrong points we get the wrong ref seq back
                # This mapper does call functions to "fill" and "replace" ref but
                # the first runs without effect on intronic inputs, and the output
                # genomic ref is used, so leave in the more complex form for now
                self.no_norm_evm = AssemblyMapper(
                        validator.hdp,
                        assembly_name=self.build,
                        alt_aln_method=validator.alt_aln_method,
                        normalize=False,
                        replace_reference=True
                        )
            var_n = hgvs_delins_parts_to_hgvs_obj(
                    self.reference,
                    self.prefix,
                    self.variant_position,
                    '',
                    f"{self.repeat_sequence * int(self.copy_number)}")
            self.original_position = copy.deepcopy(self.variant_position)
            intronic_genomic_variant = self.no_norm_evm.t_to_g(var_n)
            self.genomic_conversion = intronic_genomic_variant
            self.intronic_g_reference = intronic_genomic_variant.ac
            pos = intronic_genomic_variant.posedit.pos
            # untested ref!
            self.reference_sequence_bases = self.repeat_sequence * int(
                (pos.end - pos.start + 1) / len(self.repeat_sequence))
        elif self.prefix == "c":
            # don't do c->n mapping on t_to_g targets
            # if tx is coding and var is not c we get an error
            self._get_c_tx_info(validator)
            pos = self.convert_c_to_n_coordinates()
            # untested ref!
            self.reference_sequence_bases = self.repeat_sequence * int(
                (pos.end - pos.start + 1) / len(self.repeat_sequence))

        self.variant_position = pos


    def get_range_from_single_or_start_pos(self, validator):
        """
        Gets full range of the variant if this is needed,
        Used both when a single start position is supplied
        and to rebuild ranges for validation purposes.

        Currently this is the only place that intronic g locations are derived

        Uses: self.variant_position
              validator (a VariantValidator object for data fetch, intronic
              genomic mappings, etc.)
        Sets: self.variant_position, to the start position of the repeat (if it
              is not already done). In the case of intronic variants with -1
              strand mapping this will be the end of the within transcript
              position pair.
        Returns: The n based coordinate span of the repeat region, as found

        """
        # due to the nature of hgvs object pos we always have end but it may == start
        start_pos = self.variant_position.start
        end_pos = self.variant_position.end
        if isinstance(self.variant_position, BaseOffsetInterval) and (
                start_pos.offset or end_pos.offset):
            logger.info(
                "Re-fetching the range using adaptions for exon handling " +
                f"using the range {str(self.variant_position)} with a " +
                f"copy number of {self.copy_number} and repeat of "+
                self.repeat_sequence
            )
        elif end_pos:
            logger.info(
                "Re-fetching the range from the start position of a given " +
                f"range using {start_pos} from {str(self.variant_position)} with a" +
                f" copy number of {self.copy_number} and repeat of "+
                self.repeat_sequence
            )
            self.variant_position.end = copy.copy(end_pos)

        # Get the full reference sequence range of the variant
        # the within_ref_pos here is the start position in 0 based coordinates

        self.original_position = copy.deepcopy(self.variant_position)
        if not (isinstance(self.variant_position, BaseOffsetInterval) and (
                start_pos.offset or end_pos.offset)):
            if not self.prefix == "c":
                within_ref_pos = self.variant_position.start.base - 1
            else:
                self._get_c_tx_info(validator)
                pos = self.convert_c_to_n_coordinates()
                within_ref_pos = pos.start.base - 1
        else:
            if not self.no_norm_evm:
                # if we normalise at the wrong points we get the wrong ref seq back
                self.no_norm_evm = AssemblyMapper(
                        validator.hdp,
                        assembly_name=self.build,
                        alt_aln_method=validator.alt_aln_method,
                        normalize=False,
                        replace_reference=True
                        )
            # Handle a range variant input, given current function's purpose we
            # just dump range and re-build from the first base

            if self.variant_position.start != self.variant_position.end:
                seq_check = hgvs_delins_parts_to_hgvs_obj(
                        self.reference,
                        self.prefix,
                        self.variant_position,
                        "",
                        f"{self.repeat_sequence * int(self.copy_number)}")

                self.original_position = copy.deepcopy(self.variant_position)
                intronic_genomic_variant = self.no_norm_evm.t_to_g(seq_check)
                self.genomic_conversion = intronic_genomic_variant
                # get ref length for complex cases (needed for checking)
                ref_len = len(intronic_genomic_variant.posedit.edit.ref)
                ref_copies = int(ref_len/len(self.repeat_sequence))
                # Check the exon boundaries by using normalising mapper
                self.check_exon_boundaries(validator)
                seq_check = self.evm.g_to_t(intronic_genomic_variant, self.reference)
                regex = self.build_regex(self.repeat_sequence)

                # Check for occurrences of the repeat sequence using regex
                matches = list(regex.finditer(seq_check.posedit.edit.ref))
                repeat_count = len(matches)
                if matches[0].start() != 0:
                    # Create a new string by repeating query_str enough times
                    repeated_str = self.repeat_sequence * (len(seq_check.posedit.edit.ref)
                                                           // len(self.repeat_sequence))

                    raise RepeatSyntaxError(
                        "RepeatSyntaxError: The repeat sequence does not match the reference "
                        f"sequence at the given position {str(self.variant_position)}, expected"
                        f" {repeated_str} but the reference is {seq_check.posedit.edit.ref} at "
                        "the specified position")
                if repeat_count != ref_copies:
                    raise RepeatSyntaxError(
                        "RepeatSyntaxError: The repeat sequence does not match the expected "
                        f"copy number at position {str(self.variant_position)}. Expected "
                        f"{int(ref_copies)} occurrences of '{self.repeat_sequence}', but"
                        f" found {repeat_count} in reference sequence "
                        f"'{seq_check.posedit.edit.ref}'.")

                self.g_strand = validator.hdp.get_tx_exons(
                        self.reference, intronic_genomic_variant.ac,
                        validator.alt_aln_method)[0][3]
                # this is re-expanded later from the first base, to handle truncated input
                if  self.g_strand == -1:
                    self.variant_position.start = copy.copy(self.variant_position.end)
                else:
                    self.variant_position.end = copy.copy(self.variant_position.start)
            else:
                # Create a variant for mapping to the genome containing the whole repeat, we
                # used to use only the first base of the repeat but this breaks on -1 mapping
                # transcripts with multi-base repeats
                length = len(self.repeat_sequence)
                pos = copy.copy(self.variant_position)
                if length == 1:
                    pos = self.variant_position
                elif isinstance(self.variant_position, BaseOffsetInterval) and \
                        self.variant_position.start.offset:
                    pos.end.offset = pos.end.offset + length - 1

                intronic_variant = hgvs_delins_parts_to_hgvs_obj(
                        self.reference,
                        self.prefix,
                        pos,
                        self.repeat_sequence,
                        self.repeat_sequence,
                        )
                intronic_genomic_variant = self.no_norm_evm.t_to_g(intronic_variant)
                self.g_strand = validator.hdp.get_tx_exons(
                        intronic_variant.ac, intronic_genomic_variant.ac,
                        validator.alt_aln_method)[0][3]
                self.variant_position = intronic_genomic_variant.posedit.pos

            self.intronic_g_reference = intronic_genomic_variant.ac
            self.genomic_conversion = intronic_genomic_variant

            # Check exon boundaries
            self.check_exon_boundaries(validator)
            within_ref_pos = intronic_genomic_variant.posedit.pos.start.base - 1
            if self.g_strand == -1:
                self.repeat_sequence = self.reverse_complement(self.repeat_sequence)
        # validate that the existing range matches the repeat given
        # then (re-)expand to full length and store in n type 1 based coordinates
        self.check_reference_sequence(validator, within_ref_pos)
        ref_start_position, ref_end_position = self.get_reference_range(validator, within_ref_pos)
        ref_start_position = ref_start_position + 1
        full_range = copy.copy(self.variant_position)
        full_range.start.base = ref_start_position
        full_range.end.base = ref_end_position

        # Map the full genomic range in the description
        if self.genomic_conversion is not None:
            self.genomic_conversion.posedit.edit.ref = ""
            self.genomic_conversion.posedit.edit.alt = ""
            self.genomic_conversion.posedit.pos.start.base = ref_start_position
            self.genomic_conversion.posedit.pos.end.base = ref_end_position
        return full_range

    def build_regex(self, sequence):
        """Convert an IUPAC-coded sequence into a regex pattern."""
        regex_pattern = "".join(self._wobble_bases_map.get(base, base) for base in sequence)
        return re.compile(regex_pattern)

    def check_reference_sequence(self, validator, within_ref_pos):
        """
        Check that the current within_ref_pos is a valid start for the given repeat.
        return is unused, we raise a RepeatSyntaxError if the match fails
        """
        start = within_ref_pos
        end = within_ref_pos + len(self.repeat_sequence)
        ref = self.reference
        if self.intronic_g_reference:
            ref = self.intronic_g_reference
        requested_sequence = validator.sf.fetch_seq(ref, start, end)
        logger.info(f"Requested sequence: {requested_sequence} from {ref} at {start}-{end}")

        # Critical, do not use cached regex
        regex = self.build_regex(self.repeat_sequence)
        match = regex.search(requested_sequence)
        try:
            match.group()
            return
        except AttributeError:
            raise RepeatSyntaxError(
                f"RepeatSyntaxError: The provided repeat sequence {self.repeat_sequence } does not "
                f"match the reference sequence {requested_sequence} at the given position "
                f"{within_ref_pos+1}_{within_ref_pos + len(self.repeat_sequence)} of "
                f"reference sequence {ref}")

    def get_reference_range(self, validator, within_ref_pos):
        """
        Get the full range of the variant, starting with a within reference 0
        based start position that should line up with the start of the repeat
        sequence.
        """
        # get sequence ref
        ref = self.reference
        if self.intronic_g_reference:
            ref = self.intronic_g_reference
        logger.info(
            f"Getting the full range of the variant: "
            f"get_reference_range({ref}, {str(self.variant_position)})"
        )

        # Get the full range of the reference repeat sequence
        start_position = None
        end_position = None
        current_position = None
        while start_position is None:
            if current_position is None:
                current_position = within_ref_pos
            last_position = copy.copy(current_position)
            current_position -= len(self.repeat_sequence)
            requested_sequence = validator.sf.fetch_seq(
                ref, current_position,
                current_position + len(self.repeat_sequence))
            regex = self.build_regex(self.repeat_sequence)
            match = regex.search(requested_sequence)
            try:
                match.group()
                continue
            except AttributeError:
                start_position = last_position

        current_position = None
        while end_position is None:
            if current_position is None:
                current_position = within_ref_pos
            last_position = copy.copy(current_position)
            current_position += len(self.repeat_sequence)
            requested_sequence = validator.sf.fetch_seq(
                ref, current_position,
                current_position + len(self.repeat_sequence))
            regex = self.build_regex(self.repeat_sequence)
            match = regex.search(requested_sequence)
            try:
                match.group()
                continue
            except AttributeError:
                end_position = last_position + len(self.repeat_sequence)
        return start_position, end_position

    def reformat(self, validator):
        """Reformats and returns final formatted variant as a string"""
        logger.info(
            f"Reformatting variant: reformat({self.repeat_sequence}, "\
            f"{self.after_the_bracket}, "\
            f"{self.prefix}, {str(self.variant_position)}, {self.copy_number})"
        )

        if not self.copy_number.isdecimal():
            raise RepeatSyntaxError(
                "RepeatSyntaxError: The number of repeat units included between"
                " square brackets must be numeric")

        # Update the repeated sequence to be upper case
        self.repeat_sequence = self.repeat_sequence.upper()
        # test for non-matching chars
        if re.search("[^ACTGUMRYKSWHBVDN]", self.repeat_sequence):
            raise RepeatSyntaxError(
                "RepeatSyntaxError: Please ensure the repeated sequence"
                " includes only Aa, Cc, Tt, Gg, Uu or a valid IUPAC nucleotide code from "
                "https://genome.ucsc.edu/goldenPath/help/iupac.html")

        if self.after_the_bracket != "":
            raise RepeatSyntaxError(
                f"No information should be included after "
                f"the number of repeat units. "
                f"Currently '{self.after_the_bracket}'' is included. "
            )

        if isinstance(self.variant_position, BaseOffsetInterval) and (
            self.variant_position.start.offset or self.variant_position.end.offset):

            # Create easy variant mapper (over variant mapper) and splign locked evm
            self.evm = AssemblyMapper(validator.hdp,
                                      assembly_name=self.build,
                                      alt_aln_method=validator.alt_aln_method,
                                      normalize=True,
                                      replace_reference=True
                                      )

        # Add coordinates from start position (need range for intronic fall-back)
        self.variant_position = self.get_range_from_single_or_start_pos(validator)
        self.check_positions_given(validator)
        # Map back to transcript if intronic
        if self.genomic_conversion is not None:
            transcript_variant = self.evm.g_to_t(self.genomic_conversion, self.reference)
            self.reference = transcript_variant.ac
            self.variant_position = transcript_variant.posedit.pos
            if self.g_strand == -1:
                self.repeat_sequence = self.reverse_complement(self.repeat_sequence)
                self.reference_sequence_bases =  self.reverse_complement(
                        self.reference_sequence_bases)
        else:
            self._get_c_tx_info(validator)
            self.variant_position = self.convert_n_to_c_coordinates()
        # set inside of get_range_from_single_or_start_pos or equivalent
        ref = self.reference_sequence_bases

        final_hgvs = hgvs_delins_parts_to_hgvs_obj(
                self.reference,
                self.prefix,
                self.variant_position,
                ref,
                f"{self.repeat_sequence * int(self.copy_number)}")
        final_hgvs.posedit.expanded_rep = self.repeat_sequence

        return final_hgvs

    def _get_c_tx_info(self,validator):
        """
        Get the transcript info needed for turning n type from start of
        transcript coordinates into c type CDS relative coordinates.
        Only calls out to the database once per variant.
        Parameters
        ----------
        validator a variant validator object to use for data fetch
        Returns
        -------
        None, data is stored in self.cds_start and self.cds_end instead
        """
        if not self.prefix == "c":
            return
        if self._c_to_n_tx_maper is not None:
            return
        self._c_to_n_tx_maper = AlignmentMapper(
                validator.hdp, self.reference, self.reference, "transcript")


    def convert_n_to_c_coordinates(self):
        """
        Applies CDS based c type offset to a n type variant position
        used for c transcripts only!
        _get_c_tx_info must be called before this function is
        used.
        """
        logger.info(
            "Applying c type offset to n type coordinates: " +
            f"convert_n_to_c_coordinates({str(self.variant_position)})"
        )
        if not self.prefix == 'c':
            return self.variant_position
        return self._c_to_n_tx_maper.n_to_c(self.variant_position)

    def convert_c_to_n_coordinates(self, pos=None):
        """
        Removes the offset from c type offset variant positions
        _get_c_tx_info must be called before this function is
        used.
        """
        if pos is None:
            pos = self.variant_position
        logger.info(
            f"Removing offset: remove_offset({str(self.variant_position)})"
        )
        if not self.prefix == 'c':
            return pos
        return self._c_to_n_tx_maper.c_to_n(pos)

    def reverse_complement(self, dna_seq):
        """
        Reverse complement a DNA string using the
        :param dna_seq:
        :return: (str)
        """
        reverse_seq = dna_seq[::-1]
        return "".join([{"G": "C", "T": "A", "A": "T", "C": "G",
                                         "g": "c", "t": "a", "a": "t", "c": "g"}[base]
                                        for base in list(reverse_seq)])

    def check_exon_boundaries(self,validator):
        """
        Check the boundaries of intronic variants are correctly stated
        """
        logger.info(
            "Checking intronic variant boundaries: "+
            f"check_exon_boundaries({str(self.original_position)})"
        )
        # Return without error if no boundaries used
        if not isinstance(self.original_position, BaseOffsetInterval) or (
            isinstance(self.original_position, BaseOffsetInterval) and not (
                self.original_position.start.offset or
                self.original_position.end.offset)):
            return

        exon_data = validator.hdp.get_tx_exons(self.reference,
                                               self.genomic_conversion.ac,
                                               validator.alt_aln_method)
        transcript_exon_pos = []
        if  self.prefix == 'c':
            self._get_c_tx_info(validator)
            pos = self.convert_c_to_n_coordinates(pos=self.original_position)
        else:
            pos = self.original_position

        for exon in exon_data:
            transcript_exon_pos.append(exon['tx_end_i'])
        def check_exon_pos(exon_pos):
            if not exon_pos.offset:
                return True
            bad_exon = False
            if exon_pos.offset < 0 and exon_pos.base -1 not in transcript_exon_pos:
                bad_exon = True
            elif exon_pos.offset > 0 and exon_pos.base not in transcript_exon_pos:
                bad_exon = True
            if bad_exon:
                raise RepeatSyntaxError(
                    f"ExonBoundaryError: Position {str(self.original_position)} " +
                    "does not correspond with an exon boundary for transcript "
                    + self.reference)
            return True
        check_exon_pos(pos.start)
        # skip "end" if we don't have a range
        if pos.end == pos.start:
            return
        check_exon_pos(pos.end)


def convert_tandem(variant, validator, build, my_all):
    "convenience function to encapsulate TandemRepeats->VV integration"
    try:
        expanded_variant = TandemRepeats.parse_repeat_variant(
                variant.quibble, build, my_all, validator)
    except AttributeError:
        expanded_variant = TandemRepeats.parse_repeat_variant(
                variant, build, my_all, validator)

    if expanded_variant is False:
        return False
    expanded_var_hgvs_obj = expanded_variant.reformat(validator)

    try:
        variant.expanded_repeat = {
                "variant": expanded_var_hgvs_obj,
                "position": expanded_variant.variant_position,
                "copy_number": expanded_variant.copy_number,
                "repeat_sequence": expanded_variant.repeat_sequence,
                "reference": expanded_variant.reference,
                "prefix": expanded_variant.prefix,
                "reference_sequence_bases": expanded_variant.reference_sequence_bases}
        logger.info(f"variant.expanded_repeat: {variant.expanded_repeat}")
    except AttributeError:
        expanded_repeat = {"variant": expanded_var_hgvs_obj,
                           "position": expanded_variant.variant_position,
                           "copy_number": expanded_variant.copy_number,
                           "repeat_sequence": expanded_variant.repeat_sequence,
                           "reference": expanded_variant.reference,
                           "prefix": expanded_variant.prefix,
                           "reference_sequence_bases": expanded_variant.reference_sequence_bases}
        logger.info(f"expanded_repeat: {expanded_repeat}")
        return expanded_repeat
    return True


# <LICENSE>
# Copyright (C) 2016-2025 VariantValidator Contributors
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
