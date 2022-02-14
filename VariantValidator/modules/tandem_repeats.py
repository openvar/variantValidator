
"""
This module aims to expand the functionality of VV into expanded repeat
variants.
It achieves this by using the ex_repeat_var class with associated functions.

Authors: Robert Wilson (@RSWilson1) and Rebecca Locke (@).
"""
# Importing Modules

import re
import logging
import os

# Get path the script is run in
CURRENT_DIR = os.path.abspath(os.getcwd())

# Create and configure logger
LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"

# Set level to debug, format with date and time and re-write file each time
logging.basicConfig(
    filename=f'{CURRENT_DIR}/expanded_repeats.log',
    level=logging.DEBUG,
    format=LOG_FORMAT,
    filemode='w')

log = logging.getLogger("main log")

# Established class for expanded repeats


class ex_repeat_var:
    """
    Class used for expanded repeat variants.
    These related functions split the variant string into
    its constiuents for further analysis.
    """

    def __init__(self, variant_string, build) -> None:
        self.variant_string = variant_string.strip()  # remove whitespace
        self.build = build
        self.type = None
        self.prefix = None
        self.suffix = None
        self.no_repeats = None
        self.after_brac = None

    def check_expanded_repeat_diverging(self):
        """
        This takes a variant string and breaks it into its constituents.
        This isolates the constituents with regex.

        Paramaters
        ----------
        my_variant:str
            (Variant string i.e. LRG_199:g.1ACT[20])
        Returns
        -------
        Prints out constituents and assigns them to variables for further processing.
        Or, returns False if my_variant is not compatible.
        """
        print("Running Diverging version")
        if "[" or "]" in self.variant_string:
            # Check which transcript type is present.
            if bool(re.search("^LRG", self.variant_string)):
                self.type = "LRG"
            elif bool(re.search("^ENSG", self.variant_string)):
                self.type = "Ensembl"
            elif bool(re.search("^NM", self.variant_string)):
                self.type = "RefSeq"
            if ":" not in self.variant_string:
                print("Unable to identify a colon (:) in the variant description.\
                      A colon is required in HGVS variant")
                log.info("Unable to identify a colon (:) in the variant description.\
                      A colon is required in HGVS variant")
            else:
                self.prefix = self.variant_string.split(":")[0]
                print(f'Variant prefix: {self.prefix}')
                log.info("Variant OK. Splitting Variant based on transcript type.")
                if self.type == "LRG":
                    # Check if underscore after LRG is included
                    if "_" not in self.prefix:
                        # Add in underscore between LRG and number
                        self.prefix = re.sub(
                            r"(?i)(?<=[a-z])(?=\d)", '_', self.prefix)
                        print(f'Updated prefix: {self.prefix}')
                self.suffix = ":" + self.variant_string.split(":")[1]
                # Find whether genomic or coding
                variant_type = re.search(':(.*?)\.', self.suffix)
                print(f'Variant type: {variant_type.group(1)}')
                # Get g or c position
                var_position = re.search('\.(.*?)[ACTG]', self.suffix)
                print(f'Variant position: {var_position.group(1)}')
                if "_" in var_position.group(1):
                    start_range, end_range = var_position.group(1).split("_")
                    print(start_range)
                    print(end_range)
                    rep_seq = re.search('\.[0-9]+_[0-9]+(.*?)\[', self.suffix)
                    print(f'Repeated sequence: {rep_seq.group(1)}')
                else:
                    rep_seq = re.search('\.[0-9]+(.*?)\[', self.suffix)
                    print(f'Repeat seq without range: {rep_seq.group(1)}')
                # Get number of unit repeats
                self.no_repeats = re.search('\[(.*?)\]', self.suffix).group(1)
                print(f'Number of unit repeats: {self.no_repeats}')
                # Get anything after ] to check if extra unsupported information.
                # Or to process further for supporting other syntaxes.
                self.after_brac = re.search('\](.*)', self.suffix).group(1)
                print(f'Anything after bracket: {self.after_brac}')
        else:
            print("No expanded repeat present.")
            log.info("No Expanded repeat present, the presence of [] is\
                     essential of classifying expanded repeats\
                     please check syntax on HGVS website and try again.")

    def get_transcript_type(self):
        """
        [Find transcript type and run relevant function for processing,
        that transcript type. N.B. Future development could instead store
        the transcript and replace it with refseq.]

        Parameters
        ----------
        self.variant_string:str
                (Variant string i.e. LRG_199:g.1ACT[20])

        Returns
        -------
        Returns results from relevant check_expanded_repeat_X
        where X is the transcript type.

        Raises:
            NameError: (Error for unknown transcript type.)
        """
        if bool(re.search("^LRG", self.variant_string)):
            self.type = "LRG"
            return self.type
        elif bool(re.search("^ENSG", self.variant_string)):
            self.type = "Ensembl"
            return self.type
        elif bool(re.search("NM", self.variant_string)):
            self.type = "RefSeq"
            return self.type
        else:
            raise NameError('Unknown transcript type present. \
                            Try RefSeq transcript ID')

    def split_var_string(self):
        """
        Splits the string into two parts divided by the colon (:).
        Parameters
        ----------
        self.variant_string:str
                (Variant string i.e. LRG_199:g.1ACT[20])

        Returns
        -------
        self.prefix:str
            String for the transcript of the variant. I.e. NM_40091.5
        self.suffix:str
            String for the remaining variant string for further processing.
        """
        self.prefix = self.variant_string.split(":")[0]
        self.suffix = ":" + self.variant_string.split(":")[1]
        return self.prefix, self.suffix

    def get_transcript_name(self):
        """
        Parameters
        ----------
        self.variant_string:str
                Variant string submitted i.e. LRG_199:g.1ACT[20]

        Returns
        -------
        self.variant_string:str
            Variant string submitted i.e. LRG_199:g.1ACT[20]
        self.type:str
            The transcript type, i.e. NCBI RefSeq, Ensembl, or LRG.
        self.prefix:str
            String for the transcript of the variant. I.e. NM_40091.5
        self.suffix:str
            String for the remaining variant string for further processing.
        """
        self.prefix = self.variant_string.split(":")[0]
        self.suffix = ":" + self.variant_string.split(":")[1]
        print(f'Variant prefix: {self.prefix}')
        if self.type == "LRG":
            # Check if underscore after LRG is included
            if "_" not in self.prefix:
                # Add in underscore between LRG and number
                self.prefix = re.sub(r"(?i)(?<=[a-z])(?=\d)", '_', self.prefix)
                print(f'Updated prefix: {self.prefix}')
        print(f'Variant string suffix: {self.suffix}')
        return self.variant_string, self.type, self.prefix, self.suffix

    def get_variant_location(self):
        """
        Aim
        ----------


        Parameters
        ----------
        self.variant_string:str
                Variant string submitted i.e. LRG_199:g.1ACT[20]

        Returns
        -------
        self.variant_string:str
            Variant string submitted i.e. LRG_199:g.1ACT[20]
        self.type:str
            The transcript type, i.e. NCBI RefSeq, Ensembl, or LRG.
        self.prefix:str
            String for the transcript of the variant. I.e. NM_40091.5
        self.suffix:str
            String for the remaining variant string for further processing.
        self.no_repeats:int
            The number of expanded repeats present in the variant.
        self.after_brac
            string for any remaining str after the final bracket
            in the variant_string.
        """
        print(self.variant_string)
        print(self.type)
        print(self.prefix)
        # Get g or c position
        variant_type = re.search(':(.*?)\.', self.suffix).group(1)
        print(f'Variant type: {variant_type}')
        var_position = re.search('\.(.*?)[ACTG]', self.suffix).group(1)
        print(f'Variant position: {var_position}')
        if "_" in var_position:
            start_range, end_range = var_position.split("_")
            print(start_range)
            print(end_range)
            rep_seq = re.search('\.[0-9]+_[0-9]+(.*?)\[', self.suffix).group(1)
            print(f'Repeated sequence: {rep_seq}')
        else:
            rep_seq = re.search('\.[0-9]+(.*?)\[', self.suffix).group(1)
            print(f'Repeat seq without range: {rep_seq}')
        # Get number of unit repeats
        self.no_repeats = re.search('\[(.*?)\]', self.suffix).group(1)
        print(f'Number of unit repeats: {self.no_repeats}')
        # Get anything after ] to check
        self.after_brac = re.search('\](.*)', self.suffix).group(1)
        print(f'Anything after bracket: {self.after_brac}')
        return self.variant_string, self.type, self.prefix, self.suffix, self.no_repeats, self.after_brac


# Hard-coded variant for testing while building.
# List of variants to check format and split into constituents.
VARIANT1 = "LRG_199:g.1ACT[20]"
VARIANT2 = "LRG_199:g.1ACT[20]A"
VARIANT3 = "LRG_199:g.1AC[20]"
VARIANT4 = "LRG_199t1:c.1_3ACT[20]"
VARIANT5 = "LRG_199t1:c.1AC[20]"
VARIANT6 = "LRG_199t1:c.1act[20]"
VARIANT7 = "LRG199c.1A[1_2]"
VARIANT8 = "LRG_199:g.13ACT[20]"
VARIANT9 = "LRG_199:g.13_25ACTG[1]"
VARIANT10 = "LRG199:g.13_125ACTG[1]"
# Other types not LRG
VARIANT11 = "ENSG00000198947.15:g.1ACT[10]"
VARIANT12 = "ENST00000357033.8:c.13AC[22]"
# Missing information accepted
VARIANT13 = "LRG_199t1:c.1_ACT[20]"
# change * to + to only allow variants with range or single location


def main():
    """
    Main function for testing the functions in the script in the
    ex_repeat_var class.
    """
    t_var = ex_repeat_var(VARIANT4, "GRch37")
    results = ex_repeat_var.check_expanded_repeat_diverging(t_var)
    print(t_var.variant_string)
    print(t_var.build)
    print(t_var.type)
    print(t_var.prefix)
    print(t_var.suffix)
    print(t_var.no_repeats)
    print(t_var.after_brac)
    # print(results)
    # print(t_var.get_transcript_type())
    # print(t_var1.split_var_string())
    # print(t_var1.get_transcript_name())
    # input = ex_repeat_var.get_transcript_name(variant11,
    # ex_repeat_var.get_transcript_type(variant4))
    # print(ex_repeat_var.get_variant_location(input[0], input[1], input[2], input[3]))


# This allows the script to be run by itself or imported as a package.
if __name__ == "__main__":
    log.info('--------- Starting Script ---------')
    main()
    log.info('--------- End Script ---------')
