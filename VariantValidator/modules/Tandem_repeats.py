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
            elif re.search(pattern=REGEX2, string=raw_s):
                print(f"{raw_s}-True 2")
            else:
                print("False")  # change to raise error


#  List of variants to check format and split into constituents.
variant1 = "LRG_199:g.1ACT[20]"
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


def check_transcript_type(my_variant):
    """
    [Find transcript type and run relevant function for processing,
    that transcript type. N.B. Future development could instead store
    the transcript and replace it with refseq.]

    Args:
        my_variant ([type]): [description]

    Returns: results from relevant chech_expanded_repeat_X
    where X is the transcript type.

    Raises:
        NameError: [Error for unknown transcript type.]
    """
    if bool(re.search("^LRG", my_variant)):
        check_expanded_repeat_LRG(my_variant)
    elif bool(re.search("^ENSG", my_variant)):
        check_expanded_repeat_ENSG(my_variant)
    elif bool(re.search("NM", my_variant)):
        check_expanded_repeat_RefSeq(my_variant)
    else:
        raise NameError('Unknown transcript type present. \
                        Try RefSeq transcript ID')


def check_expanded_repeat_RefSeq(my_variant):
    """
    Summary:
    This takes a variant string and breaks it into its constituents.
    This isolates the constituents with regex.

    Args:
        my_variant ([type]): (Variant string i.e. LRG_199:g.1ACT[20])
    Returns:
    Prints out constituents and assigns them to variables for further processing.
    """
    print("Running RefSeq version")
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
            suffix = ":" + my_variant.split(":")[1]
            # Find whether genomic or coding
            variant_type = re.search(':(.*?)\.', suffix)
            print(f'Variant type: {variant_type.group(1)}')
            # Get g or c position
            var_position = re.search('\.(.*?)[ACTG]', suffix)
            print(f'Variant position: {var_position.group(1)}')
            if "_" in var_position.group(1):
                start_range, end_range = var_position.group(1).split("_")
                print(start_range)
                print(end_range)
                rep_seq = re.search('\.[0-9]+_[0-9]+(.*?)\[', suffix)
                print(f'Repeated sequence: {rep_seq.group(1)}')
            else:
                rep_seq = re.search('\.[0-9]+(.*?)\[', suffix)
                print(f'Repeat seq without range: {rep_seq.group(1)}')
            # Get number of unit repeats
            no_repeats = re.search('\[(.*?)\]', suffix)
            print(f'Number of unit repeats: {no_repeats.group(1)}')
            # Get anything after ] to check
            after_brac = re.search('\](.*)', suffix)
            print(f'Anything after bracket: {after_brac.group(1)}')
    else:
        print("No expanded repeat present.")
        return False


def check_expanded_repeat_ENSG(my_variant):
    """
    Summary:
    This takes a variant string and breaks it into its constituents.
    This isolates the constituents with regex.

    Args:
        my_variant ([type]): (Variant string i.e. LRG_199:g.1ACT[20])
    Returns:
    Prints out constituents and assigns them to variables for further processing.
    """
    print("Running ENSG version")
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
            suffix = ":" + my_variant.split(":")[1]
            # Find whether genomic or coding
            variant_type = re.search(':(.*?)\.', suffix)
            print(f'Variant type: {variant_type.group(1)}')
            # Get g or c position
            var_position = re.search('\.(.*?)[ACTG]', suffix)
            print(f'Variant position: {var_position.group(1)}')
            if "_" in var_position.group(1):
                start_range, end_range = var_position.group(1).split("_")
                print(start_range)
                print(end_range)
                rep_seq = re.search('\.[0-9]+_[0-9]+(.*?)\[', suffix)
                print(f'Repeated sequence: {rep_seq.group(1)}')
            else:
                rep_seq = re.search('\.[0-9]+(.*?)\[', suffix)
                print(f'Repeat seq without range: {rep_seq.group(1)}')
            # Get number of unit repeats
            no_repeats = re.search('\[(.*?)\]', suffix)
            print(f'Number of unit repeats: {no_repeats.group(1)}')
            # Get anything after ] to check
            after_brac = re.search('\](.*)', suffix)
            print(f'Anything after bracket: {after_brac.group(1)}')
    else:
        print("No expanded repeat present.")
        return False


def check_expanded_repeat_LRG(my_variant):
    """
    Summary:
    This takes a variant string and breaks it into its constituents.
    This isolates the constituents with regex.

    Args:
        my_variant ([type]): (Variant string i.e. LRG_199:g.1ACT[20])
    Returns:
    Prints out constituents and assigns them to variables for further processing.
    """
    print("Running LRG version")
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
            suffix = ":" + my_variant.split(":")[1]
            # Find whether genomic or coding
            variant_type = re.search(':(.*?)\.', suffix)
            print(f'Variant type: {variant_type.group(1)}')
            # Get g or c position
            var_position = re.search('\.(.*?)[ACTG]', suffix)
            print(f'Variant position: {var_position.group(1)}')
            if "_" in var_position.group(1):
                start_range, end_range = var_position.group(1).split("_")
                print(start_range)
                print(end_range)
                rep_seq = re.search('\.[0-9]+_[0-9]+(.*?)\[', suffix)
                print(f'Repeated sequence: {rep_seq.group(1)}')
            else:
                rep_seq = re.search('\.[0-9]+(.*?)\[', suffix)
                print(f'Repeat seq without range: {rep_seq.group(1)}')
            # Get number of unit repeats
            no_repeats = re.search('\[(.*?)\]', suffix)
            print(f'Number of unit repeats: {no_repeats.group(1)}')
            # Get anything after ] to check
            after_brac = re.search('\](.*)', suffix)
            print(f'Anything after bracket: {after_brac.group(1)}')
    else:
        print("No expanded repeat present.")
        return False























# incorporate diverging based on transcript type

def check_expanded_repeat_diverging(my_variant):
    """
    Summary:
    This takes a variant string and breaks it into its constituents.
    This isolates the constituents with regex.

    Args:
        my_variant ([type]): (Variant string i.e. LRG_199:g.1ACT[20])
    Returns:
    Prints out constituents and assigns them to variables for further processing.
    """
    print("Running Diverging version")
    if "[" or "]" in my_variant:
        if bool(re.search("^LRG", my_variant)):
            type = "LRG"
        elif bool(re.search("^ENSG", my_variant)):
            type = "Ensembl"
        elif bool(re.search("^NM", my_variant)):
            type = "RefSeq"
        if ":" not in my_variant:
            print("Unable to identify a colon (:) in the variant description. A colon is required in HGVS variant")
        else:
            prefix = my_variant.split(":")[0]
            print(f'Variant prefix: {prefix}')
            if type=="LRG":
                # Check if underscore after LRG is included
                if "_" not in prefix:
                    # Add in underscore between LRG and number
                    prefix = re.sub(r"(?i)(?<=[a-z])(?=\d)",'_', prefix)
                    print(f'Updated prefix: {prefix}')
            suffix = ":" + my_variant.split(":")[1]
            # Find whether genomic or coding
            variant_type = re.search(':(.*?)\.', suffix)
            print(f'Variant type: {variant_type.group(1)}')
            # Get g or c position
            var_position = re.search('\.(.*?)[ACTG]', suffix)
            print(f'Variant position: {var_position.group(1)}')
            if "_" in var_position.group(1):
                start_range, end_range = var_position.group(1).split("_")
                print(start_range)
                print(end_range)
                rep_seq = re.search('\.[0-9]+_[0-9]+(.*?)\[', suffix)
                print(f'Repeated sequence: {rep_seq.group(1)}')
            else:
                rep_seq = re.search('\.[0-9]+(.*?)\[', suffix)
                print(f'Repeat seq without range: {rep_seq.group(1)}')
            # Get number of unit repeats
            no_repeats = re.search('\[(.*?)\]', suffix)
            print(f'Number of unit repeats: {no_repeats.group(1)}')
            # Get anything after ] to check
            after_brac = re.search('\](.*)', suffix)
            print(f'Anything after bracket: {after_brac.group(1)}')
    else:
        print("No expanded repeat present.")
        return False











## Try to break this into many different methods and define class

class Tandem_variant:

    def __init__(self, variant_string, build) -> None:
        self.variant_string = variant_string
        self.build = build
        self.type = None
        self.prefix = None
        self.suffix = None
        self.no_repeats = None
        self.after_brac = None # could use empty string?


    def check_expanded_repeat_diverging(self):
        """
        Summary:
        This takes a variant string and breaks it into its constituents.
        This isolates the constituents with regex.

        Args:
            self.variant_string ([type]): (Variant string i.e. LRG_199:g.1ACT[20])
        Returns:
        Prints out constituents and assigns them to variables for further processing.
        """
        print("Running Diverging version")
        if "[" or "]" in self.variant_string:
            if bool(re.search("^LRG", self.variant_string)):
                self.type = "LRG"
            elif bool(re.search("^ENSG", self.variant_string)):
                self.type = "Ensembl"
            elif bool(re.search("^NM", self.variant_string)):
                self.type = "RefSeq"
            if ":" not in self.variant_string:
                print("Unable to identify a colon (:) in the variant description. A colon is required in HGVS variant")
            else:
                self.prefix = self.variant_string.split(":")[0]
                print(f'Variant prefix: {self.prefix}')
                if self.type=="LRG":
                    # Check if underscore after LRG is included
                    if "_" not in self.prefix:
                        # Add in underscore between LRG and number
                        self.prefix = re.sub(r"(?i)(?<=[a-z])(?=\d)",'_', self.prefix)
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
                # Get anything after ] to check
                self.after_brac = re.search('\](.*)', self.suffix).group(1)
                print(f'Anything after bracket: {self.after_brac}')
        else:
            print("No expanded repeat present.")
            return False


    def get_transcript_type(self):
        """
        [Find transcript type and run relevant function for processing,
        that transcript type. N.B. Future development could instead store
        the transcript and replace it with refseq.]

        Args:
            self.variant_string ([type]): [description]

        Returns: results from relevant chech_expanded_repeat_X
        where X is the transcript type.

        Raises:
            NameError: [Error for unknown transcript type.]
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
        self.prefix = self.variant_string.split(":")[0]
        self.suffix = ":" + self.variant_string.split(":")[1]
        return self.prefix, self.suffix


    def get_transcript_name(self):
        self.prefix = self.variant_string.split(":")[0]
        self.suffix = ":" + self.variant_string.split(":")[1]
        print(f'Variant prefix: {self.prefix}')
        if self.type=="LRG":
            # Check if underscore after LRG is included
            if "_" not in self.prefix:
                # Add in underscore between LRG and number
                self.prefix = re.sub(r"(?i)(?<=[a-z])(?=\d)",'_', self.prefix)
                print(f'Updated prefix: {self.prefix}')
        print(f'Variant string suffix: {self.suffix}')
        return self.variant_string, self.type, self.prefix, self.suffix


    def get_variant_location(self):
        print(self.variant_string)
        print(self.type)
        print(self.prefix)
        # Get g or c position
        variant_type = re.search(':(.*?)\.', self.suffix)
        print(f'Variant type: {variant_type.group(1)}')
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


def main():
    t_var1 = Tandem_variant(variant4, "GRch37")
    print(t_var1.variant_string)
    print(t_var1)
    print(t_var1.get_transcript_type())
    #print(t_var1.split_var_string())
    print(t_var1.get_transcript_name())
    #input = Tandem_variant.get_transcript_name(variant11, Tandem_variant.get_transcript_type(variant4))
    #print(Tandem_variant.get_variant_location(input[0], input[1], input[2], input[3]))


if __name__ == "__main__":
    main()
