# VariantValidator - tandem repeats Module

[![codecov](https://codecov.io/gh/openvar/variantValidator/branch/restructuring_py3/graph/badge.svg)](https://codecov.io/gh/openvar/variantValidator)

## About

VariantValidator is a user-friendly software tool designed to validate the syntax and
parameters of DNA variant descriptions according to the HGVS Sequence Variant
Nomenclature.

The aim of the tandem repeats module is to add support of tandem repeat variants from HGVS (see recommendations here: https://varnomen.hgvs.org/recommendations/DNA/variant/repeated/).

Our module add a syntax checker that takes the input for tandem repeat variants, checks the syntax and returns a HGVS formatted variant.


## Features

The basic functionality of https://variantvalidator.org/ and VariantValidator is documented here https://www.ncbi.nlm.nih.gov/pubmed/28967166

This adds further functionality by checking the syntax and reformatting tandem repeat variants. 

Function: 

```init()```

This initialises an instance of a variant with a variant string, genome build (i.e. GRCh37) and selected transcripts.

Function: 

```parse_repeat_variant()```

This is a class method which takes a provided variant, genome build (e.g. GRCh37) and selected transcripts and splits the variant into its components if it is an expanded repeat; this is determined by the presence of either square brackets ([]). The components are then assigned to class attributes so that they can be used in downstream functions. 

Function:

```check_transcript_type()```

This checks the reference provided and assigns it to the attribute ref_type.

Function:

```reformat_reference()```

This checks for common mistakes in the reference provided.

Function:

```check_genomic_or_coding()```

This checks that the reference type provided is consistent with the prefix.

Function:

```check_positions_given()```

If a range is given for the position, this is checked with this function compared to the length of the repeat sequence given and the number of repeat units.

Function:

```get_range_from_single_pos()```

This function converts a single start position into the full range of a tandem repeat. 

Function:

```reformat_reference()```

This reformats the reference part of the variant to address some simple mistakes.

Function:

```reformat_not_multiple_of_three()```

HGVS nomenclature states that if a variant is a tandem repeat, is coding and is not a multiple of three, it should be described as a duplication (if length 1) or an insertion (length 2 or more). This function reformats those variants.

Function:

```reformat()```

This function does the final formatting by calling other functions and returns the final formatted variant string. If this is a dup or ins, the variant can be put back into the main VariantValidator function. 

## Test coverage for this module

<img src="https://user-images.githubusercontent.com/30113563/154497332-514419a6-1ab2-4492-829f-3286be2db45f.png">

## How to contribute
Please refer to [CONTRIBUTING.md](https://github.com/openvar/variantValidator/blob/master/CONTRIBUTING.md)

## Acknowledgements

**VariantValidator was developed at the University of Leicester. It is now maintained and developed by the University of Manchester and is hosted (with ongoing development contributions) by the University of Leicester**

<img src="https://github.com/i3hsInnovation/resources/blob/master/images/UoM_logo.jpg?raw=true" width="40%" align="left"/>
<img src="https://github.com/i3hsInnovation/resources/blob/master/images/UoL-Logo-Full-Colour.png?raw=true" width="40%" align="right" />
