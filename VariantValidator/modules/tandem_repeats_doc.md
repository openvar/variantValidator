# VariantValidator - tandem repeats Module

[![codecov](https://codecov.io/gh/openvar/variantValidator/branch/restructuring_py3/graph/badge.svg)](https://codecov.io/gh/openvar/variantValidator)

## About

VariantValidator is a user-friendly software tool designed to validate the syntax and
parameters of DNA variant descriptions according to the HGVS Sequence Variant
Nomenclature.

The aim of the tandem repeats module is to add support of tandem repeat variants from HGVS (see recommendations here: https://varnomen.hgvs.org/recommendations/DNA/variant/repeated/).

Our module add a syntax checker that takes the input for tandem repeat variants, checks the syntax and returns a HGVS formatted variant.


## Features

The basic functionality of https://variantvalidator.org/ and VarinantValidator is documented here https://www.ncbi.nlm.nih.gov/pubmed/28967166

This adds further functionality by handling tandem repeats.

Function: init() --
This take an initialises an instance of a variant with HGVS variant string and genome build (i.e. grch37).

Function: parse_repeat_variant() --
This take an instance of a variant and splits it into it's components if it is an expanded repeat, this is determined by the presence of either square brackets ([]).

## Test coverage for this module

<img src="https://github.com/i3hsInnovation/resources/blob/master/images/UoM_logo.jpg?raw=true">

## How to contribute
Please refer to [CONTRIBUTING.md](https://github.com/openvar/variantValidator/blob/master/CONTRIBUTING.md)

## Acknowledgements

**VariantValidator was developed at the University of Leicester. It is now maintained and developed by the University of Manchester and is hosted (with ongoing development contributions) by the University of Leicester**

<img src="https://github.com/i3hsInnovation/resources/blob/master/images/UoM_logo.jpg?raw=true" width="40%" align="left"/>
<img src="https://github.com/i3hsInnovation/resources/blob/master/images/UoL-Logo-Full-Colour.png?raw=true" width="40%" align="right" />
