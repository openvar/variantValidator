[<img src="/static/img/logos/VV_logo.png?raw=true" width="20%" />](https://variantvalidator.org/)

# Validator

The Validator is the primary tool within the VariantValidator Web Interface (VVweb). It validates individual sequence variant descriptions, identifies formatting or nomenclature issues, maps variants between reference sequences, and returns HGVS-compliant variant descriptions together with supporting annotation.

The Validator accepts all VariantValidator supported input formats. For a complete description of accepted formats, see [Supported Input Formats](../user-manual/reference/supported_inputs.md).

---

# Opening the Validator

The Validator can be accessed from the VVweb home page by:

- selecting **Try it out** on the **Validator** card; 

<img src="/static/img/Validator.png" width="40%" align="left"/>
<br clear="both"/>

or

- selecting **Tools → Validator** from the navigation bar.

---

# Entering a variant

Enter a single variant description into the **Variant Description** field.

The Validator accepts a wide range of HGVS-compliant and commonly encountered non-HGVS variant description formats. Where possible, VariantValidator will automatically recognise common formatting issues and convert them into HGVS-compliant variant descriptions while reporting any corrections that have been made.

If you are unsure whether your variant description is supported, see [Supported Input Formats](../user-manual/reference/supported_inputs.md).

---

## Selecting transcripts

The Validator provides two methods for selecting transcripts.

### Transcript selection options

The **Select transcript options** menu provides predefined transcript selection strategies, including commonly used transcript collections such as RefSeq, Ensembl, MANE Select and MANE Plus Clinical.

These options allow VariantValidator to automatically determine which transcript(s) should be used during validation.

For guidance on selecting the most appropriate transcript strategy for your application, see [Transcript Selection](../user-manual/reference/transcript_selection.md).

### Specifying transcript accessions

Alternatively, one or more transcript accessions may be entered directly into the transcript input field.

Multiple transcript accessions should be separated using the `|` character.

This option is useful when validation should be restricted to one or more specific transcript reference sequences rather than using one of the predefined transcript selection strategies.

For guidance on selecting transcript reference sequences, see [Transcript Selection](../user-manual/reference/transcript_selection.md).

---

## Selecting the genome build

Choose the genome assembly appropriate for your variant:

- **GRCh38 (hg38)**
- **GRCh37 (hg19)**

This option is required when validating genomic variant descriptions where the genome assembly cannot be determined directly from the submitted variant description.

Once all options have been selected, click **Submit**.

*Insert screenshot showing a completed submission form.*

---

# Understanding the results

The Validator returns a comprehensive report describing the validated variant together with equivalent variant representations and supporting annotation.

*Insert screenshot of the validation results.*

If the submitted variant maps to multiple transcripts, a transcript selection table is displayed near the top of the report.

Selecting a transcript from this table immediately adds a new results card for that transcript without requiring the variant to be revalidated.

The transcript table also provides useful information about each transcript, including whether it is:

- the latest transcript version;
- the MANE Select transcript; or
- a MANE Plus Clinical transcript, where applicable.

This allows you to compare transcript-specific representations of the same variant and select the transcript most appropriate for your application.

---

# Understanding the validation report

The validation report is organised into several sections describing different aspects of the validated variant.

Depending on the submitted variant and available annotation, these sections may include:

- HGVS-compliant variant descriptions;
- genomic variant representations;
- recommended variant descriptions;
- transcript and protein descriptions;
- exon and intron positions;
- gene information;
- projection of genomic variants onto alternative transcripts; and
- links to external genomic resources.

Many entries within the report include direct links to external reference databases, allowing rapid access to supporting information.

The structure and contents of the report are described in more detail in [Output Formats](../user-manual/reference/output_formats.md).

---

## Validation messages

During validation, VariantValidator may report errors, warnings or informational messages.

These messages explain any problems encountered during validation, highlight automatic corrections, or provide additional guidance on the submitted variant description.

For a complete description of all validation messages and recommended actions, see [Errors and Error Codes](../user-manual/reference/errors_and_error_codes.md).

---

# Printing a validation report

Selecting **Print PDF** generates a simplified PDF version of the validation report.

The printable report is designed for inclusion in:

- publication appendices;
- clinical reports;
- laboratory documentation;
- validation records; and
- supporting information accompanying manuscripts.

The PDF presents the essential validation results in a clean, printer-friendly format while preserving the key information required for reporting and record keeping.

---

# Further Reading

The following documentation may also be useful:

- [Supported Input Formats](../user-manual/reference/supported_inputs.md)
- [Transcript Selection](../user-manual/reference/transcript_selection.md)
- [Output Formats](../user-manual/reference/output_formats.md)
- [Errors and Error Codes](../user-manual/reference/errors_and_error_codes.md)

---

# How to cite VariantValidator

If VariantValidator contributes to your research or publication, please cite the appropriate VariantValidator publication(s).

The latest citation information is available in the project repository:

https://github.com/openvar/VariantValidator#cite-us

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="/static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="/static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>