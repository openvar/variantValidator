# User Manual

Welcome to the VariantValidator User Manual.

This manual describes how to use the VariantValidator software suite. It is intended for researchers, bioinformaticians, clinical scientists and software developers who wish to validate, normalise, map and format genetic sequence variant descriptions using the VariantValidator software suite.

The documentation is organised into task-focused sections, allowing users to quickly find information relevant to their workflow.

---

## Getting Started

If you have not yet installed VariantValidator, begin with the [Installation Guide](../installation/index.md).

---

## Choose an Interface

VariantValidator can be accessed through several interfaces, allowing you to choose the workflow that best suits your requirements.

### Local Installation

A local installation provides access to the complete VariantValidator software suite, including:

- VariantValidator
- VariantFormatter
- gene2transcripts
- hgvs2reference
- Python API
- Command-line interface (CLI)

Installation instructions are provided in the [Installation Guide](../installation/index.md).

### Web Interface

The VariantValidator web interface provides interactive validation and formatting of sequence variant descriptions directly from your web browser.

Access the web interface at:

- https://variantvalidator.org/

Documentation for the web interface will be added in a future release.

### REST API

The VariantValidator REST API provides programmatic access to the VariantValidator software suite over HTTP and is intended for software integration, automated workflows and high-throughput pipelines.

Access the REST API at:

- https://rest.variantvalidator.org/

Documentation for the REST API will be added in a future release.

---

## Using VariantValidator

VariantValidator validates, normalises and maps genetic sequence variant descriptions.

Topics include:

- [VariantValidator CLI](cli/variantvalidator_cli.md)
- [VariantValidator Python API](python-api/variantvalidator_python.md)
- [Supported Input Formats](reference/supported_inputs.md)
- [Output Formats](reference/output_formats.md)
- [Transcript Selection](reference/transcript_selection.md)
- [Errors and Error Codes](reference/errors_and_error_codes.md)

---

## Using VariantFormatter

VariantFormatter formats and converts genetic sequence variant descriptions between supported representations.

Topics include:

- [VariantFormatter CLI](cli/variantformatter_cli.md)
- [VariantFormatter Python API](python-api/variantformatter_python.md)
- [Supported Input Formats](reference/supported_inputs.md)
- [Output Formats](reference/output_formats.md)
- [Transcript Selection](reference/transcript_selection.md)

---

## Additional Tools

VariantValidator includes several utility tools that support common workflows.

### gene2transcripts

The `gene2transcripts` tool retrieves transcript information for genes.

Topics include:

- [gene2transcripts CLI](cli/gene2transcripts_cli.md)
- [gene2transcripts Python API](python-api/gene2transcripts_python.md)
- [Transcript Selection](reference/transcript_selection.md)

---

### hgvs2reference

The `hgvs2reference` tool retrieves reference sequences corresponding to supported HGVS sequence variant descriptions.

Topics include:

- [hgvs2reference Python API](python-api/hgvs2reference_python.md)
- [Supported Input Formats](reference/supported_inputs.md)

---

## Reference Documentation

The reference documentation contains information shared across the VariantValidator software suite.

Topics include:

- [Supported Input Formats](reference/supported_inputs.md)
- [Output Formats](reference/output_formats.md)
- [Transcript Selection](reference/transcript_selection.md)
- [Errors and Error Codes](reference/errors_and_error_codes.md)

---

## Need Help?

If you encounter problems:

- Verify that VariantValidator has been installed and configured correctly by following the [Installation Guide](../installation/index.md).
- Review the documentation for the relevant tool.
- Search the project's GitHub Issues for known problems.
- Open a GitHub Issue to report a bug or request a feature.
- Contact the VariantValidator support team via the support page.

---

## Documentation Structure

| Section | Description |
|----------|-------------|
| [Installation Guide](../installation/index.md) | Installing, configuring and verifying VariantValidator. |
| Web Interface | Interactive validation and formatting using the hosted VariantValidator website. Documentation coming soon. |
| REST API | Programmatic access to VariantValidator over HTTP. Documentation coming soon. |
| [VariantValidator CLI](cli/variantvalidator_cli.md) | Validate, normalise and map sequence variant descriptions from the command line. |
| [VariantValidator Python API](python-api/variantvalidator_python.md) | Validate, normalise and map sequence variant descriptions from Python. |
| [VariantFormatter CLI](cli/variantformatter_cli.md) | Format and convert supported sequence variant descriptions from the command line. |
| [VariantFormatter Python API](python-api/variantformatter_python.md) | Format and convert supported sequence variant descriptions from Python. |
| [gene2transcripts CLI](cli/gene2transcripts_cli.md) | Retrieve transcript information for genes from the command line. |
| [gene2transcripts Python API](python-api/gene2transcripts_python.md) | Retrieve transcript information for genes from Python. |
| [hgvs2reference Python API](python-api/hgvs2reference_python.md) | Retrieve reference sequences from supported HGVS sequence variant descriptions. |
| [Supported Input Formats](reference/supported_inputs.md) | Supported HGVS and non-HGVS input formats accepted throughout the software suite. |
| [Output Formats](reference/output_formats.md) | Output formats returned by the VariantValidator software suite. |
| [Transcript Selection](reference/transcript_selection.md) | Transcript selection strategies and transcript models. |
| [Errors and Error Codes](reference/errors_and_error_codes.md) | Reference guide to VariantValidator errors, warnings and diagnostic messages. |

---

As the VariantValidator software suite continues to evolve, this documentation will be updated to reflect new features, interfaces and best practices. The aim is to provide clear, practical guidance for both new and experienced users while supporting reproducible and reliable variant analysis workflows.