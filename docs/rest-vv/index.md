# REST VariantValidator

The REST VariantValidator documentation describes the REST API interfaces available for integrating VariantValidator into external applications, automated pipelines and bioinformatics workflows.

VariantValidator currently provides two REST API implementations:

- [REST VariantValidator](rest_VariantValidator.md) — the original Flask-based REST API.
- [REST VariantValidator2 (SHAIP)](rest_VariantValidator2.md) — the current high-performance FastAPI implementation.

Both APIs expose VariantValidator functionality through HTTP endpoints, enabling automated validation, transcript mapping and variant formatting without requiring direct use of the Python API.

> **Note**
>
> The REST APIs are intended for software developers and automated workflows. If you wish to validate variants interactively through a web browser, see the [VariantValidator Web Interface (VVweb)](../vvweb/index.md).

---

# Choosing a REST API

## REST VariantValidator

REST VariantValidator is the original Flask-based REST implementation for VariantValidator.

It provides HTTP access to the VariantValidator framework and remains available to support existing deployments and workflows.

The repository is **not publicly accessible**. Access may be requested by contacting the VariantValidator development team through our [support page](https://variantvalidator.org/help/contact/).

See:

- [REST VariantValidator](rest_VariantValidator.md)

---

## REST VariantValidator2 (SHAIP)

REST VariantValidator2 (SHAIP) is the current generation REST API for VariantValidator.

Built using FastAPI and Python 3.12, REST VariantValidator2 provides improved performance, scalability and support for modern containerised deployments while maintaining compatibility with existing VariantValidator workflows.

The REST VariantValidator2 source repository is **not publicly accessible**.

Complete endpoint documentation, deployment guidance and installation instructions are provided through the dedicated REST VariantValidator2 documentation.

See:

- [REST VariantValidator2 (SHAIP)](rest_VariantValidator2.md)

---

# Common capabilities

Both REST APIs provide programmatic access to VariantValidator functionality, including:

- validation of HGVS variant descriptions;
- mapping between genomic, transcript and protein reference sequences;
- variant formatting;
- transcript selection;
- support for RefSeq and Ensembl transcript collections;
- support for HGVS, VCF and pseudo-VCF input formats; and
- high-throughput batch processing suitable for automated workflows.

The following reference documentation may also be useful when working with the REST APIs:

- [Supported Input Formats](../user-manual/reference/supported_inputs.md)
- [Transcript Selection](../user-manual/reference/transcript_selection.md)
- [Output Formats](../user-manual/reference/output_formats.md)
- [Errors and Error Codes](../user-manual/reference/errors_and_error_codes.md)

---

# Getting help

If you require assistance using either REST API, please contact the VariantValidator development team through our [support page](https://variantvalidator.org/help/contact/).

Access requests for the original REST VariantValidator repository may also be made through the support page.

---

# How to cite VariantValidator

If you use VariantValidator in your research, please [cite the appropriate VariantValidator publication(s)](https://github.com/openvar/VariantValidator#cite-us).

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="../static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="../static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>