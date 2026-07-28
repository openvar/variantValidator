<img src="../static/img/logos/VV_logo.png" width="20%" />

# REST VariantValidator

REST VariantValidator is the original REST API implementation for VariantValidator. It provides programmatic access to VariantValidator, VariantFormatter, Gene2Transcripts and associated tools through a REST interface, allowing integration into automated bioinformatics workflows, laboratory information management systems (LIMS) and external software applications.

The REST API exposes VariantValidator functionality through standard HTTP endpoints and uses the same VariantValidator validation and formatting infrastructure available through the other VariantValidator interfaces.

---

# Accessing REST VariantValidator

REST VariantValidator is provided as a hosted service.

The source repository is **not publicly accessible**. Access to the service may be requested by contacting the VariantValidator development team.

To request access:

1. Follow the account request process described in the REST VariantValidator access guide:

   https://github.com/openvar/rest_variantValidator/blob/master/docs/Account.md

2. Once your request has been submitted, email **admin@variantvalidator.org** so that your account can be reviewed and activated.

Once your account has been approved, you will receive the credentials required to access the REST VariantValidator service.

---

# Interactive API documentation

REST VariantValidator provides an interactive Swagger interface describing the available endpoints, required parameters and example requests.

The Swagger interface also reports the currently deployed software versions together with the transcript annotation and sequence repository releases used by the service.

---

# VariantValidator endpoints

The VariantValidator endpoints provide programmatic access to variant validation, transcript retrieval and reference sequence utilities.

Available endpoints include:

- `GET /VariantValidator/variantvalidator/{genome_build}/{variant_description}/{select_transcripts}`
- `GET /VariantValidator/variantvalidator_ensembl/{genome_build}/{variant_description}/{select_transcripts}`
- `GET /VariantValidator/tools/gene2transcripts_v2/{gene_query}/{limit_transcripts}/{transcript_set}/{genome_build}`
- `GET /VariantValidator/tools/hgvs2reference/{hgvs_description}`

These endpoints provide programmatic access to VariantValidator validation, transcript retrieval and reference sequence functionality.

See also:

- [Validator](../vvweb/validator.md)
- [Gene2Transcripts](../vvweb/gene2transcripts.md)
- [Supported Input Formats](../user-manual/reference/supported_inputs.md)
- [Transcript Selection](../user-manual/reference/transcript_selection.md)
- [Errors and Error Codes](../user-manual/reference/errors_and_error_codes.md)

---

# VariantFormatter endpoint

VariantFormatter is also available through the REST API.

Available endpoint:

- `GET /VariantFormatter/variantformatter/{genome_build}/{variant_description}/{transcript_model}/{select_transcripts}/{checkonly}`

This endpoint provides programmatic access to VariantFormatter, allowing automated formatting of genomic variant descriptions and mapping to corresponding transcript and protein representations.

VariantFormatter accepts the same supported genomic input formats as VariantValidator. However, VariantFormatter processing begins with a genomic variant and does not accept transcript, RNA or protein descriptions as starting inputs.

See also:

- [VariantFormatter CLI](../user-manual/cli/variantformatter_cli.md)
- [VariantFormatter Python API](../user-manual/python-api/variantformatter_python.md)
- [Supported Input Formats](../user-manual/reference/supported_inputs.md)
- [Transcript Selection](../user-manual/reference/transcript_selection.md)
- [Output Formats](../user-manual/reference/output_formats.md)

---

# LOVD endpoint

The LOVD endpoint provides a more granular interface to the VariantFormatter framework, exposing additional parameters for workflows requiring greater control over formatting behaviour.

Available endpoint:

- `GET /LOVD/lovd/{genome_build}/{variant_description}/{transcript_model}/{select_transcripts}/{checkonly}/{liftover}`

The endpoint is particularly useful for automated workflows that require explicit control over transcript models, validation behaviour and liftover operations.

As with VariantFormatter, processing begins with a genomic variant. Supported genomic inputs are interpreted by the underlying VariantValidator and VariantFormatter infrastructure before transcript and protein representations are generated.

See also:

- [VariantFormatter CLI](../user-manual/cli/variantformatter_cli.md)
- [VariantFormatter Python API](../user-manual/python-api/variantformatter_python.md)
- [Supported Input Formats](../user-manual/reference/supported_inputs.md)
- [Transcript Selection](../user-manual/reference/transcript_selection.md)
- [Output Formats](../user-manual/reference/output_formats.md)

---

# Service monitoring endpoints

The following endpoints are available for monitoring and testing the service:

- `GET /hello/`
- `GET /hello/limit`
- `GET /hello/trigger_error/{error_code}`

These endpoints allow administrators and automated workflows to verify that the service is operational and inspect information reported by the deployed service.

---

# Reference documentation

The following documentation may also be useful when developing applications that use REST VariantValidator:

- [Validator](../vvweb/validator.md)
- [Gene2Transcripts](../vvweb/gene2transcripts.md)
- [VariantFormatter CLI](../user-manual/cli/variantformatter_cli.md)
- [VariantFormatter Python API](../user-manual/python-api/variantformatter_python.md)
- [Supported Input Formats](../user-manual/reference/supported_inputs.md)
- [Transcript Selection](../user-manual/reference/transcript_selection.md)
- [Output Formats](../user-manual/reference/output_formats.md)
- [Errors and Error Codes](../user-manual/reference/errors_and_error_codes.md)

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="../static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="../static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>