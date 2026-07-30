<img src="../static/img/logos/VV_logo.png" width="20%" />

# REST VariantValidator

REST VariantValidator provides programmatic access to VariantValidator, VariantFormatter, Gene2Transcripts and associated tools through a REST interface, allowing integration into automated bioinformatics workflows, laboratory information management systems (LIMS) and external software applications.

The REST API exposes VariantValidator functionality through standard HTTP endpoints and uses the same VariantValidator validation and formatting infrastructure available through the other VariantValidator interfaces.

---

# Legacy API implementation

The original REST VariantValidator implementation is **no longer maintained**.

It has been replaced by a new implementation:

**https://github.com/openvar/legacy_rest_variantValidator**

The new implementation preserves the existing REST VariantValidator interface, including the original endpoints, URL structure and request parameters, allowing existing applications to continue operating without modification.

Internally, however, the service no longer performs VariantValidator or VariantFormatter processing itself. Instead, it acts as a lightweight compatibility layer that translates legacy REST VariantValidator requests into requests for REST VariantValidator2 before returning responses in the format expected by existing clients.

For users of the service, the functionality and API remain unchanged.

---

# Accessing REST VariantValidator

The legacy REST VariantValidator service is available at:

**https://rest.variantvalidator.org**

The service requires authentication before requests can be made.

The source repository is **not publicly accessible**. Access to the hosted service may be requested by contacting the VariantValidator development team.

To request access:

1. Follow the account request process described in the REST VariantValidator access guide:

   https://github.com/openvar/legacy_rest_variantValidator/blob/master/docs/Account.md

2. Once your request has been submitted, email **admin@variantvalidator.org** so that your account can be reviewed and activated.

Once your account has been approved, you will receive the credentials required to access the REST VariantValidator service.

---

# Original REST VariantValidator repository

The original REST VariantValidator repository remains publicly available for reference and for users wishing to reproduce historic deployments:

**https://github.com/openvar/rest_variantValidator**

This repository is **no longer maintained** and has been archived in its final supported state.

The implementation is permanently tied to:

- **VariantValidator 4.0.0**
- **VariantFormatter 4.0.0**

It is **not compatible** with subsequent releases of either VariantValidator or VariantFormatter and will not receive updates to support newer software versions, transcript databases or reference data.

Users requiring compatibility with current or future releases should use the maintained compatibility implementation:

**https://github.com/openvar/legacy_rest_variantValidator**

---

# Interactive API documentation

REST VariantValidator provides an interactive Swagger interface describing the available endpoints, required parameters and example requests.

The Swagger interface also reports the currently deployed software versions together with the transcript annotation and sequence repository releases used by the service.

---

# Service limits

REST VariantValidator enforces a maximum request processing time of **300 seconds**.

Requests exceeding this limit will terminate with a timeout response. This timeout helps ensure fair resource usage and prevents exceptionally long-running requests from affecting the availability of the service for other users.

---

# VariantValidator endpoints

The VariantValidator endpoints provide programmatic access to variant validation, transcript retrieval and reference sequence utilities.

Available endpoints include:

- `GET /VariantValidator/variantvalidator/{genome_build}/{variant_description}/{select_transcripts}`
- `GET /VariantValidator/variantvalidator_ensembl/{genome_build}/{variant_description}/{select_transcripts}`
- `GET /VariantValidator/tools/gene2transcripts_v2/{gene_query}/{limit_transcripts}/{transcript_set}/{genome_build}`
- `GET /VariantValidator/tools/hgvs2reference/{hgvs_description}`

> **Note**
>
> The original `gene2transcripts` endpoint has been fully deprecated and is no longer available. As the number of transcripts within the VariantValidator database has grown, the original implementation became unsustainable in terms of performance and scalability. Users should migrate to the `gene2transcripts_v2` endpoint, which provides significantly improved performance and supports the current transcript database.

These endpoints provide programmatic access to VariantValidator validation, transcript retrieval and reference sequence functionality.

See also:

- [Validator](../vvweb/validator.md)
- [Gene2Transcripts](../vvweb/gene2transcripts.md)
- [Supported Input Formats](../user-manual/reference/supported_inputs.md)
- [Transcript Selection](../user-manual/reference/transcript_selection.md)
- [Errors and Error Codes](../user-manual/reference/errors_and_error_codes.md)

---

# VariantFormatter endpoint

VariantFormatter is available through the REST API.

Available endpoint:

- `GET /VariantFormatter/variantformatter/{genome_build}/{variant_description}/{transcript_model}/{select_transcripts}/{checkonly}`

Historically, the VariantFormatter endpoint exposed a simplified interface compared with the LOVD endpoint. The current implementation has been enhanced so that VariantFormatter now provides the same level of control over formatting behaviour as the LOVD endpoint while preserving complete backwards compatibility with existing clients.

VariantFormatter accepts the same supported genomic input formats as VariantValidator. Processing begins with a genomic variant and generates equivalent transcript and protein representations.

### Legacy response format

By default, the VariantFormatter endpoint returns responses using the legacy REST VariantValidator genomic response structure to preserve compatibility with existing client software.

This behaviour can be controlled using the optional `legacy_genomic_structure` parameter:

- `True` (default) — Return the legacy REST VariantValidator response structure.
- `False` — Return the native REST VariantValidator2 genomic response structure which is identical to the VariantValidator structure.

Existing applications do not need to specify this parameter, as the default behaviour preserves the original REST VariantValidator interface.

See also:

- [VariantFormatter CLI](../user-manual/cli/variantformatter_cli.md)
- [VariantFormatter Python API](../user-manual/python-api/variantformatter_python.md)
- [Supported Input Formats](../user-manual/reference/supported_inputs.md)
- [Transcript Selection](../user-manual/reference/transcript_selection.md)
- [Output Formats](../user-manual/reference/output_formats.md)

---

# LOVD endpoint

The LOVD endpoint remains available for backwards compatibility with existing applications.

Available endpoint:

- `GET /LOVD/lovd/{genome_build}/{variant_description}/{transcript_model}/{select_transcripts}/{checkonly}/{liftover}`

The LOVD endpoint now provides the same functionality as the VariantFormatter endpoint. Both endpoints are translated into equivalent REST VariantValidator2 requests and produce consistent results. Existing applications can therefore continue to use either interface without modification.

### Legacy response format

By default, the VariantFormatter endpoint returns responses using the legacy REST VariantValidator genomic response structure to preserve compatibility with existing client software.

This behaviour can be controlled using the optional `legacy_genomic_structure` parameter:

- `True` (default) — Return the legacy REST VariantValidator response structure.
- `False` — Return the native REST VariantValidator2 genomic response structure which is identical to the VariantValidator structure.

Existing applications do not need to specify this parameter, as the default behaviour preserves the original REST VariantValidator interface.

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