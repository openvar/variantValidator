<img src="../../static/img/logos/VV_logo.png" width="20%" />

# hgvs2reference Python API

The `hgvs2ref()` method retrieves the reference sequence corresponding to an HGVS sequence variant description.

It is intended for applications that require the reference sequence underlying a submitted variant, for example to extract the reference allele or surrounding sequence.

## See also

- [VariantValidator Python API](variantvalidator_python.md) — Validate variant descriptions directly from Python.
- [Supported Input Formats](../reference/supported_inputs.md) — Supported HGVS and other variant description formats.
- [Errors and Error Codes](../reference/errors_and_error_codes.md) — Error messages and troubleshooting guidance.

---

# Importing VariantValidator

Import the `Validator` class:

```python
from VariantValidator import Validator
```

Create a `Validator` object:

```python
vv = Validator()
```

The `Validator` object can be reused for multiple queries within the same Python session.

---

# Basic Usage

Call `hgvs2ref()` with an HGVS sequence variant description:

```python
result = vv.hgvs2ref("NM_000088.4:c.589G>T")
```

The method returns a Python dictionary containing the retrieved reference sequence and associated metadata.

---

# Return Value

The returned dictionary contains the following fields.

| Field | Description |
| --- | --- |
| `variant` | The submitted HGVS sequence variant description. |
| `start_position` | The HGVS start position. |
| `end_position` | The HGVS end position. |
| `sequence` | The reference sequence corresponding to the variant coordinates. |
| `warning` | Any non-fatal warning generated during processing. |
| `error` | Error message if sequence retrieval failed. |

---

# Supported Input

`hgvs2ref()` accepts the HGVS sequence variant types described in the [Supported Input Formats](../reference/supported_inputs.md) guide.

RNA (`r.`), protein (`p.`) and mitochondrial (`m.`) sequence variants are not currently supported.

Fully intronic transcript variants cannot currently be resolved because HGVS transcript descriptions do not explicitly define the genomic reference sequence used for transcript alignment. A future HGVS nomenclature update is expected to address this limitation.

---

# Examples

## Retrieve the reference sequence for a transcript variant

```python
from VariantValidator import Validator

vv = Validator()

result = vv.hgvs2ref("NM_000088.4:c.589G>T")

print(result["sequence"])
```

---

## Retrieve the reference sequence for a genomic variant

```python
from VariantValidator import Validator

vv = Validator()

result = vv.hgvs2ref("NC_000017.11:g.50198002C>A")

print(result["sequence"])
```

---

## Inspect the complete result

The complete result can be inspected directly:

```python
from VariantValidator import Validator

vv = Validator()

result = vv.hgvs2ref("NC_000017.11:g.50198002C>A")

print(result)
```

The returned dictionary contains the sequence together with the variant coordinates and any warnings or errors generated during processing.

---

# Errors and Warnings

If the submitted variant cannot be parsed or the reference sequence cannot be retrieved, an error message is returned in the `error` field.

Warnings are returned in the `warning` field where sequence retrieval is only partially possible, such as transcript variants spanning intron boundaries.

For additional information about VariantValidator errors and warnings, see the [Errors and Error Codes](../reference/errors_and_error_codes.md) guide.

---

# Getting help

VariantValidator has been developed to support a wide range of users, from those new to HGVS nomenclature to experienced clinical scientists and bioinformaticians. If you encounter difficulties using the `hgvs2ref()` method or interpreting the returned reference sequence, we encourage you to seek assistance.

Before contacting the development team, you may find the following documentation helpful:

- [VariantValidator Python API](variantvalidator_python.md)
- [Supported Input Formats](../reference/supported_inputs.md)
- [Errors and Error Codes](../reference/errors_and_error_codes.md)

If you still require assistance, you can contact the VariantValidator team using our [contact form](https://variantvalidator.org/help/contact/).

Software bugs and feature requests can be reported through the [VariantValidator GitHub issue tracker](https://github.com/openvar/VariantValidator/issues).

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="../../static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="../../static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>
