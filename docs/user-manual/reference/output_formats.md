<img src="../../static/img/logos/VV_logo.png" width="20%" />

# Output Formats

The tools within the VariantValidator software suite return results in several formats depending on the interface and operation being used.

This guide describes the output formats provided by:

- VariantValidator;
- VariantFormatter;
- Gene2Transcripts; and
- hgvs2reference.

VariantValidator and VariantFormatter return structured validation or formatting results that can be represented as Python dictionaries, JSON or tabular output.

Gene2Transcripts and hgvs2reference return structured data suitable for direct use within Python applications or conversion to JSON.

## See also

- [Supported Input Formats](supported_inputs.md) — Input formats accepted by VariantValidator and related tools.
- [Transcript Selection](transcript_selection.md) — Available transcript selection strategies.
- [Errors and Error Codes](errors_and_error_codes.md) — Validation errors, warnings and informational messages.
- [VariantValidator Python API](../python-api/variantvalidator_python.md) — Validate variants directly from Python.
- [VariantFormatter Python API](../python-api/variantformatter_python.md) — Format genomic variants directly from Python.
- [Gene2Transcripts Python API](../python-api/gene2transcripts_python.md) — Retrieve transcript information directly from Python.
- [hgvs2reference Python API](../python-api/hgvs2reference_python.md) — Retrieve reference sequence for HGVS variants.

---

# VariantValidator

The VariantValidator `validate()` method returns a `ValOutput` object containing the results of validation.

For example:

```python
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select"
)
```

The returned `ValOutput` object provides methods for representing the validation results in different formats.

---

## Python Dictionary

Use:

```python
validation.format_as_dict()
```

to return the validation results as a Python dictionary.

For example:

```python
results = validation.format_as_dict()
```

This format is particularly useful when VariantValidator is being integrated directly into another Python application because individual fields can be accessed without serialising and parsing the results.

Metadata can be included using:

```python
results = validation.format_as_dict(with_meta=True)
```

---

## JSON

Use:

```python
validation.format_as_json()
```

to return the validation results as JSON.

For example:

```python
results = validation.format_as_json()
print(results)
```

Metadata can be included using:

```python
results = validation.format_as_json(with_meta=True)
```

JSON output is useful for:

- application integration;
- data exchange;
- storage of validation results; and
- downstream processing by software written in languages other than Python.

---

## Tabular Output

Use:

```python
validation.format_as_table()
```

to return the validation results in tabular form.

For example:

```python
results = validation.format_as_table()
print(results)
```

Metadata can be included using:

```python
results = validation.format_as_table(with_meta=True)
```

Tabular output is particularly useful for:

- batch validation;
- spreadsheet import;
- downstream command-line processing; and
- manual review of larger collections of variants.

---

# Validation Metadata

VariantValidator can include metadata describing the environment in which validation was performed.

For the Python API, metadata can be requested using the `with_meta` argument:

```python
validation.format_as_dict(with_meta=True)
```

```python
validation.format_as_json(with_meta=True)
```

```python
validation.format_as_table(with_meta=True)
```

The metadata provides information about the validation environment and may include details such as:

- the VariantValidator software version;
- database and reference data versions;
- transcript selection information; and
- other information required to interpret or reproduce the validation results.

Where validation results are being archived, shared or used as part of a reproducible analysis, retaining the associated metadata is recommended.

---

# Validation Messages

VariantValidator results may contain validation messages describing errors, warnings, automatic corrections or other information generated during validation.

These messages use standardised identifiers such as:

```text
ReferenceMismatchError
```

```text
TranscriptVersionWarning
```

```text
VariantNormalizationWarning
```

Applications processing VariantValidator output programmatically should use these stable identifiers when classifying validation outcomes rather than depending on the complete human-readable message text.

See [Errors and Error Codes](errors_and_error_codes.md) for details.

---

# VariantFormatter

VariantFormatter converts genomic variants into transcript and protein representations and returns structured formatting results.

The precise representation depends on the interface being used, but the output is designed for programmatic processing and contains the genomic variant together with its corresponding transcript-level representations.

VariantFormatter output may include:

- the validated genomic HGVS description;
- transcript HGVS descriptions;
- predicted protein consequences where applicable;
- transcript identifiers;
- gene information;
- genomic mappings; and
- warnings or other information generated during formatting.

The transcripts included in the output are controlled by the selected transcript selection strategy.

See [Transcript Selection](transcript_selection.md) for details.

---

# Gene2Transcripts

Gene2Transcripts returns structured transcript information for a supplied:

- gene symbol;
- HGNC identifier;
- RefSeq transcript accession; or
- Ensembl transcript accession.

Results are returned as structured data and can be represented as JSON.

For example, the command-line interface writes JSON directly:

```bash
gene2transcripts \
    --query COL1A1
```

Results can also be written to a file:

```bash
gene2transcripts \
    --query COL1A1 \
    --output results.json
```

When Gene2Transcripts is used through the Python API, the returned Python object can be processed directly by the calling application.

---

## Gene2Transcripts Content

Depending on the supplied query and selected options, Gene2Transcripts output may contain information describing:

- the submitted query;
- gene symbol;
- HGNC identifier;
- transcript accessions;
- transcript versions;
- transcript selection information;
- genomic spans;
- transcript annotations; and
- associated reference sequence information.

The precise transcript set returned depends on the selected transcript model and transcript selection strategy.

See [Transcript Selection](transcript_selection.md) for details.

---

# hgvs2reference

`hgvs2reference` returns the reference sequence corresponding to a submitted genomic (`g.`) or coding DNA (`c.`) HGVS sequence variant description.

For example:

```python
from VariantValidator import Validator

vv = Validator()

result = vv.hgvs2ref(
    "NC_000017.11:g.50198002C>A"
)
```

The method returns a Python dictionary.

---

## hgvs2reference Fields

The returned dictionary contains the following fields.

| Field | Description |
| --- | --- |
| `variant` | The submitted HGVS sequence variant description. |
| `start_position` | The resolved start position. |
| `end_position` | The resolved end position. |
| `sequence` | The reference sequence corresponding to the resolved coordinates. |
| `warning` | Any non-fatal warning generated during processing. |
| `error` | An error message if the reference sequence could not be retrieved. |

For example:

```python
result = vv.hgvs2ref(
    "NC_000017.11:g.50198002C>A"
)

print(result["sequence"])
```

Intronic coding DNA variants are supported when the genomic reference sequence used for the transcript alignment is explicitly supplied using compound HGVS notation.

For example:

```text
NC_000017.11(NM_000088.4):c.589+1G>T
```

See the [hgvs2reference Python API](../python-api/hgvs2reference_python.md) for details.

---

# Choosing an Output Format

The appropriate output format depends on how the results will be used.

| Output | Recommended use |
| --- | --- |
| Python dictionary | Direct integration into Python applications and pipelines. |
| JSON | Data exchange, APIs, storage and language-independent processing. |
| Tabular output | Batch analysis, spreadsheets and command-line processing. |

When working directly in Python, dictionaries generally provide the most convenient representation because individual fields can be accessed directly.

JSON is preferable when results need to be transferred between applications or stored in a portable structured format.

Tabular output is useful when results need to be reviewed manually or processed using spreadsheet and command-line tools.

---

# Batch Validator Output

The web Batch Validator returns results as a plain text, tab-delimited file.

The file is designed to be both human-readable and suitable for import into spreadsheet software such as Microsoft Excel or LibreOffice Calc.

Validation messages are reported in the **Warnings** column and may contain:

- validation errors;
- warnings;
- informational messages; and
- descriptions of automatic corrections performed during validation.

The Batch Validator also includes metadata describing the validation environment.

Users should retain this metadata when sharing or archiving validation results.

See the [Batch Validator](../../vvweb/batch_validator.md) documentation for details of using and interpreting Batch Validator output.

---

# Programmatic Processing

When VariantValidator output is being consumed by software, structured output should be preferred over parsing human-readable text.

In particular:

- use dictionary output when working directly in Python;
- use JSON when exchanging results between applications;
- use standard error and warning identifiers when interpreting validation messages; and
- retain metadata where reproducibility is important.

Human-readable warning and error descriptions may be expanded or clarified over time. Software should therefore avoid depending on exact message wording where a standardised error or warning identifier is available.

---

# Getting help

VariantValidator has been developed to support a wide range of users, from those new to HGVS nomenclature to experienced clinical scientists and bioinformaticians. If you encounter difficulties interpreting VariantValidator output, we encourage you to review the relevant documentation before contacting the development team.

You may find the following documentation helpful:

- [Supported Input Formats](supported_inputs.md)
- [Transcript Selection](transcript_selection.md)
- [Errors and Error Codes](errors_and_error_codes.md)
- [VariantValidator Python API](../python-api/variantvalidator_python.md)
- [VariantFormatter Python API](../python-api/variantformatter_python.md)
- [Gene2Transcripts Python API](../python-api/gene2transcripts_python.md)
- [hgvs2reference Python API](../python-api/hgvs2reference_python.md)

If you still require assistance, you can contact the VariantValidator team using our [contact form](https://variantvalidator.org/help/contact/).

Software bugs and feature requests can be reported through the [VariantValidator GitHub issue tracker](https://github.com/openvar/VariantValidator/issues).

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="../../static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="../../static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>
