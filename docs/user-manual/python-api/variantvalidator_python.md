# VariantValidator Python API

The VariantValidator Python API provides direct access to the VariantValidator validation engine from within Python. It 
is suitable for integrating variant validation into bioinformatics pipelines, analysis workflows, web applications and custom software.

The Python API offers access to the same validation functionality as the VariantValidator command-line interface while 
providing a programmatic interface for automated processing and downstream analysis.

For users who prefer not to write Python code:

- The [VariantValidator website](https://variantvalidator.org) provides a user-friendly alternative for validating variant descriptions.
- The [VariantValidator REST API](https://rest.variantvalidator.org) allows programmatic access to the validation services without requiring local installation.
- The [VariantValidator Command Line Interface](../cli/variantvalidator_cli.md) provides a command-line interface for validating variants locally.

---

## Basic Usage

Begin by importing the VariantValidator package and creating a `Validator` object.

```python
import json

import VariantValidator

vval = VariantValidator.Validator()
```

The `Validator` object manages access to the VariantValidator validation engine and can be reused to validate multiple 
variants within the same Python session.

Once a `Validator` object has been created, variants can be validated using the `validate()` method.

```python
import json

import VariantValidator

vval = VariantValidator.Validator()

variant = "NC_000017.11:g.50198002C>A"
genome_build = "GRCh38"

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select"
)

results = validation.format_as_dict(with_meta=True)

print(json.dumps(results, indent=4, sort_keys=True))
```

This validates the supplied variant and returns the results as a Python dictionary, which can then be processed directly 
or converted to JSON.

---

## Method Signature

Variant validation is performed using the `validate()` method of the `Validator` object.

```python
validate(
    variant=None,
    genome=None,
    select_transcripts="all",
    transcript_set=None,
    liftover_level=True,
    lovd_syntax_check=False,
    shorthand_vcf=False
)
```

The method returns a `ValOutput` object containing the validation results. These results can be accessed in a number of 
formats, including Python dictionaries and JSON.

---

## Required Arguments

The following arguments are required when calling the `validate()` method.

| Argument | Description |
|----------|-------------|
| `variant` | A single variant, multiple variants as a JSON array, or a filename containing variants to validate. |
| `genome` | The reference genome assembly (e.g. `GRCh37` or `GRCh38`). |

---

## Optional Arguments

The following optional arguments control validation behaviour.

| Argument | Default | Description                                                                          |
|----------|---------|--------------------------------------------------------------------------------------|
| `transcripts` | `"all"` | Transcript selection strategy or a JSON array of transcript identifiers.             |
| `transcript_set` | `None`  | Select the transcript database (`refseq` or `ensembl`). `None` defaults to `refseq`. |
| `liftover_level` | `True`  | Enable additional liftover to available genome builds.                               |
| `lovd_syntax_check` | `False` | Enable LOVD HGVS syntax checking.                                                    |
| `shorthand_vcf` | `False` | Enable shorthand VCF parsing.                                                        |


Usage of liftover_level

| Parameter | Type            | Required | Description |
|----------|-----------------|----------|-------------|
| liftover_level | string or bool  | No | Controls genomic liftover. `True` performs full liftover, `primary` excludes alternative scaffolds, and `False` disables liftover. Defaults to `True`. |

---

## Default Behaviour

Unless otherwise specified, the VariantValidator Python API uses the following defaults.

| Setting | Default |
|---------|---------|
| Genome assembly | User supplied (required) |
| Transcript selection | `"all"` |
| Transcript database | `refseq` |
| Liftover | Disabled |
| LOVD syntax checker | Disabled |
| Shorthand VCF parsing | Disabled |
| Return type | `ValOutput` object |
| Metadata | Disabled |

The `validate()` method returns a `ValOutput` object containing the validation results.

For example:

```python
validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select"
)
```

The `ValOutput` object provides several methods for accessing the validation results in different formats.

The optional `with_meta` argument controls whether validation metadata is included in the output. Setting `with_meta=True` 
includes additional information such as the software version, database versions, transcript selection strategy and other 
validation metadata. Setting `with_meta=False` returns only the validation results.

Return the results as a Python dictionary:

```python
results = validation.format_as_dict(with_meta=True)
```

Return the results as formatted JSON:

```python
results = validation.format_as_json(with_meta=True)
```

Return the results as tabular output:

```python
results = validation.format_as_table(with_meta=True)
```

A detailed description of each output format, including the available formatting options, is provided in the
[Output Formats](../reference/output_formats.md) guide.

---

## Supported Input Formats

The VariantValidator Python API accepts the same input formats as the command-line interface.

Supported variant descriptions include:

- Genomic HGVS (g. notation)
- Coding HGVS (c. notation)
- Non-coding HGVS (n. notation)
- RNA HGVS (r. notation)
- Protein HGVS (p. notation, including single-letter and three-letter amino acid codes - basic validation, not recommended)
- Pseudo-VCF/Chromosome coordinate notation (e.g. `17-50198002-C-A` or `17:50198002:C:A`)
- VCF notation (i.e. full VCF lines with chromosome, position, reference and alternate alleles)

The `variant` argument accepts:

- A single variant description.
- Multiple variant descriptions supplied as a JSON array.
- A text file containing one variant description per line.

See the [Supported Input Formats](../reference/supported_inputs.md) guide for a complete description of supported input formats and examples.

---

## Transcript Selection

VariantValidator supports multiple transcript selection strategies.

These include:

- MANE Select transcripts
- MANE Select and Plus Clinical transcripts
- All transcripts overlapping a genomic variant at their latest version
- All transcripts overlapping a genomic variant at all versions
- User-specified transcript lists

See the [Transcript Selection](../reference/transcript_selection.md) guide for details.

---

## Examples

### Validate a genomic variant

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select"
)

print(json.dumps(validation.format_as_dict(with_meta=True), indent=4, sort_keys=True))
```

---

### Validate a genomic variant using the Ensembl transcript set

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select",
    transcript_set="ensembl"
)

print(json.dumps(validation.format_as_dict(with_meta=True), indent=4, sort_keys=True))
```

---

### Validate a transcript variant

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NM_000088.4:c.589G>T",
    genome="GRCh38",
    select_transcripts="mane_select"
)

print(json.dumps(validation.format_as_dict(with_meta=True), indent=4, sort_keys=True))
```

---

### Validate an Ensembl transcript variant

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="ENST00000225964.10:c.589G>T",
    genome="GRCh38",
    select_transcripts="mane_select",
    transcript_set="ensembl"
)

print(json.dumps(validation.format_as_dict(with_meta=True), indent=4, sort_keys=True))
```

---

### Validate multiple variants

VariantValidator accepts multiple variants using a JSON array.

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant='["NC_000017.11:g.50198002C>A","NM_000088.4:c.589G>T"]',
    genome="GRCh38",
    select_transcripts="mane_select"
)

print(json.dumps(validation.format_as_dict(with_meta=True), indent=4, sort_keys=True))
```

Each variant is validated independently, and the results are returned in the order in which the variants were supplied.

> **Note:** RefSeq and Ensembl variant descriptions must **not** be mixed within the same validation request. Submit RefSeq and Ensembl variants in separate validation requests.

---

### Using the `transcripts` argument

> **Note:** The `transcripts` argument only affects genomic variants. It is ignored when validating transcript variants because the transcript is already explicitly defined by the input variant.

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NM_000088.3:c.589G>T",
    genome="GRCh38",
    select_transcripts="mane_select"
)

print(json.dumps(validation.format_as_dict(with_meta=True), indent=4, sort_keys=True))
```

The transcript specified by the input variant is preserved and is **not** replaced by the MANE Select transcript (`NM_000088.4`).

RefSeq and Ensembl transcript identifiers must **not** be mixed when using the `transcripts` argument.

---

### Restrict the output to MANE Select transcripts

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select"
)

print(json.dumps(validation.format_as_dict(with_meta=True), indent=4, sort_keys=True))
```

---

### Restrict the output to a single specified transcript

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts'["NM_000088.4"]'
)

print(json.dumps(validation.format_as_dict(with_meta=True), indent=4, sort_keys=True))
```

---

### Restrict the output to multiple specified transcripts

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts'["NM_000088.3","NM_000088.4"]'
)

print(json.dumps(validation.format_as_dict(with_meta=True), indent=4, sort_keys=True))
```

---

### Restrict the output to multiple user-specified transcripts

The `transcripts` argument accepts a JSON array of transcript identifiers.

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts'["NM_000088.3","NM_000088.4"]'
)

print(json.dumps(validation.format_as_dict(with_meta=True), indent=4, sort_keys=True))
```

---

### Enable additional liftover

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select",
    liftover_level=True
)

print(json.dumps(validation.format_as_dict(with_meta=True), indent=4, sort_keys=True))
```

---

### Enable the LOVD syntax checker

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select",
    lovd_syntax_check=True
)

print(json.dumps(validation.format_as_dict(with_meta=True), indent=4, sort_keys=True))
```

---

### Enable shorthand VCF parsing

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="17-50198002-C-A",
    genome="GRCh38",
    select_transcripts="mane_select",
    shorthand_vcf=True
)

print(json.dumps(validation.format_as_dict(with_meta=True), indent=4, sort_keys=True))
```

---

### Validate variants from an input file

VariantValidator can validate multiple variants from a text file.

Each line of the input file should contain a single supported variant description.

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="variants.txt",
    genome="GRCh38",
    select_transcripts="mane_select"
)

print(json.dumps(validation.format_as_dict(with_meta=True), indent=4, sort_keys=True))
```

---

### Write the validation results to a JSON file

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select"
)

with open("results.json", "w") as fh:
    json.dump(validation.format_as_dict(with_meta=True), fh, indent=4, sort_keys=True)
```

---

### Generate formatted JSON

```python
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select"
)

print(validation.format_as_json(with_meta=True))
```

---

### Generate tabular output

```python
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select"
)

print(validation.format_as_table(with_meta=True))
```

---

## Common Errors

Common problems include:

- Invalid HGVS syntax.
- Unsupported reference sequences.
- Missing or invalid required method arguments.
- Unable to connect to the VariantValidator databases.
- Missing or incorrect configuration file.
- Invalid transcript selection strategy.
- Invalid JSON arrays supplied to `variant` or `select_transcripts`.

Most exceptions include an explanatory message describing the cause of the problem.

For a complete description of Python exceptions, validation errors and troubleshooting guidance, see the
[Errors and Error Codes](../reference/errors_and_error_codes.md) guide.

---

## Related Documentation

- [VariantValidator Command Line Interface](../cli/variantvalidator_cli.md)
- [Supported Input Formats](../reference/supported_inputs.md)
- [Output Formats](../reference/output_formats.md)
- [Transcript Selection](../reference/transcript_selection.md)
- [Errors and Error Codes](../reference/errors_and_error_codes.md)
