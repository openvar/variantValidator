<img src="../../static/img/logos/VV_logo.png" width="20%" />

# VariantValidator Python API

The VariantValidator Python API provides direct access to the VariantValidator validation engine from within Python. It is suitable for integrating variant validation into bioinformatics pipelines, analysis workflows, web applications and custom software.

The Python API offers access to the same validation functionality as the VariantValidator command-line interface while providing a programmatic interface for automated processing and downstream analysis.

For users who prefer not to write Python code:

- The [VariantValidator website](https://variantvalidator.org) provides a user-friendly alternative for validating variant descriptions.
- The [VariantValidator REST API](https://rest.variantvalidator.org) allows programmatic access to the validation services without requiring local installation.
- The [VariantValidator Command Line Interface](../cli/variantvalidator_cli.md) provides a command-line interface for validating variants locally.

## See also

- [VariantValidator Command Line Interface](../cli/variantvalidator_cli.md) — Validate variant descriptions from the command line.
- [Gene2Transcripts Python API](gene2transcripts_python.md) — Retrieve transcript information directly from Python.
- [hgvs2reference Python API](hgvs2reference_python.md) — Retrieve reference sequence for HGVS variants.
- [Supported Input Formats](../reference/supported_inputs.md) — Supported variant description formats.
- [Output Formats](../reference/output_formats.md) — Available output formats and returned data.
- [Transcript Selection](../reference/transcript_selection.md) — Available transcript selection strategies.
- [Errors and Error Codes](../reference/errors_and_error_codes.md) — Error messages and troubleshooting guidance.

---

# Basic Usage

Begin by importing the VariantValidator package and creating a `Validator` object.

```python
import VariantValidator

vval = VariantValidator.Validator()
```

The `Validator` object manages access to the VariantValidator validation engine and can be reused to validate multiple variants within the same Python session.

Once a `Validator` object has been created, variants can be validated using the `validate()` method.

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select"
)

results = validation.format_as_dict(with_meta=True)

print(json.dumps(results, indent=4, sort_keys=True))
```

This validates the supplied variant and returns a `ValOutput` object. The validation results can then be returned as a Python dictionary, JSON or tabular output.

---

# Method Signature

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

The method returns a `ValOutput` object containing the validation results.

---

# Required Arguments

The following arguments are required when calling the `validate()` method.

| Argument | Description |
| --- | --- |
| `variant` | A single variant, multiple variants as a JSON array, or a filename containing variants to validate. |
| `genome` | The reference genome assembly, for example `GRCh37` or `GRCh38`. |

---

# Optional Arguments

The following optional arguments control validation behaviour.

| Argument | Default | Description |
| --- | --- | --- |
| `select_transcripts` | `"all"` | Transcript selection strategy or a JSON array of transcript identifiers. |
| `transcript_set` | `None` | Select the transcript database (`refseq` or `ensembl`). `None` defaults to `refseq`. |
| `liftover_level` | `True` | Controls generation of genomic representations on additional genome builds. |
| `lovd_syntax_check` | `False` | Enable LOVD HGVS syntax checking. |
| `shorthand_vcf` | `False` | Enable shorthand VCF parsing. |

## `liftover_level`

The `liftover_level` argument controls genomic liftover.

| Value | Description |
| --- | --- |
| `True` | Perform full liftover. |
| `"primary"` | Perform liftover while excluding alternative scaffolds. |
| `False` | Disable liftover. |

---

# Default Behaviour

Unless otherwise specified, the VariantValidator Python API uses the following defaults.

| Setting | Default |
| --- | --- |
| Genome assembly | User supplied |
| Transcript selection | `"all"` |
| Transcript database | `refseq` |
| Liftover | Enabled |
| LOVD syntax checker | Disabled |
| Shorthand VCF parsing | Disabled |
| Return type | `ValOutput` object |

For example:

```python
validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38"
)
```

The `ValOutput` object provides several methods for accessing the validation results in different formats.

The optional `with_meta` argument used by the output formatting methods controls whether validation metadata is included in the returned output. Setting `with_meta=True` includes additional information such as the software version, database versions, transcript selection strategy and other validation metadata.

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

A detailed description of each output format, including the available formatting options, is provided in the [Output Formats](../reference/output_formats.md) guide.

---

# Supported Input Formats

The VariantValidator Python API accepts the same input formats as the command-line interface.

Supported variant descriptions include:

- Genomic HGVS (`g.` notation)
- Coding HGVS (`c.` notation)
- Non-coding HGVS (`n.` notation)
- RNA HGVS (`r.` notation)
- Protein HGVS (`p.` notation, including single-letter and three-letter amino acid codes — basic validation, not recommended)
- Pseudo-VCF/chromosome coordinate notation, for example `17-50198002-C-A` or `17:50198002:C:A`
- VCF notation, including full VCF lines containing chromosome, position, reference and alternate alleles

The `variant` argument accepts:

- A single variant description.
- Multiple variant descriptions supplied as a JSON array.
- A text file containing one variant description per line.

See the [Supported Input Formats](../reference/supported_inputs.md) guide for a complete description of supported input formats and examples.

---

# Transcript Selection

VariantValidator supports multiple transcript selection strategies.

These include:

- MANE Select transcripts
- MANE Select and Plus Clinical transcripts
- All transcripts overlapping a genomic variant at their latest version
- All transcripts overlapping a genomic variant at all versions
- User-specified transcript lists

The `select_transcripts` argument controls transcript selection.

See the [Transcript Selection](../reference/transcript_selection.md) guide for details.

---

# Examples

## Validate a genomic variant

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select"
)

print(json.dumps(
    validation.format_as_dict(with_meta=True),
    indent=4,
    sort_keys=True
))
```

---

## Validate a genomic variant using the Ensembl transcript set

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

print(json.dumps(
    validation.format_as_dict(with_meta=True),
    indent=4,
    sort_keys=True
))
```

---

## Validate a transcript variant

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NM_000088.4:c.589G>T",
    genome="GRCh38",
    select_transcripts="mane_select"
)

print(json.dumps(
    validation.format_as_dict(with_meta=True),
    indent=4,
    sort_keys=True
))
```

---

## Validate an Ensembl transcript variant

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

print(json.dumps(
    validation.format_as_dict(with_meta=True),
    indent=4,
    sort_keys=True
))
```

---

## Validate multiple variants

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

print(json.dumps(
    validation.format_as_dict(with_meta=True),
    indent=4,
    sort_keys=True
))
```

Each variant is validated independently, and the results are returned in the order in which the variants were supplied.

> **Note:** RefSeq and Ensembl variant descriptions must **not** be mixed within the same validation request. Submit RefSeq and Ensembl variants in separate validation requests.

---

# Using `select_transcripts`

The `select_transcripts` argument controls which transcripts are returned when validating genomic variants.

> **Note:** `select_transcripts` only affects genomic variants. It is ignored when validating transcript variants because the transcript is already explicitly defined by the input variant.

For example:

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="NM_000088.3:c.589G>T",
    genome="GRCh38",
    select_transcripts="mane_select"
)

print(json.dumps(
    validation.format_as_dict(with_meta=True),
    indent=4,
    sort_keys=True
))
```

The transcript specified by the input variant is preserved and is **not** replaced by the MANE Select transcript (`NM_000088.4`).

RefSeq and Ensembl transcript identifiers must **not** be mixed when using `select_transcripts`.

---

## Restrict the output to MANE Select transcripts

```python
validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select"
)
```

---

## Restrict the output to a single specified transcript

```python
validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts='["NM_000088.4"]'
)
```

---

## Restrict the output to multiple specified transcripts

```python
validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts='["NM_000088.3","NM_000088.4"]'
)
```

RefSeq and Ensembl transcript identifiers must not be mixed in the same explicit transcript list.

---

# Liftover

VariantValidator can generate genomic representations on additional genome assemblies.

The `liftover_level` argument controls this behaviour.

| Value | Description |
| --- | --- |
| `True` | Perform full liftover. |
| `"primary"` | Perform liftover while excluding alternative scaffolds. |
| `False` | Disable liftover. |

For example:

```python
validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select",
    liftover_level=True
)
```

To disable liftover:

```python
validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select",
    liftover_level=False
)
```

---

# Enable the LOVD Syntax Checker

LOVD HGVS syntax checking can be enabled using `lovd_syntax_check=True`.

```python
validation = vval.validate(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select",
    lovd_syntax_check=True
)
```

---

# Enable Shorthand VCF Parsing

Shorthand VCF parsing can be enabled using `shorthand_vcf=True`.

```python
validation = vval.validate(
    variant="17-50198002-C-A",
    genome="GRCh38",
    select_transcripts="mane_select",
    shorthand_vcf=True
)
```

---

# Validate Variants from an Input File

VariantValidator can validate multiple variants from a text file.

Each line of the input file should contain a single supported variant description.

```python
import VariantValidator

vval = VariantValidator.Validator()

validation = vval.validate(
    variant="variants.txt",
    genome="GRCh38",
    select_transcripts="mane_select"
)
```

---

# Write Validation Results to a JSON File

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
    json.dump(
        validation.format_as_dict(with_meta=True),
        fh,
        indent=4,
        sort_keys=True
    )
```

---

# Generate Formatted JSON

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

# Generate Tabular Output

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

# Common Errors

Common problems include:

- Invalid HGVS syntax.
- Unsupported reference sequences.
- Missing or invalid required method arguments.
- Unable to connect to the VariantValidator databases.
- Missing or incorrect configuration file.
- Invalid transcript selection strategy.
- Invalid JSON arrays supplied to `variant` or `select_transcripts`.

Most errors include an explanatory message describing the cause of the problem.

For a complete description of validation errors and troubleshooting guidance, see the [Errors and Error Codes](../reference/errors_and_error_codes.md) guide.

---

# Getting help

VariantValidator has been developed to support a wide range of users, from those new to HGVS nomenclature to experienced clinical scientists and bioinformaticians. If you encounter difficulties using the VariantValidator Python API or interpreting validation results, we encourage you to seek assistance.

Before contacting the development team, you may find the following documentation helpful:

- [VariantValidator Command Line Interface](../cli/variantvalidator_cli.md)
- [Supported Input Formats](../reference/supported_inputs.md)
- [Output Formats](../reference/output_formats.md)
- [Transcript Selection](../reference/transcript_selection.md)
- [Errors and Error Codes](../reference/errors_and_error_codes.md)

If you still require assistance, you can contact the VariantValidator team using our [contact form](https://variantvalidator.org/help/contact/).

Software bugs and feature requests can be reported through the [VariantValidator GitHub issue tracker](https://github.com/openvar/VariantValidator/issues).

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="../../static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="../../static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>
