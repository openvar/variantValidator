# VariantFormatter Python API

The VariantFormatter Python API provides direct access to the VariantFormatter formatting engine from within Python. It is suitable for integrating variant formatting into bioinformatics pipelines, analysis workflows, web applications and custom software.

The Python API offers access to the same formatting functionality as the VariantFormatter command-line interface while providing a programmatic interface for automated processing and downstream analysis.

For users who prefer not to write Python code:

- The [VariantValidator website](https://variantvalidator.org) provides a user-friendly interface for formatting and validating variant descriptions.
- The [VariantValidator REST API](https://rest.variantvalidator.org) allows programmatic access to the formatting services without requiring local installation.
- The [VariantFormatter Command Line Interface](../cli/variantformatter_cli.md) provides a command-line interface for formatting variants locally.

---

## Basic Usage

Begin by importing the VariantFormatter package and creating a `SimpleVariantFormatter` object.

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()
```

The `SimpleVariantFormatter` object manages access to the VariantFormatter formatting engine and can be reused to format multiple variants within the same Python session.

Once a `SimpleVariantFormatter` object has been created, variants can be formatted using the `format()` method.

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome_build="GRCh38",
    select_transcripts="mane_select"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

This formats the supplied variant and returns the results as a Python dictionary that can be processed directly or converted to JSON.

---

## Method Signature

Variant formatting is performed using the `format()` method.

```python
format(
    variant=None,
    genome_build=None,
    transcript_model=None,
    select_transcripts=None,
    checkOnly=False,
    liftover=False
)
```

---

## Required Arguments

The following arguments are required when calling the `format()` method.

| Argument | Description |
|----------|-------------|
| `variant` | A single variant, multiple variants as a JSON array, or a filename containing variants to format. |
| `genome_build` | The reference genome assembly (e.g. `GRCh37` or `GRCh38`). |

---

## Optional Arguments

The following optional arguments control formatting behaviour.

| Argument | Default | Description |
|----------|---------|-------------|
| `transcript_model` | `None` | Select the transcript database (`refseq`, `ensembl` or `all`). |
| `select_transcripts` | `None` | Restrict the returned transcript representations. |
| `checkOnly` | `False` | Validate genomic HGVS syntax only without transcript or protein mapping. |
| `liftover` | `False` | Generate equivalent genomic representations on compatible genome assemblies. |

Usage of liftover_level

| Parameter | Type            | Required | Description |
|----------|-----------------|----------|-------------|
| liftover_level | string or bool  | No | Controls genomic liftover. `True` performs full liftover, `primary` excludes alternative scaffolds, and `False` disables liftover. Defaults to `True`. |

---

## Default Behaviour

Unless otherwise specified, the VariantFormatter Python API uses the following defaults.

| Setting | Default |
|---------|---------|
| Genome assembly | User supplied (required) |
| Transcript selection | All compatible transcripts |
| Transcript database | `refseq` |
| Genomic syntax checking only | Disabled |
| Liftover | Disabled |
| Output format | Python dictionary |

The `format()` method returns a Python dictionary containing the formatted variant representations.

For example,

```python
results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome_build="GRCh38"
)
```

The returned dictionary can be processed directly or converted to formatted JSON.

```python
print(json.dumps(results, indent=4, sort_keys=True))
```

A detailed description of the output format is provided in the
[Output Formats](../reference/output_formats.md) guide.

---

## Supported Input Formats

VariantFormatter accepts the same input formats as the command-line interface.

Supported variant descriptions include:

- Genomic HGVS (g. notation)
- Coding HGVS (c. notation)
- Non-coding HGVS (n. notation)
- RNA HGVS (r. notation)
- Protein HGVS (p. notation)
- Pseudo-VCF/Chromosome coordinate notation (e.g. `17-50198002-C-A` or `17:50198002:C:A`)
- VCF notation (i.e. full VCF lines with chromosome, position, reference and alternate alleles)

The `variant` argument accepts:

- A single variant description.
- Multiple variant descriptions supplied as a JSON array.
- A text file containing one variant description per line.

See the [Supported Input Formats](../reference/supported_inputs.md) guide for a complete description of supported input formats and examples.

---

## Transcript Selection

VariantFormatter supports multiple transcript selection strategies.

These include:

- MANE Select transcripts
- MANE Select and Plus Clinical transcripts
- All transcripts overlapping a genomic variant at their latest version
- All transcripts overlapping a genomic variant at all versions
- User-specified transcript lists

See the [Transcript Selection](../reference/transcript_selection.md) guide for details.

---

## Examples

### Format a genomic variant

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome_build="GRCh38"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Format a genomic variant using the Ensembl transcript database

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome_build="GRCh38",
    transcript_model="ensembl"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Format a genomic variant using all transcript databases

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome_build="GRCh38",
    transcript_model="all"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Format a transcript variant

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="NM_000088.4:c.589G>T",
    genome_build="GRCh38"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Format an Ensembl transcript variant

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="ENST00000225964.10:c.589G>T",
    genome_build="GRCh38",
    transcript_model="ensembl"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Format multiple variants

VariantFormatter accepts multiple variants using a JSON array.

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant='["NC_000017.11:g.50198002C>A","NM_000088.4:c.589G>T"]',
    genome_build="GRCh38"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

Each variant is formatted independently, and the results are returned in the order in which the variants were supplied.

**Note:** RefSeq and Ensembl variant descriptions must **not** be mixed within the same formatting request. Submit RefSeq and Ensembl variants in separate formatting requests.

---

### Using `select_transcripts`

> **Note:** The `select_transcripts` argument only affects genomic variants. It is ignored when formatting transcript variants because the transcript is already explicitly defined by the input variant.

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="NM_000088.3:c.589G>T",
    genome_build="GRCh38",
    select_transcripts="mane_select"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

The transcript specified by the input variant is preserved and is **not** replaced by the MANE Select transcript (`NM_000088.4`).

RefSeq and Ensembl transcript identifiers must **not** be mixed when using `select_transcripts`.

---

### Restrict the output to MANE Select transcripts

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome_build="GRCh38",
    select_transcripts="mane_select"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Restrict the output to a single specified transcript

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome_build="GRCh38",
    select_transcripts='["NM_000088.4"]'
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Restrict the output to multiple specified transcripts

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome_build="GRCh38",
    select_transcripts='["NM_000088.3","NM_000088.4"]'
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Restrict the output to multiple user-selected transcripts

The `select_transcripts` argument accepts a JSON array of transcript identifiers.

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome_build="GRCh38",
    select_transcripts='["NM_000088.3","NM_000088.4"]'
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Validate genomic HGVS syntax only

The `checkOnly` argument validates genomic HGVS syntax without generating transcript or protein representations.

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome_build="GRCh38",
    checkOnly=True
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Generate lifted-over genomic representations

The `liftover` argument includes equivalent genomic representations on compatible genome assemblies.

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome_build="GRCh38",
    liftover=True
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Format pseudo-VCF notation

VariantFormatter accepts pseudo-VCF chromosome coordinate notation.

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="17-50198002-C-A",
    genome_build="GRCh38"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

or

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="17:50198002:C:A",
    genome_build="GRCh38"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Format variants from an input file

VariantFormatter can format multiple variants from a text file.

Each line of the input file should contain a single supported variant description.

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="variants.txt",
    genome_build="GRCh38"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Write the formatted results to a JSON file

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome_build="GRCh38"
)

with open("results.json", "w") as fh:
    json.dump(results, fh, indent=4, sort_keys=True)
```

---

## Common Errors

Common problems include:

- Invalid HGVS syntax.
- Unsupported reference sequences.
- Missing genome build.
- Invalid transcript selection.
- Unable to connect to the VariantValidator databases.
- Missing or incorrect configuration file.

Most errors include an explanatory message describing the cause of the problem.

For a complete description of Python exceptions, error messages and troubleshooting guidance, see the
[Errors and Error Codes](../reference/errors_and_error_codes.md) guide.

---

## Related Documentation

- [VariantFormatter Command Line Interface](../cli/variantformatter_cli.md)
- [VariantValidator Python API](variantvalidator_python.md)
- [Supported Input Formats](../reference/supported_inputs.md)
- [Output Formats](../reference/output_formats.md)
- [Transcript Selection](../reference/transcript_selection.md)
- [Errors and Error Codes](../reference/errors_and_error_codes.md)


