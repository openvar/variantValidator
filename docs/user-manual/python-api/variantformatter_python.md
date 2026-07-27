# VariantFormatter Python API

The VariantFormatter Python API provides direct access to the VariantFormatter formatting engine from within Python. It is suitable for integrating genomic variant formatting into bioinformatics pipelines, analysis workflows, web applications and custom software.

VariantFormatter accepts genomic variant descriptions as input and generates corresponding genomic, transcript and protein representations where appropriate.

For users who prefer not to write Python code:

- The [VariantValidator website](https://variantvalidator.org) provides a user-friendly interface for formatting and validating variant descriptions.
- The [VariantValidator REST API](https://rest.variantvalidator.org) allows programmatic access to VariantValidator services without requiring local installation.
- The [VariantFormatter Command Line Interface](../cli/variantformatter_cli.md) provides a command-line interface for formatting variants locally.

---

## Basic Usage

Begin by importing VariantFormatter and creating a `SimpleVariantFormatter` object.

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()
```

The `SimpleVariantFormatter` object manages access to the VariantFormatter formatting engine and can be reused to format multiple variants within the same Python session.

Variants are formatted using the `format()` method.

```python
results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select",
)

print(json.dumps(results, indent=4, sort_keys=True))
```

The returned Python dictionary can be processed directly or converted to JSON.

---

## Method Signature

Variant formatting is performed using the `format()` method.

```python
format(
    variant,
    genome,
    transcript_model="refseq",
    select_transcripts="mane_select",
    checkOnly=False,
    liftover_level=True,
    legacy_genomic_structure=True,
)
```

---

## Required Arguments

| Argument | Description |
|----------|-------------|
| `variant` | A genomic variant description or multiple genomic variants supplied in a supported batch format. |
| `genome` | Reference genome assembly: `GRCh37`, `GRCh38`, `hg19` or `hg38`. |

---

## Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `transcript_model` | `refseq` | Transcript database to use: `refseq`, `ensembl` or `all`. |
| `select_transcripts` | `mane_select` | Controls which transcript representations are returned. |
| `checkOnly` | `False` | Validate and format the genomic variant without transcript or protein mapping. |
| `liftover_level` | `True` | Controls generation of genomic representations on another genome assembly. |
| `legacy_genomic_structure` | `True` | Preserve the historical VariantFormatter genomic loci structure. Set to `False` to return the VariantValidator genomic loci structure. |

---

## Default Behaviour

Unless otherwise specified, VariantFormatter uses the following behaviour:

| Setting | Default |
|---------|---------|
| Genome assembly | User supplied |
| Transcript selection | MANE Select |
| Transcript database | RefSeq |
| Genomic syntax checking only | Disabled |
| Liftover | Enabled |
| Genomic loci structure | Legacy VariantFormatter structure |
| Output format | Python dictionary |

For example:

```python
results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
)
```

A detailed description of the returned data is provided in the [Output Formats](../reference/output_formats.md) guide.

---

## Supported Input Formats

VariantFormatter accepts **genomic variants as input**.

Supported input formats include:

- Genomic HGVS (`g.` notation) using supported genomic reference sequences, including `NC_`, `NT_` and `NW_` accessions.
- Pseudo-VCF chromosome-coordinate notation, for example:
  - `17-50198002-C-A`
  - `17:50198002:C:A`

Transcript (`c.` and `n.`), RNA (`r.`) and protein (`p.`) HGVS descriptions are not accepted as VariantFormatter input.

VariantFormatter operates from a genomic variant and maps that variant to relevant transcript and protein representations.

See the [Supported Input Formats](../reference/supported_inputs.md) guide for further details.

---

## Transcript Selection

VariantFormatter maps genomic variants to overlapping transcripts.

The `select_transcripts` argument controls which transcripts are returned.

Supported transcript selection strategies include:

- `mane_select` — MANE Select transcripts.
- `mane` — MANE Select and MANE Plus Clinical transcripts.
- `all` — all relevant transcripts at their latest version.
- `raw` — all relevant transcript versions.
- Explicit user-selected transcript identifiers.

See the [Transcript Selection](../reference/transcript_selection.md) guide for further details.

---

## Examples

### Format a genomic variant

```python
import json
from VariantFormatter.simpleVariantFormatter import SimpleVariantFormatter

formatter = SimpleVariantFormatter()

results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Format a genomic variant using Ensembl transcripts

```python
results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    transcript_model="ensembl",
)
```

---

### Format a genomic variant using RefSeq and Ensembl transcripts

```python
results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    transcript_model="all",
)
```

---

### Format pseudo-VCF notation

Hyphen-delimited pseudo-VCF input:

```python
results = formatter.format(
    variant="17-50198002-C-A",
    genome="GRCh38",
)
```

Colon-delimited pseudo-VCF input:

```python
results = formatter.format(
    variant="17:50198002:C:A",
    genome="GRCh38",
)
```

---

## Format Multiple Variants

Multiple variants can be supplied as a JSON array.

```python
results = formatter.format(
    variant=(
        '["NC_000017.11:g.50198002C>A",'
        '"NC_000016.10:g.15738651_15738652inv"]'
    ),
    genome="GRCh38",
)
```

Each genomic variant is processed independently and returned in the result dictionary.

---

## Selecting Transcripts

### Restrict output to MANE Select transcripts

```python
results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane_select",
)
```

---

### Return MANE Select and MANE Plus Clinical transcripts

```python
results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="mane",
)
```

---

### Return all latest transcript versions

```python
results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="all",
)
```

---

### Return all transcript versions

```python
results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts="raw",
)
```

---

### Restrict output to a single specified transcript

```python
results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts='["NM_000088.4"]',
)
```

---

### Restrict output to multiple specified transcripts

```python
results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    select_transcripts='["NM_000088.3","NM_000088.4"]',
)
```

RefSeq and Ensembl transcript identifiers must not be mixed in the same explicit transcript list.

---

## Validate Genomic HGVS Only

The `checkOnly` argument validates and formats the genomic variant without generating transcript or protein mappings.

```python
results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    checkOnly=True,
)
```

---

## Liftover

VariantFormatter can generate equivalent genomic representations on another genome assembly.

The `liftover_level` argument controls this behaviour.

| Value | Description |
|-------|-------------|
| `True` | Perform full liftover. |
| `"primary"` | Perform liftover while excluding alternative scaffolds. |
| `False` | Disable liftover. |

For example:

```python
results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    liftover_level=True,
)
```

To disable liftover:

```python
results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
    liftover_level=False,
)
```

---

## Genomic Loci Output Structure

VariantFormatter historically uses a genomic loci structure in which each genome build contains an additional accession-keyed level.

This structure remains the default for backwards compatibility:

```python
results = formatter.format(
    variant="NC_000016.10:g.15738651_15738652inv",
    genome="GRCh38",
    legacy_genomic_structure=True,
)
```

For example, `primary_assembly_loci` has the form:

```python
{
    "grch38": {
        "NC_000016.10": {
            "hgvs_genomic_description":
                "NC_000016.10:g.15738651_15738652inv",
            "vcf": {
                "chr": "16",
                "pos": "15738651",
                "ref": "GT",
                "alt": "AC",
            },
        }
    }
}
```

### VariantValidator genomic structure

Set `legacy_genomic_structure=False` to return genomic loci using the VariantValidator structure:

```python
results = formatter.format(
    variant="NC_000016.10:g.15738651_15738652inv",
    genome="GRCh38",
    legacy_genomic_structure=False,
)
```

The additional accession-keyed level is removed:

```python
{
    "grch38": {
        "hgvs_genomic_description":
            "NC_000016.10:g.15738651_15738652inv",
        "vcf": {
            "chr": "16",
            "pos": "15738651",
            "ref": "GT",
            "alt": "AC",
        },
    }
}
```

This option affects the structure used to return genomic loci; it does not change the underlying variant mapping.

---

## Write Results to a JSON File

The returned dictionary can be written directly to JSON.

```python
import json

results = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
)

with open("results.json", "w") as fh:
    json.dump(
        results,
        fh,
        indent=4,
        sort_keys=True,
    )
```

---

## Reusing the Formatter

A `SimpleVariantFormatter` instance can be reused for multiple formatting requests.

```python
formatter = SimpleVariantFormatter()

result_1 = formatter.format(
    variant="NC_000017.11:g.50198002C>A",
    genome="GRCh38",
)

result_2 = formatter.format(
    variant="NC_000016.10:g.15738651_15738652inv",
    genome="GRCh38",
)
```

Reusing the formatter avoids unnecessarily recreating the VariantFormatter environment for each request.

---

## Common Errors

Common problems include:

- invalid genomic HGVS syntax;
- unsupported input types;
- unsupported reference sequences;
- a reference sequence that does not correspond to the selected genome build;
- invalid pseudo-VCF input;
- invalid transcript selection;
- invalid transcript model selection;
- coordinates outside the reference sequence;
- inability to connect to the VariantValidator databases;
- a missing or incorrect VariantValidator configuration file.

VariantFormatter returns explanatory warnings or errors where possible.

For further information, see the [Errors and Error Codes](../reference/errors_and_error_codes.md) guide.

---

## Related Documentation

- [VariantFormatter Command Line Interface](../cli/variantformatter_cli.md)
- [VariantValidator Python API](variantvalidator_python.md)
- [Supported Input Formats](../reference/supported_inputs.md)
- [Output Formats](../reference/output_formats.md)
- [Transcript Selection](../reference/transcript_selection.md)
- [Errors and Error Codes](../reference/errors_and_error_codes.md)