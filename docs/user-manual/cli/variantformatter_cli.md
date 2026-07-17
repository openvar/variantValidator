# VariantFormatter Command Line Interface

The VariantFormatter Command Line Interface (CLI) provides a simple way to convert genetic variant descriptions between genomic, transcript and protein representations directly from the command line. It is suitable for formatting individual variants, processing batches of variants and generating structured JSON output for downstream analysis.

The CLI is intended for users who wish to use VariantFormatter without writing Python code.

For users who are not familiar with command-line tools or Python programming:

- The [VariantValidator website](https://variantvalidator.org) provides a user-friendly interface for formatting and validating variant descriptions.
- The [VariantValidator REST API](https://rest.variantvalidator.org) allows programmatic access to the formatting services without requiring local installation.

---

## Basic Usage

The simplest way to format a variant is to provide the variant description and the genome assembly.

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38
```

VariantFormatter converts the supplied variant into alternative genomic, transcript and protein representations where appropriate and returns the results as JSON.

---

## Command Syntax

```text
variantformatter [OPTIONS]
```

To display the complete list of available options:

```bash
variantformatter --help
```

---

## Required Arguments

The following argument is always required.

| Argument | Description |
|----------|-------------|
| `-v`, `--variant` | The variant description(s) to format. |

---

## Common Options

Commonly used command-line options include:

| Option                       | Description                                                                 |
|------------------------------|-----------------------------------------------------------------------------|
| `-g`, `--genome`             | Specify the reference genome assembly (e.g. `GRCh37` or `GRCh38`).          |
| `-t`, `--select-transcripts` | Restrict the returned transcript representations.                           |
| `--transcript-model`         | Select the transcript database (`refseq`, `ensembl` or `all`).              |
| `--check-only`               | Validate genomic HGVS syntax without transcript or protein mapping.         |
| `-l`, `--liftover-level`      | Generate equivalent genomic representations on alternate genome assemblies. |
| `-o`, `--output`             | Write the results to a file.                                                |
| `--help`                     | Display the command help message.                                           |

---

## Default Behaviour

Unless otherwise specified, VariantFormatter uses the following defaults.

| Setting | Default                    |
|---------|----------------------------|
| Genome assembly | `GRCh38`                   |
| Transcript selection | `mane_select`              |
| Transcript database | `refseq`                   |
| Genomic syntax checking only | Disabled                   |
| Liftover | Enabled                    |
| Output format | JSON                       |
| Output destination | Standard output (`stdout`) |

These defaults can be overridden using the appropriate command-line options described below and demonstrated in the examples.

---

## Supported Input Formats

VariantFormatter accepts a range of commonly used variant descriptions, including:

- Genomic HGVS (g. notation)
- Coding HGVS (c. notation)
- Non-coding HGVS (n. notation)
- RNA HGVS (r. notation)
- Protein HGVS (p. notation)
- Pseudo-VCF/Chromosome coordinate notation (e.g. `17-50198002-C-A` or `17:50198002:C:A`)
- VCF notation (i.e. full VCF lines with chromosome, position, reference and alternate alleles)

See the [Supported Input Formats](../reference/supported_inputs.md) guide for a complete list of supported variant representations.

---

## Output Format

VariantFormatter returns formatted JSON output.

Output can be written directly to the terminal or saved to a file using the `--output` option.

A detailed description of the JSON output is provided in the [Output Formats](../reference/output_formats.md) guide.

---

## Transcript Selection

VariantFormatter supports the same transcript selection strategies as VariantValidator.

These include:

- MANE Select transcripts
- MANE Select and Plus Clinical transcripts
- All transcripts overlapping a genomic variant at their latest version
- All transcripts overlapping a genomic variant at all versions
- User-specified transcript lists

See the [Transcript Selection](../reference/transcript_selection.md) guide for complete details.

---

## Examples

### Format a genomic variant

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38
```

---

### Format a genomic variant using the Ensembl transcript database

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --transcript-model ensembl
```

---

### Format a genomic variant using all transcript databases

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --transcript-model all
```

---

### Format a transcript variant

```bash
variantformatter \
    --variant "NM_000088.4:c.589G>T" \
    --genome GRCh38
```

---

### Format an Ensembl transcript variant

```bash
variantformatter \
    --variant "ENST00000225964.10:c.589G>T" \
    --genome GRCh38 \
    --transcript-model ensembl
```

---

### Format multiple variants

VariantFormatter accepts multiple variants using a JSON array.

```bash
variantformatter \
    --variant '["NC_000017.11:g.50198002C>A","NM_000088.4:c.589G>T"]' \
    --genome GRCh38
```

Each variant is formatted independently, and the results are returned in the order in which the variants were supplied.

**Note:** RefSeq and Ensembl variant descriptions must **not** be mixed within the same formatting request. Submit RefSeq and Ensembl variants in separate commands.

---

### Using the `--select-transcripts` option

> **Note:** The `--select-transcripts` option only affects genomic variants. It is ignored when formatting transcript variants, as the transcript is already explicitly defined in the input.

For example,

```bash
variantformatter \
    --variant "NM_000088.3:c.589G>T" \
    --genome GRCh38 \
    --select-transcripts mane_select
```

returns the formatting for `NM_000088.3:c.589G>T`.

The transcript specified by the input variant is preserved and is **not** replaced by the MANE Select transcript (`NM_000088.4`).

RefSeq and Ensembl transcript identifiers must **not** be mixed when using the `--select-transcripts` option.

---

### Restrict the output to MANE Select transcripts

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts mane_select
```

---

### Restrict the output to a single specified transcript

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts '["NM_000088.4"]'
```

Only the specified transcript is included in the output.

---

### Restrict the output to multiple specified transcripts

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts '["NM_000088.3","NM_000088.4"]'
```

Only the specified transcript versions are included in the output.

---

### Restrict the output to multiple user-selected transcripts

The `--select-transcripts` option accepts a JSON array of transcript identifiers.

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts '["NM_000088.3","NM_000088.4"]'
```

---

### Validate genomic HGVS syntax only

The `--check-only` option validates genomic HGVS syntax without generating transcript or protein representations.

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --check-only
```

---

### Generate lifted-over genomic representations

The `--liftover-level` option includes equivalent genomic representations on alternate genome assemblies.

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --liftover-level
```
Usage

| Parameter | Type            | Required | Description |
|----------|-----------------|----------|-------------|
| liftover_level | string or bool  | No | Controls genomic liftover. `True` performs full liftover, `primary` excludes alternative scaffolds, and `False` disables liftover. Defaults to `True`. |

---

### Format pseudo-VCF notation

VariantFormatter accepts pseudo-VCF chromosome coordinate notation.

```bash
variantformatter \
    --variant "17-50198002-C-A" \
    --genome GRCh38
```

or

```bash
variantformatter \
    --variant "17:50198002:C:A" \
    --genome GRCh38
```

---

### Format variants from an input file

VariantFormatter can format multiple variants from a text file.

Each line of the input file should contain a single supported variant description.

```bash
variantformatter \
    --variant variants.txt \
    --genome GRCh38
```

---

### Write the formatted results to a JSON file

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --output results.json
```

---

### Display the command help

```bash
variantformatter --help
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

For a complete description of command-line error messages, exit codes and troubleshooting guidance, see the
[Errors and Error Codes](../reference/errors_and_error_codes.md) guide.

---

## Related Documentation

- [VariantFormatter Python API](../python-api/variantformatter_python.md)
- [VariantValidator Command Line Interface](variantvalidator_cli.md)
- [Supported Input Formats](../reference/supported_inputs.md)
- [Output Formats](../reference/output_formats.md)
- [Transcript Selection](../reference/transcript_selection.md)
- [Errors and Error Codes](../reference/errors_and_error_codes.md)


