# VariantValidator Command Line Interface

The VariantValidator Command Line Interface (CLI) provides a simple way to validate and normalise genetic variant
descriptions directly from the command line. It is suitable for validating individual variants, processing batches of
variants, and generating structured output for downstream analysis.

The CLI is intended for users who wish to use VariantValidator without writing Python code.

For users who are not familiar with command-line tools or Python programming:

- The [VariantValidator website](https://variantvalidator.org) provides a user-friendly alternative for validating variant descriptions.
- The [VariantValidator REST API](https://rest.variantvalidator.org) allows programmatic access to the validation services without requiring local installation.

---

## Basic Usage

The simplest way to validate a variant is to provide the variant description and the genome assembly.

```bash
variantvalidator \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38
```

VariantValidator validates the supplied variant, performs normalisation where appropriate, and returns the validation results.

---

## Command Syntax

```text
variantvalidator [OPTIONS]
```

To display the complete list of available options:

```bash
variantvalidator --help
```

---

## Required Arguments

The following argument is always required:

| Argument | Description |
|----------|-------------|
| `-v`, `--variant` | The variant description(s) to validate. |

---

## Common Options

Commonly used command-line options include:

| Option                       | Description |
|------------------------------|-------------|
| `-g`, `--genome`             | Specify the reference genome assembly (e.g. `GRCh37` or `GRCh38`). |
| `-t`, `--select-transcripts` | Restrict the output to selected transcripts. |
| `--transcript-set`           | Select a predefined transcript set. |
| `-f`, `--output-format`      | Specify the output format. |
| `-l`, `--liftover-level`      | Generate equivalent genomic representations on alternate genome assemblies. |
| `-o`, `--output`             | Write the results to a file. |
| `-m`, `--meta`               | Include metadata in the output. |
| `--help`                     | Display the command help message. |

---

## Default Behaviour

Unless otherwise specified, VariantValidator uses the following defaults:

| Setting              | Default                    |
|----------------------|----------------------------|
| Genome assembly      | `GRCh38`                   |
| Transcript selection | `mane_select`              |
| Transcript database  | `refseq`                   |
| Liftover             | Enabled                    |
| LOVD syntax checker  | Enabled                    |
| Output format        | JSON                       |
| Output destination   | Standard output (`stdout`) |
| Metadata             | Disabled                   |

These defaults can be overridden using the appropriate command-line options described above and demonstrated in the examples below.

---

## Supported Input Formats

VariantValidator accepts a range of commonly used variant descriptions, including:

- Genomic HGVS (g. notation)
- Coding HGVS (c. notation)
- Non-coding HGVS (n. notation)
- RNA HGVS (r. notation)
- Protein HGVS (p. notation, including single-letter and three-letter amino acid codes - basic validation, not recommended)
- Pseudo-VCF/Chromosome coordinate notation (e.g. `17-50198002-C-A` or `17:50198002:C:A`)
- VCF notation (i.e. full VCF lines with chromosome, position, reference and alternate alleles)

See the [Supported Input Formats](../reference/supported_inputs.md) guide for a complete list of supported variant representations.

---

## Output Formats

VariantValidator supports two output formats suitable for interactive use and downstream processing.

Available output formats include:

- JSON
- Tabular output

Output can be written directly to the terminal or saved to a file using the `--output` option.

A detailed description of each output format, including the JSON schema and tabular output columns, is provided in the
[Output Formats](../reference/output_formats.md) guide.

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

```bash
variantvalidator \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38
```

---

### Validate a genomic variant using the Ensembl transcript set

```bash
variantvalidator \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --transcript-set ensembl
```

---

### Validate a transcript variant

```bash
variantvalidator \
    --variant "NM_000088.4:c.589G>T" \
    --genome GRCh38
```

---

### Validate an Ensembl transcript variant

```bash
variantvalidator \
    --variant "ENST00000225964.10:c.589G>T" \
    --genome GRCh38 \
    --transcript-set ensembl
```

---

### Validate multiple variants

VariantValidator accepts multiple variants using a JSON array.

```bash
variantvalidator \
    --variant '["NC_000017.11:g.50198002C>A","NM_000088.4:c.589G>T"]' \
    --genome GRCh38
```

Each variant is validated independently, and the results are returned in the order in which the variants were supplied.

**Note:** RefSeq and Ensembl variant descriptions must **not** be mixed within the same validation request. Submit RefSeq and Ensembl variants in separate commands.

---

### Using the `--select-transcripts` option

> **Note:** The `--select-transcripts` option only affects genomic variants. It is ignored when validating transcript variants, as the transcript is already explicitly defined in the input.

For example,

```bash
variantvalidator \
    --variant "NM_000088.3:c.589G>T" \
    --genome GRCh38 \
    --select-transcripts mane_select
```

returns the validation for `NM_000088.3:c.589G>T`.

The transcript specified by the input variant is preserved and is **not** replaced by the MANE Select transcript (`NM_000088.4`).

RefSeq and Ensembl transcript identifiers must **not** be mixed when using the `--select-transcripts` option.

---

### Restrict the output to MANE Select transcripts

```bash
variantvalidator \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts mane_select
```

---

### Restrict the output to a single specified transcript

```bash
variantvalidator \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts '["NM_000088.4"]'
```

Only the specified transcript is included in the output.

---

### Restrict the output to multiple specified transcripts

```bash
variantvalidator \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts '["NM_000088.3","NM_000088.4"]'
```

Only the specified transcript versions are included in the output.

---

### Restrict the output to multiple user-selected transcripts

The `--select-transcripts` option accepts a JSON array of transcript identifiers.

```bash
variantvalidator \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts '["NM_000088.3","NM_000088.4"]'
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

| Parameter | Type            | Required | Description |
|----------|-----------------|----------|-------------|
| liftover_level | string or bool  | No | Controls genomic liftover. `True` performs full liftover, `primary` excludes alternative scaffolds, and `False` disables liftover. Defaults to `True`. |

---

### Write the results to a JSON file

```bash
variantvalidator \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --output-format json \
    --output results.json
```

---

### Display the results as a table

```bash
variantvalidator \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --output-format table
```

See the [Output Formats](../reference/output_formats.md) guide for a detailed description of the available output formats.

---

### Write the results as a table

```bash
variantvalidator \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --output-format table \
    --output results.tsv
```

---

### Validate variants from an input file

VariantValidator can validate multiple variants from a text file.

Each line of the input file should contain a single supported variant description.

```bash
variantvalidator \
    --variant variants.txt \
    --genome GRCh38
```

---

### Validate variants from an input file and write the results to a JSON file

```bash
variantvalidator \
    --variant variants.txt \
    --genome GRCh38 \
    --output-format json \
    --output results.json
```

---

### Display the command help

```bash
variantvalidator --help
```

---

## Common Errors

Common problems include:

- Invalid HGVS syntax.
- Unsupported reference sequences.
- Missing genome build.
- Unable to connect to the VariantValidator databases.
- Missing or incorrect configuration file.

Most errors include an explanatory message describing the cause of the problem.

For a complete description of command-line error messages, exit codes and troubleshooting guidance, see the
[Errors and Error Codes](../reference/errors_and_error_codes.md) guide.

---

## Related Documentation

- [VariantValidator Python API](../python-api/variantvalidator_python.md)
- [Supported Input Formats](../reference/supported_inputs.md)
- [Output Formats](../reference/output_formats.md)
- [Transcript Selection](../reference/transcript_selection.md)
- [Errors and Error Codes](../reference/errors_and_error_codes.md)