# VariantFormatter Command Line Interface

The VariantFormatter Command Line Interface (CLI) provides a simple way to format genomic variant descriptions and generate corresponding genomic, transcript and protein representations directly from the command line.

It is suitable for formatting individual variants, processing batches of variants and generating structured JSON output for downstream analysis.

The CLI is intended for users who wish to use VariantFormatter without writing Python code.

For users who are not familiar with command-line tools or Python programming:

- The [VariantValidator website](https://variantvalidator.org) provides a user-friendly interface for formatting and validating variant descriptions.
- The [VariantValidator REST API](https://rest.variantvalidator.org) allows programmatic access to the formatting services without requiring local installation.

---

## Basic Usage

The simplest way to format a variant is to provide a genomic variant description and the genome assembly.

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38
```

VariantFormatter validates and formats the supplied genomic variant, maps it to relevant transcripts and proteins where appropriate, and returns the results as JSON.

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
| `-v`, `--variant` | The genomic variant description(s) to format. |

---

## Common Options

Commonly used command-line options include:

| Option | Description |
|--------|-------------|
| `-g`, `--genome` | Specify the reference genome assembly (e.g. `GRCh37` or `GRCh38`). |
| `-t`, `--select-transcripts` | Restrict the returned transcript representations. |
| `--transcript-model` | Select the transcript database (`refseq`, `ensembl` or `all`). |
| `--check-only` | Validate genomic HGVS syntax without transcript or protein mapping. |
| `-l`, `--liftover-level` | Control generation of equivalent genomic representations on alternate genome assemblies. |
| `-o`, `--output` | Write the results to a JSON file. |
| `--help` | Display the command help message. |

---

## Default Behaviour

Unless otherwise specified, VariantFormatter uses the following defaults.

| Setting | Default |
|---------|---------|
| Genome assembly | `GRCh38` |
| Transcript selection | `mane_select` |
| Transcript database | `refseq` |
| Genomic syntax checking only | Disabled |
| Liftover | Enabled |
| Output format | JSON |
| Output destination | Standard output (`stdout`) |

These defaults can be overridden using the appropriate command-line options described below.

---

## Supported Input Formats

VariantFormatter accepts **genomic variant descriptions** as input and generates corresponding transcript and protein descriptions where appropriate.

Supported input formats include:

- Genomic HGVS (`g.` notation) using supported genomic reference sequences, including `NC_`, `NT_` and `NW_` accessions.
- Pseudo-VCF chromosome-coordinate notation, for example:
  - `17-50198002-C-A`
  - `17:50198002:C:A`

Transcript (`c.` and `n.`), RNA (`r.`) and protein (`p.`) HGVS descriptions are **not accepted as VariantFormatter input**.

VariantFormatter operates from a genomic variant and maps that variant to relevant transcript and protein representations.

See the [Supported Input Formats](../reference/supported_inputs.md) guide for further details.

---

## Output Format

VariantFormatter returns structured JSON output.

For each submitted genomic variant, the output can include:

- the formatted genomic HGVS description;
- a pseudo-VCF representation;
- transcript HGVS descriptions;
- predicted protein HGVS descriptions;
- genomic representations on other genome assemblies when liftover is enabled;
- warnings and errors associated with genomic or transcript mapping.

Output can be written directly to the terminal or saved to a file using the `--output` option.

A detailed description of the JSON output is provided in the [Output Formats](../reference/output_formats.md) guide.

---

## Transcript Selection

VariantFormatter maps genomic variants to overlapping transcripts.

The `--select-transcripts` option controls which transcripts are included in the output.

Supported transcript selection strategies include:

- MANE Select transcripts;
- MANE Select and Plus Clinical transcripts;
- all transcripts overlapping the genomic variant at their latest version;
- all transcript versions;
- user-specified transcript lists.

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

### Format a pseudo-VCF variant

VariantFormatter accepts pseudo-VCF chromosome-coordinate notation.

Hyphen-delimited input can be supplied as:

```bash
variantformatter \
    --variant "17-50198002-C-A" \
    --genome GRCh38
```

Colon-delimited input can also be supplied:

```bash
variantformatter \
    --variant "17:50198002:C:A" \
    --genome GRCh38
```

---

### Format multiple variants

Multiple variants can be supplied by repeating the `--variant` argument:

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --variant "NC_000016.10:g.15738651_15738652inv" \
    --genome GRCh38
```

Multiple variants can also be supplied as a pipe-delimited value:

```bash
variantformatter \
    --variant "17-50198002-C-A|16-15738651-GT-AC" \
    --genome GRCh38
```

or as a JSON array:

```bash
variantformatter \
    --variant '["NC_000017.11:g.50198002C>A","NC_000016.10:g.15738651_15738652inv"]' \
    --genome GRCh38
```

Each variant is processed as part of the formatting request and returned in the JSON output.

---

## Selecting Transcripts

### Restrict the output to MANE Select transcripts

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts mane_select
```

---

### Restrict the output to MANE transcripts

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts mane
```

---

### Return all latest transcript versions

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts all
```

---

### Return all transcript versions

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts raw
```

---

### Restrict the output to a single specified transcript

The `--select-transcripts` option can also be supplied with an explicit transcript identifier.

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts "NM_000088.4"
```

Only the specified transcript is requested for the genomic variant.

---

### Restrict the output to multiple specified transcripts

Multiple transcript identifiers can be supplied as a pipe-delimited value:

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts "NM_000088.3|NM_000088.4"
```

or as a JSON array:

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --select-transcripts '["NM_000088.3","NM_000088.4"]'
```

RefSeq and Ensembl transcript identifiers must not be mixed in the same explicit transcript list.

---

## Transcript Models

The `--transcript-model` option controls which transcript database is used.

Available values are:

| Value | Description |
|-------|-------------|
| `refseq` | Use RefSeq transcripts. |
| `ensembl` | Use Ensembl transcripts. |
| `all` | Use both RefSeq and Ensembl transcripts. |

For example:

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --transcript-model ensembl
```

---

## Validate Genomic HGVS Syntax Only

The `--check-only` option validates the genomic HGVS description without generating transcript or protein representations.

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --check-only
```

This can be useful when only validation and genomic formatting are required.

---

## Liftover

VariantFormatter can generate genomic representations on another genome assembly.

The `--liftover-level` option controls this behaviour.

| Value | Description |
|-------|-------------|
| `true` | Perform full liftover. |
| `primary` | Perform liftover while excluding alternative scaffolds. |
| `false` | Disable liftover. |

Liftover is enabled by default.

For example:

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --liftover-level true
```

To disable liftover:

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --liftover-level false
```

---

## Format Variants from an Input File

VariantFormatter can read variants from a text or JSON file by prefixing the filename with `@`.

### Text file

Each non-empty line should contain one variant description. Lines beginning with `#` are ignored.

For example, `variants.txt` might contain:

```text
NC_000017.11:g.50198002C>A
NC_000016.10:g.15738651_15738652inv
17-50198002-C-A
```

Run:

```bash
variantformatter \
    --variant @variants.txt \
    --genome GRCh38
```

### JSON file

A JSON input file should contain an array of variant description strings.

For example, `variants.json`:

```json
[
    "NC_000017.11:g.50198002C>A",
    "NC_000016.10:g.15738651_15738652inv"
]
```

Run:

```bash
variantformatter \
    --variant @variants.json \
    --genome GRCh38
```

---

## Write Results to a JSON File

By default, VariantFormatter writes JSON to standard output.

Use `--output` to write the results to a file:

```bash
variantformatter \
    --variant "NC_000017.11:g.50198002C>A" \
    --genome GRCh38 \
    --output results.json
```

---

## Display Command Help

To display the complete command-line help:

```bash
variantformatter --help
```

---

## Common Errors

Common problems include:

- invalid genomic HGVS syntax;
- unsupported input types;
- unsupported reference sequences;
- a reference sequence that does not correspond to the selected genome build;
- invalid transcript selection;
- invalid transcript model selection;
- invalid pseudo-VCF input;
- an input coordinate outside the reference sequence;
- inability to connect to the VariantValidator databases;
- a missing or incorrect VariantValidator configuration file.

VariantFormatter returns explanatory warnings or errors where possible.

For a complete description of command-line error messages, exit codes and troubleshooting guidance, see the [Errors and Error Codes](../reference/errors_and_error_codes.md) guide.

---

## Related Documentation

- [VariantFormatter Python API](../python-api/variantformatter_python.md)
- [VariantValidator Command Line Interface](variantvalidator_cli.md)
- [Supported Input Formats](../reference/supported_inputs.md)
- [Output Formats](../reference/output_formats.md)
- [Transcript Selection](../reference/transcript_selection.md)
- [Errors and Error Codes](../reference/errors_and_error_codes.md)
