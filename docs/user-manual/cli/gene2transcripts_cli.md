# Gene2Transcripts Command Line Interface

The Gene2Transcripts Command Line Interface (CLI) provides a simple way to retrieve transcript information for genes,
transcript accessions and HGNC identifiers directly from the command line. It is suitable for interactive queries,
batch processing and integration into bioinformatics workflows.

The CLI is intended for users who wish to use Gene2Transcripts without writing Python code.

For users who are not familiar with command-line tools or Python programming:

- The [VariantValidator website](https://variantvalidator.org) provides a user-friendly interface for exploring transcript information.
- The [VariantValidator REST API](https://rest.variantvalidator.org) allows programmatic access without requiring local installation.

---

## Basic Usage

The simplest way to retrieve transcript information is to provide a gene symbol.

```bash
gene2transcripts \
    --gene COL1A1
```

Gene2Transcripts retrieves transcript information associated with the supplied query and returns the results as JSON.

---

## Command Syntax

```text
gene2transcripts [OPTIONS]
```

To display the complete list of available options:

```bash
gene2transcripts --help
```

---

## Required Arguments

The following argument is always required.

| Argument | Description |
|----------|-------------|
| `-q`, `--query` | A gene symbol, transcript accession, HGNC identifier, multiple queries, or an input file. |

---

## Common Options

Commonly used command-line options include:

| Option | Description |
|---------|-------------|
| `-g`, `--genome` | Specify the reference genome assembly (`GRCh37` or `GRCh38`). |
| `-t`, `--select-transcripts` | Restrict the returned transcript set. |
| `--transcript-model` | Select the transcript database (`refseq` or `ensembl`). |
| `-o`, `--output` | Write the results to a file. |
| `--no-web-searches` | Disable HGNC web lookups. |
| `--no-genomic-spans` | Omit genomic span information. |
| `--lovd-syntax-check` | Enable LOVD syntax checking. |
| `--help` | Display the command help message. |

---

## Default Behaviour

Unless otherwise specified, Gene2Transcripts uses the following defaults.

| Setting | Default                   |
|---------|---------------------------|
| Genome assembly | `GRCh38`                  |
| Transcript database | `refseq`                  |
| Transcript selection | All transcripts `all`      |         
| HGNC web lookups | Enabled                   |
| Genomic spans | Included                  |
| LOVD syntax checker | Disabled                  |
| Output format | JSON                      |
| Output destination | Standard output (`stdout`) |

These defaults can be overridden using the command-line options described above and demonstrated in the examples below.

---

## Supported Query Types

Gene2Transcripts accepts several different query types, including:

- Gene symbols
- Transcript accessions
- HGNC identifiers
- Multiple queries supplied as a JSON array
- Multiple queries supplied as a pipe-delimited list
- Text files containing one query per line
- JSON files containing an array of queries

See the [Supported Input Formats](../reference/supported_inputs.md) guide for additional examples.

---

## Output Formats

Gene2Transcripts returns JSON output.

Output can be written directly to the terminal or saved to a file using the `--output` option.

A detailed description of the output format is provided in the
[Output Formats](../reference/output_formats.md) guide.

---

## Transcript Selection

Gene2Transcripts supports the same transcript selection strategies as VariantValidator.

These include:

- MANE Select transcripts
- MANE Select and Plus Clinical transcripts
- All transcripts overlapping a genomic variant at their latest version
- All transcripts overlapping a genomic variant at all versions
- User-specified transcript lists

See the [Transcript Selection](../reference/transcript_selection.md) guide for complete details.

---

## Examples

### Retrieve transcripts for a gene symbol

```bash
gene2transcripts \
    --query COL1A1
```

---

### Retrieve transcripts using the Ensembl transcript database

```bash
gene2transcripts \
    --query COL1A1 \
    --transcript-model ensembl
```

---

### Query a RefSeq transcript accession

```bash
gene2transcripts \
    --query NM_000088.4
```

---

### Query an Ensembl transcript accession

```bash
gene2transcripts \
    --query ENST00000225964.10 \
    --transcript-model ensembl
```

---

### Query an HGNC identifier

```bash
gene2transcripts \
    --query HGNC:2197
```

---

### Query multiple entries using a JSON array

```bash
gene2transcripts \
    --query '["COL1A1","COL1A2"]'
```

Each query is processed independently and returned in the order supplied.

**Note:** RefSeq and Ensembl transcript accessions should not be mixed within the same query.

---

### Query multiple entries using a pipe-delimited list

```bash
gene2transcripts \
    --query "COL1A1|COL1A2|COL3A1"
```

---

### Restrict the output to MANE Select transcripts

```bash
gene2transcripts \
    --query COL1A1 \
    --select-transcripts mane_select
```

---

### Restrict the output to a single specified transcript

```bash
gene2transcripts \
    --query COL1A1 \
    --select-transcripts '["NM_000088.4"]'
```

---

### Restrict the output to multiple specified transcripts

```bash
gene2transcripts \
    --query COL1A1 \
    --select-transcripts '["NM_000088.3","NM_000088.4"]'
```

RefSeq and Ensembl transcript identifiers must **not** be mixed when using `--select-transcripts`.

---

### Specify the genome assembly

```bash
gene2transcripts \
    --query COL1A1 \
    --genome GRCh37
```

---

### Disable HGNC web lookups

```bash
gene2transcripts \
    --query COL1A1 \
    --no-web-searches
```

---

### Omit genomic span information

```bash
gene2transcripts \
    --query COL1A1 \
    --no-genomic-spans
```

---

### Enable LOVD syntax checking

```bash
gene2transcripts \
    --query COL1A1 \
    --lovd-syntax-check
```

---

### Write the results to a JSON file

```bash
gene2transcripts \
    --query COL1A1 \
    --output results.json
```

---

### Query from an input file

Each line of the input file should contain a single supported query.

```bash
gene2transcripts \
    --query queries.txt
```

---

### Query from an input file and write the results to a JSON file

```bash
gene2transcripts \
    --query queries.txt \
    --output results.json
```

---

### Display the command help

```bash
gene2transcripts --help
```

---

## Common Errors

Common problems include:

- Unknown gene symbol or transcript accession.
- Unsupported query format.
- Invalid transcript selection.
- Unable to connect to external services.
- Missing or incorrect configuration file.

Most errors include an explanatory message describing the cause of the problem.

For a complete description of command-line error messages, exit codes and troubleshooting guidance, see the
[Errors and Error Codes](../reference/errors_and_error_codes.md) guide.

---

## Related Documentation

- [Gene2Transcripts Python API](../python-api/gene2transcripts_python.md)
- [Supported Input Formats](../reference/supported_inputs.md)
- [Output Formats](../reference/output_formats.md)
- [Transcript Selection](../reference/transcript_selection.md)
- [Errors and Error Codes](../reference/errors_and_error_codes.md)