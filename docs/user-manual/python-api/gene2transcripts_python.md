# Gene2Transcripts Python API

The Gene2Transcripts Python API provides direct access to the Gene2Transcripts engine from within Python. It is suitable
for integrating transcript retrieval into bioinformatics pipelines, analysis workflows, web applications and custom
software.

The Python API offers access to the same functionality as the Gene2Transcripts command-line interface while providing a
programmatic interface for automated processing and downstream analysis.

For users who prefer not to write Python code:

- The [VariantValidator website](https://variantvalidator.org) provides a user-friendly interface for exploring transcript information.
- The [VariantValidator REST API](https://rest.variantvalidator.org) allows programmatic access without requiring local installation.
- The [Gene2Transcripts Command Line Interface](../cli/gene2transcripts_cli.md) provides a command-line interface for retrieving transcript information locally.

---

## Basic Usage

Begin by importing the VariantValidator package and creating a `Validator` object.

```python
import json
import VariantValidator

vval = VariantValidator.Validator()
```

The `Validator` object manages access to the Gene2Transcripts engine and can be reused for multiple queries within the
same Python session.

Once a `Validator` object has been created, transcript information can be retrieved using the `gene2transcripts()`
method.

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(
    query="COL1A1"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

This retrieves transcript information for the supplied query and returns the results as a Python dictionary, which can
then be processed directly or converted to JSON.

---

## Method Signature

Transcript information is retrieved using the `gene2transcripts()` method of the `Validator` object.

```python
gene2transcripts(
    query=None,
    validator=None,
    bypass_web_searches=False,
    select_transcripts="all",
    transcript_set="refseq",
    genome_build="GRCh38",
    bypass_genomic_spans=False,
    lovd_syntax_check=False
)
```

---

## Required Arguments

The following argument is required when calling the `gene2transcripts()` method.

| Argument | Description |
|----------|-------------|
| `query` | A gene symbol, transcript accession, HGNC identifier, multiple queries, or an input file. |

---

## Optional Arguments

The following optional arguments control transcript retrieval behaviour.

| Argument | Default | Description                                                              |
|----------|---------|--------------------------------------------------------------------------|
| `bypass_web_searches` | `False` | Disable HGNC web lookups.                                                |
| `select_transcripts` | `"all"` | Transcript selection strategy or JSON array of transcript identifiers.   |
| `transcript_set` | `"refseq"` | Transcript database (`refseq` or `ensembl`).                             |
| `genome_build` | `"GRCh38"` | Return transcripts for a specific genome build.                          |
| `bypass_genomic_spans` | `False` | Omit genomic span and alignment information.                             |
| `lovd_syntax_check` | `False` | Enable LOVD HGVS syntax checking.                                        |

---

## Default Behaviour

Unless otherwise specified, the Gene2Transcripts Python API uses the following defaults.

| Setting | Default |
|---------|---------|
| Genome assembly | `GRCh38` |
| Transcript database | `refseq` |
| Transcript selection | `all` (all transcripts) |
| HGNC web lookups | Enabled |
| Genomic spans | Included |
| LOVD syntax checker | Disabled |
| Return type | Python `dict` |

---

## Supported Query Types

The Gene2Transcripts Python API accepts the same query types as the command-line interface.

Supported queries include:

- Gene symbols
- Transcript accessions
- HGNC identifiers
- Multiple queries supplied as a JSON array
- Multiple queries supplied as a pipe-delimited list
- Text files containing one query per line
- JSON files containing an array of queries

See the [Supported Input Formats](../reference/supported_inputs.md) guide for a complete description of supported input formats and examples.

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

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(query="COL1A1")

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Retrieve transcripts using the Ensembl transcript database

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(
    query="COL1A1",
    transcript_set="ensembl"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Query a RefSeq transcript accession

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(query="NM_000088.4")

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Query an Ensembl transcript accession

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(
    query="ENST00000225964.10",
    transcript_set="ensembl"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Query an HGNC identifier

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(query="HGNC:2197")

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Query multiple entries

Gene2Transcripts accepts multiple queries using a JSON array.

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(
    query='["COL1A1","COL1A2"]'
)

print(json.dumps(results, indent=4, sort_keys=True))
```

Each query is processed independently and returned in the order supplied.

**Note:** RefSeq and Ensembl transcript accessions should not be mixed within the same query.

---

### Using `select_transcripts`

The `select_transcripts` argument restricts the transcripts returned for a query.

RefSeq and Ensembl transcript identifiers must **not** be mixed when using `select_transcripts`.

---

### Restrict the output to MANE Select transcripts

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(
    query="COL1A1",
    select_transcripts="mane_select"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Restrict the output to a single specified transcript

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(
    query="COL1A1",
    select_transcripts='["NM_000088.4"]'
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Restrict the output to multiple specified transcripts

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(
    query="COL1A1",
    select_transcripts='["NM_000088.3","NM_000088.4"]'
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Specify the genome assembly

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(
    query="COL1A1",
    genome_build="GRCh37"
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Disable HGNC web lookups

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(
    query="COL1A1",
    bypass_web_searches=True
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Omit genomic span information

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(
    query="COL1A1",
    bypass_genomic_spans=True
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Enable LOVD syntax checking

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(
    query="COL1A1",
    lovd_syntax_check=True
)

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Query from an input file

Each line of the input file should contain a single supported query.

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(query="queries.txt")

print(json.dumps(results, indent=4, sort_keys=True))
```

---

### Write the results to a JSON file

```python
import json
import VariantValidator

vval = VariantValidator.Validator()

results = vval.gene2transcripts(query="COL1A1")

with open("results.json", "w") as fh:
    json.dump(results, fh, indent=4, sort_keys=True)
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

- [Gene2Transcripts Command Line Interface](../cli/gene2transcripts_cli.md)
- [Supported Input Formats](../reference/supported_inputs.md)
- [Output Formats](../reference/output_formats.md)
- [Transcript Selection](../reference/transcript_selection.md)
- [Errors and Error Codes](../reference/errors_and_error_codes.md)