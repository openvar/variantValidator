<img src="../../static/img/logos/VV_logo.png" width="20%" />

# HGVS2Reference Command Line Interface

The HGVS2Reference Command Line Interface (CLI) provides a simple way to retrieve the reference sequence corresponding to an HGVS variant description directly from the command line.

The CLI supports genomic (`g.`), coding (`c.`) and non-coding transcript (`n.`) HGVS variants, including transcript variants with explicit genomic context (`NC_(NM_)` and `NC_(NR_)`).

It is intended for users who wish to retrieve reference sequence information without writing Python code.

For users who are not familiar with command-line tools or Python programming:

- The [VariantValidator website](https://variantvalidator.org) provides a user-friendly interface for validating HGVS variant descriptions.
- The [VariantValidator REST API](https://rest.variantvalidator.org) provides programmatic access to VariantValidator over HTTP.

## See also

- [HGVS2Reference Python API](../python-api/hgvs2reference_python.md) — Access HGVS2Reference directly from Python.
- [Supported Input Formats](../reference/supported_inputs.md) — Supported HGVS variant formats.
- [Reference Sequences](../reference/reference_sequences.md) — Understanding retrieved reference sequence data.
- [Errors and Error Codes](../reference/errors_and_error_codes.md) — Error messages and troubleshooting guidance.

---

# Basic Usage

The simplest way to retrieve a reference sequence is to provide an HGVS variant description.

```bash
hgvs2reference \
    --query NM_000546.6:c.215C>G
```

HGVS2Reference retrieves the reference sequence corresponding to the supplied HGVS variant and returns the results as JSON.

---

# Command Syntax

```text
hgvs2reference [OPTIONS]
```

To display the complete list of available options:

```bash
hgvs2reference --help
```

---

# Required Arguments

The following argument is always required.

| Argument | Description |
| --- | --- |
| `-q`, `--query` | A supported HGVS variant description. |

---

# Common Options

| Option | Description |
| --- | --- |
| `-o`, `--output` | Write the results to a JSON file. |
| `--log-level` | Set the console logging level. |
| `--help` | Display the command help message. |
| `--version` | Display the installed version. |

---

# Default Behaviour

Unless otherwise specified, HGVS2Reference uses the following defaults.

| Setting | Default |
| --- | --- |
| Output format | JSON |
| Output destination | Standard output (`stdout`) |
| Logging level | `WARNING` |

These defaults can be overridden using the command-line options described below.

---

# Supported Input Formats

HGVS2Reference accepts a single HGVS variant description.

Supported coordinate systems include:

- Genomic (`g.`)
- Coding DNA (`c.`)
- Non-coding transcript (`n.`)
- Transcript variants with explicit genomic context (`NC_(NM_)` and `NC_(NR_)`)

Intronic transcript variants are supported provided that genomic context is supplied.

Examples of supported queries include:

```text
NC_000017.11:g.7676594C>T

NM_000546.6:c.215C>G

NR_002196.3:n.145G>A

NC_000017.11(NM_000546.6):c.375+1G>A
```

---

# Output

HGVS2Reference returns JSON describing the retrieved reference sequence together with positional information.

Typical output includes:

- Submitted variant
- Start position
- End position
- Retrieved reference sequence
- Warning messages (if applicable)
- Error messages (if applicable)

Output can be written directly to the terminal or saved to a file using the `--output` option.

---

# Examples

## Retrieve the reference sequence for a genomic variant

```bash
hgvs2reference \
    --query NC_000017.11:g.7676594C>T
```

---

## Retrieve the reference sequence for a coding variant

```bash
hgvs2reference \
    --query NM_000546.6:c.215C>G
```

---

## Retrieve the reference sequence for a non-coding transcript variant

```bash
hgvs2reference \
    --query NR_002196.3:n.145G>A
```

---

## Retrieve the reference sequence for an intronic transcript variant

```bash
hgvs2reference \
    --query NC_000017.11(NM_000546.6):c.375+1G>A
```

Transcript variants with intronic coordinates require explicit genomic context.

---

## Write the results to a JSON file

```bash
hgvs2reference \
    --query NM_000546.6:c.215C>G \
    --output results.json
```

---

## Increase logging output

```bash
hgvs2reference \
    --query NM_000546.6:c.215C>G \
    --log-level INFO
```

---

## Display the command help

```bash
hgvs2reference --help
```

---

## Display the installed version

```bash
hgvs2reference --version
```

---

# Common Errors

Common problems include:

- Unsupported HGVS coordinate system.
- Invalid HGVS syntax.
- Intronic transcript variants supplied without genomic context.
- Reference sequence unavailable.
- Missing or incorrect VariantValidator configuration.

Most errors include an explanatory message describing the cause of the problem.

For a complete description of command-line error messages, exit codes and troubleshooting guidance, see the [Errors and Error Codes](../reference/errors_and_error_codes.md) guide.

---

# Getting Help

VariantValidator has been developed to support a wide range of users, from those new to HGVS nomenclature to experienced clinical scientists and bioinformaticians. If you encounter difficulties using the HGVS2Reference command-line interface or interpreting the returned results, you may find the following documentation helpful:

- [HGVS2Reference Python API](../python-api/hgvs2reference_python.md)
- [Supported Input Formats](../reference/supported_inputs.md)
- [Reference Sequences](../reference/reference_sequences.md)
- [Errors and Error Codes](../reference/errors_and_error_codes.md)

If you still require assistance, you can contact the VariantValidator team using our [contact form](https://variantvalidator.org/help/contact/).

Software bugs and feature requests can be reported through the [VariantValidator GitHub issue tracker](https://github.com/openvar/VariantValidator/issues).

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="../../static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="../../static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>
