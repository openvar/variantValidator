<img src="../../static/img/logos/VV_logo.png" width="20%" />

# LOVD HGVS Syntax Checker Command-Line Interface

The LOVD HGVS Syntax Checker CLI provides access to the **LOVD HGVS Syntax Checker** from the command line using the VariantValidator Python package.

The CLI performs syntax checking of HGVS sequence variant descriptions and gene symbols using the locally installed LOVD HGVS Syntax Checker. If the local installation is unavailable, the CLI automatically falls back to the LOVD web API where supported.

Unlike the VariantValidator validation engine, the LOVD HGVS Syntax Checker does **not** require the VariantValidator databases or configuration to be installed.

The official LOVD HGVS Syntax Checker project is maintained by the LOVD team:

- GitHub: https://github.com/LOVDnl/HGVS-syntax-checker
- Documentation: https://lovd.nl/HGVS/

---

# Installation

The LOVD HGVS Syntax Checker is distributed as part of the VariantValidator Python package.

Create and activate the recommended conda environment before installing VariantValidator.

```bash
conda env create -f environment.yml
conda activate vvenv
```

Install VariantValidator into the environment.

```bash
pip install ./packaging/variantvalidator
```

The LOVD HGVS Syntax Checker itself must then be downloaded and configured locally.

```bash
python -m VariantValidator.bin.setup_lovd_syntax_checker
```

This command:

- downloads the latest LOVD HGVS Syntax Checker release;
- downloads the required support files;
- builds the local syntax checker cache; and
- prepares the command-line interface for use.

The LOVD project is actively developed. Running this command periodically is recommended to ensure that the latest version of the syntax checker and supporting reference data are installed.

> **Note**
>
> The LOVD HGVS Syntax Checker does **not** require the VariantValidator databases or configuration file.
>
> If you wish to use the full VariantValidator validation engine, including transcript mapping and variant validation, see the [Installation Guide](../../installation/installation.md).

---

# Running the LOVD HGVS Syntax Checker

The command-line tool is:

```bash
lovd-hgvs-syntax-checker
```

Results are written as formatted JSON.

---

## Checking an HGVS variant description

Submit a single HGVS sequence variant description using the `--query` option.

```bash
lovd-hgvs-syntax-checker \
    --query "NM_000059.4:c.7790G>A"
```

or

```bash
lovd-hgvs-syntax-checker \
    -q "NM_000059.4:c.7790G>A"
```

---

## Checking a genomic variant

Genomic HGVS descriptions are handled in the same way.

```bash
lovd-hgvs-syntax-checker \
    -q "NC_000013.11:g.32355250G>A"
```

---

## Checking a gene symbol

Gene symbols are checked using the `--gene` option.

```bash
lovd-hgvs-syntax-checker \
    -q BRCA2 \
    --gene
```

---

## Writing results to a file

JSON output may be written directly to a file.

```bash
lovd-hgvs-syntax-checker \
    -q "NM_000059.4:c.7790G>A" \
    --output result.json
```

---

# Output

The LOVD HGVS Syntax Checker always returns formatted JSON.

The returned information originates from the LOVD HGVS Syntax Checker and includes syntax validation results together with metadata describing the version of the checker used to perform the analysis.

---

# Relationship to VariantValidator

The LOVD HGVS Syntax Checker performs **syntax checking only**.

If you require:

- HGVS validation;
- transcript mapping;
- genomic projection;
- protein consequence prediction;
- transcript selection; or
- variant formatting,

use the VariantValidator command-line interface instead.

See:

- [VariantValidator CLI](variantvalidator_cli.md)

---

# Further Reading

The official documentation for the LOVD HGVS Syntax Checker is maintained by the LOVD project.

- Documentation: https://lovd.nl/HGVS/
- GitHub repository: https://github.com/LOVDnl/HGVS-syntax-checker

VariantValidator documentation:

- [VariantValidator CLI](variantvalidator_cli.md)
- [Python API](../python-api/variantvalidator_python.md)
- [Supported Input Formats](../reference/supported_inputs.md)
- [Errors and Error Codes](../reference/errors_and_error_codes.md)

---

# How to cite

If the LOVD HGVS Syntax Checker contributes to your research or publication, please cite:

> Fokkema IF, Kroon M, López Hernández JA, Asscheman D, Lugtenburg I, Hoogenboom J, den Dunnen JT. **The LOVD3 platform: efficient genome-wide sharing of genetic variants.** *European Journal of Human Genetics* (2021). https://doi.org/10.1038/s41431-021-00959-x. :contentReference[oaicite:0]{index=0}

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="../static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="../static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>
