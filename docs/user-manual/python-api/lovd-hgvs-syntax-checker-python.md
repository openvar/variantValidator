# LOVD HGVS Syntax Checker Python API

The LOVD HGVS Syntax Checker Python API provides convenient access to the **LOVD HGVS Syntax Checker** from within Python applications.

The API performs syntax checking of HGVS sequence variant descriptions and gene symbols using the locally installed LOVD HGVS Syntax Checker. If the local installation is unavailable, the API automatically falls back to the LOVD web API where supported.

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
pip install .
```

The LOVD HGVS Syntax Checker itself must then be downloaded and configured locally.

```bash
python -m VariantValidator.bin.setup_lovd_syntax_checker
```

This command:

- downloads the latest LOVD HGVS Syntax Checker release;
- downloads the required support files;
- builds the local syntax checker cache; and
- prepares the Python API for use.

The LOVD project is actively developed. Running this command periodically is recommended to ensure that the latest version of the syntax checker and supporting reference data are installed.

> **Note**
>
> The LOVD HGVS Syntax Checker does **not** require the VariantValidator databases or configuration file.
>
> If you wish to use the full VariantValidator validation engine, including transcript mapping and variant validation, see the [Installation Guide](../../installation/installation.md).

---

# Importing the API

Import the LOVD HGVS Syntax Checker module.

```python
from VariantValidator.modules import lovd_api
```

The primary entry point is:

```python
lovd_api.lovd_syntax_check()
```

This function automatically performs syntax checking using the locally installed LOVD HGVS Syntax Checker and falls back to the LOVD web API if necessary.

---

# Checking an HGVS variant description

Pass an HGVS sequence variant description directly to the API.

```python
from VariantValidator.modules import lovd_api

result = lovd_api.lovd_syntax_check(
    "NM_000059.4:c.7790G>A"
)

print(result)
```

---

# Checking a genomic variant

Genomic HGVS descriptions are handled in the same way.

```python
from VariantValidator.modules import lovd_api

result = lovd_api.lovd_syntax_check(
    "NC_000013.11:g.32355250G>A"
)
```

---

# Checking a gene symbol

Gene symbols are checked by enabling the `is_a_gene` option.

```python
from VariantValidator.modules import lovd_api

result = lovd_api.lovd_syntax_check(
    "BRCA2",
    is_a_gene=True,
)
```

---

# Returned data

The API returns a Python dictionary.

For successful queries the returned dictionary contains the syntax checking results together with metadata describing the version of the LOVD HGVS Syntax Checker used to generate the results.

For example:

```python
{
    "data": [...],
    "url": "...",
    "version": "3.x.x"
}
```

If an error occurs, the returned dictionary contains a `lovd_api_error` entry describing the problem.

---

# Local execution and automatic fallback

The API first attempts to execute the locally installed LOVD HGVS Syntax Checker.

If the local installation is unavailable or an unexpected error occurs, the API automatically falls back to the LOVD web API where supported.

This behaviour is automatic and requires no additional configuration.

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

use the VariantValidator Python API instead.

See:

- [VariantValidator Python API](variantvalidator_python.md)

---

# Further Reading

The official documentation for the LOVD HGVS Syntax Checker is maintained by the LOVD project.

- Documentation: https://lovd.nl/HGVS/
- GitHub repository: https://github.com/LOVDnl/HGVS-syntax-checker

VariantValidator documentation:

- [VariantValidator Python API](variantvalidator_python.md)
- [Supported Input Formats](../reference/supported_inputs.md)
- [Errors and Error Codes](../reference/errors_and_error_codes.md)
