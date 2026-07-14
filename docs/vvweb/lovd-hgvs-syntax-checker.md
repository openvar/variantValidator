# LOVD HGVS Syntax Checker

The LOVD HGVS Syntax Checker web interface is developed and maintained by the LOVD project.

It provides an interactive web interface for checking the syntax of HGVS sequence variant descriptions and gene symbols without requiring local software installation. The tool is intended for users wishing to perform individual syntax checks through a web browser.

VariantValidator incorporates the LOVD HGVS Syntax Checker within several of its own interfaces, but the web interface described here is provided directly by the LOVD project.

---

# Accessing the Web Interface

The LOVD HGVS Syntax Checker is available at:

https://lovd.nl/HGVS/

No installation or account is required.

---

# Using the Web Interface

Enter either:

- an HGVS sequence variant description; or
- a gene symbol.

The LOVD HGVS Syntax Checker analyses the submitted query and reports any syntax errors, warnings or recommendations according to the HGVS nomenclature.

For guidance on interpreting the results and the behaviour of the syntax checker, refer to the official LOVD documentation.

---

# Relationship to VariantValidator

VariantValidator uses the LOVD HGVS Syntax Checker to provide optional syntax checking within several tools.

VariantValidator interfaces include:

- [LOVD HGVS Syntax Checker CLI](../user-manual/cli/lovd-hgvs-syntax-checker-cli.md)
- [LOVD HGVS Syntax Checker Python API](../user-manual/python-api/lovd-hgvs-syntax-checker-python.md)
- [LOVD HGVS Syntax Checker REST API](../rest-vv/lovd-hgvs-syntax-checker.md)

If you require full HGVS validation, transcript mapping, variant formatting or genomic projection, use the VariantValidator tools instead.

---

# Documentation

The LOVD HGVS Syntax Checker is documented by the LOVD project.

Official resources include:

- Website: https://lovd.nl/HGVS/
- REST API (Swagger): https://api.lovd.nl/swagger/
- GitHub repository: https://github.com/LOVDnl/HGVS-syntax-checker

---

# How to cite

If the LOVD HGVS Syntax Checker contributes to your research or publication, please cite:

> Fokkema IF, Kroon M, López Hernández JA, Asscheman D, Lugtenburg I, Hoogenboom J, den Dunnen JT. **The LOVD3 platform: efficient genome-wide sharing of genetic variants.** *European Journal of Human Genetics* (2021). https://doi.org/10.1038/s41431-021-00959-x. :contentReference[oaicite:0]{index=0}

If VariantValidator also contributes to your work through its integration of the LOVD HGVS Syntax Checker or other VariantValidator functionality, please also cite the appropriate VariantValidator publication(s).

The latest VariantValidator citation information is available in the project repository:

https://github.com/openvar/VariantValidator#cite-us

