# LOVD HGVS Syntax Checker REST API

The LOVD HGVS Syntax Checker REST API is developed and maintained by the LOVD project.

It provides programmatic access to the LOVD HGVS Syntax Checker, allowing applications to perform HGVS syntax checking using simple HTTP requests that return JSON responses.

VariantValidator uses this API as a fallback when a local installation of the LOVD HGVS Syntax Checker is unavailable.

---

# API Documentation

The official interactive API documentation is available via the Swagger interface:

[https://api.lovd.nl/swagger/](https://api.lovd.nl/swagger/)

The Swagger documentation describes:

- available endpoints;
- request parameters;
- response formats;
- error responses; and
- example requests.

As the API is maintained by the LOVD project, users should refer to the official documentation for the most up-to-date information.

---

# Returned Data

The LOVD HGVS Syntax Checker REST API returns JSON responses describing the syntax analysis of the submitted variant or gene symbol.

The precise structure of the returned data is documented in the official Swagger documentation.

---

# Relationship to VariantValidator

VariantValidator provides several interfaces that utilise the LOVD HGVS Syntax Checker.

These include:

- [LOVD HGVS Syntax Checker CLI](../user-manual/cli/lovd-hgvs-syntax-checker-cli.md)
- [LOVD HGVS Syntax Checker Python API](../user-manual/python-api/lovd-hgvs-syntax-checker-python.md)

The VariantValidator Python API automatically uses a locally installed LOVD HGVS Syntax Checker where available and falls back to the LOVD REST API when required.

---

# Further Reading

Official LOVD resources:

- REST API (Swagger): [https://api.lovd.nl/swagger/](https://api.lovd.nl/swagger/)
- Website: [https://lovd.nl/HGVS/](https://lovd.nl/HGVS/)
- GitHub repository: [https://github.com/LOVDnl/HGVS-syntax-checker](https://github.com/LOVDnl/HGVS-syntax-checker)

---

# How to cite

If the LOVD HGVS Syntax Checker contributes to your research or publication, please cite:

> Fokkema IF, Kroon M, López Hernández JA, Asscheman D, Lugtenburg I, Hoogenboom J, den Dunnen JT. **The LOVD3 platform: efficient genome-wide sharing of genetic variants.** *European Journal of Human Genetics* (2021). https://doi.org/10.1038/s41431-021-00959-x. :contentReference[oaicite:0]{index=0}

# How to cite VariantValidator

If you use VariantValidator in your research, please [cite the appropriate VariantValidator publication(s)](https://github.com/openvar/VariantValidator#cite-us).

