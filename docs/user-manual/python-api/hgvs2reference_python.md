# hgvs2reference Python API

The `hgvs2ref()` method retrieves the reference sequence corresponding to an HGVS sequence variant description.

It is intended for applications that require the reference sequence underlying a submitted variant, for example to extract the reference allele or surrounding sequence.

---

# Importing VariantValidator

```python
from VariantValidator import Validator
```

Create a validator instance:

```python
vv = Validator()
```

---

# Basic usage

Call `hgvs2ref()` with an HGVS sequence variant description.

```python
result = vv.hgvs2ref("NM_000088.4:c.589G>T")
```

The method returns a dictionary containing the retrieved reference sequence and associated metadata.

---

# Return value

The returned dictionary contains the following fields.

| Field | Description |
|-------|-------------|
| `variant` | The submitted HGVS sequence variant description. |
| `start_position` | The HGVS start position. |
| `end_position` | The HGVS end position. |
| `sequence` | The reference sequence corresponding to the variant coordinates. |
| `warning` | Any non-fatal warning generated during processing. |
| `error` | Error message if sequence retrieval failed. |

---

# Supported input

`hgvs2ref()` accepts the HGVS sequence variant types described in the [Supported Input Formats](../reference/supported_inputs.md) guide.

RNA (`r.`), protein (`p.`), and mitochondrial (`m.`) sequence variants are not currently supported.

Fully intronic transcript variants cannot currently be resolved because HGVS transcript descriptions do not explicitly define the genomic reference sequence used for transcript alignment. A future HGVS nomenclature update is expected to address this limitation.

---

# Example

```python
from VariantValidator import Validator

vv = Validator()

result = vv.hgvs2ref("NC_000017.11:g.50198002C>A")

print(result["sequence"])
```

---

# Errors and warnings

If the submitted variant cannot be parsed or the reference sequence cannot be retrieved, an error message is returned in the `error` field.

Warnings are returned in the `warning` field where sequence retrieval is only partially possible, such as transcript variants spanning intron boundaries.