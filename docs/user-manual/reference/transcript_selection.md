# Transcript Selection

Many genes produce multiple transcript isoforms, each with a different transcript reference sequence. Unless a transcript is explicitly specified in the submitted variant description, the VariantValidator software suite must determine which transcript(s) should be used during validation, formatting or annotation.

Transcript selection controls which transcript reference sequences are considered when reporting transcript-level variants.

This guide describes the transcript selection strategies shared throughout the VariantValidator software suite.

---

# Why transcript selection matters

A single genomic variant may overlap multiple transcripts. Consequently, different transcript reference sequences may produce different HGVS transcript descriptions, protein consequences, exon numbering and coding positions.

For example, a single genomic variant may produce transcript descriptions for:

```text
NM_000088.4
NM_000088.5
NM_001278074.2
ENST00000371953.9
```

Transcript selection determines which of these transcript reference sequences are returned.

---

# Specifying an explicit transcript

The most precise method is to specify the desired transcript reference sequence directly.

Examples include:

```text
NM_000088.4
```

```text
NR_110010.2
```

```text
ENST00000371953.9
```

When an explicit transcript accession is supplied, only that transcript is processed.

---

# Transcript selection terms

Instead of specifying an individual transcript accession, VariantValidator provides several predefined transcript selection strategies.

## `mane`

Returns a single MANE transcript.

Where available, the MANE Select transcript is returned. If no MANE Select transcript exists but a MANE Plus Clinical transcript is available, the MANE Plus Clinical transcript is returned.

This option always attempts to return a single representative MANE transcript.

---

## `mane_select`

Returns only the MANE Select transcript.

If a MANE Select transcript is not available for the submitted gene, no transcript is returned.

---

## `select`

Returns the preferred transcript for the selected transcript collection.

For RefSeq transcript alignments this corresponds to the RefSeq Select transcript.

For Ensembl transcript alignments this corresponds to the Ensembl Select transcript where available.

---

## `all`

Returns every supported transcript associated with the submitted variant.

This option provides the most comprehensive transcript annotation and is useful for research, annotation pipelines and transcript comparison.

---

## `raw`

Returns all available transcript mappings without transcript filtering.

This option is primarily intended for advanced users and downstream software.

---

# RefSeq and Ensembl transcript collections

The transcript identifiers returned depend upon the selected transcript alignment method.

RefSeq transcript alignments return RefSeq transcript identifiers, for example:

```text
NM_000088.4
NR_110010.2
```

Ensembl transcript alignments return Ensembl transcript identifiers, for example:

```text
ENST00000371953.9
```

The selected transcript alignment method therefore determines both the available transcript accessions and the transcript selection terms that are applied.

---

# Gene symbol input

Transcript selection is particularly important when a gene symbol is supplied instead of a transcript reference sequence.

For example:

```text
COL1A1:c.589G>T
```

Because a gene may have multiple transcript isoforms, a gene symbol alone does not uniquely identify the intended transcript.

If a transcript selection strategy is supplied, VariantValidator automatically determines the appropriate transcript reference sequence and continues processing.

If no transcript selection strategy is supplied, VariantValidator cannot determine which transcript was intended. Validation stops and returns the compatible transcript reference sequences that may be used to resubmit the variant.

---

# Recommendations

For routine clinical interpretation and reporting, `mane_select` is generally recommended where available.

When reproducing published work or validating a known transcript variant, the transcript accession should be specified explicitly.

Use `all` when comprehensive transcript annotation is required.

Use `raw` only when unfiltered transcript mappings are required for downstream processing.
