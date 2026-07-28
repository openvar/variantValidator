<img src="../../static/img/logos/VV_logo.png" width="20%" />

# Supported Input Formats

This guide describes the input formats accepted by tools within the VariantValidator software suite.

The following tools accept different classes of input depending on their purpose:

- [VariantValidator](#variantvalidator)
- [VariantFormatter](#variantformatter)
- [gene2transcripts](#gene2transcripts)
- [hgvs2reference](#hgvs2reference)

VariantValidator accepts a broad range of sequence variant descriptions, including fully compliant HGVS descriptions, commonly encountered non-HGVS representations, legacy formats and common formatting mistakes.

VariantFormatter uses the same genomic input processing infrastructure and therefore accepts the same supported **genomic input formats** as VariantValidator. However, VariantFormatter is specifically a genomic-to-transcript/protein formatting tool and does not accept transcript, RNA or protein variants as its starting input.

Where a particular tool imposes additional restrictions, these are described in the relevant section below.

## See also

- [Transcript Selection](transcript_selection.md) — Transcript selection strategies such as `mane_select`, `mane`, `select` and `all`.
- [Output Formats](output_formats.md) — Results returned by each tool.
- [Errors and Error Codes](errors_and_error_codes.md) — Validation errors, warnings and informational messages.
- [VariantValidator Python API](../python-api/variantvalidator_python.md) — Validate variants directly from Python.
- [VariantFormatter Python API](../python-api/variantformatter_python.md) — Format genomic variants directly from Python.

---

# Scope of Supported Input Formats

The examples described in this guide are representative rather than exhaustive.

VariantValidator has been designed to recognise and interpret a wide range of real-world sequence variant descriptions. In addition to fully compliant HGVS sequence variant descriptions, it accepts many commonly encountered non-HGVS formats, legacy representations and common formatting mistakes originating from clinical laboratories, research pipelines, databases and published literature.

Where possible, non-standard inputs are recognised automatically and converted into valid HGVS sequence variant descriptions before validation. When an unambiguous correction cannot be made, informative warnings or error messages are returned to help the user generate a valid description.

VariantFormatter shares VariantValidator's genomic input processing and therefore accepts the same supported genomic representations. Once a genomic input has been interpreted, VariantFormatter maps the genomic variant to corresponding transcript and protein representations where appropriate.

Support for additional input formats continues to evolve as new real-world examples are encountered. If you regularly encounter a sequence variant description that is not currently recognised, we encourage you to submit a feature request through the project's contact page so that support can be considered for a future release.

---

# VariantValidator

VariantValidator accepts HGVS sequence variant descriptions at genomic, transcript, RNA, protein and mitochondrial levels, together with a broad range of non-HGVS and recoverable input formats.

---

## HGVS Sequence Variant Types

VariantValidator supports the HGVS sequence variant types defined by the Human Genome Variation Society (HGVS).

### Genomic variants (`g.`)

Genomic variants describe sequence changes relative to a genomic reference sequence.

Example:

```text
NC_000017.11:g.50198002C>A
```

---

### Coding DNA variants (`c.`)

Coding DNA variants describe sequence changes relative to a coding transcript.

Example:

```text
NM_000088.4:c.589G>T
```

Coding variants are validated directly against the specified transcript reference sequence.

---

### Non-coding transcript variants (`n.`)

Non-coding transcript variants are recognised and validated.

Example:

```text
NR_023343.1:n.245G>A
```

---

### RNA variants (`r.`)

RNA sequence variants are recognised and validated.

Example:

```text
NM_000088.4:r.589g>u
```

RNA variants must use the IUPAC RNA alphabet, for example **U** rather than **T**. Where possible, informative error messages are returned for incorrect RNA syntax.

---

### Protein variants (`p.`)

Protein sequence variants are accepted.

Example:

```text
NP_000079.2:p.Gly197Val
```

Protein-level validation is necessarily more limited than genomic or transcript-level validation because multiple nucleotide variants can produce the same protein consequence. Whenever possible, genomic or transcript variants should be used.

---

### Mitochondrial variants (`m.`)

Mitochondrial variants are recognised and validated.

Example:

```text
NC_012920.1:m.3243A>G
```

The software automatically recognises mitochondrial reference sequences and checks that the appropriate HGVS sequence type (`m.`) is used.

---

### Allele descriptions

VariantValidator supports HGVS allele descriptions containing one or more sequence variants.

Examples include:

```text
NM_000088.4:c.[589G>T;642+1G>A]
```

```text
NC_000017.11:g.[50198002C>A;50198015del]
```

Allele descriptions are recognised automatically and decomposed into their component variants for validation. Each variant is validated independently, and guidance is provided to assist reconstruction of the complete allele description where appropriate.

---

### Expanded repeat descriptions

VariantValidator supports HGVS expanded repeat descriptions.

Example:

```text
NC_000004.12:g.3074877CAG[42]
```

Expanded repeat descriptions are interpreted according to the current HGVS recommendations. Where appropriate, equivalent sequence variant descriptions are generated for downstream validation and normalisation.

---

### Uncertain and fuzzy positions

VariantValidator supports HGVS notation describing uncertain or imprecisely defined variant locations.

Examples include:

```text
NC_000017.11:g.(50198000_50198005)_(50198020_50198025)del
```

```text
NC_000017.11:g.(?_50198002)_(50198020_?)del
```

Uncertain intervals, fuzzy breakpoints and unknown variant boundaries are recognised where supported by the HGVS recommendations.

If a submitted description cannot be interpreted unambiguously, VariantValidator returns an informative validation message describing the problem.

---

### LRG reference sequences

Legacy LRG reference sequences are recognised.

Examples include:

```text
LRG_199:c.589G>T
```

```text
LRG_199t1:c.589G>T
```

Where appropriate, LRG identifiers are converted to their equivalent RefSeq reference sequences before further processing.

---

### Intronic and compound reference sequence variants

VariantValidator supports intronic HGVS sequence variant descriptions.

Examples include:

```text
NM_000088.4:c.589+1G>T
```

```text
NM_000088.4:c.690-2A>G
```

HGVS also permits transcript variants to specify the genomic reference sequence used for transcript alignment.

Examples include:

```text
NG_007400.1(NM_000088.4):c.589+1G>T
```

```text
NC_000017.11(NM_000088.4):c.589G>T
```

```text
NW_012345678.9(NM_000088.4):c.589G>T
```

```text
NT_012345678.9(NM_000088.4):c.589G>T
```

Compound reference sequence descriptions define the genomic sequence against which the transcript is aligned. This allows transcript variants to be interpreted relative to a RefSeqGene record, chromosome, scaffold or other supported genomic reference sequence.

VariantValidator recognises these descriptions and performs transcript mapping using the specified genomic alignment.

---

### Predicted variants

Predicted sequence variant descriptions are supported.

Examples include:

```text
NM_000088.4:r.(589G>T)
```

```text
NP_000079.2:p.(Gly197Val)
```

Predicted variants are recognised according to the HGVS recommendations and are retained as predicted descriptions throughout processing.

---

## Common Non-HGVS Input Formats

In addition to fully compliant HGVS sequence variant descriptions, VariantValidator accepts commonly encountered non-HGVS variant representations.

These formats frequently originate from variant calling pipelines, databases, spreadsheets, clinical reports, publications and legacy software.

Where possible, VariantValidator automatically converts these representations into valid HGVS sequence variant descriptions before validation.

---

### Gene symbols used as reference sequence identifiers

Gene symbols are frequently, but incorrectly, used in place of transcript reference sequence identifiers.

For example:

```text
COL1A1:c.589G>T
```

This is **not valid HGVS nomenclature**, because HGVS requires a reference sequence identifier rather than a gene symbol.

If a transcript selection strategy is supplied, VariantValidator can substitute an appropriate transcript reference sequence and continue validation.

If the intended transcript cannot be determined unambiguously, validation stops and guidance is returned.

---

### Pseudo-VCF chromosome notation

VariantValidator accepts simplified chromosome-coordinate formats commonly used in spreadsheets, databases and analysis pipelines.

Examples include:

```text
17-50198002-C-A
```

```text
17:50198002:C:A
```

```text
chr17-50198002-C-A
```

```text
chr17:50198002:C:A
```

These descriptions are automatically recognised and converted into the corresponding genomic HGVS sequence variant description before validation.

Pseudo-VCF descriptions containing multiple alternate alleles are also recognised.

For example:

```text
17-50198002-C-A,G,T
```

```text
17:50198002:C:A,G,T
```

Each alternate allele is decomposed into an independent sequence variant description before conversion to HGVS.

---

### Genome assembly prefixes

Some pipelines include the genome assembly as part of the chromosome-coordinate description.

Examples include:

```text
GRCh38-17-50198002-C-A
```

```text
GRCh38:17:50198002:C:A
```

```text
hg19-17-50198002-C-A
```

```text
hg38:17:50198002:C:A
```

VariantValidator recognises these formats, extracts the genome assembly, converts the description into HGVS format and validates the resulting variant.

If the embedded genome assembly conflicts with the selected genome assembly, an informative error is returned.

---

### VCF/HGVS hybrid formats

VariantValidator accepts hybrid formats that combine HGVS reference sequence identifiers with VCF-style coordinate or allele notation.

Examples include:

```text
NC_000017.11:50198002:C:A
```

```text
NC_000017.11-50198002-C-A
```

```text
NC_000017.11:g.50198002:C:A
```

```text
NC_000017.11:g.50198002-C-A
```

```text
NC_000017.11(GRCh38):g.50198002C>A
```

```text
NC_000017.11(hg38):g.50198002C>A
```

```text
Chr17(GRCh38):g.50198002C>A
```

```text
Chr17(hg19):g.48275363C>A
```

These representations are recognised, converted into genomic HGVS and then validated.

---

### Chromosome identifiers

Chromosome identifiers are recognised in commonly used forms, including:

```text
17
```

```text
chr17
```

```text
Chr17
```

Where a chromosome identifier forms part of a genomic variant description, it is mapped to the appropriate genomic reference sequence accession for the selected genome assembly.

---

### Variant Call Format (VCF)

VariantValidator accepts Variant Call Format representations and converts the chromosome, position, reference allele and alternate allele into the corresponding genomic HGVS sequence variant description.

For example:

```text
17    50198002    .    C    A
```

Multi-allelic representations can be decomposed into individual variants and processed independently.

---

## Additional Supported Formats

The examples presented above illustrate commonly encountered HGVS and non-HGVS input formats accepted by VariantValidator. They are **not** intended to be an exhaustive list of every supported syntax.

VariantValidator has been developed to recognise the broad range of sequence variant descriptions encountered in real-world clinical and research workflows.

Where possible, VariantValidator automatically converts recoverable representations into valid HGVS sequence variant descriptions before validation. When an unambiguous correction cannot be made, the software returns informative warnings or errors describing the problem.

Together with the integrated LOVD HGVS Syntax Checker, VariantValidator supports a broad range of HGVS syntax, legacy representations and commonly encountered input formats.

If you encounter a sequence variant description that is not recognised, or would like support for an additional input format, please submit a feature request through the project's contact page.

---

# VariantFormatter

VariantFormatter is designed to map **genomic variants** to corresponding transcript and protein representations for automated bioinformatics pipelines, software integration and other programmatic workflows.

VariantFormatter accepts the same supported **genomic input formats** as VariantValidator. This includes compliant genomic HGVS descriptions as well as genomic non-HGVS, pseudo-VCF, VCF-style, hybrid and recoverable representations recognised by the VariantValidator genomic input processing pipeline.

The important distinction between VariantValidator and VariantFormatter is therefore not the range of supported genomic syntax.

Instead:

- **VariantValidator** accepts genomic, transcript, RNA, protein and other supported sequence variant types and performs comprehensive validation and recovery.
- **VariantFormatter** accepts genomic variant input and maps that genomic variant to transcript and protein representations.

Transcript (`c.` and `n.`), RNA (`r.`) and protein (`p.`) descriptions are therefore not accepted as starting inputs to VariantFormatter.

---

## HGVS Genomic Sequence Variants

Genomic HGVS (`g.`) descriptions are accepted using supported genomic reference sequences.

Examples include:

```text
NC_000017.11:g.50198002C>A
```

```text
NT_187361.1:g.1000A>G
```

```text
NW_012345678.9:g.1000A>G
```

Accepted genomic descriptions are validated and normalised before being mapped to relevant transcript and protein representations.

---

## Pseudo-VCF Chromosome Notation

VariantFormatter accepts the same supported genomic pseudo-VCF representations as VariantValidator.

Examples include:

```text
17-50198002-C-A
```

```text
17:50198002:C:A
```

```text
chr17-50198002-C-A
```

```text
chr17:50198002:C:A
```

These descriptions are converted into genomic HGVS before transcript mapping.

---

## Genome Assembly Prefixes

Genomic pseudo-VCF descriptions may include an assembly identifier.

Examples include:

```text
GRCh38-17-50198002-C-A
```

```text
GRCh38:17:50198002:C:A
```

```text
hg19-17-50198002-C-A
```

```text
hg38:17:50198002:C:A
```

The assembly information is interpreted during genomic input processing.

---

## Genomic VCF/HGVS Hybrid Formats

VariantFormatter accepts the genomic hybrid representations recognised by VariantValidator.

Examples include:

```text
NC_000017.11:50198002:C:A
```

```text
NC_000017.11-50198002-C-A
```

```text
NC_000017.11:g.50198002:C:A
```

```text
NC_000017.11:g.50198002-C-A
```

```text
NC_000017.11(GRCh38):g.50198002C>A
```

```text
NC_000017.11(hg38):g.50198002C>A
```

```text
Chr17(GRCh38):g.50198002C>A
```

```text
Chr17(hg19):g.48275363C>A
```

These descriptions are resolved to an appropriate genomic reference sequence and converted into genomic HGVS before formatting continues.

---

## Multiple Alternate Alleles

Genomic pseudo-VCF input containing multiple alternate alleles is recognised.

Examples include:

```text
17-50198002-C-A,G,T
```

```text
17:50198002:C:A,G,T
```

Each alternate allele is processed as an independent genomic variant.

---

## Variant Call Format Representations

VariantFormatter accepts supported genomic VCF representations processed by the shared VariantValidator genomic input conversion pipeline.

For example:

```text
17    50198002    .    C    A
```

Where multiple alternate alleles are supplied, they can be decomposed into individual genomic variants for formatting.

---

## Invalid VariantFormatter Starting Types

VariantFormatter is not a general HGVS-to-HGVS converter.

The following are not accepted as starting inputs:

```text
NM_000088.4:c.589G>T
```

```text
NR_023343.1:n.245G>A
```

```text
NM_000088.4:r.589g>u
```

```text
NP_000079.2:p.Gly197Val
```

These representations may be produced as part of VariantFormatter output, but VariantFormatter processing begins with a genomic variant.

Users who need to validate transcript, RNA, protein or other non-genomic HGVS descriptions should use VariantValidator.

---

# gene2transcripts

Unlike VariantValidator and VariantFormatter, `gene2transcripts` accepts gene and transcript identifiers rather than sequence variant descriptions.

The following input formats are supported.

---

### HGNC gene symbols

Examples:

```text
COL1A1
```

```text
BRCA1
```

```text
TP53
```

---

### HGNC identifiers

Examples:

```text
HGNC:2197
```

```text
HGNC:1100
```

---

### RefSeq transcript identifiers

Examples:

```text
NM_000088.4
```

```text
NR_023343.1
```

---

### Ensembl transcript identifiers

Examples:

```text
ENST00000225964.10
```

```text
ENST00000357654.9
```

---

# hgvs2reference

`hgvs2reference` accepts genomic (`g.`) and coding DNA (`c.`) HGVS sequence variant descriptions and returns the corresponding reference sequence.

---

### Genomic variants (`g.`)

Example:

```text
NC_000017.11:g.50198002C>A
```

Genomic variants are resolved directly against the specified genomic reference sequence.

---

### Coding DNA variants (`c.`)

Example:

```text
NM_000088.4:c.589G>T
```

Coding DNA variants are resolved against the specified transcript reference sequence.

Coding DNA descriptions may also specify intronic positions. For intronic `c.` variants, the genomic reference sequence used for the transcript alignment must be explicitly specified using the HGVS compound reference sequence format.

For example:

```text
NC_000017.11(NM_000088.4):c.589+1G>T
```

The genomic `NC_` reference sequence explicitly identifies the genomic sequence against which the transcript is aligned. This allows `hgvs2reference` to resolve the intronic coordinate and retrieve the corresponding reference sequence.

An intronic coding DNA description containing only the transcript reference sequence, for example:

```text
NM_000088.4:c.589+1G>T
```

does not explicitly identify the genomic reference sequence used for the transcript alignment and therefore cannot be resolved by `hgvs2reference`.

---

### Limitations

The current implementation does not support:

- Non-coding transcript variants (`n.`)
- RNA variants (`r.`)
- Protein variants (`p.`)
- Mitochondrial variants (`m.`)

Intronic `c.` variants are supported when the genomic reference sequence is explicitly specified using the compound `NC_(NM)` HGVS format.

---

# Getting help

VariantValidator has been developed to support a wide range of users, from those new to HGVS nomenclature to experienced clinical scientists and bioinformaticians. If you encounter difficulties preparing a supported input or determining which input format should be used, we encourage you to seek assistance.

Before contacting the development team, you may find the following documentation helpful:

- [Transcript Selection](transcript_selection.md)
- [Output Formats](output_formats.md)
- [Errors and Error Codes](errors_and_error_codes.md)
- [VariantValidator Python API](../python-api/variantvalidator_python.md)
- [VariantFormatter Python API](../python-api/variantformatter_python.md)

If you still require assistance, you can contact the VariantValidator team using our [contact form](https://variantvalidator.org/help/contact/).

Software bugs and feature requests can be reported through the [VariantValidator GitHub issue tracker](https://github.com/openvar/VariantValidator/issues).

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="../../static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="../../static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>
