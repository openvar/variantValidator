# Supported Input Formats

This guide describes the input formats accepted throughout the VariantValidator software suite.

The supported input formats are common across the following tools:

- [VariantValidator](#variantvalidator)
- [VariantFormatter](#variantformatter)
- [gene2transcripts](#gene2transcripts)
- [hgvs2reference](#hgvs2reference)

These tools share the same parsing and normalisation engine. Consequently, many of the input formats described in this guide are recognised consistently throughout the software suite unless explicitly stated otherwise.

Where a particular tool supports additional input formats or imposes specific restrictions, these are described in the relevant section of this guide and in that tool's individual documentation.

VariantValidator accepts a wide range of sequence variant descriptions. In addition to fully compliant HGVS sequence variant descriptions, it recognises many commonly encountered non-HGVS representations, legacy formats and common formatting mistakes, automatically converting them into valid HGVS sequence variant descriptions where possible.

The parser automatically detects the submitted format, performs any necessary preprocessing, and then validates or formats the resulting variant. Where an unambiguous correction cannot be made, informative warnings or error messages are returned to assist the user.

!!! tip "Related documentation"

    - See the [Transcript Selection](transcript_selection.md) guide for information on transcript selection strategies such as `mane_select`, `mane`, `select` and `all`.
    - See the [Output Formats](output_formats.md) guide to understand the results returned by each tool.
    - See the [Errors and Error Codes](errors_and_error_codes.md) reference for explanations of validation errors and warnings.

The following sections describe the supported input formats accepted by each tool within the VariantValidator software suite.

---

## Scope of supported input formats

The examples described in this guide are representative rather than exhaustive.

VariantValidator has been designed to recognise and interpret a wide range of real-world sequence variant descriptions. In addition to fully compliant HGVS sequence variant descriptions, it accepts many commonly encountered non-HGVS formats, legacy representations and common formatting mistakes originating from clinical laboratories, research pipelines, databases and published literature.

Together, VariantValidator and the integrated LOVD HGVS Syntax Checker recognise and process a broader range of HGVS syntax, legacy representations and common user errors than any other currently available software. Many non-standard inputs are recognised automatically and, where possible, converted into valid HGVS sequence variant descriptions before validation. When automatic correction is not possible, VariantValidator provides informative error messages and guidance to help users generate a valid HGVS sequence variant description.

Support for additional input formats continues to evolve as new real-world examples are encountered. If you regularly encounter a sequence variant description that is not currently recognised, we encourage you to submit a feature request through the project's contact page so that support can be considered for a future release.

---

# VariantValidator

---

## HGVS Sequence Variant Types

VariantValidator supports the HGVS sequence variant types defined by the Human Genome Variation Society (HGVS).

### Genomic variants (`g.`)

Genomic variants describe sequence changes relative to a genomic reference sequence.

Example:

```text
NC_000017.11:g.50198002C>A
```

This is the most commonly submitted variant type and forms the primary input for VariantValidator.

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

RNA variants must use the IUPAC RNA alphabet (for example, **U** rather than **T**). Where possible, informative error messages are returned for incorrect RNA syntax.

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

Uncertain intervals, fuzzy breakpoints, and unknown variant boundaries are recognised where supported by the HGVS recommendations. If a submitted description cannot be interpreted unambiguously, VariantValidator returns an informative validation message describing the problem.

---

### LRG reference sequences

Legacy LRG reference sequences are recognised throughout the software suite.

Examples include:

```text
LRG_199:c.589G>T
```

```text
LRG_199t1:c.589G>T
```

Where appropriate, LRG identifiers are automatically converted to their equivalent RefSeq reference sequences before further processing.

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

Compound reference sequence descriptions define the genomic sequence against which the transcript is aligned. This allows transcript variants to be interpreted relative to an alternative genomic reference sequence, such as a RefSeqGene record, chromosome, scaffold or other supported genomic reference sequence.

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

In addition to fully compliant HGVS sequence variant descriptions, VariantValidator accepts several commonly encountered non-HGVS variant representations.

These formats commonly originate from variant calling pipelines, databases, spreadsheets and legacy software. Where possible, VariantValidator automatically converts these representations into valid HGVS sequence variant descriptions before validation.

The following non-HGVS input formats are recognised.


### Gene symbols used as reference sequence identifiers

Gene symbols are frequently—but incorrectly—used in place of transcript reference sequence identifiers.

For example,

```text
COL1A1:c.589G>T
```

This is **not valid HGVS nomenclature**, as HGVS requires a transcript reference sequence (for example, `NM_000088.4`) rather than a gene symbol. Because many genes have multiple transcript reference sequences, a gene symbol alone is insufficient to uniquely identify an HGVS sequence variant description.

If a transcript selection strategy is supplied (for example `mane_select`, `mane`, `select`, or a user-specified transcript accession), VariantValidator automatically substitutes the appropriate transcript reference sequence and continues validation.

If no transcript selection strategy is provided, VariantValidator cannot determine which transcript was intended. Validation stops and a warning is returned together with a list of compatible transcript reference sequences that may be used to resubmit the variant.

Providing a transcript selection strategy or explicitly specifying the desired transcript is therefore recommended whenever gene symbols are used as input.

---

### Pseudo-VCF chromosome notation

VariantValidator accepts several simplified chromosome coordinate formats that are commonly used in spreadsheets, databases and analysis pipelines.

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

VariantValidator also accepts pseudo-VCF descriptions containing multiple alternate alleles.

For example,

```text
17-50198002-C-A,G,T
```

```text
17:50198002:C:A,G,T
```

HGVS requires each alternate allele to be represented as an independent sequence variant description. VariantValidator therefore decomposes pseudo-VCF descriptions containing multiple alternate alleles into individual variants before converting each into HGVS format and validating them independently.

---

### Genome assembly prefixes

Some pipelines include the genome assembly as part of the chromosome coordinate description.

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

VariantValidator recognises these formats, extracts the genome assembly, converts the description into HGVS format and validates the resulting variant. If the embedded genome assembly conflicts with the selected genome assembly, validation fails with an informative error.

---

### VCF/HGVS hybrid formats

VariantValidator accepts several hybrid formats that combine HGVS reference sequence identifiers with VCF-style coordinate or allele notation.

These representations commonly arise following partial or incomplete conversion of VCF records into HGVS sequence variant descriptions.

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

VariantValidator recognises these hybrid representations, extracts any embedded genome assembly information where present, resolves chromosome identifiers to the appropriate reference sequence accession, converts the description into valid HGVS syntax, and then performs standard validation.

If the embedded genome assembly conflicts with the selected genome assembly, validation fails with an informative error.

---

### Chromosome identifiers

Chromosome identifiers are accepted in several commonly used forms.

Examples include:

```text
17
```

```text
chr17
```

```text
Chr17
```

These identifiers are automatically mapped to the appropriate genomic reference sequence accession for the selected genome assembly.

---

### Variant Call Format (VCF)

VariantValidator accepts complete Variant Call Format (VCF) records copied directly from VCF files.

For example,

```text
17    50198002    .    C    A
```

The chromosome, position, reference allele and alternate allele fields are used to construct the corresponding HGVS sequence variant description. Additional VCF columns, including quality scores, filter status and genotype information, are ignored.

Both single-record VCF input and multi-allelic VCF records are recognised. Where appropriate, multiple alternate alleles are decomposed into individual HGVS sequence variant descriptions and validated independently.

---

The following sections describe any additional input conventions or restrictions that apply to the remaining tools within the VariantValidator software suite.

---

## Additional supported formats

The examples presented in this document illustrate the most commonly encountered HGVS and non-HGVS input formats accepted by VariantValidator. They are **not** intended to be an exhaustive list of every supported syntax.

VariantValidator has been developed to recognise the broad range of sequence variant descriptions encountered in real-world clinical and research workflows. In addition to supporting the current HGVS recommendations, it recognises many legacy representations, common formatting mistakes and hybrid notations that frequently occur in publications, databases, spreadsheets and bioinformatics pipelines.

Where possible, VariantValidator automatically converts these representations into valid HGVS sequence variant descriptions before validation. When an unambiguous correction cannot be made, the software returns informative warnings or error messages describing the problem and, where appropriate, guidance on how the submitted description can be corrected.

Together with the integrated LOVD HGVS Syntax Checker, VariantValidator supports a broader range of HGVS syntax and commonly encountered input formats than any other currently available validation tool. New input formats and common user mistakes are continually incorporated as they are encountered to improve compatibility and user experience.

If you encounter a sequence variant description that is not recognised, or would like support for an additional input format, please submit a feature request through the project's contact page.

---

# VariantFormatter

VariantFormatter is designed for automated bioinformatics pipelines and software integration.

Unlike VariantValidator, which is intended to recognise and validate a very broad range of HGVS and non-HGVS sequence variant descriptions, VariantFormatter assumes that submitted variants are already of good quality. It performs minimal input recovery before converting accepted variant descriptions into equivalent genomic, transcript and protein representations.

This design makes VariantFormatter well suited to production workflows where variant descriptions have already been validated or originate from trusted sources.

VariantFormatter supports the following input formats:

- Genomic HGVS (`g.`) sequence variant descriptions.
- Pseudo-VCF chromosome coordinate formats.
- Pseudo-VCF descriptions containing multiple alternate alleles.
- Complete Variant Call Format (VCF) records.

Unlike VariantValidator, VariantFormatter does not attempt to recognise the extensive range of legacy HGVS syntax, hybrid formats or common user formatting errors accepted by the VariantValidator parser. Users wishing to validate or recover imperfect sequence variant descriptions should use VariantValidator before formatting.

---

## HGVS genomic sequence variants

Genomic HGVS (`g.`) sequence variant descriptions are the primary input accepted by VariantFormatter.

Example:

```text
NC_000017.11:g.50198002C>A
```

These descriptions are formatted into equivalent transcript and protein representations where possible.

---

## Pseudo-VCF chromosome notation

VariantFormatter accepts the simplified chromosome coordinate formats commonly produced by variant calling pipelines.

Examples include:

```text
17-50198002-C-A
```

```text
17:50198002:C:A
```

```text
GRCh38-17-50198002-C-A
```

```text
GRCh38:17:50198002:C:A
```

These descriptions are converted into genomic HGVS sequence variant descriptions before formatting.

---

## Multiple alternate alleles

Pseudo-VCF descriptions containing multiple alternate alleles are supported.

Examples include:

```text
17-50198002-C-A,G,T
```

```text
17:50198002:C:A,G,T
```

Each alternate allele is decomposed into an independent HGVS sequence variant description before formatting.

---

## Variant Call Format (VCF)

VariantFormatter accepts complete Variant Call Format (VCF) records.

For example,

```text
17    50198002    .    C    A
```

The chromosome, position, reference allele and alternate allele fields are extracted and converted into genomic HGVS before formatting. Additional VCF columns, such as quality metrics, filters and genotype information, are ignored.

Multi-allelic VCF records are also supported and are decomposed into individual HGVS sequence variant descriptions before formatting.

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

`hgvs2reference` accepts HGVS sequence variant descriptions and returns the corresponding reference sequence.

The following input formats are supported.

---

### Genomic variants (`g.`)

Example:

```text
NC_000017.11:g.50198002C>A
```

---

### Coding DNA variants (`c.`)

Example:

```text
NM_000088.4:c.589G>T
```

Coding variants are automatically converted to their corresponding non-coding transcript coordinates before the reference sequence is retrieved.

---

### Non-coding transcript variants (`n.`)

Example:

```text
NR_023343.1:n.245G>A
```

---

### Limitations

The current implementation does not support:

- RNA variants (`r.`)
- Protein variants (`p.`)
- Mitochondrial variants (`m.`)
- Compound genomic/transcript reference sequence descriptions (for example `NG_(NM_):c.` or `NC_(NM_):c.`)
- Fully intronic transcript variants (a warning is returned requesting the use of a genomic reference sequence instead)
