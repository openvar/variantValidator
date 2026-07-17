# Error Messages and Warning Codes

VariantValidator validates sequence variant descriptions against the Human Genome Variation Society (HGVS) recommendations and a range of supporting reference datasets. During validation, one or more **errors**, **warnings**, or **informational messages** may be returned.

Each diagnostic begins with a standardised **error code** or **warning code**, followed by a descriptive message.

For example:

```text
ReferenceMismatchError: Variant reference (C) does not agree with reference sequence (G)
```

The purpose of these codes is to:

- provide stable identifiers for common validation outcomes;
- make programmatic interpretation of validation results easier;
- provide concise identifiers that can be referenced in documentation and support requests; and
- distinguish between validation failures and advisory information.

This page explains the different categories of messages produced by VariantValidator and provides guidance on interpreting them.

---

# Understanding Validation Messages

Not all validation messages have the same meaning or require the same action.

VariantValidator classifies diagnostics into three broad categories.

## Errors

Errors indicate that the submitted variant description cannot be validated as provided.

In most cases, validation cannot continue until the underlying problem has been corrected by the user.

Typical causes include:

- invalid HGVS syntax;
- incorrect reference sequences;
- incompatible coordinate systems;
- invalid variant coordinates;
- reference sequence mismatches; or
- missing information required to perform validation.

Errors should always be investigated and corrected before reporting or interpreting a variant description.

---

## Corrective Warnings

Corrective warnings indicate that VariantValidator has been able to continue processing the submitted variant, but the submitted description is not fully HGVS compliant or has required automatic correction.

Examples include:

- automatically normalised variant descriptions;
- automatic mapping between equivalent reference sequences;
- use of deprecated descriptions;
- transcript version updates; or
- variant descriptions that are technically valid but not recommended.

Although validation has succeeded, users should normally report the corrected HGVS description rather than the original submitted description.

---

## Informational Warnings

Informational warnings do not indicate that the submitted variant description is incorrect.

Instead, they provide additional information about:

- limitations of the available reference data;
- current support within VariantValidator;
- recommended best practice; or
- additional context that may be useful when interpreting the results.

These messages generally do not require any changes to the submitted variant description.

---

# Frequently Encountered Errors

The following errors are among the most commonly encountered by VariantValidator users.

Although many additional error codes exist, these are the messages most frequently raised during routine validation.

---

## ReferenceMismatchError

**Example**

```text
ReferenceMismatchError:
Variant reference (C) does not agree with reference sequence (G)
```

### What does this mean?

The submitted variant description specifies a reference nucleotide or amino acid that does not match the selected reference sequence.

For example,

```text
NM_000088.4:c.472C>T
```

states that the reference nucleotide at position 472 is **C**, whereas the selected reference sequence actually contains **G**.

### How do I correct it?

Check that:

- the correct reference sequence has been selected;
- the correct transcript version is being used;
- the variant coordinates are correct; and
- the reference allele has been transcribed correctly.

Once corrected, the variant should be resubmitted for validation.

---

## TranscriptCoordinateError

**Example**

```text
TranscriptCoordinateError:
Using a transcript reference sequence to specify a variant position that lies outside of the reference sequence.

Instead re-submit:

NC_000006.12:g.32038297C>T
```

### What does this mean?

HGVS nomenclature does not permit transcript variant descriptions whose coordinates extend beyond the transcript reference sequence.

In these situations VariantValidator determines the equivalent genomic position and recommends validating the genomic variant instead.

### How do I correct it?

Either:

- validate the suggested genomic variant; or
- choose an alternative transcript that spans the affected position.

---

## OutOfBoundsError

**Example**

```text
OutOfBoundsError:
The specified coordinate is outside the boundaries of reference sequence NC_000017.11
```

### What does this mean?

The submitted coordinate lies outside the selected reference sequence.

For example, a chromosome reference sequence may have been selected but the submitted coordinate lies beyond the end of that chromosome.

### How do I correct it?

Verify both:

- the reference sequence accession; and
- the submitted coordinates.

The variant description should then be corrected and resubmitted.

---

## ReferenceSequenceError

**Example**

```text
ReferenceSequenceError:
HGVS variant nomenclature does not allow the use of a gene symbol (COL5A1) in place of a valid reference sequence.

Re-submit using one of:

NM_000093.4
NM_001278074.1
NM_000093.3
```

### What does this mean?

HGVS variant descriptions must use a recognised reference sequence accession.

Gene symbols such as **COL5A1** cannot be used directly as HGVS reference sequences.

VariantValidator will often recommend one or more suitable transcript accessions.

### How do I correct it?

Select one of the recommended transcript reference sequences and resubmit the variant.

---

## UncertainPositionWarning

**Example**

```text
UncertainPositionWarning:
Uncertain positions are not fully supported.
```

### What does this mean?

The submitted variant contains one or more uncertain or undefined coordinates, for example:

```text
c.1576+1_?del
```

Such descriptions are permitted by HGVS nomenclature, but cannot usually be fully validated because the precise genomic location is unknown.

### How do I correct it?

If the exact variant position is known, submit the fully defined HGVS description.

Otherwise this warning simply reflects the current limitations of computational validation for uncertain positions.

> **Note**
>
> Although this message is reported as a warning, in practice it often prevents full validation because the variant coordinates are undefined.

---

# Frequently Encountered Warnings

Unlike errors, warnings do not necessarily indicate that the submitted variant description is invalid.

Many warnings simply inform users that VariantValidator has improved the submitted description, identified a more appropriate reference sequence, or wishes to draw attention to information that may affect interpretation or reporting.

Where possible, users are encouraged to report the corrected HGVS description suggested by VariantValidator.

---

## TranscriptVersionWarning

**Example**

```text
TranscriptVersionWarning:
A more recent version of the selected reference sequence NM_000022.2 is available (NM_000022.3)
```

### What does this mean?

RefSeq and Ensembl reference sequences are periodically updated to correct errors, extend untranslated regions, or improve transcript annotation.

VariantValidator reports when a newer version of the submitted transcript is available.

### Do I need to correct this?

Not necessarily.

HGVS recommendations do not require the latest transcript version to be used. However, using the most recent reference sequence is generally considered best practice and may be required by journals, laboratories or clinical reporting guidelines.

---

## TranscriptDataWarning

**Example**

```text
TranscriptDataWarning:
No transcript definition for NM_003919.3
```

### What does this mean?

VariantValidator validates transcript variants using the VariantValidator Transcript Archive (VVTA).

This warning indicates that the requested transcript version is not currently available within the installed transcript archive.

This most commonly occurs because:

- the transcript is very new;
- the transcript has been withdrawn; or
- the local transcript archive requires updating.

### Do I need to correct this?

Usually not.

If a newer transcript version is available, consider validating against that version instead.

---

## TranscriptRangeWarning

**Example**

```text
TranscriptRangeWarning:
No transcripts found that fully overlap the described variation in the genomic sequence.
```

### What does this mean?

The submitted genomic variant does not completely overlap any known transcript.

This commonly occurs when:

- the variant lies outside all annotated genes;
- the variant overlaps only part of a transcript; or
- the affected region contains no transcript annotations.

### Do I need to correct this?

No.

This warning simply indicates that transcript-level projections cannot be generated.

---

## LrgPendingWarning

**Example**

```text
LrgPendingWarning:
The current status of LRG_XXX is pending therefore changes may be made to the LRG reference sequence.
```

### What does this mean?

Locus Reference Genomic (LRG) records may be released before they have reached their final stable form.

Records designated as **Pending** may still change before becoming public reference sequences.

### Do I need to correct this?

No.

This warning is provided to make users aware that future releases of the LRG record may differ.

---

## IntronSpanningWarning

**Example**

```text
IntronSpanningWarning:
This coding sequence variant description spans at least one intron; use of the corresponding genomic sequence variant descriptions may be invalid.
```

### What does this mean?

Some transcript variants, particularly those identified by RNA sequencing, appear to span one or more introns.

Although such variants can sometimes be projected to genomic coordinates, the genomic representation may not represent a single underlying DNA mutational event.

### Do I need to correct this?

Usually not.

Instead, users should interpret projected genomic variants with appropriate caution.

---

## ProteinSupportWarning

**Example**

```text
ProteinSupportWarning:
Protein level variant descriptions are not fully supported due to redundancy in the genetic code.
```

### What does this mean?

VariantValidator can verify:

- HGVS protein syntax;
- reference amino acids; and
- protein-level nomenclature.

However, because multiple DNA sequence changes can produce the same protein alteration, it is generally impossible to determine a unique genomic or transcript variant from a protein description alone.

### Do I need to correct this?

No.

This warning simply describes a fundamental limitation of protein-level variant interpretation.

---

## RefSeqGeneWarning

**Example**

```text
RefSeqGeneWarning:
RefSeqGene record not available.
```

### What does this mean?

Not every gene has an associated RefSeqGene (NG_) reference sequence.

This warning indicates that no RefSeqGene record currently exists for the gene being validated.

### Do I need to correct this?

No.

Validation continues normally using the available transcript and genomic reference sequences.

---

## VariantNormalizationWarning

**Example**

```text
VariantNormalizationWarning:
Removing redundant reference bases from variant description.
```

### What does this mean?

Some HGVS descriptions include sequence that is already completely defined by the coordinate positions.

For example:

```text
NM_000088.4:c.1_3delATG
```

contains the deleted bases (**ATG**) even though they are already uniquely specified by positions **1–3**.

VariantValidator automatically removes redundant sequence to produce the canonical HGVS representation.

### Do I need to correct this?

No.

This warning simply informs you that the submitted description has been normalised.

---

## VariantMappingWarning

**Example**

```text
VariantMappingWarning:
NM_000088.3:c.589GG>CT automapped to NM_000088.3:c.589_590delGGinsCT
```

### What does this mean?

The submitted variant is not fully HGVS compliant.

VariantValidator has automatically converted it into the canonical HGVS representation.

### Do I need to correct this?

Yes.

The corrected HGVS description returned by VariantValidator should normally be used for reporting.

---

## LrgMappingWarning

**Example**

```text
LrgMappingWarning:
LRG_1:g.8638G>T automapped to NG_007400.1:g.8638G>T
```

or

```text
LRG_1t1:c.589G>T automapped to NM_000088.3:c.589G>T
```

### What does this mean?

LRG reference sequences have direct equivalents within the RefSeq reference sequence collection.

VariantValidator automatically reports the corresponding RefSeq representation to improve interoperability with other databases and software.

### Do I need to correct this?

No.

The original LRG description remains valid.

The mapped RefSeq representation is provided as additional information and is generally recommended for interoperability with other resources.

---

# Complete Error Code Reference

The following table summarises every error code currently recognised by VariantValidator.

Errors indicate that the submitted variant description could not be validated as provided. Unless otherwise stated, the submitted variant should be corrected and resubmitted.

| Error Code | Description |
|------------|-------------|
| AccessionVersionError | A required reference sequence version has not been supplied. |
| AlleleMergeError | Multiple variants should be merged into a single HGVS allele description. |
| AlleleSyntaxError | The submitted allele description is not valid HGVS syntax. |
| AlleleVariantError | The submitted allele description is biologically or structurally invalid. |
| AminoAcidError | An invalid amino acid code has been supplied. |
| AminoMismatchError | The supplied reference amino acid does not match the selected reference sequence. |
| BaseOffsetError | An invalid intronic or offset coordinate has been supplied. |
| CDSBoundaryError | The supplied coding coordinate lies outside the coding sequence. |
| CharacterCaseError | Uppercase and lowercase characters have been used incorrectly. |
| CircularReferenceError | The submitted variant spans the origin of a circular reference sequence and cannot be represented as submitted. |
| CodingTranscriptError | A coding transcript has been used where a non-coding transcript is required. |
| ConversionNormalizationError | The requested coordinate conversion could not be normalised. |
| ExonBoundaryError | The submitted variant defines an invalid exon boundary. |
| ExpandedRepeatError | The submitted expanded repeat description is invalid. |
| FuzzyPositionError | The submitted position is undefined or uncertain. |
| FuzzyRangeError | The submitted coordinate interval is undefined or uncertain. |
| GenomeBuildError | The selected genome build is invalid or unsupported. |
| GenomeMismatchError | The submitted reference sequence is incompatible with the selected genome build. |
| GenomeReferenceError | The submitted genomic reference sequence is invalid for the selected assembly. |
| HgvsParserError | The submitted HGVS description could not be parsed. |
| HgvsSyntaxError | The submitted HGVS description does not follow HGVS syntax. |
| IncompatibleTypeError | The submitted coordinate type is incompatible with the reference sequence. |
| IncompleteVariantError | The submitted variant description is incomplete. |
| InsertionCoordinateError | Required insertion coordinates are missing. |
| InsertionLengthError | The insertion length is invalid. |
| InsertionPositionError | An insertion must be described between two adjacent positions. |
| InsertionSequenceError | No inserted sequence has been supplied. |
| InternalValidationError | An unexpected internal validation error occurred. |
| IntervalOrderError | The variant interval is invalid because the coordinates are in the wrong order. |
| IntronicAlleleError | Allelic descriptions cannot contain intronic coordinates. |
| IntronicVariantError | The submitted intronic variant cannot be validated in the current context. |
| InvalidRangeError | The submitted coordinate range is invalid. |
| InvalidReferenceError | The submitted reference sequence identifier is invalid. |
| InvalidVariantError | The submitted variant description is invalid. |
| LegacySyntaxError | The submitted variant uses obsolete HGVS nomenclature. |
| MitochondrialBuildError | The mitochondrial reference sequence is incompatible with the selected genome build. |
| MitochondrialReferenceError | A mitochondrial reference sequence has been described using the wrong coordinate system. |
| MultipleAlignmentError | Multiple incompatible alignments were identified. |
| NonCodingTranscriptError | A non-coding transcript has been used where a coding transcript is required. |
| NumericInputError | Numeric input alone is not a valid variant description. |
| OutOfBoundsError | The supplied coordinate lies outside the selected reference sequence. |
| OverlappingPositionError | The submitted variant positions overlap or are otherwise invalid. |
| ProteinNormalizationError | Protein-level variants cannot be normalised. |
| ProteinReferenceError | A protein reference sequence has been used for a nucleotide variant. |
| ProteinRequiredError | A protein reference sequence is required for this operation. |
| ReferenceMismatchError | The supplied reference nucleotide or amino acid does not match the selected reference sequence. |
| ReferenceSequenceError | The submitted reference sequence is invalid or unavailable. |
| ReferenceTypeError | The submitted reference sequence type is incompatible with the variant description. |
| RepeatSyntaxError | The submitted expanded repeat syntax is invalid. |
| RnaAlphabetError | RNA variants must use the RNA alphabet (U instead of T). |
| SequenceGenerationError | The variant sequence could not be generated. |
| SequenceLookupError | The required reference sequence could not be retrieved. |
| TranscriptCoordinateError | The submitted transcript coordinate lies outside the transcript reference sequence. |
| TranscriptDataError | Required transcript annotation is unavailable. |
| TranscriptIdentificationError | A suitable transcript could not be identified. |
| TranscriptMappingError | The submitted variant could not be mapped between reference sequences. |
| TranscriptMissingError | The requested transcript could not be found. |
| TranscriptReferenceError | A transcript reference sequence is required. |
| TranscriptSelectionError | The requested transcript selection is invalid. |
| TranscriptTypeError | The selected transcript type is incompatible with the submitted variant. |
| UncertainConversionError | Validation could not continue because an uncertain position could not be converted. |
| UncertainSequenceError | The submitted sequence contains unresolved uncertainty. |
| UndefinedSequenceError | Required reference sequence information is unavailable. |
| UnsupportedEditError | The requested edit operation is not currently supported. |
| ValidatorSubmissionError | The submitted validation request is incomplete or invalid. |
| VariantSyntaxError | The submitted variant contains an HGVS syntax error. |
| VariantTypeError | The submitted variant type is incompatible with the requested operation. |
| VcfFormatError | The submitted VCF description is incomplete or malformed. |

> **Note**
>
> Some errors include additional explanatory text, suggested transcript accessions, or recommended corrected variant descriptions. Where available, these suggestions should be followed when correcting and resubmitting the variant.

# Complete Warning Code Reference

The following table summarises all warning and informational codes currently recognised by VariantValidator.

Unlike errors, warnings do not necessarily indicate that the submitted variant description is invalid. Some warnings identify HGVS recommendations or best practice, whereas others simply provide additional information about the validation process.

Users should always review any warnings that accompany a validation result.

| Warning Code | Description |
|--------------|-------------|
| AlignmentDataWarning | Alignment data are incomplete or unavailable for the requested operation. |
| AlignmentGapWarning | Differences exist between transcript and genomic reference sequences. |
| AlleleExtractionWarning | An allele description has been extracted automatically from the submitted variant. |
| AlleleValidationWarning | Individual alleles have been validated independently. |
| AlternativeRepresentationWarning | An alternative HGVS representation is recommended. |
| AmbiguousVcfWarning | The submitted VCF representation is ambiguous. |
| AutoCorrectionWarning | VariantValidator has automatically corrected the submitted variant description. |
| FrameRestoreWarning | Combining neighbouring variants restores the reading frame. |
| GenomeMappingWarning | The submitted variant has been mapped to an equivalent genomic representation. |
| GenomeReferenceWarning | The selected genomic reference sequence should be reviewed. |
| GenomeMismatchWarning | The submitted reference sequence may not correspond to the selected genome build. |
| GenomeSupportWarning | The selected genome assembly has limited support. |
| InitiationCodonWarning | The variant affects the translation initiation codon. |
| IntronBaseMismatchWarning | Reference bases differ after transcript-to-genome mapping. |
| IntronSpanningWarning | The submitted variant spans one or more introns. |
| IntronicValidationWarning | Validation of intronic variants is limited. |
| LrgMappingWarning | The submitted LRG reference has been mapped to the equivalent RefSeq reference sequence. |
| LrgPendingWarning | The submitted LRG record is still pending publication and may change. |
| MultipleAlleleWarning | Multiple alternate alleles have been detected. |
| ParRegionWarning | The variant lies within or near the pseudoautosomal region. |
| ProteinLimitationWarning | Protein-level validation is inherently limited by the redundancy of the genetic code. |
| ProteinSupportWarning | Protein-level validation is only partially supported. |
| ProteinTerminationWarning | The variant affects the protein termination codon. |
| ProteinTranslationWarning | Protein translation could not be completed exactly as expected. |
| PseudoautosomalRegionWarning | The submitted variant lies within a pseudoautosomal region. |
| ReferenceUpdateWarning | A newer equivalent reference sequence is available. |
| RefSeqGeneWarning | No RefSeqGene record is available for the affected gene. |
| SingleTranscriptWarning | Validation has been restricted to a single transcript. |
| TranscriptBoundsWarning | The supplied transcript coordinates extend beyond the available transcript annotation. |
| TranscriptRangeWarning | No transcript completely spans the submitted variant. |
| TranscriptSelectionWarning | Transcript selection has been modified automatically. |
| TranscriptSuggestionWarning | One or more alternative transcript reference sequences are available. |
| TranscriptVersionWarning | A newer version of the selected transcript is available. |
| UncertainPositionWarning | Validation of uncertain variant positions is limited. |
| UncertainVariantWarning | The submitted variant contains unresolved uncertainty. |
| VariantConversionWarning | The submitted variant has been converted into HGVS-compliant form. |
| VariantMappingWarning | The submitted variant has been automatically mapped to an equivalent HGVS representation. |
| VariantNormalizationWarning | The submitted variant has been normalised to the canonical HGVS representation. |
| VcfConversionWarning | The submitted VCF description has been converted into HGVS nomenclature. |

## Informational Codes

The following informational codes are provided for completeness. These messages do not indicate errors or warnings, but provide additional context regarding validation or external libraries.

| Information Code | Description |
|------------------|-------------|
| LovdSyntaxcheckInvalid | Invalid syntax reported by the integrated LOVD syntax checker. |
| LovdSyntaxcheckLibraryVersion | Version of the LOVD syntax checking library used during validation. |
| LovdSyntaxcheckSource | Source of the LOVD syntax checking results. |
| LovdSyntaxcheckValid | Confirmation that the submitted syntax is valid according to the LOVD syntax checker. |
| LovdSyntaxcheckWarning | Additional information reported by the LOVD syntax checker. |
| ProteinTranslationInfo | Additional information relating to predicted protein translation. |

---

# Notes for API Users

Every validation message begins with a stable machine-readable code followed by a human-readable description.

For example:

```text
ReferenceMismatchError: Variant reference (C) does not agree with reference sequence (G)
```

Applications should use the code (`ReferenceMismatchError`) when filtering, categorising or processing validation results programmatically, rather than relying on the wording of the descriptive text.

The descriptive text may be expanded over time to improve clarity, while the corresponding error or warning code will remain stable wherever possible.

Multiple errors and warnings may be returned for a single submitted variant. Users are encouraged to review all reported diagnostics before interpreting or reporting a variant description.