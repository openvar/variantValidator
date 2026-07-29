# -*- coding: utf-8 -*-

"""
Provide the functional and object-oriented VariantFormatter APIs.
"""

import collections
import json

import VariantFormatter
import VariantFormatter.variantformatter as vf
import VariantValidator
from VariantValidator.modules import vcf_to_pvcf

import logging

logger = logging.getLogger(
    f"VariantValidator.VariantFormatter.{__name__.removeprefix('VariantFormatter.')}"
)


class FormatterSubmissionError(Exception):
    pass


_GLOBAL_VFO = None
_METADATA = None


def _get_global_validator():
    """Return the lazily created validator used by the legacy API."""
    global _GLOBAL_VFO

    if _GLOBAL_VFO is None:
        _GLOBAL_VFO = VariantValidator.Validator()

    return _GLOBAL_VFO


def _get_metadata():
    """Return cached VariantValidator and VariantFormatter metadata."""
    global _METADATA

    if _METADATA is None:
        metadata = _get_global_validator().my_config()
        metadata["variantformatter_version"] = VariantFormatter.__version__

        sr_root, sr_version = metadata["vvseqrepo_db"].split("/")[-2:]
        metadata["vvseqrepo_db"] = f"{sr_root}/{sr_version}"

        _METADATA = metadata

    return _METADATA


def _normalise_liftover_level(liftover_level):
    """Normalise and validate the requested liftover level."""
    if liftover_level in ("True", True, 1):
        return True

    if liftover_level in (False, "False", 0, None):
        return False

    if liftover_level == "primary":
        return "primary"

    raise FormatterSubmissionError(
        f"liftover_level '{liftover_level}' is not supported. "
        "Use True, False, None or 'primary'."
    )


def _validate_submission(variant, genome, liftover_level):
    """Validate formatter input and return the normalised liftover level."""
    if genome is None:
        raise FormatterSubmissionError("Genome build is required")

    if variant is None:
        raise FormatterSubmissionError("Variant is required")

    return _normalise_liftover_level(liftover_level)


def _normalise_transcript_selection(specify_transcripts):
    """Normalise JSON-style transcript-selection arguments."""
    mapping = {
        '["all"]': "all",
        '["raw"]': "raw",
        '["mane"]': "mane",
        '["mane_select"]': "mane_select",
        '["select"]': "select",
    }
    return mapping.get(specify_transcripts, specify_transcripts)


def _normalise_batch_input(batch_input):
    """Return formatter input as a list of variants."""
    if isinstance(batch_input, list):
        return batch_input

    try:
        parsed_input = json.loads(batch_input)
    except (json.decoder.JSONDecodeError, TypeError):
        return [batch_input]

    if isinstance(parsed_input, list):
        return parsed_input

    return [parsed_input]


def _contains_hgvs_type(variant):
    """Return whether the string contains an HGVS sequence type marker."""
    return any(
        marker in variant
        for marker in ("g.", "c.", "r.", "n.", "m.", "o.")
    )


def _process_vcf_input(variant):
    """Convert a tab-delimited VCF record to VariantFormatter shorthand."""
    warnings = []

    if "\t" not in variant or _contains_hgvs_type(variant):
        return variant, warnings

    try:
        converted_variant = vcf_to_pvcf.vcf_to_shorthand(variant)
    except vcf_to_pvcf.VcfConversionError:
        return variant, warnings

    variant = converted_variant
    warnings.append(
        f"VcfConversionWarning: VCF line identified and converted to {variant}"
    )

    vcf_data = variant.replace(":", "-").split("-")

    if len(vcf_data) < 4:
        return variant, warnings

    position = vcf_data[2]
    edit = vcf_data[3]
    edit_lower = edit.lower()

    if (
        position.isdigit()
        and ("del" in edit_lower or "inv" in edit_lower)
        and not any(
            _contains_hgvs_type(item)
            for item in vcf_data
        )
    ):
        variant = (
            f"{vcf_data[0]}:{vcf_data[1]}_"
            f"{position}{edit_lower}"
        )
        warnings.append(
            f"VcfConversionWarning: CNV identified, and mapped to {variant}"
        )

    return variant, warnings


def _looks_like_pseudo_vcf(variant):
    """Return whether input appears to use pseudo-VCF shorthand."""
    if variant.startswith("LRG"):
        return False

    separator_positions = [
        position
        for position in (variant.find(":"), variant.find("-"))
        if position > 0
    ]

    if not separator_positions:
        return False

    prefix = variant[:min(separator_positions)]
    return prefix.replace("_", "").isalnum()


def _recover_incomplete_pseudo_vcf(
        variant,
        pseudo_vcf,
        genome_build,
        validator,
        output,
):
    """Attempt to recover incomplete pseudo-VCF input through validation."""
    try:
        result = validator.validate(
            variant,
            genome_build,
            "check_only",
        ).format_as_dict(test=True)

        hgvs = result["intergenic_variant_1"][
            "primary_assembly_loci"
        ][genome_build.lower()]["hgvs_genomic_description"]

        if "NC_" not in hgvs:
            raise KeyError

    except Exception:
        output["errors"].append(
            f"{pseudo_vcf} is an unsupported format: "
            "For assistance, submit variant description "
            "to https://rest.variantvalidator.org"
        )
        output["flag"] = "submission_warning"
        return None

    output["errors"].append(
        f"{pseudo_vcf} is not HGVS compliant because a valid "
        f"reference sequence has not been provided. "
        f"Updating to {hgvs}"
    )

    return hgvs


def _expand_pseudo_vcf(variant, genome_build, validator, output):
    """Return the variants represented by pseudo-VCF input."""
    if not _looks_like_pseudo_vcf(variant):
        return [variant]

    delimiter = ":" if ":" in variant else "-"
    vcf_list = variant.split(delimiter)

    if len(vcf_list) != 4:
        recovered_variant = _recover_incomplete_pseudo_vcf(
            variant,
            variant,
            genome_build,
            validator,
            output,
        )

        if recovered_variant is None:
            return []

        return [recovered_variant]

    alternate = vcf_list[-1]

    if "," not in alternate:
        return [variant]

    prefix = vcf_list[:3]
    return [
        delimiter.join(prefix + [alt])
        for alt in alternate.split(",")
    ]


def _format_variant(
        needs_formatting,
        genome_build,
        validator,
        transcript_model,
        specify_transcripts,
        check_only,
        liftover,
        legacy_genomic_structure=True,
):
    """Format one variant and return its formatter object and data."""
    result = vf.FormatVariant(
        needs_formatting,
        genome_build,
        validator,
        transcript_model,
        specify_transcripts,
        check_only,
        liftover,
        legacy_genomic_structure=legacy_genomic_structure,
    )

    structured_data = result.stucture_data()

    return result, structured_data[needs_formatting]


def _format_impl(
        batch_input,
        genome_build,
        transcript_model=None,
        specify_transcripts=None,
        check_only=False,
        liftover=False,
        validator=None,
        testing=None,
        legacy_genomic_structure=True,
):
    """Shared implementation for the functional and object-oriented APIs."""
    validator.testing = bool(testing)

    specify_transcripts = _normalise_transcript_selection(
        specify_transcripts
    )
    validator.select_transcripts = specify_transcripts

    formatter_transcripts = (
        None
        if specify_transcripts == "all"
        else specify_transcripts
    )

    formatted_variants = collections.OrderedDict()

    for submitted_variant in _normalise_batch_input(batch_input):
        variant = submitted_variant.strip()

        variant, vcf_processing_warnings = _process_vcf_input(variant)
        variant = "".join(variant.split())

        output = collections.OrderedDict(
            errors=[],
            flag=None,
        )
        formatted_variants[variant] = output

        format_these = _expand_pseudo_vcf(
            variant,
            genome_build,
            validator,
            output,
        )

        for needs_formatting in format_these:
            result, formatted_data = _format_variant(
                needs_formatting,
                genome_build,
                validator,
                transcript_model,
                formatter_transcripts,
                check_only,
                liftover,
                legacy_genomic_structure=legacy_genomic_structure,
            )

            output["flag"] = result.warning_level
            output[needs_formatting] = formatted_data

            if vcf_processing_warnings:
                formatted_data["genomic_variant_warnings"] = (
                    vcf_processing_warnings
                )

    formatted_variants["metadata"] = _get_metadata()

    return formatted_variants


def format(
        variant=None,
        genome=None,
        transcript_model=None,
        select_transcripts=None,
        checkOnly=False,
        liftover_level=True,
        validator=None,
        testing=None,
        legacy_genomic_structure=True,
):
    """
    Format one or more variants using the legacy functional API.

    Parameters retain their historical names for backwards compatibility.
    """
    liftover_level = _validate_submission(
        variant,
        genome,
        liftover_level,
    )

    if validator is None:
        validator = _get_global_validator()

    return _format_impl(
        batch_input=variant,
        genome_build=genome,
        transcript_model=transcript_model,
        specify_transcripts=select_transcripts,
        check_only=checkOnly,
        liftover=liftover_level,
        validator=validator,
        testing=testing,
        legacy_genomic_structure=legacy_genomic_structure,
    )


class SimpleVariantFormatter:
    """
    Object-oriented VariantFormatter API.

    Each instance owns a VariantValidator instance and can therefore be
    independently pooled or reused.
    """

    def __init__(self, *, testing=False, legacy_genomic_structure=True):
        self.validator = VariantValidator.Validator()
        self.testing = testing
        self.legacy_genomic_structure = legacy_genomic_structure

    def format(
            self,
            variant=None,
            genome=None,
            transcript_model=None,
            select_transcripts=None,
            checkOnly=False,
            liftover_level=True,
            legacy_genomic_structure=None,
    ):
        """Format one or more variants."""
        liftover_level = _validate_submission(
            variant,
            genome,
            liftover_level,
        )

        if legacy_genomic_structure is None:
            legacy_genomic_structure = self.legacy_genomic_structure

        return _format_impl(
            batch_input=variant,
            genome_build=genome,
            transcript_model=transcript_model,
            specify_transcripts=select_transcripts,
            check_only=checkOnly,
            liftover=liftover_level,
            validator=self.validator,
            testing=self.testing,
            legacy_genomic_structure=legacy_genomic_structure,
        )


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
