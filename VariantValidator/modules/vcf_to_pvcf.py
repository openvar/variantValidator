class VcfConversionError(Exception):
    """Custom exception raised when a VCF line cannot be converted to shorthand."""
    pass


def split_vcf_line(vcf_line):
    """
    Split a VCF line by detecting its delimiter.

    Raises VcfConversionError if no supported delimiter is found.
    """
    line = vcf_line.strip()

    if "\t" in line:
        return line.split("\t")

    fields = line.split()
    if len(fields) > 1:
        return fields

    raise VcfConversionError(
        "Unable to detect delimiter. Expected tab ('\\t') or "
        "whitespace-separated values."
    )


def vcf_to_shorthand(vcf_line):
    """
    Convert a single VCF line into chr-start-end-TYPE shorthand.

    Handles:
      - CNVs with <DEL>/<DUP>/<INV>
      - Simple SNVs/indels (chr-pos-ref-alt)
      - Optional CN field

    Raises:
      VcfConversionError with descriptive message if conversion fails.
    """

    # Skip headers
    if vcf_line.startswith("#"):
        raise VcfConversionError(
            "Header line cannot be converted. "
            "Please provide a variant record line."
        )

    fields = split_vcf_line(vcf_line)

    if len(fields) < 5:
        raise VcfConversionError(
            f"VCF line has insufficient columns "
            f"(found {len(fields)}, expected ≥5). "
            "Ensure the line contains at least CHROM, POS, ID, REF, ALT."
        )

    chrom, pos, _id, ref, alt = fields[:5]
    info = fields[7] if len(fields) > 7 else ""

    # Validate position
    try:
        pos = int(pos)
    except ValueError:
        raise VcfConversionError(
            f"Invalid POS field: '{pos}' is not an integer."
        )

    # Structural variants
    if alt in ("<DEL>", "<DUP>", "<INV>", "DEL", "DUP", "INV"):
        end = None
        cn = None

        for entry in info.split(";"):
            if entry.startswith("END="):
                try:
                    end = int(entry.split("=", 1)[1])
                except ValueError:
                    raise VcfConversionError(
                        f"Invalid END value in INFO field: '{entry}'."
                    )

            elif entry.startswith("SVLEN=") and end is None:
                try:
                    end = pos + abs(
                        int(entry.split("=", 1)[1])
                    )
                except ValueError:
                    raise VcfConversionError(
                        f"Invalid SVLEN value in INFO field: '{entry}'."
                    )

            elif entry.startswith("CN="):
                cn = entry.split("=", 1)[1]

        if end is None:
            raise VcfConversionError(
                "Cannot determine end position for structural variant. "
                "INFO field must contain END= or SVLEN=."
            )

        alt_clean = (
            alt[1:-1]
            if alt.startswith("<") and alt.endswith(">")
            else alt
        )

        shorthand = f"{chrom}-{pos}-{end}-{alt_clean}"

        if cn:
            shorthand += f"[CN{cn}]"

        return shorthand

    # Simple SNV/indel
    if not ref or not alt:
        raise VcfConversionError(
            "Missing REF or ALT allele. Cannot convert variant."
        )

    return f"{chrom}-{pos}-{ref}-{alt}"


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
