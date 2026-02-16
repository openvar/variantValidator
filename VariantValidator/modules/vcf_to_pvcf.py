import re

class VcfConversionError(Exception):
    """Custom exception raised when a VCF line cannot be converted to shorthand."""
    pass


def split_vcf_line(vcf_line):
    """
    Split a VCF line by detecting delimiter (tab).
    Raises VcfConversionError if no supported delimiter is found.
    """
    line = vcf_line.strip()

    if "\t" in line:
        return line.split("\t")

    elif re.search("\s+", line):
        return re.split(r"\s+", line.strip())

    raise VcfConversionError(
        "Unable to detect delimiter. Expected tab ('\\t') values."
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
            "Header line cannot be converted. Please provide a variant record line."
        )

    fields = split_vcf_line(vcf_line)

    if len(fields) < 5:
        raise VcfConversionError(
            f"VCF line has insufficient columns (found {len(fields)}, expected ≥5). "
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
                    end = int(entry.split("=")[1])
                except ValueError:
                    raise VcfConversionError(
                        f"Invalid END value in INFO field: '{entry}'."
                    )

            elif entry.startswith("SVLEN=") and end is None:
                try:
                    end = pos + abs(int(entry.split("=")[1]))
                except ValueError:
                    raise VcfConversionError(
                        f"Invalid SVLEN value in INFO field: '{entry}'."
                    )

            elif entry.startswith("CN="):
                cn = entry.split("=")[1]

        if end is None:
            raise VcfConversionError(
                "Cannot determine end position for structural variant. "
                "INFO field must contain END= or SVLEN=."
            )

        # Remove angle brackets
        if ">" in alt:
            alt_clean = alt[1:-1]
        else:
            alt_clean = alt

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

# <LICENSE>
# Copyright (C) 2016-2026 VariantValidator Contributors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>
