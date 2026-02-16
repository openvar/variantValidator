import pytest

from VariantValidator.modules import vcf_to_pvcf
from VariantValidator.modules.vcf_to_pvcf import VcfConversionError


# ============================================================
# Simple SNVs
# ============================================================

def test_snv_tab_delimited():
    line = "chr2\t1500000\t.\tA\tT\t.\tPASS\t."
    assert vcf_to_pvcf.vcf_to_shorthand(line) == "chr2-1500000-A-T"


# ============================================================
# Sequence-based indels (standard VCF representation)
# ============================================================

def test_sequence_deletion():
    # deletion of "T" at position 1500001
    line = "chr3\t1500000\t.\tAT\tA\t.\tPASS\t."
    assert vcf_to_pvcf.vcf_to_shorthand(line) == "chr3-1500000-AT-A"


def test_sequence_insertion():
    # insertion of "T" after A
    line = "chr3\t1500000\t.\tA\tAT\t.\tPASS\t."
    assert vcf_to_pvcf.vcf_to_shorthand(line) == "chr3-1500000-A-AT"


def test_larger_sequence_deletion():
    line = "chr5\t2000000\t.\tATGC\tA\t.\tPASS\t."
    assert vcf_to_pvcf.vcf_to_shorthand(line) == "chr5-2000000-ATGC-A"


# ============================================================
# Symbolic structural variants
# ============================================================

def test_symbolic_del_with_end():
    line = "chr1\t1000000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=1005000"
    assert vcf_to_pvcf.vcf_to_shorthand(line) == "chr1-1000000-1005000-DEL"


def test_symbolic_dup_with_svlen():
    line = "chr1\t2000000\t.\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;SVLEN=10000"
    assert vcf_to_pvcf.vcf_to_shorthand(line) == "chr1-2000000-2010000-DUP"


def test_symbolic_inv_with_end():
    line = "chr1\t3000000\t.\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END=3005000"
    assert vcf_to_pvcf.vcf_to_shorthand(line) == "chr1-3000000-3005000-INV"


def test_symbolic_dup_with_copy_number():
    line = "chr1\t2000000\t.\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;SVLEN=10000;CN=4"
    assert vcf_to_pvcf.vcf_to_shorthand(line) == "chr1-2000000-2010000-DUP[CN4]"


# ============================================================
# Error handling
# ============================================================

def test_header_line_raises():
    with pytest.raises(VcfConversionError):
        vcf_to_pvcf.vcf_to_shorthand("#CHROM\tPOS\tID\tREF\tALT")


def test_insufficient_columns():
    with pytest.raises(VcfConversionError):
        vcf_to_pvcf.vcf_to_shorthand("chr1\t100")


def test_invalid_pos():
    with pytest.raises(VcfConversionError):
        vcf_to_pvcf.vcf_to_shorthand("chr1\tNotAnInt\t.\tA\tT\t.\tPASS\t.")


def test_missing_end_for_symbolic_sv():
    with pytest.raises(VcfConversionError):
        vcf_to_pvcf.vcf_to_shorthand("chr1\t1000000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL")


def test_unknown_delimiter():
    with pytest.raises(VcfConversionError):
        vcf_to_pvcf.vcf_to_shorthand("chr1|1000000|.|A|T")


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


