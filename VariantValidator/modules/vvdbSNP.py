import os
import glob
import pysam
from configparser import ConfigParser

# ---- Configuration ----
CONFIG_DIR = os.path.join(os.path.expanduser('~'), '.variantvalidator')

# Load config
config = ConfigParser()
config.read(CONFIG_DIR)

try:
    BCF_DIR = config['dbSNP']['location']  # directory where BCFs are stored
except KeyError:
    raise KeyError('dbSNP location not found in config file')

# ---- Automatically register the latest BCFs ----
BCF_FILES = {}

# GRCh37
files_37 = sorted(glob.glob(os.path.join(BCF_DIR, "vvdbSNP_GRCh37_*.bcf")))
if not files_37:
    raise FileNotFoundError(f"No BCF files matching 'vvdbSNP_GRCh37_*.bcf' found in {BCF_DIR}")
BCF_FILES["GRCh37"] = files_37[-1]

# GRCh38
files_38 = sorted(glob.glob(os.path.join(BCF_DIR, "vvdbSNP_GRCh38_*.bcf")))
if not files_38:
    raise FileNotFoundError(f"No BCF files matching 'vvdbSNP_GRCh38_*.bcf' found in {BCF_DIR}")
BCF_FILES["GRCh38"] = files_38[-1]

# ---- Lookup functions ----

def chrpos_to_rsid(chrom: str, pos: str, ref: str, alt: str, genome_build: str):
    """
    Lookup dbSNP rsID from chromosome, position, ref, alt for a given genome build.
    pos is passed as a string but converted internally.
    """
    pos_int = int(pos)
    bcf_path = BCF_FILES.get(genome_build)
    if not bcf_path:
        raise ValueError(f"No BCF registered for genome_build: {genome_build}")

    with pysam.VariantFile(bcf_path) as bcf:
        for rec in bcf.fetch(chrom, pos_int - 1, pos_int):
            if rec.ref == ref and alt in rec.alts:
                return rec.id
    return None

def rsid_to_chrpos(rsid: str, genome_build: str):
    """
    Lookup chromosome, position, ref, alt for a dbSNP rsID for a given genome build.
    Returns a list of tuples: (chrom, pos, ref, alt)
    """
    bcf_path = BCF_FILES.get(genome_build)
    if not bcf_path:
        raise ValueError(f"No BCF registered for genome_build: {genome_build}")

    results = []
    with pysam.VariantFile(bcf_path) as bcf:
        for rec in bcf.fetch():
            if rec.id == rsid:
                for alt in rec.alts:
                    results.append((rec.chrom, rec.pos, rec.ref, alt))
    return results

# ---- Test / Example usage ----
# ---- Test / Example usage ----
if __name__ == "__main__":
    import time

    print("Registered BCF files:")
    for gb, path in BCF_FILES.items():
        print(f"  {gb}: {path}")
    print("Time:", time.strftime("%Y-%m-%d %H:%M:%S"), "\n")

    # ---- Test GRCh37 using actual variant we have ----
    chrom = "NC_000001.10"
    pos = "10001"
    ref = "T"
    alt = "A"
    genome_build = "GRCh37"

    start = time.time()
    rsid = chrpos_to_rsid(chrom, pos, ref, alt, genome_build)
    print(f"chrpos_to_rsid({chrom}, {pos}, {ref}>{alt}, {genome_build}) -> {rsid}")
    print("Time:", time.strftime("%Y-%m-%d %H:%M:%S"), f"Elapsed: {time.time() - start:.2f}s\n")

    # rsID -> chr/pos
    test_rsid = "rs1570391677"
    genome_build = "GRCh37"

    start = time.time()
    results = rsid_to_chrpos(test_rsid, genome_build)
    print(f"rsid_to_chrpos({test_rsid}, {genome_build}) -> {results}")
    print("Time:", time.strftime("%Y-%m-%d %H:%M:%S"), f"Elapsed: {time.time() - start:.2f}s")


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
