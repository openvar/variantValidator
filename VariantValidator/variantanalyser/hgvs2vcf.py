# -*- coding: utf-8 -*-
"""
hgvs2vcf.py

A variety of functions that convert parder hgvs objects into VCF component parts
Each function has a slightly difference emphasis

1. hgvs2vcf
Simple conversionwhich ensures identity is as 5 prime as possible by adding an extra 5
prime base. Necessary for most gap handling situations

2. report_hgvs2vcf
Used to report the Most true representation of the VCF i.e. 5 prime normalized but no
additional bases added. NOTE: no gap handling capabilities

3. pos_lock_hgvs2vcf
No normalization at all. No additional bases added. Simply returns an in-situ VCF

4. hard_right_hgvs2vcf and hard_left_hgvs2vcf
Designed specifically for gap handling.
hard left pushes as 5 prime as possible and adds additional bases
hard right pushes as 3 prime as possible and adds additional bases
"""

# Import  modules
import re
import copy
import supported_chromosome_builds as supportedChromosomeBuilds

# Import Biopython modules
from Bio.Seq import Seq


# Database connections and hgvs objects are now passed from VariantValidator.py

def hgvs2vcf(hgvs_genomic, primary_assembly, reverse_normalizer, sf):
    hgvs_genomic_variant = hgvs_genomic
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic_variant)
    # hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

    # Chr
    chr = supportedChromosomeBuilds.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac, primary_assembly)
    if chr is not None:
        pass
    else:
        chr = reverse_normalized_hgvs_genomic.ac

    if re.search('[GATC]+\=', str(reverse_normalized_hgvs_genomic.posedit)):
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif (re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('del', str(
            reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        alt_start = start - 1  #
        # Recover sequences
        ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), alt_start, end - 1)
        ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        # Assemble
        pos = start
        ref = ref_seq
        alt = ref_seq + ins_seq

        # Substitutions
    elif re.search('>', str(reverse_normalized_hgvs_genomic.posedit)):
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('ins', str(
            reverse_normalized_hgvs_genomic.posedit)):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
                # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, start)
        # Assemble
        pos = str(start)
        ref = pre_base + hgvs_del_seq
        alt = pre_base


    # inv
    elif re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
                # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        bs = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start - 1, adj_start)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = ins_seq
        if re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())


    # Delins
    elif (re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and re.search('ins', str(
            reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
                # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = vcf_del_seq[:1] + ins_seq

        # Duplications
    elif (re.search('dup', str(reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)  #
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2  #
        start = start - 1  #
        # Recover sequences
        dup_seq = reverse_normalized_hgvs_genomic.posedit.edit.ref
        vcf_ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_ref_seq
        alt = vcf_ref_seq + dup_seq
    else:
        chr = ''
        ref = ''
        alt = ''
        pos = ''

    # ensure as 5' as possible
    if chr != '' and pos != '' and ref != '' and alt != '':
        if len(ref) > 1:
            rsb = list(str(ref))
            if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity':
                pos = int(pos) - 1
                prev = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), pos - 1, pos)
                pos = str(pos)
                ref = prev + ref
                alt = prev + alt

    # Dictionary the VCF
    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': reverse_normalized_hgvs_genomic}
    return vcf_dict


def report_hgvs2vcf(hgvs_genomic, primary_assembly, reverse_normalizer, sf):
    hgvs_genomic_variant = hgvs_genomic

    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic_variant)
    # hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

    # UCSC Chr
    ucsc_chr = supportedChromosomeBuilds.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac, primary_assembly)
    if ucsc_chr is not None:
        pass
    else:
        ucsc_chr = reverse_normalized_hgvs_genomic.ac

    # GRC Chr
    grc_chr = supportedChromosomeBuilds.to_chr_num_refseq(reverse_normalized_hgvs_genomic.ac, primary_assembly)
    if grc_chr is not None:
        pass
    else:
        grc_chr = reverse_normalized_hgvs_genomic.ac

    # GRC Chr
    grc_chr = supportedChromosomeBuilds.to_chr_num_refseq(reverse_normalized_hgvs_genomic.ac, primary_assembly)
    if grc_chr is not None:
        pass
    else:
        grc_chr = reverse_normalized_hgvs_genomic.ac

    if re.search('[GATC]+\=', str(reverse_normalized_hgvs_genomic.posedit)):
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif (re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('del', str(
            reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        alt_start = start - 1  #
        # Recover sequences
        ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), alt_start, end - 1)
        ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        # Assemble
        pos = start
        ref = ref_seq
        alt = ref_seq + ins_seq

    # Substitutions
    elif re.search('>', str(reverse_normalized_hgvs_genomic.posedit)):
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('ins', str(
            reverse_normalized_hgvs_genomic.posedit)):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, start)
        # Assemble
        pos = str(start)
        ref = pre_base + hgvs_del_seq
        alt = pre_base


    # inv
    elif re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        bs = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start - 1, adj_start)
        # Assemble
        pos = str(start)
        # pos = str(start-1)
        # ref = bs + vcf_del_seq
        ref = vcf_del_seq
        alt = ins_seq
        if re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif (re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and re.search('ins', str(
            reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        # pos = str(start)
        # ref = vcf_del_seq
        # alt = vcf_del_seq[:1] + ins_seq
        pos = str(start + 1)
        ref = vcf_del_seq[1:]
        alt = ins_seq

    # Duplications
    elif (re.search('dup', str(reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)  #
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2  #
        start = start - 1  #
        # Recover sequences
        dup_seq = reverse_normalized_hgvs_genomic.posedit.edit.ref
        vcf_ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start + 1)
        ref = vcf_ref_seq[1:]
        alt = vcf_ref_seq[1:] + dup_seq
    else:
        chr = ''
        ref = ''
        alt = ''
        pos = ''

    # Dictionary the VCF
    vcf_dict = {'pos': str(pos), 'ref': ref, 'alt': alt, 'ucsc_chr': ucsc_chr, 'grc_chr': grc_chr,
                'normalized_hgvs': reverse_normalized_hgvs_genomic}
    return vcf_dict


def pos_lock_hgvs2vcf(hgvs_genomic, primary_assembly, reverse_normalizer, sf):
    # Replace reference manually
    if hgvs_genomic.posedit.edit.ref == '':
        hgvs_genomic.posedit.edit.ref = sf.fetch_seq(str(hgvs_genomic.ac), hgvs_genomic.posedit.pos.start.base - 1,
                                                     hgvs_genomic.posedit.pos.end.base)

    reverse_normalized_hgvs_genomic = hgvs_genomic
    if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity' and len(
            reverse_normalized_hgvs_genomic.posedit.edit.ref) == 0:
        reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(reverse_normalized_hgvs_genomic)

    # hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

    # Chr
    chr = supportedChromosomeBuilds.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac, primary_assembly)
    if chr is not None:
        pass
    else:
        chr = reverse_normalized_hgvs_genomic.ac

    if re.search('[GATC]+\=', str(reverse_normalized_hgvs_genomic.posedit)):
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif (re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('del', str(
            reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        alt_start = start - 1  #
        # Recover sequences
        ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), alt_start, end - 1)
        ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        # Assemble
        pos = start
        ref = ref_seq
        alt = ref_seq + ins_seq

    # Substitutions
    elif re.search('>', str(reverse_normalized_hgvs_genomic.posedit)):
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('ins', str(
            reverse_normalized_hgvs_genomic.posedit)):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, start)
        # Assemble
        pos = str(start)
        ref = pre_base + hgvs_del_seq
        alt = pre_base


    # inv
    elif re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        bs = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start - 1, adj_start)
        # Assemble
        pos = str(start)
        # pos = str(start-1)
        # ref = bs + vcf_del_seq
        ref = vcf_del_seq
        alt = ins_seq
        if re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif (re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and re.search('ins', str(
            reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = vcf_del_seq[:1] + ins_seq


    # Duplications
    elif (re.search('dup', str(reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)  #
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2  #
        start = start - 1  #
        # Recover sequences
        dup_seq = reverse_normalized_hgvs_genomic.posedit.edit.ref
        vcf_ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_ref_seq
        alt = vcf_ref_seq + dup_seq
    else:
        chr = ''
        ref = ''
        alt = ''
        pos = ''

    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': reverse_normalized_hgvs_genomic}
    return vcf_dict


def hard_right_hgvs2vcf(hgvs_genomic, primary_assembly, hn, sf):
    hgvs_genomic_variant = hgvs_genomic
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    normalized_hgvs_genomic = hn.normalize(hgvs_genomic_variant)

    # Chr
    chr = supportedChromosomeBuilds.to_chr_num_ucsc(normalized_hgvs_genomic.ac, primary_assembly)
    if chr is not None:
        pass
    else:
        chr = normalized_hgvs_genomic.ac

    if re.search('[GATC]+\=', str(normalized_hgvs_genomic.posedit)):
        pos = str(normalized_hgvs_genomic.posedit.pos.start)
        ref = normalized_hgvs_genomic.posedit.edit.ref
        alt = normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif (re.search('ins', str(normalized_hgvs_genomic.posedit)) and not re.search('del', str(
            normalized_hgvs_genomic.posedit))):
        end = int(normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(normalized_hgvs_genomic.posedit.pos.start.base)
        alt_start = start - 1  #
        # Recover sequences
        ref_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), alt_start, end - 1)
        ins_seq = normalized_hgvs_genomic.posedit.edit.alt
        # Assemble
        pos = start
        ref = ref_seq
        alt = ref_seq + ins_seq

    # Substitutions
    elif re.search('>', str(normalized_hgvs_genomic.posedit)):
        ref = normalized_hgvs_genomic.posedit.edit.ref
        alt = normalized_hgvs_genomic.posedit.edit.alt
        pos = str(normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif re.search('del', str(normalized_hgvs_genomic.posedit)) and not re.search('ins',
                                                                                  str(normalized_hgvs_genomic.posedit)):
        end = int(normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        try:
            ins_seq = normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), start, end)
        pre_base = sf.fetch_seq(str(normalized_hgvs_genomic.ac), adj_start, start)
        # Assemble
        pos = str(start)
        ref = pre_base + hgvs_del_seq
        alt = pre_base

    # inv
    elif re.search('inv', str(normalized_hgvs_genomic.posedit)):
        end = int(normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 1
        start = start
        try:
            ins_seq = normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), start, end)
        vcf_del_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), adj_start, end)
        bs = sf.fetch_seq(str(normalized_hgvs_genomic.ac), adj_start - 1, adj_start)
        # Assemble
        pos = str(start)
        # pos = str(start-1)
        # ref = bs + vcf_del_seq
        ref = vcf_del_seq
        alt = ins_seq
        if re.search('inv', str(normalized_hgvs_genomic.posedit)):
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif (re.search('del', str(normalized_hgvs_genomic.posedit)) and re.search('ins',
                                                                               str(normalized_hgvs_genomic.posedit))):
        end = int(normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(normalized_hgvs_genomic.posedit.pos.start.base - 1)
        adj_start = start - 1
        start = start
        try:
            ins_seq = normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), start, end)
        vcf_del_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = vcf_del_seq[:1] + ins_seq


    # Duplications
    elif (re.search('dup', str(normalized_hgvs_genomic.posedit))):
        end = int(normalized_hgvs_genomic.posedit.pos.end.base)  #
        start = int(normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2  #
        start = start - 1  #
        # Recover sequences
        dup_seq = normalized_hgvs_genomic.posedit.edit.ref
        vcf_ref_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_ref_seq
        alt = vcf_ref_seq + dup_seq
    else:
        chr = ''
        ref = ''
        alt = ''
        pos = ''

    # ADD SURROUNDING BASES
    if chr != '' and pos != '' and ref != '' and alt != '':
        # Add 2 post bases
        pos = int(pos)
        pre_end_pos = pos + len(ref)
        end_pos = pre_end_pos + 1
        post = sf.fetch_seq(str(normalized_hgvs_genomic.ac), pre_end_pos - 1, end_pos)
        ref = ref + post
        alt = alt + post

    # Dictionary the VCF
    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': normalized_hgvs_genomic}
    return vcf_dict


def hard_left_hgvs2vcf(hgvs_genomic, primary_assembly, reverse_normalizer, sf):
    hgvs_genomic_variant = hgvs_genomic
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic_variant)
    hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

    # Chr
    chr = supportedChromosomeBuilds.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac, primary_assembly)
    if chr is not None:
        pass
    else:
        chr = reverse_normalized_hgvs_genomic.ac

    if re.search('[GATC]+\=', str(reverse_normalized_hgvs_genomic.posedit)):
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif (re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('del', str(
            reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        alt_start = start - 1  #
        # Recover sequences
        ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), alt_start, end - 1)
        ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        # Assemble
        pos = start
        ref = ref_seq
        alt = ref_seq + ins_seq

    # Substitutions
    elif re.search('>', str(reverse_normalized_hgvs_genomic.posedit)):
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('ins', str(
            reverse_normalized_hgvs_genomic.posedit)):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, start)
        # Assemble
        pos = str(start)
        ref = pre_base + hgvs_del_seq
        alt = pre_base


    # inv
    elif re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        bs = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start - 1, adj_start)
        # Assemble
        pos = str(start)
        # pos = str(start-1)
        # ref = bs + vcf_del_seq
        ref = vcf_del_seq
        alt = ins_seq
        if re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif (re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and re.search('ins', str(
            reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = vcf_del_seq[:1] + ins_seq


    # Duplications
    elif (re.search('dup', str(reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)  #
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2  #
        start = start - 1  #
        # Recover sequences
        dup_seq = reverse_normalized_hgvs_genomic.posedit.edit.ref
        vcf_ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_ref_seq
        alt = vcf_ref_seq + dup_seq
    else:
        chr = ''
        ref = ''
        alt = ''
        pos = ''

    # ADD SURROUNDING BASES
    if chr != '' and pos != '' and ref != '' and alt != '':
        pre_pos = int(pos) - 1
        pre_pos
        prev = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), pre_pos - 1, pre_pos)
        pos = str(pre_pos)
        ref = prev + ref
        alt = prev + alt

    # Dictionary the VCF
    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': reverse_normalized_hgvs_genomic}
    return vcf_dict


def hgvs_ref_alt(hgvs_variant, sf):
    if re.search('[GATC]+\=', str(hgvs_variant.posedit)):
        ref = hgvs_variant.posedit.edit.ref
        alt = hgvs_variant.posedit.edit.ref

    # Insertions
    elif (re.search('ins', str(hgvs_variant.posedit)) and not re.search('del', str(hgvs_variant.posedit))):
        end = int(hgvs_variant.posedit.pos.end.base)
        start = int(hgvs_variant.posedit.pos.start.base)
        alt_start = start - 1  #
        # Recover sequences
        ref_seq = sf.fetch_seq(str(hgvs_variant.ac), alt_start, end)
        ins_seq = hgvs_variant.posedit.edit.alt
        # Assemble
        ref = ref_seq
        alt = ref_seq[:1] + ins_seq + ref_seq[-1:]

    # Substitutions
    elif re.search('>', str(hgvs_variant.posedit)):
        ref = hgvs_variant.posedit.edit.ref
        alt = hgvs_variant.posedit.edit.alt

    # Deletions
    elif re.search('del', str(hgvs_variant.posedit)) and not re.search('ins', str(hgvs_variant.posedit)):
        ref = hgvs_variant.posedit.edit.ref
        alt = ''

    # inv
    elif re.search('inv', str(hgvs_variant.posedit)):
        ref = hgvs_variant.posedit
        my_seq = Seq(ref)
        alt = str(my_seq.reverse_complement())

    # Delins
    elif (re.search('del', str(hgvs_variant.posedit)) and re.search('ins', str(hgvs_variant.posedit))):
        ref = hgvs_variant.posedit.edit.ref
        alt = hgvs_variant.posedit.edit.alt

    # Duplications
    elif (re.search('dup', str(hgvs_variant.posedit))):
        ref = hgvs_variant.posedit.edit.ref
        alt = hgvs_variant.posedit.edit.ref + hgvs_variant.posedit.edit.ref
    else:
        ref = ''
        alt = ''

    ref_alt_dict = {'ref': ref, 'alt': alt}
    return ref_alt_dict

# Copyright (C) 2018  Peter Causey-Freeman, University of Leicester
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