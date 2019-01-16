# -*- coding: utf-8 -*-
"""
Liftover between genome builds is most accurate when mapping via a RefSeq transcript.
For intergenic regions, the process is more complex.
Lift position > Check bases > Lift back and confirm the original position
"""

# import modules
import hgvs.exceptions
import hgvs.sequencevariant
import re
import os
import vvChromosomes
import vvHGVS
from pyliftover import LiftOver
import warnings
from Bio.Seq import Seq

# Pre compile variables
hgvs.global_config.formatting.max_ref_length = 1000000

# Determine whether a liftover directory has been added to the environment
PYLIFTOVER_DIR = os.environ.get('PYLIFTOVER_DIR')

def mystr(hgvs_nucleotide):
    hgvs_nucleotide_refless = hgvs_nucleotide.format({'max_ref_length': 0})
    return hgvs_nucleotide_refless

def liftover(hgvs_genomic, build_from, build_to, hn, vm, vr, hdp, hp, reverse_normalizer, sf, evm):

    """
    :param hgvs_genomic: hgvs_object genomic description accession NC, NT, or NW. Not NG
    :param build_from:
    :param build_to:
    :return: lifted {}
    Step 1, attempt to liftover using a common RefSeq transcript
    """

    try:
        hgvs_genomic = hp.parse_hgvs_variant(hgvs_genomic)
    except TypeError:
        pass

    # Create return dictionary
    lifted_response = {}

    # Check genome build type
    if re.match('GRC', build_from):
        from_set = 'grc_chr'
        alt_from_set = 'ucsc_chr'
        if re.search('37', build_from):
            lo_from = 'hg19'
            alt_build_from = 'hg19'
        if re.search('38', build_from):
            lo_from = 'hg38'
            alt_build_from = 'hg38'

    else:
        from_set = 'ucsc_chr'
        alt_from_set = 'grc_chr'
        if re.search('19', build_from):
            lo_from = 'hg19'
            alt_build_from = 'GRCh37'
        if re.search('38', build_from):
            lo_from = 'hg38'
            alt_build_from = 'GRCh38'

    if re.match('GRC', build_to):
        to_set = 'grc_chr'
        alt_to_set = 'ucsc_chr'
        if re.search('37', build_to):
            lo_to = 'hg19'
            alt_build_to = 'hg19'
        if re.search('38', build_to):
            lo_to = 'hg38'
            alt_build_to = 'hg38'
    else:
        to_set = 'ucsc_chr'
        alt_to_set = 'grc_chr'
        if re.search('19', build_to):
            lo_to = 'hg19'
            alt_build_to = 'GRCh37'
        if re.search('38', build_to):
            lo_to = 'hg38'
            alt_build_to = 'GRCh38'

    # populate the variant from data
    vcf = hgvs2vcf.report_hgvs2vcf(hgvs_genomic, build_from, reverse_normalizer, sf)

    # Create to and from dictionaries
    lifted_response[build_from.lower()] = {}
    lifted_response[build_from.lower()][hgvs_genomic.ac] = {'hgvs_genomic_description': mystr(hgvs_genomic),
                                   'vcf': {
                                       'chr': vcf[from_set],
                                       'pos': str(vcf['pos']),
                                       'ref': vcf['ref'],
                                       'alt': vcf['alt']}
                                   }
    lifted_response[alt_build_from.lower()] = {}
    lifted_response[alt_build_from.lower()][hgvs_genomic.ac] = {'hgvs_genomic_description': mystr(hgvs_genomic),
                                   'vcf': {
                                       'chr': vcf[alt_from_set],
                                       'pos': str(vcf['pos']),
                                       'ref': vcf['ref'],
                                       'alt': vcf['alt']}
                                   }
    # From dictionary currently blank
    lifted_response[build_to.lower()] = {}
    lifted_response[alt_build_to.lower()] = {}

    # Get a list of overlapping RefSeq transcripts
    # Note, due to 0 base positions in UTA (I think) occasionally tx will
    rts_list = hdp.get_tx_for_region(hgvs_genomic.ac, 'splign', hgvs_genomic.posedit.pos.start.base - 1,
                                     hgvs_genomic.posedit.pos.end.base - 1)
    rts_dict = {}
    tx_list = False
    for tx_dat in rts_list:
        rts_dict[tx_dat[0]] = True
    rts_list_2 = evm.relevant_transcripts(hgvs_genomic)
    for tx_dat_2 in rts_list_2:
        rts_dict[tx_dat_2] = True
    if rts_dict != {}:
        tx_list = rts_dict.keys()

    # Try to liftover
    if tx_list is not False:
        selected = []
        for tx in tx_list:
            # identify the first transcript if any
            options = hdp.get_tx_mapping_options(tx)
            for op in options:
                if re.match('NC_', op[1]):
                    if re.match('GRC', build_to):
                        sfm = vvChromosomes.to_chr_num_refseq(op[1], build_to)
                    if re.match('hg', build_to):
                        sfm = vvChromosomes.to_chr_num_ucsc(op[1], build_to)
                    if sfm is not None:
                        selected.append([op[0], op[1]])
            for op in options:
                if re.match('NT_', op[1]):
                    if re.match('GRC', build_to):
                        sfm = vvChromosomes.to_chr_num_refseq(op[1], build_to)
                    if re.match('hg', build_to):
                        sfm = vvChromosomes.to_chr_num_ucsc(op[1], build_to)
                    if sfm is not None:
                        selected.append([op[0], op[1]])
            for op in options:
                if re.match('NW_', op[1]):
                    if re.match('GRC', build_to):
                        sfm = vvChromosomes.to_chr_num_refseq(op[1], build_to)
                    if re.match('hg', build_to):
                        sfm = vvChromosomes.to_chr_num_ucsc(op[1], build_to)
                    if sfm is not None:
                        selected.append([op[0], op[1]])

        # remove duplicate chroms
        filtered_1 = {}
        if selected:
            for chroms in selected:
                if chroms[1] in filtered_1.keys():
                    pass
                else:
                    filtered_1[chroms[1]] = chroms[0]
            added_data = False
            for key, val in filtered_1.iteritems():
                try:
                    # Note, due to 0 base positions in UTA (I think) occasionally tx will
                    # be identified that cannot be mapped to.
                    # In this instance, do not mark added data as True
                    hgvs_tx = vm.g_to_t(hgvs_genomic, val)
                    hgvs_alt_genomic = vm.t_to_g(hgvs_tx, key)
                    alt_vcf = vvHGVS.report_hgvs2vcf(hgvs_alt_genomic, build_to, reverse_normalizer, sf)

                    # Add the to build dictionaries
                    lifted_response[build_to.lower()][hgvs_alt_genomic.ac] = {
                                                                    'hgvs_genomic_description': mystr(hgvs_alt_genomic),
                                                                      'vcf': {
                                                                          'chr': alt_vcf[to_set],
                                                                          'pos': str(alt_vcf['pos']),
                                                                          'ref': alt_vcf['ref'],
                                                                          'alt': alt_vcf['alt']}
                                                                      }
                    lifted_response[alt_build_to.lower()][hgvs_alt_genomic.ac] = {
                                                                    'hgvs_genomic_description': mystr(hgvs_alt_genomic),
                                                                      'vcf': {
                                                                          'chr': alt_vcf[alt_to_set],
                                                                          'pos': str(alt_vcf['pos']),
                                                                          'ref': alt_vcf['ref'],
                                                                          'alt': alt_vcf['alt']}
                                                                      }
                    added_data = True
                except hgvs.exceptions.HGVSInvalidIntervalError as e:
                    continue

            if lifted_response != {} and added_data is not False:
                return lifted_response
            else:
                pass
        else:
            # liftover has failed
            pass

    """
    Step 2, attempt to liftover using PyLiftover. 
    Lift position > Check bases > Lift back and confirm the original position
    """

    # Note: pyliftover uses the UCSC liftOver tool.
    # https://pypi.org/project/pyliftover/
    # Once validated, download the UCSC liftover files from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/

    # The structure of the following code comes from VV pymod, so need to create a list
    genome_builds = [build_to]

    # Create liftover vcf
    from_vcf = vvHGVS.report_hgvs2vcf(hgvs_genomic, lo_from, reverse_normalizer, sf)

    if PYLIFTOVER_DIR is not None:
        lo_filename_to = PYLIFTOVER_DIR + "%sTo%s.over.chain" % (lo_from, lo_to)
        lo_filename_to = str(lo_filename_to.replace('Tohg', 'ToHg'))

        lo = LiftOver(lo_filename_to)
    else:
        lo = LiftOver(lo_from, lo_to)

    # Fix the GRC CHR
    if re.match('chr', from_vcf[from_set]):
        liftover_list = lo.convert_coordinate(from_vcf[from_set], int(from_vcf['pos']))
    else:
        my_chrom = 'chr' + from_vcf[from_set]
        liftover_list = lo.convert_coordinate(my_chrom, int(from_vcf['pos']))


    # Create dictionary
    primary_genomic_dicts = {}
    for lifted in liftover_list:
        chr = lifted[0]
        pos = lifted[1]
        orientated = lifted[2]

        lifted_ref_bases = from_vcf['ref']
        lifted_alt_bases = from_vcf['alt']

        # Inverted sequence
        if orientated != '+':
            my_seq = Seq(lifted_ref_bases)
            lifted_ref_bases = my_seq.reverse_complement()
            your_seq = Seq(lifted_alt_bases)
            lifted_alt_bases = your_seq.reverse_complement()
        accession = vvChromosomes.to_accession(chr, lo_to)
        if accession is None:
            wrn = 'Unable to identify an equivalent %s chromosome ID for %s' % (str(lo_to), str(chr))
            warnings.warn(wrn)
            continue
        else:
            not_delins = accession + ':g.' + str(pos) + '_' + str(
                (pos - 1) + len(lifted_ref_bases)) + 'del' + lifted_ref_bases + 'ins' + lifted_alt_bases
            hgvs_not_delins = hp.parse_hgvs_variant(not_delins)
            try:
                vr.validate(hgvs_not_delins)
            except hgvs.exceptions.HGVSError as e:
                warnings.warn(str(e))
                # Most likely incorrect bases
                continue
            else:
                hgvs_lifted = hn.normalize(hgvs_not_delins)
                # Now try map back
                if PYLIFTOVER_DIR is not None:
                    lo_filename_from = PYLIFTOVER_DIR + "%sTo%s.over.chain" % (lo_to, lo_from)

                    lo_filename_from = str(lo_filename_from.replace('Tohg', 'ToHg'))
                    lo = LiftOver(lo_filename_from)
                else:
                    lo = LiftOver(lo_to, lo_from)

                # Lift back
                liftback_list = lo.convert_coordinate(chr, pos)

                for lifted_back in liftback_list:
                    # Pull out the good guys!
                    # Need to add chr to the from_set
                    if not re.match('chr', lifted_back[0]):
                        my_from_chr = 'chr' + lifted_back[0]
                    else:
                        my_from_chr = lifted_back[0]

                    if lifted_back[0] == from_vcf[from_set] or lifted_back[0] == my_from_chr:
                        if lifted_back[1] == int(from_vcf['pos']):
                            for build in genome_builds:
                                vcf_dict = vvHGVS.report_hgvs2vcf(hgvs_lifted, build, reverse_normalizer, sf)
                                if re.match('GRC', build):
                                    lifted_response[build_to.lower()][hgvs_lifted.ac] = {
                                        'hgvs_genomic_description': mystr(hgvs_lifted),
                                        'vcf': {'chr': vcf_dict['grc_chr'],
                                                'pos': str(vcf_dict['pos']),
                                                'ref': vcf_dict['ref'],
                                                'alt': vcf_dict['alt']
                                                }
                                    }
                                    lifted_response[alt_build_to.lower()][hgvs_lifted.ac] = {
                                        'hgvs_genomic_description': mystr(hgvs_lifted),
                                        'vcf': {'chr': vcf_dict['ucsc_chr'],
                                                'pos': str(vcf_dict['pos']),
                                                'ref': vcf_dict['ref'],
                                                'alt': vcf_dict['alt']
                                                }
                                    }
                                else:
                                    lifted_response[build_to.lower()][hgvs_lifted.ac] = {
                                        'hgvs_genomic_description': mystr(hgvs_lifted),
                                        'vcf': {'chr': vcf_dict['ucsc_chr'],
                                                'pos': str(vcf_dict['pos']),
                                                'ref': vcf_dict['ref'],
                                                'alt': vcf_dict['alt']
                                                }
                                    }
                                    lifted_response[alt_build_to.lower()][hgvs_lifted.ac] = {
                                        'hgvs_genomic_description': mystr(hgvs_lifted),
                                        'vcf': {'chr': vcf_dict['grc_chr'],
                                                'pos': str(vcf_dict['pos']),
                                                'ref': vcf_dict['ref'],
                                                'alt': vcf_dict['alt']
                                                }
                                    }
    return lifted_response

# <LICENSE>
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


