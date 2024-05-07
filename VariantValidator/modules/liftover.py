# -*- coding: utf-8 -*-
"""
Liftover between genome builds is most accurate when mapping via a RefSeq transcript.
For intergenic regions, the process is more complex.
Lift position > Check bases > Lift back and confirm the original position
"""

# import modules
import vvhgvs.exceptions
import vvhgvs.sequencevariant
import logging
from . import seq_data
from . import hgvs_utils
from pyliftover import LiftOver
from Bio.Seq import Seq
import copy

# Pre compile variables
vvhgvs.global_config.formatting.max_ref_length = 1000000

logger = logging.getLogger(__name__)


def mystr(hgvs_nucleotide):
    hgvs_nucleotide_refless = hgvs_nucleotide.format({'max_ref_length': 0})
    return hgvs_nucleotide_refless


def liftover(hgvs_genomic, build_from, build_to, hn, reverse_normalizer, evm, validator,
             specify_tx=False, liftover_level=False, g_to_g=False, gap_map=False, vfo=False,
             specified_tx_variant=False):
    """
    Step 1, attempt to liftover using a common RefSeq transcript
    Step 2, attempt to liftover using PyLiftover.
    Lift position > Check bases > Lift back and confirm the original position
    :param hgvs_genomic:
    :param build_from:
    :param build_to:
    :param hn:
    :param reverse_normalizer:
    :param evm:
    :param validator: Validator obj
    :param specify_tx: Specify a specific transcript = False or str(transcript_ID)
    :param liftover_level: False or 'primary'
    :param g_to_g: True or False
    :param gap_map: True or VariantFormatter gap_map function passed (Required for VariantFormatter methods only)
    :param vfo: False or VariantFormatter VFO object passed
    :param specified_tx_variant: False or specific HGVS transcript object
    :return:
    """
    try:
        hgvs_genomic = validator.hp.parse(hgvs_genomic)
    except TypeError as e:
        logger.debug("Except passed, %s", e)

    # Create return dictionary
    lifted_response = {}
    # Check genome build type
    if 'GRC' in build_from:
        from_set = 'grc_chr'
        alt_from_set = 'ucsc_chr'
        if '37' in build_from:
            lo_from = 'hg19'
            alt_build_from = 'hg19'
        elif '38' in build_from:
            lo_from = 'hg38'
            alt_build_from = 'hg38'
        else:
            lo_from = ''
            alt_build_from = ''

    else:
        from_set = 'ucsc_chr'
        alt_from_set = 'grc_chr'
        if '19' in build_from:
            lo_from = 'hg19'
            alt_build_from = 'GRCh37'
        elif '38' in build_from:
            lo_from = 'hg38'
            alt_build_from = 'GRCh38'
        else:
            lo_from = ''
            alt_build_from = ''

    if 'GRC' in build_to:
        to_set = 'grc_chr'
        alt_to_set = 'ucsc_chr'
        if '37' in build_to:
            lo_to = 'hg19'
            alt_build_to = 'hg19'
        elif '38' in build_to:
            lo_to = 'hg38'
            alt_build_to = 'hg38'
        else:
            lo_to = ''
            alt_build_to = ''
    else:
        to_set = 'ucsc_chr'
        alt_to_set = 'grc_chr'
        if '19' in build_to:
            lo_to = 'hg19'
            alt_build_to = 'GRCh37'
        elif '38' in build_to:
            lo_to = 'hg38'
            alt_build_to = 'GRCh38'
        else:
            lo_to = ''
            alt_build_to = ''

    # populate the variant from data
    vcf = hgvs_utils.report_hgvs2vcf(hgvs_genomic, build_from, reverse_normalizer, validator.sf)

    # Create to and from dictionaries
    lifted_response[build_from.lower()] = {}
    lifted_response[build_from.lower()][hgvs_genomic.ac] = {
        'hgvs_genomic_description': mystr(hgvs_genomic),
        'vcf': {
           'chr': vcf[from_set],
           'pos': str(vcf['pos']),
           'ref': vcf['ref'],
           'alt': vcf['alt']
        }
    }
    lifted_response[alt_build_from.lower()] = {}
    lifted_response[alt_build_from.lower()][hgvs_genomic.ac] = {
        'hgvs_genomic_description': mystr(hgvs_genomic),
        'vcf': {
           'chr': vcf[alt_from_set],
           'pos': str(vcf['pos']),
           'ref': vcf['ref'],
           'alt': vcf['alt']
        }
    }

    # From dictionary currently blank
    lifted_response[build_to.lower()] = {}
    lifted_response[alt_build_to.lower()] = {}

    # Get a list of overlapping RefSeq transcripts
    # Note, due to 0 base positions in UTA (I think) occasionally tx will
    rts_list = validator.hdp.get_tx_for_region(hgvs_genomic.ac, 'splign', hgvs_genomic.posedit.pos.start.base - 1,
                                               hgvs_genomic.posedit.pos.end.base)  # - 1)
    rts_dict = {}
    tx_list = False
    if g_to_g is True:
        pass
    else:
        for tx_dat in rts_list:
            rts_dict[tx_dat[0]] = True
        if evm is not None:
            rts_list_2 = evm.relevant_transcripts(hgvs_genomic)
        else:
            rts_list_2 = []
        for tx_dat_2 in rts_list_2:
            rts_dict[tx_dat_2] = True
        if rts_dict != {}:
            tx_list = list(rts_dict.keys())

    # Try to liftover
    if tx_list is not False:
        selected = []
        # Liftover via a specific tx if it can be done!
        if specify_tx is not False:
            tx_list = [specify_tx]
        for tx in tx_list:
            # identify the first transcript if any
            options = validator.hdp.get_tx_mapping_options(tx)
            for op in options:
                sff = None
                sft = None
                if liftover_level is None:
                    continue
                if op[1].startswith('NC_'):
                    if build_to.startswith('GRC'):
                        sft = seq_data.to_chr_num_refseq(op[1], build_to)
                    if build_to.startswith('hg'):
                        sft = seq_data.to_chr_num_ucsc(op[1], build_to)
                    if build_from.startswith('GRC'):
                        sff = seq_data.to_chr_num_refseq(op[1], build_from)
                    if build_from.startswith('hg'):
                        sff = seq_data.to_chr_num_ucsc(op[1], build_from)
                    if sff is not None or sft is None:
                        selected.append([op[0], op[1], build_from, alt_build_from])
                    elif sff is None or sft is not None:
                        selected.append([op[0], op[1], build_to, alt_build_to])
                if liftover_level == 'primary':
                    continue
                else:
                    if op[1].startswith('NT_'):
                        if build_to.startswith('GRC'):
                            sft = seq_data.to_chr_num_refseq(op[1], build_to)
                        if build_to.startswith('hg'):
                            sft = seq_data.to_chr_num_ucsc(op[1], build_to)
                        if build_from.startswith('GRC'):
                            sff = seq_data.to_chr_num_refseq(op[1], build_from)
                        if build_from.startswith('hg'):
                            sff = seq_data.to_chr_num_ucsc(op[1], build_from)
                        if sff is not None or sft is None:
                            selected.append([op[0], op[1], build_from, alt_build_from])
                        elif sff is None or sft is not None:
                            selected.append([op[0], op[1], build_to, alt_build_to])
                    if op[1].startswith('NW_'):
                        if build_to.startswith('GRC'):
                            sft = seq_data.to_chr_num_refseq(op[1], build_to)
                        if build_to.startswith('hg'):
                            sft = seq_data.to_chr_num_ucsc(op[1], build_to)
                        if build_from.startswith('GRC'):
                            sff = seq_data.to_chr_num_refseq(op[1], build_from)
                        if build_from.startswith('hg'):
                            sff = seq_data.to_chr_num_ucsc(op[1], build_from)
                        if sff is not None or sft is None:
                            selected.append([op[0], op[1], build_from, alt_build_from])
                        elif sff is None or sft is not None:
                            selected.append([op[0], op[1], build_to, alt_build_to])

        # remove duplicate chroms
        filtered_1 = {}
        if selected:
            for chroms in selected:
                if chroms[1] not in list(filtered_1.keys()):
                    filtered_1[chroms[1]] = [chroms[0], chroms[2], chroms[3]]
            added_data = False
            for key, val in list(filtered_1.items()):
                try:
                    # Note, due to 0 base positions in UTA (I think) occasionally tx will
                    # be identified that cannot be mapped to.
                    # In this instance, do not mark added data as True
                    hgvs_tx = validator.vm.g_to_t(hgvs_genomic, val[0])
                    hgvs_alt_genomic = validator.vm.t_to_g(hgvs_tx, key)

                    # Gap compensation edit for the VariantFormatter pathway
                    if gap_map is not False and build_from not in val:
                        # Set genome assembly for gap mapping
                        get_assembly = seq_data.supported_for_mapping(key, "GRCh37")
                        if get_assembly is True:
                            map_to_assembly = "GRCh37"
                        get_assembly = seq_data.supported_for_mapping(key, "GRCh38")
                        if get_assembly is True:
                            map_to_assembly = "GRCh38"
                        try:
                            am_i_gapped = gap_map(specified_tx_variant, hgvs_alt_genomic, map_to_assembly, vfo)
                        except AttributeError:
                            if specified_tx_variant is None:
                                am_i_gapped = gap_map(hgvs_tx, hgvs_alt_genomic, map_to_assembly, vfo)

                        hgvs_alt_genomic = am_i_gapped["hgvs_genomic"]

                    alt_vcf = hgvs_utils.report_hgvs2vcf(hgvs_alt_genomic, build_to, reverse_normalizer, validator.sf)
                    alt_vcf_b = hgvs_utils.report_hgvs2vcf(hgvs_alt_genomic, build_from, reverse_normalizer,
                                                           validator.sf)

                    # Handle mitochondrial liftovers
                    if 'NC_012920.1' in hgvs_alt_genomic.ac or 'NC_001807.4' in hgvs_alt_genomic.ac:
                        hgvs_alt_genomic.type = "m"

                    # Add the to build dictionaries
                    if val[1] == build_to:
                        lifted_response[build_to.lower()][hgvs_alt_genomic.ac] = {
                            'hgvs_genomic_description': mystr(hgvs_alt_genomic),
                            'vcf': {
                                'chr': alt_vcf[to_set],
                                'pos': str(alt_vcf['pos']),
                                'ref': alt_vcf['ref'],
                                'alt': alt_vcf['alt']
                            }
                        }
                    if val[2] == alt_build_to:
                        lifted_response[alt_build_to.lower()][hgvs_alt_genomic.ac] = {
                            'hgvs_genomic_description': mystr(hgvs_alt_genomic),
                            'vcf': {
                                'chr': alt_vcf[alt_to_set],
                                'pos': str(alt_vcf['pos']),
                                'ref': alt_vcf['ref'],
                                'alt': alt_vcf['alt']
                            }
                        }
                    # Overwrite build from info as PAR may require additional info
                    if val[1] == build_from:
                        lifted_response[build_from.lower()][hgvs_alt_genomic.ac] = {
                            'hgvs_genomic_description': mystr(hgvs_alt_genomic),
                            'vcf': {
                                'chr': alt_vcf_b[to_set],
                                'pos': str(alt_vcf_b['pos']),
                                'ref': alt_vcf_b['ref'],
                                'alt': alt_vcf_b['alt']
                            }
                        }
                    if val[2] == alt_build_from:
                        lifted_response[alt_build_from.lower()][hgvs_alt_genomic.ac] = {
                            'hgvs_genomic_description': mystr(hgvs_alt_genomic),
                            'vcf': {
                                'chr': alt_vcf_b[alt_to_set],
                                'pos': str(alt_vcf_b['pos']),
                                'ref': alt_vcf_b['ref'],
                                'alt': alt_vcf_b['alt']
                            }
                        }

                    # Add gap warnings if found
                    try:
                        lifted_response["am_i_gapped"] = am_i_gapped
                    except UnboundLocalError:
                        pass

                    added_data = True

                except vvhgvs.exceptions.HGVSError:
                    continue

            if lifted_response != {} and added_data is not False:
                return lifted_response

    # Note: pyliftover uses the UCSC liftOver tool.
    # https://pypi.org/project/pyliftover/
    # Once validated, download the UCSC liftover files from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/
    # The structure of the following code comes from VV pymod, so need to create a list
    genome_builds = [build_to]

    # Create liftover vcf
    from_vcf = hgvs_utils.report_hgvs2vcf(hgvs_genomic, lo_from, reverse_normalizer, validator.sf)
    lo = LiftOver(lo_from, lo_to)

    # Fix the GRC CHR
    if from_vcf[from_set].startswith('chr'):
        liftover_list = lo.convert_coordinate(from_vcf[from_set], int(from_vcf['pos']))
    else:
        my_chrom = 'chr' + from_vcf[from_set]
        liftover_list = lo.convert_coordinate(my_chrom, int(from_vcf['pos']))

    # Create dictionary
    for lifted in liftover_list:
        chrom = lifted[0]
        pos = lifted[1]
        orientated = lifted[2]

        lifted_ref_bases = from_vcf['ref']
        lifted_alt_bases = from_vcf['alt']
        if hgvs_genomic.posedit.edit.type == "dup":
            # put complete original in ref and both copies of dup in alt
            lifted_ref_bases = lifted_ref_bases + lifted_alt_bases[1:]
            lifted_alt_bases = lifted_alt_bases + lifted_alt_bases[1:]

        # Inverted sequence
        if orientated != '+':
            continue

        # Find the accession
        accession = seq_data.to_accession(chrom, lo_to)
        if accession is None:
            wrn = 'Unable to identify an equivalent %s chromosome ID for %s' % (str(lo_to), str(chrom))
            logger.info(wrn)
            continue

        else:
            # Correct 37 to GRCh38 mito liftover - Applies when lifting from GRCh37 only!
            if "38" in build_to and "GRCh37" in build_from and accession == "NC_012920.1":
                mito_correction = True
                hgvs_lifted = hgvs_genomic

                # Fix the GRC CHR
                if from_vcf[from_set].startswith('chr'):
                    chrom = from_vcf[from_set]
                    pos = int(from_vcf['pos'])
                else:
                    chrom = 'chr' + from_vcf[from_set]
                    pos = int(from_vcf['pos'])
                liftback_list = [(chrom, pos, "+", "GRCh38"), (chrom, pos, "+", "GRCh37")]

            # Correct 38 to GRCh37 mito liftover - Applies when lifting from GRCh38/hg38 only!
            elif build_to == "GRCh37" and "38" in build_from and accession == "NC_001807.4":
                mito_correction = True
                hgvs_lifted = hgvs_genomic

                # Flag lifted genome build as hg19
                lst_liftover = list(lifted)
                lst_liftover[-1] = "hg19"
                lifted = tuple(lst_liftover)

                # Fix the GRC CHR
                if from_vcf[from_set].startswith('chr'):
                    chrom = from_vcf[from_set]
                    pos = int(from_vcf['pos'])
                else:
                    chrom = 'chr' + from_vcf[from_set]
                    pos = int(from_vcf['pos'])
                liftback_list = [(chrom, pos, "+", "GRCh37"), lifted]

                # Create the necessary hg19 mito hgvs
                m19_not_delins = accession + ':g.' + str(lifted[1]) + '_' + str(
                    (int(lifted[1]) - 1) + len(lifted_ref_bases)) + 'delins' + lifted_alt_bases
                m19_not_delins = str(m19_not_delins)
                m19_hgvs_not_delins = validator.hp.parse_hgvs_variant(m19_not_delins)
                try:
                    m19_hgvs_lifted = hn.normalize(m19_hgvs_not_delins)
                except vvhgvs.exceptions.HGVSError:
                    pass

            # Correct 37 to GRCh38 mito liftover - Applies when lifting from GRCh38/hg38 only!
            elif build_to == "GRCh38" and ("37" in build_from or "19" in build_from) and hgvs_genomic.ac \
                    == "NC_001807.4":
                mito_correction = True

                # Flag lifted genome build as hg19
                m19_hgvs_lifted = hgvs_genomic
                lst_liftover = list(lifted)
                lst_liftover[-1] = "GRCh37"
                lifted = tuple(lst_liftover)

                # Fix the GRC CHR
                if from_vcf[from_set].startswith('chr'):
                    chrom = from_vcf[from_set]
                    pos = int(from_vcf['pos'])
                else:
                    chrom = 'chr' + from_vcf[from_set]
                    pos = int(from_vcf['pos'])
                liftback_list = [lifted, (chrom, pos, "+", "hg19")]

                # Create the necessary hg19 mito hgvs
                m38_not_delins = accession + ':g.' + str(lifted[1]) + '_' + str(
                    (int(lifted[1]) - 1) + len(lifted_ref_bases)) + 'delins' + lifted_alt_bases
                m38_not_delins = str(m38_not_delins)
                m38_hgvs_not_delins = validator.hp.parse_hgvs_variant(m38_not_delins)
                try:
                    hgvs_lifted = hn.normalize(m38_hgvs_not_delins)
                except vvhgvs.exceptions.HGVSError:
                     pass

            else:
                mito_correction = False
                not_delins = accession + ':g.' + str(pos) + '_' + str(
                    (pos - 1) + len(lifted_ref_bases)) + 'delins' + lifted_alt_bases
                not_delins = str(not_delins)
                hgvs_not_delins = validator.hp.parse_hgvs_variant(not_delins)

                try:
                    hgvs_lifted = hn.normalize(hgvs_not_delins)
                except vvhgvs.exceptions.HGVSDataNotAvailableError:
                    continue
                except vvhgvs.exceptions.HGVSInvalidVariantError:
                    continue

                # Now try map back
                lo = LiftOver(lo_to, lo_from)
                # Lift back
                liftback_list = lo.convert_coordinate(chrom, pos)

            for lifted_back in liftback_list:
                # for hg19 and GRCh37 mito, we need to accign the origin build
                mito_build = False
                try:
                    if mito_correction is True:
                        mito_build = lifted_back[3]
                except IndexError:
                    pass

                # Set the necessary hg19 m. data
                if mito_build == 'hg19':
                    try:
                        hgvs_lifted = m19_hgvs_lifted
                    except NameError:
                        continue

                # Pull out the good guys!
                # Need to add chr to the from_set
                if not lifted_back[0].startswith('chr'):
                    my_from_chr = 'chr' + lifted_back[0]
                else:
                    my_from_chr = lifted_back[0]

                if lifted_back[0] == from_vcf[from_set] or lifted_back[0] == my_from_chr:
                    if lifted_back[1] == int(from_vcf['pos']) or mito_build == "hg19" or mito_build == "GRCh37":
                        for build in genome_builds:
                            vcf_dict = hgvs_utils.report_hgvs2vcf(
                                hgvs_lifted, build, reverse_normalizer, validator.sf)

                            # Handle mitochondrial liftovers
                            if 'NC_012920.1' in hgvs_lifted.ac or 'NC_001807.4' in hgvs_lifted.ac:
                                hgvs_lifted.type = "m"

                            # Compile the dictionary
                            if mito_build is False:
                                if build.startswith('GRC'):
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

                            else:
                                if mito_build.startswith('GRC'):
                                    lifted_response[mito_build.lower()][hgvs_lifted.ac] = {
                                        'hgvs_genomic_description': mystr(hgvs_lifted),
                                        'vcf': {'chr': vcf_dict['grc_chr'],
                                                'pos': str(vcf_dict['pos']),
                                                'ref': vcf_dict['ref'],
                                                'alt': vcf_dict['alt']
                                                }
                                    }

                                    if lifted_response["grch38"] == {}:
                                        lifted_response["grch38"][hgvs_lifted.ac] = {
                                            'hgvs_genomic_description': mystr(hgvs_lifted),
                                            'vcf': {'chr': vcf_dict['grc_chr'],
                                                    'pos': str(vcf_dict['pos']),
                                                    'ref': vcf_dict['ref'],
                                                    'alt': vcf_dict['alt']
                                                    }
                                        }
                                        lifted_response["hg38"][hgvs_lifted.ac] = {
                                            'hgvs_genomic_description': mystr(hgvs_lifted),
                                            'vcf': {'chr': 'chrM',
                                                    'pos': str(vcf_dict['pos']),
                                                    'ref': vcf_dict['ref'],
                                                    'alt': vcf_dict['alt']
                                                    }
                                        }

                                else:
                                    lifted_response[mito_build.lower()][hgvs_lifted.ac] = {
                                        'hgvs_genomic_description': mystr(hgvs_lifted),
                                        'vcf': {'chr': vcf_dict['ucsc_chr'],
                                                'pos': str(vcf_dict['pos']),
                                                'ref': vcf_dict['ref'],
                                                'alt': vcf_dict['alt']
                                                }
                                    }

    # Remove known bad lifts
    cp_lifted_response = copy.deepcopy(lifted_response)
    for key, val in cp_lifted_response.items():
        if key == "hg19" and "NC_012920.1" in val.keys():
            lifted_response.pop(key)
        elif key == "grch37" and "NC_001807.4" in val.keys():
            lifted_response[key].pop("NC_001807.4")

    return lifted_response

# <LICENSE>
# Copyright (C) 2016-2024 VariantValidator Contributors
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
