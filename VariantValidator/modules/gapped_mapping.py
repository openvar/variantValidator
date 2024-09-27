import copy
import re
import logging
import vvhgvs.exceptions
from . import utils as fn
from . import hgvs_utils
from . import seq_data

logger = logging.getLogger(__name__)


class GapMapper(object):

    def __init__(self, variant, validator):
        """
        Sets initial values
        :param variant: variant.Variant()
        :param validator: Validator()
        """
        self.variant = variant
        self.validator = validator
        self.gapped_transcripts = ''
        self.auto_info = ''
        self.orientation = None
        self.hgvs_genomic_possibilities = []
        self.disparity_deletion_in = []
        self.hgvs_genomic_5pr = None
        self.tx_hgvs_not_delins = None

    def make_gap_warnings(self, tx_ac, gen_ac, primary_assembly, message=None):

        # Look at Cigar strings and calculate the gap size and location
        tx_exons = self.validator.hdp.get_tx_exons(tx_ac, gen_ac, alt_aln_method=self.validator.alt_aln_method)

        # Locate all the gaps
        gap_in_alignment = []
        for exon in tx_exons:
            data_required = [exon[0],
                             exon[1],
                             exon[3],
                             int(exon[5]) + 1,
                             int(exon[6]),
                             int(exon[7]),
                             int(exon[8]) + 1,
                             exon[9]]

            if "I" in data_required[-1] or "D" in data_required[-1]:
                gap_in_alignment.append(data_required)

        # Create warnings
        gap_information_dict = {"gapped_alignment_warning": "",
                                "auto_info": ""
                                }
        # Identify gaps
        if gap_in_alignment is not []:
            found_gaps = []
            for gap_loc in gap_in_alignment:
                cigar = gap_loc[-1]
                split_my_cigar = cigar.replace("=", "=:")
                split_my_cigar = split_my_cigar.replace("I", "I:")
                split_my_cigar = split_my_cigar.replace("D", "D:")
                split_my_cigar = split_my_cigar.replace("X", "X:")
                split_my_cigar = split_my_cigar.split(":")

                # Get annotation
                tx_exon_start = int(gap_loc[3])
                tx_annotation = self.validator.hdp.get_tx_identity_info(gap_loc[0])
                try:
                    cds_start = int(tx_annotation[3]) + 1
                    cds_end = int(tx_annotation[4])
                    c_tx_exon_start = tx_exon_start - cds_start
                except TypeError:
                    c_tx_exon_start = tx_exon_start

                # Get all gap locations in the transcript split into ins and del
                for gap in split_my_cigar:
                    gap = gap.replace("X", "=")

                    if "=" in gap and "I" not in gap and "D" not in gap:
                        if "NM_" in tx_ac:
                            c_tx_exon_start = c_tx_exon_start + int(gap.split("=")[0])
                        else:
                            c_tx_exon_start = c_tx_exon_start + int(gap.split("=")[0])

                    elif "D" in gap:
                        if "NM_" in tx_ac:
                            pos_n_len = ["c." + str(c_tx_exon_start) + "_"
                                         + str(c_tx_exon_start + int(gap.split("D")[0])+1),
                                         str(gap.split("D")[0]) + " extra bases"
                                         ]
                        else:
                            pos_n_len = ["n." + str(c_tx_exon_start) + "_"
                                         + str(c_tx_exon_start + int(gap.split("D")[0])+1),
                                         str(gap.split("D")[0]) + " extra bases"
                                         ]

                        found_gaps.append(pos_n_len[1] + " between " + pos_n_len[0])
                        c_tx_exon_start = c_tx_exon_start + int(gap.split("D")[0])

                    elif "I" in gap:
                        if "NM_" in tx_ac:
                            pos_n_len = ["c." + str(c_tx_exon_start) + "_"
                                         + str(c_tx_exon_start + 1),
                                         str(gap.split("I")[0]) + " fewer bases"
                                        ]

                        elif "NR_" in tx_ac:
                            pos_n_len = ["n." + str(c_tx_exon_start) + "_"
                                         + str(c_tx_exon_start + 1),
                                         str(gap.split("I")[0]) + " fewer bases"
                                        ]

                        found_gaps.append(pos_n_len[1] + " between " + pos_n_len[0])
                        c_tx_exon_start = c_tx_exon_start + int(gap.split("I")[0])

            # Correct for UTR variants
            if "NM_" in tx_ac:
                cp_found_gaps = copy.copy(found_gaps)
                found_gaps = []
                for each_found in cp_found_gaps:
                    crds = each_found.split("c.")[-1]
                    start = int(crds.split("_")[0])

                    # 3 prime UTR
                    if start+cds_start >= cds_end:
                        utr_3 = start - cds_end
                        utr_3_pos = "*%s_*%s" % (str(utr_3 + cds_start), str(utr_3+1+cds_start))
                        each_found = each_found.replace(crds, utr_3_pos)
                        found_gaps.append(each_found)

                    # 5 prime UTR
                    elif start+cds_start <= cds_start:
                        found_gaps.append(each_found)

                    # CDS gap
                    else:
                        found_gaps.append(each_found)

            # Create the warnings
            if message is not None:
                gap_string = message
            else:
                gap_string = ", and ".join(found_gaps)

            gapped_alignment_warning = """Submitted description does not represent a true variant because 
it is an artefact of aligning %s with %s (genome build %s)""" % (tx_ac, gen_ac, primary_assembly)

            auto_info = """%s contains %s than %s""" % (tx_ac, gap_string, gen_ac)

            gap_information_dict["gapped_alignment_warning"] = gapped_alignment_warning. \
                replace("\n", "")
            gap_information_dict["auto_info"] = auto_info.replace("\n", "")

        return gap_information_dict

    def gapped_g_to_c(self, rel_var, select_transcripts_dict):
        """
        Gap aware projection from g. to c.
        """
        # RefSeq or Ensembl?
        expanded_genomic_for_ensembl = False
        if self.validator.alt_aln_method == 'genebuild':
            # Expand the genomic variant to include the flanking bases as a delins
            reverse_normalized_hgvs_genomic = self.validator.reverse_hn.normalize(self.variant.hgvs_genomic)

            # VCF
            vcf_dict = hgvs_utils.hgvs2vcf(reverse_normalized_hgvs_genomic, self.variant.primary_assembly,
                                           self.variant.reverse_normalizer, self.validator.sf,
                                           extra_flank_bases=4)
            pos = vcf_dict['pos']
            ref = vcf_dict['ref']
            alt = vcf_dict['alt']

            # Generate an end position
            end = str(int(pos) + len(ref) - 1)
            pos = str(pos)
            expanded_genomic_for_ensembl = self.validator.hp.parse_hgvs_variant(reverse_normalized_hgvs_genomic.ac
                                                                                + ':' +
                                                                   reverse_normalized_hgvs_genomic.type + '.' + pos +
                                                                   '_' + end +
                                                                   'del' + ref + 'ins' + alt)

        # Set variables for problem specific warnings
        gapped_alignment_warning = ''
        corrective_action_taken = ''
        self.gapped_transcripts = ''
        self.auto_info = ''
        self.disparity_deletion_in = []

        # Create a pseudo VCF so that normalization can be applied and a delins can be generated
        hgvs_genomic_variant = self.variant.hgvs_genomic
        # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic_variant)
        self.hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

        # VCF
        vcf_dict = hgvs_utils.hgvs2vcf(reverse_normalized_hgvs_genomic, self.variant.primary_assembly,
                                       self.variant.reverse_normalizer, self.validator.sf)
        pos = vcf_dict['pos']
        ref = vcf_dict['ref']
        alt = vcf_dict['alt']

        # Generate an end position
        end = str(int(pos) + len(ref) - 1)
        pos = str(pos)

        # take a look at the input genomic variant for potential base salvage
        stash_ac = vcf_dict['chr']
        stash_pos = int(vcf_dict['pos'])
        stash_ref = vcf_dict['ref']
        stash_alt = vcf_dict['alt']
        stash_end = end
        stash_input = str(self.variant.post_format_conversion)
        # Re-Analyse genomic positions
        if 'NG_' in str(self.variant.hgvs_formatted):
            c = rel_var[0]
            if hasattr(c.posedit.edit, 'ref') and c.posedit.edit.ref is not None:
                c.posedit.edit.ref = c.posedit.edit.ref.upper()
            if hasattr(c.posedit.edit, 'alt') and c.posedit.edit.alt is not None:
                c.posedit.edit.alt = c.posedit.edit.alt.upper()
            stash_input = self.validator.myevm_t_to_g(c, self.variant.no_norm_evm,
                                                      self.variant.primary_assembly, self.variant.hn)
        if 'NC_' in str(stash_input) or 'NT_' in str(stash_input) or 'NW_' in str(stash_input):
            try:
                hgvs_stash = self.validator.hp.parse_hgvs_variant(stash_input)
            except:
                hgvs_stash = stash_input
            if hasattr(hgvs_stash.posedit.edit, 'ref') and hgvs_stash.posedit.edit.ref is not None:
                hgvs_stash.posedit.edit.ref = hgvs_stash.posedit.edit.ref.upper()
            if hasattr(hgvs_stash.posedit.edit, 'alt') and hgvs_stash.posedit.edit.alt is not None:
                hgvs_stash.posedit.edit.alt = hgvs_stash.posedit.edit.alt.upper()

            # MAKE A NO NORM HGVS2VCF
            stash_dict = hgvs_utils.pos_lock_hgvs2vcf(hgvs_stash, self.variant.primary_assembly,
                                                      self.variant.reverse_normalizer, self.validator.sf)
            stash_ac = hgvs_stash.ac
            stash_pos = int(stash_dict['pos'])
            stash_ref = stash_dict['ref']
            stash_alt = stash_dict['alt']
            # Generate an end position
            stash_end = str(stash_pos + len(stash_ref) - 1)

        # Store a not real deletion insertion
        stored_hgvs_not_delins = self.validator.hp.parse_hgvs_variant(self.hgvs_genomic_5pr.ac + ':' +
            self.hgvs_genomic_5pr.type + '.' + pos + '_' + end + 'del' + ref + 'ins' + alt)
        assert stored_hgvs_not_delins != ''
        stash_hgvs_not_delins = self.validator.hp.parse_hgvs_variant(stash_ac + ':' +
            self.hgvs_genomic_5pr.type + '.' + str(stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)
        pre_loop_stash_hgvs_not_delins = copy.copy(stash_hgvs_not_delins)


        # make an empty rel_var
        nw_rel_var = []

        # loop through rel_var and amend where required
        for var in rel_var:
            stash_hgvs_not_delins = pre_loop_stash_hgvs_not_delins

            # Blank the self.tx_hgvs_not_delins
            self.tx_hgvs_not_delins = None

            # Store the current hgvs:c. description
            try:
                saved_hgvs_coding = self.validator.hp.parse_hgvs_variant(var)
                original_var = self.validator.hp.parse_hgvs_variant(var)
            except TypeError:
                saved_hgvs_coding = var
                original_var = var
            except vvhgvs.exceptions.HGVSInvalidVariantError:
                saved_hgvs_coding = var
                original_var = var

            # Remove un-selected transcripts
            if self.validator.select_transcripts != 'all' and self.validator.select_transcripts != 'raw' and \
                    "select" not in self.validator.select_transcripts and \
                    "mane" not in self.validator.select_transcripts and "refseqgene" not in \
                    self.validator.select_transcripts:
                tx_ac = saved_hgvs_coding.ac
                # If it's in the selected tx dict, keep it
                if tx_ac.split('.')[0] in list(select_transcripts_dict.keys()):
                    pass
                # If not get rid of it!
                else:
                    continue

            # Filter for Select transcripts only
            elif self.validator.select_transcripts == "select":
                tx_ac = saved_hgvs_coding.ac
                annotation = self.validator.db.get_transcript_annotation(tx_ac)
                if '"select": "MANE"' in annotation or '"select": "RefSeq"' in annotation or \
                        '"select": "Ensembl"' in annotation:
                    pass
                else:
                    continue

            # Filter for MANE transcripts only
            elif self.validator.select_transcripts == "mane":
                tx_ac = saved_hgvs_coding.ac
                annotation = self.validator.db.get_transcript_annotation(tx_ac)
                if '"mane_select": true' in annotation or '"mane_plus_clinical": true' in annotation:
                    pass
                else:
                    continue

            # Filter for mane Select transcripts only
            elif self.validator.select_transcripts == "mane_select":
                tx_ac = saved_hgvs_coding.ac
                annotation = self.validator.db.get_transcript_annotation(tx_ac)
                if '"mane_select": true' in annotation:
                    pass
                else:
                    continue

            # Filter for RefSeq Select transcripts only
            elif self.validator.select_transcripts == "refseq_select":
                tx_ac = saved_hgvs_coding.ac
                annotation = self.validator.db.get_transcript_annotation(tx_ac)
                if '"refseq_select": true' in annotation:
                    pass
                else:
                    continue

            # Filter for ensembl Select transcripts only
            elif self.validator.select_transcripts == "ensembl_select":
                tx_ac = saved_hgvs_coding.ac
                annotation = self.validator.db.get_transcript_annotation(tx_ac)
                if '"ensembl_select": true' in annotation:
                    pass
                else:
                    continue

            # Only apply to known gapped alignment genes
            symbol = self.validator.db.get_gene_symbol_from_transcript_id(saved_hgvs_coding.ac)
            if seq_data.gap_black_list(symbol) is True:

                # Applies to ensembl only
                if expanded_genomic_for_ensembl is not False:
                    try:
                        hgvs_refreshed_variant = self.validator.vm.g_to_t(expanded_genomic_for_ensembl,
                                                                      saved_hgvs_coding.ac,
                                                                      alt_aln_method=self.validator.alt_aln_method)
                    except vvhgvs.exceptions.HGVSInvalidIntervalError as e:
                        if "start or end or both are beyond the bounds of transcript record" in str(e):
                            continue
                    try:
                        hgvs_refreshed_variant = self.validator.vm.n_to_c(hgvs_refreshed_variant)
                    except vvhgvs.exceptions.HGVSError:
                        pass

                    # Get the ref length difference
                    genomic_ref_len = len(expanded_genomic_for_ensembl.posedit.edit.ref)
                    transcript_ref_len = len(hgvs_refreshed_variant.posedit.edit.ref)
                    if genomic_ref_len != transcript_ref_len:
                        message = None
                        if genomic_ref_len > transcript_ref_len:
                            gap_length = genomic_ref_len - transcript_ref_len
                            message = f"{gap_length} fewer bases"
                        elif genomic_ref_len < transcript_ref_len:
                            gap_length = transcript_ref_len - genomic_ref_len
                            message = f"{gap_length} extra bases"
                        else:
                            message = None

                        if message is not None:
                            gap_warnings = self.make_gap_warnings(hgvs_refreshed_variant.ac,
                                                                  expanded_genomic_for_ensembl.ac,
                                                                  self.variant.primary_assembly,
                                                                  message=message)
                            gapped_alignment_warning = gap_warnings["gapped_alignment_warning"]
                            if self.auto_info != "":
                                self.auto_info = self.auto_info  # + ", and " + gap_warnings["auto_info"]
                            else:
                                self.auto_info = gap_warnings["auto_info"]

                        # Will filter out intronic variants since intronic variants will not normalize
                        try:
                            hgvs_refreshed_variant = self.validator.genebuild_normalizer_cross.normalize(
                                hgvs_refreshed_variant)
                        except vvhgvs.exceptions.HGVSError:
                            nw_rel_var.append(saved_hgvs_coding)
                        else:
                            nw_rel_var.append(hgvs_refreshed_variant)
                        continue
                    else:
                        nw_rel_var.append(saved_hgvs_coding)

                # This next section is looking for exonic gaps so cannot be applied to intronic positions
                needs_a_push = False
                merged_variant = False

                if "+" not in str(saved_hgvs_coding.posedit.pos) and "-" not in str(saved_hgvs_coding.posedit.pos):
                    """
                    Directly search for gaps using vcf hard_pushing left
                    """
                    try:
                        vcf__dict = hgvs_utils.hard_left_hgvs2vcf(saved_hgvs_coding,
                                                                  self.variant.primary_assembly,
                                                                  self.variant.hn,
                                                                  self.variant.reverse_normalizer,
                                                                  self.validator.sf,
                                                                  saved_hgvs_coding.ac,
                                                                  self.validator.hdp,
                                                                  self.validator.alt_aln_method,
                                                                  self.validator.hp,
                                                                  self.validator.vm,
                                                                  self.validator.merge_hgvs_3pr,
                                                                  genomic_ac=hgvs_genomic_variant.ac)

                        if vcf__dict['needs_a_push'] is True:
                            needs_a_push = True
                            merged_variant = vcf__dict['merged_variant']
                            if merged_variant is not False:
                                try:
                                    merged_variant = self.validator.vm.n_to_c(merged_variant)
                                except TypeError:
                                    pass
                                except vvhgvs.exceptions.HGVSInvalidVariantError:
                                    pass
                                except vvhgvs.exceptions.HGVSUsageError:
                                    pass

                    except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                        pass
                    except vvhgvs.exceptions.HGVSDataNotAvailableError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass

                    """
                    Directly search for gaps using vcf hard_pushing right
                    """
                    try:
                        vcf__dict = hgvs_utils.hard_right_hgvs2vcf(saved_hgvs_coding,
                                                                   self.variant.primary_assembly,
                                                                   self.variant.hn,
                                                                   self.variant.reverse_normalizer,
                                                                   self.validator.sf,
                                                                   saved_hgvs_coding.ac,
                                                                   self.validator.hdp,
                                                                   self.validator.alt_aln_method,
                                                                   self.validator.hp,
                                                                   self.validator.vm,
                                                                   self.validator.merge_hgvs_3pr,
                                                                   genomic_ac=hgvs_genomic_variant.ac,)

                        if vcf__dict['needs_a_push'] is True:
                            needs_a_push = True
                            merged_variant = vcf__dict['merged_variant']
                            if merged_variant is not False:
                                try:
                                    merged_variant = self.validator.vm.n_to_c(merged_variant)
                                except TypeError:
                                    pass
                                except vvhgvs.exceptions.HGVSInvalidVariantError:
                                    pass

                    except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                        pass
                    except vvhgvs.exceptions.HGVSDataNotAvailableError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass

                    # Collect the hard_pushed variant information and adjust the variants accordingly
                    if needs_a_push is not False:
                        if merged_variant is not False:
                            saved_hgvs_coding = merged_variant
                            stash_hgvs_not_delins = self.validator.vm.t_to_g(saved_hgvs_coding,
                                                                             hgvs_genomic_variant.ac,
                                                                             alt_aln_method=self.validator.alt_aln_method)

                            # The merged variant may have created an ins or a del
                            if stash_hgvs_not_delins.posedit.edit.type == "del":
                                stash_hgvs_not_delins.posedit.edit.alt = ""
                            if stash_hgvs_not_delins.posedit.edit.type == "ins":
                                get_ref = copy.deepcopy(stash_hgvs_not_delins)
                                get_ref.posedit.edit.ref = ''
                                get_ref.posedit.edit.alt = ''
                                get_ref = self.variant.hn.normalize(get_ref)
                                ref_bases = get_ref.posedit.edit.ref
                                stash_hgvs_not_delins.posedit.edit.ref = ref_bases
                                stash_hgvs_not_delins.posedit.edit.alt = ref_bases[0] \
                                                                         + stash_hgvs_not_delins.posedit.edit.alt \
                                                                         + ref_bases[1]

                # Get orientation of the gene wrt genome and a list of exons mapped to the genome
                ori = self.validator.tx_exons(tx_ac=saved_hgvs_coding.ac, alt_ac=self.hgvs_genomic_5pr.ac,
                                              alt_aln_method=self.validator.alt_aln_method)
                try:
                    self.orientation = int(ori[0]['alt_strand'])
                except TypeError:
                    continue

                # Set intronic params
                intronic_variant = 'false'
                hgvs_seek_var = self.get_hgvs_seek_var(self.variant.hgvs_genomic, saved_hgvs_coding)

                if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                        saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                        hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                        saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset):
                    pass
                else:
                    hgvs_seek_var = saved_hgvs_coding

                try:
                    self.variant.hn.normalize(hgvs_seek_var)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                    error = str(e)
                    if 'Normalization of intronic variants is not supported' in error or \
                            'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                            intronic_variant = 'hard_fail'
                        else:
                            # Double check to see whether the variant is actually intronic?
                            for exon in ori:
                                genomic_start = int(exon['alt_start_i'])
                                genomic_end = int(exon['alt_end_i'])
                                if genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end and \
                                        genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end:
                                    intronic_variant = 'false'
                                    break
                                else:
                                    intronic_variant = 'true'
                except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                    if "insertion length must be 1" in str(e):
                        pass

                if re.search(r'\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search(r'\d+-', str(
                        hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+\+', str(
                        hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+-', str(hgvs_seek_var.posedit.pos)):

                    # Double check to see whether the variant is actually intronic?
                    for exon in ori:
                        genomic_start = int(exon['alt_start_i'])
                        genomic_end = int(exon['alt_end_i'])
                        if (genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                                genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                            intronic_variant = 'false'
                            break
                        else:
                            intronic_variant = 'true'

                # If exonic, process
                if intronic_variant != 'true' and intronic_variant != 'hard_fail':
                    # Attempt to find gaps in reference sequence by catching disparity in genome length and
                    # overlapping transcript lengths
                    self.disparity_deletion_in = ['false', 'false']
                    hgvs_not_delins = ''

                    if stored_hgvs_not_delins != '':
                        # Refresh hgvs_not_delins from stored_hgvs_not_delins
                        try:
                            hgvs_not_delins = self.dup_ins_5prime_shift(stored_hgvs_not_delins, saved_hgvs_coding)
                        except vvhgvs.exceptions.HGVSInvalidIntervalError:
                            hgvs_not_delins = stored_hgvs_not_delins

                        try:
                            self.tx_hgvs_not_delins = self.variant.no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                      saved_hgvs_coding.ac)
                        except vvhgvs.exceptions.HGVSInvalidIntervalError as e:
                            if "start or end or both are beyond the bounds of transcript record" in str(e):
                                self.tx_hgvs_not_delins = saved_hgvs_coding
                            else:
                                self.tx_hgvs_not_delins = self.variant.no_norm_evm.g_to_n(self.hgvs_genomic_5pr,
                                                                                          saved_hgvs_coding.ac)
                        except vvhgvs.exceptions.HGVSError as e:
                            if str(e) == 'start or end or both are beyond the bounds of transcript record':
                                self.tx_hgvs_not_delins = saved_hgvs_coding

                        # Create normalized version of tx_hgvs_not_delins
                        rn_tx_hgvs_not_delins = copy.deepcopy(self.tx_hgvs_not_delins)

                        # Check for +ve base and adjust
                        if ('+' in str(rn_tx_hgvs_not_delins.posedit.pos.start) or '-' in
                            str(rn_tx_hgvs_not_delins.posedit.pos.start)) and (
                                '+' in str(rn_tx_hgvs_not_delins.posedit.pos.end) or '-' in
                                str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                            rn_tx_hgvs_not_delins = self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins)

                        elif '+' in str(rn_tx_hgvs_not_delins.posedit.pos.end):
                            rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_end_base_to_next_nonoffset(
                                rn_tx_hgvs_not_delins, saved_hgvs_coding, back=False)

                        elif '+' in str(rn_tx_hgvs_not_delins.posedit.pos.start):
                            rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_start_base_to_previous_nonoffset(
                                rn_tx_hgvs_not_delins, saved_hgvs_coding)

                        # Check for -ve base and adjust
                        elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.end) and '-' in \
                                str(rn_tx_hgvs_not_delins.posedit.pos.start):
                            rn_tx_hgvs_not_delins = self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins)

                        elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.end):
                            rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_end_base_to_next_nonoffset(
                                rn_tx_hgvs_not_delins, saved_hgvs_coding)

                        elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.start):
                            rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_start_base_to_previous_nonoffset(
                                rn_tx_hgvs_not_delins, saved_hgvs_coding, with_base_subtract=True)

                        # Logic
                        if len(hgvs_not_delins.posedit.edit.ref) < len(rn_tx_hgvs_not_delins.posedit.edit.ref):
                            gap_length = len(rn_tx_hgvs_not_delins.posedit.edit.ref) - len(
                                hgvs_not_delins.posedit.edit.ref)
                            self.disparity_deletion_in = ['chromosome', gap_length]
                        elif len(hgvs_not_delins.posedit.edit.ref) > len(rn_tx_hgvs_not_delins.posedit.edit.ref):
                            gap_length = len(hgvs_not_delins.posedit.edit.ref) - len(
                                rn_tx_hgvs_not_delins.posedit.edit.ref)
                            self.disparity_deletion_in = ['transcript', gap_length]
                        else:
                            # store stash_hgvs_not_delins for restorstion after error below
                            restore_stash_hgvs_not_delins = copy.copy(stash_hgvs_not_delins)
                            try:
                                hgvs_stash_t = self.validator.vm.g_to_t(stash_hgvs_not_delins, saved_hgvs_coding.ac,
                                                                        alt_aln_method=self.validator.alt_aln_method)
                            except vvhgvs.exceptions.HGVSError as e:
                                if 'bounds' in str(e):
                                    stash_hgvs_not_delins = copy.copy(stored_hgvs_not_delins)
                                    try:
                                        hgvs_stash_t = self.validator.vm.g_to_t(stash_hgvs_not_delins,
                                                                                saved_hgvs_coding.ac,
                                                                                alt_aln_method=self.validator.alt_aln_method)
                                    except vvhgvs.exceptions.HGVSError:
                                        hgvs_stash_t = saved_hgvs_coding

                            if len(stash_hgvs_not_delins.posedit.edit.ref) > len(hgvs_stash_t.posedit.edit.ref):
                                try:
                                    self.variant.hn.normalize(hgvs_stash_t)
                                except Exception as e:
                                    logger.debug("Except passed, %s", e)
                                else:
                                    gap_length = len(stash_hgvs_not_delins.posedit.edit.ref) - len(
                                        hgvs_stash_t.posedit.edit.ref)
                                    self.disparity_deletion_in = ['transcript', gap_length]
                                    try:
                                        self.tx_hgvs_not_delins = self.validator.vm.c_to_n(hgvs_stash_t)
                                    except:
                                        self.tx_hgvs_not_delins = hgvs_stash_t
                                    hgvs_not_delins = stash_hgvs_not_delins
                            elif hgvs_stash_t.posedit.pos.start.offset != 0 or hgvs_stash_t.posedit.pos.end.offset != 0:
                                self.disparity_deletion_in = ['transcript', 'Requires Analysis']
                                try:
                                    self.tx_hgvs_not_delins = self.validator.vm.c_to_n(hgvs_stash_t)
                                except:
                                    self.tx_hgvs_not_delins = hgvs_stash_t
                                hgvs_not_delins = stash_hgvs_not_delins
                                self.hgvs_genomic_5pr = stash_hgvs_not_delins
                            else:
                                try:
                                    var_a = self.variant.hn.normalize(hgvs_stash_t)
                                    var_b = self.variant.hn.normalize(original_var)
                                except vvhgvs.exceptions.HGVSError:
                                    pass
                                else:
                                    if var_a.posedit.edit.type != var_b.posedit.edit.type:
                                        gap_warnings = self.make_gap_warnings(self.tx_hgvs_not_delins.ac,
                                                                              self.hgvs_genomic_5pr.ac,
                                                                              self.variant.primary_assembly)

                                        if gap_warnings["gapped_alignment_warning"] is not None \
                                                and gap_warnings["auto_info"] is not None:
                                            gapped_alignment_warning = gap_warnings["gapped_alignment_warning"]
                                            if ("fewer" in gap_warnings["auto_info"] or
                                                    "extra" in gap_warnings["auto_info"]):
                                                self.auto_info = self.auto_info + gap_warnings["auto_info"]

                            # Restore stash_hgvs_not_delins
                            stash_hgvs_not_delins = restore_stash_hgvs_not_delins

                    # Final sanity checks
                    try:
                        self.validator.vm.g_to_t(hgvs_not_delins, self.tx_hgvs_not_delins.ac,
                                                 alt_aln_method=self.validator.alt_aln_method)
                    except Exception as e:
                        if str(e) == 'start or end or both are beyond the bounds of transcript record':
                            hgvs_not_delins = saved_hgvs_coding
                            self.disparity_deletion_in = ['false', 'false']
                        logger.info(str(e))
                    try:
                        self.variant.hn.normalize(self.tx_hgvs_not_delins)
                    except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                        error = str(e)
                        if 'Normalization of intronic variants is not supported' in error or \
                                'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                            if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                                hgvs_not_delins = saved_hgvs_coding
                                self.disparity_deletion_in = ['false', 'false']
                            elif 'Normalization of intronic variants is not supported' in error:
                                # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                                self.disparity_deletion_in = ['transcript', 'Requires Analysis']
                        logger.info(error)

                    # Pre-processing of self.tx_hgvs_not_delins
                    try:
                        if self.tx_hgvs_not_delins.posedit.edit.alt is None:
                            self.tx_hgvs_not_delins.posedit.edit.alt = ''
                    except Exception as e:
                        if str(e) == "'Dup' object has no attribute 'alt'":
                            tx_hgvs_not_delins_delins_from_dup = fn.hgvs_dup2indel(self.tx_hgvs_not_delins)
                            self.tx_hgvs_not_delins = \
                                self.validator.hp.parse_hgvs_variant(tx_hgvs_not_delins_delins_from_dup)

                    # GAP IN THE TRANSCRIPT DISPARITY DETECTED
                    if self.disparity_deletion_in[0] == 'transcript':

                        # Check for issue https://github.com/openvar/variantValidator/issues/385 where the gap is
                        # being identified but oddly the vm is not compensating, likely due to odd sequence
                        try:
                            if len(self.tx_hgvs_not_delins.posedit.edit.ref) > \
                                    len(self.tx_hgvs_not_delins.posedit.edit.alt):
                                gen_len_difference = len(hgvs_not_delins.posedit.edit.ref) - \
                                                     len(hgvs_not_delins.posedit.edit.alt)
                                tx_len_difference = len(self.tx_hgvs_not_delins.posedit.edit.ref) - \
                                                    len(self.tx_hgvs_not_delins.posedit.edit.alt)
                            else:
                                gen_len_difference = len(hgvs_not_delins.posedit.edit.alt) - \
                                                     len(hgvs_not_delins.posedit.edit.ref)
                                tx_len_difference = len(self.tx_hgvs_not_delins.posedit.edit.alt) - \
                                                    len(self.tx_hgvs_not_delins.posedit.edit.ref)

                            # The logic here. Since there is a gap in the transcript,
                            # the actual length should be == gen_len_difference - 1 not == gen_len_difference
                            if tx_len_difference - self.disparity_deletion_in[1] == gen_len_difference:
                                # So here we know we need to knock off disparity_deletion_in[1] bases
                                if len(hgvs_not_delins.posedit.edit.alt) == len(self.tx_hgvs_not_delins.posedit.edit.alt):
                                    if self.orientation == 1:
                                        self.tx_hgvs_not_delins.posedit.edit.ref = hgvs_not_delins.posedit.ref
                                    else:
                                        replace_ref_bases = self.validator.revcomp(hgvs_not_delins.posedit.edit.ref)
                                        self.tx_hgvs_not_delins.posedit.edit.ref = replace_ref_bases
                                    self.tx_hgvs_not_delins.posedit.pos.end.offset = self.disparity_deletion_in[1]

                        except TypeError:
                            pass
                        except AttributeError:
                            pass

                        gap_warnings = self.make_gap_warnings(self.tx_hgvs_not_delins.ac,
                                                              self.hgvs_genomic_5pr.ac,
                                                              self.variant.primary_assembly)

                        if gap_warnings["gapped_alignment_warning"] is not None \
                                and gap_warnings["auto_info"] is not None:
                            gapped_alignment_warning = gap_warnings["gapped_alignment_warning"]
                            if ("fewer" in gap_warnings["auto_info"] or
                                    "extra" in gap_warnings["auto_info"]):
                                self.auto_info = self.auto_info + gap_warnings["auto_info"]

                        # ANY VARIANT WHOLLY WITHIN THE GAP
                        hgvs_refreshed_variant = self.transcript_disparity(reverse_normalized_hgvs_genomic,
                                                                           stored_hgvs_not_delins,
                                                                           self.variant.hgvs_genomic, 1)
                    # GAP IN THE CHROMOSOME
                    elif self.disparity_deletion_in[0] == 'chromosome':
                        # Set warning variables
                        gap_warnings = self.make_gap_warnings(self.tx_hgvs_not_delins.ac,
                                                              self.hgvs_genomic_5pr.ac,
                                                              self.variant.primary_assembly)

                        if gap_warnings["gapped_alignment_warning"] is not None \
                                and gap_warnings["auto_info"] is not None:
                            gapped_alignment_warning = gap_warnings["gapped_alignment_warning"]
                            if ("fewer" in gap_warnings["auto_info"] or
                                    "extra" in gap_warnings["auto_info"]):
                                self.auto_info = self.auto_info + gap_warnings["auto_info"]
                        hgvs_refreshed_variant = self.tx_hgvs_not_delins

                    else:
                        # Have we already had a hard push?
                        # if needs_a_push is not False:
                        if merged_variant is not False:
                            hgvs_refreshed_variant = saved_hgvs_coding
                        else:
                            try:
                                # Try the push to see if a gap is identified
                                hgvs_stash = copy.deepcopy(stash_hgvs_not_delins)
                                stash_ac = hgvs_stash.ac

                                # Make a hard left and hard right not delins g.
                                stash_dict_right = hgvs_utils.hard_right_hgvs2vcf(hgvs_stash,
                                                                                  self.variant.primary_assembly,
                                                                                  self.variant.hn,
                                                                                  self.variant.reverse_normalizer,
                                                                                  self.validator.sf,
                                                                                  saved_hgvs_coding.ac,
                                                                                  self.validator.hdp,
                                                                                  self.validator.alt_aln_method,
                                                                                  self.validator.hp,
                                                                                  self.validator.vm,
                                                                                  self.validator.merge_hgvs_3pr)
                                stash_pos_right = int(stash_dict_right['pos'])
                                stash_ref_right = stash_dict_right['ref']
                                stash_alt_right = stash_dict_right['alt']
                                stash_end_right = str(stash_pos_right + len(stash_ref_right) - 1)
                                stash_hgvs_not_delins_right = self.validator.hp.parse_hgvs_variant(stash_ac + ':' +
                                    hgvs_stash.type + '.' + str(stash_pos_right) + '_' + stash_end_right + 'del' +
                                    stash_ref_right + 'ins' + stash_alt_right)

                                stash_dict_left = hgvs_utils.hard_left_hgvs2vcf(hgvs_stash,
                                                                                self.variant.primary_assembly,
                                                                                self.variant.hn,
                                                                                self.variant.reverse_normalizer,
                                                                                self.validator.sf,
                                                                                saved_hgvs_coding.ac,
                                                                                self.validator.hdp,
                                                                                self.validator.alt_aln_method,
                                                                                self.validator.hp,
                                                                                self.validator.vm,
                                                                                self.validator.merge_hgvs_3pr)
                                stash_pos_left = int(stash_dict_left['pos'])
                                stash_ref_left = stash_dict_left['ref']
                                stash_alt_left = stash_dict_left['alt']
                                stash_end_left = str(stash_pos_left + len(stash_ref_left) - 1)
                                stash_hgvs_not_delins_left = self.validator.hp.parse_hgvs_variant(
                                    stash_ac + ':' + hgvs_stash.type + '.' + str(
                                        stash_pos_left) + '_' + stash_end_left + 'del' + stash_ref_left + 'ins' +
                                    stash_alt_left)
                            except vvhgvs.exceptions.HGVSDataNotAvailableError:
                                continue

                            # Map in-situ to the transcript left and right
                            try:
                                tx_hard_right = self.validator.vm.g_to_t(stash_hgvs_not_delins_right,
                                                                         saved_hgvs_coding.ac,
                                                                         alt_aln_method=self.validator.alt_aln_method)
                            except Exception:
                                tx_hard_right = saved_hgvs_coding
                            else:
                                normalize_stash_right = self.variant.hn.normalize(stash_hgvs_not_delins_right)
                                if str(normalize_stash_right.posedit) == str(stash_hgvs_not_delins.posedit):
                                    tx_hard_right = saved_hgvs_coding
                            try:
                                tx_hard_left = self.validator.vm.g_to_t(stash_hgvs_not_delins_left,
                                                                        saved_hgvs_coding.ac,
                                                                        alt_aln_method=self.validator.alt_aln_method)
                            except Exception:
                                tx_hard_left = saved_hgvs_coding
                            else:
                                normalize_stash_left = self.variant.hn.normalize(stash_hgvs_not_delins_left)
                                if str(normalize_stash_left.posedit) == str(stash_hgvs_not_delins.posedit):
                                    tx_hard_left = saved_hgvs_coding

                            try:
                                # The Logic - Currently limited to genome gaps
                                if len(stash_hgvs_not_delins_right.posedit.edit.ref) < len(
                                        tx_hard_right.posedit.edit.ref) or \
                                        len(stash_hgvs_not_delins_right.posedit.edit.ref) > \
                                        len(tx_hard_right.posedit.edit.ref):
                                    tx_hard_right = self.variant.hn.normalize(tx_hard_right)

                                    gap_warnings = self.make_gap_warnings(self.tx_hgvs_not_delins.ac,
                                                                          self.hgvs_genomic_5pr.ac,
                                                                          self.variant.primary_assembly)

                                    if gap_warnings["gapped_alignment_warning"] is not None \
                                            and gap_warnings["auto_info"] is not None:
                                        gapped_alignment_warning = gap_warnings["gapped_alignment_warning"]
                                        if ("fewer" in gap_warnings["auto_info"] or
                                                "extra" in gap_warnings["auto_info"]):
                                            self.auto_info = self.auto_info + gap_warnings["auto_info"]
                                    hgvs_refreshed_variant = tx_hard_right

                                elif len(stash_hgvs_not_delins_left.posedit.edit.ref) < \
                                        len(tx_hard_left.posedit.edit.ref) or \
                                        len(stash_hgvs_not_delins_left.posedit.edit.ref) > \
                                        len(tx_hard_left.posedit.edit.ref):
                                    tx_hard_left = self.variant.hn.normalize(tx_hard_left)
                                    gap_warnings = self.make_gap_warnings(self.tx_hgvs_not_delins.ac,
                                                                          self.hgvs_genomic_5pr.ac,
                                                                          self.variant.primary_assembly)

                                    if gap_warnings["gapped_alignment_warning"] is not None \
                                            and gap_warnings["auto_info"] is not None:
                                        gapped_alignment_warning = gap_warnings["gapped_alignment_warning"]
                                        if ("fewer" in gap_warnings["auto_info"] or
                                                "extra" in gap_warnings["auto_info"]):
                                            self.auto_info = self.auto_info + gap_warnings["auto_info"]
                                    hgvs_refreshed_variant = tx_hard_left

                                else:
                                    # Keep the same by re-setting rel_var
                                    hgvs_refreshed_variant = saved_hgvs_coding

                            except TypeError:
                                # e.g. chr1:156561557G>GGGGTC (investigate at a later date)
                                hgvs_refreshed_variant = saved_hgvs_coding

                            except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                                # e.g. NG_005895.1:g.3684_44407del
                                hgvs_refreshed_variant = saved_hgvs_coding

                    # Edit the output
                    hgvs_refreshed_variant = self.edit_output(hgvs_refreshed_variant, saved_hgvs_coding)

                    # Send to empty nw_rel_var
                    if hgvs_refreshed_variant.posedit.edit.type == "delins" and \
                            hgvs_refreshed_variant.posedit.edit.alt == "":
                        correct = str(hgvs_refreshed_variant).replace("ins", "")
                        hgvs_refreshed_variant = self.validator.hp.parse(correct)
                    nw_rel_var.append(hgvs_refreshed_variant)

                # Otherwise these variants need to be set
                else:
                    corrective_action_taken = ''
                    gapped_alignment_warning = ''
                    # Send to empty nw_rel_var
                    nw_rel_var.append(saved_hgvs_coding)

            # Otherwise these variants need to be set
            else:
                corrective_action_taken = ''
                gapped_alignment_warning = ''
                # Send to empty nw_rel_var
                nw_rel_var.append(saved_hgvs_coding)

        data = {
            'gapped_alignment_warning': gapped_alignment_warning,
            'corrective_action_taken': corrective_action_taken,
            'auto_info': self.auto_info,
            'disparity_deletion_in': self.disparity_deletion_in,
            'gapped_transcripts': self.gapped_transcripts
        }

        return data, nw_rel_var

    def g_to_t_compensation(self, ori, hgvs_coding, rec_var):
        self.orientation = int(ori[0]['alt_strand'])
        self.hgvs_genomic_possibilities = []
        hgvs_genomic = self.validator.myevm_t_to_g(hgvs_coding, self.variant.no_norm_evm, self.variant.primary_assembly,
                                                   self.variant.hn)

        logger.debug('g_to_t gap code 1 active')
        rn_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic)
        self.hgvs_genomic_possibilities.append([rn_hgvs_genomic, ['false', 'false']])

        try:
            if self.orientation != -1:
                chromosome_normalized_hgvs_coding = self.variant.reverse_normalizer.normalize(hgvs_coding)
            else:
                chromosome_normalized_hgvs_coding = self.variant.hn.normalize(hgvs_coding)
        except vvhgvs.exceptions.HGVSUnsupportedOperationError:
            chromosome_normalized_hgvs_coding = hgvs_coding

        most_3pr_hgvs_genomic = self.validator.myvm_t_to_g(chromosome_normalized_hgvs_coding, hgvs_genomic.ac,
                                                           self.variant.no_norm_evm, self.variant.hn)
        self.hgvs_genomic_possibilities.append([most_3pr_hgvs_genomic, ['false', 'false']])

        # Push from side to side to try pick up odd placements
        # MAKE A NO NORM HGVS2VCF
        # First to the right
        hgvs_stash = copy.deepcopy(hgvs_coding)
        stash_tx_right = ''
        stash_tx_left = ''
        try:
            hgvs_stash = self.variant.no_norm_evm.c_to_n(hgvs_stash)
        except Exception as e:
            logger.debug("Except passed, %s", e)
        try:
            stash_ac = hgvs_stash.ac
            stash_dict = hgvs_utils.hard_right_hgvs2vcf(hgvs_stash,
                                                        self.variant.primary_assembly,
                                                        self.variant.hn,
                                                        self.variant.reverse_normalizer,
                                                        self.validator.sf,
                                                        hgvs_coding.ac,
                                                        self.validator.hdp,
                                                        self.validator.alt_aln_method,
                                                        self.validator.hp,
                                                        self.validator.vm,
                                                        self.validator.merge_hgvs_3pr,
                                                        genomic_ac=hgvs_genomic.ac)

            stash_pos = int(stash_dict['pos'])
            stash_ref = stash_dict['ref']
            stash_alt = stash_dict['alt']
            # Generate an end position
            stash_end = str(stash_pos + len(stash_ref) - 1)
            # make a not real deletion insertion
            stash_hgvs_not_delins = self.validator.hp.parse_hgvs_variant(
                stash_ac + ':' + hgvs_stash.type + '.' + str(
                    stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)
            try:
                stash_hgvs_not_delins = self.variant.no_norm_evm.n_to_c(stash_hgvs_not_delins)
            except Exception as e:
                logger.debug("Except passed, %s", e)

            test_stash_tx_right = copy.deepcopy(stash_hgvs_not_delins)
            stash_genomic = self.validator.myvm_t_to_g(test_stash_tx_right, hgvs_genomic.ac, self.variant.no_norm_evm,
                                                       self.variant.hn)
            if len(test_stash_tx_right.posedit.edit.ref) == (
                    (stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                stash_tx_right = test_stash_tx_right
                if hasattr(test_stash_tx_right.posedit.edit,
                           'alt') and test_stash_tx_right.posedit.edit.alt is not None:
                    alt = test_stash_tx_right.posedit.edit.alt
                else:
                    alt = ''
                if hasattr(stash_genomic.posedit.edit,
                           'alt') and stash_genomic.posedit.edit.alt is not None:
                    g_alt = stash_genomic.posedit.edit.alt
                else:
                    g_alt = ''
                if (len(alt) - (
                        test_stash_tx_right.posedit.pos.end.base - test_stash_tx_right.posedit.pos.start.base) + 1) != (
                        len(g_alt) - (
                        stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])
                else:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
            elif test_stash_tx_right.posedit.edit.type == 'identity':
                reform_ident = str(test_stash_tx_right).split(':')[0]
                reform_ident = reform_ident + ':c.' + str(test_stash_tx_right.posedit.pos) + 'del' + str(
                    test_stash_tx_right.posedit.edit.ref)  # + 'ins' + str(test_stash_tx_right.posedit.edit.alt)
                hgvs_reform_ident = self.validator.hp.parse_hgvs_variant(reform_ident)
                try:
                    self.variant.hn.normalize(hgvs_reform_ident)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if re.search('spanning the exon-intron boundary', error):
                        stash_tx_right = test_stash_tx_right
                        self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_right = test_stash_tx_right
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])
            else:
                try:
                    self.variant.hn.normalize(test_stash_tx_right)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_right = test_stash_tx_right
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])

            # IDENTIFY GAP AND HARD SET
            if stash_dict['needs_a_push'] is True or stash_dict['identifying_g_variant'] is not False:

                # Look for merged variant from hard push
                if stash_dict['merged_variant'] is not False:
                    merged_variant = stash_dict['merged_variant']
                    try:
                        merged_variant = self.validator.vm.n_to_c(merged_variant)
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass
                    identifying_variant = stash_dict['identifying_variant']
                    try:
                        identifying_variant = self.validator.vm.n_to_c(identifying_variant)
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass
                    stash_genomic = self.validator.myvm_t_to_g(identifying_variant, stash_genomic.ac,
                                                               self.variant.no_norm_evm,
                                                               self.variant.hn)
                    stash_hgvs_not_delins = merged_variant
                    test_stash_tx_right = merged_variant

                elif stash_dict['identifying_g_variant'] is not False:
                    stash_hgvs_not_delins = stash_dict['identifying_g_variant']
                    stash_genomic = stash_dict['identifying_g_variant']

                # Look for gap info
                normalized_stash_genomic = self.variant.hn.normalize(stash_genomic)
                stash_tx_right = test_stash_tx_right
                if stash_hgvs_not_delins.posedit.edit.type == "ins":
                    len_tx = 2
                else:
                    len_tx = len(stash_hgvs_not_delins.posedit.edit.ref)
                if stash_genomic.posedit.edit.type == "ins":
                    len_gen = 2
                else:
                    len_gen = len(stash_genomic.posedit.edit.ref)
                if len_tx > len_gen:
                    gap_in = 'chromosome'
                    gap_len = len_tx - len_gen
                else:
                    gap_in = 'transcript'
                    gap_len = len_gen - len_tx

                # Set the options to a single option based on the results of pushing
                self.hgvs_genomic_possibilities = [[normalized_stash_genomic, [gap_in,
                                                                               gap_len,
                                                                               stash_hgvs_not_delins,
                                                                               stash_genomic]]]

        except vvhgvs.exceptions.HGVSError as e:
            logger.debug("Except passed, %s", e)
        # Intronic positions not supported. Will cause a Value Error
        except ValueError as e:
            logger.debug("Except passed, %s", e)

        # Then to the left
        hgvs_stash = copy.deepcopy(hgvs_coding)
        try:
            hgvs_stash = self.variant.no_norm_evm.c_to_n(hgvs_stash)
        except Exception as e:
            logger.debug("Except passed, %s", e)
        try:
            stash_ac = hgvs_stash.ac
            stash_dict = hgvs_utils.hard_left_hgvs2vcf(hgvs_stash,
                                                       self.variant.primary_assembly,
                                                       self.variant.hn,
                                                       self.variant.reverse_normalizer,
                                                       self.validator.sf,
                                                       hgvs_coding.ac,
                                                       self.validator.hdp,
                                                       self.validator.alt_aln_method,
                                                       self.validator.hp,
                                                       self.validator.vm,
                                                       self.validator.merge_hgvs_3pr,
                                                       genomic_ac=hgvs_genomic.ac)

            stash_pos = int(stash_dict['pos'])
            stash_ref = stash_dict['ref']
            stash_alt = stash_dict['alt']
            # Generate an end position
            stash_end = str(stash_pos + len(stash_ref) - 1)
            # make a not real deletion insertion
            stash_hgvs_not_delins = self.validator.hp.parse_hgvs_variant(
                stash_ac + ':' + hgvs_stash.type + '.' + str(
                    stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)
            try:
                stash_hgvs_not_delins = self.variant.no_norm_evm.n_to_c(stash_hgvs_not_delins)
            except Exception as e:
                logger.debug("Except passed, %s", e)
                # Store a tx copy for later use
            test_stash_tx_left = copy.deepcopy(stash_hgvs_not_delins)
            stash_genomic = self.validator.myvm_t_to_g(test_stash_tx_left, hgvs_genomic.ac, self.variant.no_norm_evm,
                                                       self.variant.hn)

            if len(test_stash_tx_left.posedit.edit.ref) == ((stash_genomic.posedit.pos.end.base -
                                                             stash_genomic.posedit.pos.start.base) + 1):
                if hasattr(test_stash_tx_left.posedit.edit,
                           'alt') and test_stash_tx_left.posedit.edit.alt is not None:
                    alt = test_stash_tx_left.posedit.edit.alt
                else:
                    alt = ''
                if hasattr(stash_genomic.posedit.edit,
                           'alt') and stash_genomic.posedit.edit.alt is not None:
                    g_alt = stash_genomic.posedit.edit.alt
                else:
                    g_alt = ''

                if (len(alt) - (
                        test_stash_tx_left.posedit.pos.end.base - test_stash_tx_left.posedit.pos.start.base) + 1) != (
                        len(g_alt) - (
                        stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])
                else:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
            elif test_stash_tx_left.posedit.edit.type == 'identity':
                reform_ident = str(test_stash_tx_left).split(':')[0]
                reform_ident = reform_ident + ':c.' + str(test_stash_tx_left.posedit.pos) + 'del' + str(
                    test_stash_tx_left.posedit.edit.ref)  # + 'ins' + str(test_stash_tx_left.posedit.edit.alt)
                hgvs_reform_ident = self.validator.hp.parse_hgvs_variant(reform_ident)
                try:
                    self.variant.hn.normalize(hgvs_reform_ident)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if re.search('spanning the exon-intron boundary', error):
                        stash_tx_left = test_stash_tx_left
                        self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_left = test_stash_tx_left
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])
            else:
                try:
                    self.variant.hn.normalize(test_stash_tx_left)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_left = test_stash_tx_left
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])

            # IDENTIFY GAP AND HARD SET
            if stash_dict['needs_a_push'] is True or stash_dict['identifying_g_variant'] is not False:

                # Look for merged variant from hard push
                if stash_dict['merged_variant'] is not False:
                    merged_variant = stash_dict['merged_variant']
                    try:
                        merged_variant = self.validator.vm.n_to_c(merged_variant)
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass
                    identifying_variant = stash_dict['identifying_variant']
                    try:
                        identifying_variant = self.validator.vm.n_to_c(identifying_variant)
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass
                    stash_genomic = self.validator.myvm_t_to_g(identifying_variant, stash_genomic.ac,
                                                               self.variant.no_norm_evm,
                                                               self.variant.hn)
                    stash_hgvs_not_delins = merged_variant
                    test_stash_tx_left = merged_variant

                elif stash_dict['identifying_g_variant'] is not False:
                    stash_hgvs_not_delins = stash_dict['identifying_g_variant']
                    stash_genomic = stash_dict['identifying_g_variant']

                # Look for gap info
                normalized_stash_genomic = self.variant.hn.normalize(stash_genomic)
                stash_tx_left = test_stash_tx_left
                if stash_hgvs_not_delins.posedit.edit.type == "ins":
                    len_tx = 2
                else:
                    len_tx = len(stash_hgvs_not_delins.posedit.edit.ref)
                if stash_genomic.posedit.edit.type == "ins":
                    len_gen = 2
                else:
                    len_gen = len(stash_genomic.posedit.edit.ref)
                if len_tx > len_gen:
                    gap_in = 'chromosome'
                    gap_len = len_tx - len_gen
                else:
                    gap_in = 'transcript'
                    gap_len = len_gen - len_tx

                # Set the options to a single option based on the results of pushing
                self.hgvs_genomic_possibilities = [[normalized_stash_genomic, [gap_in,
                                                                               gap_len,
                                                                               stash_hgvs_not_delins,
                                                                               stash_genomic]]]
        except vvhgvs.exceptions.HGVSError as e:
            logger.debug("Except passed, %s", e)
        # Intronic positions not supported. Will cause a Value Error
        except ValueError as e:
            logger.debug("Except passed, %s", e)

        # direct mapping from reverse_normalized transcript insertions in the delins format
        self.rev_norm_ins(hgvs_coding, hgvs_genomic)

        # Set variables for problem specific warnings
        self.gapped_transcripts = ''
        self.auto_info = ''

        # Mark as not disparity detected
        self.disparity_deletion_in = ['false', 'false']

        # Loop through to see if a gap can be located
        # Set the variables required for corrective normalization
        possibility_counter = 0
        suppress_c_normalization = 'false'  # Applies to boundary crossing normalization

        # If hard pushing identified a variant (in theory there can be only 1) then remove all other possibilities
        hard_possibility = []
        for hard_set_check in self.hgvs_genomic_possibilities:
            if hard_set_check[1] != ['false', 'false']:
                hard_possibility.append(hard_set_check)

        # Copy a version of hgvs_genomic_possibilities
        for a_possibility in self.hgvs_genomic_possibilities:
            possibility = a_possibility[0]
            disparity_info = a_possibility[1]
            possibility_counter = possibility_counter + 1

            # Loop out stash possibilities which will not spot gaps so are empty
            if possibility == '':
                continue

            # Use VCF generation code to push hgvs_genomic as for 5 prime as possible to uncover gaps
            hgvs_genomic_variant = copy.deepcopy(possibility)
            reverse_normalized_hgvs_genomic = ''
            # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
            try:
                reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic_variant)
            except vvhgvs.exceptions.HGVSError as e:
                # Strange error caused by gap in genomic
                error = str(e)
                if 'base start position must be <= end position' in error:
                    if hgvs_genomic.posedit.edit.type == 'delins':
                        start = hgvs_genomic.posedit.pos.start.base
                        end = hgvs_genomic.posedit.pos.end.base
                        lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                        rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                        hgvs_genomic.posedit.edit.ref = lhb + rhb
                        hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                        hgvs_genomic.posedit.pos.start.base = end
                        hgvs_genomic.posedit.pos.end.base = start
                        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic)
                    if hgvs_genomic.posedit.edit.type == 'del':
                        start = hgvs_genomic.posedit.pos.start.base
                        end = hgvs_genomic.posedit.pos.end.base
                        lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                        rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                        hgvs_genomic.posedit.edit.ref = lhb + rhb
                        hgvs_genomic.posedit.edit.alt = lhb + rhb
                        hgvs_genomic.posedit.pos.start.base = end
                        hgvs_genomic.posedit.pos.end.base = start
                        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic)
                if 'insertion length must be 1' in error:
                    if hgvs_genomic.posedit.edit.type == 'ins':
                        start = hgvs_genomic.posedit.pos.start.base
                        end = hgvs_genomic.posedit.pos.end.base
                        lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                        rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start, end)
                        hgvs_genomic.posedit.edit.ref = lhb + rhb
                        hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic)

            self.hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

            # Create VCF
            vcf_dict = hgvs_utils.hgvs2vcf(reverse_normalized_hgvs_genomic, self.variant.primary_assembly,
                                           self.variant.reverse_normalizer, self.validator.sf)
            pos = vcf_dict['pos']
            ref = vcf_dict['ref']
            alt = vcf_dict['alt']

            # Generate an end position
            end = str(int(pos) + len(ref) - 1)
            pos = str(pos)

            # Store a not real deletion insertion to test for gapping
            stored_hgvs_not_delins = self.validator.hp.parse_hgvs_variant(str(self.hgvs_genomic_5pr.ac) + ':' +
                                                                          self.hgvs_genomic_5pr.type +
                                                                          '.' + pos + '_' + end + 'del' + ref +
                                                                          'ins' + alt)

            # Detect intronic variation using normalization
            intronic_variant = 'false'

            # Save a copy of current hgvs_coding
            try:
                saved_hgvs_coding = self.variant.no_norm_evm.g_to_t(stored_hgvs_not_delins,
                                                                    hgvs_coding.ac)
            except vvhgvs.exceptions.HGVSInvalidIntervalError as e:
                if str(e) == 'start or end or both are beyond the bounds of transcript record':
                    continue
                else:
                    saved_hgvs_coding = self.variant.no_norm_evm.g_to_n(reverse_normalized_hgvs_genomic,
                                                                        hgvs_coding.ac)

            # Look for normalized variant options that do not match hgvs_coding
            hgvs_seek_var = self.get_hgvs_seek_var(hgvs_genomic, hgvs_coding)

            if (
                    hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                    hgvs_coding.posedit.pos.start.base + hgvs_coding.posedit.pos.start.offset) and (
                    hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                    hgvs_coding.posedit.pos.end.base + hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
                pass
            else:
                hgvs_seek_var = saved_hgvs_coding

            try:
                self.variant.hn.normalize(hgvs_seek_var)
            except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                error = str(e)
                if 'Normalization of intronic variants is not supported' in error or \
                        'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                    if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        intronic_variant = 'hard_fail'
                    else:
                        # Double check to see whether the variant is actually intronic?
                        for exon in ori:
                            genomic_start = int(exon['alt_start_i'])
                            genomic_end = int(exon['alt_end_i'])
                            if (genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                                    genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                                intronic_variant = 'false'
                                break
                            else:
                                intronic_variant = 'true'

            if intronic_variant != 'hard_fail':
                if re.search(r'\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search(r'\d+-', str(
                        hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+\+', str(
                        hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+-', str(hgvs_seek_var.posedit.pos)):
                    # Double check to see whether the variant is actually intronic?
                    for exon in ori:
                        genomic_start = int(exon['alt_start_i'])
                        genomic_end = int(exon['alt_end_i'])
                        if (genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                                genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                            intronic_variant = 'false'
                            break
                        else:
                            intronic_variant = 'true'

            if re.search(r'\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search(r'\d+-', str(
                    hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+\+', str(
                    hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+-', str(hgvs_seek_var.posedit.pos)):
                # Double check to see whether the variant is actually intronic?
                for exon in ori:
                    genomic_start = int(exon['alt_start_i'])
                    genomic_end = int(exon['alt_end_i'])
                    if (genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                            genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                        intronic_variant = 'false'
                        break
                    else:
                        intronic_variant = 'true'

            if intronic_variant != 'true':
                # Flag RefSeqGene for amendment
                # amend_RefSeqGene = 'false'
                # Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping
                # transcript lengths
                hgvs_not_delins = ''
                if stored_hgvs_not_delins != '':
                    # Refresh hgvs_not_delins from stored_hgvs_not_delins
                    hgvs_not_delins = self.dup_ins_5prime_shift(stored_hgvs_not_delins, saved_hgvs_coding)

                    try:
                        self.tx_hgvs_not_delins = self.variant.no_norm_evm.g_to_n(hgvs_not_delins,
                                                                                  saved_hgvs_coding.ac)
                    except vvhgvs.exceptions.HGVSInvalidIntervalError:
                        self.tx_hgvs_not_delins = self.variant.no_norm_evm.g_to_n(reverse_normalized_hgvs_genomic,
                                                                                  saved_hgvs_coding.ac)
                    # Create normalized version of tx_hgvs_not_delins
                    rn_tx_hgvs_not_delins = copy.deepcopy(self.tx_hgvs_not_delins)

                    # Check for +1 base and adjust
                    if '+' in str(rn_tx_hgvs_not_delins.posedit.pos.end) and '+' in str(
                            rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins = self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins)

                    elif '+' in str(rn_tx_hgvs_not_delins.posedit.pos.end):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_end_base_to_next_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding, back=False)

                    elif '+' in str(rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_start_base_to_previous_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding)

                    # Check for -ve base and adjust
                    elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.end) and '-' in str(
                            rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins = self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins)

                    elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.end):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_end_base_to_next_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding)

                    elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_start_base_to_previous_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding, with_base_subtract=True)

                    # Logic
                    hgvs_not_delins = self.logic_check(hgvs_not_delins, rn_tx_hgvs_not_delins, hgvs_coding)

                # 'At hgvs_genomic'
                # Final sanity checks
                try:
                    self.validator.vm.g_to_t(hgvs_not_delins, self.tx_hgvs_not_delins.ac,
                                             alt_aln_method=self.validator.alt_aln_method)
                except Exception as e:
                    if str(e) == 'start or end or both are beyond the bounds of transcript record':
                        continue
                try:
                    self.variant.hn.normalize(self.tx_hgvs_not_delins)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                    error = str(e)
                    if 'Normalization of intronic variants is not supported' in error or \
                            'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                            continue
                        elif 'Normalization of intronic variants is not supported' in error:
                            # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                            self.disparity_deletion_in = ['transcript', 'Requires Analysis']

                # Recreate hgvs_genomic
                if self.disparity_deletion_in[0] == 'transcript':
                    hgvs_genomic = hgvs_not_delins

                # Find oddly placed gaps where the tx variant is encompassed in the gap
                if self.disparity_deletion_in[0] == 'false' and (possibility_counter == 3 or possibility_counter == 4):
                    rg = self.variant.reverse_normalizer.normalize(hgvs_not_delins)
                    rtx = self.validator.vm.g_to_t(rg, self.tx_hgvs_not_delins.ac,
                                                   alt_aln_method=self.validator.alt_aln_method)
                    fg = self.variant.hn.normalize(hgvs_not_delins)
                    ftx = self.validator.vm.g_to_t(fg, self.tx_hgvs_not_delins.ac,
                                                   alt_aln_method=self.validator.alt_aln_method)
                    if (rtx.posedit.pos.start.offset == 0 and rtx.posedit.pos.end.offset == 0) and (
                            ftx.posedit.pos.start.offset != 0 and ftx.posedit.pos.end.offset != 0):
                        exons = self.validator.hdp.get_tx_exons(ftx.ac, hgvs_not_delins.ac,
                                                                self.validator.alt_aln_method)
                        exonic = False
                        for ex_test in exons:
                            if ftx.posedit.pos.start.base in range(ex_test[6], ex_test[
                                    7]) and ftx.posedit.pos.end.base in range(ex_test[6], ex_test[7]):
                                exonic = True
                        if exonic is True:
                            hgvs_not_delins = fg
                            hgvs_genomic = fg
                            self.hgvs_genomic_5pr = fg
                            try:
                                self.tx_hgvs_not_delins = self.validator.vm.c_to_n(ftx)
                            except Exception:
                                self.tx_hgvs_not_delins = ftx
                            self.disparity_deletion_in = ['transcript', 'Requires Analysis']

                # Pre-processing of self.tx_hgvs_not_delins
                try:
                    if self.tx_hgvs_not_delins.posedit.edit.alt is None:
                        self.tx_hgvs_not_delins.posedit.edit.alt = ''
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        tx_hgvs_not_delins_delins_from_dup = fn.hgvs_dup2indel(self.tx_hgvs_not_delins)
                        self.tx_hgvs_not_delins = self.validator.hp.parse_hgvs_variant(
                            tx_hgvs_not_delins_delins_from_dup)

                # Has a hard set variant been identified from pushes?
                hard_set_outputs = False
                if disparity_info != ['false', 'false'] and len(disparity_info) == 4:
                    self.tx_hgvs_not_delins = disparity_info[2]
                    self.disparity_deletion_in = [disparity_info[0], disparity_info[1]]
                    hgvs_refreshed_variant = hgvs_coding
                    hgvs_genomic = possibility
                    suppress_c_normalization = 'true'
                    hard_set_outputs = True

                # GAP IN THE TRANSCRIPT DISPARITY DETECTED
                elif self.disparity_deletion_in[0] == 'transcript':
                    # Suppress intron boundary crossing due to non-intron intron based c. seq annotations
                    suppress_c_normalization = 'true'
                    # amend_RefSeqGene = 'true'
                    # ANY VARIANT WHOLLY WITHIN THE GAP
                    hgvs_refreshed_variant = self.transcript_disparity(reverse_normalized_hgvs_genomic,
                                                                       stored_hgvs_not_delins, hgvs_genomic, 2)

                # GAP IN THE CHROMOSOME
                elif self.disparity_deletion_in[0] == 'chromosome':
                    suppress_c_normalization = 'true'
                    if possibility_counter == 3:
                        hgvs_refreshed_variant = stash_tx_right
                    elif possibility_counter == 4:
                        hgvs_refreshed_variant = stash_tx_left
                    else:
                        hgvs_refreshed_variant = chromosome_normalized_hgvs_coding

                else:
                    # Keep the same by re-setting rel_var
                    hgvs_refreshed_variant = hgvs_coding
                # amend_RefSeqGene = 'false'

                # Edit the output
                if 'NM_' in str(hgvs_refreshed_variant.ac) and not 'c' in str(hgvs_refreshed_variant.type):
                    hgvs_refreshed_variant = self.variant.no_norm_evm.n_to_c(hgvs_refreshed_variant)

                try:
                    self.variant.hn.normalize(hgvs_refreshed_variant)
                except Exception as e:
                    error = str(e)

                    # Ensure the final variant is not intronic nor does it cross exon boundaries
                    if 'Normalization of intronic variants is not supported' in error or \
                            'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        hgvs_refreshed_variant = saved_hgvs_coding
                    else:
                        logger.info(error)
                        continue

                # Quick check to make sure the coding variant has not changed
                try:
                    to_test = self.variant.hn.normalize(hgvs_refreshed_variant)
                except:
                    to_test = hgvs_refreshed_variant
                if str(to_test.posedit.edit) != str(hgvs_coding.posedit.edit):
                    # Try the next available genomic option
                    if hgvs_coding.posedit.edit.type == 'identity' and to_test.posedit.edit.type == 'identity':
                        hgvs_coding = to_test
                    elif hard_set_outputs is True:
                        hgvs_coding = to_test
                    else:
                        continue

                if hard_set_outputs is not True:
                    # Update hgvs_genomic
                    hgvs_genomic = self.validator.myvm_t_to_g(hgvs_refreshed_variant, hgvs_genomic.ac,
                                                              self.variant.no_norm_evm, self.variant.hn)
                    if hgvs_genomic.posedit.edit.type == 'identity':
                        re_c = self.validator.vm.g_to_t(hgvs_genomic, hgvs_refreshed_variant.ac,
                                                        alt_aln_method=self.validator.alt_aln_method)
                        if (self.variant.hn.normalize(re_c)) != (self.variant.hn.normalize(hgvs_refreshed_variant)):
                            shuffle_left_g = copy.copy(hgvs_genomic)
                            shuffle_left_g.posedit.edit.ref = ''
                            shuffle_left_g.posedit.edit.alt = ''
                            shuffle_left_g.posedit.pos.start.base = shuffle_left_g.posedit.pos.start.base - 1
                            shuffle_left_g.posedit.pos.end.base = shuffle_left_g.posedit.pos.end.base - 1
                            shuffle_left_g = self.variant.reverse_normalizer.normalize(shuffle_left_g)
                            re_c = self.validator.vm.g_to_t(shuffle_left_g, hgvs_refreshed_variant.ac,
                                                            alt_aln_method=self.validator.alt_aln_method)
                            if (self.variant.hn.normalize(re_c)) != (self.variant.hn.normalize(hgvs_refreshed_variant)):
                                hgvs_genomic = shuffle_left_g

            # Break if gap has been detected
            if self.disparity_deletion_in[0] != 'false':
                break

        # Normailse hgvs_genomic
        try:
            hgvs_genomic = self.variant.hn.normalize(hgvs_genomic)
        except vvhgvs.exceptions.HGVSError as e:
            # Strange error caused by gap in genomic

            if 'base start position must be <= end position' in error and self.disparity_deletion_in[0] == 'chromosome':
                if hgvs_genomic.posedit.edit.type == 'delins':
                    start = hgvs_genomic.posedit.pos.start.base
                    end = hgvs_genomic.posedit.pos.end.base
                    lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                    rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                    hgvs_genomic.posedit.edit.ref = lhb + rhb
                    hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                    hgvs_genomic.posedit.pos.start.base = end
                    hgvs_genomic.posedit.pos.end.base = start
                    hgvs_genomic = self.variant.hn.normalize(hgvs_genomic)
                if hgvs_genomic.posedit.edit.type == 'del':
                    start = hgvs_genomic.posedit.pos.start.base
                    end = hgvs_genomic.posedit.pos.end.base
                    lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                    rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                    hgvs_genomic.posedit.edit.ref = lhb + rhb
                    hgvs_genomic.posedit.edit.alt = lhb + rhb
                    hgvs_genomic.posedit.pos.start.base = end
                    hgvs_genomic.posedit.pos.end.base = start
                    hgvs_genomic = self.variant.hn.normalize(hgvs_genomic)

        return hgvs_genomic, suppress_c_normalization, hgvs_coding

    def g_to_t_gapped_mapping_stage2(self, ori, hgvs_coding, hgvs_genomic):
        logger.debug('g_to_t gap code 2 active')
        hgvs_genomic_variant = hgvs_genomic
        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic_variant)
        self.hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)
        vcf_dict = hgvs_utils.hgvs2vcf(reverse_normalized_hgvs_genomic, self.variant.primary_assembly,
                                       self.variant.reverse_normalizer, self.validator.sf)
        pos = vcf_dict['pos']
        ref = vcf_dict['ref']
        alt = vcf_dict['alt']

        # DO NOT DELETE
        # Generate an end position
        end = str(int(pos) + len(ref) - 1)
        pos = str(pos)
        stored_hgvs_not_delins = self.validator.hp.parse_hgvs_variant(str(
            self.hgvs_genomic_5pr.ac) + ':' + self.hgvs_genomic_5pr.type + '.' + pos + '_' + end + 'del' + ref +
                                                                      'ins' + alt)
        self.orientation = int(ori[0]['alt_strand'])
        saved_hgvs_coding = copy.deepcopy(hgvs_coding)

        # is it in an exon?
        is_it_in_an_exon = 'no'
        for exon in ori:
            genomic_start = int(exon['alt_start_i'])
            genomic_end = int(exon['alt_end_i'])
            # Take from stored copy
            if (genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                    genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                is_it_in_an_exon = 'yes'
        if is_it_in_an_exon == 'yes':
            # map form reverse normalized g. to c.
            # Attempt to find gaps in reference sequence by catching disparity in genome length and
            # overlapping transcript lengths
            self.disparity_deletion_in = ['false', 'false']
            hgvs_not_delins = ''
            hard_fail = 'false'
            if stored_hgvs_not_delins != '':
                # Refresh hgvs_not_delins from stored_hgvs_not_delins
                hgvs_not_delins = self.dup_ins_5prime_shift(stored_hgvs_not_delins, saved_hgvs_coding)
                try:
                    self.tx_hgvs_not_delins = self.variant.no_norm_evm.g_to_n(hgvs_not_delins,
                                                                              saved_hgvs_coding.ac)
                except Exception as e:
                    if str(e) == 'start or end or both are beyond the bounds of transcript record':
                        self.tx_hgvs_not_delins = hgvs_coding
                        hard_fail = 'true'

                # Create normalized version of self.tx_hgvs_not_delins
                rn_tx_hgvs_not_delins = copy.deepcopy(self.tx_hgvs_not_delins)
                # Check for +ve base and adjust
                if '+' in str(rn_tx_hgvs_not_delins.posedit.pos.end) and '+' in \
                        str(rn_tx_hgvs_not_delins.posedit.pos.start):
                    rn_tx_hgvs_not_delins = self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins)

                elif '+' in str(rn_tx_hgvs_not_delins.posedit.pos.end):
                    rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_end_base_to_next_nonoffset(
                        rn_tx_hgvs_not_delins, saved_hgvs_coding, back=False)

                elif '+' in str(rn_tx_hgvs_not_delins.posedit.pos.start):
                    rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_start_base_to_previous_nonoffset(
                        rn_tx_hgvs_not_delins, saved_hgvs_coding)

                # Check for -ve base and adjust
                elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.end) and '-' in str(
                        rn_tx_hgvs_not_delins.posedit.pos.start):
                    rn_tx_hgvs_not_delins = self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins)

                elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.end):
                    rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_end_base_to_next_nonoffset(
                        rn_tx_hgvs_not_delins, saved_hgvs_coding)

                elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.start):
                    rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_start_base_to_previous_nonoffset(
                        rn_tx_hgvs_not_delins, saved_hgvs_coding, with_base_subtract=True)

                # Logic
                hgvs_not_delins = self.logic_check(hgvs_not_delins, rn_tx_hgvs_not_delins, hgvs_coding)

            # Final sanity checks
            try:
                self.validator.vm.g_to_t(hgvs_not_delins, self.tx_hgvs_not_delins.ac,
                                         alt_aln_method=self.validator.alt_aln_method)
            except Exception as e:
                if str(e) == 'start or end or both are beyond the bounds of transcript record':
                    logger.warning(str(e))
                    return True
            try:
                self.variant.hn.normalize(self.tx_hgvs_not_delins)
            except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                error = str(e)
                if 'Normalization of intronic variants is not supported' in error or \
                        'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                    if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        logger.warning(error)
                        return True
                    elif 'Normalization of intronic variants is not supported' in error:
                        # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                        self.disparity_deletion_in = ['transcript', 'Requires Analysis']

            if hard_fail == 'true':
                self.disparity_deletion_in = ['false', 'false']

            # Recreate hgvs_genomic
            if self.disparity_deletion_in[0] == 'transcript':
                hgvs_genomic = hgvs_not_delins

            # Pre-processing of tx_hgvs_not_delins
            try:
                if self.tx_hgvs_not_delins.posedit.edit.alt is None:
                    self.tx_hgvs_not_delins.posedit.edit.alt = ''
            except Exception as e:
                if str(e) == "'Dup' object has no attribute 'alt'":
                    tx_hgvs_not_delins_delins_from_dup = fn.hgvs_dup2indel(self.tx_hgvs_not_delins)
                    self.tx_hgvs_not_delins = self.validator.hp.parse_hgvs_variant(tx_hgvs_not_delins_delins_from_dup)

            # GAP IN THE TRANSCRIPT DISPARITY DETECTED
            if self.disparity_deletion_in[0] == 'transcript':
                # ANY VARIANT WHOLLY WITHIN THE GAP
                hgvs_refreshed_variant = self.transcript_disparity(reverse_normalized_hgvs_genomic,
                                                                   stored_hgvs_not_delins, hgvs_genomic, 3)

            # GAP IN THE CHROMOSOME
            elif self.disparity_deletion_in[0] == 'chromosome':
                hgvs_refreshed_variant = self.tx_hgvs_not_delins

            else:
                # Keep the same by re-setting rel_var
                hgvs_refreshed_variant = saved_hgvs_coding

            # Edit the output
            hgvs_refreshed_variant = self.edit_output(hgvs_refreshed_variant, saved_hgvs_coding)

            # Sort out equality to equality c. events where the code will add 2 additional bases
            if hgvs_coding.posedit.edit.type == 'identity' and hgvs_refreshed_variant.posedit.edit.type == 'identity':
                pass
            else:
                hgvs_coding = copy.deepcopy(hgvs_refreshed_variant)

        return hgvs_coding

    def g_to_t_gap_compensation_version3(self, hgvs_alt_genomic, hgvs_coding, ori, alt_chr, rec_var):

        self.orientation = int(ori[0]['alt_strand'])
        hgvs_genomic = copy.deepcopy(hgvs_alt_genomic)

        logger.debug('g_to_t gap code 3 active')
        rn_hgvs_genomic = self.variant.reverse_normalizer.normalize(hgvs_alt_genomic)
        self.hgvs_genomic_possibilities.append([rn_hgvs_genomic, ['false', 'false']])
        if self.orientation != -1:
            try:
                chromosome_normalized_hgvs_coding = self.variant.reverse_normalizer.normalize(
                    hgvs_coding)
            except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                chromosome_normalized_hgvs_coding = hgvs_coding
        else:
            try:
                chromosome_normalized_hgvs_coding = self.variant.hn.normalize(hgvs_coding)
            except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                chromosome_normalized_hgvs_coding = hgvs_coding

        most_3pr_hgvs_genomic = self.validator.myvm_t_to_g(chromosome_normalized_hgvs_coding, alt_chr,
                                                           self.variant.no_norm_evm, self.variant.hn)
        self.hgvs_genomic_possibilities.append([most_3pr_hgvs_genomic, ['false', 'false']])

        # First to the right
        hgvs_stash = copy.deepcopy(hgvs_coding)
        stash_tx_right = ''
        stash_tx_left = ''

        # Capture instances where variant merging hard-sets the outputs
        try:
            hgvs_stash = self.variant.no_norm_evm.c_to_n(hgvs_stash)
        except Exception as e:
            logger.debug("Except passed, %s", e)
        try:
            stash_ac = hgvs_stash.ac
            stash_dict = hgvs_utils.hard_right_hgvs2vcf(hgvs_stash,
                                                        self.variant.primary_assembly,
                                                        self.variant.hn,
                                                        self.variant.reverse_normalizer,
                                                        self.validator.sf,
                                                        hgvs_coding.ac,
                                                        self.validator.hdp,
                                                        self.validator.alt_aln_method,
                                                        self.validator.hp,
                                                        self.validator.vm,
                                                        self.validator.merge_hgvs_3pr,
                                                        genomic_ac=hgvs_alt_genomic.ac)

            stash_pos = int(stash_dict['pos'])
            stash_ref = stash_dict['ref']
            stash_alt = stash_dict['alt']
            # Generate an end position
            stash_end = str(stash_pos + len(stash_ref) - 1)
            # make a not real deletion insertion
            stash_hgvs_not_delins = self.validator.hp.parse_hgvs_variant(
                stash_ac + ':' + hgvs_stash.type + '.' + str(
                    stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)
            try:
                stash_hgvs_not_delins = self.variant.no_norm_evm.n_to_c(stash_hgvs_not_delins)
            except Exception as e:
                logger.debug("Except passed, %s", e)

            # Store a tx copy for later use
            test_stash_tx_right = copy.deepcopy(stash_hgvs_not_delins)
            stash_genomic = self.validator.myvm_t_to_g(test_stash_tx_right, hgvs_alt_genomic.ac,
                                                       self.variant.no_norm_evm, self.variant.hn)

            if len(test_stash_tx_right.posedit.edit.ref) == ((stash_genomic.posedit.pos.end.base -
                                                              stash_genomic.posedit.pos.start.base) + 1):
                stash_tx_right = test_stash_tx_right
                if hasattr(test_stash_tx_right.posedit.edit,
                           'alt') and test_stash_tx_right.posedit.edit.alt is not None:
                    alt = test_stash_tx_right.posedit.edit.alt
                else:
                    alt = ''
                if hasattr(stash_genomic.posedit.edit,
                           'alt') and stash_genomic.posedit.edit.alt is not None:
                    g_alt = stash_genomic.posedit.edit.alt
                else:
                    g_alt = ''
                if (len(alt) - (
                        test_stash_tx_right.posedit.pos.end.base - test_stash_tx_right.posedit.pos.start.base) + 1) != (
                        len(g_alt) - (
                        stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])
                else:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
            elif test_stash_tx_right.posedit.edit.type == 'identity':
                reform_ident = str(test_stash_tx_right).split(':')[0]
                reform_ident = reform_ident + ':c.' + str(test_stash_tx_right.posedit.pos) + 'del' + str(
                    test_stash_tx_right.posedit.edit.ref)  # + 'ins' + str(test_stash_tx_right.posedit.edit.alt)
                hgvs_reform_ident = self.validator.hp.parse_hgvs_variant(reform_ident)
                try:
                    self.variant.hn.normalize(hgvs_reform_ident)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'spanning the exon-intron boundary' in error:
                        stash_tx_right = test_stash_tx_right
                        self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_right = test_stash_tx_right
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])
            else:
                try:
                    self.variant.hn.normalize(test_stash_tx_right)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_right = test_stash_tx_right
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])

            # IDENTIFY GAP AND HARD SET

            if stash_dict['needs_a_push'] is True or stash_dict['identifying_g_variant'] is not False:

                # Look for merged variant from hard push
                if stash_dict['merged_variant'] is not False:
                    merged_variant = stash_dict['merged_variant']
                    try:
                        merged_variant = self.validator.vm.n_to_c(merged_variant)
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass
                    identifying_variant = stash_dict['identifying_variant']
                    try:
                        identifying_variant = self.validator.vm.n_to_c(identifying_variant)
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass
                    stash_genomic = self.validator.myvm_t_to_g(identifying_variant, stash_genomic.ac,
                                                               self.variant.no_norm_evm,
                                                               self.variant.hn)
                    stash_hgvs_not_delins = merged_variant
                    test_stash_tx_right = merged_variant

                elif stash_dict['identifying_g_variant'] is not False:
                    stash_hgvs_not_delins = stash_dict['identifying_g_variant']
                    stash_genomic = stash_dict['identifying_g_variant']

                # Look for gap info
                normalized_stash_genomic = self.variant.hn.normalize(stash_genomic)
                stash_tx_right = test_stash_tx_right
                if stash_hgvs_not_delins.posedit.edit.type == "ins":
                    len_tx = 2
                else:
                    len_tx = len(stash_hgvs_not_delins.posedit.edit.ref)
                if stash_genomic.posedit.edit.type == "ins":
                    len_gen = 2
                else:
                    len_gen = len(stash_genomic.posedit.edit.ref)
                if len_tx > len_gen:
                    gap_in = 'chromosome'
                    gap_len = len_tx - len_gen
                else:
                    gap_in = 'transcript'
                    gap_len = len_gen - len_tx

                # Set the options to a single option based on the results of pushing
                self.hgvs_genomic_possibilities = [[normalized_stash_genomic, [gap_in,
                                                                               gap_len,
                                                                               stash_hgvs_not_delins,
                                                                               stash_genomic]]]

        except vvhgvs.exceptions.HGVSError as e:
            logger.debug("Except passed, %s", e)
        except ValueError as e:
            logger.debug("Except passed, %s", e)

        # Then to the left
        hgvs_stash = copy.deepcopy(hgvs_coding)
        try:
            hgvs_stash = self.variant.no_norm_evm.c_to_n(hgvs_stash)
        except Exception as e:
            logger.debug("Except passed, %s", e)
        try:
            stash_ac = hgvs_stash.ac
            stash_dict = hgvs_utils.hard_left_hgvs2vcf(hgvs_stash,
                                                       self.variant.primary_assembly,
                                                       self.variant.hn,
                                                       self.variant.reverse_normalizer,
                                                       self.validator.sf,
                                                       hgvs_coding.ac,
                                                       self.validator.hdp,
                                                       self.validator.alt_aln_method,
                                                       self.validator.hp,
                                                       self.validator.vm,
                                                       self.validator.merge_hgvs_3pr,
                                                       genomic_ac=hgvs_alt_genomic.ac)

            stash_pos = int(stash_dict['pos'])
            stash_ref = stash_dict['ref']
            stash_alt = stash_dict['alt']
            # Generate an end position
            stash_end = str(stash_pos + len(stash_ref) - 1)
            # make a not real deletion insertion
            stash_hgvs_not_delins = self.validator.hp.parse_hgvs_variant(
                stash_ac + ':' + hgvs_stash.type + '.' + str(
                    stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)
            try:
                stash_hgvs_not_delins = self.variant.no_norm_evm.n_to_c(stash_hgvs_not_delins)
            except Exception as e:
                logger.debug("Except passed, %s", e)
                # Store a tx copy for later use
            test_stash_tx_left = copy.deepcopy(stash_hgvs_not_delins)
            stash_genomic = self.validator.myvm_t_to_g(test_stash_tx_left, hgvs_alt_genomic.ac,
                                                  self.variant.no_norm_evm, self.variant.hn)

            if len(test_stash_tx_left.posedit.edit.ref) == ((stash_genomic.posedit.pos.end.base -
                                                             stash_genomic.posedit.pos.start.base) + 1):
                stash_tx_left = test_stash_tx_left
                if hasattr(test_stash_tx_left.posedit.edit,
                           'alt') and test_stash_tx_left.posedit.edit.alt is not None:
                    alt = test_stash_tx_left.posedit.edit.alt
                else:
                    alt = ''
                if hasattr(stash_genomic.posedit.edit,
                           'alt') and stash_genomic.posedit.edit.alt is not None:
                    g_alt = stash_genomic.posedit.edit.alt
                else:
                    g_alt = ''
                if (len(alt) - (
                        test_stash_tx_left.posedit.pos.end.base - test_stash_tx_left.posedit.pos.start.base) + 1) != (
                        len(g_alt) - (
                        stash_genomic.posedit.pos.end.base - stash_genomic.posedit.pos.start.base) + 1):
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])
                else:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
            elif test_stash_tx_left.posedit.edit.type == 'identity':
                reform_ident = str(test_stash_tx_left).split(':')[0]
                reform_ident = reform_ident + ':c.' + str(test_stash_tx_left.posedit.pos) + 'del' + str(
                    test_stash_tx_left.posedit.edit.ref)  # + 'ins' + str(test_stash_tx_left.posedit.edit.alt)
                hgvs_reform_ident = self.validator.hp.parse_hgvs_variant(reform_ident)
                try:
                    self.variant.hn.normalize(hgvs_reform_ident)
                except vvhgvs.exceptions.HGVSError as e:
                    error = str(e)
                    if 'spanning the exon-intron boundary' in error:
                        stash_tx_left = test_stash_tx_left
                        self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_left = test_stash_tx_left
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])

            else:
                try:
                    self.variant.hn.normalize(test_stash_tx_left)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                    self.hgvs_genomic_possibilities.append(['', ['false', 'false']])
                else:
                    stash_tx_left = test_stash_tx_left
                    self.hgvs_genomic_possibilities.append([stash_genomic, ['false', 'false']])

            # IDENTIFY GAP AND HARD SET
            if stash_dict['needs_a_push'] is True or stash_dict['identifying_g_variant'] is not False:

                # Look for merged variant from hard push
                if stash_dict['merged_variant'] is not False:
                    merged_variant = stash_dict['merged_variant']
                    try:
                        merged_variant = self.validator.vm.n_to_c(merged_variant)
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass
                    identifying_variant = stash_dict['identifying_variant']
                    try:
                        identifying_variant = self.validator.vm.n_to_c(identifying_variant)
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass
                    stash_genomic = self.validator.myvm_t_to_g(identifying_variant, stash_genomic.ac,
                                                               self.variant.no_norm_evm,
                                                               self.variant.hn)
                    stash_hgvs_not_delins = merged_variant
                    test_stash_tx_left = merged_variant

                elif stash_dict['identifying_g_variant'] is not False:
                    stash_hgvs_not_delins = stash_dict['identifying_g_variant']
                    stash_genomic = stash_dict['identifying_g_variant']

                # Look for gap info
                normalized_stash_genomic = self.variant.hn.normalize(stash_genomic)
                stash_tx_left = test_stash_tx_left
                if stash_hgvs_not_delins.posedit.edit.type == "ins":
                    len_tx = 2
                else:
                    len_tx = len(stash_hgvs_not_delins.posedit.edit.ref)
                if stash_genomic.posedit.edit.type == "ins":
                    len_gen = 2
                else:
                    len_gen = len(stash_genomic.posedit.edit.ref)
                if len_tx > len_gen:
                    gap_in = 'chromosome'
                    gap_len = len_tx - len_gen
                else:
                    gap_in = 'transcript'
                    gap_len = len_gen - len_tx

                # Set the options to a single option based on the results of pushing
                self.hgvs_genomic_possibilities = [[normalized_stash_genomic, [gap_in,
                                                                               gap_len,
                                                                               stash_hgvs_not_delins,
                                                                               stash_genomic]]]

        except vvhgvs.exceptions.HGVSError as e:
            logger.debug("Except passed, %s", e)
        except ValueError as e:
            logger.debug("Except passed, %s", e)

        # direct mapping from reverse_normalized transcript insertions in the delins format
        self.rev_norm_ins(hgvs_coding, hgvs_genomic)

        # Set variables for problem specific warnings
        self.gapped_transcripts = ''
        self.auto_info = ''

        # Mark as not disparity detected
        self.disparity_deletion_in = ['false', 'false']
        # Loop through to see if a gap can be located
        possibility_counter = 0

        # If hard pushing identified a variant (in theory there can be only 1) then remove all other possibilities
        hard_possibility = []
        for hard_set_check in self.hgvs_genomic_possibilities:
            if hard_set_check[1] != ['false', 'false']:
                hard_possibility.append(hard_set_check)
        if len(hard_possibility) >= 1:
            self.hgvs_genomic_possibilities = hard_possibility

        for a_possibility in self.hgvs_genomic_possibilities:
            possibility = a_possibility[0]
            disparity_info = a_possibility[1]
            possibility_counter = possibility_counter + 1
            # Loop out stash possibilities which will not spot gaps so are empty
            if possibility == '':
                continue

            # Use VCF generation code to push hgvs_genomic as for 5 prime as possible to uncover gaps
            hgvs_genomic_variant = possibility
            reverse_normalized_hgvs_genomic = ''

            # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
            try:
                reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(
                    hgvs_genomic_variant)
            except vvhgvs.exceptions.HGVSError as e:
                # Strange error caused by gap in genomic
                error = str(e)
                if 'base start position must be <= end position' in error:
                    if hgvs_genomic.posedit.edit.type == 'delins':
                        start = hgvs_genomic.posedit.pos.start.base
                        end = hgvs_genomic.posedit.pos.end.base
                        lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                        rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                        hgvs_genomic.posedit.edit.ref = lhb + rhb
                        hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                        hgvs_genomic.posedit.pos.start.base = end
                        hgvs_genomic.posedit.pos.end.base = start
                        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(
                            hgvs_genomic)
                    if hgvs_genomic.posedit.edit.type == 'del':
                        start = hgvs_genomic.posedit.pos.start.base
                        end = hgvs_genomic.posedit.pos.end.base
                        lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), end - 1, end)
                        rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                        hgvs_genomic.posedit.edit.ref = lhb + rhb
                        hgvs_genomic.posedit.edit.alt = lhb + rhb
                        hgvs_genomic.posedit.pos.start.base = end
                        hgvs_genomic.posedit.pos.end.base = start
                        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(
                            hgvs_genomic)
                if 'insertion length must be 1' in error:
                    if hgvs_genomic.posedit.edit.type == 'ins':
                        start = hgvs_genomic.posedit.pos.start.base
                        end = hgvs_genomic.posedit.pos.end.base
                        ref_bases = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, end)
                        lhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start - 1, start)
                        rhb = self.validator.sf.fetch_seq(str(hgvs_genomic.ac), start, end)
                        hgvs_genomic.posedit.edit.ref = lhb + rhb
                        hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
                        reverse_normalized_hgvs_genomic = self.variant.reverse_normalizer.normalize(
                            hgvs_genomic)

            self.hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)
            # Store a copy for later use

            # Make VCF
            vcf_dict = hgvs_utils.hgvs2vcf(reverse_normalized_hgvs_genomic, self.variant.primary_assembly,
                                           self.variant.reverse_normalizer, self.validator.sf)
            pos = vcf_dict['pos']
            ref = vcf_dict['ref']
            alt = vcf_dict['alt']

            # Look for exonic gaps within transcript or chromosome

            # Generate an end position
            end = str(int(pos) + len(ref) - 1)
            pos = str(pos)

            # Store a not real deletion insertion to test for gapping
            stored_hgvs_not_delins = self.validator.hp.parse_hgvs_variant(str(
                self.hgvs_genomic_5pr.ac) + ':' + self.hgvs_genomic_5pr.type + '.' + pos + '_' + end + 'del' +
                                                                          ref + 'ins' + alt)
            # Save a copy of current hgvs_coding
            saved_hgvs_coding = ''
            try:
                saved_hgvs_coding = self.variant.no_norm_evm.g_to_t(stored_hgvs_not_delins,
                                                                    hgvs_coding.ac)
            except Exception as e:
                if str(e) == 'start or end or both are beyond the bounds of transcript record':
                    continue

            # Detect intronic variation using normalization
            intronic_variant = 'false'
            # Look for normalized variant options that do not match hgvs_coding
            hgvs_seek_var = self.get_hgvs_seek_var(hgvs_genomic, hgvs_coding)
            if (
                    hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                    hgvs_coding.posedit.pos.start.base + hgvs_coding.posedit.pos.start.offset) and (
                    hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                    hgvs_coding.posedit.pos.end.base + hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
                pass
            else:
                hgvs_seek_var = saved_hgvs_coding

            try:
                self.variant.hn.normalize(hgvs_seek_var)
            except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                error = str(e)
                if 'Normalization of intronic variants is not supported' in error or \
                        'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                    if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        intronic_variant = 'hard_fail'
                    else:
                        # Double check to see whether the variant is actually intronic?
                        for exon in ori:
                            genomic_start = int(exon['alt_start_i'])
                            genomic_end = int(exon['alt_end_i'])
                            if (genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                                    genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                                intronic_variant = 'false'
                                break
                            else:
                                intronic_variant = 'true'

            if intronic_variant != 'hard_fail':
                if re.search(r'\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search(r'\d+-',
                                                                                    str(
                                                                                        hgvs_seek_var.posedit.pos)) or re.search(
                    r'\*\d+\+', str(
                        hgvs_seek_var.posedit.pos)) or re.search(r'\*\d+-', str(hgvs_seek_var.posedit.pos)):
                    # Double check to see whether the variant is actually intronic?
                    for exon in ori:
                        genomic_start = int(exon['alt_start_i'])
                        genomic_end = int(exon['alt_end_i'])
                        if (genomic_start < self.hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                                genomic_start < self.hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                            intronic_variant = 'false'
                            break
                        else:
                            intronic_variant = 'true'

            if intronic_variant != 'true':
                # Flag RefSeqGene for ammendment
                # amend_RefSeqGene = 'false'
                # Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping
                # transcript lengths
                hgvs_not_delins = ''
                if stored_hgvs_not_delins != '':
                    # Refresh hgvs_not_delins from stored_hgvs_not_delins
                    hgvs_not_delins = self.dup_ins_5prime_shift(stored_hgvs_not_delins, saved_hgvs_coding)

                    self.tx_hgvs_not_delins = self.variant.no_norm_evm.g_to_n(hgvs_not_delins,
                                                                              saved_hgvs_coding.ac)
                    # Create normalized version of tx_hgvs_not_delins
                    rn_tx_hgvs_not_delins = copy.deepcopy(self.tx_hgvs_not_delins)
                    # Check for +1 base and adjust
                    if '+' in str(rn_tx_hgvs_not_delins.posedit.pos.end) and '+' in str(
                            rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins = self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins)

                    elif '+' in str(rn_tx_hgvs_not_delins.posedit.pos.end):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_end_base_to_next_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding, back=False)

                    elif '+' in str(rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_start_base_to_previous_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding)

                    # Check for -ve base and adjust
                    elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.end) and '-' in str(
                            rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins = self.remove_offsetting_to_span_gap(rn_tx_hgvs_not_delins)

                    elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.end):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_end_base_to_next_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding)

                    elif '-' in str(rn_tx_hgvs_not_delins.posedit.pos.start):
                        rn_tx_hgvs_not_delins, hgvs_not_delins = self.move_tx_start_base_to_previous_nonoffset(
                            rn_tx_hgvs_not_delins, saved_hgvs_coding, with_base_subtract=True)

                    # Logic

                    hgvs_not_delins = self.logic_check(hgvs_not_delins, rn_tx_hgvs_not_delins, hgvs_coding,
                                                       do_continue=True, offset_check=True)

                # Final sanity checks
                try:
                    self.validator.vm.g_to_t(hgvs_not_delins,
                                             self.tx_hgvs_not_delins.ac,
                                             alt_aln_method=self.validator.alt_aln_method)
                except Exception as e:
                    if str(e) == 'start or end or both are beyond the bounds of transcript record':
                        continue
                try:
                    self.variant.hn.normalize(self.tx_hgvs_not_delins)
                except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
                    error = str(e)
                    if 'Normalization of intronic variants is not supported' in error or \
                            'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        if 'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                            continue
                        elif 'Normalization of intronic variants is not supported' in error:
                            # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                            self.disparity_deletion_in = ['transcript', 'Requires Analysis']

                # Recreate hgvs_genomic
                if self.disparity_deletion_in[0] == 'transcript':
                    hgvs_genomic = hgvs_not_delins

                # Find oddly placed gaps where the tx variant is encompassed in the gap
                if self.disparity_deletion_in[0] == 'false' and (
                        possibility_counter == 3 or possibility_counter == 4):
                    rg = self.variant.reverse_normalizer.normalize(hgvs_not_delins)
                    rtx = self.validator.vm.g_to_t(rg, self.tx_hgvs_not_delins.ac,
                                                   alt_aln_method=self.validator.alt_aln_method)
                    fg = self.variant.hn.normalize(hgvs_not_delins)
                    ftx = self.validator.vm.g_to_t(fg, self.tx_hgvs_not_delins.ac,
                                                   alt_aln_method=self.validator.alt_aln_method)
                    if (rtx.posedit.pos.start.offset == 0 and rtx.posedit.pos.end.offset == 0) and (
                            ftx.posedit.pos.start.offset != 0 and ftx.posedit.pos.end.offset != 0):
                        exons = self.validator.hdp.get_tx_exons(ftx.ac, hgvs_not_delins.ac,
                                                                self.validator.alt_aln_method,
                                                                alt_aln_method=self.validator.alt_aln_method)
                        exonic = False
                        for ex_test in exons:
                            if ftx.posedit.pos.start.base in range(ex_test[6], ex_test[7]) and \
                                    ftx.posedit.pos.end.base in range(ex_test[6], ex_test[7]):
                                exonic = True
                        if exonic is True:
                            hgvs_genomic = fg
                            self.hgvs_genomic_5pr = fg
                            try:
                                self.tx_hgvs_not_delins = self.validator.vm.c_to_n(ftx)
                            except Exception:
                                self.tx_hgvs_not_delins = ftx
                            self.disparity_deletion_in = ['transcript', 'Requires Analysis']

                # Pre-processing of tx_hgvs_not_delins
                try:
                    if self.tx_hgvs_not_delins.posedit.edit.alt is None:
                        self.tx_hgvs_not_delins.posedit.edit.alt = ''
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        tx_hgvs_not_delins_delins_from_dup = fn.hgvs_dup2indel(self.tx_hgvs_not_delins)
                        self.tx_hgvs_not_delins = self.validator.hp.parse_hgvs_variant(
                            tx_hgvs_not_delins_delins_from_dup)

                # Has a hard set variant been identified from pushes?
                hard_set_outputs = False
                if disparity_info != ['false', 'false'] and len(disparity_info) == 4:
                    self.tx_hgvs_not_delins = disparity_info[2]
                    self.disparity_deletion_in = [disparity_info[0], disparity_info[1]]
                    hgvs_refreshed_variant = hgvs_coding
                    hgvs_alt_genomic = possibility
                    # self.variant.warnings.append("Caution should be used when reporting the displayed variant "
                    #                              "descriptions: If you are unsure, please contact admin")
                    # self.variant.warnings.append('The displayed variants may be artefacts of aligning '
                    #                              '' + hgvs_coding.ac + ' with genomic reference '
                    #                                                    '' + disparity_info[3].ac)
                    hard_set_outputs = True

                elif self.disparity_deletion_in[0] == 'transcript':
                    # ANY VARIANT WHOLLY WITHIN THE GAP
                    hgvs_refreshed_variant = self.transcript_disparity(reverse_normalized_hgvs_genomic,
                                                                       stored_hgvs_not_delins, hgvs_genomic, 4)

                # GAP IN THE CHROMOSOME
                elif self.disparity_deletion_in[0] == 'chromosome':
                    # amend_RefSeqGene = 'true'
                    if possibility_counter == 3:
                        hgvs_refreshed_variant = stash_tx_right
                    elif possibility_counter == 4:
                        hgvs_refreshed_variant = stash_tx_left
                    else:
                        hgvs_refreshed_variant = chromosome_normalized_hgvs_coding
                else:
                    # Keep the same by re-setting rel_var
                    hgvs_refreshed_variant = hgvs_coding

                # Edit the output
                if 'NM_' in str(hgvs_refreshed_variant.ac) and 'c' not in str(hgvs_refreshed_variant.type):
                    hgvs_refreshed_variant = self.variant.no_norm_evm.n_to_c(hgvs_refreshed_variant)

                try:
                    self.variant.hn.normalize(hgvs_refreshed_variant)
                except Exception as e:
                    error = str(e)
                    # Ensure the final variant is not intronic nor does it cross exon boundaries
                    if 'Normalization of intronic variants is not supported' in error or \
                            'Unsupported normalization of variants spanning the exon-intron boundary' in error:
                        hgvs_refreshed_variant = saved_hgvs_coding
                    else:
                        continue

                # Quick check to make sure the coding variant has not changed UNLESS HARD SET
                try:
                    to_test = self.variant.hn.normalize(hgvs_refreshed_variant)
                except:
                    to_test = hgvs_refreshed_variant
                if str(to_test.posedit.edit) != str(hgvs_coding.posedit.edit):
                    # Try the next available genomic option
                    if hgvs_coding.posedit.edit.type == 'identity' and to_test.posedit.edit.type == 'identity':
                        hgvs_coding = to_test
                    elif hard_set_outputs is True:
                        hgvs_coding = to_test
                    else:
                        continue

                if hard_set_outputs is False:
                    # Update hgvs_genomic
                    hgvs_alt_genomic = self.validator.myvm_t_to_g(hgvs_refreshed_variant, alt_chr,
                                                                  self.variant.no_norm_evm, self.variant.hn)
                    if hgvs_alt_genomic.posedit.edit.type == 'identity':
                        re_c = self.validator.vm.g_to_t(hgvs_alt_genomic, hgvs_refreshed_variant.ac,
                                                        alt_aln_method=self.validator.alt_aln_method)
                        if (self.variant.hn.normalize(re_c)) != (self.variant.hn.normalize(hgvs_refreshed_variant)):
                            shuffle_left_g = copy.copy(hgvs_alt_genomic)
                            shuffle_left_g.posedit.edit.ref = ''
                            shuffle_left_g.posedit.edit.alt = ''
                            shuffle_left_g.posedit.pos.start.base = shuffle_left_g.posedit.pos.start.base - 1
                            shuffle_left_g.posedit.pos.end.base = shuffle_left_g.posedit.pos.end.base - 1
                            shuffle_left_g = self.variant.reverse_normalizer.normalize(shuffle_left_g)
                            re_c = self.validator.vm.g_to_t(shuffle_left_g, hgvs_refreshed_variant.ac,
                                                            alt_aln_method=self.validator.alt_aln_method)
                            if (self.variant.hn.normalize(re_c)) != (self.variant.hn.normalize(hgvs_refreshed_variant)):
                                hgvs_alt_genomic = shuffle_left_g

                                # If it is intronic, these vairables will not have been set

            # Break if gap has been detected
            if self.disparity_deletion_in[0] != 'false':
                break

        # Normailse hgvs_genomic
        try:
            hgvs_alt_genomic = self.variant.hn.normalize(hgvs_alt_genomic)
        except vvhgvs.exceptions.HGVSError as e:
            # Strange error caused by gap in genomic
            error = str(e)
            if 'base start position must be <= end position' in error and self.disparity_deletion_in[0] == 'chromosome':
                if hgvs_alt_genomic.posedit.edit.type == 'delins':
                    start = hgvs_alt_genomic.posedit.pos.start.base
                    end = hgvs_alt_genomic.posedit.pos.end.base
                    lhb = self.validator.sf.fetch_seq(str(hgvs_alt_genomic.ac), end - 1, end)
                    rhb = self.validator.sf.fetch_seq(str(hgvs_alt_genomic.ac), start - 1, start)
                    hgvs_alt_genomic.posedit.edit.ref = lhb + rhb
                    hgvs_alt_genomic.posedit.edit.alt = lhb + hgvs_alt_genomic.posedit.edit.alt + rhb
                    hgvs_alt_genomic.posedit.pos.start.base = end
                    hgvs_alt_genomic.posedit.pos.end.base = start
                    hgvs_alt_genomic = self.variant.hn.normalize(hgvs_alt_genomic)
                if hgvs_alt_genomic.posedit.edit.type == 'del':
                    start = hgvs_alt_genomic.posedit.pos.start.base
                    end = hgvs_alt_genomic.posedit.pos.end.base
                    lhb = self.validator.sf.fetch_seq(str(hgvs_alt_genomic.ac), end - 1, end)
                    rhb = self.validator.sf.fetch_seq(str(hgvs_alt_genomic.ac), start - 1, start)
                    hgvs_alt_genomic.posedit.edit.ref = lhb + rhb
                    hgvs_alt_genomic.posedit.edit.alt = lhb + rhb
                    hgvs_alt_genomic.posedit.pos.start.base = end
                    hgvs_alt_genomic.posedit.pos.end.base = start
                    hgvs_alt_genomic = self.variant.hn.normalize(hgvs_alt_genomic)

        # check for flanking substitutions which should be dels due to gap in transtipt
        try:
            check_flank_genomic = self.validator.myvm_t_to_g(hgvs_refreshed_variant,
                                                             hgvs_alt_genomic.ac,
                                                             self.variant.no_norm_evm,
                                                             self.variant.hn)

            if ((hgvs_alt_genomic.posedit.edit.type == hgvs_refreshed_variant.posedit.edit.type and
                    hgvs_alt_genomic.posedit.edit.type == "sub") and
                    check_flank_genomic.posedit.edit.type == "del"):

                hgvs_alt_genomic = check_flank_genomic
        except UnboundLocalError:
            pass

        return hgvs_alt_genomic, hgvs_coding

    def dup_ins_5prime_shift(self, stored_hgvs_not_delins, saved_hgvs_coding):
        hgvs_not_delins = copy.deepcopy(stored_hgvs_not_delins)
        # This test will only occur in dup of single base, insertion or substitution
        if '_' not in str(hgvs_not_delins.posedit.pos):
            if 'dup' in self.hgvs_genomic_5pr.posedit.edit.type or 'ins' in self.hgvs_genomic_5pr.posedit.edit.type:
                # For gap in chr, map to t. - but because we have pushed to 5 prime by norm, add 1 to end pos
                plussed_hgvs_not_delins = copy.deepcopy(hgvs_not_delins)
                plussed_hgvs_not_delins.posedit.pos.end.base = plussed_hgvs_not_delins.posedit.pos.end.base + 1
                plussed_hgvs_not_delins.posedit.edit.ref = ''
                transcript_variant = self.variant.no_norm_evm.g_to_t(plussed_hgvs_not_delins,
                                                                     str(saved_hgvs_coding.ac))
                if ((transcript_variant.posedit.pos.end.base - transcript_variant.posedit.pos.start.base) > (
                        self.hgvs_genomic_5pr.posedit.pos.end.base - self.hgvs_genomic_5pr.posedit.pos.start.base)):
                    if 'dup' in str(self.hgvs_genomic_5pr.posedit.edit):
                        hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                        start = hgvs_not_delins.posedit.pos.start.base - 1
                        end = hgvs_not_delins.posedit.pos.end.base
                        ref_bases = self.validator.sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                        hgvs_not_delins.posedit.edit.ref = ref_bases
                        hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[
                                                                 1:] + ref_bases[1:]
                    elif 'ins' in str(self.hgvs_genomic_5pr.posedit.edit) and \
                            'del' in str(self.hgvs_genomic_5pr.posedit.edit):
                        hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                    elif 'ins' in str(self.hgvs_genomic_5pr.posedit.edit) and \
                            'del' not in str(self.hgvs_genomic_5pr.posedit.edit):
                        hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                        start = hgvs_not_delins.posedit.pos.start.base - 1
                        end = hgvs_not_delins.posedit.pos.end.base
                        ref_bases = self.validator.sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                        hgvs_not_delins.posedit.edit.ref = ref_bases
                        hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[
                                                                 1:] + ref_bases[1:]
                else:
                    if 'dup' in str(self.hgvs_genomic_5pr.posedit.edit):
                        hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                        start = hgvs_not_delins.posedit.pos.start.base - 1
                        end = hgvs_not_delins.posedit.pos.end.base
                        ref_bases = self.validator.sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                        hgvs_not_delins.posedit.edit.ref = ref_bases
                        hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[
                                                                 1:] + ref_bases[1:]
                    elif 'ins' in str(self.hgvs_genomic_5pr.posedit.edit) and \
                            'del' in str(self.hgvs_genomic_5pr.posedit.edit):
                        hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                    elif 'ins' in str(self.hgvs_genomic_5pr.posedit.edit) and \
                            'del' not in str(self.hgvs_genomic_5pr.posedit.edit):
                        hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                        start = hgvs_not_delins.posedit.pos.start.base - 1
                        end = hgvs_not_delins.posedit.pos.end.base
                        ref_bases = self.validator.sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                        hgvs_not_delins.posedit.edit.ref = ref_bases
                        hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[
                                                                 1:] + ref_bases[1:]

        return hgvs_not_delins

    def remove_offsetting_to_span_gap(self, rn_tx_hgvs_not_delins):
        # Remove offsetting to span the gap
        rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
        rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
        rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
        rn_tx_hgvs_not_delins.posedit.edit.ref = ''
        try:
            rn_tx_hgvs_not_delins.posedit.edit.alt = ''
        except Exception as e:
            logger.debug("Except passed, %s", e)

        return rn_tx_hgvs_not_delins

    def move_tx_end_base_to_next_nonoffset(self, rn_tx_hgvs_not_delins, saved_hgvs_coding, back=True):
        # move tx end base back to next available non-offset base
        rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
        rn_tx_hgvs_not_delins.posedit.edit.ref = ''

        if back:
            # Add the additional base to the ALT
            start = rn_tx_hgvs_not_delins.posedit.pos.end.base - 1
            end = rn_tx_hgvs_not_delins.posedit.pos.end.base
            ref_bases = self.validator.sf.fetch_seq(str(self.tx_hgvs_not_delins.ac), start, end)
            rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + ref_bases
        else:
            # move tx end base to next available non-offset base
            rn_tx_hgvs_not_delins.posedit.pos.end.base = self.tx_hgvs_not_delins.posedit.pos.end.base + 1
        if 'NM_' in str(rn_tx_hgvs_not_delins):
            test_tx_var = self.variant.no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
        else:
            test_tx_var = rn_tx_hgvs_not_delins
        # re-make genomic and tx
        hgvs_not_delins = self.validator.myevm_t_to_g(test_tx_var, self.variant.no_norm_evm,
                                                      self.variant.primary_assembly, self.variant.hn)
        rn_tx_hgvs_not_delins = self.variant.no_norm_evm.g_to_n(hgvs_not_delins,
                                                                str(saved_hgvs_coding.ac))
        return rn_tx_hgvs_not_delins, hgvs_not_delins

    def move_tx_start_base_to_previous_nonoffset(self, rn_tx_hgvs_not_delins, saved_hgvs_coding,
                                                 with_base_subtract=False):

        # Store the original variant
        store_rn_tx_hgvs_not_delins = copy.deepcopy(rn_tx_hgvs_not_delins)
        # move tx start base to previous available non-offset base
        rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
        if with_base_subtract and rn_tx_hgvs_not_delins.posedit.pos.start.base > 1:
            rn_tx_hgvs_not_delins.posedit.pos.start.base = rn_tx_hgvs_not_delins.posedit.pos.start.base - 1
        rn_tx_hgvs_not_delins.posedit.edit.ref = ''

        if 'NM_' in str(rn_tx_hgvs_not_delins):
            try:
                test_tx_var = self.variant.no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
            except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                if "Expected n. variant;" in str(e):
                    rn_tx_hgvs_not_delins = self.validator.vm.c_to_n(rn_tx_hgvs_not_delins)
                    test_tx_var = self.variant.no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)

        else:
            test_tx_var = rn_tx_hgvs_not_delins

        # re-make genomic and tx
        hgvs_not_delins = self.validator.myevm_t_to_g(test_tx_var, self.variant.no_norm_evm,
                                                      self.variant.primary_assembly, self.variant.hn)

        try:
            rn_tx_hgvs_not_delins = self.variant.no_norm_evm.g_to_n(hgvs_not_delins,
                                                                str(saved_hgvs_coding.ac))
        except vvhgvs.exceptions.HGVSInvalidIntervalError:
            rn_tx_hgvs_not_delins = test_tx_var

        if store_rn_tx_hgvs_not_delins.posedit.pos.start.offset != 0 and \
                rn_tx_hgvs_not_delins.posedit.pos.start.offset == 0 \
                and store_rn_tx_hgvs_not_delins.posedit.edit.ref == rn_tx_hgvs_not_delins.posedit.edit.ref:
            offset = store_rn_tx_hgvs_not_delins.posedit.pos.start.offset
            add_in_these_bases = store_rn_tx_hgvs_not_delins.posedit.edit.ref[0:0+offset]
            rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + add_in_these_bases
        else:
            rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0

        return rn_tx_hgvs_not_delins, hgvs_not_delins

    def c2_pos_edit(self, hgvs_genomic):
        try:
            c2 = self.validator.vm.n_to_c(self.tx_hgvs_not_delins)
        except:
            c2 = self.tx_hgvs_not_delins
        c1 = copy.deepcopy(c2)
        c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
        c1.posedit.pos.start.offset = 0
        c1.posedit.pos.end = c2.posedit.pos.start
        c1.posedit.edit.ref = ''
        c1.posedit.edit.alt = ''
        if self.orientation != -1:
            g1 = self.validator.vm.t_to_g(c1, hgvs_genomic.ac, alt_aln_method=self.validator.alt_aln_method)
            g2 = self.validator.vm.t_to_g(c2, hgvs_genomic.ac, alt_aln_method=self.validator.alt_aln_method)
            g1.posedit.edit.alt = g1.posedit.edit.ref
        else:
            g1 = self.validator.vm.t_to_g(c2, hgvs_genomic.ac, alt_aln_method=self.validator.alt_aln_method)
            g2 = self.validator.vm.t_to_g(c1, hgvs_genomic.ac, alt_aln_method=self.validator.alt_aln_method)
            g2.posedit.edit.alt = g2.posedit.edit.ref
        reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
        alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
        g3 = copy.deepcopy(g1)
        g3.posedit.pos.end.base = g2.posedit.pos.end.base
        g3.posedit.edit.ref = reference
        g3.posedit.edit.alt = alternate
        c3 = self.validator.vm.g_to_t(g3, c1.ac, alt_aln_method=self.validator.alt_aln_method)
        hgvs_refreshed_variant = c3

        return hgvs_refreshed_variant

    def c1_pos_edit(self, hgvs_genomic):
        try:
            c1 = self.validator.vm.n_to_c(self.tx_hgvs_not_delins)
        except:
            c1 = self.tx_hgvs_not_delins

        c2 = copy.deepcopy(c1)
        c2.posedit.pos.start = c1.posedit.pos.end
        c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
        c2.posedit.pos.end.offset = 0
        c2.posedit.edit.ref = ''
        c2.posedit.edit.alt = ''

        if self.orientation != -1:
            g1 = self.validator.vm.t_to_g(c1, hgvs_genomic.ac, alt_aln_method=self.validator.alt_aln_method)
            g2 = self.validator.vm.t_to_g(c2, hgvs_genomic.ac, alt_aln_method=self.validator.alt_aln_method)
            g2.posedit.edit.alt = g2.posedit.edit.ref
        else:
            g1 = self.validator.vm.t_to_g(c2, hgvs_genomic.ac, alt_aln_method=self.validator.alt_aln_method)
            g2 = self.validator.vm.t_to_g(c1, hgvs_genomic.ac, alt_aln_method=self.validator.alt_aln_method)
            g1.posedit.edit.alt = g1.posedit.edit.ref

        reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
        alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]

        g3 = copy.deepcopy(g1)
        g3.posedit.pos.end.base = g2.posedit.pos.end.base
        g3.posedit.edit.ref = reference
        g3.posedit.edit.alt = alternate
        c3 = self.validator.vm.g_to_t(g3, c1.ac, alt_aln_method=self.validator.alt_aln_method)
        hgvs_refreshed_variant = c3

        return hgvs_refreshed_variant

    def transcript_disparity(self, reverse_normalized_hgvs_genomic, stored_hgvs_not_delins, hgvs_genomic,
                             running_option):

        if ('+' in str(self.tx_hgvs_not_delins.posedit.pos.start) or '-' in str(
                self.tx_hgvs_not_delins.posedit.pos.start)) and (
                '+' in str(self.tx_hgvs_not_delins.posedit.pos.end) or '-' in str(
                self.tx_hgvs_not_delins.posedit.pos.end)):
            self.gapped_transcripts = self.gapped_transcripts + ' ' + str(self.tx_hgvs_not_delins.ac)

            # Copy the current variant
            tx_gap_fill_variant = copy.deepcopy(self.tx_hgvs_not_delins)
            try:
                if tx_gap_fill_variant.posedit.edit.alt is None:
                    tx_gap_fill_variant.posedit.edit.alt = ''
            except Exception as e:
                if str(e) == "'Dup' object has no attribute 'alt'":
                    tx_gap_fill_variant_delins_from_dup = fn.hgvs_dup2indel(tx_gap_fill_variant)
                    tx_gap_fill_variant = self.validator.hp.parse_hgvs_variant(
                        tx_gap_fill_variant_delins_from_dup)

            # Identify which half of the NOT-intron the start position of the variant is in
            if '-' in str(tx_gap_fill_variant.posedit.pos.start):
                tx_gap_fill_variant.posedit.pos.start.base = tx_gap_fill_variant.posedit.pos.start.base - 1
                tx_gap_fill_variant.posedit.pos.start.offset = int('0')  # int('+1')
                tx_gap_fill_variant.posedit.pos.end.offset = int('0')  # int('-1')
                tx_gap_fill_variant.posedit.edit.alt = ''
                tx_gap_fill_variant.posedit.edit.ref = ''
            elif '+' in str(tx_gap_fill_variant.posedit.pos.start):
                tx_gap_fill_variant.posedit.pos.start.offset = int('0')  # int('+1')
                tx_gap_fill_variant.posedit.pos.end.base = tx_gap_fill_variant.posedit.pos.end.base + 1
                tx_gap_fill_variant.posedit.pos.end.offset = int('0')  # int('-1')
                tx_gap_fill_variant.posedit.edit.alt = ''
                tx_gap_fill_variant.posedit.edit.ref = ''

            try:
                tx_gap_fill_variant = self.validator.vm.n_to_c(tx_gap_fill_variant)
            except Exception as e:
                logger.debug("Except passed, %s", e)
            genomic_gap_fill_variant = self.validator.vm.t_to_g(tx_gap_fill_variant, reverse_normalized_hgvs_genomic.ac,
                                                                alt_aln_method=self.validator.alt_aln_method)
            genomic_gap_fill_variant.posedit.edit.alt = genomic_gap_fill_variant.posedit.edit.ref

            try:
                c_tx_hgvs_not_delins = self.validator.vm.n_to_c(self.tx_hgvs_not_delins)
            except Exception:
                c_tx_hgvs_not_delins = copy.copy(self.tx_hgvs_not_delins)
            genomic_gap_fill_variant_alt = self.validator.vm.t_to_g(c_tx_hgvs_not_delins, self.hgvs_genomic_5pr.ac,
                                                                    alt_aln_method=self.validator.alt_aln_method)

            # Ensure an ALT exists
            try:
                if genomic_gap_fill_variant_alt.posedit.edit.alt is None:
                    genomic_gap_fill_variant_alt.posedit.edit.alt = 'X'
            except Exception as e:
                if str(e) == "'Dup' object has no attribute 'alt'":
                    genomic_gap_fill_variant_delins_from_dup = fn.hgvs_dup2indel(genomic_gap_fill_variant)
                    genomic_gap_fill_variant = self.validator.hp.parse_hgvs_variant(
                        genomic_gap_fill_variant_delins_from_dup)
                    genomic_gap_fill_variant_alt_delins_from_dup = fn.hgvs_dup2indel(genomic_gap_fill_variant_alt)
                    genomic_gap_fill_variant_alt = self.validator.hp.parse_hgvs_variant(
                        genomic_gap_fill_variant_alt_delins_from_dup)

            # Correct insertion alts
            if genomic_gap_fill_variant_alt.posedit.edit.type == 'ins':
                append_ref = self.validator.sf.fetch_seq(genomic_gap_fill_variant_alt.ac,
                                                         genomic_gap_fill_variant_alt.posedit.pos.start.base - 1,
                                                         genomic_gap_fill_variant_alt.posedit.pos.end.base)
                genomic_gap_fill_variant_alt.posedit.edit.alt = append_ref[0] + \
                    genomic_gap_fill_variant_alt.posedit.edit.alt + append_ref[1]

            # Split the reference and replacing alt sequence into a dictionary
            reference_bases = list(genomic_gap_fill_variant.posedit.edit.ref)
            if genomic_gap_fill_variant_alt.posedit.edit.alt is not None:
                alternate_bases = list(genomic_gap_fill_variant_alt.posedit.edit.alt)
            else:
                # Deletions with no ins
                pre_alternate_bases = list(genomic_gap_fill_variant_alt.posedit.edit.ref)
                alternate_bases = []
                for base in pre_alternate_bases:
                    alternate_bases.append('X')

            # Create the dictionaries
            ref_start = genomic_gap_fill_variant.posedit.pos.start.base
            alt_start = genomic_gap_fill_variant_alt.posedit.pos.start.base
            ref_base_dict = {}
            for base in reference_bases:
                ref_base_dict[ref_start] = str(base)
                ref_start = ref_start + 1

            alt_base_dict = {}

            # NEED TO SEARCH FOR RANGE = and replace with interval_range
            # Need to search for int and replace with integer

            # Note, all variants will be forced into the format delete insert
            # Deleted bases in the ALT will be substituted for X
            for integer in range(genomic_gap_fill_variant_alt.posedit.pos.start.base,
                                 genomic_gap_fill_variant_alt.posedit.pos.end.base + 1, 1):
                if integer == alt_start:
                    alt_base_dict[integer] = str(''.join(alternate_bases))
                else:
                    alt_base_dict[integer] = 'X'

            # Generate the alt sequence
            alternate_sequence_bases = []
            for integer in range(genomic_gap_fill_variant.posedit.pos.start.base,
                                 genomic_gap_fill_variant.posedit.pos.end.base + 1, 1):
                if integer in list(alt_base_dict.keys()):
                    alternate_sequence_bases.append(alt_base_dict[integer])
                else:
                    alternate_sequence_bases.append(ref_base_dict[integer])
            alternate_sequence = ''.join(alternate_sequence_bases)
            alternate_sequence = alternate_sequence.replace('X', '')

            # Add the new alt to the gap fill variant and generate transcript variant
            genomic_gap_fill_variant.posedit.edit.alt = alternate_sequence
            hgvs_refreshed_variant = self.validator.vm.g_to_t(genomic_gap_fill_variant,
                                                              tx_gap_fill_variant.ac,
                                                              alt_aln_method=self.validator.alt_aln_method)

            # Set warning
            gap_size = str(len(genomic_gap_fill_variant.posedit.edit.ref) - 2)
            self.disparity_deletion_in[1] = [gap_size]

        else:

            if self.tx_hgvs_not_delins.posedit.pos.start.offset == 0 and \
                    self.tx_hgvs_not_delins.posedit.pos.end.offset == 0:
                # In this instance, we have identified a transcript gap but the n. version of
                # the transcript variant but do not have a position which actually hits the gap,
                # so the variant likely spans the gap, and is not picked up by an offset.
                try:
                    c1 = self.validator.vm.n_to_c(self.tx_hgvs_not_delins)
                except:
                    c1 = self.tx_hgvs_not_delins
                g3 = self.validator.nr_vm.t_to_g(c1, hgvs_genomic.ac, alt_aln_method=self.validator.alt_aln_method)

                # Check to see if gap is already handled
                hgvs_genomic_norm = self.variant.hn.normalize(hgvs_genomic)
                if (((g3.posedit.pos.end.base - g3.posedit.pos.start.base) >
                        (hgvs_genomic_norm.posedit.pos.end.base - hgvs_genomic_norm.posedit.pos.start.base)) and
                hgvs_genomic_norm.posedit.edit.type == 'del' and
                    g3.posedit.pos.end.base == hgvs_genomic_norm.posedit.pos.end.base):
                    hgvs_refreshed_variant = self.tx_hgvs_not_delins
                    return hgvs_refreshed_variant

                g3.posedit.pos.end.base = g3.posedit.pos.start.base + (len(g3.posedit.edit.ref) - 1)
                try:
                    c2 = self.validator.vm.g_to_t(g3, c1.ac, alt_aln_method=self.validator.alt_aln_method)
                    if c2.posedit.pos.start.offset == 0 and c2.posedit.pos.end.offset == 0:
                        pass
                    else:
                        self.tx_hgvs_not_delins = c2
                        try:
                            self.tx_hgvs_not_delins = self.validator.vm.c_to_n(self.tx_hgvs_not_delins)
                        except vvhgvs.exceptions.HGVSError as e:
                            logger.debug("Except passed, %s", e)
                except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                    logger.debug("Except passed, %s", e)

            if '+' in str(self.tx_hgvs_not_delins.posedit.pos.start) and \
                    '+' not in str(self.tx_hgvs_not_delins.posedit.pos.end):
                hgvs_refreshed_variant = self.c2_pos_edit(hgvs_genomic)

            elif '+' in str(self.tx_hgvs_not_delins.posedit.pos.end) and \
                    '+' not in str(self.tx_hgvs_not_delins.posedit.pos.start):
                hgvs_genomic_norm = self.variant.hn.normalize(hgvs_genomic)
                self.auto_info = self.auto_info
                self.gapped_transcripts = self.gapped_transcripts + ' ' + str(self.tx_hgvs_not_delins.ac)
                hgvs_refreshed_variant = self.c1_pos_edit(hgvs_genomic)

            elif '-' in str(self.tx_hgvs_not_delins.posedit.pos.start) and \
                    '-' not in str(self.tx_hgvs_not_delins.posedit.pos.end):
                hgvs_refreshed_variant = self.c2_pos_edit(hgvs_genomic)

            elif '-' in str(self.tx_hgvs_not_delins.posedit.pos.end) and \
                    '-' not in str(self.tx_hgvs_not_delins.posedit.pos.start):
                self.auto_info = self.auto_info #
                self.gapped_transcripts = self.gapped_transcripts + ' ' + str(self.tx_hgvs_not_delins.ac)

                # Have variation in first copy here!
                if running_option == 1:
                    try:
                        c1 = self.validator.vm.n_to_c(self.tx_hgvs_not_delins)
                    except:
                        c1 = self.tx_hgvs_not_delins
                    c2 = copy.deepcopy(c1)
                    c2.posedit.pos.start = c1.posedit.pos.end
                    c2.posedit.pos.end.base = c1.posedit.pos.end.base
                    c2.posedit.pos.end.offset = 0
                    c2.posedit.edit.ref = ''
                    c2.posedit.edit.alt = ''
                    g2 = self.validator.vm.t_to_g(c2, self.variant.hgvs_genomic.ac,
                                                  alt_aln_method=self.validator.alt_aln_method)
                    c2 = self.validator.vm.g_to_t(g2, c2.ac, alt_aln_method=self.validator.alt_aln_method)
                    # reference = c1.posedit.edit.ref + c2.posedit.edit.ref[1:]
                    alternate = c1.posedit.edit.alt + c2.posedit.edit.ref[1:]
                    c3 = copy.deepcopy(c1)
                    c3.posedit.pos.end = c2.posedit.pos.end
                    c3.posedit.edit.ref = ''  # reference
                    c3.posedit.edit.alt = alternate
                    hgvs_refreshed_variant = c3
                else:
                    hgvs_refreshed_variant = self.c1_pos_edit(hgvs_genomic)

            else:
                # Have variation in second copy here!
                if running_option == 2:
                    self.tx_hgvs_not_delins.posedit.pos.end.base = self.tx_hgvs_not_delins.posedit.pos.start.base + len(
                        self.tx_hgvs_not_delins.posedit.edit.ref) - 1
                elif running_option != 4:
                    self.gapped_transcripts = self.gapped_transcripts + ' ' + str(self.tx_hgvs_not_delins.ac)

                hgvs_refreshed_variant = self.tx_hgvs_not_delins

        return hgvs_refreshed_variant

    def edit_output(self, hgvs_refreshed_variant, saved_hgvs_coding):
        if 'NM_' in str(hgvs_refreshed_variant.ac) and 'c' not in str(hgvs_refreshed_variant.type):
            hgvs_refreshed_variant = self.variant.evm.n_to_c(hgvs_refreshed_variant)

        try:
            hgvs_refreshed_variant = self.variant.hn.normalize(hgvs_refreshed_variant)
            if hgvs_refreshed_variant.posedit.edit.type == 'delins' and \
                    hgvs_refreshed_variant.posedit.edit.ref[-1] == \
                    hgvs_refreshed_variant.posedit.edit.alt[-1]:
                hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[
                                                          0:-1]
                hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[
                                                          0:-1]
                hgvs_refreshed_variant.posedit.pos.end.base = hgvs_refreshed_variant.posedit.pos.end.base - 1
                hgvs_refreshed_variant = self.variant.hn.normalize(hgvs_refreshed_variant)

            elif hgvs_refreshed_variant.posedit.edit.type == 'delins' and \
                    hgvs_refreshed_variant.posedit.edit.ref[0] == \
                    hgvs_refreshed_variant.posedit.edit.alt[0]:
                hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[
                                                          1:]
                hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[
                                                          1:]
                hgvs_refreshed_variant.posedit.pos.start.base = hgvs_refreshed_variant.posedit.pos.start.base + 1
                hgvs_refreshed_variant = self.variant.hn.normalize(hgvs_refreshed_variant)

        except Exception as e:
            error = str(e)
            # Ensure the final variant is not intronic nor does it cross exon boundaries
            if 'Normalization of intronic variants is not supported' in error or \
                    'Unsupported normalization of variants spanning the exon-intron boundary' in error or \
                    "Unsupported normalization of variants spanning the UTR-exon boundary" in error:
                hgvs_refreshed_variant = saved_hgvs_coding

                if hgvs_refreshed_variant.posedit.edit.type == 'delins' and \
                        hgvs_refreshed_variant.posedit.edit.alt == "":
                    hgvs_refreshed_variant = str(hgvs_refreshed_variant).replace("ins", "")
                    hgvs_refreshed_variant = self.validator.hp.parse(hgvs_refreshed_variant)

        return hgvs_refreshed_variant

    def logic_check(self, hgvs_not_delins, rn_tx_hgvs_not_delins, hgvs_coding, do_continue=False, offset_check=False):
        # Logic
        if hgvs_not_delins.posedit.edit.ref is None:
            hgvs_not_delins.posedit.edit.ref = ''
        if rn_tx_hgvs_not_delins.posedit.edit.ref is None:
            rn_tx_hgvs_not_delins.posedit.edit.ref = ''
        if len(hgvs_not_delins.posedit.edit.ref) < len(rn_tx_hgvs_not_delins.posedit.edit.ref):
            gap_length = len(rn_tx_hgvs_not_delins.posedit.edit.ref) - len(hgvs_not_delins.posedit.edit.ref)
            self.disparity_deletion_in = ['chromosome', gap_length]
        elif len(hgvs_not_delins.posedit.edit.ref) > len(rn_tx_hgvs_not_delins.posedit.edit.ref):
            gap_length = len(hgvs_not_delins.posedit.edit.ref) - len(rn_tx_hgvs_not_delins.posedit.edit.ref)
            self.disparity_deletion_in = ['transcript', gap_length]
        else:
            re_capture_tx_variant = []
            for an_internal_possibility in self.hgvs_genomic_possibilities:

                # Set variables from list formats
                try:
                    internal_possibility = an_internal_possibility[1][3]
                except IndexError:
                    internal_possibility = an_internal_possibility[0]

                # Continue
                if internal_possibility == '':
                    continue
                hgvs_t_possibility = self.validator.vm.g_to_t(internal_possibility,
                                                              hgvs_coding.ac,
                                                              alt_aln_method=self.validator.alt_aln_method)
                if hgvs_t_possibility.posedit.edit.type == 'ins':
                    try:
                        hgvs_t_possibility = self.validator.vm.c_to_n(hgvs_t_possibility)
                    except Exception as e:
                        if do_continue:
                            continue
                        logger.debug("Except passed, %s", e)
                    if offset_check:
                        if hgvs_t_possibility.posedit.pos.start.offset != 0 or \
                                hgvs_t_possibility.posedit.pos.end.offset != 0:
                            continue
                    ins_ref = self.validator.sf.fetch_seq(hgvs_t_possibility.ac,
                                                          hgvs_t_possibility.posedit.pos.start.base - 1,
                                                          hgvs_t_possibility.posedit.pos.start.base + 1)
                    try:
                        hgvs_t_possibility = self.validator.vm.n_to_c(hgvs_t_possibility)
                    except Exception as e:
                        if do_continue:
                            continue
                        logger.debug("Except passed, %s", e)
                    hgvs_t_possibility.posedit.edit.ref = ins_ref
                    hgvs_t_possibility.posedit.edit.alt = ins_ref[
                                                              0] + hgvs_t_possibility.posedit.edit.alt + ins_ref[1]
                if internal_possibility.posedit.edit.type == 'ins':
                    ins_ref = self.validator.sf.fetch_seq(internal_possibility.ac,
                                                          internal_possibility.posedit.pos.start.base - 1,
                                                          internal_possibility.posedit.pos.end.base)
                    internal_possibility.posedit.edit.ref = ins_ref
                    internal_possibility.posedit.edit.alt = ins_ref[
                                                                0] + internal_possibility.posedit.edit.alt + ins_ref[1]

                if len(hgvs_t_possibility.posedit.edit.ref) < len(internal_possibility.posedit.edit.ref):
                    gap_length = len(internal_possibility.posedit.edit.ref) - len(hgvs_t_possibility.posedit.edit.ref)
                    re_capture_tx_variant = ['transcript', gap_length, hgvs_t_possibility]
                    hgvs_not_delins = internal_possibility
                    self.hgvs_genomic_5pr = internal_possibility
                    break

            if re_capture_tx_variant:
                try:
                    self.tx_hgvs_not_delins = self.validator.vm.c_to_n(re_capture_tx_variant[2])
                except:
                    self.tx_hgvs_not_delins = re_capture_tx_variant[2]
                self.disparity_deletion_in = re_capture_tx_variant[0:-1]

        return hgvs_not_delins

    def get_hgvs_seek_var(self, hgvs_genomic, hgvs_coding, ori=None, with_query_genomic=False):
        if not ori:
            ori = self.orientation

        if ori == -1:
            try:
                query_genomic = self.variant.reverse_normalizer.normalize(hgvs_genomic)
            except:
                query_genomic = hgvs_genomic
        else:
            # position genomic at its most 3 prime position
            try:
                query_genomic = self.variant.hn.normalize(hgvs_genomic)
            except:
                query_genomic = hgvs_genomic

        # Normalise intronic, if called with query_genomic
        if with_query_genomic:
            if hgvs_coding.posedit.pos.start.offset != 0:
                try:
                    hgvs_coding = self.variant.evm.g_to_t(query_genomic, hgvs_coding.ac)
                except vvhgvs.exceptions.HGVSInvalidIntervalError:
                    pass

        # Map to the transcript and test for movement
        try:
            hgvs_seek_var = self.variant.evm.g_to_t(query_genomic, hgvs_coding.ac)
        except vvhgvs.exceptions.HGVSError:
            hgvs_seek_var = hgvs_coding

        if with_query_genomic:
            return hgvs_seek_var, query_genomic, hgvs_coding

        return hgvs_seek_var

    def rev_norm_ins(self, hgvs_coding, hgvs_genomic):
        # direct mapping from reverse_normalized transcript insertions in the delins format
        try:
            if hgvs_coding.posedit.edit.type == 'ins':
                most_5pr_hgvs_transcript_variant = copy.deepcopy(hgvs_coding)
                most_3pr_hgvs_transcript_variant = self.variant.reverse_normalizer.normalize(hgvs_coding)
                try:
                    n_3pr = self.validator.vm.c_to_n(most_3pr_hgvs_transcript_variant)
                    n_5pr = self.validator.vm.c_to_n(most_5pr_hgvs_transcript_variant)
                except:
                    n_3pr = most_3pr_hgvs_transcript_variant
                    n_5pr = most_5pr_hgvs_transcript_variant
                # Make into a delins by adding the ref bases to the variant ref and alt
                pr3_ref = self.validator.sf.fetch_seq(hgvs_coding.ac, n_3pr.posedit.pos.start.base - 1,
                                                      n_3pr.posedit.pos.end.base)
                pr5_ref = self.validator.sf.fetch_seq(hgvs_coding.ac, n_5pr.posedit.pos.start.base - 1,
                                                      n_5pr.posedit.pos.end.base)
                most_3pr_hgvs_transcript_variant.posedit.edit.ref = pr3_ref
                most_5pr_hgvs_transcript_variant.posedit.edit.ref = pr5_ref
                most_3pr_hgvs_transcript_variant.posedit.edit.alt = pr3_ref[0] + \
                    most_3pr_hgvs_transcript_variant.posedit.edit.alt + pr3_ref[1]
                most_5pr_hgvs_transcript_variant.posedit.edit.alt = pr5_ref[0] + \
                    most_5pr_hgvs_transcript_variant.posedit.edit.alt + pr5_ref[1]
                # Map to the genome
                genomic_from_most_3pr_hgvs_transcript_variant = self.validator.vm.t_to_g(
                    most_3pr_hgvs_transcript_variant, hgvs_genomic.ac, alt_aln_method=self.validator.alt_aln_method)
                genomic_from_most_5pr_hgvs_transcript_variant = self.validator.vm.t_to_g(
                    most_5pr_hgvs_transcript_variant, hgvs_genomic.ac, alt_aln_method=self.validator.alt_aln_method)

                # Normalize - If the variant spans a gap it should then form a static genomic variant
                try:
                    genomic_from_most_3pr_hgvs_transcript_variant = self.variant.hn.normalize(
                        genomic_from_most_3pr_hgvs_transcript_variant)
                except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                    error = str(e)
                    if error == 'base start position must be <= end position':
                        start = genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.start.base
                        end = genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.end.base
                        genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.start.base = end
                        genomic_from_most_3pr_hgvs_transcript_variant.posedit.pos.end.base = start
                        genomic_from_most_3pr_hgvs_transcript_variant = self.variant.hn.normalize(
                            genomic_from_most_3pr_hgvs_transcript_variant)
                try:
                    genomic_from_most_5pr_hgvs_transcript_variant = self.variant.hn.normalize(
                        genomic_from_most_5pr_hgvs_transcript_variant)
                except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                    error = str(e)
                    if error == 'base start position must be <= end position':
                        start = genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.start.base
                        end = genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.end.base
                        genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.start.base = end
                        genomic_from_most_5pr_hgvs_transcript_variant.posedit.pos.end.base = start
                        genomic_from_most_5pr_hgvs_transcript_variant = self.variant.hn.normalize(
                            genomic_from_most_5pr_hgvs_transcript_variant)

                try:
                    if genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.alt is None:
                        genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.alt = ''
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        genomic_from_most_3pr_hgvs_transcript_variant_delins_from_dup = fn.hgvs_dup2indel(
                            genomic_from_most_3pr_hgvs_transcript_variant)
                        genomic_from_most_3pr_hgvs_transcript_variant = self.validator.hp.parse_hgvs_variant(
                            genomic_from_most_3pr_hgvs_transcript_variant_delins_from_dup)

                try:
                    if most_3pr_hgvs_transcript_variant.posedit.edit.alt is None:
                        most_3pr_hgvs_transcript_variant.posedit.edit.alt = ''
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        most_3pr_hgvs_transcript_variant_delins_from_dup = fn.hgvs_dup2indel(
                            most_3pr_hgvs_transcript_variant)
                        most_3pr_hgvs_transcript_variant = self.validator.hp.parse_hgvs_variant(
                            most_3pr_hgvs_transcript_variant_delins_from_dup)

                try:
                    if genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.alt is None:
                        genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.alt = ''
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        genomic_from_most_5pr_hgvs_transcript_variant_delins_from_dup = fn.hgvs_dup2indel(
                            genomic_from_most_5pr_hgvs_transcript_variant)
                        genomic_from_most_5pr_hgvs_transcript_variant = self.validator.hp.parse_hgvs_variant(
                            genomic_from_most_5pr_hgvs_transcript_variant_delins_from_dup)

                try:
                    if most_5pr_hgvs_transcript_variant.posedit.edit.alt is None:
                        most_5pr_hgvs_transcript_variant.posedit.edit.alt = ''
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        most_5pr_hgvs_transcript_variant_delins_from_dup = fn.hgvs_dup2indel(
                            most_5pr_hgvs_transcript_variant)
                        most_5pr_hgvs_transcript_variant = self.validator.hp.parse_hgvs_variant(
                            most_5pr_hgvs_transcript_variant_delins_from_dup)

                if len(genomic_from_most_3pr_hgvs_transcript_variant.posedit.edit.alt) < len(
                        most_3pr_hgvs_transcript_variant.posedit.edit.alt):
                    self.hgvs_genomic_possibilities.append([genomic_from_most_3pr_hgvs_transcript_variant,
                                                            ['false', 'false']])
                if len(genomic_from_most_5pr_hgvs_transcript_variant.posedit.edit.alt) < len(
                        most_5pr_hgvs_transcript_variant.posedit.edit.alt):
                    self.hgvs_genomic_possibilities.append([genomic_from_most_5pr_hgvs_transcript_variant,
                                                            ['false', 'false']])

        except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
            logger.debug("Except passed, %s", e)

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
