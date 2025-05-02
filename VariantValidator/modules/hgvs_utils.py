"""
A variety of functions that convert parser hgvs objects into VCF component parts
Each function has a slightly difference emphasis
"""

# Import modules
import re
import copy
from . import seq_data
from . import utils

# Import Biopython modules
from Bio.Seq import Seq
import vvhgvs
import vvhgvs.exceptions
from vvhgvs.location import BaseOffsetInterval, Interval
# used to set coordinate origin point i.e. seq start vs CDS start/end
from vvhgvs.enums import Datum
from vvhgvs.location import AAPosition
from vvhgvs.edit import AASub, AARefAlt, Dup, NARefAlt
from vvhgvs.posedit import PosEdit

# Database connections and hgvs objects are now passed from VariantValidator.py

# Error handling
class PseudoVCF2HGVSError(Exception):
    pass

class VVPosEdit(PosEdit):
    "override class for posedit to get VV specific formatting"
    met_variation = None
    def __init__(
            self,pos,edit,
            uncertain=False,
            nucleotide_not_equal = None,
            met_variation = None):
        PosEdit.__init__(self,pos=pos,edit=edit,uncertain=uncertain)
        # used to append methylation variation to the output
        self.met_variation = met_variation
        # used for formatting prot consequence of nuc input (for Ter/*)
        self.nucleotide_not_equal = nucleotide_not_equal

    def __eq__(self, other):
        "make VVPosEdit == PosEdit when appropriate"
        if not isinstance(other, PosEdit):
            return NotImplemented
        return self.pos == other.pos and self.edit == other.edit and \
                self.uncertain == other.uncertain

    def __hash__(self):
        return hash((self.pos,self.edit,self.uncertain))

    def format(self, conf=None):
        """Formatting the string of PosEdit with vv edits
        Handles brackets for predicted variants slightly differently
        Also accounts for single base =
        """
        edit = str(self.edit.format(conf))
        if self.pos is None:
            return edit
        elif self.pos.start is None and self.pos.end is None and self.uncertain:
            return edit + "?"
        # do Ter<aa NO>Ter/*<aa NO>* as just Ter=/*= and other AA specific edits
        if type(self.pos.start) is AAPosition:
            if (type(self.edit) is AASub and self.pos.start.aa == '*' and self.edit.alt == '*') or (
                type(self.edit) is AARefAlt and self.pos.start.aa == '*' and
                self.edit.ref == self.pos.start.aa and self.edit.ref==self.edit.alt):
                if not self.nucleotide_not_equal:
                    three_base_convert = False
                    if conf and "p_3_letter" in conf and conf["p_3_letter"] is not None:
                        three_base_convert = conf["p_3_letter"]
                        force_Ter_star = False
                    if conf and "p_term_asterisk" in conf and conf["p_term_asterisk"] is not None:
                        force_Ter_star = conf["p_term_asterisk"]
                    if three_base_convert and not force_Ter_star:
                        formatted_str = f"Ter="
                    else:
                        formatted_str = f"*="
                else: # force coordinates on nucleotide change
                    formatted_str = f"{self.pos.format(conf)}="
            else:
                formatted_str = f"{self.pos.format(conf)}{edit}"
        elif type(self.edit) in [NARefAlt, Dup]:
            if self.edit.ref and 'N' in self.edit.ref or \
                    type(self.edit) is NARefAlt and self.edit.alt and 'N' in self.edit.alt:
                edit = str(self.edit.format()) # ignore instruction to remove ref in N case
                formatted_str = f"{self.pos.format(conf)}{edit}"
            elif type(self.edit) is Dup:
                if self.edit.ref and len(self.edit.ref) == 1:
                    # do not add ref for single base dup (also apply to other single base changes?)
                    formatted_str = f"{self.pos.format(conf)}dup"
                else:
                    formatted_str = f"{self.pos.format(conf)}{edit}"
            else:
                 formatted_str = f"{self.pos.format(conf)}{edit}"
            if type(self.edit) is NARefAlt and self.met_variation:
                if self.edit.ref == self.edit.alt:
                    formatted_str = formatted_str[:-1] + self.met_variation
                #Do we want to do more, if met annotation not valid?
        else:
            formatted_str = f"{self.pos.format(conf)}{edit}"

        if self.uncertain:
            if self.edit in ["0", ""]:
                return f"({formatted_str}?)"
            else:
                return f"({formatted_str})"
        return formatted_str

    __str__ = format


def vcfcp_to_hgvsstr(vcf_dict, start_hgvs):
    """
    converts  vcf components to a string hgvs variant
    :param vcf_dict:
    :return: str(hgvs_variant) with no normalization
    """
    pos = int(vcf_dict['pos'])
    ref = vcf_dict['ref']
    alt = vcf_dict['alt']
    # Generate an end position
    end = str(pos + len(ref) - 1)
    str_hgvs = "%s:%s.%s_%sdel%sins%s" % (start_hgvs.ac, start_hgvs.type, pos, end, ref, alt)
    return str_hgvs


def vcfcp_to_hgvs_obj(vcf_dict, start_hgvs):
    """
    Converts updated vcf components and an original hgvs variant into a hgvs
    object.
    params:
     vcf_dict: VCF dict with updated data
     start_hgvs: orig hgvs, or other object with .ac and .type variables
    returns: hgvs variant object (not normalised)
    """
    pos = int(vcf_dict['pos'])
    return vvhgvs.sequencevariant.SequenceVariant(
            ac=start_hgvs.ac,
            type=start_hgvs.type,
            posedit=vvhgvs.posedit.PosEdit(
                vvhgvs.location.Interval(
                    start=vvhgvs.location.SimplePosition(base=pos),
                    end=vvhgvs.location.SimplePosition(base=pos + len(vcf_dict['ref']) - 1),
                    uncertain=start_hgvs.posedit.pos.uncertain
                    ),
                vvhgvs.edit.NARefAlt(ref=vcf_dict['ref'], alt=vcf_dict['alt'])
                )
            )

def unset_hgvs_obj_ref(hgvs):
    """
    Remove/unset ref bases from hgvs object, in the manner needed for output,
    but without re-parsing from text.
    """
    edit = hgvs.posedit.edit
    if edit.type in ['inv', 'dup']:
        edit.ref = ''
    elif edit.ref is not None and edit.alt is not None:
       #if #edit.alt != edit.ref and \
        if len(edit.alt) == 1 and len(edit.ref) == 1:
            pass
        elif edit.alt == edit.ref and len(edit.ref):
            pass # edit.ref = ''
        elif 'N' in edit.ref:
            pass
        else:
            edit.ref = ''
    elif edit.ref is not None:
        edit.ref = ''
    hgvs.posedit.edit = edit
    return hgvs

def hgvs_dup_to_delins(hgvs_dup):
    """
    Simple utility shorthand for hgvs_delins_parts_to_hgvs_obj, substitutes for
    hgvs_dup2indel, but provides hgvs object output when given a hgvs object
    containing a duplication input, as opposed to hgvs_dup_to_delins's text
    .output.
    param:
     hgvs_dup: A hgvs object containg a duplication (this will not be checked)
               Required.
    returns: A hgvs variant object with a delins equivalent to the input.
    """
    return hgvs_delins_parts_to_hgvs_obj(
            hgvs_dup.ac,
            hgvs_dup.type,
            hgvs_dup.posedit.pos,
            hgvs_dup.posedit.edit.ref,
            hgvs_dup.posedit.edit.ref + hgvs_dup.posedit.edit.ref)

def _derive_hgvs_obj_coordinate_origin(starts,ref_type=None,end=None):
    """
    retrun a pair of coordinate_origin types in the from of Datum.SEQ_START
    Datum.CDS_START or Datum.CDS_END, when given input coordinates/type
    also can return a start stop pair truncated ('*' removed) for int
    conversion
    param:
     starts: Start location must be int or str, required.
     ref_type: Type of ref given, if not 'c' return is simplified, optional.
     stop location: Like start, but not required, if missing Datum pair are
                    identical.
    returns: A pair of Datum enum results for the given input and a coordinate
             pair as a list, trimmed if needed for int conversion.
    """
    if ref_type != 'c':
        return Datum.SEQ_START, Datum.SEQ_START, [starts,end]
    loc_start = starts
    loc_end = end
    coordinate_origin = Datum.CDS_START
    if type(starts) is str and starts[0] == '*':
        loc_start = starts[1::]
        coordinate_origin = Datum.CDS_END
    end_coordinate_origin = coordinate_origin
    if end:
        end_coordinate_origin = Datum.CDS_START
        if type(end) is str and end[0] == '*':
            loc_end = end[1::]
            end_coordinate_origin = Datum.CDS_END
    return coordinate_origin,end_coordinate_origin,[loc_start,loc_end]

def _hgvs_offset_pos_from_str_in(starts,length,ref_type=None,end=None):
    """
    Handle offset positions a bit better.
    This does not fully handle complex stop position issues, and thus will
    always return end as as intronic + if input start is intronic +.
    We also assume that if we are given an end specifically but not a span
    then we got non offset coordinates (these do need to be stored as a
    BaseOffsetPosition for validation purposes)
    """
    sep = ''
    # set up the start point for the coordinates if c type input is given
    # also strip '*' if present
    coordinate_origin, end_coordinate_origin, loc = \
            _derive_hgvs_obj_coordinate_origin(
                    starts,ref_type=ref_type,end=end)
    loc_start = loc[0]
    loc_end = loc[1]
    assert type(loc_start) is int or '_' not in loc_start
    # set end pos if given
    end_pos = None
    if end and type(end) is str:
        if '-' in loc_end[1:]:
            prefix, sep, offset = loc_end[1:].partition('-')
            prefix = loc_end[0] + prefix
            end_pos = vvhgvs.location.BaseOffsetPosition(
                    base=int(prefix),
                    offset=-int(offset),
                    datum=end_coordinate_origin)
        elif '+' in loc_end:
            prefix, sep, offset = loc_end.partition('+')
            end_pos = vvhgvs.location.BaseOffsetPosition(
                    base=int(prefix),
                    offset=int(offset),
                    datum=end_coordinate_origin)
        else:
            end_pos = vvhgvs.location.BaseOffsetPosition(
                    base = int(loc_end),
                    datum=end_coordinate_origin)
    elif end:
        end_pos = vvhgvs.location.BaseOffsetPosition(
                base = int(end),
                datum=end_coordinate_origin)
    # set start deriving end if needed
    if type(starts) is str and '-' in starts[1:]:
        prefix, sep, offset = loc_start[1:].partition('-')
        prefix = loc_start[0] + prefix
        end_offset = int(offset) - (length -1)
        start_pos = vvhgvs.location.BaseOffsetPosition(
                base=int(prefix),
                offset=-int(offset),
                datum=coordinate_origin)
        if end_pos:# already set
            pass
        elif end_offset < 0:
            prefix = int(prefix) - end_offset
            end_pos = vvhgvs.location.BaseOffsetPosition(
                    base=int(prefix),
                    offset=0,
                    datum=end_coordinate_origin)
        else:
            end_pos = vvhgvs.location.BaseOffsetPosition(
                    base=int(prefix),
                    offset=-end_offset,
                    datum=end_coordinate_origin)
    elif type(starts) is str and '+' in starts:
        # no simple way to know whether stop exceeds end of exon
        pos = starts
        prefix, sep,offset = loc_start.partition('+')
        start_pos = vvhgvs.location.BaseOffsetPosition(
                base=int(prefix),
                offset=int(offset),
                datum=coordinate_origin)
        if not end_pos:
            end_pos = vvhgvs.location.BaseOffsetPosition(
                    base=int(prefix),
                    offset=int(offset) + length -1,
                    datum=end_coordinate_origin)
    else: # we got int input for start, or similar so do simple position
        start_pos = vvhgvs.location.BaseOffsetPosition(
                base=int(starts),
                datum=coordinate_origin)
        if not end_pos:
            end_pos = vvhgvs.location.BaseOffsetPosition(
                    base=int(starts)+ length -1,
                    datum=end_coordinate_origin)

    return start_pos, end_pos

def to_vv_hgvs(hgvs):
    "Simple recreate as VV PosEdit vis shim on hgvs_obj_from_existing_edit"
    hgvs = hgvs_obj_from_existing_edit(
            hgvs.ac,
            hgvs.type,
            hgvs.posedit.pos,
            hgvs.posedit.edit,
            unc_posedit=hgvs.posedit.uncertain)
    return hgvs

def hgvs_obj_from_existing_edit(ref_ac,ref_type, starts, edit,
                                end=None, offset_pos=False,
                                unc_posedit=None):
    """
    Converts a set of inputs, including a valid edit from an existing hgvs
    object into a new hgvs object
    params:
     ref_ac: ref accession for output hgvs, required!
     ref_type: ref type eg. g or c for hgvs object, required!
     starts: The location where the coordinates for the delins start, or an
             existing span (as a BaseOffsetInterval which is compatible with
             the hgvs object code), required!
     edit: A valid hgvs object edit, required!

     end: The end location, optional, used to avoid recalculating an already
          known end, and may be also used to test predicted end, though this
          requires a later validate. Unused if a span is given for "starts".
     offset_pos: Are the locations simple or do they need to be the more complex
                 BaseOffsetPosition type? Flag, optional. Unused if span given
                 for "starts".
    returns: hgvs variant object (not normalised)
    """
    if isinstance(starts, Interval):
        return vvhgvs.sequencevariant.SequenceVariant(
                ac=ref_ac,
                type=ref_type,
                posedit=VVPosEdit(
                    starts,
                    edit,
                    uncertain=unc_posedit)
                )
    if offset_pos or ref_type in ['c', 'n']:
        # if we got a null ref and no ref end presume ins i.e. bases either side
        # of insertion variation give 2bp len to substitute for now
        if edit.ref is None:
            length = 2
        else:
            length = len(edit.ref)
        start_pos, end_pos = _hgvs_offset_pos_from_str_in(
                starts,length,ref_type=ref_type,end=end)
        return vvhgvs.sequencevariant.SequenceVariant(
                ac=ref_ac,
                type=ref_type,
                posedit=VVPosEdit(
                    vvhgvs.location.BaseOffsetInterval(
                        start=start_pos,
                        end=end_pos),
                    edit,
                    uncertain=unc_posedit
                    )
                )
    if end:
        ends = int(end)
    else:
        ends = int(starts) + len(edit.ref), -1
    return vvhgvs.sequencevariant.SequenceVariant(
            ac=ref_ac,
            type=ref_type,
            posedit=VVPosEdit(
                vvhgvs.location.Interval(
                    start=vvhgvs.location.SimplePosition(base=int(starts)),
                    end=vvhgvs.location.SimplePosition(base=ends),
                    ),
                edit,
                uncertain=unc_posedit
                )
            )


def hgvs_delins_parts_to_hgvs_obj(ref_ac,ref_type, starts, delete, insert,end=None,offset_pos=False):
    """
    Converts a set of inputs, usually partially from a hgvs object but with
    updates into a new hgvs delins object
    params:
     ref_ac: ref accession for output hgvs, required!
     ref_type: ref type eg. g or c for hgvs object, required!
     starts: The location where the coordinates for the delins start, or an
             existing span (as a BaseOffsetInterval which is compatible with
             the hgvs object code), required!
     delete: The reference sequence over the affected span, required!
     insert: The replacement non ref sequence over the affected span, required!

     end: The end location, optional, used to avoid recalculating an already
          known end, and may be also used to test predicted end, though this
          requires a later validate. unused if a span is given for "starts".
     offset_pos: Are the locations simple or do they need to be the more complex
                 BaseOffsetPosition type? Flag, optional. Unused if span given
                 for "starts".
    returns: hgvs variant object (not normalised)
    """
    if type(starts) in [BaseOffsetInterval, Interval]:
        return vvhgvs.sequencevariant.SequenceVariant(
                ac=ref_ac,
                type=ref_type,
                posedit=VVPosEdit(
                    starts,
                    vvhgvs.edit.NARefAlt(ref=delete, alt=insert)
                    )
                )

    if offset_pos or ref_type in ['c', 'n']:
        start_pos, end_pos = _hgvs_offset_pos_from_str_in(starts,len(delete),ref_type=ref_type,end=end)
        return vvhgvs.sequencevariant.SequenceVariant(
                ac=ref_ac,
                type=ref_type,
                posedit=VVPosEdit(
                    vvhgvs.location.BaseOffsetInterval(start=start_pos,end=end_pos),
                    vvhgvs.edit.NARefAlt(ref=delete, alt=insert)
                    )
                )
    pos = int(starts)
    if end:
        ends = int(end)
    else:
        ends = pos + len(delete) - 1
    return vvhgvs.sequencevariant.SequenceVariant(
            ac=ref_ac,
            type=ref_type,
            posedit=vvhgvs.posedit.PosEdit(
                vvhgvs.location.Interval(
                    start=vvhgvs.location.SimplePosition(base=pos),
                    end=vvhgvs.location.SimplePosition(base=ends),
                    ),
                vvhgvs.edit.NARefAlt(ref=delete, alt=insert)
                )
            )


def hgvs_to_delins_hgvs(hgvs_object, hp, hn, allow_fix=False):
    """
    :param hgvs_object: parsed hgvs string
    :param hp: hgvs_parser
    :param hn: hgvs_normalizer (check function for hn vs reverse hn rules)
    :return: hgvs_object in delins format, see if statements for the details
    """
    if hgvs_object.posedit.edit.type == "delins":
        return hgvs_object
    # Duplications (alt = ref + ref)
    if hgvs_object.posedit.edit.type == "dup":
        v_pos = hgvs_object.posedit.pos.start.base
        v_ref = hgvs_object.posedit.edit.ref
        v_alt = v_ref + v_ref

    # Insertions (Generate the ref, then alt = ref[0] + insertion + ref[1]
    if hgvs_object.posedit.edit.type == "ins":

        # Handle incorrectly formatted ins
        ref_not_two = False
        try:
            if (hgvs_object.posedit.pos.end.base - hgvs_object.posedit.pos.start.base) > 1:
                ref_not_two = True
        except vvhgvs.exceptions.HGVSError:
            pass

        alt_bs = hgvs_object.posedit.edit.alt
        hgvs_object.posedit.edit.alt = ""
        hgvs_object.posedit.edit.ref = ""
        hgvs_object = hn.normalize(hgvs_object)

        if ref_not_two is False or allow_fix is False:
            hgvs_object.posedit.edit.alt = \
                hgvs_object.posedit.edit.ref[0] + \
                alt_bs + \
                hgvs_object.posedit.edit.ref[-1]
        else:
            hgvs_object.posedit.edit.alt = \
                hgvs_object.posedit.edit.ref[0] + \
                alt_bs

        # No stringing needed, return directly
        return hgvs_object

    # Deletions (Handles simple conversion by making alt = "")
    if hgvs_object.posedit.edit.type == "del":
        hgvs_object.posedit.edit.alt = ""
        return hgvs_object

    # Create the object directly via vcfcp_to_hgvs_obj
    return vcfcp_to_hgvs_obj({"pos": v_pos, "ref": v_ref, "alt": v_alt}, hgvs_object)

def pvcf_to_hgvs(query, selected_assembly, normalization_direction, reverse_normalizer, validator):
    """
    :param query: pseudo_vcf string
    :param selected_assembly:
    :param normalization_direction: normalization direction an integer, 5 or 3.
    :param reverse_normalizer:
    :param validator:
    :return:
    """
    # Set normalizer
    selected_normalizer = None
    if normalization_direction == 3:
        selected_normalizer = validator.hn
    if normalization_direction == 5:
        selected_normalizer = reverse_normalizer

    # Gel stye pVCF
    query = query.replace(':', '-')
    pre_input = copy.deepcopy(query)
    vcf_elements = pre_input.split('-')

    # VCF type 1
    if re.search(r'-\d+-[GATC]+-[GATC]+', query):
        query = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], vcf_elements[3])
    elif re.search(r'-\d+-[GATC]+-', query):
        query = '%s:%s%s>%s' % (vcf_elements[0], vcf_elements[1], vcf_elements[2], vcf_elements[2])
    else:
        raise PseudoVCF2HGVSError('Unsupported format: VCF specification 4.1 or later')

    # Chr16:2099572TC>T
    try:
        input_list = query.split(':')
        position_and_edit = input_list[1]
        if not re.match(r'N[CGWT]_', query) and not re.match(r'LRG_\d+$', query):
            chr_num = str(input_list[0])
            chr_num = chr_num.upper()
            chr_num = chr_num.strip()
            if re.match(r'CHR', chr_num):
                chr_num = chr_num.replace('CHR', '')
            # Use selected assembly
            accession = seq_data.to_accession(chr_num, selected_assembly)
            if accession is None:
                error = chr_num + ' is not part of genome build ' + selected_assembly + ' or is not supported'
                raise PseudoVCF2HGVSError(error)
        else:
            accession = input_list[0]

        # Assign reference sequence type
        ref_type = ':g.'
        if 'LRG_' in accession:
            accession = validator.db.get_refseq_id_from_lrg_id(accession)

        # Reformat the variant
        query = str(accession) + ref_type + str(position_and_edit)
    except Exception as e:
        error = str(e)
        raise PseudoVCF2HGVSError(error)

    # Find not_sub type in input e.g. GGGG>G
    not_sub = copy.deepcopy(query)
    not_sub_find = re.compile(r"([GATCgatc]+)>([GATCgatc]+)")
    if not_sub_find.search(not_sub):
        try:
            # If the length of either side of the substitution delimer (>) is >1
            matches = not_sub_find.search(not_sub)
            if len(matches.group(1)) > 1 or len(matches.group(2)) > 1 or re.search(
                    r"([GATCgatc]+)>([GATCgatc]+),([GATCgatc]+)", query):
                # Search for and remove range
                range = re.compile(r"([0-9]+)_([0-9]+)")
                if range.search(not_sub):
                    m = not_sub_find.search(not_sub)
                    start = m.group(1)
                    delete = m.group(2)
                    beginning_string, middle_string = not_sub.split(':')
                    middle_string = middle_string.split('_')[0]
                    end_string = start + '>' + delete
                    not_sub = beginning_string + ':' + middle_string + end_string
                # Split description
                ref_ac, _sep, remainder = not_sub.partition(':')
                ref_type, _sep, posedit = remainder.partition('.')
                pos_ref, _sep, insert = posedit.partition('>')
                # Split remainder using matches
                r = re.compile(r"([0-9]+)([GATCgatc]+)")
                try:
                    m = r.search(pos_ref)
                    delete = m.group(2)
                    starts = posedit.split(delete)[0]
                    hgvs_re_try = hgvs_delins_parts_to_hgvs_obj(
                            ref_ac,
                            ref_type,
                            starts, delete[0], insert)
                    hgvs_re_try.posedit.edit.ref = delete
                    start_pos = str(hgvs_re_try.posedit.pos.start)
                    end_pos = None
                    if '-' in start_pos:
                        base, offset = start_pos.split('-')
                        new_offset = 0 - int(offset) + (len(delete))
                        end_pos = base + '-' + str(new_offset)
                    elif '+' in start_pos:
                        base, offset = start_pos.split('+')
                        new_offset = 0 + int(offset) + (len(delete) - 1)
                        end_pos = base + '+' + str(new_offset)
                    else:
                        end_pos = int(start_pos) + (len(delete) - 1)
                except:
                    not_delins = not_sub
                # Parse into hgvs object
                try:
                    hgvs_not_delins = hgvs_delins_parts_to_hgvs_obj(
                            ref_ac, ref_type, start_pos,
                            delete, insert,
                            ends=end_pos)
                except vvhgvs.exceptions.HGVSError as e:
                    # Sort out multiple ALTS from VCF inputs
                    if re.search("([GATCgatc]+)>([GATCgatc]+),([GATCgatc]+)", not_delins):
                        error = 'Multiple ALTs not supported by this function'
                        raise PseudoVCF2HGVSError(error)
                    else:
                        error = str(e)
                        raise PseudoVCF2HGVSError(error)

                # HGVS will deal with the errors
                hgvs_object = hgvs_not_delins
            else:
                # we know that this should be a sub type variant, and so ends with R>A,
                # where R and A is 1 base of ref or alt respectivly (since the second
                # match did not trigger).
                start = position_and_edit[:-4]
                if '_' in position_and_edit[:-3]:
                    start, _sep, end = position_and_edit.partition('_')
                hgvs_object =  hgvs_delins_parts_to_hgvs_obj(
                        str(accession),
                        ref_type,
                        int(start),
                        position_and_edit[-3],
                        position_and_edit[-1])

        except Exception as e:
            error = str(e)
            raise PseudoVCF2HGVSError(error)
    else:
        # we should not get here! if we can we need to handle it
        raise PseudoVCF2HGVSError('Unsupported format: VCF specification 4.1 or later!')

    # Normalize
    hgvs_object = selected_normalizer.normalize(hgvs_object)
    # return
    return hgvs_object


def hgvs2vcf(hgvs_genomic, primary_assembly, reverse_normalizer, sf, extra_flank_bases=0):
    """
    Simple conversion which ensures identity is as 5 prime as possible by adding an extra 5
    prime base. Necessary for most gap handling situations

    :param hgvs_genomic:
    :param primary_assembly:
    :param reverse_normalizer:
    :param sf:
    :return:
    """
    hgvs_genomic_variant = hgvs_genomic
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    if reverse_normalizer is None:
        reverse_normalized_hgvs_genomic = hgvs_genomic_variant
    else:
        reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic_variant)
    # hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

    # Chr
    chr = seq_data.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac, primary_assembly)
    if chr is not None:
        pass
    else:
        chr = reverse_normalized_hgvs_genomic.ac

    # Identity
    if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity':
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'ins':
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
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'sub':
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'del':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        # Recover sequences
        hgvs_del_seq_w_pre_base = sf.fetch_seq(
                str(reverse_normalized_hgvs_genomic.ac),
                adj_start, end)
        pos = str(start)
        ref = hgvs_del_seq_w_pre_base
        alt = hgvs_del_seq_w_pre_base[0]

    # inv
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'inv':
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
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = ins_seq
        if 'inv' in str(reverse_normalized_hgvs_genomic.posedit):
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'delins':
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
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = vcf_del_seq[:1] + ins_seq

    # Duplications
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'dup':
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
            if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity':
                pos = int(pos) - 1
                prev = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), pos - 1, pos)
                pos = str(pos)
                ref = prev + ref
                alt = prev + alt

    # Add flank bases if requested
    if extra_flank_bases > 0:
        original_pos = pos
        pos = str(int(pos) - extra_flank_bases)
        left_flank = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), int(pos) - 1, int(original_pos) - 1)
        right_flank = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), int(original_pos) + len(ref) - 1,
                                   int(original_pos) + len(ref) - 1 + extra_flank_bases)
        ref = left_flank + ref + right_flank
        alt = left_flank + alt + right_flank

    # Dictionary the VCF
    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': reverse_normalized_hgvs_genomic}
    return vcf_dict


def report_hgvs2vcf(hgvs_genomic, primary_assembly, reverse_normalizer, sf):
    """
    Used to report the Most true representation of the VCF i.e. 5 prime normalized but no
    additional bases added. NOTE: no gap handling capabilities

    :param hgvs_genomic:
    :param primary_assembly:
    :param reverse_normalizer:
    :param sf:
    :return:
    """

    import time
    hgvs_genomic_variant = hgvs_genomic

    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic_variant)

    ucsc_pa = ''
    grc_pa = ''
    ucsc_chr = ''
    grc_chr = ''
    chrs = {}
    # Sort the primary assemblies or go through all valid assemblies
    if primary_assembly == 'All':
        # return all valid genome builds on our report output list
        gen_name_map = {
            'GRCh37':'grch37',
            'hg19':'hg19',
            'GRCh38':'grch38',
            'hg38':'hg38'}

        genomes = ['GRCh37','hg19','GRCh38','hg38']
        for genome in genomes:
            if not seq_data.supported_for_mapping(hgvs_genomic_variant.ac, genome):
                continue
            if genome.startswith('GRC'):
                chrom = seq_data.to_chr_num_refseq(
                        reverse_normalized_hgvs_genomic.ac,
                        genome)
            else:
                chrom = seq_data.to_chr_num_ucsc(
                        reverse_normalized_hgvs_genomic.ac,
                        genome)
            if chrom is None:
                chrom = hgvs_genomic_variant.ac
            chrs[gen_name_map[genome]]=chrom
    else:
        if 'GRC' in primary_assembly:
            if '37' in primary_assembly:
                ucsc_pa = 'hg19'
                grc_pa = primary_assembly
            if '38' in primary_assembly:
                ucsc_pa = 'hg38'
                grc_pa = primary_assembly
        else:
            if '19' in primary_assembly:
                ucsc_pa = primary_assembly
                grc_pa = 'GRCh37'
            if '38' in primary_assembly:
                ucsc_pa = primary_assembly
                grc_pa = 'GRCh38'
        # UCSC Chr
        ucsc_chr = seq_data.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac, ucsc_pa)
        if ucsc_chr is not None:
            pass
        else:
            ucsc_chr = reverse_normalized_hgvs_genomic.ac

        # GRC Chr
        grc_chr = seq_data.to_chr_num_refseq(reverse_normalized_hgvs_genomic.ac, grc_pa)
        if grc_chr is not None:
            pass
        else:
            grc_chr = reverse_normalized_hgvs_genomic.ac

    # Identity
    if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity':
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'ins':
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
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'sub':
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'del':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        # Recover sequences
        if adj_start >= 0:
            hgvs_del_seq_w_pre_base = sf.fetch_seq(
                    str(reverse_normalized_hgvs_genomic.ac),
                    adj_start, end)
            ref = hgvs_del_seq_w_pre_base
            alt = hgvs_del_seq_w_pre_base[0]
            pos = str(start)
        else:
            hgvs_del_seq_w_post_base = sf.fetch_seq(
                    str(reverse_normalized_hgvs_genomic.ac),
                    start, end + 1)
            ref = hgvs_del_seq_w_post_base
            alt = hgvs_del_seq_w_post_base[-1]
            pos = "1"

    # inv
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'inv':
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
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        # pos = str(start-1)
        # ref = bs + vcf_del_seq
        ref = vcf_del_seq
        alt = ins_seq
        if reverse_normalized_hgvs_genomic.posedit.edit.type == 'inv':
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'delins':
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
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        # pos = str(start)
        # ref = vcf_del_seq
        # alt = vcf_del_seq[:1] + ins_seq
        pos = str(start + 1)
        ref = vcf_del_seq[1:]
        alt = ins_seq

    # Duplications
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'dup':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)  #
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2  #
        start = start - 1  #
        # Recover sequences
        vcf_ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_ref_seq[0]
        alt = vcf_ref_seq
    else:
        ref = ''
        alt = ''
        pos = ''

    # Dictionary the VCF
    vcf_dict = {'pos': str(pos), 'ref': ref, 'alt': alt, 'ucsc_chr': ucsc_chr, 'grc_chr': grc_chr,
                'normalized_hgvs': reverse_normalized_hgvs_genomic,'chrs_by_genome':chrs}
    return vcf_dict


def pos_lock_hgvs2vcf(hgvs_genomic, primary_assembly, reverse_normalizer, sf):
    """
    No normalization at all. No additional bases added. Simply returns an in-situ VCF
    :param hgvs_genomic:
    :param primary_assembly:
    :param reverse_normalizer:
    :param sf:
    :return:
    """
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
    chr = seq_data.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac, primary_assembly)
    if chr is not None:
        pass
    else:
        chr = reverse_normalized_hgvs_genomic.ac

    # Identity
    if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity':
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'ins':
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
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'sub':
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'del':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        # Recover sequences
        hgvs_del_seq_w_pre_base = sf.fetch_seq(
                str(reverse_normalized_hgvs_genomic.ac),
                adj_start, end)
        # Assemble
        pos = str(start)
        ref = hgvs_del_seq_w_pre_base
        alt = hgvs_del_seq_w_pre_base[0]

    # inv
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'inv':
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
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        # pos = str(start-1)
        # ref = bs + vcf_del_seq
        ref = vcf_del_seq
        alt = ins_seq
        if 'inv' in str(reverse_normalized_hgvs_genomic.posedit):
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'delins':
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
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = vcf_del_seq[:1] + ins_seq

    # Duplications
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'dup':
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

def pre_push_vcf_tx_g_map_fix(norm_hgvs_transcript, un_norm_hgvs, genomic_ac,var_mapper,normaliser,seq_fetcher,mapped_hgvs):
    """
    Pre-fix transcript inputs for the left/right pushing vcf functions, to make
    sure they map correctly. Requires n type pre-normalised/reverse normalised,
    input.
    :pram hgvs_transcript: transcript type variant to check for mapping
    :pram un_norm_hgvs: un-normalised version of the transcript type variant
                        needed to check for del status (???->delins is needed
                        for some, non del, post normalisation, descriptions)
    :pram genomic_ac: genomic ac for mapping to be checked
    :pram mapped_hgvs: mapped_hgvs a genomic mapping derived from
                       hgvs_transcript, to allow shared mapping step when
                       rev-norm == norm
    :pram var_mapper: mapper for n->g mapping unused if pre-mapped
    :pram normaliser: hgvs variant normaliser
    :pram seq_fetcher: seq fetcher from hgvs (currently a SeqRepo wrapper)
    :pram mapped_hgvs: mapped variant (if present we do not use the mapper)
    """


    # Variants in/on a genomic gap that cause issues need sorting by making them span so coordinates do not reverse
    if mapped_hgvs:
        hgvs_genomic_g = mapped_hgvs
    else:
        hgvs_genomic_g = var_mapper.n_to_g(norm_hgvs_transcript, genomic_ac)

    try:
        normaliser.normalize(hgvs_genomic_g)
    except vvhgvs.exceptions.HGVSInvalidVariantError as e:
        if "base start position must be <= end position" in str(e):
            hgvs_genomic_g_identity = copy.deepcopy(hgvs_genomic_g)
            # First, get ins and del into delins
            if hgvs_genomic_g_identity.posedit.edit.type == "dup":
                # Duplications have no "alt" in the object structure so convert to delins
                stb = hgvs_genomic_g_identity.posedit.pos.end.base
                edb = hgvs_genomic_g_identity.posedit.pos.start.base
                hgvs_genomic_g_identity.posedit.pos.start.base = stb
                hgvs_genomic_g_identity.posedit.pos.end.base = edb
                hgvs_genomic_g_identity.posedit.edit.ref = seq_fetcher.fetch_seq(hgvs_genomic_g_identity.ac, stb - 1, edb)
                hgvs_genomic_g_identity = hgvs_to_delins_hgvs(hgvs_genomic_g_identity, None, normaliser)
            else:
                adjust_s = hgvs_genomic_g_identity.posedit.pos.end.base
                adjust_e = hgvs_genomic_g_identity.posedit.pos.start.base
                hgvs_genomic_g_identity.posedit.pos.start.base = adjust_s
                hgvs_genomic_g_identity.posedit.pos.end.base = adjust_e

            hgvs_genomic_g_identity.posedit.edit.ref = ''
            hgvs_genomic_g_identity.posedit.edit.alt = ''
            hgvs_genomic_g_identity = normaliser.normalize(hgvs_genomic_g_identity)

            hgvs_genomic_n_gap = var_mapper.g_to_n(hgvs_genomic_g_identity, norm_hgvs_transcript.ac)

            hgvs_genomic_n_identity = copy.deepcopy(hgvs_genomic_n_gap)
            hgvs_genomic_n_identity.posedit.edit.ref = ""
            hgvs_genomic_n_identity.posedit.edit.alt = ""
            hgvs_genomic_n_identity = normaliser.normalize(hgvs_genomic_n_identity)
            """
            At this stage we have the gap position at g. in hgvs_genomic_g_identity
            # The impact in hgvs_genomic_n_gap
            # And the range the variant needs to span to cover the gap in hgvs_genomic_n_identity
            """
            # First, get ins and del into delins this has to also check the un-normalised version to avoid
            # problems on right push
            if un_norm_hgvs.posedit.edit.type == "del" or norm_hgvs_transcript.posedit.edit.type == "del":
                norm_hgvs_transcript = hgvs_to_delins_hgvs(norm_hgvs_transcript, None, normaliser)
            elif norm_hgvs_transcript.posedit.edit.type == "ins":
                norm_hgvs_transcript = hgvs_to_delins_hgvs(norm_hgvs_transcript, None, normaliser)
            elif norm_hgvs_transcript.posedit.edit.type == "dup":
                norm_hgvs_transcript = hgvs_to_delins_hgvs(norm_hgvs_transcript, None, normaliser)

            # Now create the variant
            if hgvs_genomic_n_identity.posedit.pos.start.base < \
                    norm_hgvs_transcript.posedit.pos.start.base:
                v1 = copy.deepcopy(hgvs_genomic_n_identity)
                v1.posedit.pos.end.base = norm_hgvs_transcript.posedit.pos.start.base - 1
                v1.posedit.edit.ref = ""
                v1.posedit.edit.alt = ""
                v1 = normaliser.normalize(v1)
                if (hgvs_genomic_n_identity.posedit.pos.end.base
                        > norm_hgvs_transcript.posedit.pos.end.base):
                    v3 = copy.deepcopy(hgvs_genomic_n_identity)
                    v3.posedit.pos.start.base = norm_hgvs_transcript.posedit.pos.end.base + 1
                    v3.posedit.edit.ref = ""
                    v3.posedit.edit.alt = ""
                    v3 = normaliser.normalize(v3)
                else:
                    v3 = False

                # Assemble
                hgvs_genomic_n_assembled = copy.deepcopy(hgvs_genomic_n_identity)
                if v3 is not False:
                    hgvs_genomic_n_assembled.posedit.pos.end.base = v3.posedit.pos.end.base

                    ass_ref = (v1.posedit.edit.ref +
                               norm_hgvs_transcript.posedit.edit.ref +
                               v3.posedit.edit.ref)
                    ass_alt = (v1.posedit.edit.alt +
                               norm_hgvs_transcript.posedit.edit.alt +
                               v3.posedit.edit.alt)
                else:
                    hgvs_genomic_n_assembled.posedit.pos.end.base = \
                        norm_hgvs_transcript.posedit.pos.end.base
                    ass_ref = (v1.posedit.edit.ref +
                               norm_hgvs_transcript.posedit.edit.ref)
                    ass_alt = (v1.posedit.edit.alt +
                               norm_hgvs_transcript.posedit.edit.alt)
                hgvs_genomic_n_assembled.posedit.edit.ref = ass_ref
                hgvs_genomic_n_assembled.posedit.edit.alt = ass_alt
                norm_hgvs_transcript = copy.deepcopy(hgvs_genomic_n_assembled)
            else:
                v3 = copy.deepcopy(hgvs_genomic_n_identity)
                v3.posedit.pos.end.base = norm_hgvs_transcript.posedit.pos.start.base - 1
                v3.posedit.edit.ref = ""
                v3.posedit.edit.alt = ""
                if (hgvs_genomic_n_identity.posedit.pos.end.base
                        > norm_hgvs_transcript.posedit.pos.end.base):
                    v3 = copy.deepcopy(hgvs_genomic_n_identity)
                    v3.posedit.pos.start.base = norm_hgvs_transcript.posedit.pos.end.base + 1
                    v3.posedit.edit.ref = ""
                    v3.posedit.edit.alt = ""
                    v3 = normaliser.normalize(v3)
                else:
                    v3 = False

                # Assemble
                hgvs_genomic_n_assembled = copy.deepcopy(hgvs_genomic_n_identity)
                if v3 is not False:
                    hgvs_genomic_n_assembled.posedit.pos.end.base = v3.posedit.pos.end.base
                    ass_ref = (norm_hgvs_transcript.posedit.edit.ref +
                               v3.posedit.edit.ref)
                    ass_alt = (norm_hgvs_transcript.posedit.edit.alt +
                               v3.posedit.edit.alt)
                else:
                    hgvs_genomic_n_assembled.posedit.pos.end.base = \
                        norm_hgvs_transcript.posedit.pos.end.base
                    ass_ref = norm_hgvs_transcript.posedit.edit.ref
                    ass_alt = norm_hgvs_transcript.posedit.edit.alt
                hgvs_genomic_n_assembled.posedit.edit.ref = ass_ref
                hgvs_genomic_n_assembled.posedit.edit.alt = ass_alt
                norm_hgvs_transcript = copy.deepcopy(hgvs_genomic_n_assembled)
    return norm_hgvs_transcript



def hard_right_hgvs2vcf(hgvs_genomic, primary_assembly, hn, reverse_normalizer, sf, tx_ac, hdp, alt_aln_method, hp, vm,
                        mrg, genomic_ac=False, mapped_g=False, pre_norm=False):
    """
    Designed specifically for gap handling.
    hard right pushes as 3 prime as possible and adds additional bases
    :param hgvs_genomic:
    :param primary_assembly:
    :param hn:
    :param reverse_normalizer:
    :param sf:
    :param tx_ac: Transcipt ac when genomic var is input
    :param hdp:
    :param alt_aln_method:
    :param hp:
    :param vm:
    :param genomic_ac: Genomic ac when transcirpt var is input *Must* be false for genomic var
    :return:
    """
    if hgvs_genomic.type == 'g':
        # Reverse normalize input prior to convert: NOTE will replace ref
        if pre_norm:
            normalized_hgvs_genomic = pre_norm
        else:
            normalized_hgvs_genomic = hn.normalize(hgvs_genomic)
    else:
        # c. must be in n. format
        if hgvs_genomic.type == 'c':
            hgvs_genomic = vm.c_to_n(hgvs_genomic)
        if pre_norm:
            normalized_hgvs_genomic = pre_norm
        else:
            normalized_hgvs_genomic = hn.normalize(hgvs_genomic)
        normalized_hgvs_genomic = pre_push_vcf_tx_g_map_fix(
                normalized_hgvs_genomic,
                hgvs_genomic,
                genomic_ac,
                vm,
                hn,# normaliser
                sf,# seq_fetcher
                mapped_g)

    # Chr
    if hgvs_genomic.type == 'g':
        chr = seq_data.to_chr_num_ucsc(normalized_hgvs_genomic.ac, primary_assembly)
        if chr is None:
           chr = normalized_hgvs_genomic.ac
    else:
        chr = normalized_hgvs_genomic.ac

    # identity
    if normalized_hgvs_genomic.posedit.edit.type == 'identity':
        pos = str(normalized_hgvs_genomic.posedit.pos.start)
        ref = normalized_hgvs_genomic.posedit.edit.ref
        alt = normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif normalized_hgvs_genomic.posedit.edit.type == 'ins':
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
    elif normalized_hgvs_genomic.posedit.edit.type == 'sub':
        ref = normalized_hgvs_genomic.posedit.edit.ref
        alt = normalized_hgvs_genomic.posedit.edit.alt
        pos = str(normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif normalized_hgvs_genomic.posedit.edit.type == 'del':
        end = int(normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        # Recover sequences
        hgvs_del_seq_w_adj_start = sf.fetch_seq(str(normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = hgvs_del_seq_w_adj_start
        alt = hgvs_del_seq_w_adj_start[0]

    # inv
    elif normalized_hgvs_genomic.posedit.edit.type == 'inv':
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
        vcf_del_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        # pos = str(start-1)
        # ref = bs + vcf_del_seq
        ref = vcf_del_seq
        alt = ins_seq
        if normalized_hgvs_genomic.posedit.edit.type == 'inv':
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif normalized_hgvs_genomic.posedit.edit.type == 'delins':
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
        vcf_del_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = vcf_del_seq[:1] + ins_seq

    # Duplications
    elif normalized_hgvs_genomic.posedit.edit.type == 'dup':
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
    # If possible, capture and alt variant that spans the gap
    merged_variant = False
    pre_merged_variant = False
    identifying_variant = False
    identifying_g_variant = False
    needs_a_push = False
    if chr != '' and pos != '' and ref != '' and alt != '':

        # Set exon boundary
        if genomic_ac is False:
            # Find the boundaries at the genomic level for the current exon
            exon_set = hdp.get_tx_exons(tx_ac, hgvs_genomic.ac, alt_aln_method)
            exon_end_genomic = None
            for exon in exon_set:
                if int(exon[7]) + 1 <= int(pos) <= int(exon[8]):
                    exon_end_genomic = int(exon[8])
                    break
        else:
            # Trick the system using transcript positions
            exon_set = hdp.get_tx_exons(hgvs_genomic.ac, genomic_ac, alt_aln_method)
            exon_end_genomic = None
            for exon in exon_set:
                if int(exon[5]) + 1 <= int(pos) <= int(exon[6]):
                    exon_end_genomic = int(exon[6])
                    break

        # Set loop variables for extending the push
        push_ref = ref
        push_alt = alt
        working_pos = int(pos) + len(ref)
        if genomic_ac is False:
            genomic_ac = hgvs_genomic.ac
        # Clear staging_loop
        staging_loop = 0

        # Loop and add bases - up to the range defined below - unless we go into an intron/past the transcript
        max_push_length = 50
        try:
            flank_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac), working_pos - 1, working_pos + max_push_length)
        except  vvhgvs.exceptions.HGVSDataNotAvailableError as e:
            if "ValueError: stop out of range" in str(e):
                flank_seq = False
                # this means that we went beyond the end of the seq this should be very rare but is possible
                # fall back to old behaviour here
            else:
                raise e
        for push in range(max_push_length):
            if flank_seq:
                push_ref = push_ref + flank_seq[push]
                push_alt = push_alt + flank_seq[push]
            else:
                post = sf.fetch_seq(str(normalized_hgvs_genomic.ac), working_pos - 1, working_pos)
                push_ref = push_ref + post
                push_alt = push_alt + post

            # Create a not_delins for normalisation checking
            offset_pos = True
            if hgvs_genomic.type in ['g','m']:
                offset_pos=False
            normlize_check_variant = hgvs_delins_parts_to_hgvs_obj(
                    hgvs_genomic.ac,
                    hgvs_genomic.type,
                    int(pos),
                    push_ref,
                    push_alt,
                    end=working_pos,
                    offset_pos=offset_pos)

            # Check to see of we end up spanning a gap
            try:
                if hgvs_genomic.type != "g":
                    normlize_check_mapped = vm.n_to_g(normlize_check_variant, genomic_ac)
                else:
                    normlize_check_mapped = vm.g_to_n(normlize_check_variant, tx_ac, alt_aln_method)

            # Catch out-of-bounds errors
            except vvhgvs.exceptions.HGVSInvalidIntervalError:
                needs_a_push = False
                break

            """
            Break out from loop parameters
            """
            if normlize_check_mapped.posedit.pos.start.base > normlize_check_mapped.posedit.pos.end.base:
                needs_a_push = False
                break
            if len(normlize_check_mapped.posedit.edit.ref) <= 1:
                staging_loop = staging_loop + 1

            # Check here for the gap (Has it been crossed?) Note: if gap in tx, we have the whole gap spanned
            if (((len(normlize_check_mapped.posedit.edit.ref) != len(normlize_check_variant.posedit.edit.ref) and
                  len(normlize_check_mapped.posedit.edit.ref) > 1))
                    or
                    (normlize_check_variant.posedit.edit.type == 'identity')
                    and len(normlize_check_mapped.posedit.edit.alt) != len(normlize_check_variant.posedit.edit.ref)):

                # Add the identifying variant
                identifying_variant = normlize_check_variant

                if push == 0:  # Already crossing the gap so return original vcf
                    end_seq_check_variant = copy.copy(normlize_check_variant)
                    # end_seq_check_variant.posedit.edit.alt = end_seq_check_variant.posedit.edit.ref
                else:
                    # Look to see if the gap has been identified by addition of bases in sequence
                    end_seq_check_variant = hgvs_delins_parts_to_hgvs_obj(
                            hgvs_genomic.ac,
                            hgvs_genomic.type,
                            normlize_check_variant.posedit.pos.end.base - 1 - staging_loop,
                            push_ref[-2 - staging_loop:],
                            push_ref[-2 - staging_loop],
                            end=normlize_check_variant.posedit.pos.end.base,
                            offset_pos=True)

                # Check to see of we end up spanning a gap at the last 2 bases
                if hgvs_genomic.type != "g":
                    end_seq_check_mapped = vm.n_to_g(end_seq_check_variant, genomic_ac)
                else:
                    end_seq_check_mapped = vm.g_to_n(end_seq_check_variant, tx_ac)

                # Look for flank subs that may be missed when naieve mapping c > c made a delins from a sub
                # This is a hgvs.py quirl for flanking subs in the antisense oriemntation and refers to
                # https://github.com/openvar/variantValidator/issues/651
                try:
                    normalized_end_seq_check_mapped = hn.normalize(end_seq_check_mapped)
                    normalized_end_seq_check_variant = hn.normalize(end_seq_check_variant)
                    if (normalized_end_seq_check_mapped.type == 'g' and
                            normalized_end_seq_check_variant.type == 'n' and
                            normalized_end_seq_check_mapped.posedit.edit.type == 'sub' and
                            normalized_end_seq_check_variant == normalized_hgvs_genomic):

                        # double check the original mapping
                        sub_map = vm.g_to_t(normalized_end_seq_check_mapped,  normalized_end_seq_check_variant.ac)
                        if sub_map.posedit.edit.type == 'sub':
                            needs_a_push = True
                            merged_variant = sub_map
                            identifying_g_variant = end_seq_check_mapped
                            break
                except vvhgvs.exceptions.HGVSUnsupportedOperationError:
                    pass

                # For genomic_variant mapped onto gaps, we end up with an offset
                start_offset = False
                end_offset = False
                try:
                    end_seq_check_mapped.posedit.pos.start.offset
                except AttributeError:
                    start_offset = False
                else:
                    if end_seq_check_mapped.posedit.pos.start.offset != 0:
                        start_offset = True
                try:
                    end_seq_check_mapped.posedit.pos.end.offset
                except AttributeError:
                    end_offset = False
                else:
                    if end_seq_check_mapped.posedit.pos.end.offset != 0:
                        end_offset = True
                if start_offset is True or end_offset is True:

                    # To identify the gap, we need to span it before mapping back
                    if end_offset is True:
                        end_seq_check_mapped.posedit.pos.end.base = end_seq_check_mapped.posedit.pos.start.base + 1
                        end_seq_check_mapped.posedit.pos.end.offset = 0
                        end_seq_check_mapped.posedit.edit.ref = ''
                        norml_end_seq_check_mapped = end_seq_check_mapped
                    elif start_offset is True:
                        end_seq_check_mapped.posedit.pos.start.base = end_seq_check_mapped.posedit.pos.end.base - 1
                        end_seq_check_mapped.posedit.pos.start.offset = 0
                        end_seq_check_mapped.posedit.edit.ref = ''
                        norml_end_seq_check_mapped = end_seq_check_mapped

                    # now map back onto original reference sequence
                    try:
                        norml_end_seq_check_mapped = vm.c_to_n(norml_end_seq_check_mapped)  # Need in n. context
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass
                    if hgvs_genomic.type == "g":
                        map_back = vm.n_to_g(norml_end_seq_check_mapped, genomic_ac)
                    else:
                        map_back = vm.g_to_n(norml_end_seq_check_mapped, tx_ac)

                    # Normalize variants, original and the gap induced variant (note, variant pre-normalized)
                    map_back = hn.normalize(map_back)  # gap is left so normalize right
                    map_back_rn = reverse_normalizer.normalize(map_back)
                    try:
                        map_back = vm.c_to_n(map_back)  # Need in n. context
                        map_back_rn = vm.c_to_n(map_back_rn)
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass

                    # Can the variants be normalized together
                    if ((
                            (map_back.posedit.pos.end.base >=
                             normalized_hgvs_genomic.posedit.pos.start.base - 1)
                            and
                            (map_back.posedit.pos.end.base <=
                             normalized_hgvs_genomic.posedit.pos.end.base + 1)
                    )
                            or
                            (
                                    (map_back_rn.posedit.pos.end.base >=
                                     normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.end.base <=
                                     normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back.posedit.pos.start.base >=
                                     normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back.posedit.pos.start.base <=
                                     normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back_rn.posedit.pos.start.base >=
                                     normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.start.base <=
                                     normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )):

                        # Create a variant that reflects the impact of the gap.
                        # This uses variant merging
                        # We merge the "gap" variant and the variant itself
                        v1 = hgvs_genomic
                        v2 = map_back
                        if v2.posedit.edit.type == "identity":
                            needs_a_push = True  # Return new vcf only
                            break
                        if "g" not in hgvs_genomic.type:
                            v1 = vm.n_to_g(hgvs_genomic, genomic_ac)
                            v2 = vm.n_to_g(map_back, genomic_ac)
                        try:
                            v1 = hn.normalize(v1)
                            v2 = hn.normalize(v2)
                        except vvhgvs.exceptions.HGVSInvalidVariantError:
                            needs_a_push = True  # Return new vcf only
                            break
                        else:
                            try:
                                if v1.posedit.pos.start.base < v2.posedit.pos.start.base:
                                    pre_merged_variant = mrg([v1, v2], reverse_normalizer, final_norm=False)
                                else:
                                    pre_merged_variant = mrg([v2, v1], reverse_normalizer, final_norm=False)
                                if "g" in pre_merged_variant.type:
                                    merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                else:
                                    merged_variant = pre_merged_variant
                            except utils.mergeHGVSerror as e:
                                needs_a_push = True  # Return new vcf only
                                break
                            except vvhgvs.exceptions.HGVSParseError:
                                needs_a_push = True  # Return new vcf only
                                break

                            # Ensure merged variant is not in a "non-intron" if mapped back to n.
                            if merged_variant is not False:
                                try:
                                    if (merged_variant.posedit.pos.start.offset != 0
                                            or merged_variant.posedit.pos.start.offset != 0):
                                        # Try from normalized genomic
                                        pre_merged_variant = hn.normalize(pre_merged_variant)
                                        test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                        if (test_merged_variant.posedit.pos.start.offset == 0
                                                and test_merged_variant.posedit.pos.start.offset == 0):
                                            merged_variant = pre_merged_variant
                                        else:
                                            pre_merged_variant = reverse_normalizer.normalize(pre_merged_variant)
                                            test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                            if (test_merged_variant.posedit.pos.start.offset == 0
                                                    and test_merged_variant.posedit.pos.start.offset == 0):
                                                merged_variant = pre_merged_variant
                                    # Map back to n.
                                    if "g" in merged_variant.type:
                                        merged_variant = vm.g_to_n(merged_variant, tx_ac)
                                except AttributeError:
                                    pass
                            needs_a_push = True  # Keep the new vcf
                            break
                    else:
                        needs_a_push = False  # Restore old vcf
                        break

                # Or we have identified the gap again at the expected position
                if len(end_seq_check_mapped.posedit.edit.ref) != len(end_seq_check_variant.posedit.edit.ref):

                    """
                    At this stage, we have done the following, illustrated by a  gap in transcript

                    g. NNNNNNNNNN
                    n. NNNNNNN--N

                    We forced the gap to be projected by making the end_seq_check_variant n.=

                             NN  Deletion in g.
                             |
                    g. NNNNNNNN
                    n. NNNNNNNN

                    So we need to make the g. == again before mapping back, which will make an ins in the n.
                    """

                    # Now normalize the variants to see if they meet
                    norml_end_seq_check_mapped = copy.deepcopy(end_seq_check_mapped)
                    norml_end_seq_check_mapped.posedit.edit.alt = norml_end_seq_check_mapped.posedit.edit.ref

                    # now map back onto original reference sequence
                    try:
                        norml_end_seq_check_mapped = vm.c_to_n(norml_end_seq_check_mapped)  # Need in n. context
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass
                    if hgvs_genomic.type == "g":
                        map_back = vm.n_to_g(norml_end_seq_check_mapped, genomic_ac)
                    else:
                        map_back = vm.g_to_n(norml_end_seq_check_mapped, tx_ac)

                    # In transcript gaps, this can push us fully into the gap
                    try:
                        if map_back.posedit.pos.start.offset != 0 and map_back.posedit.pos.start.offset != 0:
                            needs_a_push = False
                            break
                    except AttributeError:
                        pass

                    # Normalize variants, original and the gap induced variant (note, variant pre-normalized)
                    try:
                        map_back = hn.normalize(map_back)  # gap is left so normalize right
                        map_back_rn = reverse_normalizer.normalize(map_back)
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        needs_a_push = False  # Restore old vcf
                        break
                    try:
                        map_back = vm.c_to_n(map_back)  # Need in n. context
                        map_back_rn = vm.c_to_n(map_back_rn)
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass

                    # Is the gap variant the same as the incoming variant?
                    if normalized_hgvs_genomic == map_back:
                        needs_a_push = True
                        push_ref = end_seq_check_variant.posedit.edit.ref
                        push_alt = end_seq_check_variant.posedit.edit.alt
                        pos = end_seq_check_variant.posedit.pos.start.base
                        break

                    # Can the variants be normalized together
                    if ((
                            (map_back.posedit.pos.end.base >=
                             normalized_hgvs_genomic.posedit.pos.start.base - 1)
                            and
                            (map_back.posedit.pos.end.base <=
                             normalized_hgvs_genomic.posedit.pos.end.base + 1)
                    )
                            or
                            (
                                    (map_back_rn.posedit.pos.end.base >=
                                     normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.end.base <=
                                     normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back.posedit.pos.start.base >=
                                     normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back.posedit.pos.start.base <=
                                     normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back_rn.posedit.pos.start.base >=
                                     normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.start.base <=
                                     normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )):

                        # Create a variant that reflects the impact of the gap.
                        # This uses variant merging
                        # We merge the "gap" variant and the variant itself
                        v1 = hgvs_genomic
                        v2 = map_back

                        if v2.posedit.edit.type == "identity":
                            needs_a_push = True  # Return new vcf only
                            break
                        if "g" not in hgvs_genomic.type:
                            v1 = vm.n_to_g(hgvs_genomic, genomic_ac)
                            v2 = vm.n_to_g(map_back, genomic_ac)

                        # Known examples of incorrect formatting from vm
                        ################################################

                        # 1. vm causes an insertion length of > 1 because of the gap - issue #392
                        if "ins" in v1.posedit.edit.type and "sub" in v2.posedit.edit.type:
                            try:
                                v1 = hn.normalize(v1)
                            except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                                if "insertion length must be 1" in str(e):
                                    v1 = hgvs_to_delins_hgvs(v1, hp, hn, allow_fix=True)
                                    identifying_g_variant = v1

                        elif "ins" in v2.posedit.edit.type and "sub" in v1.posedit.edit.type:
                            try:
                                v2 = hn.normalize(v2)
                            except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                                if "insertion length must be 1" in str(e):
                                    v2 = hgvs_to_delins_hgvs(v2, hp, hn, allow_fix=True)
                                    identifying_g_variant = v2
                        try:
                            v1 = hn.normalize(v1)
                            v2 = hn.normalize(v2)
                        except vvhgvs.exceptions.HGVSInvalidVariantError:
                            needs_a_push = True  # Return new vcf only
                            break
                        else:
                            try:
                                if v1.posedit.pos.start.base < v2.posedit.pos.start.base:
                                    pre_merged_variant = mrg([v1, v2], reverse_normalizer, final_norm=False)
                                else:
                                    pre_merged_variant = mrg([v2, v1], reverse_normalizer, final_norm=False)
                                if "g" in pre_merged_variant.type:
                                    merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                else:
                                    merged_variant = pre_merged_variant
                            except utils.mergeHGVSerror:
                                try:
                                    if v1.posedit.pos.start.base < v2.posedit.pos.start.base:
                                        pre_merged_variant = mrg([v1, v2], hn, final_norm=False)
                                    else:
                                        pre_merged_variant = mrg([v2, v1], hn, final_norm=False)
                                    if "g" in pre_merged_variant.type:
                                        merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                    else:
                                        merged_variant = pre_merged_variant
                                except utils.mergeHGVSerror:
                                    needs_a_push = False  # Return new vcf only
                                    break
                                except vvhgvs.exceptions.HGVSParseError:
                                    needs_a_push = False  # Return new vcf only
                                    break

                            except vvhgvs.exceptions.HGVSParseError:
                                needs_a_push = False  # Return new vcf only
                                break

                            # Ensure merged variant is not in a "non-intron" if mapped back to n.
                            if merged_variant is not False:
                                try:
                                    if (merged_variant.posedit.pos.start.offset != 0
                                            or merged_variant.posedit.pos.end.offset != 0):
                                        # Try from normalized genomic
                                        try:
                                            pre_merged_variant = hn.normalize(pre_merged_variant)
                                        except vvhgvs.exceptions.HGVSError:
                                            pass
                                        test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)

                                        if (test_merged_variant.posedit.pos.start.offset == 0
                                                and test_merged_variant.posedit.pos.start.offset == 0):
                                            merged_variant = test_merged_variant
                                        else:
                                            pre_merged_variant = reverse_normalizer.normalize(pre_merged_variant)
                                            test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                            if (test_merged_variant.posedit.pos.start.offset == 0
                                                    and test_merged_variant.posedit.pos.start.offset == 0):
                                                merged_variant = pre_merged_variant
                                    # Map back to n.
                                    if "g" in merged_variant.type:
                                        identifying_g_variant = merged_variant
                                        merged_variant = vm.g_to_n(merged_variant, tx_ac)
                                except AttributeError:
                                    pass

                            needs_a_push = True  # Keep the new vcf
                            break
                    else:
                        needs_a_push = False  # Restore old vcf
                        break

                else:
                    # Everything missed, assume no push required
                    needs_a_push = False
                    break

            # exon boundary hit. Break before intron
            elif working_pos == exon_end_genomic:
                break

            # Continue looping
            else:
                working_pos = working_pos + 1
                continue

        # Clear staging_loop
        staging_loop = 0

        # Create vcf dict
        if needs_a_push is True:
            # Re-sep pos-ref-alt (pos remains equal)
            ref = push_ref
            alt = push_alt

    # Dictionary the VCF
    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': normalized_hgvs_genomic,
                'merged_variant': merged_variant, 'identifying_variant': identifying_variant,
                'pre_merged_variant': pre_merged_variant, 'identifying_g_variant': identifying_g_variant}
    str_hgvs = vcfcp_to_hgvsstr(vcf_dict, hgvs_genomic)
    vcf_dict['str_hgvs'] = str_hgvs
    vcf_dict['needs_a_push'] = needs_a_push
    return vcf_dict

def hard_left_hgvs2vcf(hgvs_genomic, primary_assembly, hn, reverse_normalizer, sf, tx_ac, hdp, alt_aln_method,
                       hp, vm, mrg, genomic_ac=False, mapped_g=False, pre_norm=False ):
    """
    Designed specifically for gap handling.
    hard left pushes as 5 prime as possible and adds additional bases
    :param hgvs_genomic:
    :param primary_assembly:
    :param hn:
    :param reverse_normalizer:
    :param sf:
    :param tx_ac:
    :param hdp:
    :param alt_aln_method:
    :param hp:
    :param vm:
    :param genomic_ac:
    :param mapped_g: genomic mapping, used if a transcript type is input
    :return:
    """
    if hgvs_genomic.type == 'g':
        # Reverse normalize input prior to convert: NOTE will replace ref
        if pre_norm:
            reverse_normalized_hgvs_genomic = pre_norm
        else:
            reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic)
    else:
        # c. must be in n. format
        if hgvs_genomic.type == 'c':
            hgvs_genomic = vm.c_to_n(hgvs_genomic)
        if pre_norm:
            reverse_normalized_hgvs_genomic = pre_norm
        else:
            reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic)
        reverse_normalized_hgvs_genomic = pre_push_vcf_tx_g_map_fix(
                reverse_normalized_hgvs_genomic,
                hgvs_genomic,
                genomic_ac,
                vm,
                hn,# normaliser
                sf,# seq_fetcher
                mapped_g)

    # Chr
    if hgvs_genomic.type == 'g':
        chr = seq_data.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac, primary_assembly)
        if chr is None:
           chr = reverse_normalized_hgvs_genomic.ac
    else:
        chr = reverse_normalized_hgvs_genomic.ac

    # Identity
    if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity':
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'ins':
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
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'sub':
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'del':
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        # Recover sequences
        hgvs_del_seq_w_pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = hgvs_del_seq_w_pre_base
        alt = hgvs_del_seq_w_pre_base[0]

    # inv
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'inv':
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
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        # pos = str(start-1)
        # ref = bs + vcf_del_seq
        ref = vcf_del_seq
        alt = ins_seq
        if reverse_normalized_hgvs_genomic.posedit.edit.type == 'inv':
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'delins':
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
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start)
        ref = vcf_del_seq
        alt = vcf_del_seq[:1] + ins_seq

    # Duplications
    elif reverse_normalized_hgvs_genomic.posedit.edit.type == 'dup':
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
    # If possible, capture and alt variant that spans the gap
    merged_variant = False
    pre_merged_variant = False
    identifying_variant = False
    identifying_g_variant = False
    needs_a_push = False

    if chr != '' and pos != '' and ref != '' and alt != '':
        # Set exon boundary
        if genomic_ac is False:
            # Find the boundaries at the genomic level for the current exon
            exon_set = hdp.get_tx_exons(tx_ac, hgvs_genomic.ac, alt_aln_method)
            exon_start_genomic = None
            for exon in exon_set:
                if int(exon[7]) + 1 <= int(pos) <= int(exon[8]):
                    exon_start_genomic = int(exon[7] + 1)
                    break
        else:
            # Trick the system using transcript positions
            exon_set = hdp.get_tx_exons(hgvs_genomic.ac, genomic_ac, alt_aln_method)
            exon_start_genomic = None
            for exon in exon_set:
                if int(exon[5]) + 1 <= int(pos) <= int(exon[6]):
                    exon_start_genomic = int(exon[5] + 1)
                    break

        # Set loop variables for extending the push
        push_ref = ref
        push_alt = alt
        push_pos_by = 1
        needs_a_push = False
        staging_loop = 0
        if genomic_ac is False:
            genomic_ac = hgvs_genomic.ac
        # Loop and add bases - up to the range defined below - unless we go into an intron/past the transcript
        max_push_length = 50
        pos = int(pos)
        if pos - 1 - max_push_length < 0:
            max_push_length = pos -1

        flank_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), pos - max_push_length - 1, pos -1)
        for push in range(max_push_length):
            pre_pos = int(pos) - push_pos_by
            push_ref = flank_seq[-push_pos_by] + push_ref
            push_alt = flank_seq[-push_pos_by] + push_alt

            # Create a not_delins for normalisation checking
            var_end = pre_pos + len(push_ref) - 1
            normlize_check_variant = hgvs_delins_parts_to_hgvs_obj(
                    hgvs_genomic.ac,
                    hgvs_genomic.type,
                    pre_pos, push_ref, push_alt,
                    end=var_end,
                    offset_pos=True)
            # Check to see of we end up spanning a gap
            try:
                if hgvs_genomic.type != "g":
                    normlize_check_mapped = vm.n_to_g(normlize_check_variant, genomic_ac)
                else:
                    normlize_check_mapped = vm.g_to_n(normlize_check_variant, tx_ac, alt_aln_method)
            # Catch out-of-bounds errors
            except vvhgvs.exceptions.HGVSInvalidIntervalError:
                needs_a_push = False
                break

            """
            Break out from loop parameters
            """
            if normlize_check_mapped.posedit.pos.start.base > normlize_check_mapped.posedit.pos.end.base:
                needs_a_push = False
                break
            if len(normlize_check_mapped.posedit.edit.ref) <= 1:
                staging_loop = staging_loop + 1

            # Check here for the gap (Has it been crossed?) Note: if gap in tx, we have the whole gap spanned
            if (((len(normlize_check_mapped.posedit.edit.ref) != len(normlize_check_variant.posedit.edit.ref) and
                  len(normlize_check_mapped.posedit.edit.ref) > 1))
                    or
                    (normlize_check_variant.posedit.edit.type == 'identity')
                    and len(normlize_check_mapped.posedit.edit.alt) != len(normlize_check_variant.posedit.edit.ref)):

                # Add the identifying variant
                identifying_variant = normlize_check_variant
                if push == 0:  # Already crossing the gap so return original vcf
                    end_seq_check_variant = copy.deepcopy(normlize_check_variant)
                    # end_seq_check_variant.posedit.edit.alt = end_seq_check_variant.posedit.edit.ref

                else:
                    # Look to see if the gap has been identified by addition of bases in sequence
                    end_seq_check_variant = hgvs_delins_parts_to_hgvs_obj(
                            hgvs_genomic.ac,
                            hgvs_genomic.type,
                            normlize_check_variant.posedit.pos.start.base,
                            push_ref[0:2 + staging_loop], push_ref[0:2 + staging_loop],
                            end=normlize_check_variant.posedit.pos.start.base + 1 + staging_loop,
                            offset_pos=True)

                # Check to see of we end up spanning a gap at the last 2 bases
                if hgvs_genomic.type != "g":
                    end_seq_check_mapped = vm.n_to_g(end_seq_check_variant, genomic_ac)
                else:
                    end_seq_check_mapped = vm.g_to_n(end_seq_check_variant, tx_ac)

                # For genomic_variant mapped onto gapps, we end up with an offset
                start_offset = False
                end_offset = False
                try:
                    end_seq_check_mapped.posedit.pos.start.offset
                except AttributeError:
                    start_offset = False
                else:
                    if end_seq_check_mapped.posedit.pos.start.offset != 0:
                        start_offset = True
                try:
                    end_seq_check_mapped.posedit.pos.end.offset
                except AttributeError:
                    end_offset = False
                else:
                    if end_seq_check_mapped.posedit.pos.end.offset != 0:
                        end_offset = True
                if start_offset is True or end_offset is True:
                    # To identify the gap, we need to span it before mapping back
                    if end_offset is True:
                        end_seq_check_mapped.posedit.pos.end.base = end_seq_check_mapped.posedit.pos.start.base + 1
                        end_seq_check_mapped.posedit.pos.end.offset = 0
                        end_seq_check_mapped.posedit.edit.ref = ''
                        norml_end_seq_check_mapped = end_seq_check_mapped
                    elif start_offset is True:
                        end_seq_check_mapped.posedit.pos.start.base = end_seq_check_mapped.posedit.pos.end.base - 1
                        end_seq_check_mapped.posedit.pos.start.offset = 0
                        end_seq_check_mapped.posedit.edit.ref = ''
                        norml_end_seq_check_mapped = end_seq_check_mapped

                    # now map back onto original reference sequence
                    try:
                        norml_end_seq_check_mapped = vm.c_to_n(norml_end_seq_check_mapped)  # Need in n. context
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass
                    if hgvs_genomic.type == "g":
                        map_back = vm.n_to_g(norml_end_seq_check_mapped, genomic_ac)
                    else:
                        map_back = vm.g_to_n(norml_end_seq_check_mapped, tx_ac)

                    # Normalize variants, original and the gap induced variant (note, variant pre-normalized)
                    if map_back.posedit.pos.start.base > map_back.posedit.pos.end.base:
                        needs_a_push = False
                        break
                    map_back = hn.normalize(map_back)  # gap is left so normalize right
                    map_back_rn = reverse_normalizer.normalize(map_back)
                    try:
                        map_back = vm.c_to_n(map_back)  # Need in n. context
                        map_back_rn = vm.c_to_n(map_back_rn)
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass

                    # Can the variants be normalized together
                    if ((
                            (map_back.posedit.pos.end.base >=
                             reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                            and
                            (map_back.posedit.pos.end.base <=
                             reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                    )
                            or
                            (
                                    (map_back_rn.posedit.pos.end.base >=
                                     reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.end.base <=
                                     reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back.posedit.pos.start.base >=
                                     reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back.posedit.pos.start.base <=
                                     reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back_rn.posedit.pos.start.base >=
                                     reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.start.base <=
                                     reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )):

                        # Create a variant that reflects the impact of the gap.
                        # This uses variant merging
                        # We merge the "gap" variant and the variant itself
                        v1 = hgvs_genomic
                        v2 = map_back

                        if "g" not in hgvs_genomic.type:
                            v1 = vm.n_to_g(hgvs_genomic, genomic_ac)
                            v2 = vm.n_to_g(map_back, genomic_ac)
                        try:
                            v1 = reverse_normalizer.normalize(v1)
                            v2 = reverse_normalizer.normalize(v2)
                        except vvhgvs.exceptions.HGVSInvalidVariantError:
                            needs_a_push = True  # Restore old vcf
                            push_pos_by = push_pos_by + 1
                            break
                        else:
                            try:
                                if v1.posedit.pos.start.base < v2.posedit.pos.start.base:
                                    pre_merged_variant = mrg([v1, v2], reverse_normalizer, final_norm=False)
                                else:
                                    pre_merged_variant = mrg([v2, v1], reverse_normalizer, final_norm=False)
                                if "g" in pre_merged_variant.type:
                                    merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                else:
                                    merged_variant = pre_merged_variant
                            except utils.mergeHGVSerror as e:
                                needs_a_push = True  # Return new vcf only
                                push_pos_by = push_pos_by + 1
                                break
                            except vvhgvs.exceptions.HGVSParseError as e:
                                needs_a_push = True  # Return new vcf only
                                push_pos_by = push_pos_by + 1
                                break

                            # Ensure merged variant is not in a "non-intron" if mapped back to n.
                            if merged_variant is not False:
                                try:
                                    if (merged_variant.posedit.pos.start.offset != 0
                                            or merged_variant.posedit.pos.end.offset != 0):
                                        # Try from normalized genomic
                                        try:
                                            pre_merged_variant = hn.normalize(pre_merged_variant)
                                        except vvhgvs.exceptions.HGVSError:
                                            pass
                                        test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                        if (test_merged_variant.posedit.pos.start.offset == 0
                                                and test_merged_variant.posedit.pos.start.offset == 0):
                                            merged_variant = test_merged_variant
                                        else:
                                            pre_merged_variant = reverse_normalizer.normalize(pre_merged_variant)
                                            test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                            if (test_merged_variant.posedit.pos.start.offset == 0
                                                    and test_merged_variant.posedit.pos.start.offset == 0):
                                                merged_variant = pre_merged_variant
                                    # Map back to n.
                                    if "g" in merged_variant.type:
                                        merged_variant = vm.g_to_n(merged_variant, tx_ac)
                                except AttributeError:
                                    pass

                            needs_a_push = True  # Keep the new vcf
                            push_pos_by = push_pos_by + 1
                            break
                    else:
                        needs_a_push = False  # Restore old vcf
                        push_pos_by = push_pos_by + 1
                        break

                # Or we have identified the gap again at the expected position
                if len(end_seq_check_mapped.posedit.edit.ref) != len(end_seq_check_variant.posedit.edit.ref):
                    """
                    At this stage, we have done the following, illustrated by a  gap in transcript

                    g. NNNNNNNNNN
                    n. NNNNNNN--N

                    We forced the gap to be projected by making the end_seq_check_variant n.=

                             NN  Deletion in g.
                             |
                    g. NNNNNNNN
                    n. NNNNNNNN                   

                    So we need to make the g. == again before mapping back, which will make an ins in the n.
                    """

                    # Now normalize the variants to see if they meet
                    norml_end_seq_check_mapped = copy.deepcopy(end_seq_check_mapped)
                    norml_end_seq_check_mapped.posedit.edit.alt = norml_end_seq_check_mapped.posedit.edit.ref

                    # now map back onto original reference sequence
                    try:
                        norml_end_seq_check_mapped = vm.c_to_n(norml_end_seq_check_mapped)  # Need in n. context
                    except TypeError:
                        pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass

                    if hgvs_genomic.type == "g":
                        map_back = vm.n_to_g(norml_end_seq_check_mapped, genomic_ac)
                    else:
                        map_back = vm.g_to_n(norml_end_seq_check_mapped, tx_ac)

                    # In transcript gaps, this can push us fully into the gap
                    try:
                        if map_back.posedit.pos.start.offset != 0 and map_back.posedit.pos.start.offset != 0:
                            needs_a_push = False
                            push_pos_by = push_pos_by + 1
                            break
                    except AttributeError:
                        pass

                    # Normalize variants, original and the gap induced variant (note, variant pre-normalized)
                    if map_back.posedit.pos.start.base > map_back.posedit.pos.end.base:
                        needs_a_push = False
                        break
                    try:
                        map_back = hn.normalize(map_back)  # gap is left so normalize right
                        map_back_rn = reverse_normalizer.normalize(map_back)
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        needs_a_push = False  # Restore old vcf
                        push_pos_by = push_pos_by + 1
                        break
                    try:
                        map_back = vm.c_to_n(map_back)  # Need in n. context
                        map_back_rn = vm.c_to_n(map_back_rn)
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass

                    # Is the gap variant the same as the incoming variant?
                    if reverse_normalized_hgvs_genomic == map_back_rn:
                        needs_a_push = True
                        push_pos_by = 1
                        push_ref = end_seq_check_variant.posedit.edit.ref
                        push_alt = end_seq_check_variant.posedit.edit.alt
                        pos = end_seq_check_variant.posedit.pos.start.base
                        break

                    # Can the variants be normalized together
                    if ((
                            (map_back.posedit.pos.end.base >=
                             reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                            and
                            (map_back.posedit.pos.end.base <=
                             reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                    )
                            or
                            (
                                    (map_back_rn.posedit.pos.end.base >=
                                     reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.end.base <=
                                     reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back.posedit.pos.start.base >=
                                     reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back.posedit.pos.start.base <=
                                     reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back_rn.posedit.pos.start.base >=
                                     reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back_rn.posedit.pos.start.base <=
                                     reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )
                            or
                            (
                                    (map_back.posedit.pos.start.base <=
                                     reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
                                    and
                                    (map_back.posedit.pos.end.base >=
                                     reverse_normalized_hgvs_genomic.posedit.pos.end.base + 1)
                            )):

                        # Create a variant that reflects the impact of the gap.
                        # This uses variant merging
                        # We merge the "gap" variant and the variant itself
                        v1 = hgvs_genomic
                        v2 = map_back

                        if "g" not in hgvs_genomic.type:
                            if (hgvs_genomic.posedit.edit.type == "dup"
                                    and map_back.posedit.edit.type == "del"
                                    and hgvs_genomic.posedit.edit.ref == map_back.posedit.edit.ref):
                                v1 = vm.n_to_g(hgvs_genomic, genomic_ac)
                                v2 = vm.n_to_g(map_back, genomic_ac)

                            elif (hgvs_genomic.posedit.edit.type == "del"
                                    and map_back.posedit.edit.type == "dup"
                                    and hgvs_genomic.posedit.edit.ref == map_back.posedit.edit.ref):
                                v1 = vm.n_to_g(hgvs_genomic, genomic_ac)
                                v2 = vm.n_to_g(map_back, genomic_ac)

                            elif ((map_back.posedit.edit.type == "dup" or map_back.posedit.edit.type == "del") and
                                  hgvs_genomic.posedit.pos.start.base > map_back.posedit.pos.end.base + 1):
                                v1 = vm.n_to_g(hgvs_genomic, genomic_ac)
                                v3 = hgvs_delins_parts_to_hgvs_obj(
                                        v2.ac,
                                        v2.type,
                                        v2.posedit.pos,
                                        v2.posedit.edit.ref,
                                        v2.posedit.edit.ref)
                                v2 = vm.n_to_g(v3, genomic_ac)

                            else:
                                v1 = vm.n_to_g(hgvs_genomic, genomic_ac)
                                v2 = vm.n_to_g(map_back, genomic_ac)


                        # Known examples of incorrect formatting from vm
                        ################################################

                        # 1. vm causes an insertion length of > 1 because of the gap - issue #392
                        if "ins" in v1.posedit.edit.type and "sub" in v2.posedit.edit.type:
                            try:
                                v1 = hn.normalize(v1)
                            except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                                if "insertion length must be 1" in str(e):
                                    v1 = hgvs_to_delins_hgvs(v1, hp, hn, allow_fix=True)
                                    identifying_g_variant = v1

                        elif "ins" in v2.posedit.edit.type and "sub" in v1.posedit.edit.type:
                            try:
                                v2 = hn.normalize(v2)
                            except vvhgvs.exceptions.HGVSInvalidVariantError as e:
                                if "insertion length must be 1" in str(e):
                                    v2 = hgvs_to_delins_hgvs(v2, hp, hn, allow_fix=True)
                                    identifying_g_variant = v2

                        try:
                            v1 = reverse_normalizer.normalize(v1)
                            v2 = reverse_normalizer.normalize(v2)
                        except vvhgvs.exceptions.HGVSInvalidVariantError:
                            needs_a_push = True  # Restore old vcf
                            push_pos_by = push_pos_by + 1
                            break
                        else:
                            try:
                                if v1.posedit.pos.start.base < v2.posedit.pos.start.base:
                                    pre_merged_variant = mrg([v1, v2], reverse_normalizer, final_norm=False)
                                else:
                                    pre_merged_variant = mrg([v2, v1], reverse_normalizer, final_norm=False)
                                if "g" in pre_merged_variant.type:
                                    # identifying_g_variant = pre_merged_variant
                                    merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                else:
                                    merged_variant = pre_merged_variant
                            except utils.mergeHGVSerror as e:
                                needs_a_push = True  # Return new vcf only
                                push_pos_by = push_pos_by + 1
                                break
                            except vvhgvs.exceptions.HGVSParseError as e:
                                needs_a_push = True  # Return new vcf only
                                break

                            # Ensure merged variant is not in a "non-intron" if mapped back to n.
                            if merged_variant is not False:
                                try:
                                    if (merged_variant.posedit.pos.start.offset != 0
                                            or merged_variant.posedit.pos.start.offset != 0):
                                        # Try from normalized genomic
                                        pre_merged_variant = hn.normalize(pre_merged_variant)
                                        test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                        if (test_merged_variant.posedit.pos.start.offset == 0
                                                and test_merged_variant.posedit.pos.start.offset == 0):
                                            merged_variant = pre_merged_variant
                                        else:
                                            pre_merged_variant = reverse_normalizer.normalize(pre_merged_variant)
                                            test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                            if (test_merged_variant.posedit.pos.start.offset == 0
                                                    and test_merged_variant.posedit.pos.start.offset == 0):
                                                merged_variant = pre_merged_variant

                                    # Map back to n.
                                    if "g" in merged_variant.type:
                                        merged_variant = vm.g_to_n(merged_variant, tx_ac)
                                except AttributeError:
                                    pass

                            needs_a_push = True  # Keep the new vcf
                            push_pos_by = push_pos_by + 1
                            break
                    else:
                        needs_a_push = False  # Restore old vcf
                        push_pos_by = push_pos_by + 1
                        break

                else:
                    # Everything missed, assume no push required
                    needs_a_push = False
                    push_pos_by = push_pos_by + 1
                    break

            # exon boundary hit. Break before intron
            elif pre_pos == exon_start_genomic:
                push_pos_by = push_pos_by + 1
                break

            # Continue looping
            else:
                push_pos_by = push_pos_by + 1
                continue

        # Clear staging_loop
        staging_loop = 0

        # Populate vcf dict
        if needs_a_push is True:
            # Re-sep pos-ref-alt
            pos = (int(pos) - (push_pos_by - 1))
            ref = push_ref
            alt = push_alt

    # Dictionary the VCF
    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': reverse_normalized_hgvs_genomic,
                'merged_variant': merged_variant, 'identifying_variant': identifying_variant,
                'pre_merged_variant': pre_merged_variant, 'identifying_g_variant': identifying_g_variant}
    str_hgvs = vcfcp_to_hgvsstr(vcf_dict, hgvs_genomic)
    vcf_dict['str_hgvs'] = str_hgvs
    vcf_dict['needs_a_push'] = needs_a_push
    return vcf_dict


def hgvs_ref_alt(hgvs_variant, sf):
    if re.search(r'[GATC]+=', str(hgvs_variant.posedit)):
        ref = hgvs_variant.posedit.edit.ref
        alt = hgvs_variant.posedit.edit.ref

    # Insertions
    elif 'ins' in str(hgvs_variant.posedit) and 'del' not in str(hgvs_variant.posedit):
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
    elif '>' in str(hgvs_variant.posedit):
        ref = hgvs_variant.posedit.edit.ref
        alt = hgvs_variant.posedit.edit.alt

    # Deletions
    elif 'del' in str(hgvs_variant.posedit) and 'ins' not in str(hgvs_variant.posedit):
        ref = hgvs_variant.posedit.edit.ref
        alt = ''

    # inv
    elif 'inv' in str(hgvs_variant.posedit):
        ref = hgvs_variant.posedit.edit.ref
        my_seq = Seq(ref)
        alt = str(my_seq.reverse_complement())

    # Delins
    elif 'del' in str(hgvs_variant.posedit) and 'ins' in str(hgvs_variant.posedit):
        ref = hgvs_variant.posedit.edit.ref
        alt = hgvs_variant.posedit.edit.alt

    # Duplications
    elif 'dup' in str(hgvs_variant.posedit):
        ref = hgvs_variant.posedit.edit.ref
        alt = hgvs_variant.posedit.edit.ref + hgvs_variant.posedit.edit.ref
    else:
        ref = ''
        alt = ''

    ref_alt_dict = {'ref': ref, 'alt': alt}
    return ref_alt_dict


def incomplete_alignment_mapping_t_to_g(validator, variant):
    output = None
    mapping_options = validator.hdp.get_tx_mapping_options(variant.input_parses.ac)
    for option in mapping_options:
        if option[2] == validator.alt_aln_method and "NC_" not in option[1]:
            in_assembly = seq_data.to_chr_num_refseq(option[1], variant.primary_assembly)
            if in_assembly is not None:
                try:
                    output = validator.vm.t_to_g(variant.input_parses, option[1])
                    if variant.input_parses.posedit.edit.type == "identity":
                        output.posedit.edit.alt = output.posedit.edit.ref
                except vvhgvs.exceptions.HGVSError:
                    pass
    return output

# <LICENSE>
# Copyright (C) 2016-2025 VariantValidator Contributors
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
