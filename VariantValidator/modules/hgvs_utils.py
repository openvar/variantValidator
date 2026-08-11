import copy, re
from . import seq_data
from . import utils
from . import hgvs_position_utils

# Import vv_hgvs modules
import vvhgvs
import vvhgvs.exceptions
from vvhgvs.location import BaseOffsetInterval, BaseOffsetPosition, Interval, SimplePosition

# used to set coordinate origin point i.e. seq start vs CDS start/end
from vvhgvs.enums import Datum
from vvhgvs.location import AAPosition
from vvhgvs.edit import AASub, AARefAlt, Dup, NARefAlt
from vvhgvs.posedit import PosEdit
import logging

logger = logging.getLogger(__name__)

# Custom error handling
class PseudoVCF2HGVSError(Exception):
    pass  # Pass exception

class VVPosEdit(PosEdit):
    "override class for posedit to get VV specific formatting"
    met_variation = None
    expanded_rep = None
    def __init__(
            self,pos,edit,
            uncertain=False,
            nucleotide_not_equal = None,
            met_variation = None,
            expanded_rep = None):
        PosEdit.__init__(self,pos=pos,edit=edit,uncertain=uncertain)
        # used to append methylation variation to the output
        self.met_variation = met_variation
        # used for formatting prot consequence of nuc input (for Ter/*)
        self.nucleotide_not_equal = nucleotide_not_equal
        # used for expanded repeat type variant output (do not normalise!)
        self.expanded_rep = expanded_rep

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
            else: # must be NARefAlt, not dup
                if self.expanded_rep:
                    # dup can be converted into this but should be stated as NARefAlt, this is
                    # not robust to mapping and/or normalisation (may norm down to exclude part
                    # of expanded repeat's range)
                    rep = int(len(self.edit.alt)/len(self.expanded_rep))
                    formatted_str = f"{self.pos.format(conf)}{self.expanded_rep}[{rep}]"
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

def unset_hgvs_obj_ref(hgvs, vf_mode=False):
    """
    Remove/unset ref bases from an HGVS object where appropriate for output,
    without re-parsing from text.
    """
    edit = hgvs.posedit.edit

    # Nothing to alter for unknown/empty protein edits.
    if hgvs.type == "p" and vf_mode:
        if edit == "?":
            return hgvs

        if (
            getattr(edit, "ref", None) is None
            and getattr(edit, "alt", None) is None
        ):
            return hgvs

    logger.info(
        "unset_hgvs_obj_ref ENTER: hgvs=%s, edit_type=%s, ref=%r, alt=%r",
        hgvs,
        edit.type,
        getattr(edit, "ref", None),
        getattr(edit, "alt", None),
    )

    # Identity variants must remain identities. In VariantFormatter mode,
    # multi-base identities should be rendered without explicit sequence.
    if edit.type == "identity":
        if vf_mode and len(edit.ref) > 1:
            logger.info(
                "unset_hgvs_obj_ref IDENTITY: hgvs=%s, ref=%r, alt=%r",
                hgvs,
                edit.ref,
                edit.alt,
            )

            edit.ref = ""
            edit.alt = ""
            hgvs.posedit.edit = edit

            logger.info(
                "unset_hgvs_obj_ref EXIT: hgvs=%s, edit_type=%s, ref=%r, alt=%r",
                hgvs,
                edit.type,
                edit.ref,
                edit.alt,
            )

        return hgvs

    if edit.type in ["inv", "dup"]:
        edit.ref = ""

    elif edit.ref is not None and edit.alt is not None:
        if len(edit.alt) == 1 and len(edit.ref) == 1:
            return hgvs

        if "N" in edit.ref:
            return hgvs

        edit.ref = ""

    elif edit.ref is not None:
        edit.ref = ""

    hgvs.posedit.edit = edit

    logger.info(
        "unset_hgvs_obj_ref EXIT: hgvs=%s, edit_type=%s, ref=%r, alt=%r",
        hgvs,
        edit.type,
        getattr(edit, "ref", None),
        getattr(edit, "alt", None),
    )

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
        start_pos = vvhgvs.location.BaseOffsetPosition(
                base=int(prefix),
                offset=-int(offset),
                datum=coordinate_origin)
        if not end_pos:
            end_offset = int(offset) - (length -1)
            if end_offset < 0:
                prefix = int(prefix) - end_offset
                end_pos = vvhgvs.location.BaseOffsetPosition(
                        base=int(prefix),
                        offset=0,
                        datum=end_coordinate_origin)
            else :
                end_offset = int(offset) - (length -1)
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
                base=int(loc_start),
                datum=coordinate_origin)
        if not end_pos:
            end_pos = vvhgvs.location.BaseOffsetPosition(
                    base=int(loc_start)+ length -1,
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

def hgvs_obj_from_existing_edit(ref_ac, ref_type, starts, edit,
                                end=None, offset_pos=False,
                                unc_posedit=None):
    """
    Build an HGVS variant object using an existing HGVS edit.

    starts may be:
      - an existing Interval/BaseOffsetInterval
      - an existing SimplePosition/BaseOffsetPosition
      - a string/integer position

    end may likewise be an existing position or a string/integer position.
    """
    # Existing interval: use directly.
    if isinstance(starts, (BaseOffsetInterval, Interval)):
        position = starts

    # Existing HGVS position objects: construct the appropriate interval
    # directly and preserve the object lifecycle.
    elif isinstance(starts, BaseOffsetPosition):
        end_pos = end if end is not None else starts

        if not isinstance(end_pos, BaseOffsetPosition):
            _, end_pos = _hgvs_offset_pos_from_str_in(
                starts,
                len(edit.ref) if edit.ref is not None else 2,
                ref_type=ref_type,
                end=end
            )

        position = BaseOffsetInterval(
            start=starts,
            end=end_pos
        )

    elif isinstance(starts, SimplePosition):
        if end is None:
            if edit.ref is None:
                end_pos = SimplePosition(base=starts.base + 1)
            else:
                end_pos = SimplePosition(
                    base=starts.base + len(edit.ref) - 1
                )
        elif isinstance(end, SimplePosition):
            end_pos = end
        else:
            end_pos = SimplePosition(base=int(end))

        position = Interval(
            start=starts,
            end=end_pos
        )

    # String/numeric coding positions.
    elif offset_pos or ref_type in ('c', 'n'):
        length = 2 if edit.ref is None else len(edit.ref)

        start_pos, end_pos = _hgvs_offset_pos_from_str_in(
            starts,
            length,
            ref_type=ref_type,
            end=end
        )

        position = BaseOffsetInterval(
            start=start_pos,
            end=end_pos
        )

    # String/numeric simple positions.
    else:
        start_pos = int(starts)

        if end is not None:
            end_pos = int(end)
        else:
            end_pos = start_pos + len(edit.ref) - 1

        position = Interval(
            start=SimplePosition(base=start_pos),
            end=SimplePosition(base=end_pos)
        )

    return vvhgvs.sequencevariant.SequenceVariant(
        ac=ref_ac,
        type=ref_type,
        posedit=VVPosEdit(
            position,
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
            posedit=VVPosEdit(
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

def _select_pvcf_normalizer(normalization_direction, reverse_normalizer, validator):
    """Return the normalizer for the requested VCF normalisation direction."""
    return {3: validator.hn, 5: reverse_normalizer}[normalization_direction]

def _pvcf_to_hgvs_input(query):
    """Convert pseudo-VCF input to an HGVS-like substitution description."""
    query = query.replace(":", "-")
    vcf_elements = query.split("-")
    if re.search(r"-\d+-[GATC]+-[GATC]+", query):
        return f"{vcf_elements[0]}:{vcf_elements[1]}{vcf_elements[2]}>{vcf_elements[3]}"
    if re.search(r"-\d+-[GATC]+-", query):
        return f"{vcf_elements[0]}:{vcf_elements[1]}{vcf_elements[2]}>{vcf_elements[2]}"
    raise PseudoVCF2HGVSError("Unsupported format: VCF specification 4.1 or later")

def _resolve_pvcf_accession(accession, selected_assembly, validator):
    """Resolve a pseudo-VCF chromosome/LRG identifier to an accession."""
    if accession.startswith(("NC_", "NG_", "NW_", "NT_")):
        return accession
    if re.fullmatch(r"LRG_\d+", accession):
        return validator.db.get_refseq_id_from_lrg_id(accession)
    chr_num = accession.strip().upper()
    if chr_num.startswith("CHR"):
        chr_num = chr_num[3:]
    accession = seq_data.get_accession(chr_num, selected_assembly)
    if accession is None: # Accession is not set
        raise PseudoVCF2HGVSError(
            f"{chr_num} is not part of genome build {selected_assembly} or is not supported"
        )
    return accession

def _pvcf_get_alleles(position_and_edit):
    """Extract reference and alternate alleles from a pseudo-VCF edit."""
    match = re.search(r"([GATCgatc]+)>([GATCgatc]+)", position_and_edit)
    if match is None:
        raise PseudoVCF2HGVSError("Unsupported format: VCF specification 4.1 or later!")
    return match.groups()

def _pvcf_build_simple_hgvs(accession, ref_type, position_and_edit):
    """Build an HGVS object for a single-base pseudo-VCF substitution."""
    match = re.fullmatch(r"(\d+)(?:_(\d+))?([GATCgatc])>([GATCgatc])", position_and_edit)
    if match is None:
        raise PseudoVCF2HGVSError(
            f"Unable to parse pseudo-VCF substitution: {position_and_edit}"
        )
    start, end, ref, alt = match.groups()
    return hgvs_delins_parts_to_hgvs_obj(
        accession, ref_type, int(start), ref, alt,
        end=int(end) if end is not None else None,
    )

def _pvcf_build_multibase_hgvs(accession, ref_type, position_and_edit):
    """Build an HGVS object for a multi-base pseudo-VCF edit."""
    ref, alt = _pvcf_get_alleles(position_and_edit)
    not_sub = f"{accession}{ref_type}{position_and_edit}"
    if re.search(r"[0-9]+_[0-9]+", not_sub):
        beginning_string, middle_string = not_sub.split(":", 1)
        middle_string = middle_string.split("_", 1)[0]
        not_sub = f"{beginning_string}:{middle_string}{ref}>{alt}"

    ref_ac, _, remainder = not_sub.partition(":")
    hgvs_ref_type, _, posedit = remainder.partition(".")
    pos_ref, _, insert = posedit.partition(">")
    match = re.search(r"([0-9]+)([GATCgatc]+)", pos_ref)
    if match is None:
        raise PseudoVCF2HGVSError(f"Unable to parse reference sequence from {not_sub}")

    delete = match.group(2)
    starts = posedit.split(delete, 1)[0]
    temporary = hgvs_delins_parts_to_hgvs_obj(
        ref_ac, hgvs_ref_type, starts, delete[0], insert
    )
    temporary.posedit.edit.ref = delete
    start = temporary.posedit.pos.start

    if isinstance(start, BaseOffsetPosition):
        if start.offset < 0:
            end = BaseOffsetPosition(
                base=start.base, offset=-start.offset + len(delete), datum=start.datum
            )
        else: # Make base offset position
            end = BaseOffsetPosition(
                base=start.base, offset=start.offset + len(delete) - 1, datum=start.datum
            )
    else: # Make simple position
        end = SimplePosition(base=start.base + len(delete) - 1)

    return hgvs_obj_from_existing_edit(
        ref_ac, hgvs_ref_type, start,
        vvhgvs.edit.NARefAlt(ref=delete, alt=insert), end=end,
    )

def _pvcf_build_hgvs_object(accession, ref_type, position_and_edit):
    """Build a pseudo-VCF HGVS object while preserving HGVS position objects."""
    ref, alt = _pvcf_get_alleles(position_and_edit)
    if len(ref) == 1 and len(alt) == 1 and "," not in position_and_edit:
        return _pvcf_build_simple_hgvs(accession, ref_type, position_and_edit)
    return _pvcf_build_multibase_hgvs(accession, ref_type, position_and_edit)

def pvcf_to_hgvs(query, selected_assembly, normalization_direction, reverse_normalizer, validator):
    """Convert a pseudo-VCF description to an HGVS object."""
    selected_normalizer = _select_pvcf_normalizer(
        normalization_direction, reverse_normalizer, validator
    )
    query = _pvcf_to_hgvs_input(query)
    accession, position_and_edit = query.split(":", 1)
    accession = _resolve_pvcf_accession(accession, selected_assembly, validator)
    hgvs_object = _pvcf_build_hgvs_object(accession, ":g.", position_and_edit)
    return selected_normalizer.normalize(hgvs_object)

def _hgvs_vcf_sequence(hgvs, sf, report_mode=False):
    """Convert an HGVS edit to VCF position/ref/alt components."""
    edit = hgvs.posedit.edit
    position = hgvs.posedit.pos
    edit_type = edit.type

    if edit_type == "identity":
        return str(position.start), edit.ref, edit.ref
    if edit_type == "ins":
        end = int(position.end.base)
        start = int(position.start.base)
        ref_seq = sf.fetch_seq(hgvs.ac, start - 1, end - 1)
        return start, ref_seq, ref_seq + edit.alt
    if edit_type == "sub":
        return str(position), edit.ref, edit.alt
    if edit_type == "del":
        end = int(position.end.base)
        start = int(position.start.base)
        adj_start = start - 2
        if report_mode and adj_start < 0:
            ref_seq = sf.fetch_seq(hgvs.ac, start, end + 1)
            return "1", ref_seq, ref_seq[-1]
        ref_seq = sf.fetch_seq(hgvs.ac, adj_start, end)
        return str(start - 1), ref_seq, ref_seq[0]
    if edit_type == "inv":
        start = int(position.start.base)
        end = int(position.end.base)
        ref_seq = getattr(edit, "ref", None)
        if not ref_seq:
            ref_seq = sf.fetch_seq(hgvs.ac, start - 1, end)
        return str(start), ref_seq, utils.simple_dna_revcomp(ref_seq)
    if edit_type == "delins":
        start = int(position.start.base)
        end = int(position.end.base)
        ins_seq = edit.alt or ""
        if report_mode:
            ref_seq = sf.fetch_seq(hgvs.ac, start - 1, end)
            return str(start), ref_seq, ins_seq
        ref_seq = sf.fetch_seq(hgvs.ac, start - 2, end)
        return str(start - 1), ref_seq, ref_seq[0] + ins_seq
    if edit_type == "dup":
        end = int(position.end.base)
        start = int(position.start.base)
        ref_seq = sf.fetch_seq(hgvs.ac, start - 2, end)
        if report_mode:
            return str(start - 1), ref_seq[0], ref_seq
        return str(start - 1), ref_seq, ref_seq + edit.ref
    return "", "", ""

def _hgvs2vcf_chromosome(hgvs, primary_assembly):
    return seq_data.get_chr_num_ucsc(hgvs.ac, primary_assembly) or hgvs.ac

def _report_vcf_chromosomes(hgvs, primary_assembly):
    if primary_assembly == "All":
        gen_name_map = {"GRCh37": "grch37", "hg19": "hg19", "GRCh38": "grch38", "hg38": "hg38"}
        chrs = {}
        for genome, output_name in gen_name_map.items():
            if not seq_data.is_supported_for_mapping(hgvs.ac, genome):
                continue
            chrom = (
                seq_data.get_chr_num_refseq(hgvs.ac, genome)
                if genome.startswith("GRC")
                else seq_data.get_chr_num_ucsc(hgvs.ac, genome)
            )
            chrs[output_name] = chrom or hgvs.ac
        return "", "", chrs

    ucsc_pa = ""
    grc_pa = ""
    if "GRC" in primary_assembly:
        if "37" in primary_assembly:
            ucsc_pa = "hg19"
            grc_pa = primary_assembly # inherits
        if "38" in primary_assembly:
            ucsc_pa = "hg38"
            grc_pa = primary_assembly # inherits
    else: # When hg formart us used rather than GRCh
        if "19" in primary_assembly:
            ucsc_pa = primary_assembly # inherits
            grc_pa = "GRCh37"
        if "38" in primary_assembly:
            ucsc_pa = primary_assembly # inherits
            grc_pa = "GRCh38"

    return (
        seq_data.get_chr_num_ucsc(hgvs.ac, ucsc_pa) or hgvs.ac,
        seq_data.get_chr_num_refseq(hgvs.ac, grc_pa) or hgvs.ac,
        {},
    )

def _add_vcf_flanks(hgvs, sf, pos, ref, alt, extra_flank_bases):
    if extra_flank_bases <= 0:
        return pos, ref, alt
    original_pos = pos
    pos = str(int(pos) - extra_flank_bases)
    left_flank = sf.fetch_seq(hgvs.ac, int(pos) - 1, int(original_pos) - 1)
    right_flank = sf.fetch_seq(
        hgvs.ac,
        int(original_pos) + len(ref) - 1,
        int(original_pos) + len(ref) - 1 + extra_flank_bases,
    )
    return pos, left_flank + ref + right_flank, left_flank + alt + right_flank

def hgvs2vcf(hgvs_genomic, primary_assembly, reverse_normalizer, sf, extra_flank_bases=0):
    """Convert HGVS to the standard VCF representation."""
    normalized = (
        hgvs_genomic
        if reverse_normalizer is None
        else reverse_normalizer.normalize(hgvs_genomic)
    )
    chrom = _hgvs2vcf_chromosome(normalized, primary_assembly)
    pos, ref, alt = _hgvs_vcf_sequence(normalized, sf)

    if chrom and pos and ref and alt and len(ref) > 1:
        if normalized.posedit.edit.type == "identity":
            pos_int = int(pos) - 1
            previous = sf.fetch_seq(normalized.ac, pos_int - 1, pos_int)
            pos = str(pos_int)
            ref = previous + ref
            alt = previous + alt
    pos, ref, alt = _add_vcf_flanks(normalized, sf, pos, ref, alt, extra_flank_bases)
    return {"chr": chrom, "pos": pos, "ref": ref, "alt": alt, "normalized_hgvs": normalized}

def report_hgvs2vcf(hgvs_genomic, primary_assembly, reverse_normalizer, sf):
    """Return the report VCF representation without additional flank bases."""
    normalized = reverse_normalizer.normalize(hgvs_genomic)
    ucsc_chr, grc_chr, chrs = _report_vcf_chromosomes(normalized, primary_assembly)
    pos, ref, alt = _hgvs_vcf_sequence(normalized, sf, report_mode=True)
    return {
        "pos": str(pos),
        "ref": ref,
        "alt": alt,
        "ucsc_chr": ucsc_chr,
        "grc_chr": grc_chr,
        "normalized_hgvs": normalized,
        "chrs_by_genome": chrs,
    }

def pos_lock_hgvs2vcf(hgvs_genomic,
                      primary_assembly,
                      reverse_normalizer,
                      sf):
    """Return an in-situ VCF representation without normalisation."""
    if hgvs_genomic.posedit.edit.ref == "":
        hgvs_genomic.posedit.edit.ref = sf.fetch_seq(
            hgvs_genomic.ac,
            hgvs_genomic.posedit.pos.start.base - 1,
            hgvs_genomic.posedit.pos.end.base,
        )

    normalized = hgvs_genomic
    if normalized.posedit.edit.type == "identity" and not normalized.posedit.edit.ref:
        normalized = reverse_normalizer.normalize(normalized)

    chrom = _hgvs2vcf_chromosome(normalized, primary_assembly)
    pos, ref, alt = _hgvs_vcf_sequence(normalized, sf)
    return {"chr": chrom, "pos": pos, "ref": ref, "alt": alt, "normalized_hgvs": normalized}

def pre_push_vcf_tx_g_map_fix(
        norm_hgvs_transcript,
        un_norm_hgvs,
        genomic_ac,
        var_mapper,
        normaliser,
        seq_fetcher,
        mapped_hgvs
):
    """
    Prepare a transcript HGVS variant for left/right VCF pushing when its
    transcript-to-genome mapping crosses a genomic alignment gap.

    Parameters
    ----------
    norm_hgvs_transcript
        Normalised or reverse-normalised n. HGVS variant.
    un_norm_hgvs
        Unnormalised transcript HGVS variant. This is retained so that the
        original edit type can be inspected where normalisation has changed
        the representation.
    genomic_ac
        Genomic accession to which the transcript variant is mapped.
    var_mapper
        Variant mapper used for n-to-g and g-to-n mapping.
    normaliser
        HGVS normaliser.
    seq_fetcher
        Sequence fetcher used when reconstructing duplication alleles.
    mapped_hgvs
        Existing genomic mapping. If supplied, the n-to-g mapping step is
        skipped.

    Returns
    -------
    SequenceVariant
        The input transcript HGVS variant, or an expanded transcript variant
        that spans an alignment gap where required.
    """

    # Reuse an existing genomic mapping where available.
    if mapped_hgvs:
        hgvs_genomic_g = mapped_hgvs
    else:
        hgvs_genomic_g = var_mapper.n_to_g(
            norm_hgvs_transcript,
            genomic_ac
        )

    # Most mappings require no correction. The special handling below is
    # required when an alignment gap produces reversed genomic coordinates.
    try:
        normaliser.normalize(hgvs_genomic_g)
        return norm_hgvs_transcript
    except vvhgvs.exceptions.HGVSInvalidVariantError as e:
        if "base start position must be <= end position" not in str(e):
            return norm_hgvs_transcript

    hgvs_genomic_g_identity = copy.deepcopy(hgvs_genomic_g)

    # Reverse the genomic interval so that it can be represented and
    # normalised as an identity spanning the alignment gap.
    start_base = hgvs_genomic_g_identity.posedit.pos.end.base
    end_base = hgvs_genomic_g_identity.posedit.pos.start.base

    hgvs_genomic_g_identity.posedit.pos.start.base = start_base
    hgvs_genomic_g_identity.posedit.pos.end.base = end_base

    # Duplications have no alt sequence in the HGVS object, so populate the
    # reference sequence before converting the edit to delins.
    if hgvs_genomic_g_identity.posedit.edit.type == "dup":
        hgvs_genomic_g_identity.posedit.edit.ref = seq_fetcher.fetch_seq(
            hgvs_genomic_g_identity.ac,
            start_base - 1,
            end_base
        )
        hgvs_genomic_g_identity = hgvs_to_delins_hgvs(
            hgvs_genomic_g_identity,
            None,
            normaliser
        )

    hgvs_genomic_g_identity.posedit.edit.ref = ""
    hgvs_genomic_g_identity.posedit.edit.alt = ""
    hgvs_genomic_g_identity = normaliser.normalize(
        hgvs_genomic_g_identity
    )

    # Map the genomic gap back onto the transcript to determine the
    # transcript interval that must be spanned.
    hgvs_genomic_n_gap = var_mapper.g_to_n(
        hgvs_genomic_g_identity,
        norm_hgvs_transcript.ac
    )

    hgvs_genomic_n_identity = copy.deepcopy(hgvs_genomic_n_gap)
    hgvs_genomic_n_identity.posedit.edit.ref = ""
    hgvs_genomic_n_identity.posedit.edit.alt = ""
    hgvs_genomic_n_identity = normaliser.normalize(
        hgvs_genomic_n_identity
    )

    # Convert edits that need explicit ref/alt alleles before assembling the
    # transcript variant across the gap. Check the unnormalised variant as
    # well because right normalisation can alter a deletion representation.
    if (
        un_norm_hgvs.posedit.edit.type == "del"
        or norm_hgvs_transcript.posedit.edit.type in {"del", "ins", "dup"}
    ):
        norm_hgvs_transcript = hgvs_to_delins_hgvs(
            norm_hgvs_transcript,
            None,
            normaliser
        )

    gap_start = hgvs_genomic_n_identity.posedit.pos.start.base
    gap_end = hgvs_genomic_n_identity.posedit.pos.end.base
    variant_start = norm_hgvs_transcript.posedit.pos.start.base
    variant_end = norm_hgvs_transcript.posedit.pos.end.base

    # Build sequence preceding the submitted variant when the transcript gap
    # begins before the variant.
    if gap_start < variant_start:
        v1 = copy.deepcopy(hgvs_genomic_n_identity)
        v1.posedit.pos.end.base = variant_start - 1
        v1.posedit.edit.ref = ""
        v1.posedit.edit.alt = ""
        v1 = normaliser.normalize(v1)

        # Build sequence following the variant when the gap extends beyond it.
        if gap_end > variant_end:
            v3 = copy.deepcopy(hgvs_genomic_n_identity)
            v3.posedit.pos.start.base = variant_end + 1
            v3.posedit.edit.ref = ""
            v3.posedit.edit.alt = ""
            v3 = normaliser.normalize(v3)
        else:
            v3 = None

        hgvs_genomic_n_assembled = copy.deepcopy(
            hgvs_genomic_n_identity
        )

        if v3 is not None:
            hgvs_genomic_n_assembled.posedit.pos.end.base = (
                v3.posedit.pos.end.base
            )
            ass_ref = (
                v1.posedit.edit.ref
                + norm_hgvs_transcript.posedit.edit.ref
                + v3.posedit.edit.ref
            )
            ass_alt = (
                v1.posedit.edit.alt
                + norm_hgvs_transcript.posedit.edit.alt
                + v3.posedit.edit.alt
            )
        else:
            hgvs_genomic_n_assembled.posedit.pos.end.base = variant_end
            ass_ref = (
                v1.posedit.edit.ref
                + norm_hgvs_transcript.posedit.edit.ref
            )
            ass_alt = (
                v1.posedit.edit.alt
                + norm_hgvs_transcript.posedit.edit.alt
            )

    else:
        # No sequence is required before the variant. Only append sequence
        # following it if the transcript gap extends beyond the variant.
        if gap_end > variant_end:
            v3 = copy.deepcopy(hgvs_genomic_n_identity)
            v3.posedit.pos.start.base = variant_end + 1
            v3.posedit.edit.ref = ""
            v3.posedit.edit.alt = ""
            v3 = normaliser.normalize(v3)
        else:
            v3 = None

        hgvs_genomic_n_assembled = copy.deepcopy(
            hgvs_genomic_n_identity
        )

        if v3 is not None:
            hgvs_genomic_n_assembled.posedit.pos.end.base = (
                v3.posedit.pos.end.base
            )
            ass_ref = (
                norm_hgvs_transcript.posedit.edit.ref
                + v3.posedit.edit.ref
            )
            ass_alt = (
                norm_hgvs_transcript.posedit.edit.alt
                + v3.posedit.edit.alt
            )
        else:
            hgvs_genomic_n_assembled.posedit.pos.end.base = variant_end
            ass_ref = norm_hgvs_transcript.posedit.edit.ref
            ass_alt = norm_hgvs_transcript.posedit.edit.alt

    hgvs_genomic_n_assembled.posedit.edit.ref = ass_ref
    hgvs_genomic_n_assembled.posedit.edit.alt = ass_alt

    return hgvs_genomic_n_assembled


def _prepare_hard_hgvs(
        hgvs_genomic,
        primary_assembly,
        normalizer,
        hn,
        sf,
        vm,
        tx_ac,
        map_dat,
        alt_aln_method,
        genomic_ac,
        mapped_g,
        pre_norm,
):
    """Prepare the HGVS object and initial VCF components for hard pushing.

    The left and right push algorithms deliberately remain separate. This
    helper only handles their identical input preparation and VCF conversion.
    """
    if hgvs_genomic.type == "c":
        hgvs_genomic = vm.c_to_n(hgvs_genomic)

    if pre_norm:
        normalized_hgvs_genomic = pre_norm
    else:
        normalized_hgvs_genomic = normalizer.normalize(hgvs_genomic)

    if hgvs_genomic.type != "g":
        normalized_hgvs_genomic = pre_push_vcf_tx_g_map_fix(
            normalized_hgvs_genomic,
            hgvs_genomic,
            genomic_ac,
            vm,
            hn,
            sf,
            mapped_g,
        )

    if hgvs_genomic.type == "g":
        chrom = seq_data.get_chr_num_ucsc(
            normalized_hgvs_genomic.ac,
            primary_assembly,
        ) or normalized_hgvs_genomic.ac
    else:
        chrom = normalized_hgvs_genomic.ac

    pos, ref, alt = _hgvs_vcf_sequence(normalized_hgvs_genomic, sf)
    if not (pos and ref and alt):
        chrom = ""

    return (
        hgvs_genomic,
        normalized_hgvs_genomic,
        chrom,
        pos,
        ref,
        alt,
    )


def _hard_exon_boundary(
        map_dat,
        tx_ac,
        hgvs_ac,
        genomic_ac,
        alt_aln_method,
        pos,
        direction,
):
    """Return the exon boundary used by a hard push.

    ``direction`` is ``right`` for the 3-prime boundary and ``left`` for
    the 5-prime boundary. The mapping column selection is unchanged from
    the original hard push implementations.
    """
    if genomic_ac is False:
        exon_set = map_dat.mapped_exons(
            tx_ac,
            hgvs_ac,
            alt_aln_method=alt_aln_method,
        )
        start_column, end_column = 7, 8
    else:
        exon_set = map_dat.mapped_exons(
            hgvs_ac,
            genomic_ac,
            alt_aln_method=alt_aln_method,
        )
        start_column, end_column = 5, 6

    for exon in exon_set:
        if int(exon[start_column]) + 1 <= int(pos) <= int(exon[end_column]):
            if direction == "right":
                return int(exon[end_column])
            return int(exon[start_column] + 1)

    return None


def hard_right_hgvs2vcf(hgvs_genomic, primary_assembly, hn, reverse_normalizer, sf, tx_ac, map_dat, alt_aln_method, hp, vm,
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
    :param map_dat: cached fetcher/store for transcript mapping data
    :param alt_aln_method:
    :param hp:
    :param vm:
    :param mrg:
    :param genomic_ac: Genomic ac when transcirpt var is input *Must* be false for genomic var
    :return:
    """
    (
        hgvs_genomic,
        normalized_hgvs_genomic,
        chr,
        pos,
        ref,
        alt,
    ) = _prepare_hard_hgvs(
        hgvs_genomic,
        primary_assembly,
        hn,
        hn,
        sf,
        vm,
        tx_ac,
        map_dat,
        alt_aln_method,
        genomic_ac,
        mapped_g,
        pre_norm,
    )

    # Add surrounding bases
    # If possible, capture and alt variant that spans the gap
    merged_variant = False
    pre_merged_variant = False
    identifying_variant = False
    identifying_g_variant = False
    needs_a_push = False

    if chr != "" and pos != "" and ref != "" and alt != "":

        # Set exon boundary.
        exon_end_genomic = _hard_exon_boundary(
            map_dat,
            tx_ac,
            hgvs_genomic.ac,
            genomic_ac,
            alt_aln_method,
            pos,
            "right",
        )

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
            flank_seq = sf.fetch_seq(normalized_hgvs_genomic.ac, working_pos - 1, working_pos + max_push_length)
        except  vvhgvs.exceptions.HGVSDataNotAvailableError as e:
            if "ValueError: stop out of range" in str(e):
                flank_seq = False
                # this means that we went beyond the end of the seq this should be very rare but is possible
                # fall back to old behaviour here
            else:
                raise e
        for push in range(max_push_length):
            try:
                if flank_seq:
                    push_ref = push_ref + flank_seq[push]
                    push_alt = push_alt + flank_seq[push]
                else:
                    post = sf.fetch_seq(normalized_hgvs_genomic.ac, working_pos - 1, working_pos)
                    push_ref = push_ref + post
                    push_alt = push_alt + post
            except IndexError:
                needs_a_push = False
                break

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
                    normlize_check_mapped = vm.g_to_n(normlize_check_variant,
                                                      tx_ac, alt_aln_method)

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
            if not normlize_check_mapped.posedit.edit.ref or len(normlize_check_mapped.posedit.edit.ref) <= 1:
                staging_loop = staging_loop + 1

            # exon boundary hit. Break before intron
            if working_pos > exon_end_genomic:
                needs_a_push = False
                break

            # Check here for the gap (Has it been crossed?) Note: if gap in tx, we have the whole gap spanned
            elif (normlize_check_mapped.posedit.edit.ref and
                  ((len(normlize_check_mapped.posedit.edit.ref) != len(normlize_check_variant.posedit.edit.ref) and
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
                                    pre_merged_variant = mrg([v1, v2], reverse_normalizer, final_norm=False, map_dat=map_dat)
                                else:
                                    pre_merged_variant = mrg([v2, v1], reverse_normalizer, final_norm=False, map_dat=map_dat)
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
                                    if hgvs_position_utils.either_position_is_intronic(merged_variant):
                                        # Try from normalized genomic
                                        pre_merged_variant = hn.normalize(pre_merged_variant)
                                        test_merged_variant = vm.g_to_n(
                                            pre_merged_variant,
                                            tx_ac
                                        )

                                        if not hgvs_position_utils.either_position_is_intronic(
                                                test_merged_variant
                                        ):
                                            merged_variant = pre_merged_variant
                                        else:
                                            pre_merged_variant = reverse_normalizer.normalize(
                                                pre_merged_variant
                                            )
                                            test_merged_variant = vm.g_to_n(
                                                pre_merged_variant,
                                                tx_ac
                                            )

                                            if not hgvs_position_utils.either_position_is_intronic(
                                                    test_merged_variant
                                            ):
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
                                    pre_merged_variant = mrg([v1, v2], reverse_normalizer, final_norm=False, map_dat=map_dat)
                                else:
                                    pre_merged_variant = mrg([v2, v1], reverse_normalizer, final_norm=False, map_dat=map_dat)
                                if "g" in pre_merged_variant.type:
                                    merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                else:
                                    merged_variant = pre_merged_variant
                            except utils.mergeHGVSerror:
                                try:
                                    if v1.posedit.pos.start.base < v2.posedit.pos.start.base:
                                        pre_merged_variant = mrg([v1, v2], hn, final_norm=False, map_dat=map_dat)
                                    else:
                                        pre_merged_variant = mrg([v2, v1], hn, final_norm=False, map_dat=map_dat)
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
                                    if hgvs_position_utils.either_position_is_intronic(merged_variant):
                                        # Try from normalized genomic
                                        try:
                                            pre_merged_variant = hn.normalize(pre_merged_variant)
                                        except vvhgvs.exceptions.HGVSError:
                                            pass

                                        test_merged_variant = vm.g_to_n(
                                            pre_merged_variant,
                                            tx_ac
                                        )

                                        if not hgvs_position_utils.either_position_is_intronic(
                                                test_merged_variant
                                        ):
                                            merged_variant = pre_merged_variant
                                        else:
                                            pre_merged_variant = reverse_normalizer.normalize(
                                                pre_merged_variant
                                            )
                                            test_merged_variant = vm.g_to_n(
                                                pre_merged_variant,
                                                tx_ac
                                            )

                                            if not hgvs_position_utils.either_position_is_intronic(
                                                    test_merged_variant
                                            ):
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

    # Dictionary VCF
    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': normalized_hgvs_genomic,
                'merged_variant': merged_variant, 'identifying_variant': identifying_variant,
                'pre_merged_variant': pre_merged_variant, 'identifying_g_variant': identifying_g_variant}
    str_hgvs = vcfcp_to_hgvsstr(vcf_dict, hgvs_genomic)
    vcf_dict['str_hgvs'] = str_hgvs
    vcf_dict['needs_a_push'] = needs_a_push
    return vcf_dict # Return dict

def hard_left_hgvs2vcf(hgvs_genomic, primary_assembly, hn, reverse_normalizer, sf, tx_ac, map_dat, alt_aln_method,
                       hp, vm, mrg, genomic_ac=False, mapped_g=False, pre_norm=False):
    """
    Designed specifically for gap handling - hard left pushes as 5 prime as possible and adds additional bases
    :param hgvs_genomic
    :param primary_assembly
    :param hn
    :param reverse_normalizer
    :param sf
    :param tx_ac
    :param map_dat: cached fetcher/store for transcript mapping data
    :param alt_aln_method
    :param hp
    :param vm
    :param mrg
    :param genomic_ac
    :param mapped_g: genomic mapping, used if a transcript type is input
    :return
    """
    (
        hgvs_genomic,
        reverse_normalized_hgvs_genomic,
        chr,
        pos,
        ref,
        alt,
    ) = _prepare_hard_hgvs(
        hgvs_genomic,
        primary_assembly,
        reverse_normalizer,
        hn,
        sf,
        vm,
        tx_ac,
        map_dat,
        alt_aln_method,
        genomic_ac,
        mapped_g,
        pre_norm,
    )

    # Add surrounding bases
    # If possible, capture and alt variant that spans the gap
    merged_variant = False
    pre_merged_variant = False
    identifying_variant = False
    identifying_g_variant = False
    needs_a_push = False

    if (chr != ''
            and pos != ''
            and ref != ''
            and alt != ''):
        # Set exon boundary.
        exon_start_genomic = _hard_exon_boundary(
            map_dat,
            tx_ac,
            hgvs_genomic.ac,
            genomic_ac,
            alt_aln_method,
            pos,
            "left",
        )

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
                    normlize_check_mapped = vm.g_to_n(normlize_check_variant,
                                                      tx_ac, alt_aln_method)
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

            if not normlize_check_mapped.posedit.edit.ref or len(normlize_check_mapped.posedit.edit.ref) <= 1:
                staging_loop = staging_loop + 1

            # Check here for the gap (Has it been crossed?) Note: if gap in tx, we have the whole gap spanned
            #
            if ((normlize_check_mapped.posedit.edit.ref and
                (len(normlize_check_mapped.posedit.edit.ref) != len(normlize_check_variant.posedit.edit.ref) and
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
                                    pre_merged_variant = mrg([v1, v2], reverse_normalizer, final_norm=False, map_dat=map_dat)
                                else:
                                    pre_merged_variant = mrg([v2, v1], reverse_normalizer, final_norm=False, map_dat=map_dat)
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
                                    if hgvs_position_utils.either_position_is_intronic(merged_variant):
                                        # Try from normalized genomic
                                        try:
                                            pre_merged_variant = hn.normalize(pre_merged_variant)
                                        except vvhgvs.exceptions.HGVSError:
                                            pass

                                        test_merged_variant = vm.g_to_n(
                                            pre_merged_variant,
                                            tx_ac
                                        )

                                        if not hgvs_position_utils.either_position_is_intronic(
                                                test_merged_variant
                                        ):
                                            merged_variant = test_merged_variant
                                        else:
                                            pre_merged_variant = reverse_normalizer.normalize(pre_merged_variant)
                                            test_merged_variant = vm.g_to_n(pre_merged_variant, tx_ac)
                                            if not hgvs_position_utils.either_position_is_intronic(
                                                    test_merged_variant
                                            ):
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
                        if hgvs_position_utils.both_positions_are_intronic(map_back):
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
                                    pre_merged_variant = mrg([v1, v2], reverse_normalizer, final_norm=False, map_dat=map_dat)
                                else:
                                    pre_merged_variant = mrg([v2, v1], reverse_normalizer, final_norm=False, map_dat=map_dat)
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
                                    if hgvs_position_utils.either_position_is_intronic(merged_variant):
                                        # Try from normalized genomic
                                        pre_merged_variant = hn.normalize(pre_merged_variant)
                                        test_merged_variant = vm.g_to_n(
                                            pre_merged_variant,
                                            tx_ac
                                        )

                                        if not hgvs_position_utils.either_position_is_intronic(
                                                test_merged_variant
                                        ):
                                            merged_variant = pre_merged_variant
                                        else:
                                            pre_merged_variant = reverse_normalizer.normalize(
                                                pre_merged_variant
                                            )
                                            test_merged_variant = vm.g_to_n(
                                                pre_merged_variant,
                                                tx_ac
                                            )

                                            if not hgvs_position_utils.either_position_is_intronic(
                                                    test_merged_variant
                                            ):
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

    # Dictionary VCF
    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': reverse_normalized_hgvs_genomic,
                'merged_variant': merged_variant, 'identifying_variant': identifying_variant,
                'pre_merged_variant': pre_merged_variant, 'identifying_g_variant': identifying_g_variant}
    str_hgvs = vcfcp_to_hgvsstr(vcf_dict, hgvs_genomic)
    vcf_dict['str_hgvs'] = str_hgvs
    vcf_dict['needs_a_push'] = needs_a_push
    return vcf_dict # Return dict

def hgvs_ref_alt(hgvs_variant,
                 sf):
    edit = hgvs_variant.posedit.edit
    edit_type = edit.type

    # Identity
    if edit_type == 'identity':
        ref = edit.ref
        alt = edit.ref
    # Ins
    elif edit_type == 'ins':
        end = hgvs_variant.posedit.pos.end.base
        start = hgvs_variant.posedit.pos.start.base
        alt_start = start - 1

        # Recover sequence
        ref_seq = sf.fetch_seq(hgvs_variant.ac, alt_start, end)
        ins_seq = edit.alt

        # Assemble vcf
        ref = ref_seq # stays equivalent
        alt = (ref_seq[:1] +
               ins_seq +
               ref_seq[-1:])
    # Subs
    elif edit_type == 'sub':
        ref = edit.ref
        alt = edit.alt
    # Dels
    elif edit_type == 'del':
        ref = edit.ref
        alt = ""
    # Invs
    elif edit_type == 'inv':
        ref = edit.ref
        alt = utils.simple_dna_revcomp(ref)
    # Delins variants
    elif edit_type == 'delins':
        ref = edit.ref
        alt = edit.alt
    # Dups
    elif edit_type == 'dup':
        ref = edit.ref
        alt = edit.ref + edit.ref

    else: # Not defined
        ref = ""
        alt = ""
    return {'ref': ref, 'alt': alt}

def incomplete_alignment_mapping_t_to_g(validator, variant):
    output = None
    mapping_options = variant.map_dat.mapping_options(variant.input_parses.ac,hdp=validator.hdp)
    for option in mapping_options:
        if option[2] == validator.alt_aln_method and "NC_" not in option[1]:
            in_assembly = seq_data.get_chr_num_refseq(option[1], variant.primary_assembly)
            if in_assembly is not None:
                try:
                    output = validator.vm.t_to_g(variant.input_parses, option[1])
                    if variant.input_parses.posedit.edit.type == "identity":
                        output.posedit.edit.alt = output.posedit.edit.ref
                except vvhgvs.exceptions.HGVSError:
                    pass
    return output

# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
