import re
import copy
from vvhgvs.assemblymapper import AssemblyMapper
from VariantValidator.modules import utils as vv_utils
from VariantValidator.modules.hgvs_utils import hgvs_obj_from_existing_edit
import vvhgvs.exceptions
from vvhgvs.enums import Datum
from vvhgvs.location import Interval
# fuzzy but specified ended intervals, take a normal interval for both ends
# can take a pair of Intervals or BaseOffsetIntervals
class FEInterval(Interval):
    uncertain = True
    pass
    def validate(self):
        if self.start:
            (res, msg) = self.start.validate()
            if res != ValidationLevel.VALID:
                return (res, msg)
        if self.end:
            (res, msg) = self.end.validate()
            if res != ValidationLevel.VALID:
                return (res, msg)
        # Check start less than or equal to end
        # for now overlap is allowed so long as some is distinct
        if not self.start or not self.end:
            return (ValidationLevel.VALID, None)
        try:
            if self.start.start <= self.end.end:
                return (ValidationLevel.VALID, None)
            else:
                return (ValidationLevel.ERROR, "base start position must be <= end position")
        except HGVSUnsupportedOperationError as err:
            return (ValidationLevel.WARNING, str(err))

    def format(self, conf=None):
        if self.start is None:
            return ""
        if self.end is None or self.start == self.end:
            return self.start.format(conf)
        # always uncertain
        iv = "(" +self.start.format(conf) + ")_(" + self.end.format(conf) + ")"
        return iv



# Errors
class FuzzyPositionError(Exception):
    pass


class FuzzyRangeError(Exception):
    pass


class InvalidRangeError(Exception):
    pass


class IncompatibleTypeError(Exception):
    pass


class HgvsParseError(Exception):
    pass


def fuzzy_ends(my_variant, validator):
    """
    :param my_variant:
    :return: False if no fuzzy end detected otherwise raises Exception and provides information as to where the
    fuzzy end is located
    """
    if re.match(r"(NM_|NR_|NC_|NG_|ENST)\d+\.\d+:(g|c)\.\(\?\_\d+(\+|-)?\d*\)\_\(\d+(\+|-)?\d*_\?\)"
                r"(del|dup|inv)$", my_variant.quibble):

        parts = my_variant.quibble.split(")_(")
        num1 = parts[0].split("_")[-1]
        num2 = parts[1].split("_")[0]

        intronic_positions = []
        try:
            num1 = int(num1)
        except ValueError:
            num1 = int(re.split(r'(\+|-)', num1)[0])
            intronic_positions.append(num1)
        try:
            num2 = int(num2)
        except ValueError:
            num2 = int(re.split(r'(\+|-)', num2)[0])
            intronic_positions.append(num2)

        if re.search(r":[cnr]\.", my_variant.quibble):
            exon_boundaries = vv_utils.get_exon_boundary_list(my_variant, validator)
            for pos in intronic_positions:
                if str(pos) not in exon_boundaries:
                    raise FuzzyPositionError(f"ExonBoundaryError: Position {pos} does not correspond to an exon "
                                             f"boundary for transcript {my_variant.quibble.split(':')[0]} aligned"
                                             f" to {my_variant.primary_assembly} genomic reference "
                                             f"sequence {exon_boundaries[1]}")

        if int(num1) >= int(num2):
            my_variant.warnings.append("Uncertain positions are not fully supported, however the start position is > "
                                       "the end position")
        else:
            my_variant.warnings.append("Uncertain positions are not fully supported, however the syntax is valid")
        raise FuzzyRangeError("Fuzzy/unknown variant start and end positions "
                              "in submitted variant description")

    elif re.match(r"(NM_|NR_|NC_|NG_|ENST)\d+\.\d+:(g|c)\.\(\d+(\+|-)?\d*_\d+(\+|-)?\d*\)\_\(\d+(\+|-)?\d*_\?\)"
                  r"(del|dup|inv)$", my_variant.quibble):

        # Split the string at ")_(" and "_"
        parts = my_variant.quibble.split(")_(")
        parts[0] = parts[0].split("(")[1]
        num1, num2 = parts[0].split("_")
        num3 = parts[1].split("_")[0]

        # Convert the numbers to integers
        intronic_positions = []
        try:
            num1 = int(num1)
        except ValueError:
            num1 = int(re.split(r'(\+|-)', num1)[0])
            intronic_positions.append(num1)
        try:
            num2 = int(num2)
        except ValueError:
            num2 = int(re.split(r'(\+|-)', num2)[0])
            intronic_positions.append(num2)
        try:
            num3 = int(num3)
        except ValueError:
            num3 = int(re.split(r'(\+|-)', num3)[0])
            intronic_positions.append(num3)

        if re.search(r":[cnr]\.", my_variant.quibble):
            exon_boundaries = vv_utils.get_exon_boundary_list(my_variant, validator)
            for pos in intronic_positions:
                if str(pos) not in exon_boundaries:
                    raise FuzzyPositionError(f"ExonBoundaryError: Position {pos} does not correspond to an exon "
                                             f"boundary for transcript {my_variant.quibble.split(':')[0]} aligned"
                                             f" to {my_variant.primary_assembly} genomic reference "
                                             f"sequence {exon_boundaries[1]}")

        # Check if the numbers are in order
        if num1 <= num2 <= num3:
            my_variant.warnings.append("Uncertain positions are not fully supported, however the syntax is valid")
        else:
            my_variant.warnings.append("Uncertain positions are not fully supported, however the provided positions "
                                       "are out of order")
        raise FuzzyRangeError("Fuzzy/unknown variant end position in submitted variant description")

    elif re.match(
            r"(NM_|NR_|NC_|NG_|ENST)\d+\.\d+:(g|c)\.\(\?\_\d+(\+|-)?\d*\)\_\(\d+(\+|-)?\d*_\d+(\+|-)?\d*\)"
            r"(del|dup|inv)$", my_variant.quibble):

        # Split the string at ")_(" and "_"
        parts = my_variant.quibble.split(")_(")
        # Extract the first part of the string after "(_?" and split at "_"
        num1 = parts[0].split("_")[-1]
        # Extract the second part and split at "_"
        num2, num3 = parts[1].split(")_")[0].split("_")
        num3 = num3.split(")")[0]

        # Convert the numbers to integers
        intronic_positions = []
        try:
            num1 = int(num1)
        except ValueError:
            num1 = int(re.split(r'(\+|-)', num1)[0])
            intronic_positions.append(num1)
        try:
            num2 = int(num2)
        except ValueError:
            num2 = int(re.split(r'(\+|-)', num2)[0])
            intronic_positions.append(num2)
        try:
            num3 = int(num3)
        except ValueError:
            num3 = int(re.split(r'(\+|-)', num3)[0])
            intronic_positions.append(num3)

        if re.search(r":[cnr]\.", my_variant.quibble):
            exon_boundaries = vv_utils.get_exon_boundary_list(my_variant, validator)
            for pos in intronic_positions:
                if str(pos) not in exon_boundaries:
                    raise FuzzyPositionError(f"ExonBoundaryError: Position {pos} does not correspond to an exon "
                                             f"boundary for transcript {my_variant.quibble.split(':')[0]} aligned"
                                             f" to {my_variant.primary_assembly} genomic reference "
                                             f"sequence {exon_boundaries[1]}")

        # Check if the numbers are in order
        if num1 <= num2 <= num3:
            my_variant.warnings.append("Uncertain positions are not fully supported, however the syntax is valid")
        else:
            my_variant.warnings.append("Uncertain positions are not fully supported, however the provided positions "
                                       "are out of order")
        raise FuzzyRangeError("Fuzzy/unknown variant end position in submitted variant description")

    else:
        if "?" in str(my_variant.hgvs_formatted.posedit.pos):
            if "?" in str(my_variant.hgvs_formatted.posedit.pos.end) and "?" not in str(
                    my_variant.hgvs_formatted.posedit.pos.start):
                raise FuzzyPositionError("Fuzzy/unknown variant end position in submitted variant description")
            elif "?" in str(my_variant.hgvs_formatted.posedit.pos.start) and "?" not in str(
                    my_variant.hgvs_formatted.posedit.pos.end):
                raise FuzzyPositionError("Fuzzy/unknown variant start position in submitted variant description")
            else:
                raise FuzzyPositionError("Fuzzy/unknown variant start and end positions "
                                         "in submitted variant description")

    return False


def uncertain_positions(my_variant, validator):
    # Create evm
    evm = AssemblyMapper(validator.hdp,
                         assembly_name=my_variant.primary_assembly,
                         alt_aln_method=validator.alt_aln_method,
                         normalize=False,
                         replace_reference=True)

    # Check for uncertain positions in the correct place
    if not re.search("[gcnr].\(", my_variant.quibble):
        return

    # Formats like NC_000005.9:g.(90136803_90144453)_(90159675_90261231)dup
    if ")_(" in my_variant.quibble and "?" not in my_variant.quibble:
        accession, _sep, var_type_and_posedit = my_variant.quibble.partition(':')
        var_type, _sep ,positions_and_edit = var_type_and_posedit.partition(".(")
        position_1, _sep, position_2 = positions_and_edit.partition(")_(")
        position_2, _sep, variation = position_2.partition(")")
        try:
            if var_type == 'p':
                edit = validator.hp.parse_pro_edit(variation)
            elif var_type in ['c','t','r']:
                edit = validator.hp.parse_rna_edit(variation)
            else:
                edit = validator.hp.parse_dna_edit(variation)
        except vvhgvs.exceptions.HGVSError as e:
            raise HgvsParseError(str(e))
        position_1 = position_1.replace(")", "")
        position_2 = position_2.replace("(", "")

        my_variant.reftype = f":{var_type}."

        try:
            start, _sep, end = position_1.partition('_')
            parsed_v1 = hgvs_obj_from_existing_edit(
                    accession,
                    var_type,
                    start,
                    copy.copy(edit),
                    end=end)
        except vvhgvs.exceptions.HGVSError as e:
            raise HgvsParseError(str(e))
        try:
            validator.vr.validate(parsed_v1)
        except vvhgvs.exceptions.HGVSError as e:
            if "is not known to be compatible with variant type" in str(e):
                raise IncompatibleTypeError(str(e))
            elif "base start position must be <= end position" in str(e):
                raise InvalidRangeError(f"{str(e)} in position {str(parsed_v1.posedit.pos)}")
            elif re.search("[+-]", str(parsed_v1.posedit.pos)) or re.search("[+-]", str(parsed_v1.posedit.pos)):
                pass
            else:
                raise InvalidRangeError(f"{position_1} is an invlaid range for "
                                        f"accession {accession}")
        try:
            start, _sep, end = position_2.partition('_')
            parsed_v2 = hgvs_obj_from_existing_edit(
                    accession,
                    var_type,
                    start,
                    copy.copy(edit),
                    end=end)
        except vvhgvs.exceptions.HGVSError as e:
            raise HgvsParseError(str(e))
        try:
            validator.vr.validate(parsed_v2)
        except vvhgvs.exceptions.HGVSError as e:
            if "is not known to be compatible with variant type" in str(e):
                raise IncompatibleTypeError(str(e))
            elif "base start position must be <= end position" in str(e):
                raise InvalidRangeError(f"{str(e)} in position {str(parsed_v2.posedit.pos)}")
            elif re.search("[+-]", str(parsed_v2.posedit.pos)) or re.search("[+-]", str(parsed_v2.posedit.pos)):
                pass
            else:
                raise InvalidRangeError(f"{position_2} is an invlaid range for "
                                        f"accession {accession}")

        # Check positions are in the correct order
        if parsed_v1.posedit.pos.end.base >= parsed_v2.posedit.pos.start.base:
            raise InvalidRangeError(f"Position {parsed_v1.posedit.pos} is > or overlaps {parsed_v2.posedit.pos}")
        #TODO check del length works with range given i.e del > min and del < max
        # Mark as a reformat output
        my_variant.warnings.append("Uncertain positions are not fully supported, however the syntax is valid")
        my_variant.reformat_output = "uncertain_pos"

        # Genomic Variants
        if "NC_" in my_variant.quibble:
            my_variant.hgvs_genomic = my_variant.quibble
            start_pos = position_1.split("_")[0]
            end_pos = position_2.split("_")[1]
            parsed_v3 = hgvs_obj_from_existing_edit(
                    accession,
                    var_type,
                    start_pos,
                    edit,
                    end=end_pos)
            validator.vr.validate(parsed_v3)
            if (validator.select_transcripts != "select" and ("NM_" in validator.select_transcripts or
                                                              "ENST" in validator.select_transcripts) and
                    ("|" not in validator.select_transcripts and "[" not in validator.select_transcripts)):
                pass
            elif validator.select_transcripts != "select":
                validator.select_transcripts = "select"
                my_variant.warnings.append("Only a single transcript can be processed, updating to Select")

            # Get transcripts
            rel_var = validator.relevant_transcripts(parsed_v3, evm, validator.alt_aln_method,
                                                     my_variant.reverse_normalizer, validator.select_transcripts)
            if len(rel_var) != 0:
                # Filter for Select transcripts only if transcript not stated
                my_variant.output_type_flag = "gene"
                if validator.select_transcripts == "select":
                    for variant in rel_var:
                        annotation = validator.db.get_transcript_annotation(variant.ac)
                        if '"select": "MANE"' in annotation:
                            validator.select_transcripts = variant.ac
                            break
                        elif '"select": "RefSeq"' in annotation or \
                                '"select": "Ensembl"' in annotation:
                            validator.select_transcripts = variant.ac
                            continue

                # Map uncertain positions to transcript
                ptv1 = validator.relevant_transcripts(parsed_v1, evm, validator.alt_aln_method,
                                                      my_variant.reverse_normalizer, validator.select_transcripts)
                ptv2 = validator.relevant_transcripts(parsed_v2, evm, validator.alt_aln_method,
                                                      my_variant.reverse_normalizer, validator.select_transcripts)

                # Filter correct transcript
                ptv1 = [e for e in ptv1 if validator.select_transcripts in str(e)][0]
                ptv2 = [e for e in ptv2 if validator.select_transcripts in str(e)][0]

                if ptv1.posedit.pos.start.base < ptv2.posedit.pos.start.base:
                    t_position_1 = ptv1.posedit.pos
                    t_position_2 = ptv2.posedit.pos
                else:
                    t_position_1 = ptv2.posedit.pos
                    t_position_2 = ptv1.posedit.pos
                tx_variant = hgvs_obj_from_existing_edit(
                    ptv1.ac,
                    ptv1.type,
                    FEInterval(
                        start=t_position_1,
                        end=t_position_2),
                    edit)

                my_variant.hgvs_coding = tx_variant
                my_variant.quibble = copy.copy(tx_variant)
                my_variant.hgvs_transcript_variant = tx_variant
            else:
                my_variant.output_type_flag = "intergenic"
                my_variant.quibble = copy.copy(parsed_v3)
                my_variant.warnings.append("Selected transcript does not span the entire range "
                                           "of the genomic variation")

        # Transcript Variants
        elif "NM_" in my_variant.quibble or "NR_" in my_variant.quibble or "ENST" in my_variant.quibble:
            pgv1 = evm.t_to_g(parsed_v1)
            pgv2 = evm.t_to_g(parsed_v2)
            if pgv1.posedit.pos.start.base < pgv2.posedit.pos.start.base:
                g_position_1 = pgv1.posedit.pos
                g_position_2 = pgv2.posedit.pos
            else:
                g_position_1 = pgv2.posedit.pos
                g_position_2 = pgv1.posedit.pos
            gen_variant = hgvs_obj_from_existing_edit(
                    pgv1.ac,
                    pgv1.type,
                    FEInterval(
                        start=g_position_1,
                        end=g_position_2),
                    edit)
            my_variant.hgvs_genomic = gen_variant
            my_variant.hgvs_coding = hgvs_obj_from_existing_edit(
                    parsed_v1.ac,
                    parsed_v2.type,
                    FEInterval(start=parsed_v1.posedit.pos,
                               end=parsed_v2.posedit.pos),
                    edit
                    )
            my_variant.hgvs_transcript_variant = copy.copy(my_variant.hgvs_coding)
            my_variant.quibble = copy.copy(my_variant.hgvs_coding)
            my_variant.output_type_flag = "gene"

    elif ")_(" not in my_variant.quibble and not "?" in my_variant.quibble:
        if ")(" in my_variant.quibble:
            raise InvalidRangeError("Invalid range submitted, missing underscore between stated uncertain positions")
        accession, _sep, var_type_and_posedit = my_variant.quibble.partition(':')
        var_type, _sep ,position_and_edit = var_type_and_posedit.partition(".(")
        position_1, variation = position_and_edit.split(")")
        v1 = f"{accession}:{var_type}.{position_1}="
        my_variant.reftype = f":{var_type}."
        try:
            if var_type == 'p':
                edit = validator.hp.parse_pro_edit(variation)
                eq_edit = validator.hp.parse_pro_edit('=')
            elif var_type in ['c','t','r']:
                edit = validator.hp.parse_rna_edit(variation)
                eq_edit = validator.hp.parse_rna_edit('=')
            else:
                edit = validator.hp.parse_dna_edit(variation)
                eq_edit = validator.hp.parse_dna_edit('=')
            start, _sep, end = position_1.partition('_')
            parsed_v1 = hgvs_obj_from_existing_edit(
                    accession,
                    var_type,
                    start,
                    eq_edit,
                    end=end)
        except vvhgvs.exceptions.HGVSError as e:
            raise HgvsParseError(str(e))
        try:
            validator.vr.validate(parsed_v1)
        except vvhgvs.exceptions.HGVSError as e:
            if "is not known to be compatible with variant type" in str(e):
                raise IncompatibleTypeError(str(e))
            elif "base start position must be <= end position" in str(e):
                raise InvalidRangeError(f"{str(e)} in position {str(parsed_v1.posedit.pos)}")
            elif re.search("[+-]", str(parsed_v1.posedit.pos)) or re.search("[+-]", str(parsed_v1.posedit.pos)):
                pass
            else:
                raise InvalidRangeError(f"{position_1} is an invlaid range for "
                                        f"accession {accession}")

        # mark as reformat output
        my_variant.warnings.append("Uncertain positions are not fully supported, however the syntax is valid")
        my_variant.reformat_output = "uncertain_pos"

        # Genomic Variants
        if "NC_" in my_variant.quibble:
            my_variant.hgvs_genomic = my_variant.quibble

            # Make select_transcriopts "select" unless specified
            if (validator.select_transcripts != "select" and ("NM_" in validator.select_transcripts or
                                                              "ENST" in validator.select_transcripts) and
                    ("|" not in validator.select_transcripts and "[" not in validator.select_transcripts)):
                pass
            elif validator.select_transcripts != "select":
                validator.select_transcripts = "select"
                my_variant.warnings.append("Only a single transcript can be processed, updating to Select")

            # Get transcripts
            rel_var = validator.relevant_transcripts(parsed_v1, evm, validator.alt_aln_method,
                                                     my_variant.reverse_normalizer, validator.select_transcripts)

            # Filter transcripts
            if len(rel_var) != 0:
                # Filter for Select transcripts only if transcript not stated
                my_variant.output_type_flag = "gene"
                if validator.select_transcripts == "select":
                    for variant in rel_var:
                        annotation = validator.db.get_transcript_annotation(variant.ac)
                        if '"select": "MANE"' in annotation:
                            validator.select_transcripts = variant.ac
                            break
                        elif '"select": "RefSeq"' in annotation or \
                                '"select": "Ensembl"' in annotation:
                            validator.select_transcripts = variant.ac
                            continue

                # Extract correct description
                ptv1 = [e for e in rel_var if validator.select_transcripts in str(e)]
                tx_variant = hgvs_obj_from_existing_edit(
                    ptv1[0].ac,
                    ptv1[0].type,
                    ptv1[0].posedit.pos,
                    edit)
                tx_variant.posedit.pos.uncertain = True

                my_variant.hgvs_coding = tx_variant
                my_variant.quibble = copy.copy(tx_variant)
                my_variant.hgvs_transcript_variant = tx_variant
            else:
                my_variant.output_type_flag = "intergenic"
                my_variant.quibble = copy.copy(parsed_v1)
                my_variant.warnings.append("Selected transcript does not span the entire range "
                                           "of the genomic variation")

        elif "NM_" in my_variant.quibble or "NR_" in my_variant.quibble or "ENST" in my_variant.quibble:
            pgv1 = evm.t_to_g(parsed_v1)
            gen_variant = hgvs_obj_from_existing_edit(
                    pgv1.ac,
                    pgv1.type,
                    pgv1.posedit.pos,
                    edit)
            gen_variant.posedit.pos.uncertain = True

            my_variant.hgvs_genomic = gen_variant
            my_variant.hgvs_coding = hgvs_obj_from_existing_edit(
                    parsed_v1.ac,
                    parsed_v1.type,
                    parsed_v1.posedit.pos,
                    edit)
            my_variant.hgvs_coding.posedit.pos.uncertain = True
            my_variant.hgvs_transcript_variant = my_variant.hgvs_coding
            my_variant.quibble = copy.copy(my_variant.hgvs_coding)
            my_variant.output_type_flag = "gene"

    else:
        fuzzy_ends(my_variant, validator)
    return

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
