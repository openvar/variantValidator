import re
from vvhgvs.assemblymapper import AssemblyMapper
import vvhgvs.exceptions


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


def fuzzy_ends(my_variant):
    """
    :param my_variant:
    :return: False if no fuzzy end detected otherwise raises Exception and provides information as to where the
    fuzzy end is located
    """

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
    if ")_(" in my_variant.quibble:
        accession_and_type, positions_and_edit = my_variant.quibble.split(".(")
        position_1, position_2 = positions_and_edit.split(")_(")
        position_2, variation = position_2.split(")")
        position_1 = position_1.replace(")", "")
        position_2 = position_2.replace("(", "")
        v1 = f"{accession_and_type}.{position_1}{variation}"
        v2 = f"{accession_and_type}.{position_2}{variation}"
        my_variant.reftype = f":{accession_and_type.split(':')[1]}."
        try:
            parsed_v1 = validator.hp.parse(v1)
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
                                        f"accession {accession_and_type.split(':')[0]}")
        try:
            parsed_v2 = validator.hp.parse(v2)
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
                                        f"accession {accession_and_type.split(':')[0]}")

        # Check positions are in the correct order
        if parsed_v1.posedit.pos.end.base >= parsed_v2.posedit.pos.start.base:
            raise InvalidRangeError(f"Position {parsed_v1.posedit.pos} is > or overlaps {parsed_v2.posedit.pos}")

        # Mark as a reformat output
        my_variant.warnings.append("Uncertain positions are not fully supported, however the syntax is valid")
        my_variant.reformat_output = "uncertain_pos"

        # Genomic Variants
        if "NC_" in my_variant.quibble:
            my_variant.hgvs_genomic = my_variant.quibble
            start_pos = position_1.split("_")[0]
            end_pos = position_2.split("_")[1]
            v3 = f"{accession_and_type}.{start_pos}_{end_pos}{variation}"
            parsed_v3 = validator.hp.parse(v3)
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
                    t_position_1 = str(ptv1.posedit.pos)
                    t_position_2 = str(ptv2.posedit.pos)
                else:
                    t_position_1 = str(ptv2.posedit.pos)
                    t_position_2 = str(ptv1.posedit.pos)
                tx_variant = f"{ptv1.ac}:{ptv1.type}.({t_position_1})_({t_position_2}){ptv1.posedit.edit.type}"
                my_variant.hgvs_coding = tx_variant
                my_variant.hgvs_transcript_variant = tx_variant
            else:
                my_variant.warnings.append("Selected transcript does not span the entire range "
                                           "of the genomic variation")

        # Transcript Variants
        elif "NM_" in my_variant.quibble or "NR_" in my_variant.quibble or "ENST" in my_variant.quibble:
            pgv1 = evm.t_to_g(parsed_v1)
            pgv2 = evm.t_to_g(parsed_v2)
            if pgv1.posedit.pos.start.base < pgv2.posedit.pos.start.base:
                g_position_1 = str(pgv1.posedit.pos)
                g_position_2 = str(pgv2.posedit.pos)
            else:
                g_position_1 = str(pgv2.posedit.pos)
                g_position_2 = str(pgv1.posedit.pos)
            gen_variant = f"{pgv1.ac}:{pgv1.type}.({g_position_1})_({g_position_2}){pgv1.posedit.edit.type}"
            my_variant.hgvs_genomic = gen_variant
            my_variant.hgvs_coding = my_variant.quibble
            my_variant.hgvs_transcript_variant = my_variant.quibble

    else:
        accession_and_type, position_and_edit = my_variant.quibble.split(".(")
        position_1, variation = position_and_edit.split(")")
        v1 = f"{accession_and_type}.{position_1}="
        my_variant.reftype = f":{accession_and_type.split(':')[1]}."
        try:
            parsed_v1 = validator.hp.parse(v1)
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
                                        f"accession {accession_and_type.split(':')[0]}")

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

                tx_variant = f"{ptv1[0].ac}:{ptv1[0].type}.({ptv1[0].posedit.pos}){variation}"
                my_variant.hgvs_coding = tx_variant
                my_variant.hgvs_transcript_variant = tx_variant
            else:
                my_variant.warnings.append("Selected transcript does not span the entire range "
                                           "of the genomic variation")

        elif "NM_" in my_variant.quibble or "NR_" in my_variant.quibble or "ENST" in my_variant.quibble:
            pgv1 = evm.t_to_g(parsed_v1)
            gen_variant = f"{pgv1.ac}:{pgv1.type}.({pgv1.posedit.pos}){variation}"
            my_variant.hgvs_genomic = gen_variant
            my_variant.hgvs_coding = my_variant.quibble
            my_variant.hgvs_transcript_variant = my_variant.quibble
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
