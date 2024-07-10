# -*- coding: utf-8 -*-
import re
import copy
import vvhgvs.exceptions
import vvhgvs.assemblymapper
import vvhgvs.variantmapper
import vvhgvs.normalizer
from . import utils
import VariantFormatter.formatter as formatter  # VariantFormatter has handy lightweight functions that make this easier


# Custom Exceptions
class RnaVariantSyntaxError(Exception):
    pass


class RnaDescriptions(object):
    # Initialise and add initialisation data to the object
    def __init__(self, alt_aln_method, genome_build, validator):
        """
        Setup the object
        :param alt_aln_method: String - "splign" or "genebuild"
        :param genome_build: String - "GRCh37" or "GRCh38"
        """
        self.input = None
        self.vfo = validator
        self.alt_aln_method = alt_aln_method
        self.genome_build = genome_build
        self.cross_normalizer = vvhgvs.normalizer.Normalizer(validator.hdp,
                                                             cross_boundaries=True,
                                                             shuffle_direction=3,
                                                             alt_aln_method=self.alt_aln_method)
        self.protein_variant = None
        self.protein_variant_slr = None
        self.rna_variant = None
        self.dna_variant = None
        self.is_a_prediction = False
        self.translate = formatter.hgvs_transcript2hgvs_protein
        self.parse = formatter.parse
        self.remove_reference = formatter.remove_reference
        self.usage_warnings = ["RNA (r.) descriptions are independent of cDNA descriptions (c.)",
                               "RNA descriptions must only be used if the RNA has been sequenced and must not be "
                               "inferred from a cDNA description",
                               "c. and g. descriptions provided by VariantValidator must only be used if the DNA"
                               " sequence has been confirmed"]

    def replace_and_upper(self, variant_string):
        """
        Take a cRNA formatted variant description and convert into a cDNA description
        :param variant_string:
        :return: None
        """
        reference, edit = variant_string.split(":r.")
        edit = edit.replace("u", "t")
        edit = edit.replace("t", "T")
        edit = edit.replace("g", "G")
        edit = edit.replace("c", "C")
        edit = edit.replace("a", "A")
        edit = edit.replace("dTp", "dup")

        # vv_hgvs requires c. for processing because r. is not handled correctly, so we swap to c.
        self.dna_variant = reference + ":c." + edit

    def replace_and_lower(self, variant_string):
        """
        Take a cDNA formatted variant description and convert into a cRNA description
        :param variant_string:
        :return: None
        """
        if "c." in variant_string:
            reference, edit = variant_string.split(":c.")
        else:
            reference, edit = variant_string.split(":r.")
        edit = edit.replace("T", "u")
        edit = edit.replace("G", "g")
        edit = edit.replace("C", "c")
        edit = edit.replace("A", "a")
        if self.is_a_prediction is not False:
            self.rna_variant = reference + ":r.(" + edit + ")"
        else:
            self.rna_variant = reference + ":r." + edit

    def translate_from_rna_variant(self, hgvs_variant):
        """
        Translates from the cDNA version hgvs_object formatted.
        p. variants from r. descriptions can be taken as literal so we also remove the ()
        :param hgvs_variant: must be c.
        :return: None
        """
        protein_variant = self.translate(hgvs_variant, self.genome_build, self.vfo)
        protein_variant_slr = str(utils.single_letter_protein(protein_variant))
        protein_variant = str(protein_variant)

        if self.is_a_prediction is False:
            protein_variant = protein_variant.replace("(", "")
            protein_variant = protein_variant.replace(")", "")
            protein_variant_slr = protein_variant_slr.replace("(", "")
            protein_variant_slr = protein_variant_slr.replace(")", "")
        else:
            self.input = self.is_a_prediction
        self.protein_variant = protein_variant
        self.protein_variant_slr = protein_variant_slr

    def check_syntax(self, variant_string, variant):
        """
        Checks the syntax of the variant string which is a submitted cRNA description string and performs necessary
        transformations to populate the object
        :param variant_string:
        :return: None
        """
        # Check the content of the object
        if ":r." not in str(variant_string):
            raise RnaVariantSyntaxError("The variant type for an RNA description must be r.")
        if re.search("[GATCU]", variant_string):
            raise RnaVariantSyntaxError("RNA sequence must be lower-case")
        if re.search("t", variant_string.split("r.")[1]):
            raise RnaVariantSyntaxError("RNA sequence contains Uracil (u) not Thymine (T)")
        if re.search(":r.\(", variant_string):
            self.is_a_prediction = copy.copy(variant_string)
            variant_string = variant_string.replace(":r.(", ":r.")
            variant_string = variant_string[:-1]

        # Record if passes
        self.input = variant_string

        # Replace and uppercase base sequences
        self.replace_and_upper(self.input)

        # parse the string into hgvs object
        if "NR_" in self.input or variant.transcript_type == "n":
            raise RnaVariantSyntaxError("Invalid variant type for non-coding transcript. Instead use n.")
        else:
            hgvs_rna = self.parse(self.dna_variant, self.vfo)

        # normalize
        hgvs_rna = self.cross_normalizer.normalize(hgvs_rna)
        self.dna_variant = str(hgvs_rna)

        # translate
        self.translate_from_rna_variant(hgvs_rna)

        # Convert back to rna
        hgvs_rna = self.remove_reference(hgvs_rna)
        self.replace_and_lower(str(hgvs_rna))

    def output_dict(self):
        """
        Returns a dictionary output of the outputs from the rna submission
        :return: Dict
        """
        variant_dict = {"usage_warnings": self.usage_warnings,
                        "rna_variant": self.rna_variant,
                        "translation": str(self.protein_variant),
                        "translation_slr": str(self.protein_variant_slr)}
        return variant_dict

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
