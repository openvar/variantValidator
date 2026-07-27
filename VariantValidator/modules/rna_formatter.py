# -*- coding: utf-8 -*-
import vvhgvs.exceptions
from . import utils


# Custom Exceptions
class RnaVariantSyntaxError(Exception):
    pass


class RnaDescriptions(object):
    # Initialise and add initialisation data to the object
    def __init__(self, alt_aln_method, genome_build, validator, variant):
        """
        Setup the object
        :param alt_aln_method: String - "splign" or "genebuild"
        :param genome_build: String - "GRCh37" or "GRCh38"
        """
        self.input = None
        self.cross_normalizer = variant.cross_hn
        self.evm = variant.evm
        self.normalizer = variant.hn
        self.protein_variant = None
        self.protein_variant_slr = None
        self.rna_variant = None
        self.dna_variant = None
        self.is_a_prediction = False
        self.translate = validator.myc_to_p
        self.parse = validator.hp.parse
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
        # myc_to_p(self, hgvs_transcript, evm, re_to_p, hn):
        translation = self.translate(hgvs_variant, self.evm, re_to_p=False, hn=self.normalizer)
        protein_variant = translation["hgvs_protein"]

        if self.is_a_prediction is False:
            protein_variant.posedit.uncertain = False
        else:
            self.input = self.is_a_prediction

        self.protein_variant = protein_variant
        self.protein_variant_slr = utils.single_letter_protein(protein_variant)

    def check_syntax(self, variant_string, variant):
        """
        Checks the syntax of the variant string which is a submitted cRNA description string and performs necessary
        transformations to populate the object
        :param variant_string:
        :param variant:
        :return: None
        """
        # Check the submitted RNA description
        if ":r." not in variant_string:
            raise RnaVariantSyntaxError(
                "VariantSyntaxError: The variant type for an RNA description must be r."
            )

        if any(base in variant_string for base in "GATCU"):
            raise RnaVariantSyntaxError(
                "VariantSyntaxError: RNA sequence must be lower-case"
            )

        if "t" in variant_string.split(":r.", 1)[1]:
            raise RnaVariantSyntaxError(
                "VariantSyntaxError: RNA sequence contains Uracil (u) not Thymine (T)"
            )

        if ":r.(" in variant_string:
            self.is_a_prediction = variant_string
            variant_string = variant_string.replace(":r.(", ":r.", 1)
            variant_string = variant_string[:-1]

        # Record if passes
        self.input = variant_string

        # Replace and uppercase base sequences
        self.replace_and_upper(self.input)

        # Parse the string into HGVS object
        if "NR_" in self.input or variant.transcript_type == "n":
            raise RnaVariantSyntaxError("VariantSyntaxError: Invalid variant type for non-coding transcript. "
                                        "Instead use n.")

        hgvs_rna = self.parse(self.dna_variant)

        # Normalize
        try:
            hgvs_rna = self.cross_normalizer.normalize(hgvs_rna)
        except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
            if "Normalization of intronic variants is not supported" in str(e):
                raise RnaVariantSyntaxError(
                    "VariantSyntaxError: Intronic descriptions are only valid in the context "
                    "of a c. description"
                )
            raise

        self.dna_variant = hgvs_rna

        # Translate
        self.translate_from_rna_variant(hgvs_rna)

        # Convert back to RNA
        self.replace_and_lower(utils.valstr(hgvs_rna))

    def output_dict(self):
        """
        Returns a dictionary output of the outputs from the rna submission
        :return: Dict
        """
        variant_dict = {"usage_warnings": self.usage_warnings,
                        "rna_variant": self.rna_variant,
                        "translation": str(self.protein_variant),
                        "translation_slr": self.protein_variant_slr}
        return variant_dict

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
