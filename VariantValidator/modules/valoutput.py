import logging
import json

from VariantValidator.modules import lovd_api


logger = logging.getLogger(__name__)


class ValOutput(object):
    """
    Hold final validated Variant objects and provide methods to return
    them in a number of output formats, with or without metadata.
    """

    def __init__(self, outputlist, validator):
        self.output_list = outputlist
        self.validator = validator

    def format_as_dict(self, with_meta=True, test=False):
        validation_output = {'flag': 'warning'}

        validation_error_counter = 0
        validation_obsolete_counter = 0
        validation_warning_counter = 0
        validation_intergenic_counter = 0

        if not self.output_list:
            logger.info("No variants available to output")
            validation_output['flag'] = 'empty_result'

        for variant in self.output_list:
            output_type = variant.output_type_flag

            if output_type == 'gene':
                validation_output['flag'] = 'gene_variant'

                if variant.warnings == ['Validation error']:
                    validation_error_counter += 1
                    identification_key = (
                        f'validation_error_{validation_error_counter}'
                    )

                elif (
                        variant.is_obsolete()
                        and variant.hgvs_transcript_variant == ''
                ):
                    validation_obsolete_counter += 1
                    identification_key = (
                        f'obsolete_record_{validation_obsolete_counter}'
                    )

                else:
                    identification_key = str(
                        variant.hgvs_transcript_variant
                    )

                validation_output[identification_key] = (
                    variant.output_dict(test=test)
                )

            elif output_type == 'warning':
                if variant.warnings == ['Validation error']:
                    validation_error_counter += 1
                    identification_key = (
                        f'validation_error_{validation_error_counter}'
                    )

                elif variant.is_obsolete():
                    validation_obsolete_counter += 1
                    identification_key = (
                        f'obsolete_record_{validation_obsolete_counter}'
                    )

                else:
                    validation_warning_counter += 1

                    self.lovd_syntax_check(variant)

                    identification_key = (
                        f'validation_warning_{validation_warning_counter}'
                    )

                validation_output[identification_key] = (
                    variant.output_dict(test=test)
                )

            elif output_type == 'mitochondrial':
                validation_output['flag'] = 'mitochondrial'

                if variant.warnings == ['Validation error']:
                    validation_error_counter += 1
                    identification_key = (
                        f'validation_error_{validation_error_counter}'
                    )

                elif variant.is_obsolete():
                    validation_obsolete_counter += 1
                    identification_key = (
                        f'obsolete_record_{validation_obsolete_counter}'
                    )

                else:
                    validation_warning_counter += 1
                    identification_key = (
                        f'mitochondrial_variant_{validation_warning_counter}'
                    )

                validation_output[identification_key] = (
                    variant.output_dict(test=test)
                )

            elif output_type == 'intergenic':
                validation_output['flag'] = 'intergenic'
                validation_intergenic_counter += 1

                identification_key = (
                    f'intergenic_variant_{validation_intergenic_counter}'
                )

                validation_output[identification_key] = (
                    variant.output_dict(test=test)
                )

        if with_meta:
            validation_output['metadata'] = self.add_meta()

        return validation_output

    def format_as_json(self, with_meta=True):
        return json.dumps(
            self.format_as_dict(with_meta)
        )

    def format_as_table(self, with_meta=True):
        """
        Return all validation results in table format.
        """
        outputstrings = []

        if with_meta:
            outputstrings.append(
                '# Metadata: '
                + ', '.join(
                    f'{key}: {value}'
                    for key, value in self.add_meta().items()
                )
            )

        outputstrings.append([
            'Input',
            'Warnings',
            'Select transcript',
            'HGVS_transcript',
            'HGVS_intronic_chr_context',
            'HGVS_intronic_rsg_context',
            'HGVS_RefSeqGene',
            'HGVS_LRG',
            'HGVS_LRG_transcript',
            'HGVS_Predicted_Protein',
            'HGVS_Genomic_GRCh37',
            'GRCh37_CHR',
            'GRCh37_POS',
            'GRCh37_ID',
            'GRCh37_REF',
            'GRCh37_ALT',
            'HGVS_Genomic_GRCh38',
            'GRCh38_CHR',
            'GRCh38_POS',
            'GRCh38_ID',
            'GRCh38_REF',
            'GRCh38_ALT',
            'Gene_Symbol',
            'HGNC_Gene_ID',
            'Transcript_description',
            'Alt_genomic_loci',
        ])

        empty_vcf = {
            'chr': '',
            'pos': '',
            'id': '',
            'ref': '',
            'alt': '',
        }

        for variant in self.output_list:
            if (
                    variant.output_type_flag == 'warning'
                    or variant.warnings == ['Validation error']
            ):
                self.lovd_syntax_check(variant)

            prot = ''

            if variant.hgvs_predicted_protein_consequence is not None:
                prot = (
                    variant
                    .hgvs_predicted_protein_consequence['tlr']
                )

            if variant.rna_data is not None:
                prot = variant.rna_data['translation']

            primary_loci = variant.primary_assembly_loci or {}

            grch37_data = primary_loci.get('grch37')
            grch38_data = primary_loci.get('grch38')

            grch37 = ''
            grch37_vcf = empty_vcf

            if grch37_data:
                grch37 = (
                    grch37_data['hgvs_genomic_description']
                )
                grch37_vcf = grch37_data['vcf']

            grch38 = ''
            grch38_vcf = empty_vcf

            if grch38_data:
                grch38 = (
                    grch38_data['hgvs_genomic_description']
                )
                grch38_vcf = grch38_data['vcf']

            alt_genomic = []

            if variant.alt_genomic_loci:
                for alt in variant.alt_genomic_loci:
                    for assembly in ('grch37', 'grch38'):
                        if assembly in alt:
                            alt_genomic.append(
                                alt[assembly][
                                    'hgvs_genomic_description'
                                ]
                            )

            gene_id = ''

            if variant.stable_gene_ids:
                gene_id = variant.stable_gene_ids.get(
                    'hgnc_id',
                    ''
                )

            select_tx = None

            if variant.annotations:
                select_tx = (
                    variant.annotations
                    .get('db_xref', {})
                    .get('select')
                )

            if variant.rna_data is None:
                outputstrings.append([
                    variant.original,
                    '|'.join(
                        variant.process_warnings(
                            string_all=True
                        )
                    ),
                    select_tx,
                    variant.hgvs_transcript_variant,
                    variant.genome_context_intronic_sequence,
                    variant.refseqgene_context_intronic_sequence,
                    variant.hgvs_refseqgene_variant,
                    variant.hgvs_lrg_variant,
                    variant.hgvs_lrg_transcript_variant,
                    prot,
                    grch37,
                    grch37_vcf.get('chr', ''),
                    grch37_vcf.get('pos', ''),
                    grch37_vcf.get('id', ''),
                    grch37_vcf.get('ref', ''),
                    grch37_vcf.get('alt', ''),
                    grch38,
                    grch38_vcf.get('chr', ''),
                    grch38_vcf.get('pos', ''),
                    grch38_vcf.get('id', ''),
                    grch38_vcf.get('ref', ''),
                    grch38_vcf.get('alt', ''),
                    variant.gene_symbol,
                    gene_id,
                    variant.description,
                    '|'.join(alt_genomic),
                ])

            else:
                rna_data = variant.rna_data

                outputstrings.append([
                    variant.original,
                    '|'.join(rna_data['usage_warnings']),
                    select_tx,
                    rna_data['rna_variant'],
                    '',
                    '',
                    '',
                    '',
                    '',
                    prot,
                    '',
                    '',
                    '',
                    '',
                    '',
                    '',
                    '',
                    '',
                    '',
                    '',
                    '',
                    '',
                    variant.gene_symbol,
                    gene_id,
                    variant.description,
                    '',
                ])

        return outputstrings

    def add_meta(self):
        """
        Return VariantValidator metadata.
        """
        return {
            'variantvalidator_version':
                self.validator.version,
            'variantvalidator_hgvs_version':
                self.validator.hgvsVersion,
            'vvta_version':
                self.validator.utaSchema,
            'vvseqrepo_db':
                self.validator.seqrepoVersion,
            'vvdb_version':
                self.validator.vvdbVersion,
        }

    def lovd_syntax_check(self, variant):
        """
        Add LOVD syntax-check results to the variant warnings.
        """
        if variant.lovd_syntax_check is None:
            variant.lovd_syntax_check = (
                lovd_api.lovd_syntax_check(
                    variant.original.strip(),
                    do_lovd_check=(
                        self.validator.lovd_syntax_check
                    ),
                )
            )

        lovd_result = variant.lovd_syntax_check

        if 'lovd_api_error' in lovd_result:
            return

        lovd_messages = {}
        lovd_corrections = {}

        try:
            data = lovd_result['data'][0]
        except (KeyError, IndexError, TypeError):
            data = {}

        corrected_values = data.get('corrected_values')

        if corrected_values:
            try:
                correction_items = corrected_values.items()
            except AttributeError:
                correction_items = ()

            invalid_warning_added = any(
                'LovdSyntaxcheckInvalid:' in warning
                for warning in variant.warnings
            )

            for key, val in correction_items:
                suggestion = (
                    f"LovdSyntaxcheckSuggestions: "
                    f"[suggestion = {key}, "
                    f"probability = {round(val, 2)}]"
                )

                if val == 1:
                    if key == variant.original:
                        variant.warnings.append(
                            f"LovdSyntaxcheckValid: "
                            f"{variant.original} is "
                            f"syntactically correct"
                        )

                    else:
                        if not invalid_warning_added:
                            variant.warnings.append(
                                f"LovdSyntaxcheckInvalid: "
                                f"{variant.original} is not "
                                f"syntactically correct, see "
                                f"LovdSyntaxcheckSuggestions "
                                f"for details"
                            )
                            invalid_warning_added = True

                        variant.warnings.append(suggestion)

                else:
                    if not invalid_warning_added:
                        variant.warnings.append(
                            f"LovdSyntaxcheckInvalid: "
                            f"{variant.original} is not "
                            f"syntactically correct, see "
                            f"LovdSyntaxcheckSuggestions "
                            f"for details"
                        )
                        invalid_warning_added = True

                    variant.warnings.append(suggestion)

                lovd_corrections[key] = val

        warnings = data.get('warnings')

        if warnings:
            try:
                warning_items = warnings.items()
            except AttributeError:
                warning_items = ()

            for key, val in warning_items:
                variant.warnings.append(
                    f"LovdSyntaxcheckWarning: {val}"
                )
                lovd_messages[key] = val

        errors = data.get('errors')

        if errors:
            try:
                error_items = errors.items()
            except AttributeError:
                error_items = ()

            for key, val in error_items:
                variant.warnings.append(
                    f"LovdSyntaxcheckError: {val}"
                )
                lovd_messages[key] = val

        source = lovd_result['url']
        version = lovd_result['version']

        variant.warnings.append(
            f"LovdSyntaxcheckSource: {source}"
        )
        lovd_messages['ISOURCE'] = source

        variant.warnings.append(
            f"LovdSyntaxcheckLibraryVersion: {version}"
        )
        lovd_messages['LIBRARYVERSION'] = version

        variant.lovd_messages = lovd_messages
        variant.lovd_corrections = lovd_corrections


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
