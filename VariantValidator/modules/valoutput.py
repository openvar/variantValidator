import logging
import json

logger = logging.getLogger(__name__)


class ValOutput(object):
    """This object will hold the all final, validated outputs (Variant objects) and provide methods to return this
    into a number of formats, with or without meta data"""

    def __init__(self, outputlist, validator):
        self.output_list = outputlist
        self.validator = validator

    def format_as_dict(self, with_meta=True, test=False):
        validation_output = {'flag': 'warning'}

        validation_error_counter = 0
        validation_obsolete_counter = 0
        validation_warning_counter = 0
        validation_intergenic_counter = 0

        if len(self.output_list) == 0:
            logger.warning("No variants available to output")
            validation_output['flag'] = 'empty_result'

        for variant in self.output_list:
            # For gene outputs, i.e. those that hit transcripts
            if variant.output_type_flag == 'gene':
                validation_output['flag'] = 'gene_variant'
                if variant.warnings == ['Validation error']:
                    validation_error_counter = validation_error_counter + 1
                    identification_key = 'validation_error_%s' % validation_error_counter
                else:
                    if variant.is_obsolete() and variant.hgvs_transcript_variant == '':
                        validation_obsolete_counter += 1
                        identification_key = 'obsolete_record_%s' % validation_obsolete_counter
                    else:
                        identification_key = '%s' % variant.hgvs_transcript_variant

                # if identification_key not in validation_output.keys():
                validation_output[identification_key] = variant.output_dict(test=test)

            # For warning only outputs
            # Should only ever be 1 output as an error or a warning of the following types
            # Gene symbol as reference sequence
            # Gene as transcript reference sequence

            # Note, currently there are no NM_ mito transcripts. This is expected to change. For now mito will
            # be handled here
            if variant.output_type_flag == 'warning':
                if variant.warnings == ['Validation error']:
                    validation_error_counter = validation_error_counter + 1
                    identification_key = 'validation_error_%s' % validation_error_counter
                elif variant.is_obsolete():
                    validation_obsolete_counter += 1
                    identification_key = 'obsolete_record_%s' % validation_obsolete_counter
                else:
                    validation_warning_counter = validation_warning_counter + 1
                    identification_key = 'validation_warning_%s' % validation_warning_counter
                validation_output[identification_key] = variant.output_dict(test=test)

            elif variant.output_type_flag == 'mitochondrial':
                validation_output['flag'] = 'mitochondrial'
                if variant.warnings == ['Validation error']:
                    validation_error_counter = validation_error_counter + 1
                    identification_key = 'validation_error_%s' % validation_error_counter
                elif variant.is_obsolete():
                    validation_obsolete_counter += 1
                    identification_key = 'obsolete_record_%s' % validation_obsolete_counter
                else:
                    validation_warning_counter = validation_warning_counter + 1
                    identification_key = 'mitochondrial_variant_%s' % validation_warning_counter
                validation_output[identification_key] = variant.output_dict(test=test)

            # Intergenic variants
            if variant.output_type_flag == 'intergenic':
                validation_output['flag'] = 'intergenic'
                validation_intergenic_counter = validation_intergenic_counter + 1
                identification_key = 'intergenic_variant_%s' % validation_intergenic_counter

                # Finalise the output dictionary
                validation_output[identification_key] = variant.output_dict(test=test)

        if with_meta:
            validation_output["metadata"] = self.add_meta()

        # return batch_out
        return validation_output

    def format_as_json(self, with_meta=True):
        dictionary_output = self.format_as_dict(with_meta)
        return json.dumps(dictionary_output)

    def format_as_table(self, with_meta=True):
        """
        The table format will output all results.
        :param with_meta:
        :return:
        """
        outputstrings = []
        if with_meta:
            outputstrings.append('# Metadata: ' + ', '.join(['%s: %s' % (k, v) for k, v in self.add_meta().items()]))

        outputstrings.append(['Input', 'Warnings', 'Select transcript', 'HGVS_transcript', 'HGVS_intronic_chr_context',
                              'HGVS_intronic_rsg_context', 'HGVS_RefSeqGene', 'HGVS_LRG',
                              'HGVS_LRG_transcript', 'HGVS_Predicted_Protein', 'HGVS_Genomic_GRCh37', 'GRCh37_CHR',
                              'GRCh37_POS', 'GRCh37_ID', 'GRCh37_REF', 'GRCh37_ALT', 'HGVS_Genomic_GRCh38',
                              'GRCh38_CHR', 'GRCh38_POS', 'GRCh38_ID', 'GRCh38_REF', 'GRCh38_ALT',
                              'Gene_Symbol', 'HGNC_Gene_ID', 'Transcript_description', 'Alt_genomic_loci'])
        for variant in self.output_list:
            prot = ''
            if variant.hgvs_predicted_protein_consequence is not None:
                prot = variant.hgvs_predicted_protein_consequence['tlr']
            grch37 = ''
            grch37_vcf = {'chr': '', 'pos': '', 'ref': '', 'alt': '', 'id': ''}
            if variant.primary_assembly_loci and 'grch37' in variant.primary_assembly_loci:
                grch37 = variant.primary_assembly_loci['grch37']['hgvs_genomic_description']
                grch37_vcf = variant.primary_assembly_loci['grch37']['vcf']
                grch37_vcf['id'] = '.'
            grch38 = ''
            grch38_vcf = {'chr': '', 'pos': '', 'ref': '', 'alt': '', 'id': ''}
            if variant.primary_assembly_loci and 'grch38' in variant.primary_assembly_loci:
                grch38 = variant.primary_assembly_loci['grch38']['hgvs_genomic_description']
                grch38_vcf = variant.primary_assembly_loci['grch38']['vcf']
                grch38_vcf['id'] = '.'
            alt_genomic = []
            if variant.alt_genomic_loci:
                for alt in variant.alt_genomic_loci:
                    for k, v in alt.items():
                        if k == 'grch37' or k == 'grch38':
                            alt_genomic.append(v['hgvs_genomic_description'])
            gene_id = ''
            if variant.stable_gene_ids:
                if 'hgnc_id' in variant.stable_gene_ids:
                    gene_id = variant.stable_gene_ids['hgnc_id']

            select_tx = None
            try:
                select_tx = variant.annotations['db_xref']['select']
            except TypeError:
                pass
            except KeyError:
                pass

            outputstrings.append([
                variant.original,
                '|'.join(variant.process_warnings()),
                select_tx,
                variant.hgvs_transcript_variant,
                variant.genome_context_intronic_sequence,
                variant.refseqgene_context_intronic_sequence,
                variant.hgvs_refseqgene_variant,
                variant.hgvs_lrg_variant,
                variant.hgvs_lrg_transcript_variant,
                prot,
                grch37,
                grch37_vcf['chr'],
                grch37_vcf['pos'],
                grch37_vcf['id'],
                grch37_vcf['ref'],
                grch37_vcf['alt'],
                grch38,
                grch38_vcf['chr'],
                grch38_vcf['pos'],
                grch38_vcf['id'],
                grch38_vcf['ref'],
                grch38_vcf['alt'],
                variant.gene_symbol,
                gene_id,
                variant.description,
                '|'.join(alt_genomic)
            ])
        return outputstrings

    def add_meta(self):
        """
        Returns dictionary of metadata
        :return:
        """
        metadata = {}
        metadata['variantvalidator_version'] = self.validator.version
        metadata['variantvalidator_hgvs_version'] = self.validator.hgvsVersion
        metadata['vvta_version'] = self.validator.utaSchema
        metadata['vvseqrepo_db'] = self.validator.seqrepoVersion
        metadata['vvdb_version'] = self.validator.vvdbVersion
        return metadata

# <LICENSE>
# Copyright (C) 2016-2021 VariantValidator Contributors
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
