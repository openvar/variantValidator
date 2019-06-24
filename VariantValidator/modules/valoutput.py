import os
from .logger import Logger
import json


class ValOutput(object):
    """This object will hold the all final, validated outputs (Variant objects) and provide methods to return this
    into a number of formats, with or without meta data"""

    def __init__(self, outputlist, validator):
        self.output_list = outputlist
        self.validator = validator

    def format_as_dict(self, with_meta=True):
        validation_output = {'flag': None}

        validation_error_counter = 0
        validation_obsolete_counter = 0
        validation_warning_counter = 0
        validation_intergenic_counter = 0

        if len(self.output_list) == 0:
            validation_output['flag'] = 'empty_result'

        for variant in self.output_list:
            # For gene outputs, i.e. those that hit transcripts
            if variant.output_type_flag == 'gene':
                validation_output['flag'] = 'gene_variant'
                if variant.warnings == ['Validation error']:
                    validation_error_counter = validation_error_counter + 1
                    identification_key = 'Validation_Error_%s' % validation_error_counter
                else:
                    if variant.is_obsolete() and variant.hgvs_transcript_variant == '':
                        validation_obsolete_counter += 1
                        identification_key = 'obsolete_record_%s' % validation_obsolete_counter
                    else:
                        identification_key = '%s' % variant.hgvs_transcript_variant

                # if identification_key not in validation_output.keys():
                validation_output[identification_key] = variant.output_dict()
                # else:
                # dotter = dotter + ' '
                # validation_output[identification_key + dotter] = valid_v

            # For warning only outputs
            # Should only ever be 1 output as an error or a warning of the following types
            # Gene symbol as reference sequence
            # Gene as transcript reference sequence
            if variant.output_type_flag == 'warning':
                validation_output['flag'] = 'warning'
                if variant.warnings == ['Validation error']:
                    validation_error_counter = validation_error_counter + 1
                    identification_key = 'validation_error_%s' % validation_error_counter
                elif variant.is_obsolete():
                    validation_obsolete_counter += 1
                    identification_key = 'obsolete_record_%s' % validation_obsolete_counter
                else:
                    validation_warning_counter = validation_warning_counter + 1
                    identification_key = 'validation_warning_%s' % validation_warning_counter
                validation_output[identification_key] = variant.output_dict()

            # Intergenic variants
            if variant.output_type_flag == 'intergenic':
                validation_output['flag'] = 'intergenic'
                validation_intergenic_counter = validation_intergenic_counter + 1
                identification_key = 'intergenic_variant_%s' % validation_intergenic_counter

                # Finalise the output dictionary
                validation_output[identification_key] = variant.output_dict()

        if with_meta:
            validation_output["metadata"] = self.add_meta()

        # return batch_out
        return validation_output

    def format_as_json(self, with_meta=True):
        dictionary_output = self.format_as_dict(with_meta)
        return json.dumps(dictionary_output)

    def format_as_table(self, with_meta=True):
        """
        Currently the table format will only output correctly validated results, all warnings and obsolete records will
        be squashed.
        :param with_meta:
        :return:
        """
        outputstrings = []
        if with_meta:
            outputstrings.append('#' + str(self.add_meta()))

        outputstrings.append(['Input', 'HGVS_transcript', 'HGVS_RefSeqGene', 'HGVS_LRG', 'HGVS_LRG_transcript',
                              'Gene_Symbol', 'Transcript_description'])
        for variant in self.output_list:
            if variant.output_type_flag == 'gene':
                if variant.warnings == ['Validation error'] or (variant.is_obsolete() and
                                                                variant.hgvs_transcript_variant == ''):
                    continue
                else:
                    outputstrings.append([
                        variant.original,
                        variant.hgvs_transcript_variant,
                        variant.hgvs_refseqgene_variant,
                        variant.hgvs_lrg_variant,
                        variant.hgvs_lrg_transcript_variant,
                        variant.gene_symbol,
                        variant.description
                    ])
        return outputstrings

    def add_meta(self):
        """
        Returns dictionary of metadata
        :return:
        """
        metadata = {}

        if os.environ.get("ADD_LOGS") == "True":
            logs = []
            for l in Logger.getString().split("\n"):
                logs.append(l)
            metadata["logs"] = logs

        # metadata["variant"] = batch_variant  # original input string to validate function
        # metadata["assembly"] = selected_assembly
        # metadata["transcripts"] = select_transcripts
        # metadata['seqrepo_directory'] = self.seqrepoPath
        # metadata['uta_url'] = self.utaPath
        # metadata['py_liftover_directory'] = self.liftoverPath
        # metadata['variantvalidator_data_url'] = self.db.path
        # metadata['entrez_id'] = self.entrezID
        metadata['variantvalidator_version'] = self.validator.version
        metadata['variantvalidator_hgvs_version'] = self.validator.hgvsVersion
        metadata['uta_schema'] = self.validator.utaSchema
        metadata['seqrepo_db'] = self.validator.seqrepoVersion
        return metadata
