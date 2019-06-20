import os
from .liftover import liftover
from .logger import Logger


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

                # Attempt to liftover between genome builds
                # Note: pyliftover uses the UCSC liftOver tool.
                # https://pypi.org/project/pyliftover/
                genomic_position_info = variant.primary_assembly_loci
                for g_p_key in list(genomic_position_info.keys()):
                    build_to = ''
                    build_from = ''

                    # Identify the current build and hgvs_genomic descripsion
                    if 'hg' in g_p_key:
                        # incoming_vcf = genomic_position_info[g_p_key]['vcf']
                        # set builds
                        if g_p_key == 'hg38':
                            build_to = 'hg19'
                            build_from = 'hg38'
                        if g_p_key == 'hg19':
                            build_to = 'hg38'
                            build_from = 'hg19'
                    elif 'grc' in g_p_key:
                        # incoming_vcf = genomic_position_info[g_p_key]['vcf']
                        # set builds
                        if g_p_key == 'grch38':
                            build_to = 'GRCh37'
                            build_from = 'GRCh38'
                        if g_p_key == 'grch37':
                            build_to = 'GRCh38'
                            build_from = 'GRCh37'

                    # Liftover
                    lifted_response = liftover(genomic_position_info[g_p_key]['hgvs_genomic_description'], build_from,
                                               build_to, variant.hn, variant.reverse_normalizer,
                                               variant.evm, self.validator)

                    # Sort the respomse into primary assembly and ALT
                    primary_assembly_loci = {}
                    alt_genomic_loci = []
                    for build_key, accession_dict in list(lifted_response.items()):
                        try:
                            accession_key = list(accession_dict.keys())[0]
                            if 'NC_' in accession_dict[accession_key]['hgvs_genomic_description']:
                                primary_assembly_loci[build_key.lower()] = accession_dict[accession_key]
                            else:
                                alt_genomic_loci.append({build_key.lower(): accession_dict[accession_key]})

                        # KeyError if the dicts are empty
                        except KeyError:
                            continue
                        except IndexError:
                            continue

                    # Add the dictionaries from lifted response to the output
                    if primary_assembly_loci != {}:
                        variant.primary_assembly_loci = primary_assembly_loci
                    if alt_genomic_loci:
                        variant.alt_genomic_loci = alt_genomic_loci

                # Finalise the output dictionary
                validation_output[identification_key] = variant.output_dict()

        if with_meta:
            validation_output["metadata"] = self.add_meta()

        # return batch_out
        return validation_output

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
