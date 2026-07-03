"""
A set of functions for handling VRS style data, for now focusing on creating
VRS output from hgvs input.
The first function makes a complete VRS output for a given completed
variantValidator output variant, for not converting data into "aleles"
(hgvs variants) and Cis-Phased Blocks, protin is not yet handled (this probably
needs Terminus) Derivative Molecule and Adjacency are not yet handled (and would
probably be inputted and processed as hgvs for the latter, or not at all for the
former, anyway.

Each of the following are then used to produce the VRS components.
"""
import base64
import hashlib
import os
import copy
from configparser import ConfigParser
from biocommons.seqrepo import SeqRepo
from VariantValidator.modules.complex_descriptions import FEInterval
try:
    from VariantValidator.settings import CONFIG_DIR
except ImportError:
    CONFIG_DIR = None
#from vvhgvs.location import Interval

class HGVS_to_VRS:
    """
    A module for creating VRS output from HGVS input
    To avoid overcomplicating the core VV object this code is kept independent
    from VV, however it should normally only be used by VV instances. We also
    provide facilities for setting alternative ID and seq data/length sources.
    """
    def __init__(self,seq_repo = False,id_fetch = False,
                 ref_fetch_blocksize = 10, cache = False):
        """
        Module init
        set "seq_repo" to use a pre-existing seqrepo connection, for now this
            is somewhat direct and is used for the id_fetch too by default, but
            we are set up to change this in the future if needed
        set id_fetch to use a substitute other than seqrepo as a fetch location
            for the ID this could be the VVTA or VV database if needed, though
            currently they only work for transcripts
        set ref_fetch_blocksize to change the seqrepo fetch amount, given that
            the source data is block compressed the time savings from smaller
            than ten are probably not worth it even for data which nearly never
            "rolls" the input seq
        set cache to true or a dict to cache the results for this object
            instance, (or use the existing dict as a cache).  This cache is
            intended to be used for a single user input set, which will usually
            share lots of genomic & prot mappings. It does not currently limit
            the result set size or anything else. If cache is set to false
            variant_validator_output_set_to_vrs will partly overrule this to
            create a per output-set cache to limit the re-normalisation of
            already normalised genomic mappings.
        """
        if seq_repo:
            self.seq_repo = seq_repo
        elif CONFIG_DIR:
            config = ConfigParser()
            config.read(CONFIG_DIR)
            seq_repo_path = os.path.join(config["seqrepo"]["location"],
                                   config["seqrepo"]["version"])
            self.seq_repo = SeqRepo(seq_repo_path)
        else:
            raise ValueError("HGVS_to_VRS conversion eitehr needs a VV "
                             "CONFIG_DIR or a already set up seq_repo")
        if id_fetch:
            self.id_fetch = id_fetch
        else:
            self.id_fetch = self.get_vrs_id_for_seq
        self.ref_fetch_blocksize = ref_fetch_blocksize
        if isinstance(cache,dict):
            self.cache = cache
        elif cache:
            self.cache = {}
        else:
            self.cache = False

        self.sorted_vrs_digest_keys = {
                "SequenceLocation":["end", "sequenceReference","start","type"],
                "Allele":["location","state","type"],
                "LiteralSequenceExpression":['sequence', 'type'],
                "ReferenceLengthExpression":[
                    'length','repeatSubunitLength','type'],
                "LengthExpression":['length','type'],
                "SequenceReference":['refgetAccession', 'type'],
                }
        self.vrs_digest_prefixes = {
                "SequenceLocation":'ga4gh:SL.',
                "Allele":'ga4gh:VA.',
                }

    def variant_validator_output_set_to_vrs(self,variant):
        """"
        Convert a set of VV equivalent results into a VRS format, dataset.

        VV has the following fields that may contain differing HGVS data:
            hgvs_transcript_variant,
            genome_context_intronic_sequence,
            refseqgene_context_intronic_sequence,
            hgvs_refseqgene_variant,
            hgvs_lrg_transcript_variant,
            hgvs_lrg_variant,
            alt_genomic_loci,
            primary_assembly_loci,
            exonic_positions,

        rna_data gets used instead of hgvs_transcript_variant when mapping n/a
        hgvs_predicted_protein_consequence is used for valid transcripts
        To avoid un-needed churn when variant==variant and ID/VRS ID==ID/VRS ID
        we skip/reuse.
        The current output structure is:
        {
        "selected_assembly': variant.selected_assembly,
        "submitted_variant': variant.original,

        "gene_ids":{}
        "hgvs_transcript_variant":VRS (if available and not intronic)
        "hgvs_refseqgene_variant": VRS likewise (also tagged with redundant
                                            "hgvs_lrg_transcript_variant" ?)
        "hgvs_lrg_variant": VRS#
        "predicted prot consequence":VRS?
        "primary_assembly_loci":{VRS set}
        "alt_genomic_loci":[VRS list]
        "exonic_positions": variant.exonic_positions, # if wrt transcript
        "warnings":[warn list] # Set of VariantValidator warnings:
                               # 'lovd_messages' 'lovd_corrections'
                               # include bonus "VRSConversionWarning:" for:
                               # "intronic type variant descriptions not allowed
                               # in the VRS standard"
                               # also add "VRS variant data for ### between
                               # transcript locations so  gene data has not been
                               # attached"  for r type?
        }

        Add annot_level to get annotation skipped or used for partial filtering?
        """
        output = {
            'selected_assembly': variant.selected_assembly,
            'submitted_variant': variant.original,
            'warnings_and_messages': {
                'validation_warnings': variant.process_warnings(),
                'lovd_messages': variant.lovd_messages,
                'lovd_corrections': variant.lovd_corrections,
                },
            # Do we want to skip these 2 even as empty fields for intergenic
            # transcripts?
            "gene_ids":variant.stable_gene_ids,
            'transcript_description': variant.description,
        }
        if variant.gene_symbol:
            if not output["gene_ids"]:
                output["gene_ids"] = {}
            if isinstance(output["gene_ids"], dict):
                output["gene_ids"]['current_symbol'] = variant.gene_symbol
        # annotations is already JSON (eg below), and relatively redundant with
        # the other data, should add here if we want it but skip for now.
        # e.g. annotation set: '{"db_xref": {"select": false, "ncbigene":
        # "7840", "ensemblgene": null, "hgnc": "HGNC:428", "CCDS": null},
        # "chromosome": "2", "map": "2p13.1", "note":
        # "ALMS1 centrosome and basal body associated protein", "variant": "1",
        # "refseq_select": false, "mane_select": false, "ensembl_select": false,
        # "mane_plus_clinical": false}

        # abort on bad output_type_flag, or other failure flag can be one of:
        # 'gene' 'warning' 'mitochondrial' or 'intergenic'
        if output['warnings_and_messages']['validation_warnings'] == [""] or \
                variant.output_type_flag == 'warning':
            return output
        pos_warnings = \
            "VRSUncertainPositionWarning:VRS description conversion of hgvs "\
            "data with unknown locations is not fully supported, normalisation"\
            " is disabled and other data, e.g. ins vs delins data, may be lost."

        if variant.rna_data:
            assert not variant.hgvs_transcript_variant
            assert not variant.primary_assembly_loci
            # :r. RNA versions are independent of :c. versions, supposed to be
            # derived directly from RNA sequencing data, and as such skip
            # genomic and standard cDNA variant based validation steps.
            # We add a warning message about this. For now the default one, but
            # fold the output into hgvs_transcript_variant.
            # content :"usage_warnings", "rna_variant","translation" and
            # 'translation_slr'
            # tlr translation is skipped due to vrs code relying on single
            # letter alphabet
            output['vrs_transcript_variant'] = self.hgvs_single_var_to_vrs(
                    variant.rna_data['rna_variant'],variant)
            if variant.rna_data['translation_slr'].posedit:
                output['vrs_predicted_protein_consequence'] = \
                    self.hgvs_single_var_to_vrs(
                            variant.rna_data['translation_slr'])
            output['warnings_and_messages']['r_type_input'] = \
                    variant.rna_data['usage_warnings']
            if isinstance(variant.rna_data['rna_variant'], FEInterval):
                output['warnings_and_messages']['vrs_output_warnings'] = [
                        pos_warnings]
            return output

        # now check for intergenic and/or intronic
        reset_cache = False
        if self.cache is False:
            self.cache = {}
            reset_cache = True
        # use output flag instead?
        if variant.hgvs_transcript_variant and not (
                variant.hgvs_transcript_variant.type in ['c','n'] and
                self._intronic(variant.hgvs_transcript_variant)):
            output['vrs_transcript_variant'] = self.hgvs_single_var_to_vrs(
                    variant.hgvs_transcript_variant,variant)
            if variant.hgvs_refseqgene_variant:
                output['vrs_refseqgene_variant'] = self.hgvs_single_var_to_vrs(
                        variant.hgvs_refseqgene_variant,variant)
            if variant.hgvs_predicted_protein_consequence and \
                variant.hgvs_predicted_protein_consequence['prot'] and \
                variant.hgvs_predicted_protein_consequence['prot'].posedit and \
                getattr(
                    variant.hgvs_predicted_protein_consequence['prot'].posedit,
                    'edit',False):
                output['vrs_predicted_protein_consequence'] = \
                    self.hgvs_single_var_to_vrs(
                        variant.hgvs_predicted_protein_consequence['prot'])
            # LRG is redundant due to VRS checksum based seq id, but do we want
            # to add it as a alias in a comment for the output object?
            # the RSG/LRG vrs_predicted_protein_consequence is also skipped, do
            # we want to add this too?
            if isinstance(variant.hgvs_transcript_variant, FEInterval):
                output['warnings_and_messages']['vrs_output_warnings'] = [
                        pos_warnings]
        elif variant.hgvs_transcript_variant:
            # add warning for intronic data
            output['warnings_and_messages']['vrs_output_warnings'] = [ (
                "VRSIntronWarning: VRS does not handle mappings shared between"
                " genomic and transcript reference sequences. As such only the "
                "genomic mappings for this hgvs intronic transcript "
                "variant are preserved.")]
            if isinstance(variant.hgvs_transcript_variant, FEInterval):
                output['warnings_and_messages']['vrs_output_warnings'].append(
                        pos_warnings)

            if variant.hgvs_genomic:
                output['vrs_intronic_genomic_variant'] = \
                        self.hgvs_single_var_to_vrs(variant.hgvs_genomic)
            # RSG data
            if variant.hgvs_refseqgene_variant:
                output['vrs_refseqgene_variant'] = self.hgvs_single_var_to_vrs(
                        variant.hgvs_refseqgene_variant)
        elif variant.hgvs_genomic or variant.hgvs_refseqgene_variant:
            if variant.hgvs_genomic:
                output['vrs_intergenic_genomic_variant'] = \
                        self.hgvs_single_var_to_vrs(variant.hgvs_genomic)
            # RSG data
            if variant.hgvs_refseqgene_variant:
                output['vrs_refseqgene_variant'] = self.hgvs_single_var_to_vrs(
                        variant.hgvs_refseqgene_variant)
            if (
                    variant.hgvs_genomic and
                    isinstance(variant.hgvs_genomic, FEInterval)) or (
                    variant.hgvs_refseqgene_variant and
                    isinstance(variant.hgvs_refseqgene_variant, FEInterval)):
                output['warnings_and_messages']['vrs_output_warnings'] = [
                        pos_warnings]


        # continue for universally (non :r. included data)
        primary_assembly_loci = {}
        if variant.primary_assembly_loci:
            primary_assembly_loci = variant.primary_assembly_loci
            output['VRS_mappings_for_primary_assemblies'] = {}
        for gen in primary_assembly_loci:
            if gen in ['hg19', 'hg38']:
                continue
            output['VRS_mappings_for_primary_assemblies'][gen] = \
                    self.hgvs_single_var_to_vrs(
                            primary_assembly_loci[gen][
                                'hgvs_genomic_description'])
        output['VRS_mappings_for_alt_genomic_loci'] = []
        alt_genomic_loci = []

        if variant.alt_genomic_loci:
            alt_genomic_loci = variant.alt_genomic_loci
            output['VRS_mappings_for_alt_genomic_loci'] = []
        for loci in alt_genomic_loci:
            gen_key = list(loci.keys())[0]
            if gen_key in ['hg19', 'hg38']:
                continue
            output['VRS_mappings_for_alt_genomic_loci'].append(
                    self.hgvs_single_var_to_vrs(
                        loci[gen_key]['hgvs_genomic_description']
                        ))
        # Crudely re-set cache, relies on us not caching this object in odd
        # ways later, but if we do that we probably need a full caching system
        # not something this simple.
        if reset_cache:
            self.cache = False

        return output

    def _intronic(self,hgvs):
        if isinstance(hgvs.posedit.pos,FEInterval):
            if  (hgvs.posedit.pos.start.start.offset or
                 hgvs.posedit.pos.start.end.offset or
                 hgvs.posedit.pos.end.start.offset or
                 hgvs.posedit.pos.end.end.offset):
                return True
        elif (hgvs.posedit.pos.start.offset or hgvs.posedit.pos.end.offset):
            return True
        return False

    def hgvs_single_var_to_vrs(self, hgvs_variant,variant = False):
        "convert a hgvs with a single location to VRS, if possible"
        # first detect intronic variants and abort for them
        # VRS only does the genomic sequnce in this case
        if not self.cache is False:
            hgvs_key = str(hgvs_variant)
            if hgvs_key in self.cache:
                return self.cache[hgvs_key]
        if hgvs_variant.type in ['c','n','t'] and self._intronic(hgvs_variant):
            return None
        if hgvs_variant.type == 'c':
            # fix transcript coordinates to n
            if not variant:
                raise ValueError(
                    "Meta-variant with a variant mapper needed for c->n case")
            if hgvs_variant.posedit.pos.uncertain:
                ref_seq = hgvs_variant.posedit.edit.ref
                # handle paired ambig ends
                if isinstance(hgvs_variant.posedit.pos, FEInterval):
                    hgvs1 = copy.deepcopy(hgvs_variant)
                    hgvs1.posedit.pos = hgvs1.posedit.pos.start
                    hgvs1 = variant.lose_vm.c_to_n(hgvs1)
                    hgvs2 = copy.deepcopy(hgvs_variant)
                    hgvs2.posedit.pos = hgvs2.posedit.pos.end
                    hgvs2 = variant.lose_vm.c_to_n(hgvs2)
                    hgvs_variant.type = 'n'
                    hgvs_variant.posedit.pos.start = hgvs1.posedit.pos
                    hgvs_variant.posedit.pos.end = hgvs2.posedit.pos
                else:
                    hgvs_variant = variant.lose_vm.c_to_n(hgvs_variant)
                # prevent ref for unc pos from being reset (eg NNN>TTTGGGTTTGGG)
                hgvs_variant.posedit.edit.ref = ref_seq
            else:
                hgvs_variant = variant.lose_vm.c_to_n(hgvs_variant)

        # From this point genomic and transcript VRS should handle the same
        # For objects IDs using child objects as variables the "canonical"
        # set of content is used for the ID

        # fetch vrs id for sequence
        # Derive "VRS refrence" using ID e.g.
        # {
        #  "type": "SequenceReference",
        #  "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
        #  "label": "NC_000007.14"
        #}
        #GA4GH Digest Prefix:None Inherent:[‘refgetAccession’, ‘type’]
        vrs_ref_seq = {
                "type": "SequenceReference",
                "refgetAccession":self.id_fetch(hgvs_variant.ac),
                "label": str(hgvs_variant.ac)}
        # Handle sequences with start-stop ranges, normalisation skipped
        if hgvs_variant.posedit.pos.uncertain:
            # could check for FEInterval but since this sets uncertain don't
            # for now. Currently we do not attempt to handle VRS types more
            # complex than LiteralSequenceExpression.
            # Ranges are specified but in VRS 2.1, but not
            # normalisation for range input so we don't for now.
            if not isinstance(hgvs_variant.posedit.pos, FEInterval):
                vrs_start = [
                    int(hgvs_variant.posedit.pos.start.base)-1,
                    int(hgvs_variant.posedit.pos.end.base)]
                vrs_end = [
                    int(hgvs_variant.posedit.pos.start.base)-1,
                    int(hgvs_variant.posedit.pos.end.base)]
            else:
                vrs_start = [
                    int(hgvs_variant.posedit.pos.start.start.base)-1,
                    int(hgvs_variant.posedit.pos.start.end.base)]
                vrs_end = [
                    int(hgvs_variant.posedit.pos.end.start.base)-1,
                    int(hgvs_variant.posedit.pos.end.end.base)]
            vrs_loc = {
                "type": "SequenceLocation",
                "sequenceReference":vrs_ref_seq,
                "start":vrs_start,
                "end":vrs_end,
                }
            # if del set length expr, if dupe likewise, else handle as
            # LiteralSequenceExpression
            if hgvs_variant.posedit.edit.type in 'dup':
                min_len = 2* max(vrs_end[0] - vrs_start[1],0)
                max_len = 2*(vrs_end[1] - vrs_start[0])
                vrs_state = {
                    "type": "LengthExpression",
                    "length": [min_len,max_len]}
            elif hgvs_variant.posedit.edit.alt:
                # ins, sub, inv, or delins
                vrs_state = {
                    "type": "LiteralSequenceExpression",
                    "sequence": hgvs_variant.posedit.edit.alt,
                    }
            else:
                # del bases represent the deleted ref length
                # so length is length - ref length
                min_len = max(
                        vrs_end[0] - vrs_start[1] - \
                            len(hgvs_variant.posedit.edit.ref),
                        0)
                max_len = vrs_end[1] - vrs_start[0] - \
                    len(hgvs_variant.posedit.edit.ref)
                vrs_state = {
                    "type": "LengthExpression",
                    "length":[min_len,max_len]}
            vrs_allele = {
                    "type": "Allele",
                    "location":vrs_loc,
                    "state":vrs_state,}
            vrs_allele['id'] = self.vrs_flatten(vrs_allele,digest=True)
            if not self.cache is False:
                self.cache[hgvs_key] = vrs_allele
            return vrs_allele
        # Derive normalised start and stop
        # Derive VRS sequence location from VRS sequence and normalised location
        # start by converting hgvs (1 based ins special) to vrs (0 based)
        # coordinates, also handle empty refs, this probably should not be an
        # issue for us in normal use but still needs to be handled for the
        # general case.
        ref_seq = ''
        alt_seq = ''
        hgvs_type = hgvs_variant.posedit.edit.type
        if hgvs_type == 'ins':
            #may need to switch to {"type": "Number", "value": num} post VRS 2.1
            vrs_start = int(hgvs_variant.posedit.pos.start.base)
            vrs_end = int(hgvs_variant.posedit.pos.end.base) -1
            ref_seq = ''
            alt_seq = hgvs_variant.posedit.edit.alt
        else:
            vrs_start = int(hgvs_variant.posedit.pos.start.base) -1
            vrs_end = int(hgvs_variant.posedit.pos.end.base)
            ref_seq = hgvs_variant.posedit.edit.ref
            if not ref_seq:
                ref_seq = self.seq_repo[hgvs_variant.ac][vrs_start:vrs_end]
            if hgvs_type == 'dup':
                alt_seq = ref_seq + ref_seq
            elif hgvs_type == 'identity':
                alt_seq = ref_seq
            else:
                alt_seq = hgvs_variant.posedit.edit.alt

        # To finish the creation of a VRS allele we now need to derive the
        # state, which could be either a ReferenceLengthExpression or a
        # LiteralSequenceExpression

        # If both the ref and alt disappears on down-normalisation (or variant
        # type is otherwise is equal) ReferenceLengthExpression is used with
        # repeatSubunitLength set to the total length

        # cheat shortcut for = variants
        if ref_seq == alt_seq:
            #{ "type": "Number", "value": 55}
            vrs_loc = {
                    "type": "SequenceLocation",
                    "sequenceReference":vrs_ref_seq,
                    "start":vrs_start,
                    "end":vrs_end,
                    }
            vrs_loc["id"] = self.vrs_flatten(vrs_loc,digest=True)
            vrs_state = {
                    "type": "ReferenceLengthExpression",
                    "length": vrs_end - vrs_start,
                    "repeatSubunitLength": vrs_end - vrs_start
                    }
            vrs_allele = {
                    "type": "Allele",
                    "location":vrs_loc,
                    "state":vrs_state,}
            vrs_allele["id"] = self.vrs_flatten(vrs_allele,digest=True)
            if not self.cache is False:
                self.cache[hgvs_key] = vrs_allele
            return vrs_allele
        # need to pre zero base normalise start/end
        full_ref, full_alt,min_ref,min_alt,vrs_start,vrs_end,derived = \
                self.normalise_outwards(
                        hgvs_variant.ac,
                        ref_seq,alt_seq,
                        vrs_start,
                        vrs_end)
        vrs_loc = {
                "type": "SequenceLocation",
                "sequenceReference":vrs_ref_seq,
                "start":vrs_start,
                "end":vrs_end,
                }
        vrs_loc["id"] =  self.vrs_flatten(vrs_loc,digest=True)

        #Do we want to "assert min_alt or min_ref"? This should only happen on
        # == and be disallowed by the above code, so skip for now.
        # If both ref and alt are full post down-normalisation
        # LiteralSequenceExpression is the chosen type.
        # if we have a indel even minimalised then this must be a
        # LiteralSequenceExpression
        # If we have a pure (un-ambiguous location) ins post up-normalisation
        # LiteralSequenceExpression is the chosen output type,
        # for the other extreme ReferenceLengthExpressions are used for del's
        if full_alt and not full_ref:#ins not norm outwards
            vrs_state = {
                    "type": "LiteralSequenceExpression",
                    "sequence": full_alt,
                    }
        elif full_ref and not full_alt: # del not norm outwards
            vrs_state = {
                    "type": "ReferenceLengthExpression",
                    "length": 0,
                    "repeatSubunitLength": vrs_end-vrs_start
                    }
        # For other del or ins after down-normalisation, that expand to delins
        # after up-normalisation, if the ref or alt can fully be recreated from
        # the flank we use ReferenceLengthExpressions
        elif derived:#full_ref and full_alt is implicit beyond this point
            vrs_state = {
                    "type": "ReferenceLengthExpression",
                    "length": len(full_alt),
                    # one of these must be 0 length, so + works as well as an if
                    "repeatSubunitLength":len(min_ref+min_alt)
                    }
        # otherwise we use a LiteralSequenceExpression
        elif not derived:
            vrs_state = {
                    "type": "LiteralSequenceExpression",
                    "sequence": full_alt,
                    }
        vrs_allele = {
                "type": "Allele",
                "location":vrs_loc,
                "state":vrs_state,}
        vrs_allele['id'] = self.vrs_flatten(vrs_allele,digest=True)
        if not self.cache is False:
            self.cache[hgvs_key] = vrs_allele
        return vrs_allele

    def vrs_flatten(self, vrs_dict, dict_keys = None, digest = False):
        """
        Derive a canonical VRS object or it's digest ID from a dictionaried set
        of VRS type inputs.
        VRS demands not the contents of the fields concatenated somehow, but
        instead a RFC 8785 Canonicalised JSON object containing just said fields
        (with IDs), because we can't have simple things.
        Keys must be sorted alphabetically, but we trust that any given
        dict_keys are pre-sorted.
        We also presume that the object type is accurate, and all numbers are
        integers, this fits our current use case.
        The official VRS code uses 'canoicaljson' as a dependency, but given the
        limits on the allowed input, this is overkill.
        """
        if dict_keys:
            # force specific dict keys, used in development,
            # may be removed later
            keys = dict_keys
        else:
            keys = self.sorted_vrs_digest_keys[vrs_dict['type']]
        json = '{'
        for key in keys:
            if isinstance(vrs_dict[key],int):
                json = json + f'"{key}":{vrs_dict[key]},'
            elif isinstance(vrs_dict[key],str):
                json = json + f'"{key}":"{vrs_dict[key]}",'
            elif isinstance(vrs_dict[key],list):
                assert len(vrs_dict[key]) == 2
                # 1 must be int, other may be None, start<stop could also check?
                # for now as all input is internally sourced we don't
                start, stop = vrs_dict[key]
                if start is None:
                    start = 'null'
                if stop is None:
                    stop = 'null'
                json = json + f'"{key}":"[{start},{stop}]",'
            elif isinstance(vrs_dict[key],dict):
                if 'type' in vrs_dict[key] and \
                        vrs_dict[key]['type'] in self.vrs_digest_prefixes:
                    if 'id' not in vrs_dict[key]:
                        vrs_dict[key]['id'] = self.vrs_flatten(
                                vrs_dict[key], digest = True)
                    json = json + f'"{key}":"{vrs_dict[key]["id"][9:]}",'
                else:
                    json = json + f'"{key}":{self.vrs_flatten(vrs_dict[key])},'
            else:
                raise TypeError("Unexpected type of input for VRS conversion")
        json = json[:-1] + '}'
        if not digest:
            return json
        hash_val = base64.urlsafe_b64encode(
                hashlib.sha512(
                    json.encode('utf-8')
                    ).digest()[:24]
                ).decode("utf-8")
        if vrs_dict['type'] in self.vrs_digest_prefixes:
            return self.vrs_digest_prefixes[vrs_dict['type']] + hash_val
        raise ValueError("VRS prefix not set for VRS ids of current type")

    def get_vrs_id_for_seq(self,hgvs_seq_id):
        """Convert a hgvs id into a  ga4gh truncated hash as used by RefGet,
        which is the standard used by VRS.
        We should be able to calculate this ourselves, however instead this code
        uses (and needs) access to a full seqrepo instance instead, which will
        contain a pre-generated valid id, returned by the translate_alias
        function.
        # for python the code for the seq -> id looks like
        def ga4gh_digest(seq):
            digest = hashlib.sha512(seq.encode('utf-8')).digest()
            url_safe_di = base64.urlsafe_b64encode(digest[:24]).decode("utf-8")
            return f"ga4gh:SQ.{url_safe_di}"
        """
        # we could also use VVTA for non genomic if it is faster, but this may
        # need a different translation for the url safe encoding
        out = self.seq_repo.translate_alias(hgvs_seq_id)
        for seq_id in out:
            if seq_id.startswith('ga4gh:SQ.'):
                return seq_id[6:]
        return False
        #Failure to find ID for sequence should trigger a KeyError in SeqRepo

    def _trim(self, ref_seq, alt_seq):
        """
        Trim a VRS compliant ref and alt seq to match the VRS standard
        Trim suffix (trim right) first then prefix (left), this is specified as
        the correct trim order in the VRS standard
        """
        # Early return if we can't trim (only ins or del)
        if not ref_seq or not alt_seq:
            return ref_seq, alt_seq, 0, 0
        # Early return for ==, is not classically normalised
        if alt_seq == ref_seq:
            # return empty to flag this as ==
            return '', '', 0, 0

        loop_len = len(ref_seq)
        # Early return shortcut for most common single base transition
        if loop_len == 1 and len(alt_seq) == 1:
            return ref_seq, alt_seq, 0, 0

        # Finally start to trim possible, hgvs delins type
        #Trim Right
        loop_len = min(loop_len,len(alt_seq))
        alt_len = 0
        for i in range(-1,-loop_len-1,-1):
            if ref_seq[i] == alt_seq[i]:
                alt_len = i
            else:
                break
        if alt_len:
            ref_seq = ref_seq[:alt_len]
            alt_seq = alt_seq[:alt_len]
        sufix = alt_len

        # Early return if fully trimmed i.e. one seq is a subset of the other
        loop_len = min(len(ref_seq),len(alt_seq))
        if not loop_len:
            return ref_seq, alt_seq, 0, sufix

        #Trim Left
        alt_len = 0
        for i in  range(loop_len):
            if ref_seq[i] == alt_seq[i]:
                alt_len = i + 1
            else:
                break
        if alt_len is not None:
            ref_seq = ref_seq[alt_len:]
            alt_seq = alt_seq[alt_len:]

        return ref_seq, alt_seq, alt_len, sufix

    def _push_r(self,target_ac,start,ref_flank,roll_seq):
        """
        Function for "pushing" a minimised differing sequence along a reference
        this function acts by "rolling" the diff sequence right. Used to VRS
        normalise sequence locations outwards to the full possible location
        span.
        """
        edge_found = False
        curr_roll_pos = 0
        ref_chunk_start = 0
        roll_len = len(roll_seq)
        fully_derived = False
        flank_section = ''
        nomalised_edge = start
        curr_flank_idx = 0
        while ref_flank:
            for curr_flank_idx, flank_char in enumerate(ref_flank):
                if flank_char != roll_seq[curr_roll_pos]:
                    edge_found = True
                    break
                curr_roll_pos = curr_roll_pos + 1
                if curr_roll_pos >= roll_len:
                    fully_derived = True
                    curr_roll_pos = 0
            if edge_found:
                nomalised_edge =  start + ref_chunk_start + curr_flank_idx
                flank_section = flank_section + \
                        ref_flank[:curr_flank_idx]
                break
            ref_chunk_start = ref_chunk_start + len(ref_flank)
            flank_section = flank_section + ref_flank
            ref_flank = self.seq_repo[target_ac][
                    start + ref_chunk_start:
                    start + ref_chunk_start + self.ref_fetch_blocksize]
        if fully_derived:
            bases_non_derived = 0
        else:
            bases_non_derived = roll_len-curr_roll_pos
        if not edge_found:
            if not ref_flank:
                curr_flank_idx = 0
            nomalised_edge = start + ref_chunk_start + curr_flank_idx

        return flank_section,nomalised_edge,bases_non_derived

    def _push_l(self,target_ac,start,ref_flank,roll_seq):
        """
        Function for "pushing" a minimised differing sequence along a reference
        this function acts by "rolling" the diff sequence left. Used to VRS
        normalise sequence locations outwards to the full possible location
        span.
        """
        edge_found = False
        curr_roll_pos = 0
        ref_chunk_start = 0
        roll_len = len(roll_seq)
        fully_derived = False
        curr_flank_idx = 0
        roll_seq = roll_seq[::-1]
        flank_section = ''
        nomalised_edge = 0 # if we stop before an edge is found then we hit 0
        while ref_flank:
            curr_flank_idx = 0
            for curr_flank_idx, flank_char in enumerate(reversed(ref_flank)):
                if flank_char != roll_seq[curr_roll_pos]:
                    edge_found = True
                    break
                curr_roll_pos = curr_roll_pos + 1
                if curr_roll_pos >= roll_len:
                    fully_derived = True
                    curr_roll_pos = 0
            if edge_found:
                nomalised_edge = start - ref_chunk_start - curr_flank_idx
                if curr_flank_idx :
                    flank_section = ref_flank[-curr_flank_idx:] + flank_section
                break
            ref_chunk_start = ref_chunk_start + len(ref_flank)
            flank_section = ref_flank + flank_section
            # existing end was truncated due to sequence length
            if start - ref_chunk_start > 0:
                ref_flank = self.seq_repo[target_ac][
                        max(start-ref_chunk_start-self.ref_fetch_blocksize,0):
                        start - ref_chunk_start]
            else:
                ref_flank = ''
        if fully_derived:
            bases_non_derived = 0
        else:
            bases_non_derived = roll_len-curr_roll_pos
        return flank_section,nomalised_edge,bases_non_derived

    def normalise_outwards(
            self,
            target_ac,
            ref_seq,alt_seq,
            start,end):
        """
        Normalise a variant, defined in 0 based coordinates, outwards to the
        limits of its ambiguity e.g. a ins GC in GCGC becomes ref GCGC alt
        GCGCGC.

        This is not 'gap/mapping normalisation' like what VV does inside the
        hard left/right hgvs2vcf functions, and as such ignores things like
        intron-exon boundaries, which are not mentioned in the VRS standard for
        normalisation, thus if we want to handle these then they need to be
        handled external to this normalisation step.

        This should otherwise be equivalent to fully pushing right and left, for
        ambiguous locations, at the same time.
        """
        if not ref_seq:
            assert start == end
        else:
            assert len(ref_seq) == end - start, ref_seq+str(end)+' '+str(start)
        if not alt_seq:
            alt_seq = '' # for del alt seq may start as None
        # We first minimise then re-normalise back to full, minimalisation is
        # is needed to allow "rolling" for variants like del C ins CGC in a GC
        # stretch. We may be able to skip this step with hgvs pre-normalise d
        # input.
        orig_ref = ref_seq
        orig_alt = alt_seq
        ref_seq, alt_seq, prefix, suffix  = self._trim(ref_seq, alt_seq)
        # VRS does not do any normalisation on == (which trim to '')
        if not ref_seq and not alt_seq:
            return orig_ref,orig_alt,orig_ref,orig_alt,start,end,0
        # suffix is negative so handle a such
        start = start + prefix
        end = end + suffix
        # VRS also does not push/expand normalise indels (or base transitions)
        if ref_seq and alt_seq:
            return ref_seq, alt_seq, ref_seq, alt_seq,\
                    start, end,0
        # we must now either have an ins or a del but not a delins

        # First grab seq
        ref_left_flank = ''
        ref_right_flank = ''
        if len(ref_seq) < self.ref_fetch_blocksize*2:
            # Even if we go past the end we only get an error if we exceed the
            # limits of seqrepo's index, otherwise we silently truncate!
            seq = self.seq_repo[target_ac][
                        max(start - self.ref_fetch_blocksize,0):
                        end + self.ref_fetch_blocksize]
            if start - self.ref_fetch_blocksize >= 0:
                ref_left_flank = seq[:self.ref_fetch_blocksize]
                ref_right_flank = seq[self.ref_fetch_blocksize+len(ref_seq):]
            else:
                # this is simpler since we can run relative to the start
                ref_left_flank = seq[:start]
                ref_right_flank = seq[end:]
        else:
            ref_left_flank = self.seq_repo[target_ac][
                    max(start - self.ref_fetch_blocksize,0):start]
            ref_right_flank = self.seq_repo[target_ac][
                    end:end + self.ref_fetch_blocksize]

        # Roll in both directions, this uses a roll_length and roll_seq that is
        # from either the seq and length of the ins, or the seq and length of
        # the del.
        roll_seq = ref_seq
        if alt_seq:
            roll_seq = alt_seq
        ref_left_flank, start, bases_non_derived_l = self._push_l(
                target_ac,
                start,
                ref_left_flank,
                roll_seq)
        ref_right_flank, end, bases_non_derived_r = self._push_r(
                target_ac,
                end,
                ref_right_flank,
                roll_seq)
        fully_derived = True
        if bases_non_derived_r + bases_non_derived_l > len(roll_seq):
            fully_derived = False

        #rf = ref_left_flank + ref_seq + ref_right_flank
        #lf = ref_left_flank + alt_seq + ref_right_flank
        return ref_left_flank + ref_seq + ref_right_flank,\
                ref_left_flank + alt_seq + ref_right_flank,\
                ref_seq, alt_seq,\
                start, end, fully_derived

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
