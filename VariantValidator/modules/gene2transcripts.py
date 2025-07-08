import re
import copy
import json
import vvhgvs.exceptions
import vvhgvs.sequencevariant
import logging
from . import utils as fn
from . import seq_data, lovd_api


# Pre compile variables
vvhgvs.global_config.formatting.max_ref_length = 1000000

logger = logging.getLogger(__name__)


def gene2transcripts(g2t, query, validator=False, bypass_web_searches=False, select_transcripts=None,
                     transcript_set=None, genome_build=None, bypass_genomic_spans=False, lovd_syntax_check=False):
    """
    Generates a list of transcript (UTA supported) and transcript names from a gene symbol or RefSeq transcript ID
    :param g2t: variant object
    :param query: string gene symbol or RefSeq ID (e.g. NANOG or NM_024865.3) or if used internally, variant object
    :param validator: Validator object
    :param bypass_web_searches: bool  Shortens the output by looping out code not needed for internal processing
    :param select_transcripts: bool False or string of transcript IDs "|" delimited
    :param transcript_set: String that defines all, refseq or ensembl
    :param genome_build: String GRCh37 or GRCh38
    :param bypass_genomic_spans: bool
    :return: dictionary of transcript information
    """
    # Set LOVD data
    lovd_corrections = None
    lovd_messages = None

    # submitted
    submitted = copy.copy(query)

    # List of transcripts
    sel_tx_lst = False
    if select_transcripts is not None:
        try:
            sel_tx_lst = json.loads(select_transcripts)
        except json.decoder.JSONDecodeError:
            sel_tx_lst = [select_transcripts]

    if bypass_web_searches is True:
        pass
    else:
        # Remove whitespace
        query = ''.join(query.split())

        try:
            query = int(query)
        except ValueError:
            pass
        if isinstance(query, int):
            query = f"HGNC:{str(query)}"

        # Search by gene IDs
        if "HGNC:" in query:
            store_query = query
            query = query.upper()
            if store_query != query:
                if lovd_syntax_check is True:  # Try LOVD syntax checker
                    lovd_messages, lovd_corrections = lovd_syntax_check_g2t(submitted, lovd_syntax_check)
            query = g2t.db.get_stable_gene_id_from_hgnc_id(query)[1]
            if query == "No data":
                try:
                    query = g2t.db.get_transcripts_from_annotations(store_query)
                    if "none" not in query[0]:
                        for tx in query:
                            if tx[5] != "unassigned":
                                query = tx[5]
                                break
                    elif lovd_syntax_check is True: # Try LOVD syntax checker
                        error = f"Unable to recognise {submitted}. Please provide a gene symbol"
                        lovd_messages, lovd_corrections = lovd_syntax_check_g2t(submitted, lovd_syntax_check)
                        if lovd_messages:
                            return {
                                'error': error,
                                "requested_symbol": submitted,
                                "lovd_messages": lovd_messages,
                                "lovd_corrections": lovd_corrections}
                        else:
                            return {'error': error, "requested_symbol": submitted }
                    else:
                        return {'error': f"Unable to recognise {submitted}. Please provide a gene symbol",
                                "requested_symbol": submitted}
                except TypeError:
                    pass

        query = query.upper()
        if re.search(r'\d+ORF\d+', query):
            query = query.replace('ORF', 'orf')

        if query != submitted:
            if lovd_syntax_check is True:  # Try LOVD syntax checker
                lovd_messages, lovd_corrections = lovd_syntax_check_g2t(submitted, lovd_syntax_check)

        # Quick check for LRG
        elif 'LRG' in query:
            lrg_id = query.split('T')[0]
            lrg_to_hgnc = g2t.db.get_lrg_data_from_lrg_id(lrg_id)
            if lrg_to_hgnc and lrg_to_hgnc[0] != 'none':
                query = lrg_to_hgnc[2]

        # Quick check for blank form
        if query == '':
            return {'error': 'Please enter HGNC gene name or transcript identifier (NM_, NR_, or ENST)'
                    , "requested_symbol": submitted}

    # Gather transcript information lists
    if bypass_web_searches is True:
        tx_for_gene = []
        tx_info = g2t.hdp.get_tx_identity_info(query.hgvs_coding.ac)

        # Add primary assembly queries
        for builds in query.primary_assembly_loci.keys():
            if "grc" in builds:
                tx_for_gene.append([query.gene_symbol,
                                    tx_info[3],
                                    0,
                                    query.hgvs_coding.ac,
                                    query.primary_assembly_loci[builds]['hgvs_genomic_description'].split(":")[0],
                                    validator.alt_aln_method])

        # Add refseqgene if available
        if query.hgvs_refseqgene_variant and 'NG_' in query.hgvs_refseqgene_variant.ac:
            tx_for_gene.append([query.gene_symbol,
                                tx_info[3],
                                0,
                                query.hgvs_coding.ac,
                                query.hgvs_refseqgene_variant.ac,
                                validator.alt_aln_method])

    else:
        # Search for gene symbol on Transcript inputs
        hgnc = query
        if 'NM_' in hgnc or 'NR_' in hgnc or "ENST" in hgnc:

            # Remove version
            if '.' in hgnc:
                hgnc = hgnc.split('.')[0]

            # Find latest version in UTA
            found_res = False
            for version in range(25):
                refresh_hgnc = hgnc + '.' + str(version)
                try:
                    g2t.hdp.get_tx_identity_info(refresh_hgnc)
                    tx_found = refresh_hgnc
                    found_res = True
                except vvhgvs.exceptions.HGVSError as e:
                    logger.debug("Except passed, %s", e)
            if not found_res:
                return {'error': 'No transcript definition for (tx_ac=' + hgnc + ')',
                        "requested_symbol": submitted}

            try:
                tx_info = g2t.hdp.get_tx_identity_info(tx_found)
            except vvhgvs.exceptions.HGVSError as e:
                return {'error': str(e),
                        "requested_symbol": submitted}

            hgnc = tx_info[6]
            hgnc2 = g2t.db.get_hgnc_symbol(hgnc)

            if re.match("LOC", hgnc2) and not re.match("LOC", hgnc):
                hgnc = hgnc
            else:
                hgnc = hgnc2

        # First perform a search against the input gene symbol or the symbol inferred from UTA
        symbol_identified = False
        vvta_record = g2t.hdp.get_gene_info(hgnc)

        # Check for a record
        if vvta_record is not None:
            current_sym = hgnc
            gene_name = vvta_record[3]
            hgnc_id = vvta_record[0]
            previous_sym = vvta_record[5]
            symbol_identified = True

        # No record found, is it a previous symbol?
        else:
            # Look up current name
            vvta_record = g2t.hdp.get_gene_info_by_alias(hgnc)
            if vvta_record is not None:
                if len(vvta_record) == 1:
                    current_sym = vvta_record[0][1]
                    gene_name = vvta_record[0][3]
                    hgnc_id = vvta_record[0][0]
                    previous_sym = hgnc
                    symbol_identified = True
                elif len(vvta_record) > 1:
                    return {'error': '%s is a previous symbol for %s genes. '
                                     'Refer to https://www.genenames.org/' % (hgnc, str(len(vvta_record))),
                            "requested_symbol": submitted}
                else:
                    # Is it an updated symbol?
                    try:
                        old_symbol = g2t.db.get_uta_symbol(hgnc)
                    except Exception:
                        if lovd_syntax_check is True:
                            lovd_messages, lovd_corrections = lovd_syntax_check_g2t(submitted, lovd_syntax_check)
                            if len(lovd_corrections) == 1 and list(lovd_corrections.values())[0] == 1:
                                corrected_symbol = list(lovd_corrections.keys())[0]
                                vvta_record = g2t.hdp.get_gene_info(corrected_symbol)
                                if vvta_record is not None:
                                    current_sym = corrected_symbol
                                    gene_name = vvta_record[3]
                                    hgnc_id = vvta_record[0]
                                    previous_sym = ""
                                    symbol_identified = True
                                else:
                                    return {'error': 'Unable to recognise gene symbol %s' % hgnc,
                                            "requested_symbol": submitted,
                                            "lovd_messages": lovd_messages,
                                            "lovd_corrections": lovd_corrections}
                            else:
                                return {'error': 'Invalid characters identified. Please provide a valid submission',
                                    "requested_symbol": submitted,
                                    "lovd_messages": lovd_messages,
                                    "lovd_corrections": lovd_corrections}
                        else:
                            return {'error': 'Invalid characters identified. Please provide a valid submission',
                                    "requested_symbol": submitted}
                    if old_symbol is not None:
                        vvta_record = g2t.hdp.get_gene_info(old_symbol)
                        if vvta_record is not None:
                            current_sym = hgnc
                            gene_name = vvta_record[3]
                            hgnc_id = vvta_record[0]
                            previous_sym = old_symbol
                            symbol_identified = True

        if symbol_identified is False:
            if lovd_syntax_check is True:
                lovd_messages, lovd_corrections = lovd_syntax_check_g2t(submitted, lovd_syntax_check)
                if len(lovd_corrections) == 1 and list(lovd_corrections.values())[0] == 1:
                    corrected_symbol = list(lovd_corrections.keys())[0]
                    vvta_record = g2t.hdp.get_gene_info(corrected_symbol)
                    if vvta_record is not None:
                        current_sym = corrected_symbol
                        gene_name = vvta_record[3]
                        hgnc_id = vvta_record[0]
                        previous_sym = ""
                    else:
                        return {'error': 'Unable to recognise gene symbol %s' % hgnc,
                                "requested_symbol": submitted,
                                "lovd_messages": lovd_messages,
                                "lovd_corrections": lovd_corrections}
                else:
                    return {'error': 'Unable to recognise gene symbol %s' % hgnc,
                            "requested_symbol": submitted,
                            "lovd_messages": lovd_messages,
                            "lovd_corrections": lovd_corrections}
            else:
                return {'error': 'Unable to recognise gene symbol %s' % hgnc,
                        "requested_symbol": submitted}

        # Get transcripts
        tx_for_gene = g2t.hdp.get_tx_for_gene(current_sym)
        if len(tx_for_gene) == 0:
            tx_for_gene = g2t.hdp.get_tx_for_gene(previous_sym)
        if len(tx_for_gene) == 0:
            return {'error': 'Unable to retrieve data from the VVTA, please contact admin',
                    "requested_symbol": submitted}

    # Loop through each transcript and get the relevant transcript description
    genes_and_tx = []
    recovered = []

    # Remove un-selected transcripts
    if sel_tx_lst is not False:
        kept_tx = []
        for tx in tx_for_gene:
            annotation = g2t.db.get_transcript_annotation(tx[3])
            if tx[3] in sel_tx_lst:
                kept_tx.append(tx)
            elif "mane_select" in sel_tx_lst:
                if '"mane_select": true' in annotation:
                    kept_tx.append(tx)
            elif "mane" in sel_tx_lst:
                if '"mane_select": true' in annotation or '"mane_plus_clinical": true' in annotation:
                    kept_tx.append(tx)
            elif "select" in sel_tx_lst:
                if '"mane_select": true' in annotation or '"refseq_select": true' in annotation \
                        or '"ensembl_select": true' in annotation:
                    kept_tx.append(tx)
            elif "all" in sel_tx_lst or None in sel_tx_lst:
                kept_tx.append(tx)

        tx_for_gene = kept_tx

    for line in tx_for_gene:
        # Remove unrequested transcript_sets
        if transcript_set == "refseq":
            if "ENST" in line[3]:
                continue
        elif transcript_set == "ensembl":
            if re.match("N[MR]_", line[3]):
                continue

        if (line[3].startswith('NM_') or line[3].startswith('NR_') or line[3].startswith('ENST')) and \
                '..' not in line[3] and \
                '_NG_' not in line[3] and \
                "~" not in line[3]:

            # Filter for only requested genome build
            if genome_build is not None:
                chr_num = seq_data.to_chr_num_refseq(line[4], genome_build)
                if chr_num is None:
                    continue

            # Transcript ID
            tx = line[3]

            # Protein id
            prot_id = g2t.hdp.get_pro_ac_for_tx_ac(tx)

            # Get additional tx_ information
            try:
                tx_exons = g2t.hdp.get_tx_exons(tx, line[4], line[5])
            except vvhgvs.exceptions.HGVSDataNotAvailableError:
                continue
            tx_orientation = tx_exons[0]['alt_strand']

            # Fetch the sequence details, length is the first item
            tx_anno = g2t.hdp.get_tx_seq_anno(tx)
            tx_len = tx_anno[0]

            # Exon Set
            # Collect genomic span for the transcript against known genomic/gene reference sequences
            gen_start_pos = None
            gen_end_pos = None
            exon_set = []
            # get total exons
            total_exons = len(tx_exons)
            # Set exon counter for current exon
            if tx_orientation == 1:
                current_exon_number = 0
            else:
                current_exon_number = total_exons + 1
            for tx_pos in tx_exons:
                if tx_orientation == 1:
                    current_exon_number = current_exon_number + 1
                else:
                    current_exon_number = current_exon_number - 1
                # Collect the exon_set information
                """
                tx_exons have the following attributes::
                            {
                                'tes_exon_set_id' : 98390
                                'aes_exon_set_id' : 298679
                                'tx_ac'           : 'NM_199425.2'
                                'alt_ac'          : 'NC_000020.10'
                                'alt_strand'      : -1
                                'alt_aln_method'  : 'splign'
                                'ord'             : 2
                                'tx_exon_id'      : 936834
                                'alt_exon_id'     : 2999028
                                'tx_start_i'      : 786
                                'tx_end_i'        : 1196
                                'alt_start_i'     : 25059178
                                'alt_end_i'       : 25059588
                                'cigar'           : '410='
                            }                    
                """
                current_exon = {"transcript_start": tx_pos['tx_start_i'] + 1,
                                "transcript_end": tx_pos['tx_end_i'],
                                "genomic_start": tx_pos['alt_start_i'] + 1,
                                "genomic_end": tx_pos['alt_end_i'],
                                "cigar": tx_pos['cigar'],
                                "exon_number": current_exon_number
                                }
                exon_set.append(current_exon)
                start_pos = tx_pos['alt_start_i']
                end_pos = tx_pos['alt_end_i']
                if gen_start_pos is None:
                    gen_start_pos = start_pos
                else:
                    if int(start_pos) < int(gen_start_pos):
                        gen_start_pos = int(start_pos)
                if gen_end_pos is None:
                    gen_end_pos = end_pos
                else:
                    if int(end_pos) > int(gen_end_pos):
                        gen_end_pos = int(end_pos)

            # reverse the exon_set to maintain gene and not genome orientation if gene is -1 orientated
            if tx_orientation == -1:
                exon_set.reverse()
            if bypass_genomic_spans is True:
                gen_span = False
            elif ('NG_' in line[4] or 'NC_0' in line[4]) and line[5] != 'blat':
                gen_span = True
            else:
                gen_span = False

            tx_description = g2t.db.get_transcript_description(tx)

            if tx_description == 'none':
                try:
                    g2t.db.update_transcript_info_record(tx, g2t)
                except fn.DatabaseConnectionError as e:
                    error = 'Currently unable to update gene_ids or transcript information records because ' \
                            'VariantValidator %s' % str(e)
                    # my_variant.warnings.append(error)
                    logger.warning(error)
                tx_description = g2t.db.get_transcript_description(tx)

            # Get annotation
            try:
                tx_annotation = g2t.db.get_transcript_annotation(tx)
                tx_annotation = json.loads(tx_annotation)
            # Missing annotation data
            except json.decoder.JSONDecodeError:
                continue

            # Check for duplicates
            if tx not in recovered:
                recovered.append(tx)
                if len(line) >= 3 and isinstance(line[1], int):
                    genes_and_tx.append({'reference': tx,
                                         'description': tx_description,
                                         'annotations': tx_annotation,
                                         'translation': prot_id,
                                         'length': tx_len,
                                         'coding_start': line[1] + 1,
                                         'coding_end': line[2],
                                         # 'orientation': tx_orientation,
                                         'genomic_spans': {}
                                         })
                else:
                    genes_and_tx.append({'reference': tx,
                                         'description': tx_description,
                                         'annotations': tx_annotation,
                                         'translation': prot_id,
                                         'length': tx_len,
                                         'coding_start': None,
                                         'coding_end': None,
                                         # 'orientation': tx_orientation,
                                         'genomic_spans': {}
                                         })
                # LRG information
                lrg_transcript = g2t.db.get_lrg_transcript_id_from_refseq_transcript_id(tx)
                if lrg_transcript != 'none':
                    if line[1] is None:
                        genes_and_tx.append({'reference': tx,
                                             'description': tx_description,
                                             'annotations': tx_annotation,
                                             'translation': prot_id,
                                             'length': tx_len,
                                             'coding_start': None,
                                             'coding_end': None,
                                             # 'orientation': tx_orientation,
                                             'genomic_spans': {}
                                             })
                    elif sel_tx_lst is False:
                        genes_and_tx.append({'reference': lrg_transcript,
                                             'description': tx_description,
                                             'annotations': tx_annotation,
                                             'length': tx_len,
                                             'translation': lrg_transcript.replace('t', 'p'),
                                             'coding_start': line[1] + 1,
                                             'coding_end': line[2],
                                             'genomic_spans': {}
                                             })
                    else:
                        if lrg_transcript in sel_tx_lst:
                            genes_and_tx.append({'reference': lrg_transcript,
                                                 'description': tx_description,
                                                 'annotations': tx_annotation,
                                                 'length': tx_len,
                                                 'translation': lrg_transcript.replace('t', 'p'),
                                                 'coding_start': line[1] + 1,
                                                 'coding_end': line[2],
                                                 'genomic_spans': {}
                                                 })

            # Add the genomic span information
            if gen_span is True:
                for check_tx in genes_and_tx:
                    lrg_transcript = g2t.db.get_lrg_transcript_id_from_refseq_transcript_id(tx)
                    if check_tx['reference'] == tx:
                        if gen_start_pos < gen_end_pos:
                            check_tx['genomic_spans'][line[4]] = {'start_position': gen_start_pos + 1,
                                                                  'end_position': gen_end_pos,
                                                                  'orientation': tx_orientation,
                                                                  'exon_structure': exon_set,
                                                                  "total_exons": total_exons}
                        else:
                            check_tx['genomic_spans'][line[4]] = {'start_position': gen_end_pos + 1,
                                                                  'end_position': gen_start_pos,
                                                                  'orientation': tx_orientation,
                                                                  'exon_structure': exon_set,
                                                                  "total_exons": total_exons}
                    if lrg_transcript != 'none':
                        if check_tx['reference'] == lrg_transcript:
                            if 'NG_' in line[4]:
                                lrg_id = g2t.db.get_lrg_id_from_refseq_gene_id(line[4])
                                if lrg_id[0] in lrg_transcript:
                                    check_tx['genomic_spans'][line[4]] = {'start_position': gen_start_pos + 1,
                                                                          'end_position': gen_end_pos,
                                                                          'orientation': 1,
                                                                          'exon_structure': exon_set,
                                                                          "total_exons": total_exons}

                                    check_tx['genomic_spans'][lrg_id[0]] = {'start_position': gen_start_pos + 1,
                                                                            'end_position': gen_end_pos,
                                                                            'orientation': 1,
                                                                            'exon_structure': exon_set,
                                                                            "total_exons": total_exons}

    # Return data dict
    if bypass_web_searches is True:
        g2d_data = {'transcripts': genes_and_tx}
    else:
        g2d_data = {'current_symbol': current_sym,
                    'previous_symbol': previous_sym,
                    'current_name': gene_name,
                    'lovd_corrections': lovd_corrections,
                    'lovd_messages': lovd_messages,
                    'hgnc': hgnc_id,
                    'transcripts': genes_and_tx
                    }
    g2d_data["requested_symbol"] = submitted

    return g2d_data

def lovd_syntax_check_g2t(query, lovd_syntax_check):
    # Get additional warnings
    check_with_lovd = lovd_api.lovd_syntax_check(query, do_lovd_check=lovd_syntax_check, is_a_gene=True)

    if "lovd_api_error" in check_with_lovd:
        return {}, {}

    lovd_messages = {}
    lovd_corrections = {}
    data = check_with_lovd["data"][0]
    if data["warnings"]:
        for key, val in data["warnings"].items():
            lovd_messages[key] = val

    if data["errors"]:
        for key, val in data["errors"].items():
            lovd_messages[key] = val

    if data['corrected_values']:
        for key, val in data['corrected_values'].items():
            lovd_corrections[key] = val

    lovd_messages["ISOURCE"] = check_with_lovd['url']
    lovd_messages["LIBRARYVERSION"] = check_with_lovd['version']

    return lovd_messages, lovd_corrections

# <LICENSE>
# Copyright (C) 2016-2025 VariantValidator Contributors
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
